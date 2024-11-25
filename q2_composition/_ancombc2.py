# ----------------------------------------------------------------------------
# Copyright (c) 2016-2024, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import biom
import formulaic
from formulaic.parser.types import Token
from numpy import c_, hypot
import pandas as pd
from rpy2.robjects.conversion import Converter
import rpy2.robjects.conversion as conversion
from rpy2.robjects.packages import importr
from rpy2.robjects.methods import RS4
from rpy2.robjects import pandas2ri, default_converter
import rpy2.robjects as ro
from rpy2.rinterface import NULL as RNULL

import qiime2
from qiime2.metadata import NumericMetadataColumn, CategoricalMetadataColumn

from q2_composition._format import (
    ANCOMBC2ModelStatistics, ANCOMBC2StructuralZeros, ANCOMBC2OutputDirFmt
)

r_base = importr('base')
r_phyloseq = importr('phyloseq')
r_ancombc2 = importr('ANCOMBC')


def ancombc2(
    table: biom.Table,
    metadata: qiime2.Metadata,
    fixed_effects_formula: str,
    random_effects_formula: str | None = None,
    reference_levels: list[str] | None = None,
    p_adjust_method: str = 'holm',
    prevalence_cutoff: float = 0.1,
    group: str | None = None,
    structural_zeros: bool = False,
    asymptotic_cutoff: bool = False,
    alpha: float = 0.05,
    num_processes: int = 1,
) -> ANCOMBC2OutputDirFmt:

    if structural_zeros and group is None:
        msg = (
           'The structurual zeros option was enabled but no group variable '
           'was provided. Strucutural zeros are detected according to some '
           'grouping, so a group variable must be provided.'
        )
        raise ValueError(msg)

    # process formulae
    fixed_effects_formula = _process_formula(fixed_effects_formula, metadata)
    if random_effects_formula is not None:
        random_effects_formula = _process_formula(
            random_effects_formula, metadata
        )

    # construct phyloseq object
    r_phyloseq_object = _create_phyloseq_object(
        table, metadata, reference_levels
    )

    # call ancombc2
    if group is not None:
        group = r_base.make_names(group)[0]

    with conversion.localconverter(default_converter + _get_none_converter()):
        output = r_ancombc2.ancombc2(
            data=r_phyloseq_object,
            fix_formula=fixed_effects_formula,
            rand_formula=random_effects_formula,
            p_adj_method=p_adjust_method,
            prv_cut=prevalence_cutoff,
            group=group,
            struc_zero=structural_zeros,
            alpha=alpha,
            neg_lb=asymptotic_cutoff,
            n_cl=num_processes,
            verbose=True,
        )

    # get data of interest from returned R list, put in output format
    output_format = ANCOMBC2OutputDirFmt()

    model_statistics = output[output.names.index('res')]
    with (ro.default_converter + pandas2ri.converter).context():
        model_statistics_df = ro.conversion.get_conversion().rpy2py(
            model_statistics
        )
        output_format.statistics.write_data(model_statistics_df, pd.DataFrame)

    structural_zeros = output[output.names.index('zero_ind')]
    if structural_zeros != RNULL:
        with (ro.default_converter + pandas2ri.converter).context():
            structural_zeros_df = ro.conversion.get_conversion().rpy2py(
                structural_zeros
            )
            output_format.structural_zeros.write_data(
                structural_zeros_df, pd.DataFrame
            )

    return output_format


def _process_formula(formula: str, metadata: qiime2.Metadata) -> str:
    '''
    Processes and validates `formula`. Variables in `formula` are ensured to
    be present in `metadata` and are converted to valid R identifiers as
    needed.

    Parameters
    ----------
    formula : str
        The formula to be processed.
    metadata : qiime2.Metadata
        The per-sample metadata.

    Returns
    -------
    str
        The validated formula with variables renamed to valid R identifiers
        as needed.

    Raises
    ------
    ValueError
        If an unexpected token type is encountered while parsing `formula`.
    '''
    # handle hyphens in formula
    formula, renamed_variables = _handle_hyphens(formula, metadata)

    # parse formula into tokens
    parser = formulaic.parser.parser.DefaultFormulaParser(
        include_intercept=False
    )
    tokens = list(parser.get_tokens(formula))

    # validate formula
    _validate_formula(tokens, metadata, renamed_variables)

    # convert names
    renamed_tokens = []
    for token in tokens:
        if token.kind == Token.Kind.NAME:
            r_identifier = r_base.make_names(str(token))[0]
            renamed_tokens.append(str(r_identifier))
        else:
            renamed_tokens.append(str(token))

    # reconstruct formula
    # NOTE: whitespace may differ from the original formula, but this should
    # not affect validity or semantics
    processed_formula = ' '.join(renamed_tokens)

    return processed_formula


def _handle_hyphens(
    formula: str,
    metadata: qiime2.Metadata
) -> tuple[str, list[str]]:
    '''
    Handle hyphens that may be present in metadata variables included in
    `formula`. Replace all hyphens with periods, and track renamed variables
    so that they can be searched against metadata properly later on. Note that
    only variables present in `metadata` are searched for and corrected in
    `formula`.

    Parameters
    ----------
    formula : str
        The formula string.
    metadata : qiime2.Metadata
        The per-sample metadata.

    Returns
    -------
    tuple[str, list[str]]
        A tuple containing the updated formula, and the list of renamed
        variables.
    '''
    renamed_variables: list[str] = []
    all_terms: list[str] = []

    for column in metadata.columns:
        if '-' in column and column in formula:
            renamed_variables.append(column)
            renamed_variable = column.replace('-', '.')
            formula = formula.replace(column, renamed_variable)

    return formula, renamed_variables


def _validate_formula(
    tokens: list[Token],
    metadata: qiime2.Metadata,
    renamed_variables: list[str]
) -> None:
    '''
    Asserts that the formula variables in `tokens` are present in the
    metadata. Also ensures that an independent variable is not present in the
    formula.

    Parameters
    ----------
    tokens : list[Token]
        A list of token parsed from the formula.
    metadata : qiime2.Metadata
        The per-sample metadata.
    renamed_terms : list[str]
        Terms that have had hyphens replaced with periods. Tracking these
        allows us to search for the original term in the metadata.

    Raises
    ------
    ValueError
        If one of the variables in `variables` is not a column in `metadata`.
    ValueError
        If an independent variable is specified in the formula.
    '''
    if '~' in tokens:
        msg = (
            f'A "~" symbol was detected in your formula. The '
            'independent variable should not be included in the formula.'
        )
        raise ValueError(msg)

    variables = [t for t in tokens if t.kind == Token.Kind.NAME]
    unique_variables = list(set(variables))
    for variable in unique_variables:
        variable = str(variable)

        # check for original variable name if it had its hyphens replaced
        if variable.replace('.', '-') in renamed_variables:
            variable = variable.replace('.', '-')

        try:
            metadata.get_column(variable)
        except ValueError:
            msg = (
                f'The "{variable}" variable was not found in your metadata. '
            )
            raise ValueError(msg)


def _create_phyloseq_object(
    table: biom.Table,
    metadata: qiime2.Metadata,
    reference_levels: list[str] | None = None
) -> RS4:
    '''
    Create an R phyloseq object that contains the feature table and metadata.

    Parameters
    ----------
    table : biom.Table
        The feature table to be wrapped in a phyloseq object.
    metadata : qiime2.Metadata
        The metadata to be wrapped in a phyloseq object.
    reference_levels : list[str] or None
        The desried reference levels of each of the categorical metadata
        variables included in the ANCOMBC2 formula. Specified as a list of
        "column_name::column_value" where "column_value" is the desired
        reference level of the "column_name" column.

    Returns
    -------
    RS4
        An R phyloseq object containing the feature table and metadata.
    '''
    # convert feature table to R matrix with rownames=features, colnames=samples
    table_df = table.to_dataframe(dense=True)

    with (ro.default_converter + pandas2ri.converter).context():
        r_table_df = ro.conversion.get_conversion().py2rpy(table_df)

    r_matrix = r_base.as_matrix(r_table_df)

    # construct phyloseq OTU table
    r_phyloseq_otu_table = r_phyloseq.otu_table(r_matrix, taxa_are_rows=True)

    # construct phyloseq sample (meta)data
    r_metadata_df = _convert_metadata(metadata, reference_levels)
    r_phyloseq_sample_data = r_phyloseq.sample_data(r_metadata_df)

    # wrap table and metadata into phyloseq container object
    return r_phyloseq.phyloseq(r_phyloseq_otu_table, r_phyloseq_sample_data)


def _convert_metadata(
    metadata: qiime2.Metadata,
    reference_levels: list[str] | None = None
) -> RS4:
    '''
    Converts a `qiime2.Metadata` object into an R dataframe. The
    qiime2 metadata column types are converted into their R equivalents and
    categorical columns are releveled according to `reference_levels`.

    Parameters
    ----------
    metadata : qiime2.Metadata
        The sample metadata.
    reference_levels : list[str] or None
        The desried reference levels of each of the categorical metadata
        variables included in the ANCOMBC2 formula. Specified as a list of
        "column_name::column_value" where "column_value" is the desired
        reference level of the "column_name" column.

    Returns
    -------
    RS4
        The metadata in an R dataframe.
    '''
    # convert metadata to R dataframe
    df = metadata.to_dataframe()
    with (ro.default_converter + pandas2ri.converter).context():
        r_df = ro.conversion.get_conversion().py2rpy(df)

    # ensure columns are renamed to R's fancy
    r_base.colnames(r_df)[:] = r_base.make_names(r_base.colnames(r_df))

    # convert column types
    for column_name in metadata.columns:
        column = metadata.get_column(column_name)
        r_column_name = r_base.make_names(column_name)[0]
        col_index = r_df.colnames.index(r_column_name)

        if isinstance(column, CategoricalMetadataColumn):
            r_df[col_index] = r_base.as_factor(r_df.rx2[r_column_name])
        elif isinstance(column, NumericMetadataColumn):
            r_df[col_index] = r_base.as_numeric(r_df.rx2[r_column_name])
        else:
            msg = (
                f'An unrecognized metadata column type ({type(column)}) was '
                'encountered while translating the metadata into R.'
            )
            raise ValueError(msg)

    # relevel categorical metadata columns
    if reference_levels is not None:
        columns_seen = []
        for level_spec in reference_levels:
            column_name, reference_level = _extract_column_reference_level(
                level_spec
            )

            if column_name not in metadata.columns:
                msg = (
                    f'The "{column_name}" column that was specified in the '
                    'column-reference pair was not found in the metadata.'
                )
                raise ValueError(msg)

            column = metadata.get_column(column_name)
            if isinstance(column, NumericMetadataColumn):
                msg = (
                    'Can not specify a reference level for the numeric '
                    f'metadata column {column_name}.'
                )
                raise ValueError(msg)

            if column_name in columns_seen:
                msg = (
                    'Only specify a reference level for a given categorical '
                    f'column once. The "{column_name}" column was specified '
                    'multiple times.'
                )
                raise ValueError(msg)

            columns_seen.append(column_name)

            r_column_name = r_base.make_names(column_name)
            r_column = r_df.rx2(r_column_name)
            r_column = r_base.relevel(r_column, ref=reference_level)

    return r_df


def _extract_column_reference_level(
    level_specification: str
) -> tuple[str, str]:
    '''
    Extract the column name and reference level value from a
    column name/reference level pair. Such a pair is in the format
    "column_name::column_value".

    Parameters
    ----------
    level_specification : str
        The column name/column value pair indicating the desired reference
        level. Looks like `column_name::column_value`.

    Returns
    -------
    tuple[str]
        A tuple containing (column_name, reference_level) values that are
        extracted from the level specification.
    '''
    level_spec_fields = tuple(level_specification.split('::'))

    if len(level_spec_fields) < 2:
        msg = (
            f'No reference level was detected for the {level_specification} '
            'column-reference pair. Make sure to separate the column name and '
            'reference level with "::".'
        )
        raise ValueError(msg)
    elif len(level_spec_fields) > 2:
        msg = (
            'More than one reference level were detected for the '
            f'{level_specification} column-reference pair. The expected format '
            'is "column_name::reference_level". Make sure that "::" occurs '
            'only once in each column-reference pair.'
        )
        raise ValueError(msg)

    return level_spec_fields


def _get_none_converter() -> Converter:
    '''
    Creates an rpy2.robjects.conversion.Converter object that will convert
    any python `None` to an R `NULL` used within its context handler.

    Returns
    -------
    Converter
        The rpy2 converter object.
    '''
    converter = Converter("None Converter")
    converter.py2rpy.register(type(None), lambda _: RNULL)

    return converter

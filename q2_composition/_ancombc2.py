# ----------------------------------------------------------------------------
# Copyright (c) 2016-2024, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import biom
from rpy2.robjects.packages import importr
from rpy2.robjects.methods import RS4
from rpy2.robjects import pandas2ri
import rpy2.robjects as ro

from typing import Optional

import qiime2
from qiime2.metadata import NumericMetadataColumn, CategoricalMetadataColumn

r_base = importr('base')
r_phyloseq = importr('phyloseq')


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

    # convert column types
    for column_name in metadata.columns:
        column = metadata.get_column(column_name)
        r_column_name = r_base.make_names(column_name)

        if isinstance(column, CategoricalMetadataColumn):
            r_column = r_df.rx2(r_column_name)
            r_column = r_base.as_factor(r_df.rx2(r_column_name))
        elif isinstance(column, NumericMetadataColumn):
            r_column = r_df.rx2(r_column_name)
            r_column = r_base.as_numeric(r_df.rx2(r_column_name))
        else:
            msg = (
                f'An unrecognized metadata column type ({type(column)}) was '
                'encountered while translating the metadata into R.'
            )
            raise ValueError(msg)

    # relevel
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

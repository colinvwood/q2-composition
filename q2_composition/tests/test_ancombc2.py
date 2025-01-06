# ----------------------------------------------------------------------------
# Copyright (c) 2016-2024, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
from biom import load_table
import pandas as pd
from pandas.testing import assert_frame_equal

from pathlib import Path
import unittest

import qiime2
from qiime2.metadata import NumericMetadataColumn, CategoricalMetadataColumn
from qiime2.plugin.util import transform

from q2_composition._ancombc2 import (
    r_base, ancombc2, _process_formula, _convert_metadata, _split_into_slices
)
from q2_composition._format import ANCOMBC2SliceMapping


class TestANCOMBC2Base(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.test_data_fp = Path(__file__).parent / 'data' / 'ancombc2'

        table_fp = cls.test_data_fp / 'feature-table.biom'
        cls.biom_table = load_table(table_fp)

        metadata_fp = cls.test_data_fp / 'metadata.tsv'
        cls.metadata = qiime2.Metadata.load(metadata_fp)


class TestANCOMBC2(TestANCOMBC2Base):
    def _slices_to_single_df(
        self, slices: ANCOMBC2SliceMapping
    ) -> pd.DataFrame:
        df = pd.DataFrame()
        for slice_name, slice_df in slices.items():
            if slice_name == 'structural_zeros':
                continue

            slice_df = slice_df.rename(
                lambda n: n if n == 'taxon' else f'{slice_name}_{n}',
                axis='columns'
            )

            if df.empty:
                df = slice_df
            else:
                df = pd.merge(df, slice_df, on='taxon', how='inner')

        return df

    def test_wrapped_ancombc2(self):
        '''
        Assert that ancombc2 called through qiime2 results in the same output
        as when it is called in R. The `r-model-statistics.tsv` and
        `r-structural-zeros.tsv` files were obtained by running ANCOMBC2 in R
        using the moving pictures tutorial data.
        '''
        model_stats_fp = self.test_data_fp / 'r-model-statistics.tsv'
        ground_truth_model_stats = pd.read_csv(model_stats_fp, sep='\t')
        structural_zeros_fp = self.test_data_fp / 'r-structural-zeros.tsv'
        ground_truth_struc_zeros = pd.read_csv(structural_zeros_fp, sep='\t')

        output_format = ancombc2(
            table=self.biom_table,
            metadata=self.metadata,
            fixed_effects_formula='body-site + year',
            group='body-site',
            structural_zeros=True
        )

        slices = transform(data=output_format, to_type=ANCOMBC2SliceMapping)
        model_stats = self._slices_to_single_df(slices)

        struc_zeros = output_format.structural_zeros.view(pd.DataFrame)

        print(ground_truth_model_stats.columns.difference(model_stats.columns))
        print(model_stats.columns.difference(ground_truth_model_stats.columns))
        assert_frame_equal(ground_truth_model_stats, model_stats)
        assert_frame_equal(ground_truth_struc_zeros, struc_zeros)

    def test_group_enforced_if_structural_zeros(self):
        '''
        Detecting structural zeros is only possible if a group variable is
        provided, so ensure that we error if the structural zeros flag is set
        but no group variable is passed.
        '''
        with self.assertRaisesRegex(
            ValueError, 'structural zeros option was enabled but no group'
        ):
            ancombc2(
                table=self.biom_table,
                fixed_effects_formula='body-site',
                metadata=self.metadata,
                structural_zeros=True
            )


class TestFormulaProcessing(TestANCOMBC2Base):
    def test_process_formula(self):
        '''
        Tests that invalid R identifiers in the formula are renamed to valid
        ones, that valid R identifiers are not changed, and that hyphens are
        replaced only if they are in a metadata variable.
        '''
        formula = 'year + month + day'
        processed_formula = _process_formula(formula, self.metadata)
        self.assertEqual(processed_formula, formula)

        formula = 'body-site'
        processed_formula = _process_formula(formula, self.metadata)
        self.assertEqual(processed_formula, 'body.site')

        formula = 'body-site + reported-antibiotic-usage'
        processed_formula = _process_formula(formula, self.metadata)
        self.assertEqual(
            processed_formula, 'body.site + reported.antibiotic.usage'
        )

        formula = 'year - month'
        processed_formula = _process_formula(formula, self.metadata)
        self.assertEqual(processed_formula, formula)

    def test_process_formula_errors(self):
        '''
        Tests various error conditions that should be raised when processing
        a formula. These include the presence of the dependent variable in the
        formula, variables that are not present in the metadata, and empty
        formulas.
        '''
        formula = 'abundance ~ body-site'
        with self.assertRaisesRegex(
            ValueError,
            'dependent variable should not be included in the formula'
        ):
            _process_formula(formula, self.metadata)

        formula = 'body-site + stonks'
        with self.assertRaisesRegex(
            ValueError, '"stonks" variable was not found in your metadata'
        ):
            _process_formula(formula, self.metadata)

        # NOTE: `does-not-exist` is tokenized into `[does, -, not, -, exist]`
        # and the error message will include only one of these identifiers
        # which is less than ideal, but there isn't much we can do about this
        formula = 'year + does-not-exist'
        with self.assertRaisesRegex(
            ValueError,
            r'"(does|not|exist)" variable was not found in your metadata'
        ):
            _process_formula(formula, self.metadata)

        formula = ''
        with self.assertRaisesRegex(
            ValueError, 'A formula was found to be empty'
        ):
            _process_formula(formula, self.metadata)


class TestMetadataConversion(TestANCOMBC2Base):
    def test_variable_releveling(self):
        '''
        Test the categorical variable reference levels are being enforced in
        the R metadata dataframe.
        '''
        reference_levels = ['body-site::tongue']
        releveled_metadata = _convert_metadata(self.metadata, reference_levels)
        levels = r_base.levels(releveled_metadata.rx2('body.site'))
        reference_level = levels[0]
        self.assertEqual(reference_level, 'tongue')

        reference_levels = [
            'body-site::left palm', 'reported-antibiotic-usage::Yes'
        ]
        releveled_metadata = _convert_metadata(self.metadata, reference_levels)

        body_site_levels = r_base.levels(releveled_metadata.rx2('body.site'))
        body_site_reference = body_site_levels[0]
        self.assertEqual(body_site_reference, 'left palm')

        antibiotic_levels = r_base.levels(
            releveled_metadata.rx2('reported.antibiotic.usage')
        )
        antibiotic_reference = antibiotic_levels[0]
        self.assertEqual(antibiotic_reference, 'Yes')

    def test_variable_releveling_errors(self):
        '''
        Test various error conditions that should be raised when converting
        the metadata into an R dataframe. These include reference level
        specifications of variables not found in the metadata, reference level
        specifications of numeric metadata variables, multiple reference level
        specifications for the same column, and improperly formed reference
        level specifications.
        '''
        reference_levels = ['body-site::tongue', 'fake::reference']
        with self.assertRaisesRegex(
            ValueError, 'The "fake" column.*was not found in the metadata'
        ):
            _convert_metadata(self.metadata, reference_levels)

        reference_levels = ['days-since-experiment-start::0']
        with self.assertRaisesRegex(
            ValueError,
            'Can not specify a reference level for the numeric metadata column'
        ):
            _convert_metadata(self.metadata, reference_levels)

        reference_levels = ['body-site::tongue', 'body-site::gut']
        with self.assertRaisesRegex(
            ValueError,
            'Only specify a reference level.*once'
        ):
            _convert_metadata(self.metadata, reference_levels)

        reference_levels = ['body-site::tongue::gut']
        with self.assertRaisesRegex(
            ValueError, 'More than one reference level were detected'
        ):
            _convert_metadata(self.metadata, reference_levels)

        reference_levels = ['body-site:gut']
        with self.assertRaisesRegex(
            ValueError, 'No reference level was detected'
        ):
            _convert_metadata(self.metadata, reference_levels)

        reference_levels = ['body-site']
        with self.assertRaisesRegex(
            ValueError, 'No reference level was detected'
        ):
            _convert_metadata(self.metadata, reference_levels)

    def test_variable_type_conversion(self):
        '''
        Tests that metadata variable types are properly converted to the
        corresponding types in R.
        '''
        converted_md = _convert_metadata(self.metadata)

        for column_name in self.metadata.columns:
            r_column_name = r_base.make_names(column_name)[0]
            column = self.metadata.get_column(column_name)

            if isinstance(column, CategoricalMetadataColumn):
                self.assertTrue(
                    r_base.is_factor(converted_md.rx2[r_column_name])[0]
                )
            elif isinstance(column, NumericMetadataColumn):
                self.assertTrue(
                    r_base.is_numeric(converted_md.rx2[r_column_name])[0]
                )


class TestANCOMBC2Helpers(TestANCOMBC2Base):
    def test_split_into_slices(self):
        df = pd.DataFrame({
            'taxon': ['feature1', 'feature2', 'feature3'],
            'lfc_variable.1': [0.2, 0.9, 0.1],
            'lfc_variable.2': [-0.4, 0.0, -0.3],
            'se_variable.1': [0.4, 0.02, 0.1],
            'se_variable.2': [0.04, 0.0, 0.3],
        })

        exp = ANCOMBC2SliceMapping(
            lfc=pd.DataFrame({
                'taxon': ['feature1', 'feature2', 'feature3'],
                'variable.1': [0.2, 0.9, 0.1],
                'variable.2': [-0.4, 0.0, -0.3],
            }),
            se=pd.DataFrame({
                'taxon': ['feature1', 'feature2', 'feature3'],
                'variable.1': [0.4, 0.02, 0.1],
                'variable.2': [0.04, 0.0, 0.3],
            })
        )

        obs = _split_into_slices(df)

        assert_frame_equal(exp['lfc'], obs['lfc'])
        assert_frame_equal(exp['se'], obs['se'])

    def test_split_into_slices_overlapping_prefixes(self):
        df = pd.DataFrame({
            'taxon': ['feature1', 'feature2', 'feature3'],
            'p_variable.1': [0.2, 0.9, 0.1],
            'p_variable.2': [-0.4, 0.0, -0.3],
            'passed_ss_variable.1': [0.4, 0.02, 0.1],
            'passed_ss_variable.2': [0.04, 0.0, 0.3],
        })

        exp = ANCOMBC2SliceMapping(
            p=pd.DataFrame({
                'taxon': ['feature1', 'feature2', 'feature3'],
                'variable.1': [0.2, 0.9, 0.1],
                'variable.2': [-0.4, 0.0, -0.3],
            }),
            passed_ss=pd.DataFrame({
                'taxon': ['feature1', 'feature2', 'feature3'],
                'variable.1': [0.4, 0.02, 0.1],
                'variable.2': [0.04, 0.0, 0.3],
            })
        )

        obs = _split_into_slices(df)

        assert_frame_equal(exp['p'], obs['p'])
        assert_frame_equal(exp['passed_ss'], obs['passed_ss'])

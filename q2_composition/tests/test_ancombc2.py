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

import qiime2
from qiime2.plugin.testing import TestPluginBase

from q2_composition._ancombc2 import ancombc2

class TestANCOMBC2(TestPluginBase):
    package = 'q2_composition.tests'

    def setUp(self):
        self.test_data_fp = Path(__file__).parent / 'data' / 'ancombc2'

        table_fp = self.test_data_fp / 'feature-table.biom'
        self.biom_table = load_table(table_fp)

        metadata_fp = self.test_data_fp / 'metadata.tsv'
        self.metadata = qiime2.Metadata.load(metadata_fp)

    def tearDown(self):
        pass

    def test_wrapped_ancombc2(self):
        '''
        Assert that ancombc2 called through qiime2 results in the same output
        as when it is called in R.
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

        model_stats = output_format.statistics.view(pd.DataFrame)
        struc_zeros = output_format.structural_zeros.view(pd.DataFrame)

        print('gt struc zeros', ground_truth_struc_zeros.info())
        print('struc zeros', struc_zeros.info())
        assert_frame_equal(ground_truth_model_stats, model_stats)
        assert_frame_equal(ground_truth_struc_zeros, struc_zeros)

    def test_group_enforced_if_structural_zeros(self):
        with self.assertRaisesRegex(
            ValueError, '.*structural zeros option was enabled but no group.*'
        ):
            ancombc2(
                table=self.biom_table,
                fixed_effects_formula='body-site',
                metadata=self.metadata,
                structural_zeros=True
            )

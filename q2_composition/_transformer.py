# ----------------------------------------------------------------------------
# Copyright (c) 2016-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pandas as pd

from q2_composition.plugin_setup import plugin
from q2_composition._format import (
    FrictionlessCSVFileFormat,
    ANCOMBC2ModelStatistics,
    ANCOMBC2StructuralZeros,
)


@plugin.register_transformer
def _1(obj: FrictionlessCSVFileFormat) -> pd.DataFrame:
    path = obj.view(FrictionlessCSVFileFormat)
    df = pd.read_csv(str(path))

    return df


@plugin.register_transformer
def _2(df: pd.DataFrame) -> ANCOMBC2ModelStatistics:
    format = ANCOMBC2ModelStatistics()
    with format.open() as fh:
        df.to_csv(fh, sep='\t')

    return format


@plugin.register_transformer
def _3(format: ANCOMBC2ModelStatistics) -> pd.DataFrame:
    return pd.read_csv(format.path, sep='\t')


@plugin.register_transformer
def _4(df: pd.DataFrame) -> ANCOMBC2StructuralZeros:
    format = ANCOMBC2StructuralZeros()
    with format.open() as fh:
        df.to_csv(fh, sep='\t')

    return format


@plugin.register_transformer
def _5(format: ANCOMBC2StructuralZeros) -> pd.DataFrame:
    return pd.read_csv(format.path, sep='\t')

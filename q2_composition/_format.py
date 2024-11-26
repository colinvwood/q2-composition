# ----------------------------------------------------------------------------
# Copyright (c) 2022-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from qiime2.plugin import model


class FrictionlessCSVFileFormat(model.TextFileFormat):
    """
    Format for frictionless CSV.

    More on this later.
    """
    def _validate_(self, level):
        pass


class DataPackageSchemaFileFormat(model.TextFileFormat):
    """Format for the associated metadata for each file in the DataLoaf.

    More on this later.
    """
    def _validate_(self, level):
        pass


class DataLoafPackageDirFmt(model.DirectoryFormat):
    data_slices = model.FileCollection(r'.+\.csv',
                                       format=FrictionlessCSVFileFormat)
    nutrition_facts = model.File('datapackage.json',
                                 format=DataPackageSchemaFileFormat)

    @data_slices.set_path_maker
    def _data_slices_path_maker(self, slice_name):
        return slice_name + '.csv'


class ANCOMBC2ModelStatistics(model.TextFileFormat):
    '''
    Stores the primary output table of the ANCOMBC2 method which contains
    log fold change estimates and their standard errors for the variables
    included in the mixed effects model.
    '''
    def _validate_(self, level):
        pass


class ANCOMBC2StructuralZeros(model.TextFileFormat):
    '''
    Stores the structural zeros output table of the ANCOMBC2 method.
    '''
    def _validate_(self, level):
        pass


class ANCOMBC2OutputDirFmt(model.DirectoryFormat):
    '''
    Stores the model statistics and structural zeros tables that are output
    by the ANCOMBC2 method.
    '''
    statistics = model.File(
        'ANCOMBC2-statistics.tsv',
        format=ANCOMBC2ModelStatistics
    )
    structural_zeros = model.File(
        'ANCOMBC2-structurual-zeros.tsv',
        format=ANCOMBC2StructuralZeros,
        optional=True
    )

    def _validate_(self, level):
        pass

# **************************************************************************
# *
# * Authors:    Alberto M. Parra Pérez (amparraperez@gmail.com)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

from pyworkflow.tests import *
from .test_import_small_Molecules import TestImportSmallMolecules
#from .test_list_operate import TestListOperate
from .test_export_csv import TestExportcsv
from .test_ligand_preparation import TestLigandPreparation
from .test_converter import TestConverter
from .test_protocol_consensus_docking import TestConsensusDocking
from .test_protocol_define_pockets import TestDefinePockets
from .test_protocol_score_docking import TestScoreDocking
from .test_protocol_import_setOfSequences import TestImportSequences
from .test_protocol_multiple_alignment import TestAlignSequences
from .test_protocol_import_variants import TestImportVariants
from .test_protocol_generate_variants import TestGenerateVariants
from .test_protocol_pairwise_alignment import TestPairAlignSequences
from .test_protocol_extract_sequence_rois import TestExtractROIs
from .test_protocol_map_sequence_structure_rois import TestMapROIs
from .test_protocol_define_sequence_rois import TestDefineROIs

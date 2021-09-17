# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors: Daniel Del Hoyo (ddelhoyo@cnb.csic.es)
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
# *  e-mail address 'you@yourinstitution.email'
# *
# **************************************************************************


"""
This protocol is used to import a set of pockets (of fpocket, p2rank, autoligand) from some files

"""
from pyworkflow.protocol import params
from pwem.protocols import EMProtocol
from pyworkflow.utils import Message, weakImport
from pyworkflow.utils.path import copyFile
from pwchem.objects import SetOfPockets
from pwem.objects.data import AtomStruct
from ..constants import *
from pwchem.utils import *
import os, glob

try:
  from fpocket.objects import FpocketPocket
except:
  print('Could not find fpocket plugin')
try:
    from p2rank.objects import P2RankPocket
except:
    print('Could not find p2rank plugin')
try:
    from autodock.objects import AutoLigandPocket
except:
    print('Could not find autodock plugin')


FPOCKET, P2RANK, AUTOLIGAND = 0, 1, 2

class ImportPockets(EMProtocol):
    """
    Executes the import pockets
    """
    _label = 'import pockets'
    typeChoices = ['FPocket', 'P2Rank', 'AutoLigand']

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ """
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('pocketType', params.EnumParam, default=FPOCKET,
                      label='Type of pocket', choices=self.typeChoices,
                      help='Software used to extract the described pockets')
        form.addParam('proteinFile', params.PathParam,
                      allowsNull=False, label="AtomStruct file where pockets are described",
                      help='Select the pdb/cif file where the pockets are described')
        form.addParam('filesPath', params.PathParam,
                      label="Files directory",
                      help="Directory with the files you want to import.\n\n"
                           "The path can also contain wildcards to select"
                           "from several folders. \n\n"
                           "Examples:\n"
                           "  ~/project/data/day??_files/\n"
                           "Each '?' represents one unknown character\n\n"
                           "  ~/project/data/day*_files/\n"
                           "'*' represents any number of unknown characters\n\n"
                           "  ~/project/data/day##_files/\n"
                           "'##' represents two digits that will be used as "
                           "file ID\n\n"
                           "NOTE: wildcard characters ('*', '?', '#') "
                           "cannot appear in the actual path.)")
        form.addParam('filesPattern', params.StringParam,
                      label='Pocket files pattern', default="*.pdb",
                      help="Pattern of the files to be imported.\n\n"
                           "The pattern can contain standard wildcards such as\n"
                           "*, ?, etc, or special ones like ### to mark some\n"
                           "digits in the filename as ID.\n\n"
                           "NOTE: wildcards and special characters "
                           "('*', '?', '#', ':', '%') cannot appear in the "
                           "actual path.\n\n"
                           "You may create pockets from '.pdb' (fpocket, autoligand), .csv (p2rank)")

        form.addParam('pqrPattern', params.StringParam,
                      label='PQR pattern', default="*.pqr", condition='pocketType=={}'.format(FPOCKET),
                      help="Pattern of the files to be imported.\n\n"
                           "The pattern can contain standard wildcards such as\n"
                           "*, ?, etc, or special ones like ### to mark some\n"
                           "digits in the filename as ID.\n\n"
                           "NOTE: wildcards and special characters "
                           "('*', '?', '#', ':', '%') cannot appear in the "
                           "actual path.\n\n"
                           "The .pqr files contain parameters like the radius of the fpocket alpha-spheres")

        form.addParam('resultsFile', params.PathParam, condition='pocketType=={}'.format(AUTOLIGAND),
                      allowsNull=False, label="Results file from autoligand",
                      help='Results file from autoligand with their parameters (common for all pockets)')



    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        # Insert processing steps
        self._insertFunctionStep('importStep')

    def importStep(self):
        proteinFile = self.proteinFile.get()
        pocketFileNames = self.getMatchingFiles(self.filesPath.get(), self.filesPattern.get())
        pocketFileNames = self.copyFiles(pocketFileNames, self._getExtraPath())

        outputPockets = SetOfPockets(filename=self._getPath('pockets.sqlite'))
        if self.pocketType.get() == FPOCKET:
            print('1')
            pqrFileNames = self.getMatchingFiles(self.filesPath.get(), self.pqrPattern.get())
            pqrFileNames = self.copyFiles(pqrFileNames, self._getExtraPath())
            for pocketFn in pocketFileNames:
                print(pocketFn)
                pqrFn = pocketFn.replace('atm.pdb', 'vert.pqr')
                print(pqrFn)
                if pqrFn in pqrFileNames:
                    print(2)
                    outputPockets.append(FpocketPocket(pocketFn, proteinFile, pqrFn))

        elif self.pocketType.get() == P2RANK:
            for pocketFn in pocketFileNames:
                outputPockets.append(P2RankPocket(pocketFn, proteinFile))

        elif self.pocketType.get() == AUTOLIGAND:
            if os.path.exists(self.resultsFile.get()):
                resultsFile = os.path.abspath(self.resultsFile.get())
                for pocketFn in pocketFileNames:
                    outputPockets.append(AutoLigandPocket(pocketFn, proteinFile, resultsFile))

        self._defineOutputs(outputPockets=outputPockets)

    # --------------------------- INFO functions -----------------------------------
    def _summary(self):
        summary = []
        return summary

    def _methods(self):
        methods = []
        return methods

    def _validate(self):
        errors = []
        if self.pocketType.get() == FPOCKET:
            if not tryImportFpocket():
              errors.append('FPocket plugin not found')

        elif self.pocketType.get() == P2RANK:
            if not tryImportP2Rank():
              errors.append('P2Rank plugin not found')

        elif self.pocketType.get() == AUTOLIGAND:
            if not tryImportAutoligand():
              errors.append('Autodock plugin not found')

        return errors

    # --------------------------- UTILS functions -----------------------------------

    def getMatchingFiles(self, dir, pattern):
        fileNames = []
        for fn in glob.glob(os.path.join(dir, pattern)):
            fileNames.append(fn)
        return fileNames

    def copyFiles(self, filenames, dstDir):
        newFns = []
        dstDir = os.path.abspath(dstDir)
        for fn in filenames:
            newFn = os.path.join(dstDir, os.path.split(fn)[1])
            copyFile(fn, newFn)
            newFns.append(newFn)
        return newFns


def tryImportFpocket():
  with weakImport("fpocket"):
    from fpocket.objects import FpocketPocket
    return True
  return False

def tryImportP2Rank():
  with weakImport("p2rank"):
    from p2rank.objects import P2RankPocket
    return True
  return False

def tryImportAutoligand():
  with weakImport("autodock"):
    from autodock.objects import AutoLigandPocket
    return True
  return False


# **************************************************************************
# *
# * Authors:     Carlos Oscar Sorzano (coss@cnb.csic.es)
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

import lxml.etree as ET
import os
import sys
import urllib.request

from pwem.objects import Sequence
from pwem.protocols import EMProtocol
import pyworkflow.object as pwobj
from pyworkflow.protocol.params import *
from pwchem.objects import DatabaseID, SetOfDatabaseID, Variants, SequenceFasta, SequenceVariants

class ProtChemUniprotSequenceVariants(EMProtocol):
# class ProtChemUniprotCrossRef(EMProtocol):
    """Extract natural variants from uniprot"""
    _label = 'extract variants from uniprot'

    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('fromID', BooleanParam,
                      label='Import from Id:', default = True,
                      help="Where from Uniprot or sequence is downloaded")
        form.addParam('inputUniProtKB', StringParam, condition='fromID',
                       label='Uniprot Id:', allowsNull=False,
                       help="UniProtKB ID for the variants query")
        form.addParam('inputSequence', PointerParam, pointerClass='Sequence',
                      condition='not fromID',
                      label='Input Sequence:', allowsNull=False,
                      help="Original Sequence")

    # --------------------------- INSERT steps functions --------------------
    def _insertAllSteps(self):
        if self.fromID:
            self._insertFunctionStep('downLoadFasta')
        self._insertFunctionStep('extractStep')

    def downLoadFasta(self):
        uniprotId = self.inputUniProtKB.get()
        urlFasta = "https://www.uniprot.org/uniprot/%s.fasta" % uniprotId
        fnFasta = self._getPath("%s.fasta" % uniprotId)

        if not os.path.exists(fnFasta):
            print("Fetching uniprot: %s"%urlFasta)

            try:
                urllib.request.urlretrieve(urlFasta,fnFasta)
            except: # The library raises an exception when the web is not found
                pass

        self.downloadedSequence = Sequence(id=uniprotId)
        self.downloadedSequence.importFromFile(os.path.abspath(fnFasta))


    def extractStep(self):
        if self.fromID:
            uniprotId = self.inputUniProtKB.get()
            sequence = self.downloadedSequence.getSequence()
        else:
            uniprotId = self.inputSequence.get().getId()
            sequence = self.inputSequence.get().getSequence()

        # uniprotId = item._uniprotId.get()
        print("Processing %s"%uniprotId)

        urlId = "https://www.uniprot.org/uniprot/%s.xml" % uniprotId

        fnXML=self._getExtraPath("%s.xml"%uniprotId)

        fnAll = self._getPath("NaturalVariantsUniprot.txt")
        file = open(fnAll, "w")
        file.close()
        if not os.path.exists(fnXML):
            print("Fetching uniprot: %s"%urlId)
            # for i in range(3):
            try:
                urllib.request.urlretrieve(urlId,fnXML)
            except: # The library raises an exception when the web is not found
                pass
        if os.path.exists(fnXML):
            fnAll = self.parseXML(fnXML, fnAll)

        varsSeq = SequenceVariants(sequence=sequence, id=uniprotId)
        varsSeq.setVariantsFileName(fnAll)

        self._defineOutputs(outputVariants=varsSeq)




    def _validate(self):
        errors=[]
        return errors

    def parseXML(self, fnXML, fnOut):
        tree = ET.parse(fnXML)

        for child in tree.getroot().iter():
            positionDescription = []
            strainDescription = []
            if child.tag.endswith("feature"):

                if child.attrib['type'] == 'sequence variant':
                    for childChild in child:
                        if childChild.tag.endswith("location"):
                            # print(childChild)
                            for childChildChild in childChild:
                                if childChildChild.tag.endswith("position"):
                                    print('position: ', childChildChild.attrib['position'])
                                    positionDescription.append(childChildChild.attrib['position'])

                    if child.attrib['description']:
                        # print('description: ', child.attrib['description'])
                        strainDescription.append(child.attrib['description'])

                    for childChild in child:
                        if childChild.tag.endswith("original"):
                            # print(childChild.text)
                            positionDescription.append(childChild.text)

                    for childChild in child:
                        if childChild.tag.endswith("variation"):
                            # print(childChild.text)
                            positionDescription.append(childChild.text)
            if len(positionDescription) >= 3:
            # if len(positionDescription) != 0 and len(positionDescription) != 1:
                print('positionDescription', positionDescription)
                print('description: ', child.attrib['description'])
                descriptionLinage = child.attrib['description'].split(',')
                identifiedVariants = []
                for elementDescription in range(len(descriptionLinage)):
                    # print('element: ', descriptionLinage[elementDescription])

                    if 'In strain: ' in descriptionLinage[elementDescription]:
                        strainline = descriptionLinage[elementDescription]
                        inStrain = strainline.split(':')
                        # print('inStrain: ', inStrain[1])
                        # identifiedVariants = inStrain[1]
                        identifiedVariants.append(inStrain[1])

                    if 'strain' not in descriptionLinage[elementDescription]:
                        if ';' in descriptionLinage[elementDescription]:
                            strainline = descriptionLinage[elementDescription]
                            inStrain = strainline.split(';')
                            # print('inStrain: ', inStrain[0])
                            # identifiedVariants = identifiedVariants + inStrain[0]
                            identifiedVariants.append((inStrain[0] + '.'))
                        else:
                            inStrain = descriptionLinage[elementDescription]
                            # print('inStrain: ', inStrain)
                            # identifiedVariants = identifiedVariants + inStrain
                            identifiedVariants.append(inStrain)

                    if 'in strain' in descriptionLinage[elementDescription]:
                        strainline = descriptionLinage[elementDescription]
                        inStrain = strainline.split(' ')
                        # print('inStrain: ', inStrain[3])
                        # identifiedVariants = identifiedVariants + inStrain[3]
                        identifiedVariants.append(inStrain[3])

                # print('identifiedVariants: ', identifiedVariants)

                mutantString = identifiedVariants[0].lstrip()
                for mutantIndex in range(1, (len(identifiedVariants))):
                    mutantString = mutantString + ', ' + identifiedVariants[mutantIndex].lstrip()
                print('mutantString: ', mutantString)

                # Crear un archivo con las variantes
                file = open(fnOut, "a")
                snv = file.write(positionDescription[1] + positionDescription[0] + positionDescription[2] + ' ' + mutantString +'\n')
                # snv = file.write(positionDescription[1] + positionDescription[0] + positionDescription[2] + '\n')
                file.close()
        return fnOut

    def readFasta(self, fnFasta):
        seq = ''
        file = open(fnFasta)
        file.readline()
        for line in file:
            seq += line.strip()
        file.close()
        return seq




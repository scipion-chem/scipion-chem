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

import os
import lxml.etree as ET
import urllib.request

from pwem.objects import Sequence
from pwem.protocols import EMProtocol
from pyworkflow.protocol.params import *
from pwchem.objects import SequenceVariants


class ProtChemImportVariants(EMProtocol):
    """Extract natural variants from uniprot or associates a mutation list to a sequence """
    _label = 'Import sequence variants'

    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('fromID', BooleanParam,
                      label='Import from Uniprot Id: ', default=True,
                      help="Whether to import the sequence and variants from an UniprotID or from local files")

        form.addParam('inputUniProtKB', StringParam, condition='fromID',
                      label='Uniprot Id: ', allowsNull=False,
                      help="Uniprot entry for the sequence to import")

        form.addParam('inputSequence', PointerParam, pointerClass='Sequence',
                      condition='not fromID',
                      label='Input Sequence: ', allowsNull=False,
                      help="Input sequence to associate the variants to")

        form.addParam('inputMutaList', PathParam,
                      condition='not fromID',
                      label='Input variants file: ', allowsNull=False,
                      help="Customized list of variants, based on single substitutions. "
                           "\ne.g. leucine (L) in position 5 is replaced by phenylalanine (F): L5F. "
                           "\nVariants or lineages including several mutations might also be defined, specifying them "
                           "in a second column: L5F Iota/B.1.526 ")

    def _insertAllSteps(self):
        if self.fromID:
            self._insertFunctionStep('downLoadFasta')
        self._insertFunctionStep('extractStep')

    def downLoadFasta(self):
        uniprotId = self.inputUniProtKB.get()
        urlFasta = "https://www.uniprot.org/uniprot/%s.fasta" % uniprotId
        fnFasta = self._getPath("%s.fasta" % uniprotId)

        if not os.path.exists(fnFasta):
            print("Fetching uniprot: %s" % urlFasta)
            try:
                urllib.request.urlretrieve(urlFasta, fnFasta)
            except:  # The library raises an exception when the web is not found
                pass

        self.downloadedSequence = Sequence(id=uniprotId)
        self.downloadedSequence.importFromFile(os.path.abspath(fnFasta))

    def extractStep(self):
        if self.fromID:
            uniprotId = self.inputUniProtKB.get()
            sequence = self.downloadedSequence.getSequence()
            print("Processing %s" % uniprotId)
            urlId = "https://www.uniprot.org/uniprot/%s.xml" % uniprotId
            fnXML = self._getExtraPath("%s.xml" % uniprotId)
            fnAll = self._getPath("NaturalVariantsUniprot.txt")
            file = open(fnAll, "w")
            file.close()

            if not os.path.exists(fnXML):
                print("Fetching uniprot: %s" % urlId)
                try:
                    urllib.request.urlretrieve(urlId, fnXML)
                except:  # The library raises an exception when the web is not found
                    pass

            if os.path.exists(fnXML):
                fnAll = self.parseXML(fnXML, fnAll)
            seqObj = Sequence(sequence=sequence, id=uniprotId, name=uniprotId)
            varsSeq = SequenceVariants(filename=fnAll)
            varsSeq.setSequence(seqObj)
        else:
            uniprotId = self.inputSequence.get().getSeqName()
            sequence = self.inputSequence.get().getSequence()
            mutListCustom = self.inputMutaList.get()
            fileList = open(mutListCustom, 'r')
            typeList = ''
            for line_mutListCustom in fileList:
                line_mutListCustom_list = line_mutListCustom.rstrip().split(' ')
                length_list = len(line_mutListCustom_list)
                if length_list >= 2:
                    typeList = 'Complete'
                else:
                    typeList = 'Simple'
                break
            fileList.close()

            if typeList == 'Simple':
                fileList = open(mutListCustom, 'r')
                mutListCustomSimple = self._getPath("NaturalVariantsCustomized.txt")
                file_mutListCustomSimple = open(mutListCustomSimple, 'w')
                for line_mutListCustom in fileList:
                    newFileList = line_mutListCustom.rstrip() + ' ' + line_mutListCustom
                    file_mutListCustomSimple.write(newFileList)
                fileList.close()
                file_mutListCustomSimple.close()

                seqObj = Sequence(sequence=sequence, id=uniprotId, name=uniprotId)
                varsSeq = SequenceVariants(filename=mutListCustomSimple)
                varsSeq.setSequence(seqObj)
            else:
                seqObj = Sequence(sequence=sequence, id=uniprotId, name=uniprotId)
                varsSeq = SequenceVariants(filename=mutListCustom)
                varsSeq.setSequence(seqObj)

        self._defineOutputs(outputVariants=varsSeq)

    def _validate(self):
        errors = []
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
                            for childChildChild in childChild:
                                if childChildChild.tag.endswith("position"):
                                    positionDescription.append(childChildChild.attrib['position'])
                    if child.attrib['description']:
                        strainDescription.append(child.attrib['description'])
                    for childChild in child:
                        if childChild.tag.endswith("original"):
                            positionDescription.append(childChild.text)
                    for childChild in child:
                        if childChild.tag.endswith("variation"):
                            positionDescription.append(childChild.text)
            if len(positionDescription) >= 3:
                descriptionLinage = child.attrib['description'].split(',')
                identifiedVariants = []
                for elementDescription in range(len(descriptionLinage)):
                    if 'In strain: ' in descriptionLinage[elementDescription]:
                        strainline = descriptionLinage[elementDescription]
                        inStrain = strainline.split(':')
                        identifiedVariants.append(inStrain[1])
                    if 'strain' not in descriptionLinage[elementDescription]:
                        if ';' in descriptionLinage[elementDescription]:
                            strainline = descriptionLinage[elementDescription]
                            inStrain = strainline.split(';')
                            identifiedVariants.append((inStrain[0] + '.'))
                        else:
                            inStrain = descriptionLinage[elementDescription]
                            identifiedVariants.append(inStrain)
                    if 'in strain' in descriptionLinage[elementDescription]:
                        strainline = descriptionLinage[elementDescription]
                        inStrain = strainline.split(' ')
                        identifiedVariants.append(inStrain[3])

                mutantString = identifiedVariants[0].lstrip()
                for mutantIndex in range(1, (len(identifiedVariants))):
                    mutantString = mutantString + ', ' + identifiedVariants[mutantIndex].lstrip()

                # Make file with variants
                file = open(fnOut, "a")
                snv = file.write(positionDescription[1] + positionDescription[0] + positionDescription[
                    2] + ' ' + mutantString + '\n')
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

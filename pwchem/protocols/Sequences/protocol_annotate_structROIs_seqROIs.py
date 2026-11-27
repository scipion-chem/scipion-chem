# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors: Judith Maestro (scipionchem@cnb.csic.es)
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
This protocol annotates a set of structural ROIs with the sequence ROIs (mutations of a variant, epitopes, ...)
they contain, and reports the residues of a partner chain in contact with those sequence ROIs
"""
import json, os
from scipy import spatial

from pwem.protocols import EMProtocol
from pwem.convert.atom_struct import AtomicStructHandler
from pyworkflow.protocol import params
from pyworkflow.object import String, Integer
from pyworkflow.utils import Message

from pwchem.objects import SetOfStructROIs
from pwchem.utils import getChainIds, natural_sort
from pwchem.utils.utilsFasta import pairwiseAlign


class ProtAnnotateStructROIsWithSeqROIs(EMProtocol):
    """
    Annotates a set of structural ROIs with the sequence ROIs (e.g. the mutations of a variant) located in them.

    The sequence of the sequence ROIs is aligned to the chain of interest of the structure where the structural
    ROIs are defined, so the sequence ROIs positions are translated into the structure numbering. Each structural
    ROI is then annotated with the sequence ROIs whose residues belong to it.

    Optionally, the residues of a partner chain (e.g. an antibody or a nanobody in a complex) closer than a given
    distance to those residues are also reported.
    """
    _label = 'Annotate structural ROIs with sequence ROIs'

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_INPUT)
        group = form.addGroup('Input')
        group.addParam('inputStructROIs', params.PointerParam, pointerClass='SetOfStructROIs',
                       allowsNull=False, label="Input structural ROIs: ",
                       help='Structural ROIs (e.g. the interface clusters predicted by ARCTIC-3D) to annotate')
        group.addParam('inputSequenceROIs', params.PointerParam, pointerClass='SetOfSequenceROIs',
                       allowsNull=False, label="Input sequence ROIs: ",
                       help='Sequence ROIs (e.g. the mutations of a variant) to locate in the structural ROIs')

        group.addParam('chain_name', params.StringParam, allowsNull=False,
                       label='Chain of interest: ',
                       help='Chain of the structure of the ROIs corresponding to the input sequence')
        group.addParam('preview', params.LabelParam, condition='chain_name!=""',
                       label='Preview alignment: ',
                       help='Preview the alignment of the input sequence and the selected chain')

        group = form.addGroup('Contacts with partner')
        group.addParam('doPartner', params.BooleanParam, default=True,
                       label='Report partner residues in contact: ',
                       help='Whether to report the residues of another chain (e.g. a nanobody) in contact with '
                            'the sequence ROIs found in each structural ROI')
        group.addParam('partnerChains', params.StringParam, condition='doPartner', default='',
                       label='Partner chains: ',
                       help='Chains whose residues in contact with the sequence ROIs will be reported. '
                            'If empty, every chain but the chain of interest is used')
        group.addParam('maxDistance', params.FloatParam, condition='doPartner', default=5.0,
                       label='Maximum contact distance (A): ',
                       help='Maximum distance between any atom of a sequence ROI residue and any atom of a '
                            'partner residue for them to be considered in contact')

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep(self.alignSequencesStep)
        self._insertFunctionStep(self.annotateStep)
        self._insertFunctionStep(self.defineOutputStep)

    def alignSequencesStep(self):
        seq, seqName = self.getInputSequence()
        seqAS, seqNameAS = self.getChainSequence()

        outFile = os.path.abspath(self._getPath('pairWise.fasta'))
        pairwiseAlign(seq, seqAS, outFile, seqName1=seqName, seqName2=seqNameAS, force=True)
        pairwiseAlign(seq, seqAS, outFile.replace('.fasta', '.aln'), seqName1=seqName, seqName2=seqNameAS,
                      force=True)

    def annotateStep(self):
        """Translate every sequence ROI into the structure numbering and assign it to the structural ROIs
        containing any of its residues"""
        mapDic = self.mapResidues()
        if not mapDic:
            raise Exception(f'The input sequence could not be aligned to the chain of interest, check the '
                            f'alignment in {self._getPath("pairWise.aln")}')
        self.info(f'{len(mapDic)} residues of the input sequence mapped to the chain of interest')

        seqROIsDic = {}
        for seqROI in self.inputSequenceROIs.get():
            resIds = self.getSeqROIResidues(seqROI, mapDic)
            if resIds:
                seqROIsDic[str(seqROI.getObjId())] = {'name': self.getSeqROIName(seqROI), 'residues': resIds}
            else:
                self.info(f'{self.getSeqROIName(seqROI)} could not be mapped to the chain of interest')

        annotDic, resStructDic = {}, {}
        for structROI in self.inputStructROIs.get():
            roiId = str(structROI.getObjId())
            roiResidues = set(structROI.getDecodedCResidues())
            matched = [sId for sId in seqROIsDic if set(seqROIsDic[sId]['residues']) & roiResidues]
            annotDic[roiId] = matched
            for sId in matched:
                resStructDic.setdefault(sId, []).append(roiId)

        partnerDic = self.getPartnerResidues(annotDic, seqROIsDic) if self.doPartner.get() else {}

        self.saveAnnotation({'seqROIs': seqROIsDic, 'structROIs': annotDic, 'inStructROIs': resStructDic,
                             'partner': partnerDic})
        self.writeReport(seqROIsDic, annotDic, partnerDic)

    def defineOutputStep(self):
        annot = self.readAnnotation()
        seqROIsDic, annotDic, partnerDic = annot['seqROIs'], annot['structROIs'], annot['partner']

        outROIs = SetOfStructROIs(filename=self._getPath('StructROIs.sqlite'))
        for structROI in self.inputStructROIs.get():
            roiId = str(structROI.getObjId())
            matched = annotDic.get(roiId, [])
            names = [seqROIsDic[sId]['name'] for sId in matched]
            resIds = [resId for sId in matched for resId in seqROIsDic[sId]['residues']]

            structROI._seqROIs = String(','.join(names) if names else None)
            structROI._nSeqROIs = Integer(len(names))
            structROI._seqROIResidues = String('-'.join(natural_sort(resIds)) if resIds else None)
            if self.doPartner.get():
                partnerRes = partnerDic.get(roiId, [])
                structROI._partnerResidues = String('-'.join(partnerRes) if partnerRes else None)
                structROI._nPartnerResidues = Integer(len(partnerRes))
            outROIs.append(structROI)

        inSeqROIs = self.inputSequenceROIs.get()
        outSeqROIs = inSeqROIs.createCopy(self._getPath(), copyInfo=True)
        for seqROI in inSeqROIs:
            sId = str(seqROI.getObjId())
            seqROI._structROIs = String(','.join(map(str, annot['inStructROIs'].get(sId, []))) or None)
            seqROI._structResidues = String('-'.join(seqROIsDic[sId]['residues']) if sId in seqROIsDic else None)
            outSeqROIs.append(seqROI)

        if len(outROIs) > 0:
            outROIs.buildPDBhetatmFile()
            self._defineOutputs(outputStructROIs=outROIs, outputSequenceROIs=outSeqROIs)
            self._defineSourceRelation(self.inputStructROIs, outROIs)
            self._defineSourceRelation(self.inputSequenceROIs, outSeqROIs)

    # --------------------------- INFO functions -----------------------------------
    def _summary(self):
        summary = []
        if os.path.exists(self.getReportFile()):
            with open(self.getReportFile()) as f:
                lines = f.read().strip().split('\n')
            summary += lines[:11]
            if len(lines) > 11:
                summary.append(f'... see the whole annotation in {self.getReportFile()}')
        else:
            summary.append('The protocol has not been run yet')
        return summary

    def _validate(self):
        errors = []
        if self.chain_name.get() and self.doPartner.get() and self.partnerChains.get().strip():
            if self.getChainName() in self.getPartnerChainIds():
                errors.append('The chain of interest cannot be one of the partner chains')
        return errors

    # --------------------------- UTILS functions -----------------------------------
    def getStructFile(self):
        return os.path.abspath(self.inputStructROIs.get().getProteinFile())

    def getChainName(self):
        return json.loads(self.chain_name.get())['chain'].upper().strip()

    def getModelId(self):
        return int(json.loads(self.chain_name.get()).get('model', 0))

    def getPartnerChainIds(self):
        chainsStr = self.partnerChains.get().strip()
        return getChainIds(chainsStr) if chainsStr else []

    def getInputSequence(self):
        inputObj = self.inputSequenceROIs.get()
        seqName = inputObj.getSequenceObj().getId() or inputObj.getSeqName()
        return inputObj.getSequence(), seqName

    def getChainSequence(self):
        asFile = self.getStructFile()
        handler = AtomicStructHandler(asFile)
        seq = str(handler.getSequenceFromChain(modelID=self.getModelId(), chainID=self.getChainName()))
        return seq, os.path.basename(asFile)

    def getStructModel(self):
        return AtomicStructHandler(self.getStructFile()).getStructure()[self.getModelId()]

    def getChainResidues(self, chainId=None):
        chainId = chainId if chainId else self.getChainName()
        chain = self.getStructModel()[chainId]
        return [res for res in chain.get_residues() if res.get_id()[0] == ' ']

    def parseSequences(self):
        seq, seqAS, first = '', '', True
        with open(self._getPath('pairWise.fasta')) as f:
            f.readline()
            for line in f:
                if line.startswith('>'):
                    first = False
                elif first:
                    seq += line.strip()
                else:
                    seqAS += line.strip()
        return seq, seqAS

    def mapResidues(self):
        """Map the indexes of the input sequence to the residue numbers of the chain of interest,
        as {seqIdx: structResNumber}"""
        seq, seqAS = self.parseSequences()
        resIdxs = [res.get_id()[1] for res in self.getChainResidues()]

        mapDic, i, j = {}, 0, 0
        for k in range(len(seq)):
            if seq[k] == '-':
                i += 1
            if seqAS[k] == '-':
                j += 1
            if seq[k] != '-' and seqAS[k] != '-' and k - j < len(resIdxs):
                mapDic[k - i + 1] = resIdxs[k - j]
        return mapDic

    def getSeqROIResidues(self, seqROI, mapDic):
        """Residues of a sequence ROI, as [chainId_resNumber] in the structure numbering"""
        chainId, resIds = self.getChainName(), []
        for idx in range(seqROI.getROIIdx(), seqROI.getROIIdx2() + 1):
            if idx in mapDic:
                resIds.append(f'{chainId}_{mapDic[idx]}')
        return resIds

    def getSeqROIName(self, seqROI):
        """Name of the sequence ROI: its description (e.g. the mutation it originates from), its name or its id"""
        for name in [seqROI._ROISequence.getDescription(), seqROI.getROIName(), seqROI.getROIId()]:
            if name and str(name).strip() and str(name) != 'None':
                return str(name).strip()
        return f'ROI_{seqROI.getROIIdx()}-{seqROI.getROIIdx2()}'

    def getPartnerResidues(self, annotDic, seqROIsDic):
        """Residues of the partner chains closer than the maximum distance to the sequence ROIs residues
        found in each structural ROI"""
        maxDist, model = self.maxDistance.get(), self.getStructModel()
        chainIds = self.getPartnerChainIds()
        if not chainIds:
            chainIds = [chain.id for chain in model.get_chains() if chain.id != self.getChainName()]
        self.info(f'Looking for partner residues in chains {",".join(chainIds)}')

        partnerAtoms, partnerResIds = [], []
        for chainId in chainIds:
            for res in self.getChainResidues(chainId):
                for atom in res.get_atoms():
                    partnerAtoms.append(atom.get_coord().tolist())
                    partnerResIds.append(f'{chainId}_{res.get_id()[1]}')

        resCoords = {}
        for res in self.getChainResidues():
            resCoords[f'{self.getChainName()}_{res.get_id()[1]}'] = [a.get_coord().tolist() for a in res.get_atoms()]

        partnerDic = {}
        for roiId, matched in annotDic.items():
            roiResIds = {resId for sId in matched for resId in seqROIsDic[sId]['residues']}
            coords = [coord for resId in roiResIds for coord in resCoords.get(resId, [])]
            contacts = set()
            if coords and partnerAtoms:
                dists = spatial.distance.cdist(partnerAtoms, coords)
                for i in range(len(partnerAtoms)):
                    if min(dists[i, :]) < maxDist:
                        contacts.add(partnerResIds[i])
            partnerDic[roiId] = natural_sort(list(contacts))
        return partnerDic

    def writeReport(self, seqROIsDic, annotDic, partnerDic):
        with open(self.getReportFile(), 'w') as f:
            nFound = len([roiId for roiId in annotDic if annotDic[roiId]])
            f.write(f'{len(seqROIsDic)} of the {len(self.inputSequenceROIs.get())} sequence ROIs mapped to the '
                    f'structure, found in {nFound} of the {len(self.inputStructROIs.get())} structural ROIs\n')
            for roiId, matched in annotDic.items():
                if matched:
                    names = ','.join([seqROIsDic[sId]['name'] for sId in matched])
                    line = f'ROI {roiId}: {names}'
                    if self.doPartner.get():
                        line += f' | partner residues: {",".join(partnerDic.get(roiId, [])) or "None"}'
                    f.write(line + '\n')

    def getReportFile(self):
        return self._getExtraPath('annotation.txt')

    def getAnnotationFile(self):
        return self._getExtraPath('annotation.json')

    def saveAnnotation(self, annot):
        with open(self.getAnnotationFile(), 'w') as f:
            json.dump(annot, f)

    def readAnnotation(self):
        with open(self.getAnnotationFile()) as f:
            return json.load(f)

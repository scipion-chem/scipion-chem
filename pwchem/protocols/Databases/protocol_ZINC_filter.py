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

import os, glob
from urllib.request import urlopen

from pwem.protocols import EMProtocol
from pyworkflow.protocol.params import PointerParam, EnumParam, TextParam, LabelParam, STEPS_PARALLEL

from pwchem.utils import performBatchThreading

REMOVE, KEEP = 0, 1

avaiChoices = ['not-for-sale', 'agent', 'in-stock', 'boutique', 'now', 'bb', 'on-demand', 'wait-ok', 'for-sale']
bioActChoices = ['in-cells-only', 'fda', 'world-not-fda', 'investigational-only', 'world', 'in-trials', 'in-vivo-only',
                 'in-man-only', 'in-man', 'in-vivo', 'in-cells', 'in-vitro-only', 'in-vitro']
bioGenChoices = ['nonhuman-metabolites', 'endogenous', 'metabolites', 'natural-products', 'biogenic']
otherChoices = ['aggregators', 'named']
reactChoices = ['not-anodyne', 'hot-ok', 'hot', 'reactive', 'standard', 'anodyne', 'reactive-ok', 'clean',
                'standard-ok']


class ProtChemZINCFilter(EMProtocol):
    """Filter a set of small molecules by being in all selected catalogs of ZINC.
       See https://zinc15.docking.org/substances/subsets/"""
    _label = 'ZINC filter'
    subGroups = {'Availability': avaiChoices, 'Bioactive and Drugs': bioActChoices, 'Biogenic': bioGenChoices,
                 'Other': otherChoices, 'Reactivity': reactChoices}

    def __init__(self, *args, **kwargs):
        EMProtocol.__init__(self, *args, **kwargs)
        self.stepsExecutionMode = STEPS_PARALLEL

    def _defineParams(self, form):
        form.addSection(label='Input')
        group = form.addGroup('Input')
        group.addParam('inputSet', PointerParam, pointerClass="SetOfSmallMolecules",
                        label='Set to filter:', allowsNull=False)

        group = form.addGroup('Define filter')
        group.addParam('mode', EnumParam, choices=["Remove", "Keep"], label='Mode: ', default=0,
                      help='Whether to remove or keep the entry if it is included in the selected subset')

        group.addParam('subGroup', EnumParam, choices=list(self.subGroups.keys()), label='ZINC subgroups: ', default=0,
                       help='Subset group to filter by.\n'
                            'For more info check https://zinc15.docking.org/substances/subsets/')

        for i, subKey in enumerate(self.subGroups):
            group.addParam('subset_{}'.format(subKey), EnumParam, choices=self.subGroups[subKey],
                           label='{} subset: '.format(subKey), default=0, condition='subGroup=={}'.format(i),
                           help='Filter by presence in the selected subset.\n'
                                'For more info check https://zinc15.docking.org/substances/subsets/')

        group.addParam('addFilter', LabelParam, label='Add filter expression: ',
                       help='Add filter expression to the list')

        group = form.addGroup('Filter summary')
        group.addParam('filterList', TextParam, width=70, default='', label='List of filtering expressions: ',
                       help='List of filtering expressions the molecules have to pass\n')

        form.addParallelSection(threads=4, mpi=1)

    # --------------------------- INSERT steps functions --------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('operateStep')

    def operateStep(self):
        filDic, subsets = self.parseFilter()
        nt = self.numberOfThreads.get()

        presBase, outDir = 'subsets', os.path.abspath(self._getTmpPath())
        passMols = performBatchThreading(self.checkPresence, self.inputSet.get(), nt,
                                         presBase=os.path.join(outDir, presBase))

        allPresFile = self.getSummaryFile()
        with open(allPresFile, 'w') as f:
            f.write('{}\tKEEP:\t{}\tREMOVE:\t{}\tADD\n'.
                    format('ZINC_ID', '\t'.join(filDic['Keep']), '\t'.join(filDic['Remove'])))
            for presFile in glob.glob(os.path.join(outDir, '{}_*'.format(presBase))):
                presFile = os.path.join(outDir, presFile)
                with open(presFile) as fcur:
                    f.write(fcur.read())
        
        outputSet = self.inputSet.get().create(self._getPath())
        for mol in passMols:
            outputSet.append(mol)

        if len(outputSet)>0:
            self._defineOutputs(outputSmallMolecules=outputSet)
            self._defineSourceRelation(self.inputSet, outputSet)

    def _validate(self):
        errors = []
        firstMol = self.inputSet.get().getFirstItem()
        molName = firstMol.getMolName()
        if not 'ZINC' in molName and not hasattr(firstMol, 'ZINC_ID'):
            errors.append("Cannot find the ZINC ID in the molName or a ZINC_ID attribute.\n"
                          "We suggest to use the small molecules identification protocol")

        filDic, _ = self.parseFilter()
        interSubsets = set(filDic['Keep']).intersection(set(filDic['Remove']))
        if len(interSubsets) > 0:
            errors.append("The filter cannot ask to keep and remove entries from {} subsets at the same time".
                          format(' '.join(interSubsets)))

        return errors

    def _summary(self):
        filDic, _ = self.parseFilter()
        return ['Keep molecules present in subsets: {}\nRemove molecules present in subsets: {}\n\n'
                'Check presence for each molecule in file {}'.
                    format(', '.join(filDic['Keep']), ', '.join(filDic['Remove']), self.getSummaryFile())]


    def checkPresence(self, mols, molLists, it, presBase):
        filDic, subsets = self.parseFilter()
        presFile = '{}_{}.txt'.format(presBase, it)
        with open(presFile, 'w') as f:
            pass

        for mol in mols:
            zincId = self.getZINC_ID(mol)
            if not zincId:
                continue
            url = "http://zinc15.docking.org/substances/%s" % zincId

            presentSubsets = []
            inSubsets, add = False, True
            try:
                with urlopen(url) as response:
                    mystr = response.read().decode("utf8")

                for line in mystr.split('\n'):
                    if not inSubsets:
                        if '<strong class="pull-left">In:</strong>' in line:
                            inSubsets = True

                    else:
                        if '</ul>' in line:
                            break
                        for subset in subsets:
                            if "/substances/subsets/{}/".format(subset) in line:
                                presentSubsets.append(subset)

                keepSubsets = set(filDic['Keep']).intersection(set(presentSubsets))
                removeSubsets = set(filDic['Remove']).intersection(set(presentSubsets))
                if len(keepSubsets) < len(filDic['Keep']) or len(removeSubsets) > 0:
                    add = False

                keepBools, remBools = [str(sb in presentSubsets) for sb in filDic['Keep']], \
                                      [str(sb in presentSubsets) for sb in filDic['Remove']]
                with open(presFile, 'a') as f:
                    f.write('{}\t-\t{}\t-\t{}\t{}\n'.format(zincId, '\t'.join(keepBools), '\t'.join(remBools), add))
            except:
                print("Could not retrieve {} information from ZINC. Keeping the molecule by default".format(zincId))
                with open(presFile, 'a') as f:
                    f.write('{}\t-\t{}\t-\t{}\t{}\n'.format(zincId, '\t'.join(['Not_found'] * len(filDic['Keep'])),
                                                            '\t'.join(['Not_found'] * len(filDic['Remove'])), add))

            if add:
                molLists[it].append(mol)

    def getSummaryFile(self):
        return self._getPath("summary.txt")

    def createElementLine(self):
        keep = self.getEnumText('mode')
        subsetKey = self.getEnumText('subGroup')
        subset = self.getEnumText(f'subset_{subsetKey}')
        return f'{keep} if in {subset}\n'

    def parseFilter(self):
        filDic = {'Keep': [], 'Remove': []}
        filterStr = self.filterList.get().strip()
        for fil in filterStr.split('\n'):
            filDic[fil.split()[0]].append(fil.split()[-1])

        filDic = {'Keep': list(set(filDic['Keep'])), 'Remove': list(set(filDic['Remove']))}
        subsets = set(filDic['Remove']).union(set(filDic['Keep']))
        return filDic, subsets

    def getZINC_ID(self, mol):
        zincID = None
        if hasattr(mol, 'ZINC_ID') and getattr(mol, 'ZINC_ID').get():
            zincID = getattr(mol, 'ZINC_ID').get()
        elif 'ZINC' in mol.getMolName():
            zincID = mol.getMolName()
        return zincID
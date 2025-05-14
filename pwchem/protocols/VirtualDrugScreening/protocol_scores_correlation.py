# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors: Daniel Del Hoyo (ddelhoyo@cnb.csic.es)
# *
# * Biocomputing Unit, CNB-CSIC
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

"""
"""
import os, math, scipy
import numpy as np
import matplotlib.pyplot as plt

from pyworkflow.protocol import params

from pwem.protocols import EMProtocol

DIFF, RATIO = 1, 2

class ProtScoreCorrelation(EMProtocol):
    """
    Protocol to analyze the correlation between the scores of two sets of objects
    """

    _label = 'Score correlation'
    _inLabels = ['1', '2']

    # -------------------------- DEFINE param functions ----------------------
    def _defineInput(self, form, label, allowRatio=False, ratioCondition='True'):
      form.addParam(f'inputName_{label}', params.StringParam, label="Select name for score: ", default='',
                     condition=f'{ratioCondition} and {allowRatio}', expertLevel=params.LEVEL_ADVANCED,
                     help='Select the name of the score to be plotted in the output graph')
      group = form.addGroup(f'Input Set {label}')
      group.addParam(f'inputFromFile_{label}', params.BooleanParam, default=False,
                     label='Get input scores from file: ', expertLevel=params.LEVEL_ADVANCED,
                     condition=f'{ratioCondition}',
                     help='Whether to get the input scores from a file. The file must be a comma-separated file '
                          'as "ID,Score" with the selected ID and score values.')
      group.addParam(f'filePath_{label}', params.PathParam, label='Scores file: ',
                     condition=f'{ratioCondition} and inputFromFile_{label}',
                     help='CSV file containing a set of scores to use in form "ID,Score"')

      group.addParam(f'inputSet_{label}', params.PointerParam, pointerClass="EMSet",
                     label='Input set: ', allowsNull=allowRatio,
                     condition=f'{ratioCondition} and not inputFromFile_{label}',
                     help='Input set containing the score whose score will be analyzed')
      group.addParam(f'inputID_{label}', params.StringParam, label="Select ID: ", default='',
                     condition=f'{ratioCondition} and not inputFromFile_{label}',
                     help='Select the field to be used as ID for the selected input')
      group.addParam(f'inputScore_{label}', params.StringParam, label="Select score: ", default='',
                     condition=f'{ratioCondition} and not inputFromFile_{label}',
                     help='Select a score for the selected input to be used')

      group.addParam(f'extraAction_{label}', params.EnumParam, choices=['None', 'Log'],
                     default=0, condition=f'{ratioCondition}',
                     label='Extra action to take on scores: ', expertLevel=params.LEVEL_ADVANCED,
                     help='Extra function to pass on the score')
      form.addParam(f'compareSet_{label}', params.EnumParam, choices=['None', 'Difference', 'Ratio'],
                    default=0, condition=f'{ratioCondition} and {allowRatio}',
                    label='Compare the score to other set using: ', expertLevel=params.LEVEL_ADVANCED,
                    help='Choose the force filed to perform the ligand optimization')

      if allowRatio:
        self._defineInput(form, f'{label}_ratio', allowRatio=False,
                          ratioCondition=f'compareSet_{label} in [{DIFF},{RATIO}]')

    def _defineParams(self, form):
        """ Define the input parameters that will be used.
        """
        for label in self._inLabels:
          form.addSection(label=f'Input {label}')
          self._defineInput(form, label=label, allowRatio=True)

        form.addSection(label='Correlation')
        group = form.addGroup('Correlation')
        group.addParam('corrType', params.EnumParam, choices=['Pearson', 'Spearman'], default=0,
                       label='Correlation analysis: ',
                       help='Type of correlation analysis to perform')

    # --------------------------- STEPS functions ------------------------------

    def _insertAllSteps(self):
        self._insertFunctionStep(self.correlationStep)

    def saveSummary(self):
      with open(self.getOutputTxtPath(), 'w') as f:
        f.write(f'{self.getEnumText("corrType")} correlation: {self.corResult[0]}\n'
                f'With p-value: {self.corResult[1]}')

    def correlationStep(self):
      scoreDics = {}
      for label in self._inLabels:
        scoreDics[label] = self.getScoreDic(label)
        if getattr(self, f'compareSet_{label}').get() in [DIFF, RATIO]:
          setRatioDic = self.getScoreDic(f'{label}_ratio')
          scoreDics[label] = self.operateSetDics(scoreDics[label], setRatioDic,
                                                 op=self.getEnumText(f'compareSet_{label}'))

      self.corResult = self.calculateCorrelation(scoreDics['1'], scoreDics['2'],
                                                   corrAnalysis=self.getEnumText('corrType'))
      self.saveSummary()

      self.plotScoreDistribution(scoreDics['1'], scoreDics['2'], self.corResult)

    # --------------------------- UTILS functions ------------------------------
    def _validate(self):
        """ Validate if the inputs are in mol2 or pdb format
        """
        errors = []
        return errors

    def _summary(self):
        summ = []
        if os.path.exists(self.getOutputTxtPath()):
          with open(self.getOutputTxtPath()) as f:
            summ.append(f.read())
        return summ

    def parseScoreFile(self, file):
      scoreDic = {}
      with open(file) as f:
        for line in f:
          if line.strip():
            sline = line.split(',')
            scoreDic[sline[0].strip()] = float(sline[1].strip())
      return scoreDic

    def parseScoreSet(self, inputSet, idAttr, scoreAttr):
      scoreDic = {}
      for item in inputSet:
        scoreDic[getattr(item, idAttr).get()] = float(getattr(item, scoreAttr).get())
      return scoreDic

    def getScoreDic(self, label):
      if getattr(self, f'inputFromFile_{label}'):
        scFile = getattr(self, f'filePath_{label}').get()
        scoreDic = self.parseScoreFile(scFile)
      else:
        inSet = getattr(self, f'inputSet_{label}').get()
        idAttr, scoreAttr = getattr(self, f'inputID_{label}').get(), getattr(self, f'inputScore_{label}').get()
        scoreDic = self.parseScoreSet(inSet, idAttr, scoreAttr)

      if self.getEnumText(f'extraAction_{label}') == 'Log':
        for name, sc in scoreDic.items():
          scoreDic[name] = math.log(sc)
      return scoreDic

    def operateSetDics(self, d1, d2, op='ratio'):
      nDic = {}
      for name in d1:
        sc1, sc2 = d1[name], d2[name]
        if op == 'Ratio':
          nDic[name] = sc1 / sc2
        elif op == 'Difference':
          nDic[name] = sc1 - sc2
      return nDic

    def calculateCorrelation(self, d1, d2, corrAnalysis='Pearson'):
      vs1, vs2 = [d1[name] for name in d1 if name in d2], [d2[name] for name in d1 if name in d2]
      if corrAnalysis == 'Pearson':
        return scipy.stats.pearsonr(vs1, vs2)
      elif corrAnalysis == 'Spearman':
        return scipy.stats.spearmanr(vs1, vs2)

    def plotScoreDistribution(self, d1, d2, corr):
        molNames = [name for name in d1 if name in d2]
        x, y = [d1[name] for name in molNames], [d2[name] for name in molNames]
        coef = np.polyfit(x, y, 1)
        poly1dFn = np.poly1d(coef)

        plt.title(f'Scores correlation ({self.getEnumText("corrType")}): {round(corr[0], 4)}')
        xlabel = self.inputName_1.get().strip() if self.inputName_1.get().strip() else "Input 1"
        ylabel = self.inputName_2.get().strip() if self.inputName_2.get().strip() else "Input 2"
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.plot(x, y, 'yo', x, poly1dFn(x), 'r')
        plt.grid()

        if len(molNames) < 20:
          for i, txt in enumerate(molNames):
            plt.annotate(txt, (x[i], y[i]), fontsize=6)

        plt.savefig(self.getOuputImgPath())

    def getOuputImgPath(self):
      return self._getExtraPath('scoreDistribution.png')

    def getOutputTxtPath(self):
      return self._getPath('correlation.txt')

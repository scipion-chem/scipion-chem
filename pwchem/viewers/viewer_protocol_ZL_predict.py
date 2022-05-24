# **************************************************************************
# *
# * Authors:  Carlos Oscar Sorzano (coss@cnb.csic.es)
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
from pathlib import Path

from pwchem.protocols import ProtChemZLPredict
from pwchem import Plugin
from pyworkflow.viewer import DESKTOP_TKINTER, ProtocolViewer
from pyworkflow.protocol.params import LabelParam, StringParam
import pyworkflow.utils as pwutils
from pwem.objects.data import AtomStruct

CQUARKSERVER="zhanglab.ccmb.med.umich.edu"

class ProtChemZLPredictViewer(ProtocolViewer):
    """ Visualize the output of protocol Zhang Lab """
    _label = 'viewer Zhang lab'
    _targets = [ProtChemZLPredict]
    _environments = [DESKTOP_TKINTER]

    def __init__(self, **args):
        ProtocolViewer.__init__(self, **args)

    def _defineParams(self, form):
        form.addSection(label='Visualization')
        # Select the level to show
        form.addParam('url', StringParam,
                      label="URL of Zhang Lab results",
                      help='It is only needed  only the first time.\n'
                           'Example http://zhanglab.ccmb.med.umich.edu/C-QUARK/output/QB663.')
        form.addParam('chimera', LabelParam,
                      label="See all models in chimera",
                      help='It is available when the results have first been downloaded')

    def _getVisualizeDict(self):
        return {'url': self._viewResults,
                'chimera': self._viewChimera}

    def getResultsDir(self):
        fnDir = self.protocol._getExtraPath(CQUARKSERVER)
        fnBaseDir=None
        for fn in Path(fnDir).rglob('*.html'):
            if fn.name=="index.html":
                fnBaseDir=str(fn.parent)
                break
        return fnBaseDir

    def _viewResults(self, e=None):
        import os
        views = []
        fnBaseDir = self.getResultsDir()
        if not fnBaseDir:
            if os.path.exists(self.protocol._getExtraPath(CQUARKSERVER)):
                os.system("rm -rf %s"%self.protocol._getExtraPath(CQUARKSERVER))
            url=self.url.get().strip()
            if not url.endswith("index.html"):
                url+="/index.html"
            urlDir=url.replace("/index.html","").split(".edu")[1]
            if 'I-TASSER' in urlDir:
                includeDirs = "-I jsmol,3Dmol,I-TASSER/output"
            else:
                includeDirs = "-I jsmol,3Dmol,%s" % urlDir
            os.system('cd %s; wget %s --mirror -p --convert-links -P . -r --retry-connrefused --waitretry=1 --read-timeout=20 --timeout=15 -t 5 %s'%\
                      (self.protocol._getExtraPath(),includeDirs,url))
            if 'I-TASSER' in urlDir:
                extraDownloads = ['I-TASSER/output/bin/jmol/j2s/core/package.js',
                                  'I-TASSER/output/bin/jmol/j2s/core/corescript.z.js',
                                  'I-TASSER/output/bin/jmol/j2s/core/core.z.js',
                                  'I-TASSER/output/bin/jmol/j2s/core/corebio.z.js',
                                  'I-TASSER/output/bin/jmol/j2s/JM/Resolver.js',
                                  'I-TASSER/output/bin/jmol/j2s/J/shape/Mesh.js',
                                  'I-TASSER/output/bin/jmol/j2s/J/render/MeshRenderer.js',
                                  'I-TASSER/output/bin/jmol/j2s/core/corescriptcmd.z.js',
                                  'I-TASSER/output/bin/jmol/j2s/J/thread/SpinThread.js',
                                  'I-TASSER/output/bin/jmol/j2s/J/g3d/HermiteRenderer.js',
                                  'I-TASSER/output/bin/jmol/j2s/core/coretext.z.js']
                for fn in extraDownloads:
                    os.system('cd %s; wget --mirror -p --convert-links -P . -r --retry-connrefused --waitretry=1 --read-timeout=20 --timeout=15 -t 5 https://%s/%s'%\
                              (self.protocol._getExtraPath(),CQUARKSERVER,fn))
            fnBaseDir = self.getResultsDir()
        if fnBaseDir:
            url=os.path.abspath(os.path.join(fnBaseDir,"index.html"))
            pwutils.runJob(None,"python",Plugin.getPluginHome('utils/showZLPredictResults.py')+" "+url)
            #webbrowser.open_new_tab(url)
            if not hasattr(self.protocol,"outputPdb_1"):
                for fn in Path(fnBaseDir).rglob('model*.pdb'):
                    self.constructOutput(str(fn))
        return views

    def _viewChimera(self, e=None):
        from chimera import Plugin as chimera_plugin
        args=""
        fnBaseDir = self.getResultsDir()
        if fnBaseDir:
            for fn in Path(fnBaseDir).rglob('model*.pdb'):
                args+=str(fn)+" "
        os.system("%s %s &"%(chimera_plugin.getProgram(),args))

    def constructOutput(self,fnPdb):
        fnDir, fnResults=os.path.split(fnPdb)

        pdb=AtomStruct(filename=fnPdb)
        pdbNo=fnResults[5]
        outputName = 'outputPdb_%s' % pdbNo
        if not hasattr(self.protocol,outputName):
            outputDict = {outputName: pdb}
            self.protocol._defineOutputs(**outputDict)
            self.protocol._defineSourceRelation(self.protocol.inputSeq, pdb)

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

import glob
from pathlib import Path

from bioinformatics.protocols.protocol_raptorX import ProtBioinformaticsRaptorX
from bioinformatics import Plugin
from pyworkflow.viewer import DESKTOP_TKINTER, ProtocolViewer
from pyworkflow.protocol.params import LabelParam, StringParam
import pyworkflow.utils as pwutils
from pwem.objects.data import AtomStruct

SERVERDIR = 'raptorx.uchicago.edu'

class ProtBioinformaticsRaptorXViewer(ProtocolViewer):
    """ Visualize the output of protocol Raptor X"""
    _label = 'viewer raptor X'
    _targets = [ProtBioinformaticsRaptorX]
    _environments = [DESKTOP_TKINTER]

    def __init__(self, **args):
        ProtocolViewer.__init__(self, **args)

    def _defineParams(self, form):
        form.addSection(label='Visualization')
        # Select the level to show
        form.addParam('url', StringParam,
                      label="URL of Raptor X results",
                      help='It is only needed  only the first time.\n'
                           'Example http://raptorx.uchicago.edu/StructPredV2/myjobs/19070092_548639.')

    def _getVisualizeDict(self):
        return {'url': self._viewResults}

    def getResultsDir(self):
        fnDir = self.protocol._getExtraPath(SERVERDIR)
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
            if os.path.exists(self.protocol._getExtraPath(SERVERDIR)):
                os.system("rm -rf %s"%self.protocol._getExtraPath(SERVERDIR))
            url=self.url.get().strip()
            if url.endswith('/'):
                url=url[:-1]
            includeDirs = '-I site_media/jsmol -I site_media/jsmol/j2s/core -I StructPredV2/myjobs'
            os.system('cd %s; wget %s --mirror -p --convert-links -P . -r --retry-connrefused --waitretry=1 --read-timeout=20 --timeout=15 -t 5 %s'%\
                      (self.protocol._getExtraPath(),includeDirs,url))
            fnResult = os.path.join(self.protocol._getExtraPath(SERVERDIR,'StructPredV2','myjobs',os.path.split(url)[1]))
            pwutils.createLink(fnResult,
                               os.path.join(self.protocol._getExtraPath(SERVERDIR,'StructPredV2','myjobs','index.html')))
            extraDownloads = ['site_media/jsmol/j2s/core/package.js',
                              'site_media/jsmol/j2s/core/core.z.js',
                              'site_media/jsmol/j2s/core/corescript.z.js',
                              'site_media/jsmol/j2s/core/corescript2.z.js',
                              'site_media/jsmol/j2s/core/corestate.z.js',
                              'site_media/jsmol/j2s/core/coretext.z.js',
                              'site_media/jsmol/j2s/core/corezip.z.js',
                              'site_media/jsmol/j2s/core/coremenu.z.js',
                              'site_media/jsmol/j2s/core/corebio.z.js',
                              'site_media/jsmol/j2s/J/shape/Measures.js',
                              'site_media/jsmol/j2s/J/shape/Mesh.js',
                              'site_media/jsmol/j2s/J/render/MeasuresRenderer.js',
                              'site_media/jsmol/j2s/J/render/MeshRenderer.js',
                              'site_media/jsmol/j2s/J/util/Hermite.js',
                              'site_media/jsmol/j2s/J/thread/SpinThread.js',
                              'site_media/jsmol/j2s/J/g3d/HermiteRenderer.js'
                              ]
            for fn in extraDownloads:
                os.system('cd %s; wget --mirror -p --convert-links -P . -r --retry-connrefused --waitretry=1 --read-timeout=20 --timeout=15 -t 5 http://%s/%s'%\
                          (self.protocol._getExtraPath(),SERVERDIR,fn))
            fh=open(fnResult, encoding="utf-8")
            for line in fh.readlines():
                if 'function loadSummaryPicture' in line:
                    tokens = line.split('function loadSummaryPicture')
                    otherId = tokens[1].split('(')[0]
                    fnSummary = os.path.join(self.protocol._getExtraPath(SERVERDIR,'StructPredV2','myjobs','%ssummary_data'%otherId))
                    os.system('cd %s; wget --mirror -p --convert-links -P . -r --retry-connrefused --waitretry=1 --read-timeout=20 --timeout=15 -t 5 http://%s/%s'%\
                              (self.protocol._getExtraPath(),SERVERDIR,'StructPredV2/myjobs/%ssummary_data'%otherId))
                    urlId = url.split("_")[1]
                    pwutils.createLink(fnSummary,os.path.join(self.protocol._getExtraPath(SERVERDIR,'StructPredV2','myjobs','%s.all_in_one.zip'%urlId)))
                    fnPDB = os.path.join(self.protocol._getExtraPath(SERVERDIR,'StructPredV2','myjobs','%ssummary_pdb_image'%otherId))
                    os.system('cd %s; wget --mirror -p --convert-links -P . -r --retry-connrefused --waitretry=1 --read-timeout=20 --timeout=15 -t 5 http://%s/%s'%\
                              (self.protocol._getExtraPath(),SERVERDIR,'StructPredV2/myjobs/%ssummary_pdb_image'%otherId))
                    pwutils.createLink(fnPDB,os.path.join(self.protocol._getExtraPath(SERVERDIR,'StructPredV2','myjobs','%s.pdb'%urlId)))
                    pwutils.createLink(fnPDB,self.protocol._getExtraPath('%s.pdb'%urlId))
                if ".pdb'" in line:
                    fn = line.split(' = ')[1].replace(';','').strip()
                    os.system('cd %s; wget --mirror -p --convert-links -P . -r --retry-connrefused --waitretry=1 --read-timeout=20 --timeout=15 -t 5 http://%s/%s.gz' % \
                        (self.protocol._getExtraPath(), SERVERDIR, fn))
                    fn = fn.replace("'","")[1:]
                    fnDomain = self.protocol._getExtraPath(os.path.split(fn)[1]+".gz")
                    pwutils.copyFile(self.protocol._getExtraPath(SERVERDIR+"/"+fn+".gz"),fnDomain)
                    os.system('cd %s; gunzip %s'%(self.protocol._getExtraPath(),os.path.split(fnDomain)[1]))
            fh.close()
            os.system('rm %s/wget-log*'%self.protocol._getExtraPath())

            fnBaseDir = self.getResultsDir()
        if fnBaseDir:
            url=os.path.abspath(os.path.join(fnBaseDir,"index.html"))
            pwutils.runJob(None,"python",Plugin.getPluginHome('utils/showRaptorXResults.py')+" "+url)
            if not hasattr(self.protocol,"outputPdb"):
                for fn in glob.glob(self.protocol._getExtraPath('*.pdb')):
                    suffix = ''
                    fnBase = os.path.split(fn)[1]
                    if '_' in fnBase:
                        suffix='_Domain%s'%fnBase.split('_')[1].replace('.pdb','')
                    self.constructOutput(fn, suffix)
        return views

    def constructOutput(self,fnPdb, suffix=""):
        pdb=AtomStruct(filename=fnPdb)
        outputName = 'outputPdb%s' % suffix
        if not hasattr(self.protocol,outputName):
            outputDict = {outputName: pdb}
            self.protocol._defineOutputs(**outputDict)
            self.protocol._defineSourceRelation(self.protocol.inputSeq, pdb)

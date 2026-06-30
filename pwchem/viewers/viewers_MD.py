# **************************************************************************
# *
# * Authors:     Carlos Oscar Sorzano (coss@cnb.csic.es)
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
import subprocess

import pyworkflow.viewer as pwviewer
from pyworkflow.protocol import params
from pwem.viewers import Chimera, ChimeraView

# pwchem imports
from .. import Plugin
from ..constants import MDTRAJ_DIC, OPENBABEL_DIC, TCL_MD_STR, PML_MD_STR, PML_MD_STR_AMBER, TCL_MD_LIG_STR
from ..viewers import PyMolViewer, PyMolView, VmdViewPopen
from ..objects import MDSystem

_ANAL_RMSD      = 0
_ANAL_RMSF      = 1
_ANAL_RG        = 2
_ANAL_SS        = 3
_ANAL_SASA      = 4
_ANAL_PCA       = 5
_ANAL_DISTANCE  = 6

_ANAL_CHOICES = ['RMSD', 'RMSF', 'Rg', 'Secondary Structure', 'SASA', 'PCA', 'Atom distance']

# Analyses that use the atom-selection dropdown
_USES_SEL_ATOMS   = [_ANAL_RMSD, _ANAL_RMSF]

class MDSystemViewer(pwviewer.Viewer):
  _label = 'Viewer Molecular Dynamics system'
  _environments = [pwviewer.DESKTOP_TKINTER]
  _targets = []

  def _visualize(self, obj, onlySystem=False, trjFile=None, **kwargs):
    systemFile = os.path.abspath(obj.getSystemFile())
    if not trjFile:
        trjFile = obj.hasTrajectory()

    if not trjFile or onlySystem:
        pymolV = PyMolViewer(project=self.getProject())
        return pymolV._visualize(systemFile, cwd=os.path.dirname(systemFile))

    else:
        trjFile = os.path.abspath(obj.getTrajectoryFile())
        topoFile = os.path.abspath(obj.getTopologyFile())
        if topoFile.lower().endswith(('.top', '.tpr')):
            topoFile = obj.getSystemFile() # needed for gromacs topology
        outPml = os.path.join(os.path.dirname(trjFile), 'pymolSimulation.pml')
        _, trjExt = os.path.splitext(trjFile)
        if trjExt in ['.nc', '.netcdf']:
            template = PML_MD_STR_AMBER
        else:
            template = PML_MD_STR
        with open(outPml, 'w') as f:
          f.write(template.format(os.path.abspath(topoFile),
                                    os.path.abspath(trjFile)))

        return [PyMolView(os.path.abspath(outPml), cwd=os.path.dirname(trjFile))]


class MDSystemPViewer(pwviewer.ProtocolViewer):
    """Visualize the output of a Molecular Dynamics simulation."""
    _label = 'Viewer Molecular Dynamics System'
    _targets = [MDSystem]
    _mdtrajScript   = 'mdtraj_analysis.py'
    _prolifViewScript = 'prolif_viewer.py'

    def __init__(self, **args):
        super().__init__(**args)

    # ------------------------------------------------------------------
    # Form definition
    # ------------------------------------------------------------------

    def _defineSimParams(self, form):
        group = form.addGroup('Open MD simulation')
        group.addParam('displayMdPymol', params.LabelParam,
                       label='Display trajectory with PyMol: ',
                       help='Display trajectory with PyMol.\n'
                            'Protein as NewCartoon, waters as sticks.')
        group.addParam('displayMdVMD', params.LabelParam,
                       label='Display trajectory with VMD: ',
                       help='Display trajectory with VMD.\n'
                            'Protein as NewCartoon, waters as dots.')

    def _defineVideoParams(self, form):
        form.addSection(label='Generate MD video')
        group = form.addGroup('Cinematic trajectory video',
                              help='Render a movie of the MD trajectory. The camera is auto-framed on '
                                   'the molecule and the frames are encoded into an mp4 (or gif). '
                                   'Inspired by the VisualFactory project (FindPerspective auto-camera '
                                   '+ UltimateSmoothMD smoothing).')

        group.addParam('vidEngine', params.EnumParam,
                       label='Rendering engine: ', default=0,
                       choices=['PyMol', 'ChimeraX'],
                       help='Engine used to render the movie.\n'
                            '* PyMol: head-less ray-tracing (no display needed). Robust default.\n'
                            '* ChimeraX: GPU-accelerated, VisualFactory-style cinematic quality. '
                            'Needs a graphical display (opens the ChimeraX window).')

        group.addParam('vidStyle', params.EnumParam,
                       label='Protein style: ', default=0,
                       choices=['Cartoon', 'Cartoon + sticks', 'Surface', 'Sticks', 'Ribbon'],
                       help='Representation used for the protein in the movie.')
        group.addParam('vidBg', params.EnumParam,
                       label='Background: ', default=0, choices=['White', 'Black'],
                       help='Background color of the rendered frames.')
        group.addParam('vidHighlightLig', params.BooleanParam,
                       label='Highlight ligand: ', default=True,
                       help='Show the ligand (resname "%s") as coloured sticks.'
                            % self.getMDSystem().getLigandID())

        group.addParam('vidResolution', params.EnumParam,
                       label='Resolution: ', default=1,
                       choices=['480p', '720p (HD)', '1080p (Full HD)', '2160p (4K)'],
                       help='Frame resolution. Higher resolutions render more slowly.')
        group.addParam('vidFps', params.IntParam,
                       label='Frames per second: ', default=15,
                       help='Playback speed of the resulting video.')
        group.addParam('vidStride', params.IntParam,
                       label='Frame stride: ', default=1,
                       help='Render every Nth trajectory frame (use >1 for long trajectories).')

        group.addParam('vidSmooth', params.IntParam,
                       label='Smoothing window: ', default=0, expertLevel=params.LEVEL_ADVANCED,
                       condition='vidEngine==0',
                       help='[PyMol only] Coordinate-smoothing window applied before rendering '
                            '(0 = none). Removes thermal jitter for cleaner playback '
                            '(visualization only).')
        group.addParam('vidRay', params.BooleanParam,
                       label='Cinematic ray-tracing: ', default=True, expertLevel=params.LEVEL_ADVANCED,
                       condition='vidEngine==0',
                       help='[PyMol only] Ray-trace each frame with soft shadows. '
                            'Disable for a faster, flatter preview.')
        group.addParam('vidSpin', params.BooleanParam,
                       label='Camera spin: ', default=False, expertLevel=params.LEVEL_ADVANCED,
                       help='Add a full 360 degree camera rotation across the trajectory.')
        group.addParam('vidFormat', params.EnumParam,
                       label='Output format: ', default=0, choices=['mp4', 'gif'],
                       help='Video container. mp4 needs ffmpeg; otherwise a gif is produced.')

        group.addParam('displayMDVideo', params.LabelParam,
                       label='Generate and open MD video: ',
                       help='Render the video and open it with the default system player. '
                            'The file is saved next to the trajectory as <systemName>_MDvideo.<ext>.')

    def _defineMDTrajParams(self, form):
        form.addSection(label='Trajectory analysis')
        group = form.addGroup('MDTraj analysis')

        group.addParam('mdAnalChoices', params.EnumParam,
                       label='Analysis type: ',
                       choices=_ANAL_CHOICES, default=_ANAL_RMSD,
                       help='Select the MDTraj analysis to run.\n'
                            'Relevant parameters will appear below.')

        # ── Shared: atom selection (RMSD / RMSF) ────────
        group.addParam('selAtoms', params.EnumParam,
                       label='Selection of atoms: ', default=0,
                       choices=['Protein', 'Backbone', 'CA', 'Sidechain', 'Ligand'],
                       condition=f'mdAnalChoices in {_USES_SEL_ATOMS}',
                       help='Atom selection used in the analysis.')

        # ── Shared: heavy atoms (RMSD / RMSF) ──────────────────
        group.addParam('heavyAtoms', params.BooleanParam,
                       label='Use only heavy atoms: ', default=True,
                       condition=f'mdAnalChoices in {_USES_SEL_ATOMS}',
                       help='Restrict analysis to non-hydrogen atoms.')

        # ── Secondary Structure sub-options ─────────────────────────
        group.addParam('ssDisplayType', params.EnumParam, default=0,
                       label='Display mode: ',
                       choices=['Per residue', 'Per frame', 'Heatmap'],
                       condition=f'mdAnalChoices == {_ANAL_SS}',
                       help='How to present the DSSP secondary-structure results.')

        # ── SASA sub-options ─────────────────────────────────────────
        group.addParam('sasaScope', params.EnumParam, default=0,
                       label='SASA scope: ',
                       choices=['All', 'Ligand'],
                       condition=f'mdAnalChoices == {_ANAL_SASA}',
                       help='Plot SASA for all atoms or the ligand only.\n'
                            'Computation time can be long.')

        # ── PCA sub-options ──────────────────────────────────────────
        group.addParam('pcaType', params.EnumParam, default=0,
                       label='PCA input: ',
                       choices=['Coordinates', 'Pairwise distances'],
                       condition=f'mdAnalChoices == {_ANAL_PCA}',
                       help='Run PCA on Cartesian coordinates or backbone pairwise distances.')

        # ── Distance sub-options ─────────────────────────────────────
        distLine = group.addLine('Atom indices: ',
                                 condition=f'mdAnalChoices == {_ANAL_DISTANCE}',
                                 help='Indices of the two atoms (from the PDB) whose distance '
                                      'will be tracked across frames.')
        distLine.addParam('atom1', params.IntParam, allowsNull=True, label='Atom 1: ')
        distLine.addParam('atom2', params.IntParam, allowsNull=True, label='Atom 2: ')

        # ── Run button ────────────────────────────────────────
        group.addParam('displayMDTrajAnalysis', params.LabelParam,
                       label='Display MDTraj analysis: ',
                       help='Run and display the selected analysis.')

        if self.getMDSystem().getProlifFile():
            self._defineProlifParams(form)

    def _defineProlifParams(self, form):
        form.addSection('Receptor-ligand interactions')
        group = form.addGroup('ProLIF analysis')
        group.addParam('displayFingerprint', params.LabelParam, label='Show interaction fingerprint: ',
                       help='Display interaction fingerprint generated by ProLIF')
        group.addParam('displayInterNetwork', params.LabelParam, label='Display interaction network: ',
                      help='Display ligand interaction network generated by ProLIF (web browser is opened).\n'
                           'Only the interactions that are present in more than 30% of the frames are shown.')
        group.addParam('displayProlifMatrix', params.LabelParam, label='Display similarity matrix: ',
                      help='Display Tanimoto similarity matrix calculated using ligand-target interaction fingerprint.')

    def _defineParams(self, form):
        form.addSection(label='Visualization of MD System')
        group = form.addGroup('Open MD system')
        group.addParam('displayPymol', params.LabelParam,
                       label='Open system in PyMol: ',
                       help='Display the system in the PyMol GUI.')

        if self.getMDSystem().hasTrajectory():
            self._defineSimParams(form)
            self._defineVideoParams(form)
            self._defineMDTrajParams(form)

    def getMDSystem(self, objType=MDSystem):
        if isinstance(self.protocol, objType):
            return self.protocol
        return self.protocol.outputSystem

    def _getVisualizeDict(self):
        return {
            'displayPymol':           self._showPymol,
            'displayMdPymol':         self._showMdPymol,
            'displayMdVMD':           self._showMdVMD,
            'displayMDTrajAnalysis':  self._showMDTrajAnalysis,
            'displayMDVideo':         self._showMDVideo,
            'displayFingerprint':     self._showProlifFp,
            'displayInterNetwork':    self._showProlifNetwork,
            'displayProlifMatrix':    self._showProlifMatrix,
        }

    def _showMDTrajAnalysis(self, paramName=None):
        dispatch = {
            _ANAL_RMSD:     self._showMDTrajRMSDRMSF,
            _ANAL_RMSF:     self._showMDTrajRMSDRMSF,
            _ANAL_RG:       self._showMDTrajRGAnalysis,
            _ANAL_SS:       self._showMDTrajSSAnalysis,
            _ANAL_SASA:     self._showMDTrajSASAAnalysis,
            _ANAL_PCA:      self._showMDTrajPCAAnalysis,
            _ANAL_DISTANCE: self._showDistance,
        }
        return dispatch[self.mdAnalChoices.get()]()

    # ------------------------------------------------------------------
    # Visualize methods
    # ------------------------------------------------------------------

    def _showPymol(self, paramName=None):
        system = self.getMDSystem()
        return MDSystemViewer(project=self.getProject())._visualize(system, onlySystem=True)

    def _showMdPymol(self, paramName=None):
        system = self.getMDSystem()
        return MDSystemViewer(project=self.getProject())._visualize(system)

    def writeTCL(self, outTcl, sysFile, sysExt, sysTrj, trjExt):
        system  = self.getMDSystem()
        vmdStr  = TCL_MD_STR % (sysFile, sysExt, sysTrj, trjExt)
        vmdStr += TCL_MD_LIG_STR.format(system.getLigandID())
        with open(outTcl, 'w') as f:
            f.write(vmdStr)

    def _showMdVMD(self, paramName=None):
        system = self.getMDSystem()
        outTcl = os.path.join(os.path.dirname(system.getTrajectoryFile()), 'vmdSimulation.tcl')
        sysExt = os.path.splitext(system.getFileName())[1][1:]
        trjExt = os.path.splitext(system.getTrajectoryFile())[1][1:]
        self.writeTCL(outTcl, system.getFileName(), sysExt,
                      system.getTrajectoryFile(), trjExt)
        return [VmdViewPopen('-e {}'.format(outTcl))]

    # Common option tables shared by both engines.
    _VID_STYLES  = ['cartoon', 'cartoon+sticks', 'surface', 'sticks', 'ribbon']
    _VID_BGS     = ['white', 'black']
    _VID_RESOS   = [(854, 480), (1280, 720), (1920, 1080), (3840, 2160)]
    _VID_FORMATS = ['mp4', 'gif']

    def _getVideoPaths(self):
        system  = self.getMDSystem()
        # Both engines read the structure file (.gro/.pdb) as topology, never .top/.tpr.
        structFile = os.path.abspath(system.getSystemFile())
        trjFile    = os.path.abspath(system.getTrajectoryFile())
        workDir    = os.path.dirname(trjFile)
        outBase    = '{}_MDvideo'.format(system.getSystemName())
        return system, structFile, trjFile, workDir, outBase

    def _showMDVideo(self, paramName=None):
        if self.vidEngine.get() == 1:
            return self._showMDVideoChimeraX()
        return self._showMDVideoPymol()

    def _showMDVideoPymol(self):
        system, structFile, trjFile, workDir, outBase = self._getVideoPaths()
        resos = ['480p', '720p', '1080p', '4K']

        args  = (f'-- -i "{structFile}" -t "{trjFile}" -o "{outBase}" --workdir "{workDir}" '
                 f'--style {self._VID_STYLES[self.vidStyle.get()]} --bg {self._VID_BGS[self.vidBg.get()]} '
                 f'--ligand {system.getLigandID()} --highlightLig {int(self.vidHighlightLig.get())} '
                 f'--resolution {resos[self.vidResolution.get()]} --fps {self.vidFps.get()} '
                 f'--stride {max(1, self.vidStride.get())} --smooth {max(0, self.vidSmooth.get())} '
                 f'--ray {int(self.vidRay.get())} --spin {int(self.vidSpin.get())} '
                 f'--format {self._VID_FORMATS[self.vidFormat.get()]} --open 1')

        # create_MDvideo.py lives in pwchem/scripts and runs inside PyMol (bundled in the
        # OpenBabel env): "pymol -cq create_MDvideo.py -- args".
        Plugin.runScript(self, 'create_MDvideo.py', args, env=OPENBABEL_DIC,
                         pyStr='pymol -cq', popen=True, wait=False, cwd=workDir)

    def _showMDVideoChimeraX(self):
        """Render the movie with ChimeraX (VisualFactory-style cinematic quality).

        We write a ChimeraX command file (.cxc) and launch it with the standard pwem
        ChimeraView -- the same mechanism the other chem viewers use (see
        viewers_data._viewChimera). A .cxc (NOT a Python script) is mandatory for
        movie recording: ChimeraX runs .cxc commands through its frame-aware command
        queue, so `coordset` playback advances and `wait` works. The equivalent calls
        from a Python script block the event loop and freeze ChimeraX.

        This ChimeraX also rejects the 'end' keyword in the coordset range
        (`coordset #1 1,end` -> "Expected a keyword"); only a NUMERIC range works. So
        we count the trajectory frames first (mdtraj) and bake the number in.
        It needs a graphical display (opens the ChimeraX window).
        """
        system, structFile, trjFile, workDir, outBase = self._getVideoPaths()
        width, height = self._VID_RESOS[self.vidResolution.get()]
        bg     = self._VID_BGS[self.vidBg.get()]
        style  = self._VID_STYLES[self.vidStyle.get()]
        fmt    = self._VID_FORMATS[self.vidFormat.get()]
        outFile = os.path.join(workDir, '{}.{}'.format(outBase, 'mp4' if fmt == 'gif' else fmt))

        nFrames = self._countTrajFrames(structFile, trjFile, workDir)
        cxc = self._buildChimeraXCxc(structFile, trjFile, workDir, outFile, style, bg,
                                     width, height, self.vidFps.get(), nFrames,
                                     bool(self.vidSpin.get()), system.getLigandID(),
                                     bool(self.vidHighlightLig.get()))
        cxcFile = os.path.join(workDir, '{}.cxc'.format(outBase))
        with open(cxcFile, 'w') as f:
            f.write(cxc)
        return [ChimeraView(cxcFile)]

    def _countTrajFrames(self, structFile, trjFile, workDir):
        """Count trajectory frames with mdtraj (memory-light iterload), so the
        ChimeraX coordset range can be NUMERIC. Returns an int, or None on failure."""
        counter = os.path.join(workDir, '_count_frames.py')
        with open(counter, 'w') as f:
            f.write("import sys, mdtraj as md\n"
                    "n = 0\n"
                    "for ch in md.iterload(sys.argv[1], top=sys.argv[2], chunk=500):\n"
                    "    n += ch.n_frames\n"
                    "print('NFRAMES', n)\n")
        cmd = '{} && python "{}" "{}" "{}"'.format(
            Plugin.getEnvActivationCommand(MDTRAJ_DIC), counter, trjFile, structFile)
        try:
            out = subprocess.check_output(cmd, shell=True, text=True, stderr=subprocess.STDOUT)
            for line in out.splitlines():
                if line.startswith('NFRAMES'):
                    return int(line.split()[1])
        except Exception:
            pass
        return None

    def _buildChimeraXCxc(self, structFile, trjFile, workDir, outFile, style, bg,
                          width, height, fps, nFrames, spin, ligand, highlightLig):
        """Assemble the ChimeraX .cxc command script. Silhouettes (the black
        outlines) are turned OFF and 'lighting soft' gives the clean VisualFactory
        look; 'view' auto-frames the camera (FindPerspective idea)."""
        prot = '#1 & protein'
        repCmds = []
        if style == 'surface':
            repCmds += ['hide #1 atoms', 'surface {}'.format(prot), 'color {} #6699cc'.format(prot)]
        elif style in ('cartoon', 'ribbon'):
            repCmds += ['hide #1 atoms', 'show {} cartoon'.format(prot), 'rainbow {}'.format(prot)]
        elif style == 'cartoon+sticks':
            repCmds += ['show {} cartoon'.format(prot), 'show {} atoms'.format(prot),
                        'style {} stick'.format(prot), 'rainbow {}'.format(prot)]
        elif style == 'sticks':
            repCmds += ['show #1 atoms', 'style #1 stick', 'rainbow {}'.format(prot)]

        ligCmds = []
        if highlightLig:
            lig = '#1 & :{}'.format(ligand)
            ligCmds += ['show {} atoms'.format(lig), 'style {} stick'.format(lig),
                        'color {} yellow'.format(lig), 'color {} byhetero'.format(lig)]

        # Numeric coordset range (the 'end' keyword is rejected). If the frame count
        # could not be determined, fall back to playing from the model's first frame.
        coordsetCmd = 'coordset #1 1,{}'.format(nFrames) if nFrames else 'coordset #1'
        spinCmds = ['turn y 2 180'] if spin else []   # 360 deg spin, recorded after playback
        lines = [
            'cd "{}"'.format(workDir),
            'set bgColor {}'.format(bg),
            'open "{}"'.format(structFile),
            'open "{}" structureModel #1'.format(trjFile),
            # DELETE (not hide) solvent/ions: tens of thousands of waters otherwise
            # fragment the protein surface into disconnected blobs (same fix as PyMol).
            'delete solvent',
            'delete ions',
        ]
        lines += repCmds + ligCmds
        lines += [
            'graphics silhouettes false',
            'lighting soft',
            'view',
            'movie record size {},{} supersample 3'.format(width, height),
            coordsetCmd,
        ]
        lines += spinCmds
        lines += [
            'wait',
            'movie encode "{}" framerate {} quality high'.format(outFile, fps),
        ]
        return '\n'.join(lines) + '\n'

    def _showMDTrajRMSDRMSF(self, paramName=None):
        """Handles both RMSD and RMSF (distinguished by mdAnalChoices text)."""
        system   = self.getMDSystem()
        analFlag = self.getEnumText('mdAnalChoices').lower()   # 'rmsd' or 'rmsf'
        selAtoms = self.getEnumText('selAtoms')
        args = (f'-i {system.getFileName()} -t {system.getTrajectoryFile()} '
                f'-o {system.getSystemName()} -{analFlag} -sa {selAtoms} ')
        if self.heavyAtoms.get():
            args += '-ha '
        Plugin.runScript(self, self._mdtrajScript, args, env=MDTRAJ_DIC, popen=True, wait=False)

    def _showMDTrajRGAnalysis(self, paramName=None):
        system   = self.getMDSystem()
        args = (f'-i {system.getFileName()} -t {system.getTrajectoryFile()} '
                f'-o {system.getSystemName()} -rg ')
        Plugin.runScript(self, self._mdtrajScript, args, env=MDTRAJ_DIC, popen=True, wait=False)

    def _showMDTrajSSAnalysis(self, paramName=None):
        system   = self.getMDSystem()
        ssFlags  = ['--per-residue', '--per-frame', '--heatmap']
        flag     = ssFlags[self.ssDisplayType.get()]
        args = (f'-i {system.getFileName()} -t {system.getTrajectoryFile()} {flag} ')
        Plugin.runScript(self, 'mdtraj_SS.py', args, env=MDTRAJ_DIC, popen=True, wait=False)

    def _showMDTrajSASAAnalysis(self, paramName=None):
        system   = self.getMDSystem()
        selAtoms = self.getEnumText('sasaScope')
        args = (f'-i {system.getFileName()} -t {system.getTrajectoryFile()} '
                f'-o {system.getSystemName()} -sa {selAtoms} -sasa ')
        Plugin.runScript(self, self._mdtrajScript, args, env=MDTRAJ_DIC, popen=True, wait=False)

    def _showMDTrajPCAAnalysis(self, paramName=None):
        system   = self.getMDSystem()
        pcaFlags = ['--pca-coord', '--pca-dist']
        flag     = pcaFlags[self.pcaType.get()]
        args = (f'-i {system.getFileName()} -t {system.getTrajectoryFile()} {flag} ')
        Plugin.runScript(self, 'mdtraj_PCA.py', args, env=MDTRAJ_DIC, popen=True, wait=False)

    def _showDistance(self, paramName=None):
        system       = self.getMDSystem()
        atom1, atom2 = self.atom1.get(), self.atom2.get()
        args = (f' -distance -i {system.getFileName()} -t {system.getTrajectoryFile()} '
                f'-o {system.getSystemName()} -a1 {atom1} -a2 {atom2}')
        Plugin.runScript(self, self._mdtrajScript, args, env=MDTRAJ_DIC, popen=True, wait=False)

    def _showProlifFp(self, paramName=None):
        fpPkl = self.getMDSystem().getProlifFile()
        Plugin.runScript(self, self._prolifViewScript, f'"{fpPkl}" --mode barcode',
                         env=MDTRAJ_DIC, popen=True, wait=False)

    def _showProlifNetwork(self, paramName=None):
        fpPkl = self.getMDSystem().getProlifFile()
        Plugin.runScript(self, self._prolifViewScript, f'"{fpPkl}" --mode network',
                         env=MDTRAJ_DIC, popen=True, wait=False)

    def _showProlifMatrix(self, paramName=None):
        fpPkl = self.getMDSystem().getProlifFile()
        Plugin.runScript(self, self._prolifViewScript, f'"{fpPkl}" --mode matrix',
                         env=MDTRAJ_DIC, popen=True, wait=False)



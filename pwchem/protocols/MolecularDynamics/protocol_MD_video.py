# **************************************************************************
# *
# * Authors: Joaquin Algorta (joaquin.algorta@cnb.csic.es)
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

import json
import os
import shutil

from pyworkflow.protocol import params
from pwem.protocols import EMProtocol
from pwem.objects.data import EMFile

from pwchem import Plugin
from pwchem.constants import MDTRAJ_DIC, OPENBABEL_DIC

# (UI label, script value) pairs, kept together so they can't drift out of sync.
_STYLE_OPTS = [
    ('Cartoon', 'cartoon'), ('Cartoon + sticks', 'cartoon+sticks'), ('Surface', 'surface'),
    ('Sticks', 'sticks'), ('Ribbon', 'ribbon'),
]
_BG_OPTS = [('White', 'white'), ('Black', 'black')]
_RESO_OPTS = [
    ('480p', '480p'), ('720p (HD)', '720p'), ('1080p (Full HD)', '1080p'), ('2160p (4K)', '4K'),
]
_COLORSCHEME_OPTS = [('Rainbow', 'rainbow'), ('By chain', 'chain'), ('Secondary structure', 'ss')]
_FORMAT_OPTS = [('mp4', 'mp4'), ('gif', 'gif')]
_AXIS_OPTS = [('X', 'x'), ('Y', 'y'), ('Z', 'z')]


class ProtocolMDVideo(EMProtocol):
    """Render an MD trajectory as a video (PyMOL), split into numberOfThreads
    single-threaded, resumable chunks instead of one unbounded viewer process."""
    _label = 'MD trajectory video'
    stepsExecutionMode = params.STEPS_PARALLEL

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        group = form.addGroup('Input')
        group.addParam('inputMDSystem', params.PointerParam, pointerClass='MDSystem',
                       label='MD system: ',
                       help='Molecular Dynamics System with an associated trajectory.')
        line = group.addLine('Frame range: ',
                             help='Render only this part of the trajectory (e.g. a binding event). '
                                  '0 in the last frame means "until the end".')
        line.addParam('firstFrame', params.IntParam, default=0, label='First')
        line.addParam('lastFrame', params.IntParam, default=0, label='Last')
        group.addParam('maxFrames', params.IntParam, default=0,
                       label='Max frames: ',
                       help='Cap the number of rendered frames by automatically increasing the '
                            'stride (0 = off). Use this instead of guessing a stride for long '
                            'trajectories.')

        group = form.addGroup('Representation')
        group.addParam('vidStyle', params.EnumParam, default=0,
                       label='Protein style: ', choices=self._labels(_STYLE_OPTS),
                       help='Representation used for the protein in the video.')
        group.addParam('colorScheme', params.EnumParam, default=0,
                       label='Color scheme: ', choices=self._labels(_COLORSCHEME_OPTS),
                       help='Rainbow (N->C spectrum) is unreadable for multi-chain systems; '
                            'use "By chain" or "Secondary structure" instead.')
        group.addParam('vidBg', params.EnumParam, default=0, label='Background: ',
                       choices=self._labels(_BG_OPTS), help='Background color of the rendered frames.')

        group = form.addGroup('Ligand & pocket')
        group.addParam('vidHighlightLig', params.BooleanParam, default=True,
                       label='Highlight ligand: ',
                       help='Show the ligand as coloured sticks.')
        group.addParam('showPocket', params.BooleanParam, default=False,
                       label='Show pocket residues: ', condition='vidHighlightLig',
                       help='Show residues within the pocket distance of the ligand as thin '
                            'sticks -- the usual view for binding studies.')
        group.addParam('pocketDist', params.FloatParam, default=5.0,
                       label='Pocket distance (A): ', condition='vidHighlightLig and showPocket',
                       help='Residues with an atom within this distance of the ligand are shown.')

        form.addSection(label='Video')
        group = form.addGroup('Camera')
        group.addParam('alignTraj', params.BooleanParam, default=True,
                       label='Align trajectory: ',
                       help='intra_fit on the protein backbone to remove rotation/translation '
                            'drift. Without it the molecule wanders out of the camera.')
        group.addParam('vidSpin', params.BooleanParam, default=False,
                       label='Camera spin: ',
                       help='Add a camera rotation across the trajectory.')
        group.addParam('spinAxis', params.EnumParam, default=1, choices=self._labels(_AXIS_OPTS),
                       label='Spin axis: ', condition='vidSpin')
        group.addParam('nTurns', params.FloatParam, default=1.0,
                       label='Number of turns: ', condition='vidSpin')

        group = form.addGroup('Output')
        group.addParam('vidResolution', params.EnumParam, default=1,
                       label='Resolution: ', choices=self._labels(_RESO_OPTS),
                       help='Frame resolution. Higher resolutions render more slowly.')
        group.addParam('vidFps', params.IntParam, default=15, label='Frames per second: ',
                       help='Playback speed of the resulting video.')
        group.addParam('vidStride', params.IntParam, default=1, label='Frame stride: ',
                       help='Render every Nth trajectory frame (use >1 for long trajectories).')
        group.addParam('vidFormat', params.EnumParam, default=0, label='Output format: ',
                       choices=self._labels(_FORMAT_OPTS), help='Video container. mp4 needs ffmpeg; '
                                                            'otherwise a gif is produced.')
        group.addParam('timeLabel', params.BooleanParam, default=False,
                       label='Overlay time label: ',
                       help='Overlay "t = X ns" on each frame.')
        group.addParam('vidSmooth', params.IntParam, default=0,
                       label='Smoothing window: ',
                       help='Coordinate-smoothing window applied before rendering (0 = none). '
                            'Removes thermal jitter for cleaner playback (visualization only).')
        group.addParam('vidRay', params.BooleanParam, default=True, expertLevel=params.LEVEL_ADVANCED,
                       label='Cinematic ray-tracing: ',
                       help='Ray-trace each frame with soft shadows. Disable for a faster, '
                            'flatter preview.')
        group.addParam('videoQuality', params.IntParam, default=20,
                       label='Quality (CRF): ', condition='vidFormat==0',
                       help='mp4 encoding quality (libx264 CRF). Lower = better quality, bigger file.')
        group.addParam('keepFrames', params.BooleanParam, default=False, expertLevel=params.LEVEL_ADVANCED,
                       label='Keep frame PNGs: ',
                       help='Keep the rendered PNG frames (extra/frames/) after encoding, e.g. '
                            'for figures or slides.')

        form.addParallelSection(threads=4, mpi=0)

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        pStep = self._insertFunctionStep(self.prepareStep)
        rSteps = [self._insertFunctionStep(self.renderChunkStep, i, prerequisites=[pStep])
                  for i in range(self._getNChunks())]
        eStep = self._insertFunctionStep(self.encodeStep, prerequisites=rSteps)
        self._insertFunctionStep(self.createOutputStep, prerequisites=[eStep])

    def prepareStep(self):
        """Count frames, apply first/last/stride/maxFrames, split into chunks."""
        structFile, trjFile = self._getStructTrajFiles()

        countFile = os.path.abspath(self._getExtraPath('mdvideo_count.json'))
        args = '-i "{}" -t "{}" -o "{}"'.format(structFile, trjFile, countFile)
        Plugin.runScript(self, 'count_MDvideo_frames.py', args, env=MDTRAJ_DIC)

        with open(countFile) as f:
            countInfo = json.load(f)
        nTotalTraj, timesPs = countInfo['nTotal'], countInfo['timesPs']

        first = max(0, self.firstFrame.get())
        last = min(self.lastFrame.get() or nTotalTraj, nTotalTraj)
        stride = max(1, self.vidStride.get())
        candidateIdx = list(range(first, last, stride))

        maxFrames = self.maxFrames.get()
        if maxFrames and len(candidateIdx) > maxFrames:
            extraStride = -(-len(candidateIdx) // maxFrames)   # ceil division
            candidateIdx = candidateIdx[::extraStride]

        states = [i + 1 for i in candidateIdx]   # 1-based PyMOL/mdtraj state numbers
        times = {str(gIdx): timesPs[i] / 1000.0 for gIdx, i in enumerate(candidateIdx)}

        chunks, offsets = self._splitChunks(states, self._getNChunks())

        with open(self._getPlanFile(), 'w') as f:
            json.dump({'states': states, 'chunks': chunks, 'offsets': offsets,
                      'times': times, 'nTotal': len(states)}, f)

    def renderChunkStep(self, i):
        plan = self._readPlan()
        states = plan['chunks'][i]
        if not states:
            return   # trajectory shorter than numberOfThreads: nothing for this chunk

        structFile, trjFile = self._getStructTrajFiles()
        chunkDir = os.path.abspath(self._getExtraPath('frames', 'chunk_{:02d}'.format(i)))
        os.makedirs(chunkDir, exist_ok=True)

        args = self._buildRenderArgs(structFile, trjFile, states, plan['offsets'][i],
                                     plan['nTotal'], chunkDir)
        self._runVideoScript(args)

    def encodeStep(self):
        plan = self._readPlan()
        outBase = self._getOutBase()
        extraDir = os.path.abspath(self._getExtraPath())
        framesGlob = os.path.abspath(self._getExtraPath('frames', 'chunk_*', 'frame_*.png'))
        fmt = self._selected(_FORMAT_OPTS, self.vidFormat)

        args = (f'-- -o "{outBase}" --workdir "{extraDir}" --encodeOnly --framesGlob "{framesGlob}" '
               f'--fps {self.vidFps.get()} --format {fmt} --crf {self.videoQuality.get()} '
               f'--keepFrames {int(self.keepFrames.get())}')

        if self.timeLabel.get():
            timesFile = os.path.abspath(self._getExtraPath('mdvideo_times.json'))
            with open(timesFile, 'w') as f:
                json.dump(plan['times'], f)
            args += ' --timeLabel 1 --timesFile "{}"'.format(timesFile)

        self._runVideoScript(args)

    def createOutputStep(self):
        outBase = self._getOutBase()
        fmt = self._selected(_FORMAT_OPTS, self.vidFormat)
        videoFile = self._getExtraPath('{}.{}'.format(outBase, fmt))
        if not os.path.exists(videoFile):
            # ffmpeg missing -> the script already fell back to gif
            videoFile = self._getExtraPath('{}.gif'.format(outBase))

        outVideo = EMFile(filename=os.path.abspath(videoFile))
        self._defineOutputs(outputVideo=outVideo)
        self._defineSourceRelation(self.inputMDSystem, outVideo)

    # --------------------------- UTILS functions -----------------------------------
    def _getStructTrajFiles(self):
        mdsys = self.inputMDSystem.get()
        topFile = mdsys.getTopologyFile()
        # .top/.tpr (GROMACS) aren't coordinates; other topologies (e.g. AMBER .parm7) match the trajectory.
        if not topFile or os.path.splitext(topFile)[1].lower() in ('.top', '.tpr'):
            structFile = mdsys.getSystemFile()
        else:
            structFile = topFile
        return os.path.abspath(structFile), os.path.abspath(mdsys.getTrajectoryFile())

    def _getOutBase(self):
        return '{}_MDvideo'.format(self.inputMDSystem.get().getSystemName())

    def _getNChunks(self):
        # 1 thread coordinates; the rest render chunks concurrently (1 core each).
        return max(1, self.numberOfThreads.get() - 1)

    def _getPlanFile(self):
        return self._getExtraPath('mdvideo_plan.json')

    def _readPlan(self):
        with open(self._getPlanFile()) as f:
            return json.load(f)

    def _runVideoScript(self, args):
        Plugin.runScript(self, 'create_MDvideo.py', args, env=OPENBABEL_DIC, pyStr='pymol -cq')

    def _labels(self, opts):
        return [label for label, _ in opts]

    def _values(self, opts):
        return [value for _, value in opts]

    def _selected(self, opts, enumParam):
        """The script value of an EnumParam's currently selected option."""
        return self._values(opts)[enumParam.get()]

    def _splitChunks(self, states, nChunks):
        """Split states into nChunks contiguous, near-equal groups; return (chunks, offsets)."""
        baseSize, rem = divmod(len(states), nChunks)
        chunks, offsets, pos = [], [], 0
        for i in range(nChunks):
            size = baseSize + (1 if i < rem else 0)
            offsets.append(pos)
            chunks.append(states[pos:pos + size])
            pos += size
        return chunks, offsets

    def _buildRenderArgs(self, structFile, trjFile, states, frameOffset, nTotal, chunkDir):
        # "--" marks everything after it as script args, not pymol's own.
        ligandId = self.inputMDSystem.get().getLigandID()
        style = self._selected(_STYLE_OPTS, self.vidStyle)
        bg = self._selected(_BG_OPTS, self.vidBg)
        colorScheme = self._selected(_COLORSCHEME_OPTS, self.colorScheme)
        resolution = self._selected(_RESO_OPTS, self.vidResolution)
        spinAxis = self._selected(_AXIS_OPTS, self.spinAxis)
        statesStr = ','.join(map(str, states))

        return (f'-- -i "{structFile}" -t "{trjFile}" --renderOnly --framesDir "{chunkDir}" '
               f'--states "{statesStr}" --frameOffset {frameOffset} --nTotalFrames {nTotal} '
               f'--threads 1 --style {style} --bg {bg} --colorScheme {colorScheme} '
               f'--ligand "{ligandId}" --highlightLig {int(self.vidHighlightLig.get())} '
               f'--showPocket {int(self.showPocket.get())} --pocketDist {self.pocketDist.get()} '
               f'--resolution {resolution} --smooth {max(0, self.vidSmooth.get())} '
               f'--alignTraj {int(self.alignTraj.get())} --ray {int(self.vidRay.get())} '
               f'--spin {int(self.vidSpin.get())} --spinAxis {spinAxis} --nTurns {self.nTurns.get()}')

    # --------------------------- INFO functions -----------------------------------
    def _summary(self):
        summary = []
        if hasattr(self, 'outputVideo'):
            plan = self._readPlan() if os.path.exists(self._getPlanFile()) else {}
            nFrames = plan.get('nTotal')
            fps = self.vidFps.get()
            duration = '{:.1f}'.format(nFrames / fps) if nFrames and fps else '?'
            resolution = self._selected(_RESO_OPTS, self.vidResolution)
            summary.append('Rendered {} frames at {} ({} fps, ~{} s) with PyMOL.'
                           .format(nFrames, resolution, fps, duration))
            summary.append('Output video: {}'.format(os.path.abspath(self.outputVideo.getFileName())))
        return summary

    def _validate(self):
        errors = []
        mdsys = self.inputMDSystem.get()
        if mdsys and not mdsys.hasTrajectory():
            errors.append('The input MDSystem has no associated trajectory file.')
        if self.lastFrame.get() and self.firstFrame.get() >= self.lastFrame.get():
            errors.append('First frame must be smaller than last frame.')
        if self.vidFps.get() < 1:
            errors.append('Frames per second must be at least 1.')
        if self.vidStride.get() < 1:
            errors.append('Frame stride must be at least 1.')
        if self.maxFrames.get() < 0:
            errors.append('Max frames cannot be negative.')
        if self.firstFrame.get() < 0:
            errors.append('First frame cannot be negative.')
        if self.lastFrame.get() < 0:
            errors.append('Last frame cannot be negative.')
        return errors

    def _warnings(self):
        warnings = []
        if self.vidFormat.get() == 0 and not shutil.which('ffmpeg'):
            warnings.append('ffmpeg was not found: the video will be encoded as an animated gif '
                            'instead of mp4.')
        return warnings

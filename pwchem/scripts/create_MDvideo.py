# **************************************************************************
# *
# * Authors:     Scipion-Chem team (scipionchem@cnb.csic.es)
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 3 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# **************************************************************************

"""
create_MDvideo.py  --  Cinematic Molecular Dynamics trajectory video generator.

This script renders a video (mp4 or animated gif) of an MD trajectory using the
PyMOL python interpreter in headless mode.  It is meant to be launched by the
``MDSystemPViewer`` of scipion-chem, e.g. ::

    pymol -cq create_MDvideo.py -- -i system.gro -t traj.xtc -o movie [options]

The design borrows ideas from *TheVisualHub/VisualFactory* (see
``claude/extra/VisualFactory.md``):

  * **FindPerspective** -> the camera is auto-framed on the molecular centroid
    (``orient`` + ``zoom``), so the user does not have to point the camera.
  * **UltimateSmoothMD** -> an optional coordinate-smoothing window removes the
    thermal "boiling" of the raw trajectory for clean, honest playback.  The
    smoothed coordinates are used **for visualization only**.
  * Cinematic styling -> ray-traced frames with ambient occlusion, soft
    shadows and a clean background, plus an optional slow camera spin.
  * **VideoProcessing** -> frames are muxed into a web-friendly, ``+faststart``
    H.264 mp4 with ffmpeg (falling back to an animated gif via Pillow when
    ffmpeg is not available).

Rendering always uses ``cmd.ray`` so that it works on a head-less node with no
OpenGL context.

NOTE: educational / non-commercial use, consistent with the VisualFactory and
UCSF PyMOL licensing spirit.
"""

import argparse
import glob
import os
import shutil
import subprocess
import sys

from pymol import cmd


# --------------------------------------------------------------------------- #
#  Presets
# --------------------------------------------------------------------------- #
# Resolution presets (width, height).  "4K" follows VisualFactory's 4K output.
RESOLUTIONS = {
    '480p':  (854, 480),
    '720p':  (1280, 720),
    '1080p': (1920, 1080),
    '4K':    (3840, 2160),
}

# Atoms that are ions / counter-ions and should not be treated as "the ligand".
ION_RESNAMES = ('NA', 'CL', 'SOD', 'CLA', 'K', 'MG', 'ZN', 'CA', 'POT', 'BR', 'IOD')


def parseArgs():
    p = argparse.ArgumentParser(description='Render a cinematic MD trajectory video with PyMOL.')
    p.add_argument('-i', '--inputStruct', required=True,
                   help='System / topology structure file (.pdb, .gro, .prmtop ...).')
    p.add_argument('-t', '--trajectory', required=True,
                   help='Trajectory file (.xtc, .trr, .dcd, .nc/.netcdf ...).')
    p.add_argument('-o', '--output', default='md_video',
                   help='Output basename (without extension).')
    p.add_argument('--workdir', default=None,
                   help='Directory where the video and temporary frames are written.')

    # Visual style
    p.add_argument('--style', default='cartoon',
                   choices=['cartoon', 'surface', 'sticks', 'ribbon', 'cartoon+sticks'],
                   help='Protein representation.')
    p.add_argument('--bg', default='white', choices=['white', 'black'],
                   help='Background color.')
    p.add_argument('--ligand', default='LIG',
                   help='Residue name of the ligand to highlight (default LIG).')
    p.add_argument('--highlightLig', type=int, default=1,
                   help='1: show the ligand as sticks and colour it; 0: off.')
    p.add_argument('--ray', type=int, default=1,
                   help='1: cinematic ray-tracing (shadows + ambient occlusion); '
                        '0: faster, flatter rendering.')
    p.add_argument('--spin', type=int, default=0,
                   help='1: add a full 360 deg camera spin across the trajectory.')

    # Trajectory handling
    p.add_argument('--stride', type=int, default=1,
                   help='Render every Nth frame (>=1).')
    p.add_argument('--smooth', type=int, default=0,
                   help='Coordinate smoothing window (0 = no smoothing). '
                        'Higher = smoother but less faithful (visualization only).')

    # Video parameters
    p.add_argument('--resolution', default='720p', choices=list(RESOLUTIONS.keys()),
                   help='Frame resolution.')
    p.add_argument('--fps', type=int, default=15, help='Frames per second.')
    p.add_argument('--format', default='mp4', choices=['mp4', 'gif'],
                   help='Output video container.')
    p.add_argument('--open', type=int, default=0,
                   help='1: open the finished video with the default system player.')

    # PyMOL passes the script name + everything after "--" in sys.argv.
    return p.parse_args(sys.argv[1:])


def log(msg):
    """Smart console feedback (VisualFactory style)."""
    print('[create_MDvideo] {}'.format(msg), flush=True)


def loadSystem(structFile, trajFile):
    """Load the topology and the trajectory into a single object 'mdsys'."""
    obj = 'mdsys'
    cmd.load(structFile, obj)
    nTopoStates = cmd.count_states(obj)
    # A multi-model topology would confuse load_traj: keep only the 1st state.
    if nTopoStates > 1:
        cmd.split_states(obj, 1, 1)
        cmd.delete(obj)
        cmd.set_name(cmd.get_object_list()[0], obj)

    ext = os.path.splitext(trajFile)[1].lower()
    kwargs = {}
    # AMBER ascii trajectories need an explicit format hint.
    if ext in ('.crd', '.mdcrd', '.trj'):
        kwargs['format'] = 'trj'
    cmd.load_traj(trajFile, obj, **kwargs)

    # Physically REMOVE solvent + ions (do not just hide them). In solvated MD
    # systems the tens of thousands of water atoms make PyMOL's molecular-surface
    # generation fragment the protein into disconnected "bubbles"; hiding the water
    # is not enough, it must be removed. Removing it also speeds up smoothing and
    # ray-tracing. The ligand (highlighted separately) is NOT in ION_RESNAMES.
    cmd.remove('({}) and ((solvent) or (resn {}))'.format(obj, '+'.join(ION_RESNAMES)))
    cmd.rebuild()
    return obj, cmd.count_states(obj)


def setupLighting(bg):
    """Known-good, photorealistic render settings.

    These values are deliberately conservative (the look you get from PyMOL's
    publication presets) so the output is reliably clean:
      * NO ``ray_trace_mode`` outline  -> avoids the black edges around cartoons.
      * depth-cue and ray fog OFF       -> avoids the rear of the model going dark.
      * soft shadows + mild specular    -> cinematic but not noisy.
    """
    cmd.bg_color(bg)
    # White background must be opaque in the frames (so mp4/gif show white, not black).
    cmd.set('ray_opaque_background', 1)

    cmd.set('ray_trace_mode', 0)            # photorealistic, NO outline (was the black-edge cause)
    cmd.set('antialias', 2)
    cmd.set('hash_max', 240)

    # Balanced lighting: bright enough that nothing reads as "black".
    cmd.set('ambient', 0.45)
    cmd.set('direct', 0.55)
    cmd.set('reflect', 0.40)
    cmd.set('light_count', 2)
    cmd.set('specular', 0.25)
    cmd.set('spec_count', 1)
    cmd.set('shininess', 10)
    cmd.set('ray_shadows', 1)
    cmd.set('ray_shadow_decay_factor', 0.1)

    # Disable distance darkening / fog -> the whole model stays evenly lit.
    cmd.set('depth_cue', 0)
    cmd.set('ray_trace_fog', 0)
    cmd.set('fog', 0)

    # Smooth, high quality cartoons and surfaces.
    cmd.set('cartoon_fancy_helices', 1)
    cmd.set('cartoon_smooth_loops', 1)
    cmd.set('cartoon_highlight_color', -1)
    cmd.set('surface_quality', 1)
    cmd.set('solvent_radius', 1.4)
    cmd.set('transparency_mode', 1)


def selectLigand(obj, ligand):
    """Return a selection string for the ligand (named resn, else any organic het)."""
    ligSel = '{0} and (resn {1}) and not polymer and not solvent'.format(obj, ligand)
    if cmd.count_atoms(ligSel) == 0:
        ligSel = ('{0} and not polymer and not solvent and not resn {1}'
                  .format(obj, '+'.join(ION_RESNAMES)))
    return ligSel if cmd.count_atoms(ligSel) > 0 else None


def applyStyle(obj, style, bg, ligand, highlightLig, spin=False):
    """Apply a known-good representation preset, then auto-frame the camera."""
    setupLighting(bg)
    cmd.hide('everything', obj)   # solvent/ions were already removed in loadSystem

    polymerSel = '{} and polymer'.format(obj)

    if style == 'surface':
        # Recommended preset: a clean molecular surface over the protein + the
        # ligand as bright sticks. The surface is kept OPAQUE on purpose: a
        # transparent surface forces PyMOL to ray-trace several layers per frame,
        # which is far too slow for a multi-frame movie.
        cmd.show('surface', polymerSel)
        cmd.color('skyblue', polymerSel)
        cmd.set('surface_quality', 0)      # 0 is plenty for a movie and much faster
    elif style in ('cartoon', 'cartoon+sticks', 'ribbon'):
        rep = 'ribbon' if style == 'ribbon' else 'cartoon'
        cmd.show(rep, polymerSel)
        if style == 'cartoon+sticks':
            cmd.set('cartoon_side_chain_helper', 1)
            cmd.show('sticks', polymerSel + ' and sidechain')
            cmd.set('stick_radius', 0.15, polymerSel)
            # Colour every atom (not only CA) so the sticks share the cartoon's
            # N->C rainbow instead of reading as a separate monochrome mesh.
            cmd.spectrum('count', 'rainbow', polymerSel)
        else:
            # Spectrum N->C is a classic, reliably attractive colouring.
            cmd.spectrum('count', 'rainbow', polymerSel + ' and name CA')
    elif style == 'sticks':
        cmd.show('sticks', polymerSel)
        cmd.set('stick_radius', 0.18, polymerSel)
        cmd.spectrum('count', 'rainbow', polymerSel)

    # Highlight the ligand (the part that usually matters in a binding study).
    if highlightLig:
        ligSel = selectLigand(obj, ligand)
        if ligSel:
            cmd.show('sticks', ligSel)
            cmd.set('stick_radius', 0.20, ligSel)
            cmd.color('yellow', ligSel)
            cmd.util.cnc(ligSel)               # colour heteroatoms by element
            cmd.set('stick_ball', 0)
            log('Highlighted {} ligand atoms.'.format(cmd.count_atoms(ligSel)))

    # Auto-frame the camera on the protein (+ ligand) ONLY -- never on the full
    # object, whose hidden solvent box would shrink the molecule into a corner
    # (FindPerspective auto-camera idea).
    focusSel = polymerSel if cmd.count_atoms(polymerSel) else obj
    ligSel = selectLigand(obj, ligand)
    if ligSel:
        focusSel = '({}) or ({})'.format(focusSel, ligSel)
    cmd.orient(focusSel)
    frameCamera(focusSel, surface=(style == 'surface'), spin=spin)


def frameCamera(sel, surface=False, spin=False, margin=3.0):
    """Zoom so the model fits the frame -- and *stays* in frame while spinning.

    ``cmd.zoom(complete=1)`` only fits the current orientation's bounding box, so
    when the camera spins (or for a molecular surface that bulges past the atom
    centres) parts of the model get clipped. The model rotates about its centre,
    so the worst-case projected radius is the bounding-sphere radius ``R``. When
    ``spin`` is on we add a lateral ``buffer = R - largestHalfDimension`` so the
    view is wide enough for the diagonal at every angle; when it is off we fit the
    box tightly so the molecule fills the frame. In both cases the clipping slab is
    opened to ``4R`` so the near/far planes never slice the model.
    """
    (x0, y0, z0), (x1, y1, z1) = cmd.get_extent(sel)
    dims = (x1 - x0, y1 - y0, z1 - z0)
    radius = 0.5 * (dims[0] ** 2 + dims[1] ** 2 + dims[2] ** 2) ** 0.5
    buffer = margin
    if spin:
        buffer += radius - max(dims) / 2.0     # widen so the diagonal fits at any angle
    if surface:
        buffer += 2.0                          # the surface bulges past the atom centres

    cmd.zoom(sel, buffer=buffer, complete=1)
    cmd.clip('slab', radius * 4)


def renderFramesByState(obj, nStates, args, framesDir):
    """Render explicitly per state (set state BEFORE ray-tracing)."""
    width, height = RESOLUTIONS[args.resolution]
    states = list(range(1, nStates + 1, max(1, args.stride)))
    nFrames = len(states)
    spinPerFrame = (360.0 / nFrames) if (args.spin and nFrames > 1) else 0.0

    if not args.ray:
        cmd.set('ray_shadows', 0)
        cmd.set('ambient_occlusion_mode', 0)

    log('Rendering {} frames at {}x{} (ray={}, spin={}).'
        .format(nFrames, width, height, args.ray, bool(args.spin)))

    frameFiles = []
    for idx, state in enumerate(states):
        cmd.set('state', state)
        if spinPerFrame:
            cmd.turn('y', spinPerFrame)
        framePath = os.path.join(framesDir, 'frame_{:05d}.png'.format(idx))
        cmd.ray(width, height)
        cmd.png(framePath, dpi=300 if args.ray else 150)
        frameFiles.append(framePath)
        if (idx + 1) % 10 == 0 or idx == nFrames - 1:
            log('  ...{}/{} frames'.format(idx + 1, nFrames))
    return frameFiles


def assembleVideo(frameFiles, framesDir, outBase, fps, fmt):
    """Mux PNG frames into mp4 (ffmpeg) or gif (Pillow)."""
    if not frameFiles:
        raise RuntimeError('No frames were rendered.')

    ffmpeg = shutil.which('ffmpeg')
    if fmt == 'mp4' and ffmpeg:
        outFile = outBase + '.mp4'
        pattern = os.path.join(framesDir, 'frame_%05d.png')
        cmdLine = [
            ffmpeg, '-y', '-framerate', str(fps), '-i', pattern,
            '-c:v', 'libx264', '-preset', 'slow', '-crf', '20',
            '-pix_fmt', 'yuv420p', '-movflags', '+faststart',
            # H.264 requires even dimensions.
            '-vf', 'scale=trunc(iw/2)*2:trunc(ih/2)*2',
            outFile,
        ]
        log('Encoding mp4 with ffmpeg...')
        subprocess.check_call(cmdLine)
        return outFile

    # Fallback (or explicit request): animated gif via Pillow.
    if fmt == 'mp4' and not ffmpeg:
        log('ffmpeg not found -> falling back to animated gif.')
    from PIL import Image
    outFile = outBase + '.gif'
    log('Encoding gif with Pillow...')
    frames = [Image.open(f).convert('RGB') for f in frameFiles]
    duration = int(1000.0 / max(1, fps))
    frames[0].save(outFile, save_all=True, append_images=frames[1:],
                   duration=duration, loop=0, optimize=True)
    return outFile


def main():
    args = parseArgs()

    workdir = args.workdir or os.path.dirname(os.path.abspath(args.trajectory))
    os.makedirs(workdir, exist_ok=True)
    framesDir = os.path.join(workdir, '_md_video_frames')
    if os.path.isdir(framesDir):
        shutil.rmtree(framesDir)
    os.makedirs(framesDir)
    outBase = os.path.join(workdir, args.output)

    log('System    : {}'.format(args.inputStruct))
    log('Trajectory: {}'.format(args.trajectory))

    cmd.feedback('disable', 'all', 'everything')
    obj, nStates = loadSystem(args.inputStruct, args.trajectory)
    log('Loaded {} states.'.format(nStates))

    if args.smooth > 0 and nStates > 2:
        log('Smoothing trajectory (window={}, visualization only).'.format(args.smooth))
        cmd.smooth('all', passes=1, window=args.smooth)

    applyStyle(obj, args.style, args.bg, args.ligand, bool(args.highlightLig), spin=bool(args.spin))

    frameFiles = renderFramesByState(obj, nStates, args, framesDir)
    outFile = assembleVideo(frameFiles, framesDir, outBase, args.fps, args.format)

    # Clean up the intermediate frames.
    shutil.rmtree(framesDir, ignore_errors=True)
    log('DONE. Video written to: {}'.format(outFile))

    if args.open:
        opener = shutil.which('xdg-open') or shutil.which('open')
        if opener:
            subprocess.Popen([opener, outFile])


# PyMOL execs scripts with __name__ == 'pymol' (not '__main__'), so accept both.
if __name__ in ('__main__', 'pymol'):
    main()

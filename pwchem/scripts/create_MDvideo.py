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
create_MDvideo.py -- headless PyMOL renderer for an MD trajectory video.

Launched per-chunk by ProtocolMDVideo: ``pymol -cq create_MDvideo.py -- -i sys.gro
-t traj.xtc -o movie [options]``. ``--renderOnly``/``--encodeOnly`` pick a phase;
with neither, it does both (for manual/standalone use).
"""

import argparse
import glob
import json
import os
import shutil
import subprocess
import sys

from pymol import cmd


# Resolution presets (width, height).
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
    p.add_argument('-i', '--inputStruct', default=None,
                   help='System / topology structure file (.pdb, .gro, .prmtop ...).')
    p.add_argument('-t', '--trajectory', default=None,
                   help='Trajectory file (.xtc, .trr, .dcd, .nc/.netcdf ...).')
    p.add_argument('-o', '--output', default='md_video',
                   help='Output basename (without extension).')
    p.add_argument('--workdir', default=None,
                   help='Directory where the video is written (default: trajectory dir).')
    p.add_argument('--framesDir', default=None,
                   help='Directory this chunk writes its frame PNGs to '
                        '(default: <workdir>/_md_video_frames).')
    p.add_argument('--framesGlob', default=None,
                   help='[encode] Glob pattern collecting every chunk\'s frames, '
                        'e.g. "extra/frames/chunk_*/frame_*.png".')

    # Phase selection
    p.add_argument('--renderOnly', action='store_true', help='Only render this chunk\'s frames.')
    p.add_argument('--encodeOnly', action='store_true', help='Only encode already-rendered frames.')

    # Trajectory chunk handled by this process
    p.add_argument('--states', default=None,
                   help='States to render: a comma list ("1,5,9") or "start:end:step". '
                        'Default: every --stride-th of all trajectory states.')
    p.add_argument('--frameOffset', type=int, default=0,
                   help='Global (0-based) index of the first state in --states, used for '
                        'output filenames and for an absolute, chunk-independent camera spin.')
    p.add_argument('--nTotalFrames', type=int, default=None,
                   help='Total number of frames in the whole video (spin angle denominator). '
                        'Defaults to len(states).')
    p.add_argument('--threads', type=int, default=0,
                   help='Cap PyMOL to this many threads (0 = do not restrict). Used so a '
                        'multi-chunk render does not consume all CPU cores.')

    # Visual style
    p.add_argument('--style', default='cartoon',
                   choices=['cartoon', 'surface', 'sticks', 'ribbon', 'cartoon+sticks'],
                   help='Protein representation.')
    p.add_argument('--bg', default='white', choices=['white', 'black'],
                   help='Background color.')
    p.add_argument('--colorScheme', default='rainbow', choices=['rainbow', 'chain', 'ss'],
                   help='Protein coloring: rainbow (N->C spectrum), chain, or secondary structure.')
    p.add_argument('--ligand', default='LIG',
                   help='Residue name of the ligand to highlight (default LIG).')
    p.add_argument('--highlightLig', type=int, default=1,
                   help='1: show the ligand as sticks and colour it; 0: off.')
    p.add_argument('--showPocket', type=int, default=0,
                   help='1: show pocket residues (within --pocketDist of the ligand) as thin sticks.')
    p.add_argument('--pocketDist', type=float, default=5.0,
                   help='Pocket distance cutoff in Angstrom (only used with --showPocket).')
    p.add_argument('--ray', type=int, default=1,
                   help='1: cinematic ray-tracing (shadows + ambient occlusion); '
                        '0: faster, flatter rendering.')
    p.add_argument('--spin', type=int, default=0,
                   help='1: add a camera spin across the trajectory.')
    p.add_argument('--spinAxis', default='y', choices=['x', 'y', 'z'], help='Camera spin axis.')
    p.add_argument('--nTurns', type=float, default=1.0, help='Number of full spin rotations.')

    # Trajectory handling
    p.add_argument('--stride', type=int, default=1,
                   help='[standalone] Render every Nth frame (>=1); ignored when --states is given.')
    p.add_argument('--smooth', type=int, default=0,
                   help='Coordinate smoothing window (0 = no smoothing). '
                        'Higher = smoother but less faithful (visualization only).')
    p.add_argument('--alignTraj', type=int, default=0,
                   help='1: intra_fit on the protein backbone before rendering, to remove '
                        'the overall rotation/translation drift.')

    # Video parameters
    p.add_argument('--resolution', default='720p', choices=list(RESOLUTIONS.keys()),
                   help='Frame resolution.')
    p.add_argument('--fps', type=int, default=15, help='Frames per second.')
    p.add_argument('--format', default='mp4', choices=['mp4', 'gif'],
                   help='Output video container.')
    p.add_argument('--crf', type=int, default=20, help='mp4 encoding quality (libx264 CRF; lower=better).')
    p.add_argument('--timeLabel', type=int, default=0,
                   help='1: overlay "t = X ns" on each frame (needs --timesFile).')
    p.add_argument('--timesFile', default=None,
                   help='JSON file {globalFrameIdx: timeNs} used by --timeLabel.')
    p.add_argument('--keepFrames', type=int, default=0,
                   help='1: keep the rendered PNG frames after encoding.')

    # PyMOL passes the script name + everything after "--" in sys.argv.
    return p.parse_args(sys.argv[1:])


def log(msg):
    print('[create_MDvideo] {}'.format(msg), flush=True)


def parseStatesArg(statesStr):
    """Parse --states: a comma list ("1,5,9") or a "start:end:step" range."""
    if ':' in statesStr:
        parts = [int(x) for x in statesStr.split(':')]
        start, end = parts[0], parts[1]
        step = parts[2] if len(parts) > 2 else 1
        return list(range(start, end, step))
    return [int(x) for x in statesStr.split(',') if x.strip() != '']


def resolveStates(args, nStates):
    if args.states:
        return parseStatesArg(args.states)
    return list(range(1, nStates + 1, max(1, args.stride)))


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

    # Remove (not hide) solvent/ions: otherwise waters fragment the surface.
    cmd.remove('({}) and ((solvent) or (resn {}))'.format(obj, '+'.join(ION_RESNAMES)))
    cmd.rebuild()
    return obj, cmd.count_states(obj)


def setupLighting(bg):
    """Conservative, publication-style lighting: no outline, no depth fog, soft shadows."""
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


def colorByScheme(polymerSel, scheme, style):
    """rainbow / chain / ss (secondary structure) -- rainbow alone is unreadable
    for multi-chain systems."""
    if scheme == 'chain':
        cmd.util.cbc(polymerSel)
        return
    if scheme == 'ss':
        cmd.dss(polymerSel)
        cmd.color('red', polymerSel + ' and ss H')
        cmd.color('yellow', polymerSel + ' and ss S')
        cmd.color('cyan', polymerSel + ' and not (ss H or ss S)')
        return

    # rainbow: spectrum N->C, a classic reliably attractive colouring.
    if style in ('cartoon+sticks', 'sticks'):
        cmd.spectrum('count', 'rainbow', polymerSel)
    else:
        cmd.spectrum('count', 'rainbow', polymerSel + ' and name CA')


def applyStyle(obj, args):
    """Apply a known-good representation preset, then auto-frame the camera."""
    setupLighting(args.bg)
    cmd.hide('everything', obj)   # solvent/ions were already removed in loadSystem

    polymerSel = '{} and polymer'.format(obj)

    if args.style == 'surface':
        # Opaque on purpose: a transparent surface ray-traces multiple layers per frame.
        cmd.show('surface', polymerSel)
        cmd.color('skyblue', polymerSel)
        cmd.set('surface_quality', 0)      # 0 is plenty for a movie and much faster
    elif args.style in ('cartoon', 'cartoon+sticks', 'ribbon'):
        rep = 'ribbon' if args.style == 'ribbon' else 'cartoon'
        cmd.show(rep, polymerSel)
        if args.style == 'cartoon+sticks':
            cmd.set('cartoon_side_chain_helper', 1)
            cmd.show('sticks', polymerSel + ' and sidechain')
            cmd.set('stick_radius', 0.15, polymerSel)
    elif args.style == 'sticks':
        cmd.show('sticks', polymerSel)
        cmd.set('stick_radius', 0.18, polymerSel)

    if args.style != 'surface':
        colorByScheme(polymerSel, args.colorScheme, args.style)

    # Highlight the ligand (the part that usually matters in a binding study).
    ligSel = selectLigand(obj, args.ligand)
    if args.highlightLig and ligSel:
        cmd.show('sticks', ligSel)
        cmd.set('stick_radius', 0.20, ligSel)
        cmd.color('purple', ligSel)         # yellow clashes with ss/rainbow; black is invisible on --bg black
        cmd.util.cnc(ligSel)               # colour heteroatoms by element
        cmd.set('stick_ball', 0)
        log('Highlighted {} ligand atoms.'.format(cmd.count_atoms(ligSel)))

        if args.showPocket:
            pocketSel = ('byres (({}) within {} of ({}))'
                        .format(polymerSel, args.pocketDist, ligSel))
            cmd.show('sticks', pocketSel)
            cmd.set('stick_radius', 0.10, pocketSel)
            cmd.util.cnc(pocketSel)
            log('Showing pocket residues within {} A of the ligand.'.format(args.pocketDist))

    # Frame on protein+ligand only -- the full object's hidden solvent box would shrink it.
    focusSel = polymerSel if cmd.count_atoms(polymerSel) else obj
    if ligSel:
        focusSel = '({}) or ({})'.format(focusSel, ligSel)
    cmd.orient(focusSel)
    frameCamera(focusSel, surface=(args.style == 'surface'), spin=bool(args.spin))


def frameCamera(sel, surface=False, spin=False, margin=3.0):
    """Zoom to fit, widened for spin (else the model clips as it rotates)."""
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


def renderStates(obj, states, args, framesDir):
    """Render states, named by global index. Spin is absolute (view reset + turned
    by globalIdx * anglePerFrame) so it stays continuous across chunk boundaries."""
    width, height = RESOLUTIONS[args.resolution]
    nFrames = len(states)
    nTotal = args.nTotalFrames or nFrames
    anglePerFrame = (360.0 * args.nTurns / nTotal) if (args.spin and nTotal > 1) else 0.0
    baseView = cmd.get_view() if anglePerFrame else None

    if not args.ray:
        cmd.set('ray_shadows', 0)
        cmd.set('ambient_occlusion_mode', 0)

    log('Rendering {} of {} total frames at {}x{} (ray={}, spin={}).'
        .format(nFrames, nTotal, width, height, args.ray, bool(args.spin)))

    frameFiles = []
    for localIdx, state in enumerate(states):
        globalIdx = args.frameOffset + localIdx
        cmd.set('state', state)
        if anglePerFrame:
            cmd.set_view(baseView)
            cmd.turn(args.spinAxis, globalIdx * anglePerFrame)
        framePath = os.path.join(framesDir, 'frame_{:05d}.png'.format(globalIdx))
        cmd.ray(width, height)
        cmd.png(framePath, dpi=300 if args.ray else 150)
        frameFiles.append(framePath)
        if (localIdx + 1) % 10 == 0 or localIdx == nFrames - 1:
            log('  ...{}/{} frames'.format(localIdx + 1, nFrames))
    return frameFiles


def frameGlobalIdx(framePath):
    return int(os.path.splitext(os.path.basename(framePath))[0].split('_')[-1])


def stampTimeLabels(frameFiles, timesFile):
    """Overlay "t = X ns" on each frame, from timesFile ({globalFrameIdx: timeNs})."""
    from PIL import Image, ImageDraw

    with open(timesFile) as f:
        times = json.load(f)

    font = None
    try:
        from PIL import ImageFont
        font = ImageFont.load_default()
    except Exception:
        pass

    for framePath in frameFiles:
        timeNs = times.get(str(frameGlobalIdx(framePath)))
        if timeNs is None:
            continue
        img = Image.open(framePath).convert('RGB')
        draw = ImageDraw.Draw(img)
        draw.text((12, 12), 't = {:.2f} ns'.format(timeNs), fill=(255, 0, 0), font=font)
        img.save(framePath)


def flattenFrames(frameFiles, flatDir):
    """Symlink each chunk's frames into one sequentially-named dir (ffmpeg needs that)."""
    os.makedirs(flatDir, exist_ok=True)
    flatFiles = []
    for i, src in enumerate(frameFiles):
        dst = os.path.join(flatDir, 'seq_{:05d}.png'.format(i))
        if not os.path.exists(dst):
            try:
                os.symlink(os.path.abspath(src), dst)
            except OSError:
                shutil.copy(src, dst)
        flatFiles.append(dst)
    return flatFiles


def assembleVideo(frameFiles, outBase, fps, fmt, crf=20):
    """Mux PNG frames (sorted by global frame index) into mp4 (ffmpeg) or gif (Pillow)."""
    if not frameFiles:
        raise RuntimeError('No frames were rendered.')
    frameFiles = sorted(frameFiles, key=frameGlobalIdx)

    ffmpeg = shutil.which('ffmpeg')
    if fmt == 'mp4' and ffmpeg:
        outFile = outBase + '.mp4'
        flatDir = outBase + '_ffmpeg_seq'
        flattenFrames(frameFiles, flatDir)
        pattern = os.path.join(flatDir, 'seq_%05d.png')
        cmdLine = [
            ffmpeg, '-y', '-framerate', str(fps), '-i', pattern,
            '-c:v', 'libx264', '-preset', 'slow', '-crf', str(crf),
            '-pix_fmt', 'yuv420p', '-movflags', '+faststart',
            # H.264 requires even dimensions.
            '-vf', 'scale=trunc(iw/2)*2:trunc(ih/2)*2',
            outFile,
        ]
        log('Encoding mp4 with ffmpeg...')
        subprocess.check_call(cmdLine)
        shutil.rmtree(flatDir, ignore_errors=True)
        return outFile

    # Fallback (or explicit request): animated gif via Pillow.
    if fmt == 'mp4' and not ffmpeg:
        log('ffmpeg not found -> falling back to animated gif.')
    from PIL import Image
    outFile = outBase + '.gif'
    log('Encoding gif with Pillow...')
    # A generator, not a list: holding every decoded frame at once doesn't scale.
    frameIter = (Image.open(f).convert('RGB') for f in frameFiles)
    firstFrame = next(frameIter)
    duration = int(1000.0 / max(1, fps))
    firstFrame.save(outFile, save_all=True, append_images=frameIter,
                    duration=duration, loop=0, optimize=True)
    return outFile


def main():
    args = parseArgs()
    doRender = not args.encodeOnly
    doEncode = not args.renderOnly

    workdir = args.workdir or (os.path.dirname(os.path.abspath(args.trajectory))
                               if args.trajectory else '.')
    os.makedirs(workdir, exist_ok=True)
    framesDir = args.framesDir or os.path.join(workdir, '_md_video_frames')
    outBase = os.path.join(workdir, args.output)

    frameFiles = []
    if doRender:
        os.makedirs(framesDir, exist_ok=True)
        log('System    : {}'.format(args.inputStruct))
        log('Trajectory: {}'.format(args.trajectory))

        cmd.feedback('disable', 'all', 'everything')
        if args.threads:
            cmd.set('max_threads', args.threads)

        obj, nStates = loadSystem(args.inputStruct, args.trajectory)
        log('Loaded {} states.'.format(nStates))

        if args.alignTraj:
            log('Aligning trajectory on the protein backbone (intra_fit).')
            cmd.intra_fit('{} and polymer and name CA'.format(obj))

        if args.smooth > 0 and nStates > 2:
            log('Smoothing trajectory (window={}, visualization only).'.format(args.smooth))
            cmd.smooth('all', passes=1, window=args.smooth)

        applyStyle(obj, args)
        frameFiles = renderStates(obj, resolveStates(args, nStates), args, framesDir)

    if doEncode:
        if not frameFiles:
            pattern = args.framesGlob or os.path.join(framesDir, '**', 'frame_*.png')
            frameFiles = sorted(glob.glob(pattern, recursive=True), key=frameGlobalIdx)

        if args.timeLabel and args.timesFile:
            stampTimeLabels(frameFiles, args.timesFile)

        outFile = assembleVideo(frameFiles, outBase, args.fps, args.format, args.crf)

        if not args.keepFrames:
            for frameDir in {os.path.dirname(f) for f in frameFiles}:
                shutil.rmtree(frameDir, ignore_errors=True)

        log('DONE. Video written to: {}'.format(outFile))


# PyMOL execs scripts with __name__ == 'pymol' (not '__main__'), so accept both.
if __name__ in ('__main__', 'pymol'):
    main()

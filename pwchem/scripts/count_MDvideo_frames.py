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

"""Count trajectory frames and per-frame times (mdtraj); writes {"nTotal", "timesPs"} JSON."""

import argparse
import json
import os

import mdtraj as md


def parseArgs():
    p = argparse.ArgumentParser(description='Count MD trajectory frames and their times.')
    p.add_argument('-i', '--inputStruct', required=True)
    p.add_argument('-t', '--trajectory', required=True)
    p.add_argument('-o', '--output', required=True)
    p.add_argument('-b', '--baseDir', required=True,
                   help='Every path must resolve inside this directory.')
    return p.parse_args()


def safePath(path, baseDir):
    """Canonicalise a CLI-derived path and refuse anything outside baseDir (S8707)."""
    resolved = os.path.realpath(path)
    base = os.path.realpath(baseDir)
    if resolved != base and not resolved.startswith(base + os.sep):
        raise ValueError('"{}" resolves outside the allowed directory "{}".'.format(path, base))
    return resolved


def main():
    args = parseArgs()
    inputStruct = safePath(args.inputStruct, args.baseDir)
    trajectory = safePath(args.trajectory, args.baseDir)
    output = safePath(args.output, args.baseDir)

    nTotal = 0
    timesPs = []
    for chunk in md.iterload(trajectory, top=inputStruct, chunk=500):
        nTotal += chunk.n_frames
        timesPs.extend(float(t) for t in chunk.time)

    with open(output, 'w') as f:
        json.dump({'nTotal': nTotal, 'timesPs': timesPs}, f)


if __name__ == '__main__':
    main()

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

import mdtraj as md


def parseArgs():
    p = argparse.ArgumentParser(description='Count MD trajectory frames and their times.')
    p.add_argument('-i', '--inputStruct', required=True)
    p.add_argument('-t', '--trajectory', required=True)
    p.add_argument('-o', '--output', required=True)
    return p.parse_args()


def main():
    args = parseArgs()
    nTotal = 0
    timesPs = []
    for chunk in md.iterload(args.trajectory, top=args.inputStruct, chunk=500):
        nTotal += chunk.n_frames
        timesPs.extend(float(t) for t in chunk.time)

    with open(args.output, 'w') as f:
        json.dump({'nTotal': nTotal, 'timesPs': timesPs}, f)


if __name__ == '__main__':
    main()

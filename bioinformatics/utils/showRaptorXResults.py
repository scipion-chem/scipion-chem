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

from multiprocessing import Process
import os
import sys

SERVERDIR = 'raptorx.uchicago.edu'

def serverRun(fnDir):
    import socket
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        if s.connect_ex(('localhost', 1234)) == 0:
            data = [(int(p), c) for p, c in [x.rstrip('\n').split(' ', 1) for x in os.popen('ps h -eo pid:1,command')]]
            for pid, cmd in data:
                if "python -m http.server 1234" in cmd:
                    os.system('kill -9 %d'%pid)
    print("Starting at",fnDir)
    os.system("cd %s; python -m http.server 1234"%fnDir)

if __name__ == "__main__":
    if len(sys.argv)==1:
        print("Usage: python3 showRaptorXResults.py [url]")
    url = sys.argv[1]
    fnDir=url.split(SERVERDIR)[0]
    localurl=url.split(SERVERDIR)[1]

    serverDir = os.path.join(fnDir,SERVERDIR)
    server = Process(target=serverRun, args=(serverDir,))
    server.start()
    os.system("sensible-browser http://localhost:1234/%s"%localurl)
    server.terminate()

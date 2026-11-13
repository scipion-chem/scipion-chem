import os
from subprocess import Popen

import pyworkflow.utils as pwutils
from pwem.wizards import VariableWizard

from pwchem import Plugin as pwchemPlugin
from pwchem.protocols.VirtualDrugScreening.protocol_filter_contacts import (
    ProtContactLigands
)
import json
from pwchem.constants import OPENBABEL_DIC, RDKIT_DIC

from pwchem.utils import pdbqt2other
from pwem.convert import cifToPdb

from Bio.PDB import PDBParser, MMCIFParser

class PyMol:
    """Help class to run PyMOL and manage its environment."""

    @classmethod
    def getEnviron(cls):
        environ = pwutils.Environ(os.environ)
        environ.set(
            'PATH',
            pwchemPlugin.getProgramHome(OPENBABEL_DIC, path='bin'),
            position=pwutils.Environ.BEGIN
        )
        return environ


class PreviewLigandAtomIDsWizard(VariableWizard):
    _targets, _inputs, _outputs = [], {}, {}

    def show(self, form, *params):
        import os
        import shlex
        from subprocess import Popen

        protocol = form.protocol
        ligandSet = protocol.inputLigandSet.get()

        mol_obj = next(iter(ligandSet), None)

        if mol_obj is None:
            print("ERROR: The input ligand set is empty.")
            return

        ligandFile = mol_obj.getPoseFile() or mol_obj.getFileName()

        if not ligandFile or not os.path.isfile(ligandFile):
            print("ERROR: The ligand pose file was not found.")
            return

        ligandFile = os.path.abspath(ligandFile)
        extension = os.path.splitext(ligandFile)[1].lower()
        project = protocol.getProject()

        # Convert unsupported formats to SDF using the Scipion-chem
        # Open Babel conversion script.
        if extension not in (".sdf", ".mol"):
            convertedFile = os.path.abspath(
                project.getTmpPath(
                    "ligand_atom_numbers_{}.sdf".format(
                        protocol.getObjId()
                    )
                )
            )

            os.makedirs(os.path.dirname(convertedFile), exist_ok=True)

            try:
                args = '-i "{}" -of sdf -o "{}" -od "{}" --overWrite'.format(
                    os.path.abspath(ligandFile),
                    os.path.splitext(os.path.basename(convertedFile))[0],
                    os.path.dirname(convertedFile)
                )

                pwchemPlugin.runScript(
                    protocol,
                    "obabel_IO.py",
                    args,
                    env=OPENBABEL_DIC,
                    popen=True,
                    cwd=os.path.dirname(ligandFile)
                )

                # obabel_IO.py may name the output from the input basename.
                generatedFile = os.path.join(
                    os.path.dirname(convertedFile),
                    os.path.splitext(os.path.basename(ligandFile))[0]
                    + ".sdf"
                )

                if os.path.isfile(generatedFile):
                    if generatedFile != convertedFile:
                        os.replace(generatedFile, convertedFile)

                if (
                    not os.path.isfile(convertedFile)
                    or os.path.getsize(convertedFile) == 0
                ):
                    print(
                        "ERROR: Open Babel did not generate the expected "
                        "SDF file: {}".format(convertedFile)
                    )
                    return

                ligandFile = convertedFile

            except Exception as exc:
                import traceback
                traceback.print_exc()
                print("ERROR: Could not convert ligand to SDF:", exc)
                return

        imageFile = os.path.abspath(
            project.getTmpPath(
                "ligand_atom_numbers_{}.png".format(
                    protocol.getObjId()
                )
            )
        )

        args = '{} "{}"'.format(
            shlex.quote(ligandFile),
            imageFile
        )

        try:
            pwchemPlugin.runScript(
                protocol,
                "rdkit_ligand_atom_numbers.py",
                args,
                env=RDKIT_DIC,
                popen=True,
                cwd=os.path.dirname(ligandFile)
            )
        except Exception as exc:
            import traceback
            traceback.print_exc()
            print(
                "ERROR: Could not generate the ligand atom-number figure:",
                exc
            )
            return

        if not os.path.isfile(imageFile):
            print(
                "ERROR: The RDKit script did not generate the expected "
                "image: {}".format(imageFile)
            )
            return

        print("Generated atom-number figure: {}".format(imageFile))

        try:
            Popen(["xdg-open", imageFile])
        except OSError as exc:
            print("ERROR opening ligand figure: {}".format(exc))
            print("Open manually: {}".format(imageFile))

PreviewLigandAtomIDsWizard().addTarget(
    protocol=ProtContactLigands,
    targets=['ligAtom'],
    inputs=['inputLigandSet'],
    outputs=[]
)

class PreviewReceptorAtomIDsWizard(VariableWizard):
    _targets, _inputs, _outputs = [], {}, {}

    def show(self, form, *params):
        protocol = form.protocol

        ligandSet = protocol.inputLigandSet.get()
        receptorFile = ligandSet.getProteinFile()

        if not receptorFile or not os.path.isfile(receptorFile):
            print("ERROR: Receptor file not found.")
            return

        chainData = json.loads(protocol.recChain.get())
        chainId = str(chainData["chain"])

        residueData = json.loads(protocol.recResidue.get())
        residueId = str(residueData["index"]).split("-")[0]

        project = protocol.getProject()

        # Convert PDBQT/CIF to PDB for PyMOL visualization.
        extension = os.path.splitext(receptorFile)[1].lower()

        if extension in (".pdbqt", ".cif", ".mmcif"):
            receptorPdb = os.path.abspath(
                project.getTmpPath(
                    "{}_preview.pdb".format(
                        os.path.splitext(os.path.basename(receptorFile))[0]
                    )
                )
            )

            os.makedirs(os.path.dirname(receptorPdb), exist_ok=True)

            if extension == ".pdbqt":
                pdbqt2other(protocol, receptorFile, receptorPdb)
            else:
                cifToPdb(receptorFile, receptorPdb)

        else:
            receptorPdb = os.path.abspath(receptorFile)

        scriptFile = os.path.abspath(
            project.getTmpPath(
                "previewReceptorAtomIDs_{}.pml".format(protocol.getObjId())
            )
        )

        labelScript = os.path.abspath(
            project.getTmpPath(
                "labelReceptorAtoms_{}.py".format(protocol.getObjId())
            )
        )

        with open(labelScript, "w") as f:
            f.write("from pymol import cmd\n")
            f.write(
                'selection = "receptor and chain {} and resi {} and not elem H and not elem D"\n'
                .format(chainId, residueId)
            )
            f.write("atoms = cmd.get_model(selection).atom\n")
            f.write("for i, atom in enumerate(atoms, 1):\n")
            f.write('    label = "{}:{}".format(atom.resn, atom.name)\n')
            f.write(
                '    cmd.label("receptor and index {}".format(atom.index), '
                '"\\"{}\\"".format(label))\n'
            )

        pymolReceptor = receptorPdb.replace("\\", "/")

        with open(scriptFile, "w") as f:
            f.write('load "{}", receptor\n'.format(pymolReceptor))
            f.write("hide everything\n")
            f.write("show cartoon, receptor\n")
            f.write(
                "show sticks, receptor and chain {} and resi {}\n"
                .format(chainId, residueId)
            )
            f.write(
                "color red, receptor and chain {} and resi {}\n"
                .format(chainId, residueId)
            )
            f.write(
                "zoom receptor and chain {} and resi {}\n"
                .format(chainId, residueId)
            )
            f.write("run {}\n".format(labelScript.replace("\\", "/")))
            f.write("set label_size, 18\n")
            f.write("set label_color, white\n")
            f.write(
                "zoom receptor and chain {} and resi {}\n"
                .format(chainId, residueId)
            )

        try:
            pymolBin = pwchemPlugin.getEnvPath(
                OPENBABEL_DIC, "bin/pymol"
            )
            Popen(
                [pymolBin, scriptFile],
                cwd=os.path.dirname(receptorPdb),
                env=PyMol.getEnviron()
            )
        except OSError as e:
            print("ERROR: Could not launch PyMOL: {}".format(e))


PreviewReceptorAtomIDsWizard().addTarget(
    protocol=ProtContactLigands,
    targets=["recAtom"],
    inputs=["inputLigandSet"],
    outputs=[]
)
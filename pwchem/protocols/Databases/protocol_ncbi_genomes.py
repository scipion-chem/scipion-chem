# **************************************************************************
# *
# * Authors: Noelia Auba
# *          
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307 USA
# *
# * All comments concerning this program package may be sent to the
# * e-mail address 'you@yourinstitution.email'
# *
# **************************************************************************
# **************************************************************************
# *
# * Authors: Noelia Auba
# *
# **************************************************************************

# General imports
import time
from pathlib import Path

import requests
from Bio import Entrez

# Scipion imports
from pwem.protocols import EMProtocol
from pyworkflow.protocol import params

# Plugin imports
from pwchem.objects import Genome, SetOfGenome

class ProtocolGenomeFetching(EMProtocol):
    """
    Download bacterial genome protein FASTA files from NCBI using genome accession numbers.
    The output is a SetOfGenome object containing the downloaded FASTA files.
    """

    _label = 'Genome fetching'

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """Define protocol form."""
        form.addSection(label='Input')
        group = form.addGroup('NCBI Parameters')
        group.addParam(
        'email',
        params.StringParam,
        label='Email',
        help='Email required by NCBI Entrez API policy.'
    )
        group.addParam(
            'accessions',
            params.TextParam,
            label='Accession Numbers:',
            help='Enter genome accession numbers separated by commas (NC_, NZ_, GCF_, GCA_)' \
            ' https://www.ncbi.nlm.nih.gov/datasets/genome/'
        )

        form.addParallelSection(threads=1, mpi=1)

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep(self.downloadStep)
        self._insertFunctionStep(self.createOutputStep)

    def downloadStep(self):
        Entrez.email = self.email.get()
        Entrez.tool = "GenomeDownloader"

        accessions = [a.strip() for a in self.accessions.get().split(",") if a.strip()]
        self.results = []

        outDir = Path(self._getExtraPath("genomes"))
        outDir.mkdir(exist_ok=True)

        for accession in accessions:
            self.info(f"Searching genome: {accession}")

            uid = self.searchAssembly(accession)
            if not uid:
                self.warning(f"No assembly found for {accession}")
                continue

            fasta_url, organism = self.getFasta(uid)
            if not fasta_url:
                self.warning(f"No FASTA found for {accession}")
                continue

            outFile = outDir / f"{accession}_protein.faa.gz"

            try:
                self.downloadFile(fasta_url, outFile)
                self.results.append((accession, organism, outFile))
                self.info(f"Downloaded: {accession}")

            except Exception as e:
                self.warning(f"Download failed for {accession}: {e}")

            time.sleep(2)
            resultsFile = Path(self._getExtraPath("results.tsv"))

            with open(resultsFile, "w") as f:
                for acc, org, path in self.results:
                    f.write(f"{acc}\t{org}\t{path}\n")

    def createOutputStep(self):
        genomeSet = SetOfGenome().create(
            outputPath=self._getPath(),
            suffix='genomes'
        )
        
        resultsFile = Path(self._getExtraPath("results.tsv"))

        if not resultsFile.exists():
            raise RuntimeError("Results file not found")

        with open(resultsFile) as f:
            for line in f:
                accession, organism, outFile = line.strip().split("\t")
        #for accession, organism, outFile in getattr(self, "results", []):
                genome = Genome()
                genome.setAccession(accession)
                genome.setOrganism(organism)
                genome.setFastaFile(str(outFile))
                genomeSet.append(genome)

        self._defineOutputs(outputGenomes=genomeSet)

    # --------------------------- Summary functions ------------------------------
    def _summary(self):
        if hasattr(self, "outputGenomes"):
            n = self.outputGenomes.getSize()
            return [f"{n} genome protein FASTA files were downloaded from NCBI."]
            
        return ["No genomes downloaded"]
   
    
    def _methods(self):
        if hasattr(self, "outputGenomes") and self.outputGenomes is not None:
            n = self.outputGenomes.getSize()
            return [f"{n} genome protein FASTA files were downloaded from NCBI."]
        return ["Genome fetching protocol executed."]

    def _validate(self):
        errors = []

        if not self.accessions.get().strip():
            errors.append("You must provide at least one genome accession.")

        return errors

    # --------------------------- UTILS functions ------------------------------
    def searchAssembly(self, term):
        handle = Entrez.esearch(
            db="assembly",
            term=term,
            retmax=1
        )

        record = Entrez.read(handle)
        handle.close()

        if not record["IdList"]:
            return None

        return record["IdList"][0]

    def getFasta(self, uid):
        handle = Entrez.esummary(
            db="assembly",
            id=uid,
            report="full"
        )

        summary = Entrez.read(handle)
        handle.close()

        docsum = summary["DocumentSummarySet"]["DocumentSummary"][0]
        ftp_path = docsum.get("FtpPath_RefSeq")
        organism = docsum.get("Organism")
        print(ftp_path)

        if not ftp_path:
            return None, None

        base = ftp_path.split("/")[-1]
        fasta_url = f"{ftp_path}/{base}_protein.faa.gz"
        fasta_url = fasta_url.replace("ftp://", "https://")

        return fasta_url, organism

    def downloadFile(self, url, out_path):
        r = requests.get(url, stream=True)
        r.raise_for_status()

        with open(out_path, "wb") as f:
            for chunk in r.iter_content(chunk_size=8192):
                if chunk:
                    f.write(chunk)
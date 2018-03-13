from Bio.PDB import PDBParser, PDBIO, PPBuilder
import sys

file = sys.argv[1]

io = PDBIO()
pdb = PDBParser().get_structure(file)
polipeptide = PPBuilder()

for chain in pdb.get_chains():
	io.set_structure(chain)
	io.save(pdb.get_(id) + "_" + chain.get_id() + ".pdb")

for pp in polipeptide.build_peptides(pdb):
	fasta = open(pdb.get_(id) + "_" + chain.get_id() + ".fa", "w")
	fasta.write(">" + pdb.get_(id) + "_" + chain.get_id() + "\n")
	fasta.write(pp.get_sequence())
		


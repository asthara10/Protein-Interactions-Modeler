import argparse

def SplitChain(infiles):
	"""
	blah
	"""
	chain_names = set()
	for file in options.infiles:
		name = file.split(".")[0]
		io = PDBIO()
		pdb = PDBParser().get_structure(name, file)
		polipeptide = PPBuilder()
		for chain in pdb.get_chains():
			print (chain.get_id())
			if chain.get_id() not in chain_names:
				io.set_structure(chain)
				io.save(pdb.get_id() + "_" + chain.get_id() + ".pdb")

				for pp in polipeptide.build_peptides(pdb):
					fasta = open(pdb.get_id() + "_" + chain.get_id() + ".fa", "w")
					fasta.write(">" + pdb.get_id() + "_" + chain.get_id() + "\n")
					fasta.write(str(pp.get_sequence()))
			chain_names.add(chain.get_id())

if __name__ == "__main__":

	from Bio.PDB import PDBParser, PDBIO, PPBuilder

	parser = argparse.ArgumentParser(description="blah")

	parser.add_argument('-i', '--input',
				dest= "infiles",
				action= "store",
				required = True,
				help="blah",
				nargs = "+")

	options=parser.parse_args()

	SplitChain(options.infiles)
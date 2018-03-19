import argparse

def Parse_PDB(pdbfiles):
	"""
	Parses PDB files using biopython by creating PDBParser objects.
	"""
	PDB_objects = []
	PDB_names = []
	for file in pdbfiles:
		name = file.split(".")[0]
		pdb = PDBParser().get_structure(name, file)
		PDB_objects.append(pdb)
		PDB_names.append(name)
	return (PDB_objects, PDB_names)

def SplitChain(PDB_objects):
	"""
	Splits a list of PDB files by chain creating one PDB and one FASTA file per chain.
	"""
	chain_names = set()
	File_prefix = []

	for pdb in PDB_objects:
		io = PDBIO()

		# Creates a PDB file for each chain of the original file.
		for chain in pdb.get_chains():
			if chain.get_id() not in chain_names:
				io.set_structure(chain)
				io.save(pdb.get_id() + "_" + chain.get_id() + ".pdb")
				File_prefix.append(pdb.get_id() + "_" + chain.get_id())

				# Creates a FASTA file for each chain of the original file.
				polipeptide = PPBuilder()
				for pp in polipeptide.build_peptides(pdb):
					fasta = open(pdb.get_id() + "_" + chain.get_id() + ".fa", "w")
					fasta.write(">" + pdb.get_id() + "_" + chain.get_id() + "\n")
					fasta.write(str(pp.get_sequence()))

				chain_names.add(chain.get_id())

	return File_prefix

def Run_BLAST(db_path, prefix):
	"""
	Runs a psiBLAST application.
	"""
	fasta = prefix + ".fa"
	xml = prefix + ".xml"
	pssm = prefix + "_pssm"

	Psi_query = Ncbicmd('psiblast', db = db_path, query = fasta, out = xml, evalue = 10, outfmt = 3, out_pssm = pssm)
	Psi_query() # Run BLAST.

	return xml

def Select_template(BLAST_outs):
	"""
	Selects the best templates from the BLAST output for all chains.
	"""
	Outputs = {}

	for Out in BLAST_outs:
		Blast_out = open(Out)
		First = True
		for line in Blast_out:
			if line.startswith("Sequences producing significant alignments:"):
				Blast_out.readline()
				for line in Blast_out:
					line = line.strip()
					line = line.split()
					if First:
						Min_evalue = line[len(line)-1]
						Template_evalue = [(line[0][:-2], line[len(line)-1])]
						Outputs[Out] = Template_evalue
						First = False
					else:
						if line[len(line)-1] == Min_evalue:
							Template_evalue.append((line[0][:-2], line[len(line)-1]))
							Outputs[Out] = Template_evalue
						else:
							break
	# Select all possible best templates, the ones with the minimum evalue.
	min_evalue = min(map(lambda x: min(map(lambda y: y[1], x)), Outputs.values()))
	templates = set()
	for value in Outputs.values():
		for template in value:
			if template[1] == min_evalue:
				templates.add(template[1])

	return templates

if __name__ == "__main__":

	from Bio.PDB import PDBParser, PDBIO, PPBuilder
	from Bio.PDB import PDBList
	from Bio.Blast.Applications import NcbipsiblastCommandline as Ncbicmd
	from Bio.Align.Applications import ClustalwCommandline
	import Blast_align

	parser = argparse.ArgumentParser(description="blah")

	parser.add_argument('-i', '--input',
				dest= "infiles",
				action= "store",
				required = True,
				nargs = "+",
				help="Input file names (PDB files with pairs of chains that interact), insert as many as needed.")
	parser.add_argument('-d', '--database',
				dest= "database",
				action= "store",
				required = True,
				help="Path to the pdb database in your computer.")
	

	options=parser.parse_args()

	(PDB_objects, PDB_names) = Parse_PDB(options.infiles)
	PDB_prefix = SplitChain(PDB_objects)
	BLAST_outs = []
	for prefix in PDB_prefix:
		output = Run_BLAST(options.database, prefix)
		BLAST_outs.append(output)
	Templates = Select_template(BLAST_outs)
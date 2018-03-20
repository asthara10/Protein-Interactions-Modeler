import argparse

def Parse_PDB(pdbfiles):
	"""
	Parses PDB files using biopython by creating PDBParser objects.
	"""
	PDB_objects = []
	PDB_names = []
	for file in pdbfiles:
		name = file.split(".")[0]
		pdb = PDBParser(QUIET=True).get_structure(name, file)
		PDB_objects.append(pdb)
		PDB_names.append(name)
	return (PDB_objects, PDB_names)

def SplitChain(PDB_objects):
	"""
	Splits a list of PDB files by chain creating one PDB and one FASTA file per chain.
	"""
	File_prefix = []

	for pdb in PDB_objects:
		chain_names = set()
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
		BLAST_out = open(Out)
		First = True
		for line in BLAST_out:
			if line.startswith("Sequences producing significant alignments:"):
				BLAST_out.readline()
				for line in BLAST_out:
					line = line.strip()
					line = line.split()
					if First:
						min_evalue = line[len(line)-1]
						template_evalue = [(line[0][:-2], line[len(line)-1])]
						Outputs[Out] = template_evalue
						First = False
					else:
						if line[len(line)-1] == min_evalue:
							template_evalue.append((line[0][:-2], line[len(line)-1]))
							Outputs[Out] = template_evalue
						else:
							break
	# Select all possible best templates, the ones with the minimum evalue.
	min_evalue = min(map(lambda x: min(map(lambda y: y[1], x)), Outputs.values()))
	templates = set()
	for value in Outputs.values():
		for template in value:
			if template[1] == min_evalue:
				templates.add(template[0])

	return templates

def Download_template(template):
	pdbl = PDBList()
	pdbl.retrieve_pdb_file(template, obsolete=False, pdir="./", file_format="pdb")

def Create_joined_fastas(input_PDB_objects):
	polipeptide = PPBuilder()
	first_line = True
	filename = ""

	for obj in input_PDB_objects:
		filename = filename + obj.get_id()
	filename = filename + ".fa"
	joined_fasta = open(filename, 'w')

	for obj in input_PDB_objects:
		if first_line:
			joined_fasta.write(">" + obj.get_id() + "\n")
			first_line = False
		else:
			joined_fasta.write("\n" + ">" + obj.get_id() + "\n")
		for polipep in polipeptide.build_peptides(obj):
			joined_fasta.write(str(polipep.get_sequence()))

	return filename

def Create_fasta(input_pdb):
	polipeptide = PPBuilder()
	fist_chain = True
	pdb = PDBParser().get_structure(input_pdb, "pdb" + input_pdb + ".ent")
	for pp in polipeptide.build_peptides(pdb):
		if fist_chain:
			new_fasta = open(input_pdb + "_paired.fa", "w")
			new_fasta.write(">" + input_pdb +"\n")
			new_fasta.write(str(pp.get_sequence()))
			fist_chain = False
		else:
			new_fasta.write(str(pp.get_sequence()))

def Run_clustal(fastas):
	for fasta in fastas:
		clustalw_cline = ClustalwCommandline("clustalw2", infile=fasta)
		stdout, stderr = clustalw_cline()
		with open(fasta.split(".")[0] + "_ClustalScore.txt", 'w') as scores:
			scores.write(stdout)

def Analize_clustal_score(sc_file):
	file = open(sc_file, 'r')
	for line in file:
		if line.startswith("Group "):
			line = line.strip()
			line = line.split()

if __name__ == "__main__":

	from Bio.PDB import PDBParser, PDBIO, PPBuilder
	from Bio.PDB import PDBList
	from Bio.Blast.Applications import NcbipsiblastCommandline as Ncbicmd
	from Bio.Align.Applications import ClustalwCommandline
	import copy
	#import Blast_align

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

	# Parse input files
	(PDB_input_objects, PDB_input_names) = Parse_PDB(options.infiles)
	file_prefixes = SplitChain(PDB_input_objects)
	bychain_PDBs = map(lambda x: x + ".pdb", file_prefixes)
	(PDB_bychain_objects, PDB_bychain_names) = Parse_PDB(bychain_PDBs)

	# Look for templates
	BLAST_outs = []
	for prefix in file_prefixes:
		output = Run_BLAST(options.database, prefix)
		BLAST_outs.append(output)
	Templates = Select_template(BLAST_outs)

	# Work with templates
	fasta_names = []

	for template in Templates:
		Download_template(template)
	temp_PDBs = map(lambda x: "pdb" + x + ".ent", Templates)
	print("template pdbs")
	print(temp_PDBs)
	(PDB_temp_objs, PDB_temp_names) = Parse_PDB(temp_PDBs)	
	print("template objects")
	print(PDB_temp_objs)
	template_chains = SplitChain(PDB_temp_objs)
	print("template chains prefix")
	print(template_chains)
	bychain_PDBs = map(lambda x: x + ".pdb", template_chains)
	print("template pdb by chain files")
	print(bychain_PDBs)
	(PDB_chain_temp_objs, PDB_chain_temp_names) = Parse_PDB(bychain_PDBs)
	for chain in PDB_chain_temp_objs:
		obj_list = copy.copy(PDB_bychain_objects)
		obj_list.append(chain)
		joined_file = Create_joined_fastas(obj_list)
		fasta_names.append(joined_file)
	Run_clustal(fasta_names)


	# Check if all input chains are equal
	"""joined_file = Create_joined_fastas(PDB_bychain_objects)
	Run_clustal(joined_file)
	clu_score = joined_file.split(".")[0] + "_ClustalScore.txt"
	Analize_clustal_score(clu_score)"""

	

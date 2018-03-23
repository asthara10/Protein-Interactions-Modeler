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
		filename = filename + obj.get_id() + "_"
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

"""def Create_fasta(input_pdb):
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
			new_fasta.write(str(pp.get_sequence()))"""

def Run_clustal(fastas):
	for fasta in fastas:
		clustalw_cline = ClustalwCommandline("clustalw2", infile=fasta)
		stdout, stderr = clustalw_cline()
		with open(fasta.split(".")[0] + "ClustalScore.txt", 'w') as scores:
			scores.write(stdout)

def Analize_clustal_score(sc_file, temp_name):
	file = open(sc_file, 'r')
	equivalences = {}
	aligns = []
	equiv_reg = re.compile('(Sequence )([0-9]+)(: )(\S+)')
	alig_reg = re.compile('(Sequences \()([0-9]+)(:)([0-9]+)(\) Aligned. Score:  )([0-9]+)')
	for line in file:
		if re.match(equiv_reg, line):
			match = re.match(equiv_reg, line)
			equivalences[match[2]] = match[4]
			if match[4] == temp_name:
				template = match[2]
		elif re.match(alig_reg, line):
			alig = re.match(alig_reg, line)
			if (template == alig[2]) and (alig[6] == "100"):
				aligns.append(equivalences[alig[4]]) 
			elif (template == alig[4]) and (alig[6] == "100"):
				aligns.append(equivalences[alig[2]])

	return aligns

def Find_interactions(PDB_obj):
	interact_chains = []

	chains = Selection.unfold_entities(PDB_obj, 'C')

	obj_atoms = Selection.unfold_entities(PDB_obj, 'A')
	neighbors = NeighborSearch(obj_atoms)
	for chain in chains:
		atoms = Selection.unfold_entities(chain, 'A')
		for center in atoms:
			interactions = neighbors.search(center.coord,2.5,level='C')
			ids = list(map(lambda x: x.get_id(), interactions))
			if len(ids) > 1:
				final_ids = list(filter(lambda x: x != chain.get_id(), ids))
				interact_chains.append((chain.get_id(), final_ids))
	
	return interact_chains

def Assign_query_to_temp(i, cand_list, temp_chains, Final_interactions, temp):
	j = 0
	targ = cand_list[i][0]
	
	if len(cand_list) < i :
		return True
	else:
		while j < len(cand_list[i][1]):
			#hem de comprovar que el cand_list[i][1][j] no estigui agafat per cap objectiu anterior
			if temp_chains[cand_list[i][1][j]] == None:
				temp_chains[cand_list[i][1][j]] = targ
				#si la relacio existeix en el objectiu
				for prev_targ in Final_interactions["target_interacts"].keys():
					if targ in Final_interactions["target_interacts"][prev_targ]:
						for key, val in temp_chains.items():
							if val == prev_targ:
								prev_temp = key
						#si no es compleix que la relacio existeix en el template
						if temp_chains[cand_list[i][1][j]] not in Final_interactions["temps"][temp]["temp_interact"][key]:
							#posar a none el valor del diccionari temp_chains
							temp_chains[cand_list[i][1][j]] = None
							j += 1
							return False
			else:
				j += 1
				return False
			if Assign_query_to_temp(i+1, cand_list, temp_chains, Final_interactions, temp):
				return True
			j += 1
		#posar a none el valor del diccionari temp_chains
		#temp_chains[cand_list[i][1][j]] = None
		return False

def I_Assign_query_to_temp(targ_chain_list, temp_chains, Final_interactions, temp):
	candidates = []
	temporal_cand = []
	#crear llista de  tuple (string,llistes candidats)
	for target in targ_chain_list:
		for template in Final_interactions["temps"][temp]["target_temp"].keys():
			if target in Final_interactions["temps"][temp]["target_temp"][template]:
				temporal_cand.append(template)
		candidates.append((target, temporal_cand))
		if temporal_cand != []:
			temporal_cand = []
		else:
			return

	#fer primera crida recursiva
	Assign_query_to_temp(0, candidates, temp_chains, Final_interactions, temp)

def Superimpose_chains(temp_obj, PDB_bychain_objects, temp_chains):
	i = 0
	ref_model = temp_obj[0]
	ppbuild = PPBuilder()
	template_chains = Selection.unfold_entities(temp_obj, 'C')
	min_len1 = min(list(map(lambda x: len(ppbuild.build_peptides(x)[0].get_sequence()), template_chains)))
	min_len2 = min(list(map(lambda x: len(ppbuild.build_peptides(x)[0].get_sequence()), PDB_bychain_objects)))
	min_len = min([min_len1, min_len2])
	atoms_to_be_aligned = range(2, min_len)

	for sample_structure in PDB_bychain_objects:
		sample_model = sample_structure[0]
		ref_atoms = []
		sample_atoms = []

		for ref_chain in ref_model:
			for key, val in temp_chains.items():
				if val == sample_structure.get_id():
					if key[:-2] == temp_obj.get_id():
						temp_ch = key
			if temp_obj.get_id() + "_" + ref_chain.get_id() == temp_ch:
				for ref_res in ref_chain:
					if ref_res.get_id()[1] in atoms_to_be_aligned:
						ref_atoms.append(ref_res['CA'])

		for sample_chain in sample_model:
			for sample_res in sample_chain:
				if sample_res.get_id()[1] in atoms_to_be_aligned:
					sample_atoms.append(sample_res['CA'])

		super_imposer = Superimposer()
		super_imposer.set_atoms(ref_atoms, sample_atoms)
		super_imposer.apply(sample_atoms)

		#print(super_imposer.rms)
		# possible millora si tenim temps, filtrar per rmsd

		io = PDBIO()
		io.set_structure(sample_structure)
		io.save(temp_obj.get_id() + "_" + str(i) + "_aligned.pdb", write_end = False)
		i += 1

	j = copy.copy(i)
	i = 1
	file = open(temp_obj.get_id() + "_0_aligned.pdb", 'a')
	while i < j:
		file2 = open(temp_obj.get_id() + "_" + str(i) + "_aligned.pdb")
		for line in file2:
			file.write(line)
		i += 1



if __name__ == "__main__":

	from Bio.PDB import PDBParser, PDBIO, PPBuilder, Superimposer
	from Bio.PDB import PDBList
	from Bio.PDB import NeighborSearch, Selection
	from Bio.Blast.Applications import NcbipsiblastCommandline as Ncbicmd
	from Bio.Align.Applications import ClustalwCommandline
	import copy
	import numpy
	import re
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

	Final_interactions = {}
	Final_interactions["target_interacts"] = {}
	Final_interactions["temps"] = {}

	# Parse input files
	(PDB_input_objects, PDB_input_names) = Parse_PDB(options.infiles)
	file_prefixes = SplitChain(PDB_input_objects)

	first = True
	already_added = []
	unique_chains = []
	for pref in file_prefixes:
		if first:
			unique_chains.append(pref)
			already_added.append(pref.split("_")[1])
			first = False
		if pref.split("_")[1] not in already_added:
			unique_chains.append(pref)
			already_added.append(pref.split("_")[1])

	bychain_PDBs = map(lambda x: x + ".pdb", unique_chains)
	(PDB_bychain_objects, PDB_bychain_names) = Parse_PDB(bychain_PDBs)
	for inp in PDB_input_objects: # Add data to Final_interactions dictionary
		inp_chains = inp.get_chains()
		inp_chains_ids = list(map(lambda x: x.get_id(), inp_chains))
		if inp_chains_ids[0] not in Final_interactions["target_interacts"].keys():
			Final_interactions["target_interacts"][inp_chains_ids[0]] = list(inp_chains_ids[1])
		elif (inp_chains_ids[0] in Final_interactions["target_interacts"].keys()) and (inp_chains_ids[1] not in Final_interactions["target_interacts"][inp_chains_ids[0]]):
			Final_interactions["target_interacts"][inp_chains_ids[0]].append(inp_chains_ids[1])

	# Look for templates
	BLAST_outs = []
	for prefix in file_prefixes:
		output = Run_BLAST(options.database, prefix)
		BLAST_outs.append(output)
	Templates = Select_template(BLAST_outs)

	# Work with templates
	fasta_names = []
	temp_chains = {}

	for template in Templates:
		Download_template(template)
	temp_PDBs = map(lambda x: "pdb" + x + ".ent", Templates)
	(PDB_temp_objs, PDB_temp_names) = Parse_PDB(temp_PDBs)	
	template_chains = SplitChain(PDB_temp_objs)
	bychain_PDBs = map(lambda x: x + ".pdb", template_chains)
	(PDB_chain_temp_objs, PDB_chain_temp_names) = Parse_PDB(bychain_PDBs)
	for chain in PDB_chain_temp_objs:
		obj_list = copy.copy(PDB_bychain_objects)
		obj_list.append(chain)
		joined_file = Create_joined_fastas(obj_list)
		fasta_names.append(joined_file)
	Run_clustal(fasta_names)
	first = True
	for fa_name in fasta_names:
		temp_name = fa_name.split("_")[-3] + "_" + fa_name.split("_")[-2]
		if first:
			tmp = temp_name[:-2]
			Final_interactions["temps"][temp_name[:-2]] = {}
			Final_interactions["temps"][temp_name[:-2]]["target_temp"] = {}
			first = False
		#print(temp_name)
		temp_chains[temp_name] = None
		aligns = Analize_clustal_score(fa_name.split(".")[0] + "ClustalScore.txt", temp_name)
		if tmp == temp_name[:-2]:
			Final_interactions["temps"][temp_name[:-2]]["target_temp"][temp_name] = aligns
		else:
			Final_interactions["temps"][temp_name[:-2]] = {}
			Final_interactions["temps"][temp_name[:-2]]["target_temp"] = {}
			Final_interactions["temps"][temp_name[:-2]]["target_temp"][temp_name] = aligns
			tmp = temp_name[:-2]
			
		Final_interactions["temps"][temp_name[:-2]]["temp_interact"] = {}
		temp_obj = list(filter(lambda x: x.get_id() == temp_name[:-2], PDB_temp_objs))
		list_interacts = Find_interactions(temp_obj[0])
		for interact in list_interacts:
			if interact[0] in Final_interactions["temps"][temp_name[:-2]]["temp_interact"].keys():
				Final_interactions["temps"][temp_name[:-2]]["temp_interact"][interact[0]].update(interact[1])
			else:
				Final_interactions["temps"][temp_name[:-2]]["temp_interact"][interact[0]] = set(interact[1])
	#print(Final_interactions)


	#asignation_results = Assign_query_to_temp(PDB_bychain_names, temp_chains, Final_interactions, "pdb1a0u")
	#for temp_obj in PDB_temp_objs:
		#Superimpose_chains(temp_obj, PDB_bychain_objects)
	for temp in PDB_temp_names:
		I_Assign_query_to_temp(PDB_bychain_names, temp_chains, Final_interactions, temp)

	for temp_obj in PDB_temp_objs:
		for key, val in temp_chains.items():
			if (temp_obj.get_id() == key[:-2]) and (val != None):
				#print("hola")
				Superimpose_chains(temp_obj, PDB_bychain_objects, temp_chains)
	#print(temp_chains)





	# Check if all input chains are equal
	"""joined_file = Create_joined_fastas(PDB_bychain_objects)
	Run_clustal(joined_file)
	clu_score = joined_file.split(".")[0] + "_ClustalScore.txt"
	Analize_clustal_score(clu_score)"""

	

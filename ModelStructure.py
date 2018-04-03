if __name__ == "__main__":

	from Modules import FileParsersGenerators, GeneralFunctions, ProteinWorkingFunctions, RunningAnalyzingPrograms
	#from Bio.PDB import PDBParser, PDBIO, PPBuilder, Superimposer,PDBList, NeighborSearch, Selection
	#from Bio.Blast.Applications import NcbipsiblastCommandline as Ncbicmd
	#from Bio.Align.Applications import ClustalwCommandline
	import argparse
	import copy
	#import numpy
	#import re
	import sys

	parser = argparse.ArgumentParser(description="A program to model protein structures. It uses as an imput pairs of interacting chains (in PDB format) and models the whole structure by superimposition with the best templates.")

	parser.add_argument('-i', '--input',
				dest = "infiles",
				action = "store",
				required = True,
				nargs = "+",
				help = "Input file names (PDB files with pairs of chains that interact), insert as many as needed.")
	parser.add_argument('-d', '--database',
				dest = "database",
				action = "store",
				required = True,
				help = "Path to the pdb database in your computer.")
	parser.add_argument('-v', '--verbose',
				dest = "verbose",
				action = "store",
				default = False,
				help = "Print progression log to the standard error.")
	
	options=parser.parse_args()

	# INITIALIZATING SOME VARIABLES

	# Dictionary with all interactions between chains and correspondency with target-template chians
	Final_interactions = {}
	Final_interactions["target_interacts"] = {}
	Final_interactions["temps"] = {}
	# Input chains
	already_added = []
	unique_chains = []
	# Templates
	BLAST_outs = []
	fasta_names = []
	temp_chains = {}
	#Output
	final_files = []
	correct_predictions = []

	# WORKING WITH THE INPUT

	if options.verbose:
		sys.stderr.write("Parsing input files...\n")

	# Parse input files.
	(PDB_input_objects, PDB_input_names) = FileParsersGenerators.ParsePDB(options.infiles)
	file_prefixes = FileParsersGenerators.SplitChain(PDB_input_objects)

	# Save the names of the input chains.
	first = True
	for pref in file_prefixes:
		if first:
			unique_chains.append(pref)
			already_added.append(pref.split("_")[1])
			first = False
		if pref.split("_")[1] not in already_added:
			unique_chains.append(pref)
			already_added.append(pref.split("_")[1])

	# Parse the pdb files with single chains
	bychain_PDBs = map(lambda x: x + ".pdb", unique_chains)
	(PDB_bychain_objects, PDB_bychain_names) = FileParsersGenerators.ParsePDB(bychain_PDBs)
	# Add data to Final_interactions dictionary
	for inp in PDB_input_objects:
		inp_chains = inp.get_chains()
		inp_chains_ids = list(map(lambda x: x.get_id(), inp_chains))
		if inp_chains_ids[0] not in Final_interactions["target_interacts"].keys():
			Final_interactions["target_interacts"][inp_chains_ids[0]] = list(inp_chains_ids[1])
		elif (inp_chains_ids[0] in Final_interactions["target_interacts"].keys()) and (inp_chains_ids[1] not in Final_interactions["target_interacts"][inp_chains_ids[0]]):
			Final_interactions["target_interacts"][inp_chains_ids[0]].append(inp_chains_ids[1])

	# WORKING WITH TEMPLATES

	if options.verbose:
		sys.stderr.write("Searching templates...\n")

	# Look for templates
	for prefix in file_prefixes:
		output = RunningAnalyzingPrograms.RunBLAST(options.database, prefix)
		BLAST_outs.append(output)
	Templates = ProteinWorkingFunctions.SelectTemplate(BLAST_outs)

	# Downloading, parsing and spliting by chain the templates
	for template in Templates:
		ProteinWorkingFunctions.DownloadTemplate(template)
	temp_PDBs = map(lambda x: "pdb" + x + ".ent", Templates)
	(PDB_temp_objs, PDB_temp_names) = FileParsersGenerators.ParsePDB(temp_PDBs)	
	template_chains = FileParsersGenerators.SplitChain(PDB_temp_objs)
	bychain_PDBs = map(lambda x: x + ".pdb", template_chains)
	(PDB_chain_temp_objs, PDB_chain_temp_names) = FileParsersGenerators.ParsePDB(bychain_PDBs)

	# MODELING

	if options.verbose:
		sys.stderr.write("Running ClustalW...\n")

	# Creating a fasta file for each template chain adding all the input chains
	for chain in PDB_chain_temp_objs:
		obj_list = copy.copy(PDB_bychain_objects)
		obj_list.append(chain)
		joined_file = FileParsersGenerators.CreateJoinedFastas(obj_list)
		fasta_names.append(joined_file)
	# Performing a multiple alignment with ClustalW
	RunningAnalyzingPrograms.RunClustal(fasta_names)
	# Add data to Final_interactions dictionary
	first = True
	for fa_name in fasta_names:
		temp_name = fa_name.split("_")[-3] + "_" + fa_name.split("_")[-2]
		if first:
			tmp = GeneralFunctions.GetNameWOChain(temp_name)
			Final_interactions["temps"][GeneralFunctions.GetNameWOChain(temp_name)] = {}
			Final_interactions["temps"][GeneralFunctions.GetNameWOChain(temp_name)]["target_temp"] = {}
			first = False
		temp_chains[temp_name] = None
		aligns = RunningAnalyzingPrograms.AnalizeClustalScore(fa_name.split(".")[0] + "ClustalScore.txt", temp_name, 100)
		if len(aligns) == 0:
			aligns = RunningAnalyzingPrograms.AnalizeClustalScore(fa_name.split(".")[0] + "ClustalScore.txt", temp_name, 90)
			if len(aligns) == 0:
				aligns = RunningAnalyzingPrograms.AnalizeClustalScore(fa_name.split(".")[0] + "ClustalScore.txt", temp_name, 50)
				if len(aligns) == 0:
					print("%s template is not good enough to trust its model, the program will continue but we recommend you to consider using another approach to solve your problem if this one is the only good model obtained.")
					aligns = RunningAnalyzingPrograms.AnalizeClustalScore(fa_name.split(".")[0] + "ClustalScore.txt", temp_name, 0)
		if tmp == GeneralFunctions.GetNameWOChain(temp_name):
			Final_interactions["temps"][GeneralFunctions.GetNameWOChain(temp_name)]["target_temp"][temp_name] = aligns
		else:
			Final_interactions["temps"][GeneralFunctions.GetNameWOChain(temp_name)] = {}
			Final_interactions["temps"][GeneralFunctions.GetNameWOChain(temp_name)]["target_temp"] = {}
			Final_interactions["temps"][GeneralFunctions.GetNameWOChain(temp_name)]["target_temp"][temp_name] = aligns
			tmp = GeneralFunctions.GetNameWOChain(temp_name)
			
		Final_interactions["temps"][temp_name[:-2]]["temp_interact"] = {}
		temp_obj = list(filter(lambda x: x.get_id() == GeneralFunctions.GetNameWOChain(temp_name), PDB_temp_objs))
		list_interacts = ProteinWorkingFunctions.FindInteractions(temp_obj[0], True)
		for interact in list_interacts:
			if interact[0] in Final_interactions["temps"][temp_name[:-2]]["temp_interact"].keys():
				Final_interactions["temps"][GeneralFunctions.GetNameWOChain(temp_name)]["temp_interact"][interact[0]].update(interact[1])
			else:
				Final_interactions["temps"][GeneralFunctions.GetNameWOChain(temp_name)]["temp_interact"][interact[0]] = set(interact[1])
	
	if options.verbose:
		sys.stderr.write("Assigning target-template chain relations...\n")
	
	# Using a backtracking to assign chain relations (target-template)
	for temp in PDB_temp_names:
		ProteinWorkingFunctions.I_AssignQueryToTemp(PDB_bychain_names, temp_chains, Final_interactions, temp)

	if options.verbose:
		sys.stderr.write("Superimposing chains...\n")

	# Superimposing the chains
	for temp_obj in PDB_temp_objs:
		try:
			ProteinWorkingFunctions.SuperimposeChains(final_files, temp_obj, PDB_bychain_objects, temp_chains)
		except Exception:
			pass

	# FINAL ANALYSIS OF THE OBTAINED MODELS

	if options.verbose:
		sys.stderr.write("Omitting models that contain clashes...\n")

	# Analyzing obtained models
	(PDB_final_objects, PDB_final_names) = FileParsersGenerators.ParsePDB(final_files)
	for final_obj, final_name in zip(PDB_final_objects, PDB_final_names):
		if not ProteinWorkingFunctions.FindInteractions(final_obj, False):
			correct_predictions.append(final_name)

	if options.verbose:
		sys.stderr.write("The programme finished correctly.\n")

	if len(correct_predictions) >= 1:
		print("The generated files with the models are called: ")
		for out_f in correct_predictions:
			print("%s.pdb" %(out_f))
		print("Feel free to further analyse the models and choose the one that you find more accurate.")
	else:
		print("Good models couldn't be found.")
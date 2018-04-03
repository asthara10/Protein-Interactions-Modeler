from Bio.Blast.Applications import NcbipsiblastCommandline as Ncbicmd
from Bio.Align.Applications import ClustalwCommandline
import re

def RunBLAST(db_path, prefix):
	"""
	Runs a psiBLAST application.

	Arguments:

	db_path: path to the pdb database.
	prefix: prefix of the input files to analyze and of the output files that will be generated.
	"""

	fasta = prefix + ".fa"
	xml = prefix + ".xml"
	pssm = prefix + "_pssm"

	Psi_query = Ncbicmd('psiblast', db = db_path, query = fasta, out = xml, evalue = 10, outfmt = 3, out_pssm = pssm)
	Psi_query() # Run BLAST.

	return xml

def RunClustal(fastas):
	"""
	Performs a multiple alignment running ClustalW.

	Arguments:

	fastas: list of FASTA files with all sequences to be aligned, ClustalW will be run for each one.
	"""

	for fasta in fastas:
		clustalw_cline = ClustalwCommandline("clustalw2", infile=fasta)
		stdout, stderr = clustalw_cline()
		with open(fasta.split(".")[0] + "ClustalScore.txt", 'w') as scores:
			scores.write(stdout)

def AnalizeClustalScore(sc_file, temp_name, score):
	"""
	Analyzes ClustalW output score files. 
	Selects the chains that align with an specific score or higher.

	Arguments:

	sc_file: name of the file containing the scores.
	temp_name: template name.
	score: desired score treshold.
	"""

	equivalences = {}
	aligns = []

	file = open(sc_file, 'r')

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
			if (template == alig[2]) and (int(alig[6]) >= score):
				aligns.append(equivalences[alig[4]]) 
			elif (template == alig[4]) and (int(alig[6]) >= score):
				aligns.append(equivalences[alig[2]])

	return aligns

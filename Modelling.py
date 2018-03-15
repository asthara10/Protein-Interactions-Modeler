import Blast_align
from Bio.Align.Applications import ClustalwCommandline
from Bio.PDB import PDBList

input_files = ["AB.pdb", "AC.pdb", "AD.pdb", "BC.pdb", "DC.pdb"]
Blast_outs = ["AB_A.fa.xml", "AB_B.fa.xml", "AD_D.fa.xml", "AC_C.fa.xml"]
templates = Blast_align.Select_template(Blast_outs)

templates = list(templates)
templates = [x.upper() for x in templates]

pdbl = PDBList()

for template in templates:
	pdbl.retrieve_pdb_file(template, pdir='PDB')

def Run_clustal(input_files):
	for file in input_files:
		cline = ClustalwCommandline("clustalw2", infile=file)
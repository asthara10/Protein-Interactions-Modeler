import Blast_align
from Bio.Align.Applications import ClustalwCommandline
from Bio.PDB import PDBList
from Bio.PDB import PDBParser, PDBIO, PPBuilder

input_files = ["AB.pdb", "AC.pdb", "AD.pdb", "BC.pdb", "DC.pdb"]
Blast_outs = ["AB_A.fa.xml", "AB_B.fa.xml", "AD_D.fa.xml", "AC_C.fa.xml"]
#input_fasta = ["AB_A.fa", "AB_B.fa", "AD_D.fa", "AC_C.fa"]
templates = Blast_align.Select_template(Blast_outs)

templates = list(templates)
#templates = [x.upper() for x in templates]

pdbl = PDBList()
io = PDBIO()
polipeptide = PPBuilder()
#for inputfa in input_files:

for template in templates:
	pdbl.retrieve_pdb_file(template, pdir="./", file_format='pdb', obsolete=False)#al fer maco der tot mb carpetes
	pdb = PDBParser().get_structure(template, "pdb"+template+".ent")
	for infa in input_files:
		pd1 = PDBParser().get_structure(infa.split(".")[0], infa)
		for pp in polipeptide.build_peptides(pdb):
			for qq in polipeptide.build_peptides(pd1):
				new_fasta = open(template + "_" + infa.split(".")[0]+".fa", "w")
				new_fasta.write(">" + template +"\n")
				new_fasta.write(str(pp.get_sequence()))
				new_fasta.write("\n>" + infa.split(".")[0] +"\n")
				new_fasta.write(str(qq.get_sequence()))
	



def Run_clustal(input_files):
	for file in input_files:
		cline = ClustalwCommandline("clustalw2", infile=file)

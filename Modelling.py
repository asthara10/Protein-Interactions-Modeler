import Blast_align
from Bio.Align.Applications import ClustalwCommandline
from Bio.PDB import PDBList
from Bio.PDB import PDBParser, PDBIO, PPBuilder

input_files = ["AB.pdb", "AC.pdb", "AD.pdb", "BC.pdb", "DC.pdb"]
Blast_outs = ["AB_A.fa.xml", "AB_B.fa.xml", "AD_D.fa.xml", "AC_C.fa.xml"]
input_fasta = ["AB_A.fa", "AB_B.fa", "AD_D.fa", "AC_C.fa"]
templates = Blast_align.Select_template(Blast_outs)

templates = list(templates)
#templates = [x.upper() for x in templates]
def Create_Fastas(templates, input_files):
	pdbl = PDBList()
	io = PDBIO()
	polipeptide = PPBuilder()
	#for inputfa in input_files:

	for template in templates:
		first_temp = True
		#pdbl.retrieve_pdb_file(template, pdir="./", file_format='pdb', obsolete=False)#al fer maco der tot mb carpetes
		pdb = PDBParser().get_structure(template, "./ent/pdb"+template+".ent")
		for pp in polipeptide.build_peptides(pdb):
			if first_temp:
				new_fasta = open(template + "_paired.fa", "w")
				new_fasta.write(">" + template +"\n")
				new_fasta.write(str(pp.get_sequence()))
				first_temp = False
			else:
				new_fasta.write(str(pp.get_sequence()))
		for infa in input_files:
			first_chain = True
			new_fasta = open(template + "_paired.fa", "a")
			pd1 = PDBParser().get_structure(infa.split(".")[0], infa)
			for qq in polipeptide.build_peptides(pd1):
				if first_chain:
					new_fasta.write("\n>" + infa.split(".")[0] +"\n")
					new_fasta.write(str(qq.get_sequence()))
					first_chain = False
				else:
					new_fasta.write(str(qq.get_sequence()))
	

fastas = ["1a0u_paired.fa","1a0z_paired.fa","1dxt_paired.fa","1dxu_paired.fa","1gli_paired.fa","1j7s_paired.fa","1o1l_paired.fa","1o1n_paired.fa","1y0t_paired.fa","1y0w_paired.fa"]

def Run_clustal(fastas):
	for fasta in fastas:
		#aln = open(fasta.split(".")[0]+".aln", "w")
		clustalw_cline = ClustalwCommandline("clustalw2", infile=fasta)
		clustalw_cline()

Create_Fastas(templates, input_files)
Run_clustal(fastas)
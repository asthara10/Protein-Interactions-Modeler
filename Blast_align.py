from Bio.Blast.Applications import NcbipsiblastCommandline as Ncbicmd
from Bio.Blast import NCBIXML
import argparse

def Run_BLAST():
	"""
	blah blah blah
	"""
	Psi_query = Ncbicmd('psiblast', db = db_path, query = Filename+".fasta", out = Filename+".xml", evalue = 10, outfmt = 3, out_pssm = Filename+"_pssm")

	Psi_query() #en teoria ja tenim l'output del blast i està be


def Select_template():
	"""
	blah blah blah
	"""
	Blast_out = open (Filename+".xml")
	Blast_records = NCBIXML.parse(Blast_out)
	

if __name__ == "__main__":

	parser = argparse.ArgumentParser(description="Make a psiBLAST, choose the best template and align the sequences.")

	parser.add_argument('-i', '--input',
				dest= "Filename",
				action= "store",
				default= os.getcwd(),
				help="Input FASTA file and prefix of output files.")

	options = parser.parse_args()

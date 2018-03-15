from Bio.Blast.Applications import NcbipsiblastCommandline as Ncbicmd
import argparse

def Run_BLAST(db_path, Filename):
	"""
	blah blah blah
	"""
	Psi_query = Ncbicmd('psiblast', db = db_path, query = Filename, out = Filename+".xml", evalue = 10, outfmt = 3, out_pssm = Filename+"_pssm")

	Psi_query() # Run BLAST



def Select_template(Blast_outs):
	"""
	blah blah blah
	"""
	Outputs = {}

	for Out in Blast_outs:
		Blast_out = open(Out)
		First = True
		for line in Blast_out:
			if line.startswith("Sequences producing significant alignments:"):
				Blast_out.readline()
				for line in Blast_out:
					line = line.strip()
					line = line.split()
					if First:
						Max_evalue = line[len(line)-1]
						listoftups = [(line[0][:-2], line[len(line)-1])]
						Outputs[Out] = listoftups
						First = False
					else:
						if line[len(line)-1] == Max_evalue:
							listoftups.append((line[0][:-2], line[len(line)-1]))
							Outputs[Out] = listoftups
						else:
							break

	min_evalue = min(map(lambda x: min(map(lambda y: y[1], x)), Outputs.values()))
	templates = set()
	for chain in Outputs.values():
		for tup in chain:
			if tup[1] == min_evalue:
				templates.add(tup)

	print(templates)


		
if __name__ == "__main__":

	parser = argparse.ArgumentParser(description="Make a psiBLAST, choose the best template and align the sequences.")

	parser.add_argument('-i', '--input',
				dest = "Filenames",
				action = "store",
				nargs = '+',
				required = True,
				help ="Input FASTA file and prefix of output files.")

	options = parser.parse_args()

	db_path = "./databases/pdb_seqres.txt"

	Blast_outs = []

	for Filename in options.Filenames:
		#Run_BLAST(db_path, Filename)
		Blast_outs.append(Filename+".xml")

	Select_template(Blast_outs)
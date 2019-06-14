import os
import pandas as pd  
from Bio import Entrez, SeqIO
from os.path import dirname, abspath, isfile
Entrez.email = 'td1515@ic.ac.uk'
Entrez.api_key = '41267f8592172caaa22ab00ec006c4330208'

parentdir = dirname(dirname(abspath(__file__)))

# OPEN RESULTS TABLE:
def open_results(query):
	filepath = parentdir + '/3_bidirectional-best-hits/results-Cj' + query + '-names.xlsx'
	bbh = pd.read_excel(filepath, sheet_name='BBH-eval-filter')
	bbh = bbh.loc[:, ['saccver', 'subject_gcf', 'organism_name', 'infraspecific_name', 'taxid', 'species_taxid']]
	psi = pd.read_excel(filepath, sheet_name='psiblast-eval-filter')
	psi = psi.loc[:, ['saccver', 'subject_gcf', 'organism_name', 'infraspecific_name', 'taxid', 'species_taxid']]
	return bbh, psi

def get_sequence(hit):
	handle = Entrez.efetch(db="protein", id=hit, rettype="fasta", retmode="text") # fetch the full sequence of that acc.ver
	record = SeqIO.read(handle, "fasta") # read it out
	handle.close()
	return record

def make_mfa(bbh_df, psi_df, query):
	i = 0
	bbh_hits = [i for i in bbh_df['saccver']]
	with open(query+'-BBH-proteins.mfa', 'a+') as f:
		for hit in bbh_hits:
			seq = get_sequence(hit)
			i += 1
			print(i)
			SeqIO.write(seq, f, 'fasta')
	print('made ', query, 'bbh.fa')

	psi_hits = [i for i in psi_df['saccver']]
	i = 0
	with open(query+'-psiblast-proteins.mfa', 'a+') as f:
		for hit in psi_hits:
			seq = get_sequence(hit)
			i += 1
			print(i)
			SeqIO.write(seq, f, 'fasta')
	print('made ', query, 'psi.fa')

bbh_a, psi_a = open_results('PflA')
make_mfa(bbh_a, psi_a, 'PflA')

bbh_b, psi_b = open_results('PflB')
make_mfa(bbh_b, psi_b, 'PflB')


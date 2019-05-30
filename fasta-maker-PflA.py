import os
import pandas as pd  
from Bio import Entrez, SeqIO
from os.path import dirname, abspath, isfile
Entrez.email = 'td1515@ic.ac.uk'

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# OPEN PSIBLAST-BLASTP-MATCHES AND UNIQUE-PSIBLST-HITS TABLES:
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

query = 'CjPflA'
parentdir = dirname(dirname(abspath(__file__)))
psi_nonbbh = parentdir + '/best-hits/psiblast-blastp-compare/' + query + '_psiblast_blastp_noBBHs.tsv' #path to table of psiblast-blastp hits
psi_unique = parentdir + '/best-hits/psiblast-blastp-compare/' + query + '_psiblast_novel_finds.tsv' #path to unique psiblast hits

nonbbhs = pd.read_csv(psi_nonbbh, sep='\t') #open psiblast-blastp hits as pd.df
uniques = pd.read_csv(psi_unique, sep='\t') #open unique psiblast hits as pd.df

e_cut = 1e-26
e = "{0:.0e}".format(e_cut)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# MAKE FILES TO STORE THE FASTA SEQUENCES:
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def mknew_faa(filename):
	file_exists = isfile(filename) # check if results.tsv already exists
	if file_exists: #if it existed, delete it
		os.remove(filename)
	fa = open(filename, "a+") # open (or create) results.tsv
	return(fa)

bbh_fa = mknew_faa(query + '_psi-nonbbh_evalue_' + e + '.faa') # sequences of psiblast-blastp hits will be stored here
psi_fa = mknew_faa(query + '_psi-unique_evalue_' + e + '.faa') # sequences of unique psiblast hits will be stored here

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# EXTRACT THE SACC.VERS AND APPEND:
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#first filter the hits tables for evalue
nonbbhs = nonbbhs.loc[nonbbhs['evalue'] < e_cut]
uniques = uniques.loc[uniques['evalue'] < e_cut]

#then filter to just get the sacc.ver list 
nonbbhs = nonbbhs.saccver
uniques = uniques.saccver

#function that will fetch the sequences and append them to the correct file
def fa_writer(ls, fa_file):
	for acc in ls: # top_hit_id stores the acc.ver numbers of the blast top hits
	    handle = Entrez.efetch(db="protein", id=acc, rettype="fasta", retmode="text") # fetch the full sequence of that acc.ver
	    record = SeqIO.read(handle, "fasta") # read it out
	    handle.close()
	    SeqIO.write(record, fa_file, "fasta") # save it a a file, into the appropriate subdiresctory

fa_writer(nonbbhs, bbh_fa) # take acc.vers from nonbbhs, write their sequences to bbh_fa
print('Finished writing psiblast-blastp sequences')
fa_writer(uniques, psi_fa) # take acc.vers from uniques, write their sequences to psi_fa
print('Finished writing unique psiblast sequences')

bbh_fa.close()
psi_fa.close()

	

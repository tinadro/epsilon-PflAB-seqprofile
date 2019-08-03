import pandas as pd 
import numpy as np 
from Bio import Entrez, SeqIO
Entrez.email = 'td1515@ic.ac.uk'
Entrez.api_key = '41267f8592172caaa22ab00ec006c4330208'

def get_taxid(prot_acc):
	handle = Entrez.esummary(db='protein', id=prot_acc, report='full')
	record = Entrez.read(handle)
	handle.close()
	return record[0]['TaxId']

def check_family(txid):
	handle = Entrez.efetch(db='taxonomy', id=txid)
	record = Entrez.read(handle)
	handle.close()
	return record[0]['LineageEx']

def family_fasta(pfl):
	with open(pfl+'-campylobacteraceae.mfa', 'a+') as campy, open(pfl+'-arcobacter.mfa', 'a+') as arco, open(pfl+'-helicobacteraceae.mfa', 'a+') as helico, open(pfl+'-hydrogenimonaceae.mfa', 'a+') as hydro, open(pfl+'-nautiliaceae.mfa', 'a+') as nautilia, open(pfl+'-nitratiruptor.mfa', 'a+') as nitratiruptor, open(pfl+'-sulfurovum.mfa', 'a+') as sulfurovum, open(pfl+'-thioreductor.mfa', 'a+') as thioreductor:
		count =0
		for record in SeqIO.parse(pfl+'-proteins.mfa', 'fasta'):
			print(count)
			count += 1
			if 'tblastn' in record.id:
				pass
			else:
				acc = record.id
				taxid = get_taxid(acc)
				taxinfo = check_family(taxid)
				for ele in taxinfo:
					for key in ele.keys():
						if ele[key] == '194': # Campylobacter
							SeqIO.write(record, campy, 'fasta')
						elif ele[key] == '57665': # sulfurospirillum, also campylobacteraceae
							SeqIO.write(record, campy, 'fasta')
						elif ele[key] == '2321108': # arcobacter group
							SeqIO.write(record, arco, 'fasta')
						elif ele[key] == '72293': # helicobacteraceae
							SeqIO.write(record, helico, 'fasta')
						elif ele[key] == '292630': # hydrogenimonaceae
							SeqIO.write(record, hydro, 'fasta')
						elif ele[key] == '224467': # nautiliales
							SeqIO.write(record, nautilia, 'fasta')
						elif ele[key] == '269258': # nitraruptor
							SeqIO.write(record, nitratiruptor, 'fasta')
						elif ele[key] == '265570': # sulfurovum
							SeqIO.write(record, sulfurovum, 'fasta')
						elif ele[key] == '269252': # thioreductor
							SeqIO.write(record, thioreductor, 'fasta')

family_fasta('PflA')
family_fasta('PflB')

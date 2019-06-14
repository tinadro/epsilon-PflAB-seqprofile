from Bio import SeqIO
import sys

filename = sys.argv[1] # name of file to filter
f = '.'.join(filename.split('.')[:-1]) # file name without extension

uniques = f+'_nr.mfa'

with open(uniques, 'a+') as outfile:
	ids = []
	for record in SeqIO.parse(filename, 'fasta'):
		if record.id not in ids:
			ids.append(record.id)
			SeqIO.write(record, outfile, 'fasta')


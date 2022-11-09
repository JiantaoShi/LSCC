import argparse
import sys
import re
import numpy as np
from cyvcf2 import VCF

parser = argparse.ArgumentParser()
parser.add_argument('-I', "--iFile", required=True,  help="Input maf file.")
parser.add_argument('-S', "--SNV", required=True, help="Annotated SNV file in Vcf format.")
parser.add_argument('-D', "--cScore", required=True, type=int, help="CADD score cutoff.")
parser.add_argument('-O', "--oFile", required=True, help="Output file.")
args = parser.parse_args()

keywords = ['Hugo_Symbol', 
			'Chromosome', 
			'Start_position', 
			'End_position', 
			'Reference_Allele', 
			'Tumor_Seq_Allele2', 
			'Variant_Classification', 
			'Variant_Type', 
			'COSMIC_n_overlapping_mutations', 
			'T_ref_count', 
			'T_alt_count', 
			'N_ref_count', 
			'N_alt_count', 
			'context', 
			'categ', 
			'Tumor_Sample_Barcode']

hMap = {}
header = []

with open(args.iFile) as f:
	for text in f:
		chunks  = text.split("\n")[0].split("\t")
		if chunks[0] == 'Hugo_Symbol':
			for i in range(len(chunks)):
				key = chunks[i]
				if key in keywords:
					hMap[key] = i
					header.append(key)
			break

Vcf = VCF(args.SNV)

OUT = open(args.oFile,"w")
OUT.write('\t'.join(header) + '\tCADD\teffect\n')

with open(args.iFile) as f:
	for text in f:
		chunks  = text.split("\n")[0].split("\t")
		if chunks[0] == 'Hugo_Symbol':
			continue
		if chunks[hMap['Variant_Type']] != 'SNP':
			continue
		xout = []
		for key in header:
			xout.append(chunks[hMap[key]])
		VC = chunks[hMap['Variant_Classification']]
		ALT = chunks[hMap['Tumor_Seq_Allele2']]
		IV = chunks[hMap['Chromosome']] + ':' + chunks[hMap['Start_position']] + '-' + chunks[hMap['End_position']]
		if VC == 'Silent':
			effect = 'silent'
		elif VC == '5\'Flank':
			hits = 0
			for variant in Vcf(IV):
				if variant.ALT[0] == ALT:
					hits = 1
					CADD = float(variant.INFO.get('CADD'))
					break
			if hits == 0:
				continue
			if CADD > args.cScore:
				effect = 'nonsilent'
			else:
				effect = 'noncoding'
		else:
			continue
		xout.append(str(CADD))
		xout.append(effect)
		OUT.write('\t'.join(xout) + '\n')
OUT.close()


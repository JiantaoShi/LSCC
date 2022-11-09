import argparse
import sys
import re
import numpy as np
from cyvcf2 import VCF

parser = argparse.ArgumentParser()
parser.add_argument('-I', "--iFile", required=True,  help="Input VCF file with all SNVs.")
parser.add_argument('-G', "--gFile", required=True, help="Gene table.")
parser.add_argument('-C', "--cFile", required=True, help="Mapping file between mutation and category.")
parser.add_argument('-D', "--cScore", required=True, type=int, help="CADD score cutoff.")
parser.add_argument('-O', "--oFile", required=True, help="Output file.")
args = parser.parse_args()

# load gene in to map
gMap = {}
with open(args.gFile) as f:
    for text in f:
        chunks  = text.split("\n")[0].split("\t")
        gMap[chunks[0]] = re.split(':|-', chunks[1])

cMap = {}
with open(args.cFile) as f:
    for text in f:
        chunks  = text.split("\n")[0].split("\t")
        cMap[chunks[0]] = chunks[1]

# coverage map
rowDict = {'silent':0, 'nonsilent':1, 'noncoding':2}
colDict = {'1':0, '2':1, '3':2, '4':3, '5':4, '6':5, '7':6}

Vcf = VCF(args.iFile)
OUT = open(args.oFile,"w")
OUT.write('gene\teffect\tcateg\tcoverage\n')

for gene in gMap:
	print('Running ' + gene)
	queryIR = gMap[gene]
	CHR = queryIR[0]
	gStart = int(queryIR[1])
	gEnd   = int(queryIR[2])

	summaryM = np.zeros((3, 7), dtype = 'int')
	for variant in Vcf(CHR + ':' + str(gStart) + '-' + str(gEnd)):
		if variant.INFO.get('HS') != gene:
			continue
		VC = variant.INFO.get('VC')
		if VC not in ['5\'Flank', 'Silent']:
			continue
		xBC  = variant.INFO.get('L') + '(' + variant.REF + '->' + variant.ALT[0] + ')' + variant.INFO.get('R')
		if xBC not in cMap:
			print(xBC + ' not found in provided table.')
			continue
		category = cMap[xBC]
		CADD = float(variant.INFO.get('CADD'))
		if VC == 'Silent':
			effect = 'silent'
		elif VC == '5\'Flank':
			if CADD > args.cScore:
				effect = 'nonsilent'
			else:
				effect = 'noncoding'
		summaryM[rowDict[effect], colDict[category]] += 1
	summaryM[:, 6] = np.sum(summaryM, 1)

	for effect in rowDict:
		for category in colDict:
			if category == '7':
				continue
			value = summaryM[rowDict[effect], colDict[category]]
			OUT.write(gene + '\t' + effect + '\t' + category + '\t' + str(int(round(value/3.0))) + '\n')
OUT.close()




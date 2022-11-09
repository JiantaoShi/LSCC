import argparse
import sys
from cyvcf2 import VCF

parser = argparse.ArgumentParser()
parser.add_argument('-I', "--iFile", required=True, help="Input VCF file.")
parser.add_argument('-O', "--oFile", required=True, help="Output file.")
args = parser.parse_args()

OUT = open(args.oFile,"w")

chrMap   = {}
startMap = {}
endMap   = {}
countMap = {}
blackMap = set()

for variant in VCF(args.iFile):
	CHROM = variant.CHROM
	start = variant.start
	gene  = variant.INFO.get('HS')
	VC = variant.INFO.get('VC')

	if gene == 'Unknown':
		continue
	if VC != '5\'Flank':
		continue

	if gene in chrMap:
		if chrMap[gene] == CHROM:
			countMap[gene] += 1
			if start < startMap[gene]:
				startMap[gene] = start
			elif start > endMap[gene]:
				endMap[gene] = start
		else:
			print(gene + " found in two chromosomes.")
			blackMap.add(gene)
	else:
		chrMap[gene] = CHROM
		startMap[gene] = start
		endMap[gene] = start
		countMap[gene] = 1

for gene in chrMap:
	if gene not in blackMap:
		OUT.write(gene + "\t" + chrMap[gene] + ":" + str(startMap[gene]) + "-" + str(endMap[gene]) + "\t" + str(int(round(countMap[gene]/3))) + '\n')
OUT.close()

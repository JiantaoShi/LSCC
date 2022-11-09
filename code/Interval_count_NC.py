import argparse
import sys,re
from cyvcf2 import VCF

parser = argparse.ArgumentParser()
parser.add_argument('-I', "--iFile", required=True, help="SNV BD vcf file.")
parser.add_argument('-V', "--vFile", required=True, help="Interval file.")
parser.add_argument('-G', "--gFile", required=True, help="Output file.")
parser.add_argument('-O', "--oFile", required=True, help="Output file.")
args = parser.parse_args()

# load gene in to map
gMap = {}
with open(args.gFile) as f:
    for text in f:
        chunks  = text.split("\n")[0].split("\t")
        gMap[chunks[0]] = re.split(':|-', chunks[1])

# load Interval file
IV =  []
count = []
with open(args.vFile) as f:
    for text in f:
        text  = text.split("\n")[0]
        IV.append(text)

# Query data base
Vcf = VCF(args.iFile)
for iv in IV:
	num = 0
	for variant in Vcf(iv):
		if variant.INFO.get('HS') not in gMap:
			continue
		if variant.INFO.get('VC') != '5\'Flank':
			continue
		num += 1
	count.append(num)

# write to file
OUT = open(args.oFile,"w")
for i in range(len(IV)):
	OUT.write(IV[i] + '\t' + str(count[i]/3) + '\n')
OUT.close()

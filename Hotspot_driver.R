library("GenomicRanges")

# Mutations in 5' flanking regioins
Tx = strsplit(readLines("data/LUSC.maf"), "\t")[-1]
chr = sapply(Tx, function(x) x[2])
pos = as.numeric(sapply(Tx, function(x) x[3]))
PID = sapply(Tx, function(x) x[17])
Variant_Classification = sapply(Tx, function(x) x[7])
GR  <- GRanges(seqnames = Rle(chr),  ranges = IRanges(pos,pos), PID = PID)[Variant_Classification == "5'Flank"]

# find hotspot intevals
xgr        <- GR
start(xgr) <- start(xgr) - 25
end(xgr)   <- end(xgr) + 25
hGR        <- reduce(xgr)
hGR$n      <- countOverlaps(hGR, GR, type = "any", ignore.strand = TRUE)
hGR        <- hGR[hGR$n > 2]
write.table(cbind(as.character(hGR)), file = "out/Hotspot_interval.txt", row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
# Count number of possible non-coding mutations
# python code/Interval_count_NC.py -I data/db/WPS_SNV.vcf.gz -V out/Hotspot_interval.txt -G data/db/WPS_gene_Promoter.txt -O out/Hotspot_interval_count.txt
hGR$N     <- as.numeric(read.table(file = "out/Hotspot_interval_count.txt", row.names = NULL, header = FALSE, sep = "\t")[,2])

# Background mutation rate
promoterTable = read.table(file = "data/db/WPS_gene_Promoter.txt", row.names = NULL, header = FALSE, sep = "\t")
promoterSize  = sum(as.numeric(promoterTable[,3]))
bgRate = summary(factor(GR$PID), maxsum = 1000)/promoterSize

# check each interval seperately
Ntotal = length(unique(GR$PID))
rM = matrix(NA, length(hGR), 3)

for(i in 1:length(hGR)){

	gr  = hGR[i]
	xGR = GR[countOverlaps(GR, gr, type = "any", ignore.strand = TRUE) > 0]
	meanRate = mean(bgRate[unique(xGR$PID)])*width(gr)
	Nk = length(unique(xGR$PID))
	rM[i, 1] = Nk
	rM[i, 2] = signif(pnbinom(Ntotal - Nk, Nk, meanRate, lower.tail = TRUE, log.p = FALSE), 4)
}
rM[,3] = signif(p.adjust(rM[, 2]), 4)
colnames(rM) = c("k", "pValue", "qValue")

# the nearby gene
geneGR = GRanges(as.character(promoterTable[,2]))
geneGR$gene = as.character(promoterTable[, 1])
geneGR$N    = as.numeric(promoterTable[, 3])
x = findOverlaps(hGR, geneGR, type = "any", ignore.strand = TRUE)
hGR = hGR[queryHits(x)]
hGR$gene = geneGR[subjectHits(x)]$gene
xM = rM[queryHits(x), ]
mcols(hGR) = data.frame(mcols(hGR), xM)

# save results
SID    = as.character(hGR)
window = width(hGR)
oTable = data.frame(SID, window, mcols(hGR))
write.table(oTable, file = "out/Hotspot_interval_pvalue.txt", row.names = FALSE, col.names = TRUE, sep = "\t")

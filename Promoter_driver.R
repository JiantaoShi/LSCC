library("GenomicRanges")

# Mutations in 5' flanking regioins
Tx = strsplit(readLines("data/LUSC.maf"), "\t")[-1]
chr = sapply(Tx, function(x) x[2])
pos = as.numeric(sapply(Tx, function(x) x[3]))
PID = sapply(Tx, function(x) x[17])
Variant_Classification = sapply(Tx, function(x) x[7])
GR  <- GRanges(seqnames = Rle(chr),  ranges = IRanges(pos,pos), PID = PID)[Variant_Classification == "5'Flank"]

# read promoter intervals
promoterTable = read.table(file = "data/db/WPS_gene_Promoter.txt", row.names = 1, header = FALSE, sep = "\t")
colnames(promoterTable) = c("IV", "count")

# read gene table
gM = as.matrix(read.table(file = "data/db/WPS_gene_covariates.txt", row.names = 1, header = TRUE))
for(i in 1:3){
	value  = gM[,i]
	gM[,i] = (value - mean(value))/sd(value)
}

# test each region
top = round(nrow(gM)*0.05)
promoterGR = GRanges(as.character(promoterTable[,1]))
promoterGR$count = as.numeric(promoterTable[,2])
promoterGR$gene  = rownames(promoterTable)

Ntotal = length(unique(GR$PID))
pvalue = rep(NA, length(promoterGR))

for(i in 1:length(promoterGR)){

	xgr  = promoterGR[i]
	gene = xgr$gene

	if(!gene %in% rownames(gM))
		next

	xg  = gM[gene, ]
	xGM = gM[rownames(gM) != gene, ]
	distance = rep(NA, nrow(xGM))
	for(j in 1:nrow(xGM)){
		yg = xGM[j, ]
		distance[j] = sqrt(sum((xg - yg)^2))
	}
	names(distance) = rownames(xGM)
	probes = names(sort(distance)[1:top])

	Idx  = rownames(promoterTable) %in% probes
	bgGR = promoterGR[Idx]
	bgSize = sum(bgGR$count)
	nBG    = sum(countOverlaps(GR, bgGR, type = "any", ignore.strand = TRUE) > 0)
	nFG    = sum(countOverlaps(GR, xgr,  type = "any", ignore.strand = TRUE) > 0)
	bgRate = nBG/bgSize/Ntotal
	fgRate = 1 - (1 - bgRate)^(xgr$count)

	Nobs = length(unique(GR[countOverlaps(GR, xgr,  type = "any", ignore.strand = TRUE) > 0]$PID))

	pvalue[i] = pbinom(Nobs, Ntotal, fgRate, lower.tail = FALSE, log.p = FALSE)
	if(pvalue[i] < 0.00001)
		cat(gene, Nobs, pvalue[i], "\n")
}
promoterGR$n = countOverlaps(promoterGR, GR, type = "any", ignore.strand = TRUE)
promoterGR$pvalue = pvalue
promoterGR$qvalue = p.adjust(pvalue)

# save
SID = as.character(promoterGR)
oT  = data.frame(SID, mcols(promoterGR))
write.table(oT, file = "out/Promoter_sig.txt", row.names = FALSE, col.names = TRUE, sep = "\t")

library("tximport")
library("DESeq2")

args = commandArgs(TRUE)
folder = args[1]
saveto = args[2]

samples = file.path(folder, "samples.tsv")
samples = read.table(samples)

filenames = file.path(folder, samples$filename)
names(filenames) = rownames(samples)

tx2gene = file.path(folder, "tx2gene.tsv")
tx2gene = read.csv(tx2gene, sep = "\t", header = FALSE)
colnames(tx2gene) = c("tx", "gene_id", "gene_name")
tx2gene = tx2gene[, 1:2]

txi <- tximport(filenames, type = "salmon", tx2gene = tx2gene)
assertthat::are_equal(rownames(samples), colnames(txi$counts))

dds = DESeqDataSetFromTximport(txi, colData = samples, design = ~IFNb)

dds$IFNb = relevel(dds$IFNb, ref = "untreated")
dds = DESeq(dds)
coef = "IFNb_treated_vs_untreated"

# ISGs
res = results(dds, name = coef,
              alpha = 0.01, lfcThreshold = 0.585, altHypothesis = 'greater')
res.shrunk = lfcShrink(dds, coef = coef, res = res)
res.shrunk$padj[is.na(res.shrunk$padj)] = 1

filename = file.path(saveto, "GSE118926.csv")
write.csv(res.shrunk, filename)

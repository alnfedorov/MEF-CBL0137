library("tximport")
library("DESeq2")

args = commandArgs(TRUE)
folder = args[1]
saveto = args[2]

samples = file.path(folder, "samples.tsv")
samples = read.table(samples)

samples$IFNb = relevel(factor(samples$IFNb), ref="0h")
samples$selection = relevel(factor(samples$selection), ref="input")
samples$origin = factor(samples$origin)

filenames = file.path(folder, samples$filename)
names(filenames) = rownames(samples)

tx2gene = file.path(folder, "tx2gene.tsv")
tx2gene = read.csv(tx2gene, sep = "\t", header = FALSE)
colnames(tx2gene) = c("tx", "gene_id", "gene_name")
tx2gene = tx2gene[, 1:2]

txi <- tximport(filenames, type = "salmon", tx2gene = tx2gene)
assertthat::are_equal(rownames(samples), colnames(txi$counts))

samples$all = factor(paste(samples$IFNb, samples$selection, sep="_"))
dds = DESeqDataSetFromTximport(txi, colData = samples, design = ~all)

# RIP tests
for (IFNb in c("0h", "24h", "48h")) {
  ref = paste(IFNb, "input", sep="_")
  dds$all = relevel(samples$all, ref = ref)
  deseq = DESeq(dds)

  for (treatment in c("Z22", "IgG")) {
    coef = paste("all", IFNb, treatment, "vs", ref, sep = "_")
    res = results(deseq, name = coef,
                  alpha = 0.01, lfcThreshold = 1, altHypothesis = 'greater')
    res.shrunk = lfcShrink(deseq, coef = coef, res = res)
    res.shrunk$padj[is.na(res.shrunk$padj)] = 1

    filename = paste(coef, ".csv", sep = "")
    filename = file.path(saveto, filename)
    write.csv(res.shrunk, filename)
  }
}

# Drop all RIP-related experiments to avoid messing up with the normalization and variance estimation
dds = DESeqDataSetFromTximport(txi, colData = samples, design = ~IFNb)
dds = dds[,dds$selection == "input"]

# ISGs tests / DEGs tests
treatment = c("48h", "48h", "24h")
control = c("24h", "0h", "0h")
degs = c()
for (ind in 1:3) {
  dds$IFNb = relevel(dds$IFNb, ref = control[ind])
  deseq = DESeq(dds)

  # ISGs
  coef = paste("IFNb", treatment[ind], "vs", control[ind], sep = "_")
  res = results(deseq, name = coef, alpha = 0.01, lfcThreshold = 0.585, altHypothesis = 'greater')
  res.shrunk = lfcShrink(deseq, coef = coef, res = res)
  res.shrunk$padj[is.na(res.shrunk$padj)] = 1

  filename = paste(coef, ".csv", sep = "")
  filename = file.path(saveto, filename)
  write.csv(res.shrunk, filename)

  # DEGs
  res = results(deseq, name = coef, alpha = 0.05, lfcThreshold = 0.585, altHypothesis = 'greaterAbs')
  res$padj[is.na(res$padj)] = 1
  cond.degs = res[res$padj < 0.05, ]
  degs = c(degs, row.names(cond.degs))
}

# save normalized expression
dds$IFNb = relevel(dds$IFNb, ref = "0h")
deseq = DESeq(dds)

expression = counts(deseq, normalized=TRUE)
filename = file.path(saveto, "expression.normalized.csv")
write.csv(expression, filename)

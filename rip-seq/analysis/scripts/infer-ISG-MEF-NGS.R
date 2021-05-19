library("DESeq2")

args = commandArgs(TRUE)
folder = args[1]
saveto = args[2]

samples = file.path(folder, "samples.tsv")
samples = read.table(samples, header = TRUE, row.names = 1)
samples$condition = factor(samples$condition)
samples$filename = file.path(folder, samples$filename)

counts = file.path(folder, "counts.tsv")
counts = read.table(counts, header = TRUE, row.names = 1)
assertthat::are_equal(rownames(samples), colnames(counts))

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = samples,
                              design = ~condition)
dds$condition = relevel(dds$condition, ref = "untreated")
dds = DESeq(dds)

coef = "condition_IFNb_vs_untreated"
res = results(dds, name = coef,
              alpha = 0.1, lfcThreshold = 0.585, altHypothesis = 'greater')
res.shrunk = lfcShrink(dds, coef = coef, res = res)
res.shrunk$padj[is.na(res.shrunk$padj)] = 1

write.csv(res.shrunk, saveto)

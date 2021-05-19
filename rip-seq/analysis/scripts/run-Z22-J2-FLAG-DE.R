library("tximport")
library("DESeq2")

args = commandArgs(TRUE)
folder = args[1]
saveto = args[2]

samples = file.path(folder, "samples.tsv")
samples = read.table(samples)
samples$cell = gsub("-", "_", samples$cell)

samples$cell = factor(samples$cell)
samples$IFNb = factor(samples$IFNb)
samples$antibody = factor(samples$antibody)

filenames = file.path(folder, samples$filename)
names(filenames) = rownames(samples)

txi = tximport(filenames, txIn = FALSE, type = "salmon", geneIdCol = "Name")
assertthat::are_equal(rownames(samples), colnames(txi$counts))

for (cell in c("ADAR1_WT", "ADAR1_KO")) {
  for (ifnb in c("0h", "72h")) {
    for (antibody in c("Z22", "J2", "FLAG")) {
      # no such sample -> skip
      if (sum(samples$cell == cell &
                samples$IFNb == ifnb &
                samples$antibody == antibody) == 0) {
        next
      }
      # Use only data for the corresponding cells|condition, since recommended usage of all samples leads to issues
      # with library normalization and dispersion overestimation (it is relevant only for the RIP-seq data)
      dds = DESeqDataSetFromTximport(txi, colData = samples, design = ~antibody)
      keep = dds$IFNb == ifnb & dds$cell == cell
      dds = dds[, keep]
      dds$antibody = droplevels(dds$antibody)
      dds$antibody = relevel(dds$antibody, ref = "IgG")

      keep <- rowSums(counts(dds)) >= 10
      dds <- dds[keep,]
      dds <- DESeq(dds)

      coef = paste("antibody", antibody, "vs", "IgG", sep = "_")
      res = results(dds, name = coef,
                    alpha = 0.1, lfcThreshold = 0.585, altHypothesis = 'greater')
      res.shrunk = lfcShrink(dds, coef = coef, res = res)
      res.shrunk$padj[is.na(res.shrunk$padj)] = 1

      filename = paste(cell, ifnb, antibody, "vs-IgG.csv", sep = "-")
      filename = gsub("_", "-", filename)
      filename = file.path(saveto, filename)

      write.csv(res.shrunk, filename)
    }
  }
}

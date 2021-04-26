if (!require(circlize))
    install.packages("circlize")
library(circlize)


args = commandArgs(TRUE)
z22peaks = args[1]
flagpeaks = args[2]
overlaps = args[3]


postscript(args[4],
           horizontal = FALSE, onefile = FALSE, paper = "special", height = 12, width = 12)
par(cex=1.3, mar=c(0, 0, 3, 0))

# 1. Outer circle = chromosomes
circos.par("track.height" = 0.0825, "start.degree" = 90, "cell.padding" = c(0, 0, 0, 0))
circos.initializeWithIdeogram(plotType = NULL, species="mm10")

par(ps=15)
bg.col = rep(c("#D1D1D1", "#696969"), length.out=30)
circos.track(ylim = c(0, 1), bg.col=bg.col, bg.lty = 0, panel.fun = function(x, y) {
  chr = substring(CELL_META$sector.index, 4)
  xlim = CELL_META$xlim
  ylim = CELL_META$ylim

  if (chr == "X") {
    col = "white"
  } else if (chr == "Y") {
    col = "black"
  } else if (strtoi(chr) %% 2 == 1) {
    col = "black"
  } else {
    col = "white"
  }
  circos.text(mean(xlim), mean(ylim), chr, col = col, cex=1.3,
              facing = "inside", niceFacing = FALSE)
})


# 3. Inner circle - peaks
files = c(
  z22peaks,
  flagpeaks
)
colors = c(
  "#FF4474",
  "#6AD2F4"
)

for (i in 1:length(colors)) {
  bed = read.csv(files[i], header = FALSE, sep = "\t")[,c(1, 2, 3, 9)]
  colnames(bed) = c("chr", "start", "end", 'values')
  # bed$values = rep(1, nrow(bed))
  bed$values = ifelse(bed$values > 35, yes = 35, no = bed$values)

  col = colors[i]

  circos.genomicTrack(bed, ylim=c(0, 35), bg.border = "#BABDBF",
                      panel.fun = function(region, value, ...) {
                        circos.genomicLines(region, value, type = "h",
                                            lwd=0.25, col=col, ...)

                      })
}

bed = read.csv(overlaps, header = FALSE, sep = "\t")[,c(1, 2, 3)]
colnames(bed) = c("chr", "start", "end")
bed$values = rep(1, nrow(bed))
circos.genomicTrack(bed, ylim=c(0, 1), bg.border = "#BABDBF",
                    panel.fun = function(region, value, ...) {
                      circos.genomicLines(region, value, type = "h",
                                          lwd=0.25, col="#B09ED6", ...)

                    })
colors = c(colors, "#B09ED6")

# 4. Legend

# gaps
vgap = 0.125
# start positions
xstart = -0.4
ystart = vgap * 2.25
# box height/width
boxw = 0.085
boxh = 0.05
textadj = boxh / 20


labels = c("Z22 peaks", "FLAG peaks", "Z22 & FLAG peaks")
for (i in 1:length(labels)) {
  y = ystart - vgap * i
  rect(xstart, y, xstart + boxw, y + boxh, col=colors[i], border=NA)
  text(xstart + boxw * 1.5, y + textadj, labels=labels[i], adj = c(0,0))
}


y = ystart - vgap * 4
rect(xstart, y, xstart + boxw, y + boxh, col=bg.col[2], border=FALSE)
rect(xstart + boxw, y, xstart + boxw*2, y + boxh, col=bg.col[1], border=FALSE)
text(xstart + boxw * 2.5, y + textadj, labels="# chromosome", adj = c(0,0))

firsty = ystart - vgap + boxh * 0.5
lasty = ystart - vgap * 3 + boxh * 0.25

x = xstart + boxw * 8.5
lines(c(x, x), c(firsty, lasty))

centery = (firsty + lasty) / 2
text(x + boxw * 0.25, centery, labels="CBL0137 14h", adj=c(0.5,-0.1), srt=-90)
title("Z22 & FLAG ChIP-seq (CBL0137 14h)", font.main=1)
# text(0, -1.03, "blacklisted excluded", adj=c(0.5, 1), font=4, cex=1)
dev.off()

require(MapMan2GO)
message("USAGE: Rscript path/2/MapMan2GO/exec/plotDepths.R path/2/MapMan2GO")

#' Input command line arguments:
input.args <- commandArgs(trailingOnly = TRUE)

#' Calculate median, max and mean
meanBins.GO.depth <- vector()
for (i in 1:length(mm.2.go.df$MapManBin)){
  meanBins.GO.depth[i] <- mean(mm.2.full.desc$GO.depth[grep(pattern = mm.2.go.df$MapManBin[i], mm.2.full.desc$MapManBin, fixed =T)])
}
mm.2.go.df$mean.GO.depth <- meanBins.GO.depth

medianBins.GO.depth <- vector()
for (i in 1:length(mm.2.go.df$MapManBin)){
  medianBins.GO.depth[i] <- median(mm.2.full.desc$GO.depth[grep(pattern = mm.2.go.df$MapManBin[i], mm.2.full.desc$MapManBin, fixed =T)])
}
mm.2.go.df$median.GO.depth <- medianBins.GO.depth

maxBins.GO.depth <- vector()
for (i in 1:length(mm.2.go.df$MapManBin)){
  suppressWarnings( maxBins.GO.depth[i] <- max(mm.2.full.desc$GO.depth[grep(pattern = mm.2.go.df$MapManBin[i], mm.2.full.desc$MapManBin, fixed =T)]) )
}
mm.2.go.df$max.GO.depth <- maxBins.GO.depth

#' - Boxplot of GO.depths statistics
pdf(file.path(input.args[[1]], "inst", "GO.depthMapMan-BinsDistribution.pdf"))
  boxplot(mm.2.go.df$median.GO.depth, mm.2.go.df$max.GO.depth, mm.2.go.df$mean.GO.depth, ylab = "GO.depth", names=c("Median", "Maximum", "Mean"), col = "lightgrey", pch = '-', outline=FALSE)
dev.off()

pdf(file.path(input.args[[1]], "inst", "GO.depthGOtermsDistribution.pdf"))
  def.mar <- par('mar')
  m.1 <- def.mar; m.1[[1]] <- 0
  op <- par(mfcol = 2:1, mar=m.1)
  hist(mm.2.full.desc$GO.depth, col = "lightgrey", main = "GO.depth GO terms distribution", xlab = NULL)
  m.2 <- def.mar; m.2[[3]] <- 0
  par( mar=m.2 )
  boxplot(mm.2.full.desc$GO.depth, col = "lightgrey", horizontal = TRUE, frame = FALSE, pch = '|')
  par(op)
dev.off()

#' Scatter plot Number of genes vs GO depth
pdf(file.path(input.args[[1]], "inst", "NumberOfGenesVsGODepthMean.pdf"))
  plot(mm.2.go.df$n.genes, mm.2.go.df$mean.GO.depth, xlab="Number of genes per MapMan-Bin", ylab="GO.depth", pch=19)
dev.off()

pdf(file.path(input.args[[1]], "inst", "NumberOfGenesVsGODepthMedian.pdf"))
  plot(mm.2.go.df$n.genes, mm.2.go.df$median.GO.depth, xlab="Number of genes per MapMan-Bin", ylab="GO.depth", pch=19)
dev.off()

pdf(file.path(input.args[[1]], "inst", "NumberOfGenesVsGODepthMaximum.pdf"))
  plot(mm.2.go.df$n.genes, mm.2.go.df$max.GO.depth, xlab="Number of genes per MapMan-Bin", ylab="GO.depth", pch=19)
dev.off()

save(mm.2.go, mm.2.go.df, mm.2.full.desc, mm.desc.df, file = file.path(input.args[[1]], 
    "data", "MapManBins2GO.RData"))



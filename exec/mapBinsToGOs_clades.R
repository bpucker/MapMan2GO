require(MapMan2GO)
options(mc.cores = detectCores())

message("USAGE: Rscript path/2/MapMan2GO/exec/mapBinsToGOs.R path/2/MapMan2GO path/2/MercatorRawResults.tsf /path/2/SwissProtIDs SwissProtIDs_subset_name")

input.args <- commandArgs(trailingOnly = TRUE)

#' Perform the analysis for the subset

subset.name <- input.args[[4]]
clade.ID <- unlist((fread(input.args[[3]], header = FALSE, data.table = FALSE, 
    stringsAsFactors = FALSE)), use.names = FALSE)
mmvs.work <- mm.bins.vs.sprot[which(mm.bins.vs.sprot$Swissprot.Short.ID %in% clade.ID), 
    ]

#' Map MapMan-Bins to compound Gene Ontology Term Annotations (GOA):
mm.leaf.bins.work <- unique(mmvs.work$MapManBin)
mm.2.go.work <- setNames(mclapply(mm.leaf.bins.work, compoundGoAnnotationEntropy), 
    mm.leaf.bins.work)
mm.2.go.df.work <- Reduce(rbind, mclapply(names(mm.2.go.work), function(x) {
    y <- mm.2.go.work[[x]]
    data.frame(MapManBin = x, MapManBin.GO = y[["MapManBin.GO"]], Shannon.Entropy = y[["Shannon.Entropy"]], 
        n.GO = y[["n.GO"]], n.genes = y[["n.genes"]], median.n.GO = y[["median.n.GO"]], 
        stringsAsFactors = FALSE)
}))
#' Add a full description for each MapMan-Bin's GOA including GO-Term names:
go.terms.not.in.db <- c()
mm.2.full.desc.work <- Reduce(rbind, mclapply(names(mm.2.go.work), function(m.b) {
    m.b.gos <- Reduce(intersect, mm.2.go.work[[m.b]]$genes.goa)
    Reduce(rbind, mclapply(m.b.gos, function(g.id) {
        # 
        if (g.id %in% GO.OBO$id) {
            g.name <- GO.OBO$name[[g.id]]
            g.ancestors <- GO.OBO$ancestors[[g.id]]
            g.depth <- length(g.ancestors)
            g.ont <- if (g.depth > 1) {
                GO.OBO$name[[g.ancestors[[1]]]]
            } else {
                g.name
            }
            data.frame(MapManBin = m.b, GO.Term = g.id, GO.Name = g.name, GO.Ontolgy = g.ont, 
                GO.depth = g.depth, stringsAsFactors = FALSE)
        } else {
            go.terms.not.in.db <- c(go.terms.not.in.db, g.id)
            data.frame(MapManBin = m.b, GO.Term = g.id, GO.Name = NA, GO.Ontolgy = NA, 
                GO.depth = NA, stringsAsFactors = FALSE)
        }
    }))
}))
#' Warn about missing GO Terms:
if (length(go.terms.not.in.db) > 0) {
    message("WARNING: The following GO Terms were assigned to MapManBin(s), but had no entry in the installed 'GO.db' package:\n", 
        paste(sort(unique(go.terms.not.in.db)), collapse = ", "))
}


#' Add the MapMan Bins Descriptions to the above table:
mr.full <- readMercatorResultTable(input.args[[2]], add.go.terms = FALSE)
mm.desc.df.work <- unique(mr.full[, c("BINCODE", "NAME")])
names(mm.desc.df.work) <- c("MapManBin", "Description")


value = integer(0)
mm.2.full.desc.work$Bin.Description <- as.character(unlist(mclapply(mm.2.full.desc.work$MapManBin, 
    function(x) {
        if (identical(which(mm.desc.df.work$MapManBin == x), value)) {
            mm.2.full.desc.work$Bin.Description <- NA
        } else {
            mm.desc.df.work[which(mm.desc.df.work$MapManBin == x), "Description"]
        }
        
    })))



#' Some statistics:
#' - Histogram of Number of GO Terms per MapMan-Bin GOA     # when successful, repeat this plot, it was overwritten, now the name is corrected
fl.name <- paste(ref.set.name, "NumberOfGoTermsPerMapManBinGoaHist.pdf", sep = "_")
pdf(file.path(input.args[[1]], "inst", fl.name))
plotDistAsHistAndBox(mm.2.go.df.work$n.GO, "Number of GO Terms per MapMan-Bin GOA")
dev.off()
#' - Histogram of Entropies
fl.name <- paste(ref.set.name, "MapManBinGoaEntropyHist.pdf", sep = "_")
pdf(file.path(input.args[[1]], "inst", fl.name))
plotDistAsHistAndBox(mm.2.go.df.work$Shannon.Entropy, "Shannon Entropy of compound GO Annotations per MapMan-Bin")
dev.off()
#' - Histogram of Sizes in terms of number of genes
fl.name <- paste(ref.set.name, "GenesPerMapManBinHist.pdf", sep = "_")
pdf(file.path(input.args[[1]], "inst", fl.name))
plotDistAsHistAndBox(mm.2.go.df.work$n.genes, "Number of genes per MapMan-Bin")
dev.off()
#' - Number of genes vs Number of GO Terms in the Bin-GOA:                     # Do again
fl.name <- paste(ref.set.name, "NumberOfGenesVsNumberOfGoTermsInMapManBinGOA.pdf", 
    sep = "_")
pdf(file.path(input.args[[1]], "inst", fl.name))
plot(mm.2.go.df.work$n.genes, mm.2.go.df.work$n.GO, xlab = "Number of genes per MapMan-Bin", 
    ylab = "Number of GO Terms in MapMan-Bin-GOA", pch = 20)
dev.off()
#' - Histogram of number of MapMan-Bins sharing identical GOAs
fl.name <- paste(subset.name, "NumberOfMapManBinsSharingIdentGOAsHist.pdf", sep = "_")
pdf(file.path(input.args[[1]], "inst", fl.name))
x <- as.numeric(table(mm.2.go.df.work[which(mm.2.go.df.work$MapManBin.GO != ""), "MapManBin.GO"]))
hist(x, col = "lightgrey", xlab = "Number of Bins sharing identical GOAs", main = "Histogram of number of MapManBins sharing identical GOAs")
dev.off()


#' Save results:
fl.name <- paste(ref.set.name, "MapManBins2GO.txt", sep = "_")
write.table(mm.2.full.desc.work, file.path(input.args[[1]], "inst", fl.name), sep = "\t", 
    row.names = FALSE)

#' Save binary results:
rdata.file <- (file = file.path(input.args[[1]], "data", "MapManBins2GO.RData"))
ref.set.results <- c(paste("mm.2.go", ref.set.name, sep = "."), paste("mm.2.go.df", 
    ref.set.name, sep = "."), paste("mm.2.full.desc", ref.set.name, sep = "."), 
    paste("mm.desc.df ", ref.set.name, sep = "."))
assign(ref.set.results[[1]], mm.2.go.work)
assign(ref.set.results[[2]], mm.2.go.df.work)
assign(ref.set.results[[3]], mm.2.full.desc.work)
assign(ref.set.results[[4]], mm.desc.df.work)
appendToRData(list = ref.set.results, file = rdata.file)


message("DONE")

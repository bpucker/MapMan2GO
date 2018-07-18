require(MapMan2GO)
options(mc.cores = detectCores())

message("USAGE: Rscript path/2/MapMan2GO/exec/mapBinsToGOs_batch.R path/2/MercatorRawResults.tsf start_index stop_index path_partial_result_file")

input.args <- commandArgs(trailingOnly = TRUE)

#' Unique MapMan-Bins to annotate:
mm.leaf.bins.all <- unique(mm.bins.vs.sprot$MapManBin)

#' Subset according to input arguments:
start.i <- input.args[[2]]

stop.i <- input.args[[3]]

result.file <- input.args[[4]]

mm.leaf.bins <- mm.leaf.bins.all[start.i:stop.i]
rm(mm.leaf.bins.all)


mm.2.go <- setNames(mclapply(mm.leaf.bins, compoundGoAnnotationEntropy), mm.leaf.bins)
mm.2.go.df <- Reduce(rbind, mclapply(names(mm.2.go), function(x) {
    y <- mm.2.go[[x]]
    data.frame(MapManBin = x, MapManBin.GO = y[["MapManBin.GO"]], Shannon.Entropy = y[["Shannon.Entropy"]], 
        Shannon.Entropy.BP = y[["Shannon.Entropy.BP"]], Shannon.Entropy.CC = y[["Shannon.Entropy.CC"]], 
        Shannon.Entropy.MF = y[["Shannon.Entropy.MF"]], Shannon.Entropy.not.used = y[["Shannon.Entropy.not.used"]], 
        Shannon.Entropy.not.used.BP = y[["Shannon.Entropy.not.used.BP"]], Shannon.Entropy.not.used.CC = y[["Shannon.Entropy.not.used.CC"]], 
        Shannon.Entropy.not.used.MF = y[["Shannon.Entropy.not.used.MF"]], n.GO = y[["n.GO"]], 
        n.genes = y[["n.genes"]], median.n.GO = y[["median.n.GO"]], stringsAsFactors = FALSE)
}))
#' Add a full description for each MapMan-Bin's GOA including GO-Term names:
go.terms.not.in.db <- c()
mm.2.full.desc <- Reduce(rbind, mclapply(names(mm.2.go), function(m.b) {
    m.b.gos <- Reduce(intersect, mm.2.go[[m.b]]$genes.goa)
    Reduce(rbind, mclapply(m.b.gos, function(g.id) {
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
mr.full <- as.data.table( readMercatorResultTable(input.args[[1]], add.go.terms = FALSE) )
mm.desc.df <- unique(mr.full[, c("BINCODE", "NAME")])
names(mm.desc.df) <- c("MapManBin", "Description")

value = integer(0)
mm.2.full.desc$Bin.Description <- as.character(unlist(mclapply(mm.2.full.desc$MapManBin, 
    function(x) {
        if (identical(which(mm.desc.df$MapManBin == x), value)) {
            mm.2.full.desc$Bin.Description <- NA
        } else {
            mm.desc.df[which(mm.desc.df$MapManBin == x), "Description"]
        }
        
    })))


#' Save results:
save(mm.2.go, mm.2.go.df, mm.2.full.desc, mm.desc.df, file = result.file)


message("DONE")

require(MapMan2GO)
options(mc.cores = detectCores())

message("USAGE: Rscript path/2/MapMan2GO/exec/evalGoPredictionPerformances.R path/2/mercator_result.txt path/2/blast_result_table path/2/interpro_scan_results_preprocessed path/2/preprocessed_uniprot_goa_table gold_standard_goa path/2/MapMan2GO")

input.args <- commandArgs(trailingOnly = TRUE)


#' Use the set of non SwissProt gene that have non electronically made GO term
#' annotations (see Vignette for more details):
gold.standard.ids <- sanitizeAccession(readLines(file.path(input.args[[length(input.args)]], 
    "inst", "nonSwissProtGenesWithNonElectMadeGOAs_trEMBL_IDs.txt")))
options(MapMan2GO.reference.genes = gold.standard.ids)


#' Load the input GOA tables for the predictions and the gold standard genes:
ukb.goa <- fread(input.args[[4]], sep = "\t", header = FALSE, quote = "", na.strings = "", 
    stringsAsFactors = FALSE)
gold.standard.goas <- read.table(input.args[[5]], comment.char = "!", quote = "", 
    na.strings = "", sep = "\t", header = FALSE, stringsAsFactors = FALSE)[, c("V2", 
    "V5", "V7")]
colnames(gold.standard.goas) <- c("gene", "GO", "ECO")
exclude.eco <- c("IEA", "RCA", "ISS")
gold.standard.goas.trust <- gold.standard.goas[which(!gold.standard.goas$ECO %in% 
    exclude.eco), ]
gold.standard.goas.trust.w.anc <- unique(Reduce(rbind, mclapply(gold.standard.goas.trust$gene, 
    function(gs.gene) {
        gs.goa <- gold.standard.goas.trust[which(gold.standard.goas.trust$gene == 
            gs.gene), ]
        Reduce(rbind, lapply(1:nrow(gs.goa), function(gs.goa.i) {
            gs.goa.row <- gs.goa[gs.goa.i, ]
            data.frame(gene = gs.gene, GO = addAncestors(gs.goa.row$GO), ECO = gs.goa.row$ECO, 
                stringsAsFactors = FALSE)
        }))
    })))
#' Set the gold standard reference GO annotations and the columns in which to
#' lookup gold standard proteins' identifiers and GO annotations:
message("Using gold standard without ancestral terms and predictions as they come from Best Blast and InterProScan (AsIs method).")
options(MapMan2GO.reference.gene.annos = gold.standard.goas.trust)
options(MapMan2GO.rga.gene.col = "gene")
options(MapMan2GO.rga.anno.col = "GO")


#' Compute F1-Scores and Matthews correlation coefficient of competing GO predictors.
#' ##################################################################################

#' First consider the universe of all possible GO term annotations 'as is',
#' that is without adding ancestral GO terms.

best.blast.tbl.non.sprot <- extractBestBlastHits(input.args[[2]], NULL)
blast.hit.goa <- as.data.frame(ukb.goa[ukb.goa[["V3"]] %in% best.blast.tbl.non.sprot$hit.ukb.short.id, 
    ])
best.blast.pred <- bestBlastPredictions(best.blast.tbl.non.sprot, blast.hit.goa)
#' Performance on the whole of the Gene Ontology:
options(MapMan2GO.performance.universe.annotations = GO.OBO$id)
bb.annos.no.ref.anc.non.sprot <- predictorPerformance(best.blast.pred, pa.gene.col = "query", 
    pa.anno.col = "GO", process.annos.funk = identity)
#' Performance on the sub-ontology Biological Process (BP):
options(MapMan2GO.performance.universe.annotations = GO.BP)
bb.annos.no.ref.anc.non.sprot.BP <- predictorPerformance(best.blast.pred, pa.gene.col = "query", 
    pa.anno.col = "GO", process.annos.funk = identity)
#' Performance on the sub-ontology Cellular Component (CC):
options(MapMan2GO.performance.universe.annotations = GO.CC)
bb.annos.no.ref.anc.non.sprot.CC <- predictorPerformance(best.blast.pred, pa.gene.col = "query", 
    pa.anno.col = "GO", process.annos.funk = identity)
#' Performance on the sub-ontology Molecular Function (MF):
options(MapMan2GO.performance.universe.annotations = GO.MF)
bb.annos.no.ref.anc.non.sprot.MF <- predictorPerformance(best.blast.pred, pa.gene.col = "query", 
    pa.anno.col = "GO", process.annos.funk = identity)


ipr.annos.non.sprot <- read.table(input.args[[3]], sep = "\t", header = FALSE, 
    comment.char = "", quote = "", na.strings = "", stringsAsFactors = FALSE)
#' Performance on the whole of the Gene Ontology:
options(MapMan2GO.performance.universe.annotations = GO.OBO$id)
ipr.annos.no.ref.anc.non.sprot <- predictorPerformance(ipr.annos.non.sprot, pa.gene.col = "V1", 
    pa.anno.col = "V3", process.annos.funk = identity)
#' Performance on the sub-ontology Biological Process (BP):
options(MapMan2GO.performance.universe.annotations = GO.BP)
ipr.annos.no.ref.anc.non.sprot.BP <- predictorPerformance(ipr.annos.non.sprot, 
    pa.gene.col = "V1", pa.anno.col = "V3", process.annos.funk = identity)
#' Performance on the sub-ontology Cellular Component (CC):
options(MapMan2GO.performance.universe.annotations = GO.CC)
ipr.annos.no.ref.anc.non.sprot.CC <- predictorPerformance(ipr.annos.non.sprot, 
    pa.gene.col = "V1", pa.anno.col = "V3", process.annos.funk = identity)
#' Performance on the sub-ontology Molecular Function (MF):
options(MapMan2GO.performance.universe.annotations = GO.MF)
ipr.annos.no.ref.anc.non.sprot.MF <- predictorPerformance(ipr.annos.non.sprot, 
    pa.gene.col = "V1", pa.anno.col = "V3", process.annos.funk = identity)



#' Now consider the universe of all possible GO term annotations to include
#' ancestral terms. Also extend both reference GO term annotations as well as
#' predicted GO term annotations with ancestral terms:
message("Using gold standard and predictions including ancestral GO terms (PlusAnc method).")
#' (Default process.annos.funk is MapMan2GO::addAncestors)

mercator.annos.non.sprot <- as.data.frame(readMercatorResultTable(input.args[[1]], 
    sanitize.accession = TRUE))
#' Performance on the whole of the Gene Ontology:
options(MapMan2GO.performance.universe.annotations = GO.OBO$id)
mercator.annos.perf.non.sprot <- predictorPerformance(mercator.annos.non.sprot, 
    pa.gene.col = "IDENTIFIER", pa.anno.col = "MapManBin.GO", process.predicted.annos.funk = splitMapManBinGOAs)
#' Performance on the sub-ontology Biological Process (BP):
options(MapMan2GO.performance.universe.annotations = GO.BP)
mercator.annos.perf.non.sprot.BP <- predictorPerformance(mercator.annos.non.sprot, 
    pa.gene.col = "IDENTIFIER", pa.anno.col = "MapManBin.GO", process.predicted.annos.funk = splitMapManBinGOAs)
#' Performance on the sub-ontology Cellular Component (CC):
options(MapMan2GO.performance.universe.annotations = GO.CC)
mercator.annos.perf.non.sprot.CC <- predictorPerformance(mercator.annos.non.sprot, 
    pa.gene.col = "IDENTIFIER", pa.anno.col = "MapManBin.GO", process.predicted.annos.funk = splitMapManBinGOAs)
#' Performance on the sub-ontology Molecular Function (MF):
options(MapMan2GO.performance.universe.annotations = GO.MF)
mercator.annos.perf.non.sprot.MF <- predictorPerformance(mercator.annos.non.sprot, 
    pa.gene.col = "IDENTIFIER", pa.anno.col = "MapManBin.GO", process.predicted.annos.funk = splitMapManBinGOAs)


#' Performance on the whole of the Gene Ontology:
options(MapMan2GO.performance.universe.annotations = GO.OBO$id)
ipr.annos.w.ref.anc.non.sprot <- predictorPerformance(ipr.annos.non.sprot, pa.gene.col = "V1", 
    pa.anno.col = "V3")
#' Performance on the sub-ontology Biological Process (BP):
options(MapMan2GO.performance.universe.annotations = GO.BP)
ipr.annos.w.ref.anc.non.sprot.BP <- predictorPerformance(ipr.annos.non.sprot, pa.gene.col = "V1", 
    pa.anno.col = "V3")
#' Performance on the sub-ontology Cellular Component (CC):
options(MapMan2GO.performance.universe.annotations = GO.CC)
ipr.annos.w.ref.anc.non.sprot.CC <- predictorPerformance(ipr.annos.non.sprot, pa.gene.col = "V1", 
    pa.anno.col = "V3")
#' Performance on the sub-ontology Molecular Function (MF):
options(MapMan2GO.performance.universe.annotations = GO.MF)
ipr.annos.w.ref.anc.non.sprot.MF <- predictorPerformance(ipr.annos.non.sprot, pa.gene.col = "V1", 
    pa.anno.col = "V3")

#' Performance on the whole of the Gene Ontology:
options(MapMan2GO.performance.universe.annotations = GO.OBO$id)
bb.annos.w.ref.anc.non.sprot <- predictorPerformance(best.blast.pred, pa.gene.col = "query", 
    pa.anno.col = "GO")
#' Performance on the sub-ontology Biological Process (BP):
options(MapMan2GO.performance.universe.annotations = GO.BP)
bb.annos.w.ref.anc.non.sprot.BP <- predictorPerformance(best.blast.pred, pa.gene.col = "query", 
    pa.anno.col = "GO")
#' Performance on the sub-ontology Cellular Component (CC):
options(MapMan2GO.performance.universe.annotations = GO.CC)
bb.annos.w.ref.anc.non.sprot.CC <- predictorPerformance(best.blast.pred, pa.gene.col = "query", 
    pa.anno.col = "GO")
#' Performance on the sub-ontology Molecular Function (MF):
options(MapMan2GO.performance.universe.annotations = GO.MF)
bb.annos.w.ref.anc.non.sprot.MF <- predictorPerformance(best.blast.pred, pa.gene.col = "query", 
    pa.anno.col = "GO")


#' Results in two flavours. Once including all gold standard proteins and
#' another time just those where Blast had no Hits with 100% sequence identity!
gold.standard.with.100.seq.sim.hits <- sort(unique(best.blast.tbl.non.sprot[which(best.blast.tbl.non.sprot$V3 == 
    100), "query.ukb.short.id"]))



#' Save results:
save(mercator.annos.non.sprot, mercator.annos.perf.non.sprot, mercator.annos.perf.non.sprot.BP, 
    mercator.annos.perf.non.sprot.CC, mercator.annos.perf.non.sprot.MF, ipr.annos.non.sprot, 
    ipr.annos.w.ref.anc.non.sprot, ipr.annos.w.ref.anc.non.sprot.BP, ipr.annos.w.ref.anc.non.sprot.CC, 
    ipr.annos.w.ref.anc.non.sprot.MF, ipr.annos.no.ref.anc.non.sprot, ipr.annos.no.ref.anc.non.sprot.BP, 
    ipr.annos.no.ref.anc.non.sprot.CC, ipr.annos.no.ref.anc.non.sprot.MF, best.blast.tbl.non.sprot, 
    bb.annos.no.ref.anc.non.sprot, bb.annos.no.ref.anc.non.sprot.BP, bb.annos.no.ref.anc.non.sprot.CC, 
    bb.annos.no.ref.anc.non.sprot.MF, bb.annos.w.ref.anc.non.sprot, bb.annos.w.ref.anc.non.sprot.BP, 
    bb.annos.w.ref.anc.non.sprot.CC, bb.annos.w.ref.anc.non.sprot.MF, gold.standard.with.100.seq.sim.hits, 
    file = file.path(input.args[[length(input.args)]], "data", "predictionPerformancesNonSwissProt.RData"))


message("DONE")

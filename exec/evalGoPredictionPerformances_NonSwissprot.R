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
#' Set the gold standard reference GO annotations and the columns in which to
#' lookup gold standard proteins' identifiers and GO annotations:
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
bb.annos.no.ref.anc.non.sprot <- predictorPerformance(best.blast.pred, pa.gene.col = "query", 
    pa.anno.col = "GO", process.annos.funk = identity)

ipr.annos.non.sprot <- read.table(input.args[[3]], sep = "\t", header = FALSE, 
    comment.char = "", quote = "", na.strings = "", stringsAsFactors = FALSE)
ipr.annos.no.ref.anc.non.sprot <- predictorPerformance(ipr.annos.non.sprot, pa.gene.col = "V1", 
    pa.anno.col = "V3", process.annos.funk = identity)



#' Now consider the universe of all possible GO term annotations to include
#' ancestral terms. Also extend both reference GO term annotations as well as
#' predicted GO term annotations with ancestral termsc

mercator.annos.non.sprot <- as.data.frame(readMercatorResultTable(input.args[[1]], 
    sanitize.accession = TRUE))
mercator.annos.perf.non.sprot <- predictorPerformance(mercator.annos.non.sprot, 
    pa.gene.col = "IDENTIFIER", pa.anno.col = "MapManBin.GO", process.predicted.annos.funk = splitMapManBinGOAs)

ipr.annos.w.ref.anc.non.sprot <- predictorPerformance(ipr.annos.non.sprot, pa.gene.col = "V1", 
    pa.anno.col = "V3")

bb.annos.w.ref.anc.non.sprot <- predictorPerformance(best.blast.pred, pa.gene.col = "query", 
    pa.anno.col = "GO")


#' Results in two flavours. Once including all gold standard proteins and
#' another time just those where Blast had no Hits with 100% sequence identity!
gold.standard.with.100.seq.sim.hits <- sort(unique(best.blast.tbl.non.sprot[which(best.blast.tbl.non.sprot$V3 == 
    100), "query.ukb.short.id"]))


#' Write out textual results:
capSumOut <- function(x) {
    capture.output(summary(x[, setdiff(colnames(x), "gene")]))
}
performance.summaries <- c("Mercator performance on Oryza Non SwissProt Non Elect Inferred GOA:", 
    capSumOut(mercator.annos.perf.non.sprot), "", "", "InterProScan performance on Oryza Non SwissProt Non Elect Inferred GOA:", 
    capSumOut(ipr.annos.w.ref.anc.non.sprot), "", "", "Best Blast performance on Oryza Non SwissProt Non Elect Inferred GOA:", 
    capSumOut(bb.annos.w.ref.anc.non.sprot), "", "", "Excluding those gold standard proteins for which Blast found Hits with 100 percent sequence similarity:", 
    "", "Mercator performance on Oryza Non SwissProt Non Elect Inferred GOA:", 
    capSumOut(mercator.annos.perf.non.sprot[which(!mercator.annos.perf.non.sprot$gene %in% 
        gold.standard.with.100.seq.sim.hits), ]), "", "", "InterProScan performance on Oryza Non SwissProt Non Elect Inferred GOA:", 
    capSumOut(ipr.annos.w.ref.anc.non.sprot[which(!ipr.annos.w.ref.anc.non.sprot$gene %in% 
        gold.standard.with.100.seq.sim.hits), ]), "", "", "Best Blast performance on Oryza Non SwissProt Non Elect Inferred GOA:", 
    capSumOut(bb.annos.w.ref.anc.non.sprot[which(!bb.annos.w.ref.anc.non.sprot$gene %in% 
        gold.standard.with.100.seq.sim.hits), ]))
writeLines(performance.summaries, file.path(input.args[[length(input.args)]], "inst", 
    "predictionPerformancesNonSwissProt_summaries.txt"))


#' Plot histograms of F1-Scores and Matthews correlation coefficient:
scores <- c("precision", "recall", "false.pos.rate", "f1.score", "mcc", "specificity")

#' - for mercator
for (gs.set in c("all", "non100PercSeqSimHit")) {
    for (score.i in scores) {
        plot.vals <- if (gs.set == "all") {
            mercator.annos.perf.non.sprot[, score.i]
        } else {
            mercator.annos.perf.non.sprot[which(mercator.annos.perf.non.sprot$gene %in% 
                gold.standard.with.100.seq.sim.hits), score.i]
        }
        pdf(file.path(input.args[[length(input.args)]], "inst", paste("mercator_", 
            score.i, "_", gs.set, "_non_SwissProt_Hist.pdf", sep = "")))
        plotDistAsHistAndBox(plot.vals, main = paste("Mercator", score.i, "distribution (", 
            gs.set, ")"), summary.as.title = TRUE)
        dev.off()
    }
}


#' - for InterProScan
for (gs.set in c("all", "non100PercSeqSimHit")) {
    for (score.i in scores) {
        plot.vals <- if (gs.set == "all") {
            ipr.annos.w.ref.anc.non.sprot[, score.i]
        } else {
            ipr.annos.w.ref.anc.non.sprot[which(ipr.annos.w.ref.anc.non.sprot$gene %in% 
                gold.standard.with.100.seq.sim.hits), score.i]
        }
        pdf(file.path(input.args[[length(input.args)]], "inst", paste("interProScanWithAncRefGoTerms_", 
            score.i, "_", gs.set, "_non_SwissProt_Hist.pdf", sep = "")))
        plotDistAsHistAndBox(plot.vals, main = paste("InterProScan", score.i, "distribution (", 
            gs.set, ")"), summary.as.title = TRUE)
        dev.off()
    }
}
for (gs.set in c("all", "non100PercSeqSimHit")) {
    for (score.i in scores) {
        plot.vals <- if (gs.set == "all") {
            ipr.annos.no.ref.anc.non.sprot[, score.i]
        } else {
            ipr.annos.no.ref.anc.non.sprot[which(ipr.annos.no.ref.anc.non.sprot$gene %in% 
                gold.standard.with.100.seq.sim.hits), score.i]
        }
        pdf(file.path(input.args[[length(input.args)]], "inst", paste("interProScanNoAncRefGoTerms_", 
            score.i, "_", gs.set, "_non_SwissProt_Hist.pdf", sep = "")))
        plotDistAsHistAndBox(plot.vals, main = paste("InterProScan", score.i, "distribution (", 
            gs.set, ")"), summary.as.title = TRUE)
        dev.off()
    }
}


#' - for BestBlast
for (gs.set in c("all", "non100PercSeqSimHit")) {
    for (score.i in scores) {
        plot.vals <- if (gs.set == "all") {
            bb.annos.w.ref.anc.non.sprot[, score.i]
        } else {
            bb.annos.w.ref.anc.non.sprot[which(bb.annos.w.ref.anc.non.sprot$gene %in% 
                gold.standard.with.100.seq.sim.hits), score.i]
        }
        pdf(file.path(input.args[[length(input.args)]], "inst", paste("bestBlastWithAncRefGoTerms_", 
            score.i, "_", gs.set, "_non_SwissProt_Hist.pdf", sep = "")))
        plotDistAsHistAndBox(plot.vals, main = paste("Best Blast", score.i, "distribution (", 
            gs.set, ")"), summary.as.title = TRUE)
        dev.off()
    }
}
for (gs.set in c("all", "non100PercSeqSimHit")) {
    for (score.i in scores) {
        plot.vals <- if (gs.set == "all") {
            bb.annos.no.ref.anc.non.sprot[, score.i]
        } else {
            bb.annos.no.ref.anc.non.sprot[which(bb.annos.no.ref.anc.non.sprot$gene %in% 
                gold.standard.with.100.seq.sim.hits), score.i]
        }
        pdf(file.path(input.args[[length(input.args)]], "inst", paste("bestBlastNoAncRefGoTerms_", 
            score.i, "_", gs.set, "_non_SwissProt_Hist.pdf", sep = "")))
        plotDistAsHistAndBox(plot.vals, main = paste("Best Blast", score.i, "distribution (", 
            gs.set, ")"), summary.as.title = TRUE)
        dev.off()
    }
}


#' Scatterplots of of n.truth againt false positives, true positives and n.pred:
scatter.y <- c("n.pred", "false.pos", "true.pos")

#' - for mercator
for (gs.set in c("all", "non100PercSeqSimHit")) {
    plot.df <- if (gs.set == "all") {
        mercator.annos.perf.non.sprot
    } else {
        mercator.annos.perf.non.sprot[which(!mercator.annos.perf.non.sprot$gene %in% 
            gold.standard.with.100.seq.sim.hits), ]
    }
    for (scatter.i in scatter.y) {
        pdf(file.path(input.args[[length(input.args)]], "inst", paste("mercator_", 
            scatter.i, "_", gs.set, "_non_SwissProt_Scatter.pdf", sep = "")))
        plot(plot.df$n.truth, plot.df[, scatter.i], xlab = "n.truth", ylab = scatter.i, 
            pch = 20, main = paste("Mercator", scatter.i, "distribution (", gs.set, 
                ")"))
        dev.off()
    }
}


#' - for Best Blast with ancestral GO Terms
for (gs.set in c("all", "non100PercSeqSimHit")) {
    plot.df <- if (gs.set == "all") {
        bb.annos.w.ref.anc.non.sprot
    } else {
        bb.annos.w.ref.anc.non.sprot[which(!bb.annos.w.ref.anc.non.sprot$gene %in% 
            gold.standard.with.100.seq.sim.hits), ]
    }
    for (scatter.i in scatter.y) {
        pdf(file.path(input.args[[length(input.args)]], "inst", paste("bestBlast_", 
            scatter.i, "_", gs.set, "_non_SwissProt_Scatter.pdf", sep = "")))
        plot(plot.df$n.truth, plot.df[, scatter.i], xlab = "n.truth", ylab = scatter.i, 
            pch = 20, main = paste("Best Blast", scatter.i, "distribution (", gs.set, 
                ")"))
        dev.off()
    }
}


#' - for InterProScan with ancestral GO Terms
for (gs.set in c("all", "non100PercSeqSimHit")) {
    plot.df <- if (gs.set == "all") {
        ipr.annos.w.ref.anc.non.sprot
    } else {
        ipr.annos.w.ref.anc.non.sprot[which(!ipr.annos.w.ref.anc.non.sprot$gene %in% 
            gold.standard.with.100.seq.sim.hits), ]
    }
    for (scatter.i in scatter.y) {
        pdf(file.path(input.args[[length(input.args)]], "inst", paste("interProScan_", 
            scatter.i, "_", gs.set, "_non_SwissProt_Scatter.pdf", sep = "")))
        plot(plot.df$n.truth, plot.df[, scatter.i], xlab = "n.truth", ylab = scatter.i, 
            pch = 20, main = paste("InterProScan", scatter.i, "distribution (", 
                gs.set, ")"))
        dev.off()
    }
}


#' Save results:
save(mercator.annos.non.sprot, mercator.annos.perf.non.sprot, ipr.annos.non.sprot, 
    ipr.annos.w.ref.anc.non.sprot, ipr.annos.no.ref.anc.non.sprot, best.blast.tbl.non.sprot, 
    bb.annos.no.ref.anc.non.sprot, bb.annos.w.ref.anc.non.sprot, file = file.path(input.args[[length(input.args)]], 
        "data", "predictionPerformancesNonSwissProt.RData"))


message("DONE")

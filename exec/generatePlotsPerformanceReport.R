require(MapMan2GO)

message("USAGE: Rscript path/2/MapMan2GO/exec/generatePlotsPerformanceReport.R path/2/MapMan2GO")

input.args <- commandArgs(trailingOnly=TRUE)

#' Plot histograms of Precision, Recall, False Positive Rate, Specificty, F1-Scores and Matthews correlation coefficient:
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

message("DONE")

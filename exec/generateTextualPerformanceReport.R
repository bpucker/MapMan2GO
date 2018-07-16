require(MapMan2GO)

message("USAGE: Rscript path/2/MapMan2GO/exec/generateTextualPerformanceReport.R path/2/MapMan2GO")

input.args <- commandArgs(trailingOnly = TRUE)


#' Write out textual results:
capSumOut <- function(x) {
    capture.output(summary(x[, setdiff(colnames(x), "gene")]))
}
performance.summaries <- c("Mercator performance on Oryza Non SwissProt Non Elect Inferred GOA:", 
    capSumOut(mercator.annos.perf.non.sprot), "", "", "InterProScan performance on Oryza Non SwissProt Non Elect Inferred GOA:", 
    capSumOut(ipr.annos.w.ref.anc.non.sprot), "", "", "Best Blast performance on Oryza Non SwissProt Non Elect Inferred GOA:", 
    capSumOut(bb.annos.w.ref.anc.non.sprot), "", "", "InterProScan performance on Oryza Non SwissProt Non Elect Inferred GOA excluding ancestral GO terms:", 
    capSumOut(ipr.annos.no.ref.anc.non.sprot), "", "", "Best Blast performance on Oryza Non SwissProt Non Elect Inferred GOA excluding ancestral GO terms:", 
    capSumOut(bb.annos.no.ref.anc.non.sprot), "", "", "Excluding those gold standard proteins for which Blast found Hits with 100 percent sequence similarity:", 
    "", "Mercator performance on Oryza Non SwissProt Non Elect Inferred GOA:", 
    capSumOut(mercator.annos.perf.non.sprot[which(!mercator.annos.perf.non.sprot$gene %in% 
        gold.standard.with.100.seq.sim.hits), ]), "", "", "InterProScan performance on Oryza Non SwissProt Non Elect Inferred GOA:", 
    capSumOut(ipr.annos.w.ref.anc.non.sprot[which(!ipr.annos.w.ref.anc.non.sprot$gene %in% 
        gold.standard.with.100.seq.sim.hits), ]), "", "", "Best Blast performance on Oryza Non SwissProt Non Elect Inferred GOA:", 
    capSumOut(bb.annos.w.ref.anc.non.sprot[which(!bb.annos.w.ref.anc.non.sprot$gene %in% 
        gold.standard.with.100.seq.sim.hits), ]), "", "", "InterProScan performance on Oryza Non SwissProt Non Elect Inferred GOA excluding ancestral GO terms:", 
    capSumOut(ipr.annos.no.ref.anc.non.sprot[which(!ipr.annos.no.ref.anc.non.sprot$gene %in% 
        gold.standard.with.100.seq.sim.hits), ]), "", "", "Best Blast performance on Oryza Non SwissProt Non Elect Inferred GOA excluding ancestral GO terms:", 
    capSumOut(bb.annos.no.ref.anc.non.sprot[which(!bb.annos.no.ref.anc.non.sprot$gene %in% 
        gold.standard.with.100.seq.sim.hits), ]), "", "", "", "Performance evaluation on the Biological Process category", 
    "", "", "Mercator performance on Oryza Non SwissProt Non Elect Inferred GOA for BP:", 
    capSumOut(mercator.annos.perf.non.sprot.BP), "", "", "InterProScan performance on Oryza Non SwissProt Non Elect Inferred GOA for BP:", 
    capSumOut(ipr.annos.w.ref.anc.non.sprot.BP), "", "", "Best Blast performance on Oryza Non SwissProt Non Elect Inferred GOA for BP:", 
    capSumOut(bb.annos.w.ref.anc.non.sprot.BP), "", "", "InterProScan performance on Oryza Non SwissProt Non Elect Inferred GOA for BP excluding ancestral GO terms:", 
    capSumOut(ipr.annos.no.ref.anc.non.sprot.BP), "", "", "Best Blast performance on Oryza Non SwissProt Non Elect Inferred GOA for BP excluding ancestral GO terms:", 
    capSumOut(bb.annos.no.ref.anc.non.sprot.BP), "", "", "Biological Process excluding those gold standard proteins for which Blast found Hits with 100 percent sequence similarity:", 
    "", "Mercator performance on Oryza Non SwissProt Non Elect Inferred GOA for BP:", 
    capSumOut(mercator.annos.perf.non.sprot.BP[which(!mercator.annos.perf.non.sprot.BP$gene %in% 
        gold.standard.with.100.seq.sim.hits), ]), "", "", "InterProScan performance on Oryza Non SwissProt Non Elect Inferred GOA for BP:", 
    capSumOut(ipr.annos.w.ref.anc.non.sprot.BP[which(!ipr.annos.w.ref.anc.non.sprot.BP$gene %in% 
        gold.standard.with.100.seq.sim.hits), ]), "", "", "Best Blast performance on Oryza Non SwissProt Non Elect Inferred GOA for BP:", 
    capSumOut(bb.annos.w.ref.anc.non.sprot.BP[which(!bb.annos.w.ref.anc.non.sprot.BP$gene %in% 
        gold.standard.with.100.seq.sim.hits), ]), "", "", "InterProScan performance on Oryza Non SwissProt Non Elect Inferred GOA for BP excluding ancestral GO terms:", 
    capSumOut(ipr.annos.no.ref.anc.non.sprot.BP[which(!ipr.annos.no.ref.anc.non.sprot.BP$gene %in% 
        gold.standard.with.100.seq.sim.hits), ]), "", "", "Best Blast performance on Oryza Non SwissProt Non Elect Inferred GOA for BP excluding ancestral GO terms:", 
    capSumOut(bb.annos.no.ref.anc.non.sprot.BP[which(!bb.annos.no.ref.anc.non.sprot.BP$gene %in% 
        gold.standard.with.100.seq.sim.hits), ]), "", "", "", "Performance evaluation on the Cellular Component category", 
    "", "", "Mercator performance on Oryza Non SwissProt Non Elect Inferred GOA for CC:", 
    capSumOut(mercator.annos.perf.non.sprot.CC), "", "", "InterProScan performance on Oryza Non SwissProt Non Elect Inferred GOA for CC:", 
    capSumOut(ipr.annos.w.ref.anc.non.sprot.CC), "", "", "Best Blast performance on Oryza Non SwissProt Non Elect Inferred GOA for CC:", 
    capSumOut(bb.annos.w.ref.anc.non.sprot.CC), "", "", "InterProScan performance on Oryza Non SwissProt Non Elect Inferred GOA for CC excluding ancestral GO terms:", 
    capSumOut(ipr.annos.no.ref.anc.non.sprot.CC), "", "", "Best Blast performance on Oryza Non SwissProt Non Elect Inferred GOA for CC excluding ancestral GO terms:", 
    capSumOut(bb.annos.no.ref.anc.non.sprot.CC), "", "", "Cellular Component excluding those gold standard proteins for which Blast found Hits with 100 percent sequence similarity:", 
    "", "Mercator performance on Oryza Non SwissProt Non Elect Inferred GOA for CC:", 
    capSumOut(mercator.annos.perf.non.sprot.CC[which(!mercator.annos.perf.non.sprot.CC$gene %in% 
        gold.standard.with.100.seq.sim.hits), ]), "", "", "InterProScan performance on Oryza Non SwissProt Non Elect Inferred GOA for CC:", 
    capSumOut(ipr.annos.w.ref.anc.non.sprot.CC[which(!ipr.annos.w.ref.anc.non.sprot.CC$gene %in% 
        gold.standard.with.100.seq.sim.hits), ]), "", "", "Best Blast performance on Oryza Non SwissProt Non Elect Inferred GOA for CC:", 
    capSumOut(bb.annos.w.ref.anc.non.sprot.CC[which(!bb.annos.w.ref.anc.non.sprot.CC$gene %in% 
        gold.standard.with.100.seq.sim.hits), ]), "", "", "InterProScan performance on Oryza Non SwissProt Non Elect Inferred GOA for CC excluding ancestral GO terms:", 
    capSumOut(ipr.annos.no.ref.anc.non.sprot.CC[which(!ipr.annos.no.ref.anc.non.sprot.CC$gene %in% 
        gold.standard.with.100.seq.sim.hits), ]), "", "", "Best Blast performance on Oryza Non SwissProt Non Elect Inferred GOA for CC excluding ancestral GO terms:", 
    capSumOut(bb.annos.no.ref.anc.non.sprot.CC[which(!bb.annos.no.ref.anc.non.sprot.CC$gene %in% 
        gold.standard.with.100.seq.sim.hits), ]), "", "", "", "Performance evaluation on the Molecular Function category", 
    "", "", "Mercator performance on Oryza Non SwissProt Non Elect Inferred GOA for MF:", 
    capSumOut(mercator.annos.perf.non.sprot.MF), "", "", "InterProScan performance on Oryza Non SwissProt Non Elect Inferred GOA for MF:", 
    capSumOut(ipr.annos.w.ref.anc.non.sprot.MF), "", "", "Best Blast performance on Oryza Non SwissProt Non Elect Inferred GOA for MF:", 
    capSumOut(bb.annos.w.ref.anc.non.sprot.MF), "", "", "InterProScan performance on Oryza Non SwissProt Non Elect Inferred GOA for MF excluding ancestral GO terms:", 
    capSumOut(ipr.annos.no.ref.anc.non.sprot.MF), "", "", "Best Blast performance on Oryza Non SwissProt Non Elect Inferred GOA for MF excluding ancestral GO terms:", 
    capSumOut(bb.annos.no.ref.anc.non.sprot.MF), "", "", "Molecular Functin excluding those gold standard proteins for which Blast found Hits with 100 percent sequence similarity:", 
    "", "Mercator performance on Oryza Non SwissProt Non Elect Inferred GOA for MF:", 
    capSumOut(mercator.annos.perf.non.sprot.MF[which(!mercator.annos.perf.non.sprot.MF$gene %in% 
        gold.standard.with.100.seq.sim.hits), ]), "", "", "InterProScan performance on Oryza Non SwissProt Non Elect Inferred GOA for MF:", 
    capSumOut(ipr.annos.w.ref.anc.non.sprot.MF[which(!ipr.annos.w.ref.anc.non.sprot.MF$gene %in% 
        gold.standard.with.100.seq.sim.hits), ]), "", "", "Best Blast performance on Oryza Non SwissProt Non Elect Inferred GOA for MF:", 
    capSumOut(bb.annos.w.ref.anc.non.sprot.MF[which(!bb.annos.w.ref.anc.non.sprot.MF$gene %in% 
        gold.standard.with.100.seq.sim.hits), ]), "", "", "InterProScan performance on Oryza Non SwissProt Non Elect Inferred GOA for MF excluding ancestral GO terms:", 
    capSumOut(ipr.annos.no.ref.anc.non.sprot.MF[which(!ipr.annos.no.ref.anc.non.sprot.MF$gene %in% 
        gold.standard.with.100.seq.sim.hits), ]), "", "", "Best Blast performance on Oryza Non SwissProt Non Elect Inferred GOA for MF excluding ancestral GO terms:", 
    capSumOut(bb.annos.no.ref.anc.non.sprot.MF[which(!bb.annos.no.ref.anc.non.sprot.MF$gene %in% 
        gold.standard.with.100.seq.sim.hits), ]) )
writeLines(performance.summaries, file.path(input.args[[length(input.args)]], "inst", 
    "predictionPerformancesNonSwissProt_summaries.txt"))


message("DONE")


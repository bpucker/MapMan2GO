require(MapMan2GO)

message("USAGE: Rscript path/2/MapMan2GO/exec/generateTextualPerformanceReport.R path/2/MapMan2GO")

input.args <- commandArgs(trailingOnly = TRUE)


#' Write out textual results:
capSumOut <- function(x) {
    capture.output(summary(x[, setdiff(colnames(x), "gene")]))
}
performance.summaries <- c("Mercator performance on Oryza Non SwissProt Non Elect Inferred GOA:", 
    capSumOut(mercator.annos.perf.non.sprot), "", "", "InterProScan performance on Oryza Non SwissProt Non Elect Inferred GOA:", 
    capSumOut(ipr.annos.no.ref.anc.non.sprot), "", "", "", "Performance evaluation on the Biological Process category", 
    "", "", "Mercator performance on Oryza Non SwissProt Non Elect Inferred GOA for BP:", 
    capSumOut(mercator.annos.perf.non.sprot.BP), "", "", "InterProScan performance on Oryza Non SwissProt Non Elect Inferred GOA for BP:", 
    capSumOut(ipr.annos.w.ref.anc.non.sprot.BP), "", "", "InterProScan performance on Oryza Non SwissProt Non Elect Inferred GOA for BP excluding ancestral GO terms:", 
    capSumOut(ipr.annos.no.ref.anc.non.sprot.BP),  "", "", "", "Performance evaluation on the Cellular Component category", 
    "", "", "Mercator performance on Oryza Non SwissProt Non Elect Inferred GOA for CC:", 
    capSumOut(mercator.annos.perf.non.sprot.CC), "", "", "InterProScan performance on Oryza Non SwissProt Non Elect Inferred GOA for CC:", 
    capSumOut(ipr.annos.w.ref.anc.non.sprot.CC), "", "", "InterProScan performance on Oryza Non SwissProt Non Elect Inferred GOA for CC excluding ancestral GO terms:", 
    capSumOut(ipr.annos.no.ref.anc.non.sprot.CC), "", "", "", "Performance evaluation on the Molecular Function category", 
    "", "", "Mercator performance on Oryza Non SwissProt Non Elect Inferred GOA for MF:", 
    capSumOut(mercator.annos.perf.non.sprot.MF), "", "", "InterProScan performance on Oryza Non SwissProt Non Elect Inferred GOA for MF:", 
    capSumOut(ipr.annos.w.ref.anc.non.sprot.MF), "", "", "InterProScan performance on Oryza Non SwissProt Non Elect Inferred GOA for MF excluding ancestral GO terms:", 
    capSumOut(ipr.annos.no.ref.anc.non.sprot.MF)  )
writeLines(performance.summaries, file.path(input.args[[length(input.args)]], "inst", 
    "predictionPerformancesNonSwissProt_summaries.txt"))


message("DONE")


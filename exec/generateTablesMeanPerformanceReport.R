require(MapMan2GO)

message("USAGE: Rscript path/2/MapMan2GO/exec/generateTableMeanPerformanceReport.R path/2/MapMan2GO")

input.args <- commandArgs(trailingOnly = TRUE)

#' Generate summary tables of the prediction analysis for the results: All GO terms, Biological Process, Cellular Component and Molecular Function
col.name <- c("Method", "TP", "TN", "FP", "FN")
row.name <- c("Mercator", "InterProScan")

#' All GO terms 
TP <- c(mean(mercator.annos.perf.non.sprot$recall, na.rm = T), mean(ipr.annos.w.ref.anc.non.sprot$recall, 
    na.rm = T))
TN <- c(mean(mercator.annos.perf.non.sprot$specificity, na.rm = T), mean(ipr.annos.w.ref.anc.non.sprot$specificity, 
    na.rm = T))
FP <- c(mean(mercator.annos.perf.non.sprot$false.pos.rate, na.rm = T), mean(ipr.annos.w.ref.anc.non.sprot$false.pos.rate, 
    na.rm = T))
FN <- c((1 - (mean(mercator.annos.perf.non.sprot$recall, na.rm = T))), (1 - (mean(ipr.annos.w.ref.anc.non.sprot$recall, 
    na.rm = T))))
performance.all.go <- data.frame(row.name, TP, TN, FP, FN)
colnames(performance.all.go) <- col.name
write.table(performance.all.go, file.path(input.args, "inst", "performanceEvalMeanAllGoTerms.txt"), 
    row.names = FALSE, quote = FALSE, sep = "\t")

#' Biological Process 
TP.BP <- c(mean(mercator.annos.perf.non.sprot.BP$recall, na.rm = T), mean(ipr.annos.w.ref.anc.non.sprot.BP$recall, 
    na.rm = T))
TN.BP <- c(mean(mercator.annos.perf.non.sprot.BP$specificity, na.rm = T), mean(ipr.annos.w.ref.anc.non.sprot.BP$specificity, 
    na.rm = T))
FP.BP <- c(mean(mercator.annos.perf.non.sprot.BP$false.pos.rate, na.rm = T), mean(ipr.annos.w.ref.anc.non.sprot.BP$false.pos.rate, 
    na.rm = T))
FN.BP <- c((1 - (mean(mercator.annos.perf.non.sprot.BP$recall, na.rm = T))), (1 - 
    (mean(ipr.annos.w.ref.anc.non.sprot.BP$recall, na.rm = T))))
performance.BP <- data.frame(row.name, TP.BP, TN.BP, FP.BP, FN.BP)
colnames(performance.BP) <- col.name
write.table(performance.BP, file.path(input.args, "inst", "performanceEvalMeanBP.txt"), 
    row.names = FALSE, quote = FALSE, sep = "\t")

#' Cellular Component
TP.CC <- c(mean(mercator.annos.perf.non.sprot.CC$recall, na.rm = T), mean(ipr.annos.w.ref.anc.non.sprot.CC$recall, 
    na.rm = T))
TN.CC <- c(mean(mercator.annos.perf.non.sprot.CC$specificity, na.rm = T), mean(ipr.annos.w.ref.anc.non.sprot.CC$specificity, 
    na.rm = T))
FP.CC <- c(mean(mercator.annos.perf.non.sprot.CC$false.pos.rate, na.rm = T), mean(ipr.annos.w.ref.anc.non.sprot.CC$false.pos.rate, 
    na.rm = T))
FN.CC <- c((1 - (mean(mercator.annos.perf.non.sprot.CC$recall, na.rm = T))), (1 - 
    (mean(ipr.annos.w.ref.anc.non.sprot.CC$recall, na.rm = T))))
performance.CC <- data.frame(row.name, TP.CC, TN.CC, FP.CC, FN.CC)
colnames(performance.CC) <- col.name
write.table(performance.CC, file.path(input.args, "inst", "performanceEvalMeanCC.txt"), 
    row.names = FALSE, quote = FALSE, sep = "\t")

#' Molecular Function
TP.MF <- c(mean(mercator.annos.perf.non.sprot.MF$recall, na.rm = T), mean(ipr.annos.w.ref.anc.non.sprot.MF$recall, 
    na.rm = T))
TN.MF <- c(mean(mercator.annos.perf.non.sprot.MF$specificity, na.rm = T), mean(ipr.annos.w.ref.anc.non.sprot.MF$specificity, 
    na.rm = T))
FP.MF <- c(mean(mercator.annos.perf.non.sprot.MF$false.pos.rate, na.rm = T), mean(ipr.annos.w.ref.anc.non.sprot.MF$false.pos.rate, 
    na.rm = T))
FN.MF <- c((1 - (mean(mercator.annos.perf.non.sprot.MF$recall, na.rm = T))), (1 - 
    (mean(ipr.annos.w.ref.anc.non.sprot.MF$recall, na.rm = T))))
performance.MF <- data.frame(row.name, TP.MF, TN.MF, FP.MF, FN.MF)
colnames(performance.MF) <- col.name
write.table(performance.MF, file.path(input.args, "inst", "performanceEvalMeanMF.txt"), 
    row.names = FALSE, quote = FALSE, sep = "\t")


message("DONE")

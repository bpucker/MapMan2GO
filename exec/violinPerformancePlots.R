require(MapMan2GO)

message("USAGE: RScript path/2/MapMan2Go/exec/violinPerformancePlots.R path/2/MapMan2GO")

f.n.r.m <- 1 - mercator.annos.f1$recal
f.n.r.bb <- 1 - bb.annos.w.ref.anc$recall
f.n.r.ipr <- 1 - ipr.annos.w.ref.anc$recall

pdf(file.path(input.args[[length(input.args)]], "inst", "truePositiveRate_violin.pdf"))
vioplot(mercator.annos.f1$recall[!is.na(mercator.annos.f1$recall)], bb.annos.w.ref.anc$recall[!is.na(bb.annos.w.ref.anc$recall)], 
    ipr.annos.w.ref.anc$recall[!is.na(ipr.annos.w.ref.anc$recall)], names = c("Mercator", 
        "Best Blast", "InterProScan"), col = "lightgrey")
title(main = "True Positive Rate")
dev.off()

pdf(file.path(input.args[[length(input.args)]], "inst", "falseNegativeRate_violin.pdf"))
vioplot(f.n.r.m[!is.na(f.n.r.m)], f.n.r.bb[!is.na(f.n.r.bb)], f.n.r.ipr[!is.na(f.n.r.ipr)], names = c("Mercator", "Best Blast", "InterProScan"), 
    col = "lightgrey")
title(main="False Negative Rate")
dev.off()

#' False positive rate and True negative rate can not be ploted because their quaintile values are the same, ploted as boxplot

colors <- brewer.pal(3, "Dark2") 
alpha.colors <- addAlpha( colors )

pdf(file.path(input.args[[length(input.args)]], "inst", "falsePositiveRate_boxplot.pdf"))
boxplot((mercator.annos.f1$false.pos.rate[!is.na(mercator.annos.f1$false.pos.rate)]), 
    (bb.annos.w.ref.anc$false.pos.rate[!is.na(bb.annos.w.ref.anc$false.pos.rate)]), 
    (ipr.annos.w.ref.anc$false.pos.rate[!is.na(ipr.annos.w.ref.anc$false.pos.rate)]), 
    names = c("Mercator", "Best Blast", "InterProScan"), col = alpha.colors, main = "False Positive Rate", pch = "-")
dev.off()

pdf(file.path(input.args[[length(input.args)]], "inst", "trueNegativeRate_boxplot.pdf"))
boxplot((mercator.annos.f1$specificity[!is.na(mercator.annos.f1$specificity)]), (bb.annos.w.ref.anc$specificity[!is.na(bb.annos.w.ref.anc$specificity)]), 
    (ipr.annos.w.ref.anc$specificity[!is.na(ipr.annos.w.ref.anc$specificity)]), names = c("Mercator", 
        "Best Blast", "InterProScan"), col = alpha.colors, main = "True Negative Rate", pch = "-")
dev.off()

message("DONE")

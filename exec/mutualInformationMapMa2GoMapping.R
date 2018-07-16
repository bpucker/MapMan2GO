require(MapMan2GO)
options(mc.cores = detectCores())

message("USAGE: Rscript path/2/MapMan2GO/exec/mutualInformationMapMan2GoMapping.R path/2/MapMan2GO")

input.args <- commandArgs(trailingOnly = TRUE)

cols <- lapply( mm.2.go.df$MapManBin.GO, function(mm.goa) as.integer(GO.OBO$id %in% splitMapManBinGOAs(mm.goa)) )
mm.2.go.mi.mtrx <- matrix( unlist(cols), ncol=nrow(mm.2.go.df), dimnames=list(GO.OBO$id, mm.2.go.df$MapManBin))
mm.2.go.mi <- mi.empirical(mm.2.go.mi.mtrx)
message("Mutual information between MapMan Bins and mapped GO Terms is: ", mm.2.go.mi)

#' Save results:
save(mm.2.go.mi.mtrx, mm.2.go.mi, file=file.path(input.args, "data", "mapMan2GoMutualInformation.RData"))

message("DONE")

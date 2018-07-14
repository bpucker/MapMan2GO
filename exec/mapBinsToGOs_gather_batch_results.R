require(MapMan2GO)

message("USAGE: Rscript path/2/MapMan2GO/exec/mapBinsToGOs_gather_batch_results.R path/2/batch_results_RData_dir path/2/MapMan2GO")

input.args <- commandArgs(trailingOnly = TRUE)

batch.rdata <- system(paste("ls ", input.args[[1]], "/*.RData", sep = ""), intern = TRUE)

load(batch.rdata[[1]])
for (b in batch.rdata[2:length(batch.rdata)]) {
    mm.2.full.desc.tmp <- mm.2.full.desc
    mm.2.go.tmp <- mm.2.go
    mm.2.go.df.tmp <- mm.2.go.df
    mm.desc.df.tmp <- mm.desc.df
    load(b)
    mm.2.full.desc <- rbind(mm.2.full.desc.tmp, mm.2.full.desc)
    mm.2.go <- append(mm.2.go.tmp, mm.2.go)
    mm.2.go.df <- rbind(mm.2.go.df.tmp, mm.2.go.df)
    mm.desc.df <- rbind(mm.desc.df.tmp, mm.desc.df)
}


#' Save results:
write.table(mm.2.full.desc, file.path(input.args[[length(input.args)]], "inst", 
    "MapManBins2GO.txt"), sep = "\t", row.names = FALSE)
save(mm.2.go, mm.2.go.df, mm.2.full.desc, mm.desc.df, file = file.path(input.args[[length(input.args)]], 
    "data", "MapManBins2GO.RData"))


message("DONE")

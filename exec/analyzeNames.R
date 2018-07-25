require(MapMan2GO)

message("USAGE: Rscript path/2/MapMan2GO/exec/analyzeNames.R path/2/MapMan2GO")
input.args <- commandArgs(trailingOnly = TRUE)

protID <- mercator.annos.non.sprot$IDENTIFIER[which(mercator.annos.non.sprot$IDENTIFIER != "")]
mercator.bin <- c()
name.mercator.bin <- c()
ipr.id <- c()
name.interpro <- c()
regexp.filter <- readLines(file.path(input.args[[length(input.args)]], "inst", "regular_expressions_to_filter.txt"))

for (i in 1:length(protID)) {
  name.mercator.bin[i] <- mercator.annos.non.sprot$NAME[ which(mercator.annos.non.sprot$IDENTIFIER == protID[i] ) ]
  mercator.bin[i] <- mercator.annos.non.sprot$BINCODE[ which(mercator.annos.non.sprot$IDENTIFIER == protID[i] ) ]
  name.mercator.bin[i] <- mercator.annos.non.sprot$NAME[ which(mercator.annos.non.sprot$IDENTIFIER == protID[i] ) ]
  if (length(which(ipr.annos.non.sprot$V1 == protID[i]))> 0 ) {
    ipr.id[i] <- unique(ipr.annos.non.sprot$V2[which(ipr.annos.non.sprot$V1 == protID[i])])
    name.interpro[i] <- ipr.db[[ ipr.id[i] ]]$NAME
# Need to fix. Some Swissprot proteins have more than one InterproID, i.e. protID[72]
  } else {
    ipr.id[i] <- NA
    name.interpro[i] <- NA
  }
}
name.mercator.bin.filtered <- filter(name.mercator.bin, regexp.filter)
name.interpro.filtered <- filter(name.interpro, regexp.filter)

protein.descriptions.df <- data.frame(protID, mercator.bin, name.mercator.bin.filtered, ipr.id, name.interpro.filtered)
names(protein.descriptions.df) <- c("Swissprot.ID", "MapManBin", "Bin.Description", "InterproID", "Interpro.Description")
 


require(MapMan2GO)

message("USAGE: Rscript path/2/MapMan2GO/exec/extractFromGOANonSprotNon_IEA_ISS_RCA_gene_ids.R gene_ontology_goa.txt uniprot_sprot.fasta path/2/MapMan2GO")

input.args <- commandArgs(trailingOnly = TRUE)


#' Load Gene Ontology GOA table:
go.goa <- read.table(input.args[[1]], sep = "\t", header = FALSE, quote = "", comment.char = "!", 
    na.strings = "", stringsAsFactors = FALSE)

#' Load UniprotKB SwissProt:
sprot <- read.fasta(input.args[[2]], seqtype = "AA", strip.desc = TRUE, as.string = TRUE)
sprot.short.ids <- sanitizeAccession(names(sprot))

#' Find genes IDs in above GOA table that have non electronically made GO
#' annotations and are not contained in SwissProt:
excl.evidence.codes <- c("IEA", "ISS", "RCA")
gene.ids <- sort(unique(go.goa[which(!go.goa$V7 %in% excl.evidence.codes & !go.goa$V2 %in% 
    sprot.short.ids), "V2"]))


#' Save results:
writeLines( gene.ids, file.path( input.args[[length(input.args)]], "inst", "nonSwissProtGenesWithNonElectMadeGOAs_IDs.txt" ) )


message("DONE")

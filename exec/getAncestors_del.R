require(MapMan2GO)
options(mc.cores = detectCores())
load("../data/MapManBins2GO.RData")

uniq.goas <- unique(mm.2.go.df$MapManBin.GO)
uniq.goas.counts <- vector()
uniq.goas.list <- unique(mm.2.go.df$MapManBin.GO)
for (i in 1:length(uniq.goas.list)){
  uniq.goas.counts[i] <- length(which(mm.2.go.df$MapManBin.GO  == uniq.goas.list[i]))
}
uniq.goas.list.df <- as.data.frame(uniq.goas.list)
uniq.goas.list.df$uniq.goas.counts <- uniq.goas.counts

#' Define addAncestors function and get ancestors for earch GO term
go.term.dist.to.root.funks <- list(MF = GOMFANCESTOR, BP = GOBPANCESTOR, CC = GOCCANCESTOR)
addAncestors <- function(go.terms, ancestors = getOption("MapMan2GO.term.ancestors", 
                                                         list(MF = GOMFANCESTOR, BP = GOBPANCESTOR, CC = GOCCANCESTOR)), exclude.root = TRUE, 
                         root.go = "all") {
  sort(unique(unlist(lapply(go.term, function(g.id) {
    g.t <- GOTERM[[g.id]]
    if (!is.null(g.t) && length(g.t) > 0) {
      g.name <- attr(g.t, "Term")
      g.ont <- attr(g.t, "Ontology")
      g.incl.anc <- c(g.id, get(g.id, go.term.dist.to.root.funks[[g.ont]]))
      if (exclude.root) 
        g.incl.anc <- setdiff(g.incl.anc, root.go)
      g.incl.anc
    } else {
      g.id
    }
  }))))
}

uniq.go.term <- unique(mm.2.full.desc$GO.Term)

ancestors.go <- list()
a <- list()
for (i in 1:length(uniq.go.term)){
  go.term <- uniq.go.term[i]
  a[i] <- paste(addAncestors(go.term), collapse = ",")
  ancestors.go  <- strsplit(as.character(a), ",")
}
names(ancestors.go)<-uniq.go.term[1:length(uniq.go.term)]

#' Merge ancestors.go with a list of all GOAs
uniq.goas.split<- list()
go.with.ancestors <- list()
for (i in 1:length(uniq.go.term)) {
  uniq.goas.split[i] <- strsplit(as.character(uniq.goas.list[i]), ",")
  go.with.ancestors[[i]] <- get(uniq.go.term[i], ancestors.go)
}

mm.2.go.df$MapManBin.GO <- goas.with.ancestors


#### Testing atm

duplicated.GOAs <- mm.2.go.df$MapManBin.GO[ which(duplicated(mm.2.go.df$MapManBin.GO))]
coords.duplicated.GOAs <- which(duplicated(mm.2.go.df$MapManBin.GO))
duplicated.GOAs.uniq <- unique(duplicated.GOAs)


# Look for each duplicated GOA in the data frame mm.2.go.df, must return the Bins
Bins.GOAs.ancestors <- list()
for(i in 1:length(duplicated.GOAs.ancestors.uniq)) {
  Bins.GOAs.ancestors[[i]] <- mm.2.go.df$MapManBin[grep(pattern=duplicated.GOAs.ancestors.uniq[i], mm.2.go.df$MapManBin.GO)]
}

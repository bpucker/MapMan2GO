.onLoad <- function(libname = find.package("MapMan2GO"), pkgname = "MapMan2GO") {
    data("GO_OBO", package = "MapMan2GO")
    data("MapManBinsVsSwissprotAndGOA", package = "MapMan2GO")
    data("MapManBins2GO", package = "MapMan2GO")
    data("UKB_Reference_GOA_InterPro", package = "MapMan2GO")
    data("predictionPerformances", package = "MapMan2GO")
    data("predictionPerformancesNonSwissProt", package = "MapMan2GO")
    data("EvidenceCodesOntology", package = "MapMan2GO")
    data("mapMan2GoMutualInformation", package = "MapMan2GO")
    data("interpro", package = "MapMan2GO")
}

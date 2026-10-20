#' @importFrom MetaboCoreUtils nameMapping
.onLoad <- function(libname, pkgname) {
    map <- .load_mapping()
    assign(".VARMAP", map, envir = asNamespace(pkgname))
    assign(".SV2CV", nameMapping(from = "spectra_variable", to = "cv_term",
                                 map = map), envir = asNamespace(pkgname))
    assign(".CV2SV", nameMapping(from = "cv_term", to = "spectra_variable",
                                 map = map), envir = asNamespace(pkgname))
}

#' Load the mapping for variable names from mzSpecLib to Spectra.
#'
#' @noRd
.load_mapping <- function() {
    map <- utils::read.table(system.file("txt", "mzSpecLib-mapping.txt",
                                         package = "MsBackendMzSpecLib"),
                             sep = "\t", header = TRUE)
    rownames(map) <- map[, 1L]
    map
}

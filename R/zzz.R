#' @importFrom MetaboCoreUtils nameMapping
.onLoad <- function(libname, pkgname) {
    map <- .load_mapping()
    assign(".VARMAP", map, envir = .get_envir(pkgname))
    assign(".SV2CV", nameMapping(from = "spectra_variable", to = "cv_term",
                                 map = map), envir = .get_envir(pkgname))
    assign(".CV2SV", nameMapping(from = "cv_term", to = "spectra_variable",
                                 map = map), envir = .get_envir(pkgname))
}

.get_envir <- function(x) {
    asNamespace(x)
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

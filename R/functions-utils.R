## Utility functions

#' @title Conversion between spectra variables and MS CV terms
#'
#' @description
#'
#' `spectraVariableFromCvTerm()` and `spectraVariableToCvTerm()` are helper
#' functions to *translate* between `Spectra` spectra variable names and MS
#' controlled vocabulary (CV) terms defined and used in the
#' [mzSpecLib](https://hupo-psi.github.io/mzSpecLib/specification/metadata.html)
#' format.
#'
#' The complete mapping of potential spectra variables to CV terms can be
#' retrieved with `spectraVariableMapping(MsBackendMzSpecLib())` (see examples).
#'
#' @param x `character()` with the names of spectra variables (for
#'     `spectraVariableToCvTerm()`) or CV terms (for
#'     `spectraVariableToCvTerm()`).
#'
#' @return named `character` with names being equal to `x` and values
#'     representing the mapped names/terms. `NA_character_` is returned for
#'     values in `x` that can not be mapped.
#'
#' @author Johannes Rainer
#'
#' @importFrom MetaboCoreUtils translate
#'
#' @examples
#'
#' ## Map spectra variables to CV terms used in the mzSpecLib definition
#' spectraVariableToCvTerm(c("rtime", "msLevel", "precursorMz"))
#'
#' ## For spectra variables without a matching value in the CV definition of
#' ## mzSpecLib an `NA` is returned:
#' spectraVariableToCvTerm(c("precScanNum"))
#'
#' ## Reverse mapping: from CV terms to spectra variables
#' spectraVariableFromCvTerm(c("MS:1003208", "MS:1000002"))
#'
#' ## Complete mapping from spectra variables to CV terms defined for
#' ## `MsBackendMzSpecLib` backends:
#' spectraVariableMapping(MsBackendMzSpecLib())
#'
#' @export
spectraVariableFromCvTerm <- function(x) {
    translate(x, mapping = .CV2SV, nomatch = NA_character_, warning = FALSE)
}

#' @rdname spectraVariableFromCvTerm
#'
#' @export
spectraVariableToCvTerm <- function(x) {
    translate(x, mapping = .SV2CV, nomatch = NA_character_, warning = FALSE)
}

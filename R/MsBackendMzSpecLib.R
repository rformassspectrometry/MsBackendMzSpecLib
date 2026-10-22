#' @title MsBackend for mzSpecLib files
#'
#' @name MsBackendMzSpecLib
#'
#' @description
#'
#' The `MsBackendMzSpecLib` [Spectra::MsBackend] backend provides an interface
#' to read and import mass spectrometry data (along with annotations) from data
#' files in the HUPO-PSI
#' [*mzSpecLib*](https://hupo-psi.github.io/mzSpecLib/index.html) format.
#' Through this backend, such files can then be directly used with a
#' [Spectra::Spectra] object.
#'
#' All spectra and analyte metadata information from the mzSpecLib file(s) are
#' available as [Spectra::spectraVariables()] in the `Spectra` object. The
#' mapping of the spectra attribute names (controlled vocabulary terms) to
#' spectra variable names can be listed by the
#' `spectraVariableMapping(MsBackendMzSpecLib())` function. See also the
#' [mzSpecLib documentation](https://hupo-psi.github.io/mzSpecLib/specification/metadata.html)
#' for more information on spectra attributes and the names and details for
#' the individual controlled vocabulary (CV) terms.
#'
#' @param object A `MsBackendMzSpecLib` object
#'
#' @param ... For `spectraVariableMapping()`: ignored.
#'
#' @author Johannes Rainer, Gabriele Tomè
#'
#' @examples
#'
#' ## List the mapping of spectra variables and mzSpecLib spectrum attributes
#' spectraVariableMapping(MsBackendMzSpecLib())
NULL

#' @importClassesFrom Spectra MsBackendCached
setClass(
    "MsBackendMzSpecLib",
    contains = "MsBackendCached",
    ## add slots
    prototype = prototype(
        readonly = TRUE, version = "0.1"))

#' @importMethodsFrom Spectra spectraVariableMapping
#'
#' @rdname MsBackendMzSpecLib
#'
#' @export
setMethod("spectraVariableMapping", "MsBackendMzSpecLib",
          function(object, ...) {
              .SV2CV
          })

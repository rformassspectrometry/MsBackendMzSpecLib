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
#' @section Available methods:
#'
#' The methods listed below are specifically implemented for
#' `MsBackendMzSpecLib`. The backend inherits and thus supports all other
#' methods from the [Spectra::MsBackend] class.
#'
#' - `backendInitialize()`: initializes `object` (the `MsBackendMzSpecLib`
#'   object) using mzSpecLib data files whose path is specified by `files`.
#'   This method is supposed to be called right after creating a
#'   `MsBackendMzSpecLib` object with `MsBackendMzSpecLib()` function.
#'
#' - `spectraVariables()`: returns a `character` vector with the spectra
#'   variables names of core spectra variables defined in the Spectra package
#'   and other additional variables contained in `object`.
#'   Note that also `"mz"` and `"intensity"` (which are by default not returned
#'   by the `spectraVariables,Spectra` method) are returned.
#'
#' @param BPPARAM Parameter object defining the parallel processing
#'   setup to import data in parallel. Defaults to `BPPARAM = bpparam()`.
#'   See [BiocParallel::bpparam()] for more information.
#'
#' @param files `character` specifying mzSpecLib file names.
#'
#' @param object A `MsBackendMzSpecLib` object
#'
#' @param ... For `spectraVariableMapping()`: ignored.
#'
#' @author Johannes Rainer, Gabriele Tomè
#'
#' @seealso The HUPO-PSI mzSpecLib specification,
#'   <https://www.psidev.info/mzspeclib>.
#'
#' @examples
#'
#' ## List the mapping of spectra variables and mzSpecLib spectrum attributes
#' spectraVariableMapping(MsBackendMzSpecLib())
#'
#' fl <- system.file("examples", "fetal_brain_tiny_mzspeclib.txt",
#'                  package = "MsBackendMzSpecLib")
#'
#' ## Create a MsBackend instance for that file
#' be <- backendInitialize(MsBackendMzSpecLib(), fl)
#' be
#'
NULL

#' @importClassesFrom Spectra MsBackendCached
setClass(
    "MsBackendMzSpecLib",
    contains = "MsBackendCached",
    slots = c(fileNames = "character",
              index = "data.frame",
              header = "data.frame"),
    prototype = prototype(
        fileNames = character(),
        index = data.frame(tag = character(), kind = character(),
                            attrset = logical(), key = character(),
                            start = numeric(), end = numeric(),
                            file = integer()),
        header = data.frame(group = integer(), accession = character(),
                            name = character(), value = character(),
                            file = integer()),
        readonly = TRUE,
        version = "0.1"))

#' @importMethodsFrom Spectra spectraVariableMapping
#'
#' @rdname MsBackendMzSpecLib
#'
#' @export
setMethod("spectraVariableMapping", "MsBackendMzSpecLib",
          function(object, ...) {
              .SV2CV
          })

#' @importFrom methods validObject
setValidity("MsBackendMzSpecLib", function(object) {
    msg <- NULL
    if (nrow(object@index) &&
        !all(object@index$file %in% seq_along(object@fileNames)))
        msg <- c(msg, "'index$file' must point into 'fileNames'")
    if (length(msg)) msg
    else TRUE
})

#' @rdname MsBackendMzSpecLib
#'
#' @importFrom methods callNextMethod
#'
#' @exportMethod backendInitialize
setMethod("backendInitialize", "MsBackendMzSpecLib",
          function(object, files, ..., BPPARAM = bpparam()) {
    if (missing(files))
        stop("Parameter 'files' is mandatory for 'MsBackendMzSpecLib'",
             call. = FALSE)
    if (!is.character(files))
        stop("Parameter 'files' is expected to be a character vector",
            " with the files names from where data should be",
            " imported", call. = FALSE)
    if (!length(files))
        return(object)
    files <- normalizePath(files, mustWork = FALSE)
    msg <- Spectra:::.valid_ms_backend_files_exist(files)
    if (length(msg))
        stop(msg, call. = FALSE)

    object <- .initialize(object, files, BPPARAM)

    validObject(object)
    object
})

#' @rdname MsBackendMzSpecLib
#'
#' @importFrom methods callNextMethod
#'
#' @exportMethod show
setMethod("show", "MsBackendMzSpecLib", function(object) {
    nm  <- .header_value(object@header, "MS:1003188") ## library name
    ver <- .header_value(object@header, "MS:1003190") ## library version
    lfv <- .header_value(object@header, "MS:1003186") ## library format version
    cat("MsBackendMzSpecLib with", length(object), "spectra\n")
    cat(" library format version: ", lfv, "\n")
    cat(" library name: ", ifelse(is.na(nm), "<unnamed>", nm), "\n")
    cat(ifelse(!is.na(ver), paste0(" library version: ", ver, "\n"), ""))
    cat(" files: ", object@fileNames, "\n")
})


#' @rdname MsBackendMzSpecLib
#'
#' @importFrom methods callNextMethod
#'
#' @exportMethod length
setMethod("length", "MsBackendMzSpecLib", function(x) {
    sum(x@index$kind == "Spectrum" & !x@index$attrset)
})

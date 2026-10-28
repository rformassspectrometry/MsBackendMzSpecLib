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
#' - `peaksData()`: gets the peak matrices of the spectra in the backend.
#'   Returns a `list` of `data.frame` with columns defined by parameter
#'   `columns` (which defaults to `columns = c("mz", "intensity")`. Use
#'   `peaksVariables()` to list all supported and available columns for a
#'   backend. The length of the `list` is equal to the number of spectra in
#'   `object`.
#'
#' - `peaksVariables()`: gets the supported peak variables (columns) for the
#'   backend.
#'
#' - `spectraData()`: gets spectra variables (specified by `columns`) from
#'   `object`.
#'
#' - `spectraVariables()`: returns a `character` vector with the spectra
#'   variables names of core spectra variables defined in the Spectra package
#'   and other additional variables contained in `object`.
#'   Note that also `"mz"` and `"intensity"` (which are by default not returned
#'   by the `spectraVariables,Spectra` method) are returned.
#'
#' @details
#' The `MsBackendMzSpecLib` is composed of three elements:
#' - `fileNames`: `character` vector with the associated files;
#' - `index`: `data.frame` with basic information for each section in the
#'   file. For each tag, it reports which tag it is, whether it is an
#'   AttributeSet tag, and its identifier key (if present), along with the
#'   start and end positions in the file and, if it is a spectrum, the start of
#'   the peak subsection;
#' - `header`: `data.frame` with the first parsed section of the file
#'   containing the file metadata.
#'
#'
#' @param BPPARAM Parameter object defining the parallel processing
#'   setup to import data in parallel. Defaults to `BPPARAM = bpparam()`.
#'   See [BiocParallel::bpparam()] for more information.
#'
#' @param columns For `spectraData()`: names of the spectra variables to extract
#'   from `object`. For `peaksData()`: names of the peak variables to extract.
#'   Defaults to `columns = c("mz", "intensity")`.
#'
#' @param files `character` specifying mzSpecLib file names.
#'
#' @param index `integer(1)` For `peaksVariables()`: index of the desired
#'     spectrum.
#'
#' @param object A `MsBackendMzSpecLib` object
#'
#' @param x A `MsBackendMzSpecLib` object
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
#' fl <- system.file("examples", "fetal_brain_tiny.mzspeclib.txt",
#'                  package = "MsBackendMzSpecLib")
#'
#' ## Create a MsBackend instance for that file
#' be <- backendInitialize(MsBackendMzSpecLib(), fl)
#' be
#'
#' spectraData <- spectraData(be)
#' head(spectraData)
#'
#' peaksData <- peaksData(be)
#' head(peaksData[[1]][[1]])
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

#' @importMethodsFrom Spectra spectraVariables
#'
#' @rdname MsBackendMzSpecLib
#'
#' @export
setMethod("spectraVariables", "MsBackendMzSpecLib",
          function(object, ...) {
              union(callNextMethod(object), .VARMAP$spectra_variable)
          })

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

#' @importMethodsFrom Spectra backendInitialize
#'
#' @rdname MsBackendMzSpecLib
#'
#' @importFrom methods callNextMethod
#'
#' @export
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
#' @importFrom methods show
setMethod("show", "MsBackendMzSpecLib", function(object) {
    ver <- .header_value(object@header, "MS:1003190")[1L] ## library version
    lfv <- .header_value(object@header, "MS:1003186")[1L] ## library format v
    cat("MsBackendMzSpecLib with", length(object), "spectra\n")
    cat(" library format version: ", lfv, "\n")
    cat(ifelse(!is.na(ver), paste0(" library version: ", ver, "\n"), ""))
    cat(" files: \n\t", paste(object@fileNames, collapse = "; \n\t"), "\n")
})

#' @rdname MsBackendMzSpecLib
#'
#' @importFrom methods callNextMethod
#'
#' @exportMethod length
setMethod("length", "MsBackendMzSpecLib", function(x) {
    sum(x@index$kind == "Spectrum" & !x@index$attrset)
})

#' @importMethodsFrom Spectra spectraData
#'
#' @rdname MsBackendMzSpecLib
#'
#' @importFrom Spectra rbindlistWithRownames
#'
#' @importFrom BiocParallel bplapply bpparam
#'
#' @importFrom IRanges NumericList
#'
#' @export
setMethod("spectraData", "MsBackendMzSpecLib",
        function(object, columns = spectraVariables(object),
                 BPPARAM = bpparam()) {
        idx_tbl <- object@index
        sps_list <- bplapply(seq_along(object@fileNames), function (f_i) {
            idx <- idx_tbl[idx_tbl$file == f_i, ]
            ## Get global AttributeSet = Spectrum
            set_rows <- which(idx$attrset & idx$kind == "Spectrum")
            metadata_all_spectrum <- setNames(lapply(set_rows, function(j) {
                .parse_attributes(.read_block(object@fileNames[f_i],
                                            idx$start[j], idx$end[j]))
                }), idx$key[set_rows])

            ## Get Metadata from single spectra
            rows <- which(!idx$attrset & idx$kind == "Spectrum")
            spectra_Data <- rbindlistWithRownames(
                lapply(rows, function(j) {
                    spectra_Data <- .spectrum_variables(as.integer(idx$key[j]),
                                .read_block(object@fileNames[f_i], idx$start[j],
                                            idx$peak_line[j]),
                                        metadata_all_spectrum)
                    ## Read peaks table
                    if (length(columns) != length(spectraVariables(object)) &&
                        any(.BASE_PEAK_COLS %in% columns)) {
                        p_tbl <- .peaks_table(idx$key[j], object@fileNames[f_i],
                                                idx, .BASE_PEAK_COLS)
                        if ("mz" %in% columns)
                            spectra_Data$mz <- list(p_tbl$mz)
                        if ("intensity" %in% columns)
                            spectra_Data$intensity <- list(p_tbl$intensity)
                    }
                    spectra_Data
                }), fill = TRUE)
            spectra_Data$dataOrigin <- basename(object@fileNames[f_i])
            spectra_Data$dataStorage <- object@fileNames[f_i]
            ## TODO: export optional columns by default?
            cols <- intersect(columns, names(spectra_Data))
            spectra_Data[cols]
        }, BPPARAM = BPPARAM)
        rbindlistWithRownames(sps_list, fill = TRUE)
})

#' @importMethodsFrom Spectra peaksData
#'
#' @rdname MsBackendMzSpecLib
#'
#' @importFrom Spectra rbindlistWithRownames
#'
#' @importFrom BiocParallel bplapply bpparam
#'
#' @export
setMethod("peaksData", "MsBackendMzSpecLib",
        function(object, columns = c("mz", "intensity"),
                 BPPARAM = bpparam()) {
            if (!all(columns %in% peaksVariables(object)))
                stop("Some of the requested columns are not available. Use ",
                    "peaksVariables() to list the supported columns.")
            idx_tbl <- object@index
            sps_list <- bplapply(seq_along(object@fileNames), function (f_i) {
                idx <- idx_tbl[idx_tbl$file == f_i, ]
                spectrum <- which(!idx$attrset & idx$kind == "Spectrum")
                peaks_data <- bplapply(spectrum, function(s) {
                    p_tbl <- .peaks_table(idx$key[s],
                                        object@fileNames[f_i], idx,
                                        peaksVariables(object,
                                                    as.integer(idx$key[s])))
                    p_tbl[make.names(intersect(columns,
                                        peaksVariables(object,
                                                    as.integer(idx$key[s]))))]
                }, BPPARAM = BPPARAM)
                peaks_data
            })
            sps_list
        })

#' @importMethodsFrom Spectra peaksVariables
#'
#' @rdname MsBackendMzSpecLib
#'
#' @importFrom stats na.omit
#'
#' @export
setMethod("peaksVariables", "MsBackendMzSpecLib",
        function(object, index = NULL, ...) {
            opt_columns <- spectraData(object, "peak_attribute")
            if (nrow(opt_columns)){
                if (!is.null(index))
                    opt_columns <- if (all(!is.na(opt_columns[index, ])))
                                        opt_columns[index, ]
                                    else .BASE_PEAK_COLS
                else
                    opt_columns <- na.omit(opt_columns)
                unique(unlist(lapply(opt_columns, function(o){
                    c(.BASE_PEAK_COLS, strsplit(o, "&&"))
                })))
            } else{
                .BASE_PEAK_COLS
            }
        })

#' @importFrom methods new
#'
#' @rdname MsBackendMzSpecLib
#'
#' @export
MsBackendMzSpecLib <- function() {
    new("MsBackendMzSpecLib")
}
#' @title Initialize MsBackendMzSpecLib object
#'
#' @description
#'
#' Initialize the `MsBackendMzSpecLib` object `x`:
#' - iterate over the files, for each file:
#' - get the *frame* information and build the `@indices` matrix with the
#'   frame index, scan index and file index columns.
#' - fix polarity information
#'
#' @param x `MsBackendMzSpecLib` object.
#'
#' @param file `character` specifying mzSpecLib file names.
#'
#' @param BPPARAM Parameter object defining the parallel processing
#'   setup to import data in parallel. Defaults to `BPPARAM = bpparam()`.
#'   See [BiocParallel::bpparam()] for more information.
#'
#' @author Gabriele Tomè
#'
#' @importFrom BiocParallel bplapply
#'
#' @importFrom BiocParallel bpparam
#'
#' @importFrom stats setNames
#'
#' @return initialized `MsBackendMzSpecLib` object
#'
#' @noRd
.initialize <- function(x, file = character(), BPPARAM = bpparam()){
    L <- bplapply(seq_len(length(file)), .open_mzspeclib,
                    file = file, BPPARAM = BPPARAM)

    x@fileNames <- file
    x@index <- do.call(rbind, lapply(L, "[[", 1))
    x@header <- do.call(rbind, lapply(L, "[[", 2))
    x
}

#' Import one file: tag index, header and spectra variables
#' @noRd
.open_mzspeclib <- function(i, file) {
    idx <- .index_sections(file[i])
    if (idx[1, "kind"] != "mzSpecLib")
        stop("Mandatory <mzSpecLib> section not found in ", file[i])
    idx$file <- i

    header <- .parse_attributes(.read_block(file[i], idx$start[1], idx$end[1]))

    list(index = idx[, c("file", "tag", "kind", "attrset",
                         "key", "start", "end", "peak_line")],
         header = header)
}

#' Index the section tags of an mzSpecLib file
#'
#' Streams the file in chunks (constant memory) and records the line number
#' of every section tag, without parsing the content.
#'
#' @param path `character(1)` path to the mzSpecLib file.
#'
#' @param chunk_bytes `numeric(1)` bytes read per chunk.
#'
#' @return `data.frame` with one row per tag and columns `tag`, `kind`,
#'     `attrset`, `key`, `start`, `end`, `peak_line`.
#'
#' @importFrom readr read_lines_chunked DataFrameCallback
#'
#' @noRd
.index_sections <- function(path, chunk_lines = 1e6) {
    tag_pattern <- paste0("^(?:", paste(.SECTION, collapse = "|"),
                            "|<Peaks>)\\s*$")
    idx  <- read_lines_chunked(path,
                DataFrameCallback$new(function(x, pos) {
                    i <- grepl(tag_pattern, x, perl = TRUE)
                    data.frame(tag = trimws(x[i]), line = pos + which(i) - 1L)
                }),
                chunk_size = chunk_lines, progress = FALSE)
    if (!nrow(idx))
        stop("No section tags found in the file.", call. = FALSE)
    peaks_pos <- idx[idx$tag == "<Peaks>",]
    idx <- idx[idx$tag != "<Peaks>",]

    info <- .parse_attributeset_tag(idx$tag)
    df <- data.frame(tag = idx$tag, kind = info$kind, attrset = info$attrset,
                key = info$key, start = idx$line,
                end = c(idx$line[-1L] - 1L, Inf), peak_line = NA)
    df[df$kind == "Spectrum" & !df$attrset, ]$peak_line <- peaks_pos$line
    df
}

#' Read the body of one section from the file
#'
#' @param path `character` with the path of the file to read.
#'
#' @param start `numeric` line where the section tag starts.
#'
#' @param end `numeric` line where the section ends.
#'
#' @return `character` with the non-empty lines of the section, without
#'     the tag line.
#'
#' @importFrom readr read_lines
#'
#' @noRd
.read_block <- function(path, start, end) {
    if (end-start-1 > 0) {
        body <- read_lines(path, skip = start, n_max = end-start-1,
                            progress = FALSE)
        body[nzchar(body)]
    } else
        NULL
}

#' Parse the CV terms of a section. Lines that don't match the CV term format
#' are dropped.
#'
#' @param body `character` with the lines of a section.
#'
#' @return `data.frame` with columns `group` (`integer`), `accession`,
#'     `name` and `value` (`character`).
#'
#' @noRd
.parse_attributes <- function(body) {
    if (!length(body)) {
        df <- data.frame(group = integer(), accession = character(),
                        name = character(), value = character())
    } else {
        m <- str_match(body, .TERM_REGEX)
        parsed <- !is.na(m[, 1])
        df <- data.frame(group = as.integer(m[parsed, 2]),
                        accession = m[parsed, 3],
                        name = gsub('^"|"$', "", m[parsed, 4]),
                        value = m[parsed, 5])
    }
    df
}

#' Parse a section tag to get its "type" and "key"
#'
#' @param tag `character` with tags.
#'
#' @return `list` with 3 `vectors`:
#'     - `attrset`: `logical` if the tag is an "AttributeSet",
#'     - `kind`: `character` with the section type (e.g. "Spectrum"),
#'     - `key`: `character` with the section key (`NA` if absent).
#'
#' @importFrom stringr str_match
#'
#' @noRd
.parse_attributeset_tag <- function(tag) {
    m <- str_match(tag, "^<(AttributeSet )?(\\w+)(?:=(\\w+))?>$")
    list(attrset = !is.na(m[, 2]), kind = m[, 3], key = m[, 4])
}

#' Get the value of a CV term from a header
#'
#' @param header `data.frame` as returned by `.parse_attributes()`.
#'
#' @param accession `character(1)` CV accession, e.g. `"MS:1003188"`.
#'
#' @return `character(1)` with the first matching value, or `NA` if absent.
#'
#' @noRd
.header_value <- function(header, accession) {
    v <- header$value[header$accession == accession]
    if (length(v)) v[1L] else NA_character_
}

## GLOBAL VARIABLES ------------------------------------------------------------

#' @title Regex to parse the the files
#'
#' @description
#' Piece of the regex:
#' 1. Start of the line,
#' 2. Group designator (optional),
#' 3. URI  identifier of the subject CV term (ex. `"MS:1003186"`),
#' 4. `"|"` to separate URI identifier of the CV term from the name of the CV,
#' 5. Name of the CV term (ex. `"library format version"`),
#' 6. `"="` separator between CV term and it's value,
#' 7. `"(.*)"`: value corresponding to CV term,
#' 8. End of the line.
#'
#' @noRd
.TERM_REGEX <- paste0(
    "^",
    "(?:\\[(\\d+)\\])?",
    "([A-Za-z][A-Za-z0-9_.-]*:[A-Za-z0-9_.-]+)",
    "\\|",
    "(\"[^\"]*\"|[^\"=][^=]*)",
    "=",
    "(.*)",
    "$"
)

.SECTION <- c("<mzSpecLib>",
            "<AttributeSet Spectrum=\\w+>",
            "<AttributeSet Analyte=\\w+>",
            "<AttributeSet Interpretation=\\w+>",
            "<AttributeSet Cluster=\\w+>",
            "<Cluster=\\d+>",
            "<Spectrum=\\d+>")

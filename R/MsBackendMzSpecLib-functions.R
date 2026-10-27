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
    header$file <- i

    list(index = idx[, c("file", "tag", "kind", "attrset",
                         "key", "start", "end", "peak_line")],
         header = header[, c("file", "group", "accession", "name", "value")])
}

#' Index the section tags of an mzSpecLib file
#'
#' Streams the file in chunks (constant memory) and records the line number
#' of every section tag, without parsing the content.
#'
#' @param path `character(1)` path to the mzSpecLib file.
#'
#' @param chunk_lines `numeric(1)` bytes read per chunk.
#'
#' @return `data.frame` with one row per tag and columns:
#'     - `tag`: Original tags detected in the file (e.g.,
#'           `"<AttributeSet Spectrum=1>"`).
#'     - `kind`: The type of tag detected (e.g., `"Spectrum"`).
#'     - `attrset`: `logical` indicating whether the field is an
#'           `"AttributeSet"`.
#'     - `key`: Key extracted from the tag, or `NA` if absent (e.g., `"1"`).
#'     - `start`: Line position where the section starts.
#'     - `end`: Line position where the section ends.
#'     - `peak_line`: Line position where the peak section starts if it is a
#'           spectrum section; `NA` otherwise.
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
                end = c(idx$line[-1L], Inf), peak_line = NA)
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
#' @return `data.frame` with columns:
#'     - `group`:  `integer` with the ID of the group; if not present
#'           `NA_integer_`,
#'     - `accession`: ID of the CV term (e.g. `"MS:1003188"`),
#'     - `name`: name of the CV term (e.g. `"library name"`)
#'     - `value`: `character` with the value of the CV term.
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
    if (length(v)) v else NA_character_
}

#' Spectra variables of one spectrum block (peaks are not parsed)
#'
#' @param scanIndex `integer(1)` with the scanIndex of the Spectrum.
#'
#' @param lines `character` lines of the spectrum block.
#'
#' @param mtd_all_spectrum `list` of general `<AttributeSet Spectrum=..>`
#'     tables.
#'
#' @return named `list` of scalars, one per spectra variable.
#'
#' @noRd
.spectrum_variables <- function(scanIndex, lines, mtd_all_spectrum) {
    mtd <- .merge_spectrum_mtd_section(scanIndex, lines, mtd_all_spectrum)
    ## Aggregate optional columns of peaks table in a single row
    if(length(val <- .header_value(mtd, "MS:1003254"))){
        mtd <- rbind(mtd[mtd$accession != "MS:1003254", ],
                    data.frame(group = NA, accession = "MS:1003254",
                        name = "peak attribute",
                        value = paste(val, collapse = "&&"),
                        spectra_variable = "peak_attribute", repeatable = TRUE))
    }

    spectra_Data <- as.data.frame(t(setNames(mtd$value, mtd$spectra_variable)))

    spectra_Data$scanIndex <- scanIndex
    ## Convertion to correct type
    cols <- intersect(names(spectra_Data), .VARMAP$spectra_variable)
    spectra_Data[cols] <- Map(convert_value, spectra_Data[cols],
                    .VARMAP[match(cols, .VARMAP$spectra_variable), "value"])
    ## Convert variable for Spectra
    if("rtime" %in% names(spectra_Data))
        spectra_Data <- .rtime(spectra_Data)
    if("polarity" %in% names(spectra_Data))
        spectra_Data$polarity <- .polarity(spectra_Data$polarity)
    if("spectrum_representation" %in% names(spectra_Data))
        spectra_Data$spectrum_representation <-
                .spectrum_representation(spectra_Data$spectrum_representation)
    if(all(c("isolationWindowTargetMz", "isolationWindowLowerOffset") %in%
                names(spectra_Data)))
        spectra_Data$isolationWindowLowerMz <-
                spectra_Data$isolationWindowTargetMz -
                    spectra_Data$isolationWindowLowerOffset
    if(all(c("isolationWindowTargetMz", "isolationWindowUpperOffset") %in%
                names(spectra_Data)))
        spectra_Data$isolationWindowUpperMz <-
                spectra_Data$isolationWindowTargetMz +
                    spectra_Data$isolationWindowUpperOffset

    # DataStorage, dataorigin == path file
    # Mz = NumericList
    # Intensity = NumericList
    spectra_Data
}

#' Merge global spectrum metadata with spectrum specific metadata
#'
#' @param scanIndex `integer(1)` with the scanIndex of the Spectrum.
#'
#' @param lines `character` lines of the spectrum block.
#'
#' @param mtd_all_spectrum `list` of general `<AttributeSet Spectrum=..>`
#'     tables.
#'
#' @return `data.frame` with merged and filtered spectrum metadata.
#'
#' @importFrom dplyr left_join
#'
#' @noRd
.merge_spectrum_mtd_section <- function(scanIndex, lines, mtd_all_spectrum){
    parts <- .split_by_tags(lines, .SPECTRUM_SECTION)
    mtd_spectrum <- .parse_attributes(parts[["Metadata"]])
    mtd_spectrum$spectra_variable <-
                spectraVariableFromCvTerm(mtd_spectrum$accession)

    ## Precedence: own terms > AttributeSet from last to first
    ## Attribute set name: "MS:1003212"
    ## Need to merge by spectrum Id and attribute set name
    attr_set_names <- .header_value(mtd_spectrum, "MS:1003212")
    ## Reverse order since last is more important than the first
    if (any(lapply(mtd_all_spectrum, nrow) != 0)) {
        mtd_all_spectrum <- rev(mtd_all_spectrum)
        mtd_global <- do.call(rbind, mtd_all_spectrum[grepl(
                        paste("all", scanIndex, attr_set_names, sep = "|",
                                collapse = "|"),
                        names(mtd_all_spectrum))])
        mtd_global$spectra_variable <-
                    spectraVariableFromCvTerm(mtd_global$accession)
        mtd_spectrum <- rbind(mtd_spectrum, mtd_global)
    }
    ## Parse other attribute: MS:1003275, MS:1003276
    if (length(.header_value(mtd_spectrum, "MS:1003275"))){
        mtd_spectrum <- .parse_other_attributes(mtd_spectrum)
    }
    ## Filter not allowed repeated fields
    mtd <- left_join(mtd_spectrum,
                    .VARMAP[, c("spectra_variable", "repeatable")],
                    by = "spectra_variable")
    mtd <- mtd[!duplicated(mtd$spectra_variable) |
        mtd$repeatable | is.na(mtd$repeatable), ]
    ## Update names Unit columns
    units <- which(mtd$accession == "UO:0000000")
    mtd$spectra_variable[units] <- paste(mtd$spectra_variable[units-1], "Unit",
                                        sep = "")
    mtd[is.na(mtd$spectra_variable), "spectra_variable"] <-
            mtd[is.na(mtd$spectra_variable), "name"]
    mtd
}


#' @title Split a section in its subpart
#' @noRd
.split_by_tags <- function(section, patterns) {
    tag_pattern <- paste0("^(?:", paste(patterns, collapse = "|"), ")$")
    is_tag <- grepl(tag_pattern, section, perl = TRUE)
    group  <- factor(cumsum(is_tag), levels = 0:sum(is_tag),
                     labels = c("Metadata", section[is_tag]))
    split(section[!is_tag], group[!is_tag])
}

#' "MS:1000130|positive scan" -> 1L
#' "MS:1000129|negative scan" -> 0L
#' @noRd
.polarity <- function(v) {
    if (!is.na(v)) {
        if (grep("MS:1000130", v))
            return(1L)
        else if (grep("MS:1000129", v))
            return(0L)
    }
    return(-1L)
}

#' "MS:1000127|centroid spectrum" -> TRUE
#' "MS:1000128|profile spectrum" -> FALSE
#' @noRd
.spectrum_representation <- function(v){
    if (!is.na(v)) {
        if (grep("MS:1000127", v))
            return(TRUE)
        else if (grep("MS:1000128", v))
            return(FALSE)
    }
    return(NA)
}

#' Retention time in seconds, converting based on the unit term of the same
#' group
#'
#' @noRd
.rtime <- function(spectra_Data) {
    rt_pos <- which(names(spectra_Data) == "rtime")
    if (!is.na(rt_pos)) {
        unit <- ifelse(names(spectra_Data)[rt_pos+1] == "rtimeUnit",
                        spectra_Data["rtimeUnit"], "UO:0000010|second")
        factor <- switch(sub("\\|.*$", "", unit),
                        "UO:0000031" = 60,
                        "UO:0000032" = 3600,
                        1)
        spectra_Data["rtime"] <- spectra_Data["rtime"] * factor
        spectra_Data["rtimeUnit"] <- NULL
    }
    spectra_Data
}

#' Parse `other attribute` fields and merging names of the field (`MS:1003275`)
#' and its value (`MS:1003276`)
#'
#' @return Update `data.frame` without the rows `MS:1003275` and `MS:1003276`
#'     and new rows based on their values.
#'
#' @noRd
.parse_other_attributes <- function(mtd) {
    names <- .header_value(mtd, "MS:1003275")
    values <- .header_value(mtd, "MS:1003276")
    rbind(mtd[!(mtd$accession %in% c("MS:1003275", "MS:1003276")), ],
          data.frame(group = NA, accession = "MS:0", name = names,
                     value = values, spectra_variable = names))
}

#' Convert a value to a given type
#'
#' @param x Value to convert.
#'
#' @param type `character(1)`, one of `"integer"` or `"numeric"`. Any other
#'     value returns `x` unchanged.
#'
#' @return `x` converted to the requested type.
#'
#' @noRd
convert_value <- function(x, type) {
    switch(type,
        integer = as.integer(round(as.numeric(x))),
        numeric = as.numeric(x),
        x
    )
}

#' Read the peaks of a spectrum as a data.frame
#'
#' @param s `character(1)` index of the spectrum.
#'
#' @param fileNames `character(1)` path to the file.
#'
#' @param idx `data.frame` with the `peak_line` and `end` line positions of
#'     each peak table.
#'
#' @param columns `character` column names of the peaks table.
#'
#' @return `data.frame` with the peaks, where `mz` and `intensity` are numeric.
#'
#' @noRd
.peaks_table <- function(s, fileNames, idx, columns) {
    ## Suppress warning since the spectrum can have different number of columns
    ## Read empty columns to merge later
    spectra_idx <- idx[which(idx$key == s), ]
    peaks_data <- suppressWarnings(read.delim2(text = .read_block(fileNames,
                        spectra_idx$peak_line, spectra_idx$end),
                        sep = "\t", col.names = columns))
    peaks_data$mz <- as.numeric(peaks_data$mz)
    peaks_data$intensity <- as.numeric(peaks_data$intensity)
    peaks_data
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

#' Tags of mzSpecLib file
#' @noRd
.SECTION <- c("<mzSpecLib>",
            "<AttributeSet Spectrum=\\w+>",
            "<AttributeSet Analyte=\\w+>",
            "<AttributeSet Interpretation=\\w+>",
            "<AttributeSet Cluster=\\w+>",
            "<Cluster=\\d+>",
            "<Spectrum=\\d+>")

#' Tags of a Spectrum section
#' @noRd
.SPECTRUM_SECTION <- c("<AttributeSet Spectrum=\\w+>",
                    "<Analyte=\\d+>",
                    "<Interpretation=\\d+>",
                    "<InterpretationMember=\\d+>")

.BASE_PEAK_COLS <- c("mz", "intensity", "annotation")

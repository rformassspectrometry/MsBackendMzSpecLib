test_that("MsBackendMzSpecLib works", {
    res <- MsBackendMzSpecLib()
    expect_s4_class(res, "MsBackendMzSpecLib")
})

test_that(".initialize() works", {
    x <- MsBackendMzSpecLib()
    res <- .initialize(x, file = small_file)

    expect_s4_class(res, "MsBackendMzSpecLib")
    expect_equal(res@fileNames, small_file)
    expect_equal(nrow(res@index), 4L)
    expect_equal(unique(res@index$file), 1L)
    expect_equal(nrow(res@header), 3L)
})

test_that(".open_mzspeclib() works", {
    res <- .open_mzspeclib(1, small_file)
    expect_named(res, c("index", "header"))
    expect_equal(names(res$index),
                     c("file", "tag", "kind", "attrset", "key", "start",
                       "end", "peak_line"))
    expect_true(all(res$index$file == 1L))
    expect_equal(nrow(res$index), 4L)
    expect_equal(.header_value(res$header, "MS:1003188"), "mini")
    expect_equal(nrow(res$header), 3L)

    ## Missing <mzSpecLib> or not first tag
    f <- tempfile()
    writeLines(c("<AttributeSet Spectrum=all>",
                        "MS:1000511|ms level=2",
                        "<mzSpecLib>",
                        "MS:1003188|library name=x"), f)
    expect_error(.open_mzspeclib(1, f), "Mandatory <mzSpecLib> section")

    ## Open example file
    res_test <- .open_mzspeclib(1, test_file)
    expect_equal(nrow(res_test$index), 14L)
    expect_equal(nrow(res_test$header), 15L)
    expect_equal(.header_value(res_test$header, "MS:1003188"),
                     "fetal_brain_tiny_consensus_td")
})

test_that(".index_sections() works", {
    idx <- .index_sections(small_file)

    expect_s3_class(idx, "data.frame")
    expect_equal(idx$tag, c("<mzSpecLib>", "<AttributeSet Spectrum=all>",
                                "<Spectrum=1>", "<Spectrum=2>"))
    expect_equal(idx$kind, c("mzSpecLib", "Spectrum", "Spectrum",
                                 "Spectrum"))
    expect_equal(idx$attrset, c(FALSE, TRUE, FALSE, FALSE))
    expect_equal(idx$key, c(NA, "all", "1", "2"))
    expect_equal(idx$start, c(1L, 5L, 7L, 16L))
    expect_equal(idx$end, c(5L, 7L, 16L, Inf))
    expect_equal(idx$peak_line, c(NA, NA, 13L, 18L))

    ref <- .index_sections(small_file)
    for (cs in c(1, 2, 3, 7, 1e6))
        expect_equal(.index_sections(small_file, chunk_lines = cs), ref)

    ## File without tags
    f <- tempfile()
    writeLines(c("just", "some", "text"), f)
    expect_error(.index_sections(f), "No section tags found")
})

test_that(".read_block() works", {
    expect_equal(.read_block(small_file, 1, 4),
                c("MS:1003186|library format version=1.0",
                "MS:1003188|library name=mini",
                "[1]MS:1003207|library creation software=MS:1001477|SpectraST"))

    res <- .read_block(small_file, 5, 6)
    expect_equal(res, "MS:1000511|ms level=2")

    expect_null(.read_block(small_file, 6, 6))
    expect_null(.read_block(small_file, 6, 5))

    res <- .read_block(small_file, 16, Inf)
    expect_equal(res[1], "MS:1003061|library spectrum name=OTHER/3")
    expect_equal(res[length(res)], "300.3\t30\tb4")
    expect_length(res, 3L)
})

test_that(".header_value() works", {
    res <- .open_mzspeclib(1, small_file)$header
    expect_equal(.header_value(res, "MS:1003188"), "mini")
    expect_equal(.header_value(res, "MS:1003186"), "1.0")
    expect_equal(.header_value(res, "MS:0000000"), NA_character_)
    ## empty header
    expect_equal(.header_value(.parse_attributes(character()),
                                    "MS:1003188"),
                     NA_character_)
})

test_that(".parse_attributes() works", {
    res <- .parse_attributes(character())
    expect_s3_class(res, "data.frame")
    expect_equal(nrow(res), 0L)

    res <- .parse_attributes(NULL)
    expect_s3_class(res, "data.frame")
    expect_equal(nrow(res), 0L)

    res <- .parse_attributes(small_lines[1:4])
    expect_equal(names(res), c("group", "accession", "name", "value"))

    expect_equal(nrow(res), 3L)
    expect_equal(res$group, c(NA_integer_, NA_integer_, 1L))
    expect_equal(res$accession, c("MS:1003186", "MS:1003188", "MS:1003207"))
    expect_equal(res$name, c("library format version", "library name",
                                    "library creation software"))
    expect_equal(res$value, c("1.0", "mini", "MS:1001477|SpectraST"))

    res <- .parse_attributes("MS:1003206|library creation log=a=b;c=d")
    expect_equal(res$name, "library creation log")
    expect_equal(res$value, "a=b;c=d")

    res <- .parse_attributes(c('MS:1|"quoted name"=val',
                            'MS:1|"quoted name with = equal"=val',
                            "MS:2|empty=",
                            "UO:0000000|unit=UO:0000266|electronvolt",
                            'MS:3|log=> > > fetal_brain.splib : IMPORT FROM PepXML "/data/interact.pep.xml", (Comet against "UP000005640_9606.fasta" (AA); ) [P=0.9;q=9999;n=;g=FALSE;o=FALSE;I=;_RNT=0;_RDR=100000;_DCN=0;_NAA=6;_NPK=10;_MDF=9999;_CEN=FALSE;_XAN=FALSE;_BRK=FALSE;_BRM=FALSE;_IRT=;_IRR=FALSE]'))
    expect_equal(res$name, c("quoted name", "quoted name with = equal",
                                "empty", "unit", "log"))
    expect_equal(res$value[3], "")
    expect_equal(res$accession[4], "UO:0000000")
    expect_equal(res$value[5], '> > > fetal_brain.splib : IMPORT FROM PepXML "/data/interact.pep.xml", (Comet against "UP000005640_9606.fasta" (AA); ) [P=0.9;q=9999;n=;g=FALSE;o=FALSE;I=;_RNT=0;_RDR=100000;_DCN=0;_NAA=6;_NPK=10;_MDF=9999;_CEN=FALSE;_XAN=FALSE;_BRK=FALSE;_BRM=FALSE;_IRT=;_IRR=FALSE]')

    res <- .parse_attributes(c("garbage line",
                               "MS:1000511|ms level=2",
                               "100.1\t10\tb2",
                               "no_colon|name=value"))
    expect_equal(nrow(res), 1L)
    expect_equal(res$accession, "MS:1000511")
})

test_that(".parse_attributeset_tag() works", {
    res <- .parse_attributeset_tag(c("<mzSpecLib>", "<Spectrum=1>",
                                     "<Cluster=1>"))
    expect_named(res, c("attrset", "kind", "key"))
    expect_equal(res$attrset, c(FALSE, FALSE, FALSE))
    expect_equal(res$kind, c("mzSpecLib", "Spectrum", "Cluster"))
    expect_equal(res$key, c(NA_character_, "1", "1"))

    res <- .parse_attributeset_tag(c("<AttributeSet Spectrum=all>",
                                     "<AttributeSet Analyte=HUMAN_TRYPTIC>",
                                     "<AttributeSet Interpretation=abc>"))
    expect_true(all(res$attrset))
    expect_equal(res$kind, c("Spectrum", "Analyte", "Interpretation"))
    expect_equal(res$key, c("all", "HUMAN_TRYPTIC", "abc"))

    res <- .parse_attributeset_tag(c("not a tag", "<Spectrum=>"))
    expect_equal(res$attrset, c(FALSE, FALSE))
    expect_true(all(is.na(res$kind)))
    expect_true(all(is.na(res$key)))
})

test_that(".spectrum_variables works", {
    be_test <- backendInitialize(MsBackendMzSpecLib(), small_file_opt)
    idx <- be_test@index
    rows <- which(idx$attrset & idx$kind == "Spectrum")
    global_mtd <- setNames(lapply(rows, function(j)
        .parse_attributes(.read_block(be_test@fileNames[1], idx$start[j],
                                      idx$end[j]))),
        idx$key[rows])

    ## Spectrum 1
    sv <- .spectrum_variables(1L, .read_block(be_test@fileNames[1],
                                            idx$start[3],
                                            idx$peak_line[3]),
                                global_mtd)
    expect_equal(nrow(sv), 1L)
    expect_equal(sv$scanIndex, 1L)
    expect_equal(sv$msLevel, 2L)
    expect_equal(sv$library_spectrum_name, "PEPTIDE/2")
    expect_equal(sv$peak_attribute, "optional_col")

    ## Spectrum 2
    sv2 <- .spectrum_variables(2L,
                                .read_block(be_test@fileNames[1], idx$start[4],
                                            idx$peak_line[4]),
                                global_mtd)
    expect_equal(sv2$scanIndex, 2L)
    expect_equal(sv2$library_spectrum_name, "OTHER/3")

    ## Own terms take precedence over the global AttributeSet
    own <- c("MS:1000511|ms level=1", "MS:1003061|library spectrum name=X/1",
            "MS:1000465|scan polarity=MS:1000130|positive scan",
            "MS:1000894|retention time=7",
            "MS:1000525|spectrum representation=MS:1000127|centroid spectrum")
    sv3 <- .spectrum_variables(1L, own, global_mtd)
    expect_equal(sv3$msLevel, 1L)
    expect_equal(sv3$library_spectrum_name, "X/1")

    ## Isolation window limits are derived from target and offsets
    iw <- c(own,
            "MS:1000827|isolation window target m/z=500",
            "MS:1000828|isolation window lower offset=1.5",
            "MS:1000829|isolation window upper offset=2.5")
    sv4 <- .spectrum_variables(1L, iw, global_mtd)
    expect_equal(sv4$isolationWindowTargetMz, 500)
    expect_equal(sv4$isolationWindowLowerMz, 498.5)
    expect_equal(sv4$isolationWindowUpperMz, 502.5)
})

test_that(".merge_spectrum_mtd_section works", {
    be_test <- backendInitialize(MsBackendMzSpecLib(), small_file_opt)
    idx <- be_test@index
    rows <- which(idx$attrset & idx$kind == "Spectrum")
    global_mtd <- setNames(lapply(rows, function(j)
        .parse_attributes(.read_block(be_test@fileNames[1], idx$start[j],
                                      idx$end[j]))),
        idx$key[rows])

    mtd <- .merge_spectrum_mtd_section(1L,
                                .read_block(be_test@fileNames[1], idx$start[3],
                                            idx$peak_line[3]),
                                global_mtd)
    expect_equal(colnames(mtd), c("group", "accession", "name", "value",
                                    "spectra_variable","repeatable"))
    expect_false(anyNA(mtd$spectra_variable))
    ## Own term and global term are both present
    expect_equal(mtd$value[mtd$accession == "MS:1003061"], "PEPTIDE/2")
    expect_equal(mtd$value[mtd$accession == "MS:1003254"], "optional_col")
    expect_equal(mtd$value[mtd$accession == "MS:1000511"], "2")
    ## Analyte / Interpretation sections are not part of the spectrum metadata
    expect_false("MS:1000888" %in% mtd$accession)
    expect_false("MS:1003208" %in% mtd$accession)

    ## Empty list of global tables
    mtd_no_global <- .merge_spectrum_mtd_section(1L,
                                .read_block(be_test@fileNames[1], idx$start[3],
                                            idx$peak_line[3]),
                                list())
    expect_false("MS:1000511" %in% mtd_no_global$accession)
    expect_true("MS:1003061" %in% mtd_no_global$accession)

    ## Precedence: own terms > global terms
    own <- c("MS:1000511|ms level=1")
    mtd_priority <- .merge_spectrum_mtd_section(1L, own, global_mtd)
    expect_equal(mtd_priority$value[mtd_priority$accession == "MS:1000511"],
                "1")

    ## Not-repeatable fields are only kept once
    dup <- c("MS:1000511|ms level=1")
    mtd_dup <- .merge_spectrum_mtd_section(1L, dup, global_mtd)
    expect_equal(sum(mtd_dup$accession == "MS:1000511"), 1L)
    expect_equal(mtd_dup$value[mtd_dup$accession == "MS:1000511"], "1")
})

test_that(".split_by_tags works", {
    section <- small_lines[8:12]
    res <- .split_by_tags(section, .SPECTRUM_SECTION)
    expect_type(res, "list")
    expect_equal(names(res), c("Metadata", "<Analyte=1>", "<Interpretation=1>"))
    expect_equal(res[["Metadata"]], section[1])
    expect_equal(res[["<Analyte=1>"]], section[3])
    expect_equal(res[["<Interpretation=1>"]], section[5])

    ## No tags: everything is Metadata
    res2 <- .split_by_tags(section[c(1, 3)], .SPECTRUM_SECTION)
    expect_equal(names(res2), "Metadata")
    expect_equal(res2$Metadata, section[c(1, 3)])
})

test_that(".polarity works", {
    expect_equal(.polarity("MS:1000130|positive scan"), 1L)
    expect_equal(.polarity("MS:1000129|negative scan"), 0L)
    expect_equal(.polarity(NA_character_), -1L)
    expect_equal(.polarity("MS:9999999|unknown"), -1L)
})

test_that(".spectrum_representation works", {
    expect_equal(.spectrum_representation("MS:1000127|centroid spectrum"),
                     TRUE)
    expect_equal(.spectrum_representation("MS:1000128|profile spectrum"),
                     FALSE)
    expect_equal(.spectrum_representation(NA_character_), NA)
    expect_equal(.spectrum_representation("MS:9999999|unknown"), NA)
})

test_that(".rtime works", {
    ## seconds, unit column removed
    df <- data.frame(rtime = 12.5, rtimeUnit = "UO:0000010|second", other = 1)
    res <- .rtime(df)
    expect_equal(res$rtime, 12.5)
    expect_false("rtimeUnit" %in% names(res))

    ## no unit column: assumed seconds
    df <- data.frame(rtime = 7, other = 1)
    res <- .rtime(df)
    expect_equal(res$rtime, 7)
    expect_equal(names(res), c("rtime", "other"))

    ## minutes -> seconds, unit column removed
    df <- data.frame(rtime = 2, rtimeUnit = "UO:0000031|minute", other = 1)
    res <- .rtime(df)
    expect_equal(res$rtime, 120)
    expect_equal(names(res), c("rtime", "other"))
    expect_equal(res$other, 1)

    ## hours -> seconds
    df <- data.frame(rtime = 1.5, rtimeUnit = "UO:0000032|hour", other = 1)
    expect_equal(.rtime(df)$rtime, 5400)

    ## unknown unit: default seconds
    df <- data.frame(rtime = 7, rtimeUnit = "UO:9999999|unknown", other = 1)
    expect_equal(.rtime(df)$rtime, 7)
})

test_that(".parse_other_attributes works", {
    mtd <- data.frame(
        group = c(NA, NA, NA),
        accession = c("MS:1000511", "MS:1003275", "MS:1003276"),
        name = c("ms level", "other attribute name", "other attribute value"),
        value = c("2", "myVar", "myValue"),
        spectra_variable = c("msLevel", NA, NA))
    res <- .parse_other_attributes(mtd)
    expect_equal(nrow(res), 2L)
    expect_false(any(c("MS:1003275", "MS:1003276") %in% res$accession))
    expect_equal(res$value[res$accession == "MS:1000511"], "2")
    ## new row: name of the field -> value
    new <- res[res$accession == "MS:0", ]
    expect_equal(nrow(new), 1L)
    expect_equal(new$name, "myVar")
    expect_equal(new$value, "myValue")
    expect_equal(new$spectra_variable, "myVar")
})

test_that("convert_value works", {
    expect_equal(convert_value("2", "integer"), 2L)
    expect_equal(convert_value("2.6", "integer"), 3L)
    expect_equal(convert_value(c("1", "2.4"), "integer"), c(1L, 2L))
    expect_equal(convert_value("1.5", "numeric"), 1.5)
    expect_equal(convert_value(c("1", "2.5"), "numeric"), c(1, 2.5))
    expect_equal(convert_value("abc", "character"), "abc")
    expect_equal(convert_value("abc", "other"), "abc")
    expect_equal(convert_value(NA_character_, "integer"), NA_integer_)
    expect_equal(convert_value(NA_character_, "numeric"), NA_real_)
})

test_that(".peaks_table works", {
    idx <- be@index
    cols <- c("mz", "intensity", "annotation")

    p1 <- .peaks_table("1", be@fileNames[1], idx, cols)
    expect_s3_class(p1, "data.frame")
    expect_equal(colnames(p1), cols)
    expect_equal(nrow(p1), 2L)
    expect_equal(p1$mz, c(100.1, 200.2))
    expect_equal(p1$intensity, c(10, 20))
    expect_equal(p1$annotation, c("b2", "y3"))

    ## With optional columns
    be_test <- backendInitialize(MsBackendMzSpecLib(), small_file_opt)
    idx_test <- be_test@index
    cols_test <- peaksVariables(be_test, index = 1L)
    p2 <- .peaks_table("1", be_test@fileNames[1], idx_test, cols_test)
    expect_equal(nrow(p2), 2L)
    expect_equal(p2$mz, c(100.1, 200.2))
    expect_equal(p2$intensity, c(10, 20))
    expect_equal(p2$annotation, c("b2", "y3"))
    expect_equal(p2$optional_col, c(1, 2))
})

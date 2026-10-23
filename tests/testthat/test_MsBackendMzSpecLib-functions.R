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
    expect_equal(nrow(res_test$index), 24L)
    expect_equal(nrow(res_test$header), 12L)
    expect_equal(.header_value(res_test$header, "MS:1003188"),
                     "fetal_brain_tiny")
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
    expect_equal(.read_block(small_file, 1, 5),
                c("MS:1003186|library format version=1.0",
                "MS:1003188|library name=mini",
                "[1]MS:1003207|library creation software=MS:1001477|SpectraST"))

    res <- .read_block(small_file, 5, 7)
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

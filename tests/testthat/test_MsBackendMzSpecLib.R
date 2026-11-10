test_that("spectraVariableMapping,MsBackendMzSpecLib works", {
    res <- spectraVariableMapping(MsBackendMzSpecLib())
    expect_true(is.character(res))
    expect_equal(res, MsBackendMzSpecLib:::.SV2CV)
})

test_that("spectraVariables,MsBackendMzSpecLib works", {
    sv <- spectraVariables(be)
    expect_type(sv, "character")
    expect_equal(length(sv), 109)

    expect_true(all(c("msLevel", "rtime", "dataOrigin", "dataStorage",
                      "scanIndex", "mz", "intensity") %in% sv))
    expect_true(all(.VARMAP$spectra_variable %in% sv))
})

test_that("validity of MsBackendMzSpecLib works", {
    expect_true(validObject(MsBackendMzSpecLib()))
    expect_true(validObject(be))

    bad <- be
    bad@index$file[] <- 99L
    expect_error(validObject(bad), "must point into")

    bad <- be
    bad@fileNames <- character()
    expect_error(validObject(bad), "must point into")
})

test_that("backendInitialize,MsBackendMzSpecLib works", {
    expect_error(backendInitialize(MsBackendMzSpecLib()),
                    "'files' is mandatory")
    expect_error(backendInitialize(MsBackendMzSpecLib(), 1),
                "'files' is expected to be a character vector")
    expect_error(backendInitialize(MsBackendMzSpecLib(),
                                   file.path(tempdir(), "does_not_exist.txt")),
                "not found")

    ## empty files: returns object unchanged
    be_empty <- backendInitialize(MsBackendMzSpecLib(), character())
    expect_equal(MsBackendMzSpecLib(), be_empty)

    ## Small example
    be_small <- backendInitialize(MsBackendMzSpecLib(), small_file)
    expect_s4_class(be_small, "MsBackendMzSpecLib")
    expect_true(validObject(be_small))
    expect_equal(be_small@fileNames, normalizePath(small_file))
    expect_equal(nrow(be_small@index), 4L)
    expect_true(all(be_small@index$file == 1L))
    expect_equal(nrow(be_small@header), 3L)

    ## multiple files
    be_files <- backendInitialize(MsBackendMzSpecLib(),
                                    c(small_file, test_file))
    expect_equal(length(be_files@fileNames), 2L)
    expect_equal(sort(unique(be_files@index$file)), 1:2)
    expect_equal(nrow(be_files@index), 18L)
})

test_that("show,MsBackendMzSpecLib works", {
    expect_output(show(be), "MsBackendMzSpecLib with 2 spectra")
    expect_output(show(be), "library format version:")
    expect_output(show(be), "1.0", fixed = TRUE)
    expect_output(show(be), "files:")
    expect_output(show(be), basename(small_file), fixed = TRUE)

    ## empty backend
    expect_output(show(MsBackendMzSpecLib()), "with 0 spectra")
})

test_that("length,MsBackendMzSpecLib works", {
    expect_equal(length(MsBackendMzSpecLib()), 0L)
    expect_equal(length(be), 2L)
})

test_that("spectraData,MsBackendMzSpecLib works", {
    sd <- spectraData(be, BPPARAM = BiocParallel::SerialParam())
    expect_equal(nrow(sd), 2L)
    expect_equal(sd$scanIndex, 1:2)
    expect_equal(sd$msLevel, c(2L, 2L))
    expect_equal(sd$dataOrigin, rep(basename(be@fileNames), 2))
    expect_equal(sd$dataStorage, rep(be@fileNames, 2))
    expect_equal(sd$library_spectrum_name, c("PEPTIDE/2", "OTHER/3"))

    ## subset of columns
    sd2 <- spectraData(be, c("scanIndex", "msLevel"),
                        BPPARAM = BiocParallel::SerialParam())
    expect_equal(colnames(sd2), c("scanIndex", "msLevel"))
    expect_equal(nrow(sd2), 2L)

    ## unknown columns are ignored
    sd3 <- spectraData(be, c("scanIndex", "not_a_variable"),
                        BPPARAM = BiocParallel::SerialParam())
    expect_equal(sd3, DataFrame("scanIndex" = 1:2))

    ## peaks columns
    sd4 <- spectraData(be, c("scanIndex", "mz", "intensity"),
                        BPPARAM = BiocParallel::SerialParam())
    expect_equal(colnames(sd4), c("scanIndex", "mz", "intensity"))
    expect_equal(sd4$mz[[1]], c(100.1, 200.2))
    expect_equal(sd4$mz[[2]], 300.3)
    expect_equal(sd4$intensity[[1]], c(10, 20))
    expect_equal(sd4$intensity[[2]], 30)
    sd5 <- spectraData(be, c("scanIndex", "mz"),
                        BPPARAM = BiocParallel::SerialParam())
    expect_equal(colnames(sd5), c("scanIndex", "mz"))

    ## multiple files
    be_files <- backendInitialize(MsBackendMzSpecLib(),
                                    c(small_file, small_file))
    sd6 <- spectraData(be_files, c("scanIndex", "dataOrigin", "dataStorage"),
                        BPPARAM = BiocParallel::SerialParam())
    expect_equal(nrow(sd6), 4L)
    expect_equal(sd6$scanIndex, c(1:2, 1:2))
    expect_equal(sd6$dataOrigin,
                 rep(basename(be_files@fileNames), each = 2))

    ## Read real file
    sd_test_file <- spectraData(be_test_file, c("scanIndex", "msLevel"),
                        BPPARAM = BiocParallel::SerialParam())
    expect_equal(colnames(sd_test_file), c("scanIndex", "msLevel"))
    expect_equal(nrow(sd_test_file), 10L)
})

test_that("peaksData,MsBackendMzSpecLib works", {
    pd <- peaksData(be, BPPARAM = BiocParallel::SerialParam())
    ## one element per file, each with one element per spectrum
    expect_type(pd, "list")
    expect_length(pd, 1L)
    expect_length(pd[[1]], 2L)
    p1 <- pd[[1]][[1]]
    p2 <- pd[[1]][[2]]
    expect_equal(colnames(p1), c("mz", "intensity"))
    expect_equal(p1$mz, c(100.1, 200.2))
    expect_equal(p1$intensity, c(10, 20))
    expect_equal(p2$mz, 300.3)
    expect_equal(p2$intensity, 30)

    ## subset of columns
    pd2 <- peaksData(be, columns = "mz", BPPARAM = BiocParallel::SerialParam())
    expect_equal(colnames(pd2[[1]][[1]]), "mz")
    expect_equal(colnames(pd2[[1]][[2]]), "mz")

    ## unavailable column
    expect_error(peaksData(be, columns = c("mz", "not_a_column"),
                        BPPARAM = BiocParallel::SerialParam()),
                 "not available")

    ## File with optional columns
    be_test <- backendInitialize(MsBackendMzSpecLib(), small_file_opt)
    pd3 <- peaksData(be_test,
                    c("mz", "intensity", "annotation", "optional_col"),
                    BPPARAM = BiocParallel::SerialParam())
    expect_equal(colnames(pd3[[1]][[1]]),
                    c("mz", "intensity", "annotation", "optional_col"))
    expect_equal(colnames(pd3[[1]][[2]]),
                    c("mz", "intensity", "annotation"))
    expect_equal(pd3[[1]][[1]]$optional_col, c(1, 2))

    ## Read real file
    pd_test_file <- peaksData(be_test_file,
                              BPPARAM = BiocParallel::SerialParam())
    expect_type(pd_test_file, "list")
    expect_length(pd_test_file, 1L)
    expect_length(pd_test_file[[1]], 10L)
})

test_that("peaksVariables,MsBackendMzSpecLib works", {
    pv <- peaksVariables(be)
    expect_type(pv, "character")
    expect_equal(pv, c("mz", "intensity", "annotation"))

    be_test <- backendInitialize(MsBackendMzSpecLib(), small_file_opt)
    expect_equal(peaksVariables(be_test, index = 1L),
                c("mz", "intensity", "annotation", "optional_col"))
    expect_equal(peaksVariables(be_test, index = 2L),
                c("mz", "intensity", "annotation"))
})

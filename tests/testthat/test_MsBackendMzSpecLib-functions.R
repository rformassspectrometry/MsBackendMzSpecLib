test_that("MsBackendMzSpecLib works", {
    res <- MsBackendMzSpecLib()
    expect_s4_class(res, "MsBackendMzSpecLib")
})

test_that("spectraVariableMapping,MsBackendMzSpecLib works", {
    res <- spectraVariableMapping(MsBackendMzSpecLib())
    expect_true(is.character(res))
    expect_equal(res, MsBackendMzSpecLib:::.SV2CV)
})

test_that(".onLoad works", {
    res <- callr::r(function() loadNamespace("MsBackendMzSpecLib"))
    expect_true(length(get(".VARMAP", res)) > 0)
})

test_that(".onLoad works", {
    res <- callr::r(function() loadNamespace("MsBackendMzSpecLib"))
    expect_true(length(get(".VARMAP", res)) > 0)
    ## The `asNamespace()` call in `.onLoad()` makes unit testing difficult.
    ## workaround is to use a dummy function that calls `asNamespace()`, but
    ## that can be mocked.
    my_env <- new.env()
    res <- with_mocked_bindings(
        ".get_envir" = function(x) my_env,
        code = .onLoad("MsBackendMzSpecLib", "MsBackendMzSpecLib")
    )
    expect_true(all(c(".VARMAP", ".SV2CV", ".CV2SV") %in% names(my_env)))
})

test_that(".get_envir works", {
    res <- .get_envir("MsBackendMzSpecLib")
    expect_true(any(names(res) %in% c(".VARMAP")))
})

test_that(".load_mapping works", {
    res <- .load_mapping()
    expect_true(is.data.frame(res))
    expect_true(nrow(res) > 0)
    expect_true(all(c("cv_term", "spectra_variable") %in% colnames(res)))
})

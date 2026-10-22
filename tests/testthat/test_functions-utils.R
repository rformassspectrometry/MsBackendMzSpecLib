test_that("spectraVariableFromCvTerm works", {
    res <- spectraVariableFromCvTerm(c("something"))
    expect_true(is.character(res))
    expect_equal(res, c(something = NA_character_))
    res <- spectraVariableFromCvTerm(c("MS:1003065", "no", "MS:1002476"))
    expect_equal(res, c(`MS:1003065` = "spectrum_aggregation_type",
                        no = NA_character_,
                        `MS:1002476` = "ionMobilityDriftTime"))
})

test_that("spectraVariableToCvTerm works", {
    res <- spectraVariableToCvTerm(c("rtime", "dataOrigin", "not_exists"))
    expect_equal(res, c(rtime = "MS:1000894", dataOrigin = NA_character_,
                        not_exists = NA_character_))
})

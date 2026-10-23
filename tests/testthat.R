library(testthat)

test_file <- system.file("examples", "fetal_brain_tiny.mzspeclib.txt",
                          package = "MsBackendMzSpecLib")

## Minimal fil
small_lines <- c(
    "<mzSpecLib>",
    "MS:1003186|library format version=1.0",
    "MS:1003188|library name=mini",
    "[1]MS:1003207|library creation software=MS:1001477|SpectraST",
    "<AttributeSet Spectrum=all>",
    "MS:1000511|ms level=2",
    "<Spectrum=1>",
    "MS:1003061|library spectrum name=PEPTIDE/2",
    "<Analyte=1>",
    "MS:1000888|stripped peptide sequence=PEPTIDE",
    "<Interpretation=1>",
    "MS:1003208|interpretation attribute=x",
    "<Peaks>",
    "100.1\t10\tb2",
    "200.2\t20\ty3",
    "<Spectrum=2>",
    "MS:1003061|library spectrum name=OTHER/3",
    "<Peaks>",
    "300.3\t30\tb4"
)
small_file <- tempfile()
writeLines(small_lines, small_file)

test_check("MsBackendMzSpecLib")

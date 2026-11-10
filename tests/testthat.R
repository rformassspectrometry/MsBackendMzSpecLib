library(testthat)
library(MsBackendMzSpecLib)

test_file <- system.file("examples",
                         "fetal_brain_tiny_consensus_td.mzSpecLib.txt",
                         package = "MsBackendMzSpecLib")
be_test_file <- backendInitialize(MsBackendMzSpecLib(), test_file)

## Minimal file
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

be <- backendInitialize(MsBackendMzSpecLib(), small_file)

## Minimal file with optional columns
small_lines_opt <- c(
    "<mzSpecLib>",
    "MS:1003186|library format version=1.0",
    "MS:1003188|library name=mini",
    "[1]MS:1003207|library creation software=MS:1001477|SpectraST",
    "<AttributeSet Spectrum=all>",
    "MS:1000511|ms level=2",
    "<Spectrum=1>",
    "MS:1003061|library spectrum name=PEPTIDE/2",
    "MS:1003254|peak attribute=optional_col",
    "<Analyte=1>",
    "MS:1000888|stripped peptide sequence=PEPTIDE",
    "<Interpretation=1>",
    "MS:1003208|interpretation attribute=x",
    "<Peaks>",
    "100.1\t10\tb2\t1",
    "200.2\t20\ty3\t2",
    "<Spectrum=2>",
    "MS:1003061|library spectrum name=OTHER/3",
    "<Peaks>",
    "300.3\t30\tb4"
)
small_file_opt <- tempfile()
writeLines(small_lines_opt, small_file_opt)

test_check("MsBackendMzSpecLib")

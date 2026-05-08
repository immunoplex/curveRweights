## ---- data-raw/make_example_data.R -------------------------------------------
## Run this script ONCE to create the package data files in data/
##
## In RStudio with the curveRweights project open:
##   source("data-raw/make_example_data.R")
##
## Or from the R console:
##   setwd("path/to/curveRweights")
##   source("data-raw/make_example_data.R")
## -----------------------------------------------------------------------------

library(dplyr)

# ---- Full dataset: example_assay --------------------------------------------
# 48,224 observations from a Luminex multiplex immunoassay.
# 11 antigens x 10 features x 150 subjects x 4 timepoints x 15 plates.
# Anonymised: patientid, nominal_sample_dilution, and group_b replaced
# with non-identifiable labels.

example_assay <- read.csv(
  "data-raw/dat_example_assay.csv",
  stringsAsFactors = FALSE
)

# Ensure consistent types
example_assay <- example_assay |>
  mutate(
    plate                   = as.character(plate),
    nominal_sample_dilution = as.character(nominal_sample_dilution),
    patientid               = as.character(patientid),
    group_a                 = as.character(group_a),
    group_b                 = as.character(group_b),
    antigen                 = as.character(antigen),
    feature                 = as.character(feature),
    mfi                     = as.numeric(mfi),
    predicted_concentration = as.numeric(predicted_concentration),
    se_concentration        = as.numeric(se_concentration),
    pcov                    = as.numeric(pcov)
  )

# Save full dataset
if (!dir.exists("data")) dir.create("data")
usethis::use_data(example_assay, overwrite = TRUE)

cat("Created data/example_assay.rda\n")
cat("  ", nrow(example_assay), "rows x", ncol(example_assay), "columns\n")
cat("  antigens:", paste(sort(unique(example_assay$antigen)), collapse = ", "), "\n")
cat("  features:", paste(sort(unique(example_assay$feature)), collapse = ", "), "\n")
cat("  groups_a:", paste(sort(unique(example_assay$group_a)), collapse = ", "), "\n")
cat("  groups_b:", paste(sort(unique(example_assay$group_b)), collapse = ", "), "\n")
cat("  plates:",   length(unique(example_assay$plate)), "\n")
cat("  subjects:", length(unique(example_assay$patientid)), "\n")

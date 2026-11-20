# Some tests check only the SHAPE of a result (which fields exist, which
# columns are written), so they run deliberately short MCMC chains to keep the
# suite fast. Short chains trigger Stan's low effective-sample-size warnings,
# which say nothing about what those tests check. quiet_ess() silences exactly
# those warnings and lets every other warning through. Tests that check an
# estimate run chains long enough not to need it.
quiet_ess <- function(expr) {
  withCallingHandlers(expr, warning = function(w) {
    if (grepl("Effective Samples Size (ESS) is too low", conditionMessage(w), fixed = TRUE))
      invokeRestart("muffleWarning")
  })
}

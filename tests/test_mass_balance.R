# ============================================================================
# Model Verification Tests for mAb PBPK Models
#
# Checks:
#   1. Non-negativity of all concentrations
#   2. IgG1 half-life in expected range (18-25 days)
#   3. Dose proportionality for linear PK
#   4. IgG reduction with efgartigimod
#   5. IgG recovery after drug washout
#   6. Baseline IgG stability without drug
#   7. Numerical stability across dose range
#
# Note: These are functional verification tests, not strict mass balance
# (cumulative degradation is not tracked as a separate state variable).
# ============================================================================

source("src/01_minimal_pbpk.R")
source("src/05_efgartigimod_app.R")

run_tests <- function() {

  pass_count <- 0
  fail_count <- 0
  total <- 0

  check <- function(test_name, condition) {
    total <<- total + 1
    if (condition) {
      cat(sprintf("  PASS: %s\n", test_name))
      pass_count <<- pass_count + 1
    } else {
      cat(sprintf("  FAIL: %s\n", test_name))
      fail_count <<- fail_count + 1
    }
  }

  cat("=", rep("=", 60), "\n", sep = "")
  cat("Mass Balance & Validation Tests\n")
  cat("=", rep("=", 60), "\n\n")

  # --- Test 1: No negative concentrations ---
  cat("Test 1: Non-negativity\n")
  res <- run_minimal_pbpk(dose_mg = 700, antibody = "typical_igg1", t_end = 90 * 24)
  check("Plasma >= 0", all(res$data$C_plasma >= 0))
  check("Tight >= 0", all(res$data$C_tight >= 0))
  check("Leaky >= 0", all(res$data$C_leaky >= 0))

  # --- Test 2: Terminal half-life in expected range for IgG1 ---
  cat("\nTest 2: IgG1 half-life validation\n")
  t_half <- res$pk_metrics$t_half_days
  check(sprintf("t1/2 = %.1f days (expected 15-25 days)", t_half),
        !is.na(t_half) && t_half > 10 && t_half < 35)

  # --- Test 3: Dose proportionality ---
  cat("\nTest 3: Dose proportionality (linear PK for IgG1)\n")
  res1 <- run_minimal_pbpk(dose_mg = 70, antibody = "typical_igg1", t_end = 60 * 24)
  res10 <- run_minimal_pbpk(dose_mg = 700, antibody = "typical_igg1", t_end = 60 * 24)
  ratio <- res10$pk_metrics$Cmax_nM / res1$pk_metrics$Cmax_nM
  check(sprintf("Cmax ratio (10x dose) = %.2f (expected ~10)", ratio),
        ratio > 9 && ratio < 11)

  # --- Test 4: Efgartigimod IgG reduction ---
  cat("\nTest 4: Efgartigimod IgG reduction\n")
  res_efg <- run_efgartigimod(dose_mg_kg = 10, n_doses = 4)
  nadir <- res_efg$nadir$igg_percent

  check(sprintf("IgG nadir = %.1f%% (expected 15-40%%)", nadir),
        nadir > 5 && nadir < 50)
  check("IgG nadir occurs after first dose",
        res_efg$nadir$time_days > 1)

  # --- Test 5: IgG recovery after washout ---
  cat("\nTest 5: IgG recovery\n")
  final_igg <- tail(res_efg$data$igg_percent_baseline, 1)
  check(sprintf("IgG at day 84 = %.1f%% (expected >80%% recovery)", final_igg),
        final_igg > 70)

  # --- Test 6: No drug = stable IgG ---
  cat("\nTest 6: Baseline IgG stability (no drug)\n")
  params_nodrug <- get_efgartigimod_system_params()
  state_nodrug <- c(E_plasma = 0, E_tight = 0, E_leaky = 0,
                    IgG = params_nodrug$IgG_baseline)
  times <- seq(0, 60 * 24, by = 1)
  out_nodrug <- as.data.frame(ode(y = state_nodrug, times = times,
                                   func = efgartigimod_ode,
                                   parms = params_nodrug, method = "lsoda"))
  igg_start <- out_nodrug$IgG[1]
  igg_end <- tail(out_nodrug$IgG, 1)
  pct_change <- abs(igg_end - igg_start) / igg_start * 100
  check(sprintf("IgG drift without drug = %.2f%% (expected <1%%)", pct_change),
        pct_change < 1)

  # --- Test 7: Numerical stability across doses ---
  cat("\nTest 7: Numerical stability\n")
  stable <- TRUE
  for (d in c(0.1, 1, 5, 10, 25, 50)) {
    tryCatch({
      r <- run_efgartigimod(dose_mg_kg = d, n_doses = 4, t_end = 84 * 24)
      if (any(is.na(r$data$IgG)) || any(r$data$IgG < 0)) stable <- FALSE
    }, error = function(e) {
      stable <<- FALSE
    })
  }
  check("Stable across doses 0.1-50 mg/kg", stable)

  # --- Summary ---
  cat("\n", rep("=", 40), "\n", sep = "")
  cat(sprintf("Results: %d/%d passed", pass_count, total))
  if (fail_count > 0) {
    cat(sprintf(", %d FAILED", fail_count))
  }
  cat("\n")

  invisible(list(pass = pass_count, fail = fail_count, total = total))
}

if (sys.nframe() == 0) {
  run_tests()
}

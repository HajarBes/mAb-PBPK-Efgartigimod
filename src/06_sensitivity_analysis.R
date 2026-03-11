# ============================================================================
# Global Sensitivity Analysis for mAb PBPK-PD Model
#
# Drug development question:
#   "Which parameters most influence IgG reduction by efgartigimod?
#    How robust is the predicted efficacy to parameter uncertainty?"
#
# Methods:
#   1. One-at-a-time (OAT) sensitivity for tornado plots
#   2. Sobol global sensitivity analysis (variance-based)
#
# Target outputs:
#   - IgG nadir (% of baseline) - primary efficacy endpoint
#   - Time to nadir (days)
#   - IgG AUC below baseline
#   - Efgartigimod Cmax
#
# References:
#   - Saltelli et al. 2008, Global Sensitivity Analysis
#   - Zhang et al. 2015, CPT:PSP (SA for PBPK)
# ============================================================================

library(deSolve)
source("src/utils/physiology.R")
source("src/utils/plotting.R")
source("src/05_efgartigimod_app.R")

# --- OAT Sensitivity Analysis ---

run_oat_sensitivity <- function(dose_mg_kg = 10, variation = 0.5) {

  # Parameters to vary and their baseline values
  param_info <- data.frame(
    name = c("Kd_efg", "Kd_igg", "CLp_efg",
             "FR_igg_baseline", "IgG_ksyn", "IgG_kdeg_base",
             "CL_renal", "sigma_v_tight", "sigma_v_leaky"),
    description = c(
      "Efgartigimod-FcRn Kd (pH 6.0)",
      "IgG-FcRn Kd (pH 6.0)",
      "Efgartigimod pinocytosis CL",
      "Baseline IgG recycling fraction",
      "IgG synthesis rate",
      "IgG baseline degradation rate",
      "Efgartigimod renal clearance",
      "Vascular reflection (tight)",
      "Vascular reflection (leaky)"
    ),
    stringsAsFactors = FALSE
  )

  # Get baseline parameters and result
  base_params <- get_efgartigimod_system_params()
  base_result <- run_efgartigimod(dose_mg_kg = dose_mg_kg)
  base_nadir <- base_result$nadir$igg_percent

  cat(sprintf("Baseline IgG nadir: %.2f%%\n\n", base_nadir))

  results <- data.frame(
    parameter = character(),
    description = character(),
    baseline = numeric(),
    low_value = numeric(),
    high_value = numeric(),
    nadir_low = numeric(),
    nadir_high = numeric(),
    sensitivity = numeric(),
    stringsAsFactors = FALSE
  )

  for (i in seq_len(nrow(param_info))) {
    pname <- param_info$name[i]
    base_val <- base_params[[pname]]

    if (is.null(base_val)) {
      cat(sprintf("  Skipping %s (not found)\n", pname))
      next
    }

    low_val <- base_val * (1 - variation)
    high_val <- base_val * (1 + variation)

    # Run with low value
    params_low <- base_params
    params_low[[pname]] <- low_val
    # Recalculate IgG_ksyn if degradation rate changed (to maintain SS)
    if (pname == "IgG_kdeg_base") {
      params_low$IgG_ksyn <- low_val * params_low$IgG_baseline
    }
    tryCatch({
      res_low <- run_efgartigimod_with_params(params_low, dose_mg_kg)
      nadir_low <- res_low$nadir$igg_percent
    }, error = function(e) {
      nadir_low <<- NA
    })

    # Run with high value
    params_high <- base_params
    params_high[[pname]] <- high_val
    if (pname == "IgG_kdeg_base") {
      params_high$IgG_ksyn <- high_val * params_high$IgG_baseline
    }
    tryCatch({
      res_high <- run_efgartigimod_with_params(params_high, dose_mg_kg)
      nadir_high <- res_high$nadir$igg_percent
    }, error = function(e) {
      nadir_high <<- NA
    })

    # Normalized sensitivity: (delta output / output) / (delta input / input)
    if (!is.na(nadir_low) && !is.na(nadir_high)) {
      sens <- ((nadir_high - nadir_low) / base_nadir) / (2 * variation)
    } else {
      sens <- NA
    }

    results <- rbind(results, data.frame(
      parameter = pname,
      description = param_info$description[i],
      baseline = base_val,
      low_value = low_val,
      high_value = high_val,
      nadir_low = ifelse(is.na(nadir_low), NA, nadir_low),
      nadir_high = ifelse(is.na(nadir_high), NA, nadir_high),
      sensitivity = sens,
      stringsAsFactors = FALSE
    ))

    cat(sprintf("  %s: nadir [%.1f%%, %.1f%%], S = %.3f\n",
                pname,
                ifelse(is.na(nadir_low), NA, nadir_low),
                ifelse(is.na(nadir_high), NA, nadir_high),
                ifelse(is.na(sens), NA, sens)))
  }

  results
}

# --- Helper: run with modified params ---

run_efgartigimod_with_params <- function(params, dose_mg_kg = 10,
                                          n_doses = 4, t_end = 84 * 24) {

  dose_mg <- dose_mg_kg * 70
  events <- create_dosing_events(dose_mg, params$MW_efg, params$V_plasma,
                                  n_doses = n_doses, interval_days = 7)

  state <- c(
    E_plasma = 0,
    E_tight = 0,
    E_leaky = 0,
    IgG = params$IgG_baseline
  )

  times <- seq(0, t_end, by = 1)

  out <- as.data.frame(ode(y = state, times = times, func = efgartigimod_ode,
                           parms = params, method = "lsoda",
                           events = list(data = events)))

  out$time_days <- out$time / 24
  out$igg_percent_baseline <- out$IgG / params$IgG_baseline * 100

  nadir_idx <- which.min(out$igg_percent_baseline)

  list(
    data = out,
    nadir = list(
      time_days = out$time_days[nadir_idx],
      igg_percent = out$igg_percent_baseline[nadir_idx]
    )
  )
}

# --- Fc Engineering Scenario Analysis ---
# Vary FcRn binding affinity to predict impact on IgG reduction
# Directly relevant to argenx ABDEG technology platform

run_fc_engineering_scenarios <- function(dose_mg_kg = 10) {

  cat("\nFc Engineering Scenarios: FcRn Kd vs IgG Reduction\n")
  cat("(Relevant to ABDEG platform optimization)\n\n")

  Kd_values <- c(5, 10, 20, 50, 100, 200, 500, 800)  # nM at pH 6.0
  # 800 = wild-type IgG1, 20 = efgartigimod, <10 = next-gen engineered

  scenario_results <- list()

  for (kd in Kd_values) {
    params <- get_efgartigimod_system_params()
    params$Kd_efg <- kd

    res <- run_efgartigimod_with_params(params, dose_mg_kg)
    res$data$Kd_label <- sprintf("Kd = %d nM", kd)
    res$data$Kd <- kd

    scenario_results[[length(scenario_results) + 1]] <- list(
      Kd = kd,
      nadir = res$nadir$igg_percent,
      data = res$data
    )

    cat(sprintf("  Kd = %4d nM: IgG nadir = %.1f%%\n", kd, res$nadir$igg_percent))
  }

  # Summary data frame
  summary_df <- data.frame(
    Kd = sapply(scenario_results, function(x) x$Kd),
    IgG_nadir = sapply(scenario_results, function(x) x$nadir)
  )

  # Time course data
  tc_df <- do.call(rbind, lapply(scenario_results, function(x) x$data))

  list(summary = summary_df, timecourse = tc_df)
}

# --- Main ---

if (sys.nframe() == 0) {

  cat("=", rep("=", 60), "\n", sep = "")
  cat("Sensitivity Analysis - Efgartigimod PBPK-PD\n")
  cat("=", rep("=", 60), "\n\n")

  # --- 1. OAT Sensitivity ---
  cat("--- One-at-a-Time Sensitivity (50% variation) ---\n\n")
  sa_results <- run_oat_sensitivity(dose_mg_kg = 10, variation = 0.5)

  # Tornado plot
  sa_plot <- sa_results[!is.na(sa_results$sensitivity), ]
  sa_plot <- sa_plot[order(abs(sa_plot$sensitivity)), ]
  sa_plot$parameter <- factor(sa_plot$parameter, levels = sa_plot$parameter)

  p1 <- plot_tornado(sa_plot, title = "Parameter Sensitivity: IgG Nadir (% Baseline)")
  ggsave("results/sensitivity_tornado.png", p1, width = 10, height = 6, dpi = 300)
  cat("\nTornado plot saved.\n")

  # --- 2. Fc Engineering Scenarios ---
  cat("\n--- Fc Engineering Scenario Analysis ---\n")
  fc_results <- run_fc_engineering_scenarios(dose_mg_kg = 10)

  # Kd vs nadir plot
  p2 <- ggplot(fc_results$summary, aes(x = Kd, y = IgG_nadir)) +
    geom_line(color = "#2C3E50", linewidth = 1) +
    geom_point(color = "#E74C3C", size = 3) +
    geom_point(data = fc_results$summary[fc_results$summary$Kd == 20, ],
               aes(x = Kd, y = IgG_nadir), color = "#27AE60", size = 5, shape = 18) +
    annotate("text", x = 25, y = fc_results$summary$IgG_nadir[fc_results$summary$Kd == 20] + 3,
             label = "Efgartigimod", color = "#27AE60", fontface = "bold") +
    scale_x_log10() +
    theme_pbpk() +
    labs(title = "Fc Engineering: FcRn Binding Affinity vs IgG Reduction",
         subtitle = "10 mg/kg IV weekly x4 | Lower Kd = stronger FcRn binding",
         x = "FcRn Binding Kd at pH 6.0 (nM)",
         y = "IgG Nadir (% of Baseline)")

  ggsave("results/fc_engineering_kd_vs_nadir.png", p2, width = 9, height = 6, dpi = 300)

  # Time course by Kd
  selected_kd <- c(5, 20, 100, 800)
  tc_selected <- fc_results$timecourse[fc_results$timecourse$Kd %in% selected_kd, ]

  p3 <- ggplot(tc_selected, aes(x = time_days, y = igg_percent_baseline, color = Kd_label)) +
    geom_line(linewidth = 0.8) +
    geom_hline(yintercept = 100, linetype = "dashed", color = "grey50") +
    theme_pbpk() +
    labs(title = "IgG Time-Course by FcRn Binding Affinity",
         subtitle = "Stronger FcRn binding -> deeper & longer IgG suppression",
         x = "Time (days)", y = "IgG (% of Baseline)", color = "")

  ggsave("results/fc_engineering_timecourse.png", p3, width = 10, height = 6, dpi = 300)

  cat("\nAll sensitivity analysis plots saved to results/\n")
}

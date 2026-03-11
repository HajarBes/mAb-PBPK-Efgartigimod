# ============================================================================
# Fc Engineering Scenario Analysis
#
# Drug development question:
#   "How does FcRn binding affinity engineering affect mAb half-life
#    and therapeutic window? What is the optimal Kd for FcRn-targeting
#    therapeutics like efgartigimod?"
#
# Context:
#   argenx's ABDEG technology engineers Fc domains for enhanced FcRn binding.
#   Understanding the Kd-efficacy relationship guides:
#   - Lead candidate selection
#   - Structure-activity relationship (SAR) for Fc variants
#   - Prediction of next-generation molecules
#
# Scenarios:
#   1. Vary Kd at pH 6.0: 1 nM to 1000 nM
#   2. Vary pH-dependent release: Kd(pH7.4) / Kd(pH6.0) ratio
#   3. Half-life extension (for non-FcRn-blocker mAbs with enhanced FcRn binding)
# ============================================================================

source("src/05_efgartigimod_app.R")
source("src/06_sensitivity_analysis.R")
source("src/01_minimal_pbpk.R")

# --- Scenario 1: Kd sweep for FcRn-blocking therapeutics ---
# (efgartigimod-like: drug blocks FcRn to reduce IgG)

run_kd_sweep_fcrn_blocker <- function(dose_mg_kg = 10) {

  Kd_values <- 10^seq(log10(1), log10(2000), length.out = 30)

  results <- data.frame(
    Kd = numeric(),
    IgG_nadir_pct = numeric(),
    time_to_nadir_days = numeric(),
    duration_below_50pct = numeric()
  )

  for (kd in Kd_values) {
    params <- get_efgartigimod_system_params()
    params$Kd_efg <- kd

    res <- run_efgartigimod_with_params(params, dose_mg_kg)

    # Duration below 50% baseline
    below50 <- sum(res$data$igg_percent_baseline < 50) / 24  # hours -> days

    results <- rbind(results, data.frame(
      Kd = kd,
      IgG_nadir_pct = res$nadir$igg_percent,
      time_to_nadir_days = res$nadir$time_days,
      duration_below_50pct = below50
    ))
  }

  results
}

# --- Scenario 2: Half-life extension for therapeutic mAbs ---
# (Not FcRn blockers, but mAbs with enhanced FcRn binding for longer t1/2)

run_halflife_extension <- function() {

  cat("\nHalf-Life Extension: Enhanced FcRn Binding for Therapeutic mAbs\n\n")

  # For a normal therapeutic mAb (not FcRn blocker):
  # Enhanced FcRn binding at pH 6.0 (lower Kd) -> higher FR -> longer half-life
  # BUT must release at pH 7.4 (otherwise trapped on cell surface)

  Kd_values <- c(1000, 800, 500, 200, 100, 50, 20, 10)
  # 800 = wild-type IgG1
  # <100 = Fc-engineered (YTE, LS, ABDEG variants)

  results <- data.frame(
    Kd_pH6 = numeric(),
    FR_effective = numeric(),
    predicted_halflife_days = numeric()
  )

  for (kd in Kd_values) {
    # Effective recycling at physiological IgG + therapeutic mAb
    # Therapeutic mAb at ~1000 nM (typical trough), total IgG ~67000 nM
    C_total <- 67000
    fcrn <- get_fcrn_params()
    fcrn$Kd_pH6 <- kd

    fr <- compute_effective_FR(C_total, fcrn)

    # Half-life scales with 1/(1-FR)
    # Baseline: FR=0.715, t1/2=21 days
    t_half <- 21 * (1 - 0.715) / (1 - fr)

    results <- rbind(results, data.frame(
      Kd_pH6 = kd,
      FR_effective = fr,
      predicted_halflife_days = t_half
    ))

    cat(sprintf("  Kd = %4d nM: FR = %.3f, t1/2 = %.1f days\n", kd, fr, t_half))
  }

  results
}

# --- Main ---

if (sys.nframe() == 0) {

  cat("=", rep("=", 60), "\n", sep = "")
  cat("Fc Engineering Scenario Analysis\n")
  cat("=", rep("=", 60), "\n\n")

  # --- Scenario 1: FcRn blocker Kd sweep ---
  cat("--- Scenario 1: FcRn Blocker (Efgartigimod-like) ---\n")
  kd_results <- run_kd_sweep_fcrn_blocker(dose_mg_kg = 10)

  library(ggplot2)
  library(patchwork)

  p1a <- ggplot(kd_results, aes(x = Kd, y = IgG_nadir_pct)) +
    geom_line(color = "#2C3E50", linewidth = 1) +
    geom_point(color = "#E74C3C", size = 2) +
    geom_vline(xintercept = 20, linetype = "dashed", color = "#27AE60") +
    annotate("text", x = 25, y = max(kd_results$IgG_nadir_pct) - 5,
             label = "Efgartigimod\n(Kd=20 nM)", color = "#27AE60", size = 3) +
    scale_x_log10() +
    theme_pbpk() +
    labs(title = "A) IgG Nadir vs FcRn Binding Affinity",
         x = "Kd at pH 6.0 (nM)", y = "IgG Nadir (% Baseline)")

  p1b <- ggplot(kd_results, aes(x = Kd, y = duration_below_50pct)) +
    geom_line(color = "#8E44AD", linewidth = 1) +
    geom_point(color = "#8E44AD", size = 2) +
    geom_vline(xintercept = 20, linetype = "dashed", color = "#27AE60") +
    scale_x_log10() +
    theme_pbpk() +
    labs(title = "B) Duration of >50% IgG Reduction",
         x = "Kd at pH 6.0 (nM)", y = "Days Below 50% Baseline")

  p_combined <- p1a / p1b
  ggsave("results/fc_engineering_kd_sweep.png", p_combined,
         width = 10, height = 10, dpi = 300)

  # --- Scenario 2: Half-life extension ---
  cat("\n--- Scenario 2: Half-Life Extension ---\n")
  hl_results <- run_halflife_extension()

  p2 <- ggplot(hl_results, aes(x = Kd_pH6, y = predicted_halflife_days)) +
    geom_line(color = "#2980B9", linewidth = 1) +
    geom_point(color = "#2980B9", size = 3) +
    geom_hline(yintercept = 21, linetype = "dashed", color = "grey50") +
    annotate("text", x = 500, y = 23, label = "Wild-type IgG1 (21 days)",
             color = "grey50", size = 3) +
    geom_point(data = hl_results[hl_results$Kd_pH6 == 800, ],
               color = "#E74C3C", size = 4, shape = 18) +
    annotate("text", x = 600, y = hl_results$predicted_halflife_days[hl_results$Kd_pH6 == 800] - 2,
             label = "WT IgG1", color = "#E74C3C", size = 3) +
    scale_x_log10() +
    theme_pbpk() +
    labs(title = "Predicted mAb Half-Life vs FcRn Binding Affinity",
         subtitle = "Enhanced FcRn binding at pH 6.0 extends mAb half-life",
         x = "FcRn Kd at pH 6.0 (nM)",
         y = "Predicted Half-Life (days)")

  ggsave("results/halflife_extension.png", p2, width = 9, height = 6, dpi = 300)

  cat("\nAll Fc engineering plots saved to results/\n")
}

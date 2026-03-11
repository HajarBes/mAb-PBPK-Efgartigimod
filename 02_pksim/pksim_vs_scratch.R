# ============================================================================
# PK-Sim vs Custom PBPK Model Comparison
#
# Purpose:
#   Compare efgartigimod PK predictions from our custom R-based minimal PBPK
#   model (src/05_efgartigimod_app.R) against PK-Sim's organ-resolved
#   large-molecule PBPK.
#
# Workflow:
#   1. Run our custom model to get R-based predictions
#   2. Load PK-Sim exported CSV (or generate synthetic placeholder data)
#   3. Overlay both on the same time axis
#   4. Compute quantitative comparison metrics (RMSE, fold-error, AUC ratio)
#
# Note on synthetic data:
#   Until a real PK-Sim simulation is exported, this script generates
#   synthetic PK-Sim-like data based on published PK parameters for
#   efgartigimod (Ulrichts 2018, Fc fragment t1/2 ~4.8 days). The synthetic
#   data is clearly labeled and intended only as a placeholder for
#   demonstration. Replace it with actual PK-Sim output when available.
#
# Usage:
#   source("02_pksim/pksim_vs_scratch.R")
#   (run from project root: mAb-PBPK-Efgartigimod/)
#
# ============================================================================

library(deSolve)
library(ggplot2)

source("src/utils/physiology.R")
source("src/utils/plotting.R")
source("src/03_fcrn_recycling.R")
source("src/05_efgartigimod_app.R")

# ============================================================================
# SECTION 1: Run Custom Model
# ============================================================================

cat("Running custom R-based PBPK model...\n")
res_custom <- run_efgartigimod(dose_mg_kg = 10, n_doses = 4, t_end = 84 * 24)
df_custom <- res_custom$data

# Extract relevant columns
custom_pk <- data.frame(
  time_days = df_custom$time_days,
  conc_ug_mL = df_custom$E_plasma_ug_mL,
  igg_percent = df_custom$igg_percent_baseline,
  source = "Custom R model"
)

cat(sprintf("  Custom model: Cmax = %.1f ug/mL, IgG nadir = %.1f%%\n",
            max(custom_pk$conc_ug_mL), min(custom_pk$igg_percent)))

# ============================================================================
# SECTION 2: Load PK-Sim Results
# ============================================================================

# --- Option A: Load real PK-Sim CSV ---
# Uncomment and edit the path when a real PK-Sim export is available.
# PK-Sim CSV typically has columns like:
#   Time [h], Organism|VenousBlood|Plasma|Efgartigimod|Concentration [mg/l]
#
# load_pksim_csv <- function(filepath) {
#   raw <- read.csv(filepath, check.names = FALSE)
#
#   # Identify time and concentration columns (adapt to your export format)
#   time_col <- grep("Time", names(raw), value = TRUE)[1]
#   conc_col <- grep("Concentration", names(raw), value = TRUE)[1]
#
#   df <- data.frame(
#     time_h = raw[[time_col]],
#     conc_mg_L = raw[[conc_col]]
#   )
#
#   # Convert units: mg/L = ug/mL, time h -> days
#   df$time_days <- df$time_h / 24
#   df$conc_ug_mL <- df$conc_mg_L  # mg/L = ug/mL (same unit)
#
#   df
# }
#
# pksim_raw <- load_pksim_csv("02_pksim/data/pksim_efgartigimod_results.csv")

# --- Option B: Generate synthetic PK-Sim-like data (PLACEHOLDER) ---
# This uses a two-compartment model with parameters calibrated to produce
# PK profiles consistent with what PK-Sim would generate for an Fc fragment.
# Parameters are informed by Ulrichts 2018 and typical PK-Sim behavior
# for a 54 kDa protein with no FcRn recycling.

cat("\nGenerating synthetic PK-Sim-like data (placeholder)...\n")
cat("  *** Replace with real PK-Sim export when available ***\n\n")

generate_synthetic_pksim <- function(dose_mg_kg = 10, BW = 70, n_doses = 4,
                                      interval_days = 7, t_end_days = 84) {
  # Two-compartment IV model mimicking PK-Sim organ-resolved output
  # Parameters chosen to give slightly different (but plausible) PK vs our model,
  # reflecting structural differences between lumped and organ-resolved PBPK.

  dose_mg <- dose_mg_kg * BW
  MW <- 54000  # Da

  # PK-Sim-like parameters for Fc fragment (no FcRn recycling)
  # Slightly different from our model due to organ-resolved distribution
  V1 <- 3.3     # central volume (L) - PK-Sim typically gives slightly larger
  V2 <- 5.8     # peripheral volume (L)
  CL <- 0.42    # clearance (L/day) - t1/2 ~4.5 days (PK-Sim may give slightly faster)
  Q  <- 0.85    # intercompartmental CL (L/day)

  # Two-compartment ODE
  pk_ode <- function(t, state, pars) {
    C1 <- state[1]
    C2 <- state[2]
    with(as.list(pars), {
      dC1 <- -(CL/V1) * C1 - (Q/V1) * C1 + (Q/V2) * C2
      dC2 <-  (Q/V1) * C1 - (Q/V2) * C2
      list(c(dC1, dC2))
    })
  }

  pars <- c(V1 = V1, V2 = V2, CL = CL, Q = Q)

  # Build dosing events (bolus approximation of 1h infusion)
  dose_conc <- (dose_mg * 1e3) / (MW * V1)  # ug/mL in central compartment
  # More accurately: dose in ug / V1 in mL
  dose_conc_ug_mL <- dose_mg * 1e3 / (V1 * 1e3)  # mg -> ug, L -> mL

  events <- data.frame(
    var = rep("C1", n_doses),
    time = seq(0, by = interval_days, length.out = n_doses),
    value = rep(dose_conc_ug_mL, n_doses),
    method = rep("add", n_doses)
  )

  times <- seq(0, t_end_days, by = 1/24)  # hourly
  state <- c(C1 = 0, C2 = 0)

  out <- as.data.frame(ode(y = state, times = times, func = pk_ode,
                            parms = pars, method = "lsoda",
                            events = list(data = events)))

  data.frame(
    time_days = out$time,
    conc_ug_mL = out$C1,
    source = "PK-Sim (synthetic placeholder)"
  )
}

pksim_pk <- generate_synthetic_pksim()

# Add IgG column (NA for PK-Sim - it does not natively predict IgG PD)
pksim_pk$igg_percent <- NA

cat(sprintf("  Synthetic PK-Sim: Cmax = %.1f ug/mL\n", max(pksim_pk$conc_ug_mL)))

# ============================================================================
# SECTION 3: Overlay Comparison Plots
# ============================================================================

# Combine datasets
df_compare <- rbind(
  custom_pk[, c("time_days", "conc_ug_mL", "source")],
  pksim_pk[, c("time_days", "conc_ug_mL", "source")]
)

# --- Plot 1: Linear scale PK overlay ---
p1 <- ggplot(df_compare, aes(x = time_days, y = conc_ug_mL, color = source)) +
  geom_line(linewidth = 0.8) +
  geom_vline(xintercept = c(0, 7, 14, 21), linetype = "dotted",
             color = "grey70", linewidth = 0.3) +
  scale_color_manual(values = c("Custom R model" = "#2C3E50",
                                "PK-Sim (synthetic placeholder)" = "#E74C3C")) +
  theme_pbpk() +
  labs(title = "Efgartigimod PK: Custom Model vs PK-Sim",
       subtitle = "10 mg/kg IV weekly x4, 70 kg - Linear scale",
       x = "Time (days)", y = "Plasma Concentration (ug/mL)",
       color = "Model",
       caption = "PK-Sim data is synthetic placeholder - replace with actual export")

# --- Plot 2: Semi-log PK overlay ---
p2 <- ggplot(df_compare, aes(x = time_days, y = conc_ug_mL, color = source)) +
  geom_line(linewidth = 0.8) +
  scale_y_log10(limits = c(0.01, NA)) +
  geom_vline(xintercept = c(0, 7, 14, 21), linetype = "dotted",
             color = "grey70", linewidth = 0.3) +
  scale_color_manual(values = c("Custom R model" = "#2C3E50",
                                "PK-Sim (synthetic placeholder)" = "#E74C3C")) +
  theme_pbpk() +
  labs(title = "Efgartigimod PK: Custom Model vs PK-Sim",
       subtitle = "10 mg/kg IV weekly x4, 70 kg - Semi-log scale",
       x = "Time (days)", y = "Plasma Concentration (ug/mL)",
       color = "Model",
       caption = "PK-Sim data is synthetic placeholder - replace with actual export")

# --- Plot 3: Residuals (fold-error over time) ---
# Interpolate PK-Sim onto custom model time grid
pksim_interp <- approx(pksim_pk$time_days, pksim_pk$conc_ug_mL,
                        xout = custom_pk$time_days, rule = 2)

df_resid <- data.frame(
  time_days = custom_pk$time_days,
  custom = custom_pk$conc_ug_mL,
  pksim = pksim_interp$y
)

# Fold-error: only compute where both concentrations are measurable
df_resid$fold_error <- ifelse(
  df_resid$pksim > 0.01 & df_resid$custom > 0.01,
  df_resid$custom / df_resid$pksim,
  NA
)

p3 <- ggplot(df_resid, aes(x = time_days, y = fold_error)) +
  geom_line(color = "#2C3E50", linewidth = 0.8) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = c(0.5, 2), linetype = "dotted", color = "red") +
  annotate("text", x = 70, y = 2.1, label = "2-fold boundary", color = "red", size = 3) +
  annotate("text", x = 70, y = 0.45, label = "0.5-fold boundary", color = "red", size = 3) +
  geom_vline(xintercept = c(0, 7, 14, 21), linetype = "dotted",
             color = "grey70", linewidth = 0.3) +
  coord_cartesian(ylim = c(0, 3)) +
  theme_pbpk() +
  labs(title = "Fold-Error: Custom Model / PK-Sim",
       subtitle = "Values within 0.5-2.0 indicate acceptable agreement",
       x = "Time (days)", y = "Fold-Error (Custom / PK-Sim)")

# Save plots
ggsave("results/pksim_comparison_linear.png", p1, width = 10, height = 6, dpi = 300)
ggsave("results/pksim_comparison_semilog.png", p2, width = 10, height = 6, dpi = 300)
ggsave("results/pksim_comparison_folderror.png", p3, width = 10, height = 5, dpi = 300)

cat("Comparison plots saved to results/\n\n")

# ============================================================================
# SECTION 4: Quantitative Comparison Metrics
# ============================================================================

cat("=== Quantitative Comparison Metrics ===\n\n")

# Only use time points where both models have meaningful concentrations
valid <- !is.na(df_resid$fold_error) & df_resid$custom > 0.01 & df_resid$pksim > 0.01

if (sum(valid) > 10) {

  custom_valid <- df_resid$custom[valid]
  pksim_valid <- df_resid$pksim[valid]
  fold_valid <- df_resid$fold_error[valid]

  # --- RMSE ---
  rmse <- sqrt(mean((custom_valid - pksim_valid)^2))
  cat(sprintf("RMSE (ug/mL):                     %.2f\n", rmse))

  # --- Normalized RMSE (as % of max concentration) ---
  nrmse <- rmse / max(custom_valid) * 100
  cat(sprintf("Normalized RMSE (%% of Cmax):      %.1f%%\n", nrmse))

  # --- Mean Absolute Fold-Error (MAFE) ---
  mafe <- mean(abs(log10(fold_valid)))
  cat(sprintf("Mean Absolute Fold-Error (log10):  %.3f\n", mafe))

  # --- Percentage within 2-fold ---
  within_2fold <- mean(fold_valid >= 0.5 & fold_valid <= 2.0) * 100
  cat(sprintf("Within 2-fold (%%):                %.1f%%\n", within_2fold))

  # --- Cmax comparison ---
  cmax_custom <- max(custom_pk$conc_ug_mL)
  cmax_pksim <- max(pksim_pk$conc_ug_mL)
  cat(sprintf("\nCmax custom:  %.1f ug/mL\n", cmax_custom))
  cat(sprintf("Cmax PK-Sim:  %.1f ug/mL\n", cmax_pksim))
  cat(sprintf("Cmax ratio:   %.2f\n", cmax_custom / cmax_pksim))

  # --- AUC comparison (trapezoidal rule) ---
  auc_trap <- function(time, conc) {
    n <- length(time)
    sum(diff(time) * (conc[-n] + conc[-1]) / 2)
  }

  auc_custom <- auc_trap(custom_pk$time_days, custom_pk$conc_ug_mL)
  auc_pksim <- auc_trap(pksim_pk$time_days, pksim_pk$conc_ug_mL)
  cat(sprintf("\nAUC0-84d custom:  %.1f ug*day/mL\n", auc_custom))
  cat(sprintf("AUC0-84d PK-Sim:  %.1f ug*day/mL\n", auc_pksim))
  cat(sprintf("AUC ratio:        %.2f\n", auc_custom / auc_pksim))

  # --- Terminal half-life comparison ---
  # Estimate t1/2 from terminal phase (day 30-60, after last dose)
  estimate_half_life <- function(time, conc, start_day = 30, end_day = 60) {
    idx <- time >= start_day & time <= end_day & conc > 0.01
    if (sum(idx) < 10) return(NA)
    fit <- lm(log(conc[idx]) ~ time[idx])
    kel <- -coef(fit)[2]
    log(2) / kel
  }

  t12_custom <- estimate_half_life(custom_pk$time_days, custom_pk$conc_ug_mL)
  t12_pksim <- estimate_half_life(pksim_pk$time_days, pksim_pk$conc_ug_mL)
  cat(sprintf("\nTerminal t1/2 custom:  %.1f days\n", t12_custom))
  cat(sprintf("Terminal t1/2 PK-Sim:  %.1f days\n", t12_pksim))

  # --- Summary table ---
  cat("\n=== Summary ===\n")
  cat(sprintf("%-30s %-15s %-15s %-10s\n", "Metric", "Custom", "PK-Sim", "Ratio"))
  cat(paste(rep("-", 70), collapse = ""), "\n")
  cat(sprintf("%-30s %-15.1f %-15.1f %-10.2f\n", "Cmax (ug/mL)", cmax_custom, cmax_pksim, cmax_custom/cmax_pksim))
  cat(sprintf("%-30s %-15.1f %-15.1f %-10.2f\n", "AUC 0-84d (ug*day/mL)", auc_custom, auc_pksim, auc_custom/auc_pksim))
  cat(sprintf("%-30s %-15.1f %-15.1f %-10.2f\n", "Terminal t1/2 (days)", t12_custom, t12_pksim, t12_custom/t12_pksim))
  cat(sprintf("%-30s %-15s %-15s %-10.1f%%\n", "Within 2-fold", "", "", within_2fold))
  cat(sprintf("%-30s %-15s %-15s %-10.3f\n", "MAFE (log10)", "", "", mafe))

} else {
  cat("Insufficient valid data points for comparison metrics.\n")
}

cat("\n=== Interpretation Notes ===\n")
cat("- RMSE < 20% of Cmax: good agreement for cross-platform PBPK comparison\n")
cat("- >80% within 2-fold: acceptable for PBPK predictions (Abduljalil 2022)\n")
cat("- AUC ratio 0.8-1.25: bioequivalence-level agreement\n")
cat("- Differences are expected due to structural model differences:\n")
cat("    Custom model: 3 lumped compartments, competitive FcRn inhibition\n")
cat("    PK-Sim: organ-resolved (15+ tissues), two-pore + endosomal kinetics\n")
cat("\n*** Current PK-Sim data is SYNTHETIC - metrics are illustrative only ***\n")
cat("*** Re-run after replacing with actual PK-Sim export for real comparison ***\n")

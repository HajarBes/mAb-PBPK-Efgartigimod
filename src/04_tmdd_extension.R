# ============================================================================
# Target-Mediated Drug Disposition (TMDD) Extension
#
# For therapeutic mAbs that bind a soluble or membrane-bound target,
# TMDD creates nonlinear PK: at low concentrations, target-mediated
# clearance dominates; at high concentrations, target is saturated
# and linear CL prevails.
#
# This module adds TMDD to the minimal PBPK model.
#
# TMDD equations (Mager & Jusko 2001):
#   dAb/dt = -kon*Ab*Target + koff*Complex + [PBPK terms]
#   dTarget/dt = ksyn - kdeg_T*Target - kon*Ab*Target + koff*Complex
#   dComplex/dt = kon*Ab*Target - koff*Complex - kint*Complex
#
# Application: generalized mAb with soluble target
# (efgartigimod's target is FcRn, handled separately in script 05)
#
# References:
#   - Mager & Jusko 2001, J Pharmacokinet Pharmacodyn
#   - Gibiansky et al. 2008, J Pharmacokinet Pharmacodyn (QSS approx)
# ============================================================================

library(deSolve)
source("src/utils/physiology.R")
source("src/utils/plotting.R")
source("src/01_minimal_pbpk.R")

# --- TMDD Parameters (generic example) ---

get_tmdd_params <- function() {
  list(
    # Target parameters
    R0 = 10,            # baseline target concentration (nM)
    ksyn = 1.0,         # target synthesis rate (nM/h) = kdeg_T * R0
    kdeg_T = 0.1,       # target degradation rate (1/h)
    # Binding kinetics
    kon = 0.001,        # on-rate (1/nM/h)
    koff = 0.01,        # off-rate (1/h)  -> Kd = 10 nM
    kint = 0.05         # complex internalization rate (1/h)
  )
}

# --- Minimal PBPK + TMDD ODE ---

tmdd_pbpk_ode <- function(t, state, params) {

  C_plasma <- state["C_plasma"]
  C_tight  <- state["C_tight"]
  C_leaky  <- state["C_leaky"]
  R        <- state["R"]          # free target in plasma
  RC       <- state["RC"]         # drug-target complex in plasma

  with(params, {

    # --- TMDD in plasma ---
    J_bind <- kon * C_plasma * R
    J_unbind <- koff * RC
    J_int <- kint * RC
    J_syn <- ksyn
    J_deg_T <- kdeg_T * R

    # --- PBPK transport (Shah & Betts: pinocytosis from plasma) ---
    CL_net <- CLp * (1 - FR)

    conv_to_tight <- L_tight * (1 - sigma_v_tight) * C_plasma
    conv_to_leaky <- L_leaky * (1 - sigma_v_leaky) * C_plasma
    lymph_tight <- L_tight * (1 - sigma_l) * C_tight
    lymph_leaky <- L_leaky * (1 - sigma_l) * C_leaky

    # --- Free drug in plasma ---
    dC_plasma <- (1 / V_plasma) * (
      - CL_net * C_plasma
      + lymph_tight + lymph_leaky
      - conv_to_tight - conv_to_leaky
    ) - J_bind + J_unbind

    # --- Tight tissues ---
    dC_tight <- (1 / V_tight) * (
      conv_to_tight - lymph_tight
    )

    # --- Leaky tissues ---
    dC_leaky <- (1 / V_leaky) * (
      conv_to_leaky - lymph_leaky
    )

    # --- Free target ---
    dR <- J_syn - J_deg_T - J_bind + J_unbind

    # --- Drug-target complex ---
    dRC <- J_bind - J_unbind - J_int

    list(c(
      dC_plasma = unname(dC_plasma),
      dC_tight = unname(dC_tight),
      dC_leaky = unname(dC_leaky),
      dR = unname(dR),
      dRC = unname(dRC)
    ))
  })
}

# --- Simulation ---

run_tmdd_pbpk <- function(dose_mg, t_end = 42 * 24, dt = 1) {

  pbpk <- get_minimal_pbpk_params("typical_igg1")
  tmdd <- get_tmdd_params()
  params <- c(pbpk, tmdd)

  dose_nmol <- (dose_mg * 1e6) / params$MW
  C0_plasma <- dose_nmol / params$V_plasma

  state <- c(
    C_plasma = C0_plasma,
    C_tight = 0,
    C_leaky = 0,
    R = tmdd$R0,
    RC = 0
  )

  times <- seq(0, t_end, by = dt)

  out <- as.data.frame(ode(y = state, times = times, func = tmdd_pbpk_ode,
                           parms = params, method = "lsoda"))

  out$time_days <- out$time / 24
  out$C_total <- out$C_plasma + out$RC  # total drug (free + bound)
  out$target_suppression <- (1 - out$R / tmdd$R0) * 100

  out
}

# --- Main ---

if (sys.nframe() == 0) {

  cat("=", rep("=", 60), "\n", sep = "")
  cat("TMDD Extension - Nonlinear PK from Target Binding\n")
  cat("=", rep("=", 60), "\n\n", sep = "")

  # Compare three dose levels to show nonlinear PK
  doses <- c(0.1, 1, 10) * 70  # 0.1, 1, 10 mg/kg

  results <- list()
  for (d in doses) {
    res <- run_tmdd_pbpk(dose_mg = d, t_end = 60 * 24)
    res$dose_label <- sprintf("%.1f mg/kg", d / 70)
    results[[length(results) + 1]] <- res
  }

  df_all <- do.call(rbind, results)

  # Plot: nonlinear PK
  p <- ggplot(df_all, aes(x = time_days, y = C_total, color = dose_label)) +
    geom_line(linewidth = 0.8) +
    scale_y_log10() +
    theme_pbpk() +
    labs(title = "TMDD Effect on mAb PK: Dose-Dependent Nonlinearity",
         x = "Time (days)", y = "Total Plasma Concentration (nM)",
         color = "Dose")

  ggsave("results/tmdd_nonlinear_pk.png", p, width = 9, height = 5, dpi = 300)
  cat("Plot saved to results/tmdd_nonlinear_pk.png\n")

  # Plot: target suppression
  p2 <- ggplot(df_all, aes(x = time_days, y = target_suppression, color = dose_label)) +
    geom_line(linewidth = 0.8) +
    theme_pbpk() +
    labs(title = "Target Suppression vs Dose",
         x = "Time (days)", y = "Target Suppression (%)",
         color = "Dose")

  ggsave("results/tmdd_target_suppression.png", p2, width = 9, height = 5, dpi = 300)
  cat("Plot saved to results/tmdd_target_suppression.png\n")
}

# ============================================================================
# Minimal PBPK Model for Monoclonal Antibodies
# Based on: Shah & Betts, J Pharmacokinet Pharmacodyn (2012) 39:67-86
#
# Model structure:
#   - Plasma compartment
#   - Two lumped tissue groups: "tight" (muscle, skin, adipose) and "leaky"
#     (liver, spleen, gut, kidney, heart, lung)
#   - Endosomal compartment with FcRn-mediated recycling
#   - Lymphatic return
#
# Key processes:
#   1. Convective extravasation (plasma -> tissue interstitial via lymph)
#   2. Lymphatic drainage (tissue ISF -> lymph -> plasma)
#   3. Pinocytosis by vascular endothelium (from PLASMA, not ISF)
#   4. FcRn binding in acidic endosome (pH 6.0)
#   5. FcRn-bound mAb: recycled to plasma OR transcytosed to ISF
#   6. Unbound mAb in endosome: lysosomal degradation
#
# Units: concentration (nM), volume (L), time (hours), flow (L/h)
# ============================================================================

library(deSolve)
source("src/utils/physiology.R")
source("src/utils/plotting.R")

# --- Model Parameters ---

get_minimal_pbpk_params <- function(antibody = "typical_igg1") {

  phys <- get_physiology()

  if (antibody == "typical_igg1") {
    ab <- get_igg_params()
  } else if (antibody == "efgartigimod") {
    ab <- get_efgartigimod_params()
  }

  # Lumped tissue volumes (interstitial space only - where mAb distributes)
  # Tight tissues: muscle + skin + adipose + bone + brain + rest
  V_tight <- phys$V_muscle * phys$Fi_muscle +
             phys$V_skin * phys$Fi_skin +
             phys$V_adipose * phys$Fi_adipose +
             phys$V_bone * phys$Fi_bone +
             phys$V_brain * phys$Fi_brain +
             phys$V_rest * phys$Fi_rest

  # Leaky tissues: lung + liver + spleen + gut + kidney + heart
  V_leaky <- phys$V_lung * phys$Fi_lung +
             phys$V_liver * phys$Fi_liver +
             phys$V_spleen * phys$Fi_spleen +
             phys$V_gut * phys$Fi_gut +
             phys$V_kidney * phys$Fi_kidney +
             phys$V_heart * phys$Fi_heart

  # Total lymph flow: ~2.9 L/day = 0.12 L/h (Shah & Betts 2012)
  # Distributed proportionally to ISF volume (not blood flow)
  # Using blood flow for lung gives unrealistic lymph since Q_lung = CO
  L_total <- 0.12  # L/h
  f_tight <- V_tight / (V_tight + V_leaky)
  L_tight <- L_total * f_tight
  L_leaky <- L_total * (1 - f_tight)

  # Vascular reflection coefficients
  sigma_v_tight <- phys$sigma_v_tight
  sigma_v_leaky <- phys$sigma_v_leaky

  # Lymphatic reflection coefficient
  sigma_l <- phys$sigma_l

  # FcRn recycling parameters
  FR <- ab$FR                   # fraction recycled
  CLp <- ab$CLp_minimal         # pinocytosis CL calibrated for minimal model (L/h)
  k_deg <- CLp * (1 - FR)       # net clearance rate (L/h)

  list(
    # Volumes
    V_plasma = phys$V_plasma,
    V_tight = V_tight,
    V_leaky = V_leaky,
    # Flows
    L_tight = L_tight,
    L_leaky = L_leaky,
    L_total = L_total,
    # Transport
    sigma_v_tight = sigma_v_tight,
    sigma_v_leaky = sigma_v_leaky,
    sigma_l = sigma_l,
    # FcRn
    FR = FR,
    CLp = CLp,
    k_deg = k_deg,
    # Antibody
    MW = ab$MW
  )
}

# --- ODE System ---

minimal_pbpk_ode <- function(t, state, params) {

  C_plasma <- state["C_plasma"]
  C_tight  <- state["C_tight"]
  C_leaky  <- state["C_leaky"]

  with(params, {

    # === PINOCYTOSIS FROM PLASMA (Shah & Betts 2012, Eq 1) ===
    # Vascular endothelial cells pinocytose IgG from the luminal (plasma) side.
    # CLp = total pinocytosis clearance (L/h) from plasma.
    # After endosomal sorting by FcRn:
    #   - FR fraction recycled back to PLASMA (luminal side)
    #   - (1-FR) fraction degraded in lysosomes
    # Net plasma clearance = CLp * (1 - FR)
    # Tissue distribution is ONLY via convection (lymph-driven, 2-pore)

    CL_net <- CLp * (1 - FR)   # net clearance from plasma (L/h)

    # === CONVECTIVE TRANSPORT (2-pore, lymph-driven) ===
    # Extravasation: plasma -> ISF through large pores (Shah & Betts Eq 2-3)
    conv_to_tight <- L_tight * (1 - sigma_v_tight) * C_plasma
    conv_to_leaky <- L_leaky * (1 - sigma_v_leaky) * C_plasma

    # Lymphatic drainage: ISF -> lymph -> plasma
    lymph_tight <- L_tight * (1 - sigma_l) * C_tight
    lymph_leaky <- L_leaky * (1 - sigma_l) * C_leaky

    # === PLASMA ===
    dC_plasma <- (1 / V_plasma) * (
      - CL_net * C_plasma                # FcRn-mediated net clearance
      + lymph_tight + lymph_leaky        # lymphatic return
      - conv_to_tight - conv_to_leaky    # convective extravasation
    )

    # === TIGHT TISSUES (muscle, skin, adipose) ===
    dC_tight <- (1 / V_tight) * (
      conv_to_tight       # convective extravasation from plasma
      - lymph_tight       # lymphatic drainage to plasma
    )

    # === LEAKY TISSUES (lung, liver, spleen, gut, kidney, heart) ===
    dC_leaky <- (1 / V_leaky) * (
      conv_to_leaky       # convective extravasation from plasma
      - lymph_leaky       # lymphatic drainage to plasma
    )

    list(c(dC_plasma = unname(dC_plasma),
           dC_tight = unname(dC_tight),
           dC_leaky = unname(dC_leaky)))
  })
}

# --- Simulation ---

run_minimal_pbpk <- function(dose_mg, antibody = "typical_igg1",
                              t_end = 42 * 24, dt = 1,
                              route = "iv_bolus") {

  params <- get_minimal_pbpk_params(antibody)

  # Convert dose from mg to nmol
  dose_nmol <- (dose_mg * 1e6) / params$MW  # mg -> ng -> nmol

  # Initial concentration in plasma after IV bolus
  C0_plasma <- dose_nmol / params$V_plasma

  # Initial conditions
  state <- c(C_plasma = C0_plasma, C_tight = 0, C_leaky = 0)

  # Time vector (hours)
  times <- seq(0, t_end, by = dt)

  # Solve ODE
  out <- as.data.frame(ode(y = state, times = times, func = minimal_pbpk_ode,
                           parms = params, method = "lsoda"))

  # Add derived quantities
  out$time_days <- out$time / 24
  out$C_plasma_ug_mL <- out$C_plasma * params$MW / 1e6  # nM -> ug/mL

  # PK metrics
  cmax <- max(out$C_plasma)
  auc <- sum(diff(out$time) * (head(out$C_plasma, -1) + tail(out$C_plasma, -1)) / 2)

  # Terminal half-life (from last 50% of data in log-linear phase)
  n <- nrow(out)
  idx <- which(out$C_plasma > 0.01 * cmax)
  if (length(idx) > 10) {
    tail_idx <- tail(idx, round(length(idx) * 0.5))
    fit <- lm(log(C_plasma) ~ time, data = out[tail_idx, ])
    kel <- -coef(fit)["time"]
    t_half_h <- log(2) / kel
    t_half_d <- t_half_h / 24
  } else {
    t_half_d <- NA
  }

  list(
    data = out,
    params = params,
    pk_metrics = list(
      Cmax_nM = cmax,
      AUC_nM_h = auc,
      t_half_days = unname(t_half_d)
    )
  )
}

# --- Main execution ---

if (sys.nframe() == 0) {

  cat("=" , rep("=", 60), "\n", sep = "")
  cat("Minimal PBPK Model for Monoclonal Antibodies\n")
  cat("Shah & Betts 2012 framework\n")
  cat("=", rep("=", 60), "\n\n", sep = "")

  # Simulate typical IgG1: 10 mg/kg IV bolus in 70 kg human
  dose <- 10 * 70  # mg
  cat(sprintf("Simulating: Typical IgG1, %.0f mg IV bolus (10 mg/kg)\n", dose))

  result <- run_minimal_pbpk(dose_mg = dose, antibody = "typical_igg1",
                              t_end = 60 * 24)

  cat(sprintf("\nPK Metrics:\n"))
  cat(sprintf("  Cmax (plasma):   %.1f nM (%.1f ug/mL)\n",
              result$pk_metrics$Cmax_nM,
              result$pk_metrics$Cmax_nM * 150000 / 1e6))
  cat(sprintf("  Terminal t1/2:   %.1f days\n", result$pk_metrics$t_half_days))
  cat(sprintf("  AUC(0-inf):      %.0f nM*h\n", result$pk_metrics$AUC_nM_h))
  cat(sprintf("\n  Expected t1/2 for IgG1: ~21 days\n"))

  # Plot
  sim_df <- data.frame(
    time = result$data$time_days,
    concentration = result$data$C_plasma
  )

  p <- plot_pk_profile(sim_df, title = "Minimal PBPK - Typical IgG1 (10 mg/kg IV)",
                       xlab = "Time (days)", ylab = "Plasma Concentration (nM)")

  ggsave("results/pk_profile_minimal_igg1.png", p, width = 8, height = 5, dpi = 300)
  cat("\nPlot saved to results/pk_profile_minimal_igg1.png\n")
}

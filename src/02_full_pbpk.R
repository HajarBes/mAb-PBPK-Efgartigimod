# ============================================================================
# Full Tissue PBPK Model for Monoclonal Antibodies
# Based on: Shah & Betts 2012, Ng et al. 2006, Jones et al. 2019
#
# Extension of minimal model to individual tissue compartments.
# Each tissue has:
#   - Vascular sub-compartment (rapid equilibrium with plasma)
#   - Interstitial sub-compartment (extravascular distribution)
#
# Key processes:
#   - Pinocytosis from PLASMA (vascular endothelium luminal side)
#   - FcRn recycling back to plasma
#   - Convective extravasation (2-pore, lymph-driven)
#   - Lymphatic drainage back to plasma
#
# Tissues: lung, liver, spleen, gut, kidney, heart, muscle, skin,
#          adipose, bone, brain, rest
# ============================================================================

library(deSolve)
source("src/utils/physiology.R")
source("src/utils/plotting.R")

# Tissue names
TISSUES <- c("lung", "liver", "spleen", "gut", "kidney", "heart",
             "muscle", "skin", "adipose", "bone", "brain", "rest")

# Classify tissues by capillary type
LEAKY_TISSUES <- c("lung", "liver", "spleen", "gut", "kidney", "heart")
TIGHT_TISSUES <- c("muscle", "skin", "adipose", "bone", "brain", "rest")

# --- Build parameter set ---

get_full_pbpk_params <- function(antibody = "typical_igg1") {

  phys <- get_physiology()

  if (antibody == "typical_igg1") {
    ab <- get_igg_params()
  } else if (antibody == "efgartigimod") {
    ab <- get_efgartigimod_params()
  }

  # Extract tissue-specific parameters into named vectors
  tissues <- TISSUES

  V_tissue <- setNames(
    c(phys$V_lung, phys$V_liver, phys$V_spleen, phys$V_gut,
      phys$V_kidney, phys$V_heart, phys$V_muscle, phys$V_skin,
      phys$V_adipose, phys$V_bone, phys$V_brain, phys$V_rest),
    tissues
  )

  Fv <- setNames(
    c(phys$Fv_lung, phys$Fv_liver, phys$Fv_spleen, phys$Fv_gut,
      phys$Fv_kidney, phys$Fv_heart, phys$Fv_muscle, phys$Fv_skin,
      phys$Fv_adipose, phys$Fv_bone, phys$Fv_brain, phys$Fv_rest),
    tissues
  )

  Fi <- setNames(
    c(phys$Fi_lung, phys$Fi_liver, phys$Fi_spleen, phys$Fi_gut,
      phys$Fi_kidney, phys$Fi_heart, phys$Fi_muscle, phys$Fi_skin,
      phys$Fi_adipose, phys$Fi_bone, phys$Fi_brain, phys$Fi_rest),
    tissues
  )

  Q_tissue <- setNames(
    c(phys$Q_lung, phys$Q_liver, phys$Q_spleen, phys$Q_gut,
      phys$Q_kidney, phys$Q_heart, phys$Q_muscle, phys$Q_skin,
      phys$Q_adipose, phys$Q_bone, phys$Q_brain, phys$Q_rest),
    tissues
  )

  # Vascular reflection coefficients per tissue
  sigma_v <- setNames(
    ifelse(tissues %in% LEAKY_TISSUES, phys$sigma_v_leaky, phys$sigma_v_tight),
    tissues
  )

  # Interstitial volumes
  V_is <- V_tissue * Fi

  # Vascular volumes
  V_v <- V_tissue * Fv

  # Endosomal volumes
  V_e <- V_tissue * phys$V_endo

  # Lymph flows per tissue - distributed by ISF volume (not blood flow)
  # Total lymph: 0.12 L/h (Shah & Betts 2012)
  L_total <- 0.12
  L_tissue <- L_total * V_is / sum(V_is)

  # FcRn: pinocytosis from plasma (not ISF)
  # CLp is total pinocytosis clearance from plasma
  # Net clearance = CLp * (1 - FR)
  # Use CLp_minimal (calibrated for pinocytosis-from-plasma paradigm)

  list(
    tissues = tissues,
    n_tissues = length(tissues),
    V_plasma = phys$V_plasma,
    V_tissue = V_tissue,
    V_is = V_is,
    V_v = V_v,
    V_e = V_e,
    L_tissue = L_tissue,
    sigma_v = sigma_v,
    sigma_l = phys$sigma_l,
    CLp = ab$CLp_minimal,
    FR = ab$FR,
    MW = ab$MW
  )
}

# --- ODE System ---
# State vector: C_plasma, then for each tissue: C_is_i (interstitial concentration)
# Vascular sub-compartments assumed in rapid equilibrium with plasma
# Pinocytosis from PLASMA (vascular endothelium luminal side, per Shah & Betts 2012)
# FcRn recycling returns to plasma; unbound fraction degraded in lysosomes

full_pbpk_ode <- function(t, state, params) {

  with(params, {
    C_plasma <- state[1]
    C_is <- state[2:(n_tissues + 1)]
    names(C_is) <- tissues

    # Rates for each tissue
    dC_is <- numeric(n_tissues)
    names(dC_is) <- tissues

    total_lymph_return <- 0
    total_convection_out <- 0

    for (i in seq_along(tissues)) {
      ti <- tissues[i]

      # Convective extravasation: plasma -> interstitial
      J_conv <- L_tissue[ti] * (1 - sigma_v[ti]) * C_plasma

      # Lymphatic drainage: interstitial -> lymph -> plasma
      J_lymph <- L_tissue[ti] * (1 - sigma_l) * C_is[ti]

      # Tissue interstitial: only convection in, lymph out
      # (pinocytosis is from plasma, not from ISF)
      dC_is[i] <- (1 / V_is[ti]) * (
        J_conv    # in from plasma via convection
        - J_lymph # out via lymphatic drainage
      )

      # Accumulate plasma contributions
      total_lymph_return <- total_lymph_return + J_lymph
      total_convection_out <- total_convection_out + J_conv
    }

    # Net plasma clearance: pinocytosis from plasma, FR fraction recycled back
    CL_net <- CLp * (1 - FR)

    # Plasma balance
    dC_plasma <- (1 / V_plasma) * (
      - CL_net * C_plasma          # FcRn-mediated net clearance
      + total_lymph_return         # lymphatic return from all tissues
      - total_convection_out       # convective extravasation to all tissues
    )

    list(c(dC_plasma = unname(dC_plasma), unname(dC_is)))
  })
}

# --- Simulation ---

run_full_pbpk <- function(dose_mg, antibody = "typical_igg1",
                           t_end = 42 * 24, dt = 1) {

  params <- get_full_pbpk_params(antibody)

  dose_nmol <- (dose_mg * 1e6) / params$MW
  C0_plasma <- dose_nmol / params$V_plasma

  # Initial state: drug in plasma only
  state <- c(C_plasma = C0_plasma, setNames(rep(0, params$n_tissues), params$tissues))

  times <- seq(0, t_end, by = dt)

  out <- as.data.frame(ode(y = state, times = times, func = full_pbpk_ode,
                           parms = params, method = "lsoda"))

  # Rename columns
  colnames(out) <- c("time", "C_plasma", params$tissues)

  out$time_days <- out$time / 24
  out$C_plasma_ug_mL <- out$C_plasma * params$MW / 1e6

  # Terminal half-life
  cmax <- max(out$C_plasma)
  idx <- which(out$C_plasma > 0.01 * cmax)
  if (length(idx) > 10) {
    tail_idx <- tail(idx, round(length(idx) * 0.5))
    fit <- lm(log(C_plasma) ~ time, data = out[tail_idx, ])
    kel <- -coef(fit)["time"]
    t_half_d <- unname(log(2) / kel / 24)
  } else {
    t_half_d <- NA
  }

  list(
    data = out,
    params = params,
    pk_metrics = list(
      Cmax_nM = cmax,
      t_half_days = t_half_d
    )
  )
}

# --- Main ---

if (sys.nframe() == 0) {

  cat("=", rep("=", 60), "\n", sep = "")
  cat("Full Tissue PBPK Model for Monoclonal Antibodies\n")
  cat("12-tissue model with FcRn recycling\n")
  cat("=", rep("=", 60), "\n\n", sep = "")

  dose <- 10 * 70
  cat(sprintf("Simulating: Typical IgG1, %.0f mg IV bolus (10 mg/kg)\n", dose))

  result <- run_full_pbpk(dose_mg = dose, antibody = "typical_igg1", t_end = 60 * 24)

  cat(sprintf("\nPK Metrics:\n"))
  cat(sprintf("  Cmax (plasma):   %.1f nM\n", result$pk_metrics$Cmax_nM))
  cat(sprintf("  Terminal t1/2:   %.1f days\n", result$pk_metrics$t_half_days))

  # Tissue concentration plot
  p <- plot_tissue_profiles(
    result$data,
    tissues = c("C_plasma", "muscle", "skin", "liver", "spleen", "lung"),
    time_col = "time_days",
    title = "Full PBPK - IgG1 Tissue Distribution (10 mg/kg IV)"
  )

  ggsave("results/tissue_profiles_full_igg1.png", p, width = 10, height = 6, dpi = 300)
  cat("\nPlot saved to results/tissue_profiles_full_igg1.png\n")
}

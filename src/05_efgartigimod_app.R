# ============================================================================
# Efgartigimod (VYVGART) PBPK-PD Model
#
# Drug development question:
#   "Can PBPK predict the impact of FcRn blockade by efgartigimod
#    on endogenous IgG reduction, and how does FcRn binding affinity
#    engineering affect therapeutic efficacy?"
#
# Efgartigimod mechanism:
#   - Engineered human IgG1 Fc fragment (ABDEG technology)
#   - Binds FcRn with ~40x higher affinity than wild-type IgG at pH 6.0
#   - Blocks FcRn-mediated recycling of endogenous IgG
#   - IgG diverted to lysosomal degradation -> IgG levels drop
#   - Used in generalized myasthenia gravis (gMG), CIDP
#
# Model structure:
#   1. Efgartigimod PK: minimal PBPK (3-compartment + FcRn)
#   2. Endogenous IgG: turnover model with FcRn-dependent recycling
#   3. PD link: efgartigimod occupies FcRn -> reduces IgG recycling fraction
#
# Clinical data reference:
#   - Ulrichts et al. 2018, JCI (healthy volunteers, 10 mg/kg IV)
#   - Observed: ~75% IgG reduction at nadir, recovery over 8 weeks
#
# Units: nM, L, hours
# ============================================================================

library(deSolve)
source("src/utils/physiology.R")
source("src/utils/plotting.R")
source("src/03_fcrn_recycling.R")

# --- Combined PK-PD Parameters ---

get_efgartigimod_system_params <- function() {

  phys <- get_physiology()
  fcrn <- get_fcrn_params()

  # Efgartigimod PBPK parameters (Fc fragment)
  MW_efg <- 54000
  V_plasma <- phys$V_plasma

  # Lumped tissue volumes (same as minimal PBPK)
  # Tight tissues: muscle + skin + adipose + bone + brain + rest
  V_tight <- phys$V_muscle * phys$Fi_muscle +
             phys$V_skin * phys$Fi_skin +
             phys$V_adipose * phys$Fi_adipose +
             phys$V_bone * phys$Fi_bone +
             phys$V_brain * phys$Fi_brain +
             phys$V_rest * phys$Fi_rest

  V_leaky <- phys$V_lung * phys$Fi_lung +
             phys$V_liver * phys$Fi_liver +
             phys$V_spleen * phys$Fi_spleen +
             phys$V_gut * phys$Fi_gut +
             phys$V_kidney * phys$Fi_kidney +
             phys$V_heart * phys$Fi_heart

  # Total lymph flow: ~2.9 L/day = 0.12 L/h, distributed by ISF volume
  L_total <- 0.12
  f_tight <- V_tight / (V_tight + V_leaky)
  L_tight <- L_total * f_tight
  L_leaky <- L_total * (1 - f_tight)

  # Endogenous IgG parameters
  IgG_baseline <- 66667       # ~10 mg/mL = 66,667 nM
  IgG_kdeg_base <- 0.693 / (21 * 24)  # baseline degradation rate (1/h), t1/2 = 21 days
  IgG_ksyn <- IgG_kdeg_base * IgG_baseline  # synthesis rate (nM/h) to maintain SS

  list(
    # Efgartigimod PK
    MW_efg = MW_efg,
    V_plasma = V_plasma,
    V_tight = V_tight,
    V_leaky = V_leaky,
    L_tight = L_tight,
    L_leaky = L_leaky,
    sigma_v_tight = phys$sigma_v_tight,
    sigma_v_leaky = phys$sigma_v_leaky,
    sigma_l = phys$sigma_l,
    CLp_efg = 0.015,        # pinocytosis CL for Fc fragment (L/h)
    FR_efg = 0.0,           # efgartigimod NOT recycled (binds FcRn, gets degraded)
    CL_renal = 0.004,       # small renal CL (total CL ~0.019, t1/2 ~4.8 days)

    # FcRn competition
    FcRn_total = fcrn$FcRn_total,
    Kd_efg = 20,            # efgartigimod-FcRn Kd at pH 6.0 (nM) - ABDEG engineered
    Kd_igg = fcrn$Kd_pH6,   # wild-type IgG-FcRn Kd at pH 6.0 (nM)

    # Endogenous IgG turnover
    IgG_baseline = IgG_baseline,
    IgG_ksyn = IgG_ksyn,
    IgG_kdeg_base = IgG_kdeg_base,

    # FcRn recycling at baseline
    FR_igg_baseline = 0.715,

    # Degradation parameters
    k_recycle = fcrn$k_recycle,
    k_deg_lyso = fcrn$k_deg
  )
}

# --- Combined ODE: Efgartigimod PK + IgG PD ---

efgartigimod_ode <- function(t, state, params) {

  # Efgartigimod concentrations
  E_plasma <- state["E_plasma"]   # efgartigimod in plasma
  E_tight  <- state["E_tight"]    # efgartigimod in tight tissues
  E_leaky  <- state["E_leaky"]    # efgartigimod in leaky tissues

  # Endogenous IgG
  IgG <- state["IgG"]             # total endogenous IgG (plasma, nM)

  with(params, {

    # === Efgartigimod PK (minimal PBPK, pinocytosis from plasma) ===
    # Efgartigimod is NOT recycled by FcRn (FR=0) - all pinocytosed drug is degraded
    # Net CL = CLp_efg * (1 - 0) = CLp_efg (no recycling)
    CL_net_efg <- CLp_efg + CL_renal

    conv_tight_E <- L_tight * (1 - sigma_v_tight) * E_plasma
    conv_leaky_E <- L_leaky * (1 - sigma_v_leaky) * E_plasma
    lymph_tight_E <- L_tight * (1 - sigma_l) * E_tight
    lymph_leaky_E <- L_leaky * (1 - sigma_l) * E_leaky

    dE_plasma <- (1 / V_plasma) * (
      - CL_net_efg * E_plasma
      + lymph_tight_E + lymph_leaky_E
      - conv_tight_E - conv_leaky_E
    )

    dE_tight <- (1 / V_tight) * (
      conv_tight_E - lymph_tight_E
    )

    dE_leaky <- (1 / V_leaky) * (
      conv_leaky_E - lymph_leaky_E
    )

    # === FcRn Competition (Cheng-Prusoff competitive inhibition) ===
    # Efgartigimod and endogenous IgG compete for FcRn in the endosome.
    # Since IgG is present at ~66,667 nM (>>Kd_igg=800 nM), the effective
    # IC50 for FcRn blockade is much higher than Kd_efg alone:
    #   IC50_eff = Kd_efg * (1 + IgG/Kd_igg)  [Cheng-Prusoff equation]
    #
    # Fraction of FcRn bound to IgG:
    #   f_IgG = (IgG/Kd_igg) / (1 + IgG/Kd_igg + E/Kd_efg)
    #
    # Ratio to baseline (E=0):
    #   FR_eff/FR_base = (1 + IgG/Kd_igg) / (1 + IgG/Kd_igg + E/Kd_efg)

    IgG_over_Kd <- IgG / Kd_igg         # ~83 at baseline
    E_over_Kd <- E_plasma / Kd_efg

    FR_igg_effective <- FR_igg_baseline * (1 + IgG_over_Kd) /
                        (1 + IgG_over_Kd + E_over_Kd)

    # FcRn occupancy by efgartigimod
    frac_FcRn_efg <- E_over_Kd / (1 + IgG_over_Kd + E_over_Kd)

    # === IgG Turnover ===
    # IgG clearance depends on recycling fraction:
    #   High FR -> low clearance (most IgG recycled)
    #   Low FR (FcRn blocked) -> high clearance (IgG degraded)

    # Effective IgG half-life: t1/2 = ln2 / (k_deg * (1 - FR_effective))
    # Or equivalently: CL_igg = CL_base / (1 - FR_baseline) * (1 - FR_effective)
    k_elim_igg <- IgG_kdeg_base * (1 - FR_igg_effective) / (1 - FR_igg_baseline)

    dIgG <- IgG_ksyn - k_elim_igg * IgG

    list(c(
      dE_plasma = unname(dE_plasma),
      dE_tight = unname(dE_tight),
      dE_leaky = unname(dE_leaky),
      dIgG = unname(dIgG)
    ),
    # Output variables
    FR_effective = unname(FR_igg_effective),
    FcRn_occupancy_efg = unname(frac_FcRn_efg),
    IgG_percent = unname(IgG / IgG_baseline * 100)
    )
  })
}

# --- Dosing: multiple IV infusions ---

create_dosing_events <- function(dose_mg, MW, V_plasma,
                                  n_doses = 4, interval_days = 7) {
  # Weekly dosing regimen (like Phase 2 efgartigimod trial)
  dose_nmol <- (dose_mg * 1e6) / MW
  C_add <- dose_nmol / V_plasma

  data.frame(
    var = rep("E_plasma", n_doses),
    time = seq(0, by = interval_days * 24, length.out = n_doses),
    value = rep(C_add, n_doses),
    method = rep("add", n_doses)
  )
}

# --- Run Efgartigimod Simulation ---

run_efgartigimod <- function(dose_mg_kg = 10, n_doses = 4,
                              interval_days = 7, t_end = 84 * 24) {

  params <- get_efgartigimod_system_params()
  dose_mg <- dose_mg_kg * 70  # 70 kg patient

  events <- create_dosing_events(dose_mg, params$MW_efg, params$V_plasma,
                                  n_doses = n_doses,
                                  interval_days = interval_days)

  state <- c(
    E_plasma = 0,
    E_tight = 0,
    E_leaky = 0,
    IgG = params$IgG_baseline
  )

  times <- seq(0, t_end, by = 1)  # hourly resolution

  out <- as.data.frame(ode(y = state, times = times, func = efgartigimod_ode,
                           parms = params, method = "lsoda",
                           events = list(data = events)))

  out$time_days <- out$time / 24
  out$igg_percent_baseline <- out$IgG / params$IgG_baseline * 100
  out$E_plasma_ug_mL <- out$E_plasma * params$MW_efg / 1e6

  # Find nadir
  nadir_idx <- which.min(out$igg_percent_baseline)

  list(
    data = out,
    params = params,
    nadir = list(
      time_days = out$time_days[nadir_idx],
      igg_percent = out$igg_percent_baseline[nadir_idx]
    )
  )
}

# --- Main ---

if (sys.nframe() == 0) {

  cat("=", rep("=", 60), "\n", sep = "")
  cat("Efgartigimod PBPK-PD Model\n")
  cat("Drug Development Question: Impact of FcRn blockade on IgG reduction\n")
  cat("=", rep("=", 60), "\n\n", sep = "")

  # --- Simulation 1: Standard regimen (10 mg/kg weekly x4) ---
  cat("Simulation 1: Efgartigimod 10 mg/kg IV weekly x4\n")
  res <- run_efgartigimod(dose_mg_kg = 10, n_doses = 4)

  cat(sprintf("  IgG nadir: %.1f%% of baseline at day %.1f\n",
              res$nadir$igg_percent, res$nadir$time_days))
  cat(sprintf("  Clinical target: ~75%% reduction (25%% of baseline)\n"))
  cat(sprintf("  Ulrichts 2018 observed: ~75%% reduction\n\n"))

  # --- Simulation 2: Dose-response ---
  cat("Simulation 2: Dose-response comparison\n")
  doses <- c(2, 5, 10, 25)
  dose_results <- list()

  for (d in doses) {
    r <- run_efgartigimod(dose_mg_kg = d, n_doses = 4)
    r$data$dose_label <- sprintf("%d mg/kg", d)
    dose_results[[length(dose_results) + 1]] <- r
    cat(sprintf("  %2d mg/kg: IgG nadir = %.1f%% at day %.1f\n",
                d, r$nadir$igg_percent, r$nadir$time_days))
  }

  df_doses <- do.call(rbind, lapply(dose_results, function(x) x$data))

  # --- Plot 1: IgG reduction by dose ---
  p1 <- ggplot(df_doses, aes(x = time_days, y = igg_percent_baseline, color = dose_label)) +
    geom_line(linewidth = 0.8) +
    geom_hline(yintercept = 100, linetype = "dashed", color = "grey50") +
    geom_hline(yintercept = 25, linetype = "dotted", color = "red") +
    annotate("text", x = 70, y = 28, label = "75% reduction target",
             color = "red", size = 3) +
    # Dosing indicators
    geom_vline(xintercept = c(0, 7, 14, 21), linetype = "dotted",
               color = "grey70", linewidth = 0.3) +
    theme_pbpk() +
    labs(title = "Efgartigimod Dose-Response: IgG Reduction",
         subtitle = "Weekly IV x4 doses, 70 kg patient",
         x = "Time (days)", y = "IgG (% of Baseline)", color = "Dose")

  ggsave("results/efgartigimod_igg_reduction.png", p1, width = 10, height = 6, dpi = 300)

  # --- Plot 2: Efgartigimod PK ---
  p2 <- ggplot(df_doses, aes(x = time_days, y = E_plasma, color = dose_label)) +
    geom_line(linewidth = 0.8) +
    scale_y_log10() +
    geom_vline(xintercept = c(0, 7, 14, 21), linetype = "dotted",
               color = "grey70", linewidth = 0.3) +
    theme_pbpk() +
    labs(title = "Efgartigimod Plasma PK (Multiple Dose)",
         x = "Time (days)", y = "Plasma Concentration (nM)", color = "Dose")

  ggsave("results/efgartigimod_pk_profile.png", p2, width = 10, height = 6, dpi = 300)

  # --- Plot 3: FcRn occupancy ---
  res10 <- run_efgartigimod(dose_mg_kg = 10, n_doses = 4)

  p3 <- ggplot(res10$data, aes(x = time_days, y = FcRn_occupancy_efg * 100)) +
    geom_line(color = "#9B59B6", linewidth = 1) +
    geom_vline(xintercept = c(0, 7, 14, 21), linetype = "dotted",
               color = "grey70", linewidth = 0.3) +
    theme_pbpk() +
    labs(title = "FcRn Occupancy by Efgartigimod (10 mg/kg weekly x4)",
         x = "Time (days)", y = "FcRn Occupancy (%)")

  ggsave("results/fcrn_occupancy.png", p3, width = 9, height = 5, dpi = 300)

  cat("\nPlots saved to results/\n")
}

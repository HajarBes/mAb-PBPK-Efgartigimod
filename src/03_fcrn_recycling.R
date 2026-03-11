# ============================================================================
# Detailed FcRn Recycling Module
#
# This script models the endosomal FcRn recycling pathway in detail,
# replacing the simplified "fraction recycled" (FR) parameter with
# explicit binding kinetics.
#
# Biology:
#   1. IgG is taken up from interstitial fluid by pinocytosis
#   2. Endosome acidifies to pH ~6.0
#   3. At pH 6.0, IgG binds FcRn with high affinity (Kd ~100-1000 nM)
#   4. FcRn-IgG complex is recycled to cell surface
#   5. At pH 7.4, IgG releases from FcRn (low affinity at neutral pH)
#   6. Unbound IgG in endosome is trafficked to lysosome -> degraded
#
# Key parameters:
#   - FcRn_total: total FcRn concentration in endosome
#   - Kd_pH6: dissociation constant at pH 6.0
#   - k_recycle: rate of vesicle recycling to surface
#   - k_deg: rate of lysosomal degradation
#
# References:
#   - Garg & Balthasar 2007, J Pharmacokinet Pharmacodyn
#   - Pyzik et al. 2019, Nat Rev Immunol (FcRn biology review)
#   - Datta-Mannan 2019, Drug Metab Dispos
# ============================================================================

library(deSolve)
source("src/utils/physiology.R")

# --- FcRn Endosomal Kinetics Parameters ---

get_fcrn_params <- function() {
  list(
    FcRn_total = 49800,   # total FcRn concentration in endosome (nM)
                           # Garg & Balthasar 2007 estimate
    Kd_pH6 = 800,         # IgG-FcRn Kd at pH 6.0 (nM) - typical IgG1
    Kd_pH74 = 1e6,        # IgG-FcRn Kd at pH 7.4 (nM) - very weak
    kon_pH6 = 8.9e-4,     # on-rate at pH 6.0 (1/nM/h)
    koff_pH6 = 0.712,     # off-rate at pH 6.0 (1/h) = kon * Kd
    k_uptake = 0.0366,    # pinocytosis rate constant (1/h)
    k_recycle = 0.693,    # recycling rate from endosome (1/h), t1/2 ~1h
    k_deg = 0.347,        # lysosomal degradation rate (1/h), t1/2 ~2h
    k_sort = 1.0          # endosomal sorting rate (1/h) - how fast pH drops
  )
}

# --- Endosomal Sub-model ---
# For a single tissue compartment, track:
#   C_is: interstitial free IgG
#   C_endo_free: endosomal free IgG (not bound to FcRn)
#   C_endo_bound: endosomal IgG-FcRn complex
#   FcRn_free: free FcRn in endosome
#
# NOTE: This is a detailed mechanistic module for illustration/exploration.
# The main PBPK models (01, 05) use the steady-state FR approximation instead.
# Required params: all from get_fcrn_params() plus V_is, V_endo (tissue-specific)

fcrn_endosomal_ode <- function(t, state, params) {

  C_is       <- state["C_is"]
  C_endo_free <- state["C_endo_free"]
  C_endo_bound <- state["C_endo_bound"]

  with(params, {

    # Free FcRn = total - bound
    FcRn_free <- FcRn_total - C_endo_bound
    FcRn_free <- max(FcRn_free, 0)  # prevent negative

    # --- Endosomal binding at pH 6.0 ---
    J_bind <- kon_pH6 * C_endo_free * FcRn_free
    J_unbind <- koff_pH6 * C_endo_bound

    # --- Pinocytosis: plasma -> endosome ---
    # Note: in the full PBPK context, pinocytosis is from plasma (luminal side)
    # This sub-model uses C_is as input for standalone demonstration
    J_uptake <- k_uptake * C_is * V_is  # amount/h

    # --- Recycling: bound IgG:FcRn -> surface -> release ---
    J_recycle <- k_recycle * C_endo_bound

    # --- Lysosomal degradation: free endosomal IgG ---
    J_degrade <- k_deg * C_endo_free

    # --- ODEs ---
    dC_is <- (1 / V_is) * (
      - J_uptake
      + J_recycle * V_endo
    )

    dC_endo_free <- (1 / V_endo) * J_uptake - J_bind + J_unbind - J_degrade

    dC_endo_bound <- J_bind - J_unbind - J_recycle

    list(c(
      dC_is = unname(dC_is),
      dC_endo_free = unname(dC_endo_free),
      dC_endo_bound = unname(dC_endo_bound)
    ))
  })
}

# --- Compute Effective Recycling Fraction ---
# Given FcRn parameters and IgG concentration, compute the effective FR
# This replaces the constant FR in the minimal/full PBPK

compute_effective_FR <- function(C_igg, fcrn_params = get_fcrn_params()) {
  with(fcrn_params, {
    # At steady state in endosome:
    # Bound = FcRn_total * C_igg / (Kd_pH6 + C_igg)  (Langmuir isotherm)
    # FR = Bound / (Bound + Free) = Bound / C_igg_endo

    # For simplicity, assume endosomal IgG ~ pinocytosed ISF IgG
    bound <- FcRn_total * C_igg / (Kd_pH6 + C_igg)
    free <- C_igg - bound  # approximate

    # Fraction that gets recycled vs degraded
    fr_effective <- (k_recycle * bound) / (k_recycle * bound + k_deg * max(free, 0.01))
    fr_effective
  })
}

# --- Demo: FR vs IgG concentration ---
# Shows FcRn saturation - at high IgG, recycling fraction drops

if (sys.nframe() == 0) {

  cat("FcRn Recycling Module - Effective Recycling Fraction\n\n")

  fcrn <- get_fcrn_params()

  # IgG concentrations from 1 to 100,000 nM
  C_range <- 10^seq(0, 5, length.out = 200)
  FR_range <- sapply(C_range, compute_effective_FR, fcrn_params = fcrn)

  # Physiological IgG range: ~70,000 nM (10 mg/mL)
  cat(sprintf("Physiological IgG (~10 mg/mL = 66,667 nM): FR = %.3f\n",
              compute_effective_FR(66667, fcrn)))
  cat(sprintf("Low IgG (1,000 nM):                         FR = %.3f\n",
              compute_effective_FR(1000, fcrn)))
  cat(sprintf("Very low IgG (100 nM):                      FR = %.3f\n",
              compute_effective_FR(100, fcrn)))

  # Engineered high-affinity (efgartigimod-like, Kd = 20 nM)
  fcrn_eng <- fcrn
  fcrn_eng$Kd_pH6 <- 20
  cat(sprintf("\nEngineered Fc (Kd=20nM) at 66,667 nM:      FR = %.3f\n",
              compute_effective_FR(66667, fcrn_eng)))

  # Plot
  library(ggplot2)
  source("src/utils/plotting.R")

  df <- data.frame(
    C_igg = rep(C_range, 2),
    FR = c(FR_range,
           sapply(C_range, compute_effective_FR,
                  fcrn_params = modifyList(fcrn, list(Kd_pH6 = 20)))),
    type = rep(c("Wild-type IgG1 (Kd=800 nM)", "Engineered Fc (Kd=20 nM)"),
               each = length(C_range))
  )

  p <- ggplot(df, aes(x = C_igg, y = FR, color = type)) +
    geom_line(linewidth = 1) +
    scale_x_log10(labels = scales::comma) +
    geom_vline(xintercept = 66667, linetype = "dashed", color = "grey40") +
    annotate("text", x = 66667, y = 0.1, label = "Physiological IgG",
             angle = 90, vjust = -0.5, color = "grey40") +
    theme_pbpk() +
    labs(title = "FcRn Saturation: Effective Recycling Fraction vs IgG Concentration",
         x = "IgG Concentration (nM)",
         y = "Effective Recycling Fraction (FR)",
         color = "")

  ggsave("results/fcrn_saturation_curve.png", p, width = 9, height = 5, dpi = 300)
  cat("\nPlot saved to results/fcrn_saturation_curve.png\n")
}

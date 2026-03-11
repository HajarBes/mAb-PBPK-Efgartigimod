# ============================================================================
# Special Populations Analysis
#
# Drug development question:
#   "How do physiological differences in special populations affect
#    efgartigimod PK and IgG reduction?"
#
# Populations:
#   1. Pediatric (6-12 years) - allometric scaling
#   2. Hepatic impairment - altered FcRn expression, albumin, blood flow
#   3. Renal impairment - altered GFR for Fc fragment clearance
#   4. Elderly (>65 years) - altered volumes, flows, FcRn
#   5. Obese (BMI > 35) - altered Vd, cardiac output
#
# Regulatory context:
#   PBPK-based dose adjustment predictions are accepted by FDA/EMA
#   for special populations when clinical data is limited.
#
# References:
#   - Malik & Edginton 2018, CPT:PSP (pediatric mAb PBPK)
#   - Pan et al. 2020 (FcRn and aging)
# ============================================================================

source("src/05_efgartigimod_app.R")

# --- Population-specific parameter modifications ---

get_pediatric_params <- function(age_years = 9) {
  # Allometric scaling for a ~30 kg child (age 6-12)
  params <- get_efgartigimod_system_params()
  BW_child <- 30
  BW_adult <- 70

  # Allometric scaling exponents
  # Volumes: scale by BW^1.0
  # Flows: scale by BW^0.75
  # Clearances: scale by BW^0.75

  vol_scale <- BW_child / BW_adult
  flow_scale <- (BW_child / BW_adult)^0.75

  params$V_plasma <- params$V_plasma * vol_scale
  params$V_tight <- params$V_tight * vol_scale
  params$V_leaky <- params$V_leaky * vol_scale

  params$L_tight <- params$L_tight * flow_scale
  params$L_leaky <- params$L_leaky * flow_scale
  params$CLp_efg <- params$CLp_efg * flow_scale
  params$CL_renal <- params$CL_renal * flow_scale

  # IgG synthesis scales with body size
  params$IgG_ksyn <- params$IgG_ksyn * flow_scale

  # FcRn expression may differ in children (generally similar)
  # No strong evidence for change, keep same

  params$BW <- BW_child
  params$population <- "Pediatric (6-12y, 30 kg)"
  params
}

get_hepatic_impairment_params <- function(severity = "moderate") {
  params <- get_efgartigimod_system_params()

  if (severity == "mild") {
    # Mild: 80% liver blood flow, 90% FcRn
    flow_factor <- 0.80
    fcrn_factor <- 0.90
  } else if (severity == "moderate") {
    # Moderate: 60% liver blood flow, 70% FcRn
    flow_factor <- 0.60
    fcrn_factor <- 0.70
  } else {
    # Severe: 40% liver blood flow, 50% FcRn
    flow_factor <- 0.40
    fcrn_factor <- 0.50
  }

  # Liver is in leaky compartment - reduce its contribution
  params$L_leaky <- params$L_leaky * (0.7 + 0.3 * flow_factor)
  params$FcRn_total <- params$FcRn_total * fcrn_factor

  # Reduced FcRn -> reduced baseline recycling
  params$FR_igg_baseline <- params$FR_igg_baseline * fcrn_factor

  # IgG baseline may be lower in hepatic impairment
  params$IgG_baseline <- params$IgG_baseline * 0.85

  params$population <- paste0("Hepatic impairment (", severity, ")")
  params
}

get_renal_impairment_params <- function(gfr_fraction = 0.3) {
  params <- get_efgartigimod_system_params()

  # Fc fragment has some renal clearance
  params$CL_renal <- params$CL_renal * gfr_fraction

  # Renal impairment doesn't strongly affect mAb/Fc PK
  # but GFR reduction affects small Fc fragments
  params$population <- sprintf("Renal impairment (GFR %.0f%%)", gfr_fraction * 100)
  params
}

get_elderly_params <- function(age = 75) {
  params <- get_efgartigimod_system_params()

  # Reduced cardiac output (~20%), reduced lean mass
  flow_scale <- 0.80
  params$L_tight <- params$L_tight * flow_scale
  params$L_leaky <- params$L_leaky * flow_scale

  # Slightly reduced Vd
  params$V_plasma <- params$V_plasma * 0.95
  params$V_tight <- params$V_tight * 0.90  # less muscle
  params$V_leaky <- params$V_leaky * 0.95

  # FcRn expression may slightly decrease with age
  params$FcRn_total <- params$FcRn_total * 0.90

  params$population <- sprintf("Elderly (%d years)", age)
  params
}

get_obese_params <- function(bmi = 40) {
  params <- get_efgartigimod_system_params()

  BW_obese <- 120  # ~120 kg for BMI 40
  BW_ref <- 70

  # Increased Vd (primarily adipose expansion)
  params$V_plasma <- params$V_plasma * (BW_obese / BW_ref)^0.9
  params$V_tight <- params$V_tight * (BW_obese / BW_ref)^1.2  # adipose-heavy
  params$V_leaky <- params$V_leaky * (BW_obese / BW_ref)^0.8

  # Increased cardiac output
  flow_scale <- (BW_obese / BW_ref)^0.7
  params$L_tight <- params$L_tight * flow_scale
  params$L_leaky <- params$L_leaky * flow_scale
  params$CLp_efg <- params$CLp_efg * flow_scale

  # IgG baseline slightly higher in obesity
  params$IgG_baseline <- params$IgG_baseline * 1.1

  params$BW <- BW_obese
  params$population <- sprintf("Obese (BMI %d, %d kg)", bmi, BW_obese)
  params
}

# --- Run All Populations ---

if (sys.nframe() == 0) {

  cat("=", rep("=", 60), "\n", sep = "")
  cat("Special Populations Analysis - Efgartigimod\n")
  cat("=", rep("=", 60), "\n\n")

  # Define populations and dose
  populations <- list(
    list(name = "Adult (reference)", params_fn = get_efgartigimod_system_params,
         dose_mg_kg = 10),
    list(name = "Pediatric (30 kg)", params_fn = get_pediatric_params,
         dose_mg_kg = 10),
    list(name = "Hepatic (moderate)", params_fn = get_hepatic_impairment_params,
         dose_mg_kg = 10),
    list(name = "Renal (GFR 30%)", params_fn = function() get_renal_impairment_params(0.3),
         dose_mg_kg = 10),
    list(name = "Elderly (75y)", params_fn = get_elderly_params,
         dose_mg_kg = 10),
    list(name = "Obese (BMI 40)", params_fn = get_obese_params,
         dose_mg_kg = 10)
  )

  all_results <- list()

  for (pop in populations) {
    params <- pop$params_fn()
    BW <- ifelse(!is.null(params$BW), params$BW, 70)
    dose_mg <- pop$dose_mg_kg * BW

    res <- run_efgartigimod_with_params(params, dose_mg_kg = dose_mg / 70)
    res$data$population <- pop$name

    all_results[[pop$name]] <- res

    cat(sprintf("  %-25s: IgG nadir = %5.1f%% at day %.1f\n",
                pop$name, res$nadir$igg_percent, res$nadir$time_days))
  }

  # Combined plot
  df_all <- do.call(rbind, lapply(all_results, function(x) x$data))

  p <- ggplot(df_all, aes(x = time_days, y = igg_percent_baseline, color = population)) +
    geom_line(linewidth = 0.8) +
    geom_hline(yintercept = 100, linetype = "dashed", color = "grey50") +
    geom_vline(xintercept = c(0, 7, 14, 21), linetype = "dotted",
               color = "grey70", linewidth = 0.3) +
    theme_pbpk() +
    labs(title = "Efgartigimod IgG Reduction Across Special Populations",
         subtitle = "10 mg/kg IV weekly x4 | PBPK-based predictions",
         x = "Time (days)", y = "IgG (% of Baseline)", color = "Population")

  ggsave("results/special_populations.png", p, width = 11, height = 6, dpi = 300)
  cat("\nPlot saved to results/special_populations.png\n")

  # Summary table
  cat("\n--- Summary Table ---\n")
  cat(sprintf("%-25s  %12s  %12s\n", "Population", "IgG Nadir(%)", "Nadir Day"))
  cat(rep("-", 55), "\n", sep = "")
  for (pop in populations) {
    r <- all_results[[pop$name]]
    cat(sprintf("%-25s  %12.1f  %12.1f\n",
                pop$name, r$nadir$igg_percent, r$nadir$time_days))
  }
}

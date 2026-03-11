# Human Physiological Parameters for mAb PBPK
# Sources: Shah & Betts 2012, ICRP 2002, Brown et al. 1997
# 70 kg reference human

get_physiology <- function() {

  # Body weight (kg)
  BW <- 70

  # Tissue volumes (L) - Shah & Betts 2012 Table 1
  V_plasma    <- 3.13
  V_lung      <- 0.76
  V_liver     <- 1.69
  V_spleen    <- 0.15
  V_gut       <- 1.20
  V_kidney    <- 0.28
  V_heart     <- 0.31
  V_muscle    <- 25.00
  V_skin      <- 3.30
  V_adipose   <- 10.00
  V_bone      <- 3.80
  V_brain     <- 1.40
  V_rest      <- 4.58  # remaining tissues

  # Vascular fractions (dimensionless) - Shah & Betts 2012
  Fv_lung     <- 0.55
  Fv_liver    <- 0.21
  Fv_spleen   <- 0.17
  Fv_gut      <- 0.05
  Fv_kidney   <- 0.07
  Fv_heart    <- 0.26
  Fv_muscle   <- 0.04
  Fv_skin     <- 0.05
  Fv_adipose  <- 0.05
  Fv_bone     <- 0.04
  Fv_brain    <- 0.03
  Fv_rest     <- 0.05

  # Interstitial fractions (dimensionless)
  Fi_lung     <- 0.30
  Fi_liver    <- 0.16
  Fi_spleen   <- 0.15
  Fi_gut      <- 0.14
  Fi_kidney   <- 0.13
  Fi_heart    <- 0.14
  Fi_muscle   <- 0.12
  Fi_skin     <- 0.38
  Fi_adipose  <- 0.14
  Fi_bone     <- 0.10
  Fi_brain    <- 0.09
  Fi_rest     <- 0.15

  # Plasma flows (L/h) - Shah & Betts 2012 (from cardiac output)
  CO <- 362.4  # cardiac output (L/h) = plasma flow to lungs

  Q_lung    <- CO
  Q_liver   <- 13.80  # hepatic artery only (portal from gut/spleen handled separately)
  Q_spleen  <- 5.76
  Q_gut     <- 43.56
  Q_kidney  <- 36.00
  Q_heart   <- 7.20
  Q_muscle  <- 27.72
  Q_skin    <- 18.36
  Q_adipose <- 9.36
  Q_bone    <- 5.04
  Q_brain   <- 21.60
  Q_rest    <- CO - (Q_liver + Q_spleen + Q_gut + Q_kidney +
                      Q_heart + Q_muscle + Q_skin + Q_adipose +
                      Q_bone + Q_brain)

  # Lymph flow (L/h) - approximately 0.2% of plasma flow per tissue
  L_total <- 0.002 * CO  # total lymph flow

  # Vascular reflection coefficients (sigma_v)
  # Leaky tissues: lower sigma -> more convective transport
  # Tight tissues: higher sigma -> less transport
  sigma_v_tight <- 0.95   # muscle, skin, adipose, bone, brain
  sigma_v_leaky <- 0.85   # liver, spleen, gut, kidney, heart, lung

  # Lymphatic reflection coefficient
  sigma_l <- 0.2  # for all tissues

  # Endosomal parameters (FcRn recycling)
  V_endo <- 0.005  # endosomal volume per L tissue (L/L)

  list(
    BW = BW,
    # Volumes
    V_plasma = V_plasma,
    V_lung = V_lung, V_liver = V_liver, V_spleen = V_spleen,
    V_gut = V_gut, V_kidney = V_kidney, V_heart = V_heart,
    V_muscle = V_muscle, V_skin = V_skin, V_adipose = V_adipose,
    V_bone = V_bone, V_brain = V_brain, V_rest = V_rest,
    # Vascular fractions
    Fv_lung = Fv_lung, Fv_liver = Fv_liver, Fv_spleen = Fv_spleen,
    Fv_gut = Fv_gut, Fv_kidney = Fv_kidney, Fv_heart = Fv_heart,
    Fv_muscle = Fv_muscle, Fv_skin = Fv_skin, Fv_adipose = Fv_adipose,
    Fv_bone = Fv_bone, Fv_brain = Fv_brain, Fv_rest = Fv_rest,
    # Interstitial fractions
    Fi_lung = Fi_lung, Fi_liver = Fi_liver, Fi_spleen = Fi_spleen,
    Fi_gut = Fi_gut, Fi_kidney = Fi_kidney, Fi_heart = Fi_heart,
    Fi_muscle = Fi_muscle, Fi_skin = Fi_skin, Fi_adipose = Fi_adipose,
    Fi_bone = Fi_bone, Fi_brain = Fi_brain, Fi_rest = Fi_rest,
    # Plasma flows
    CO = CO,
    Q_lung = Q_lung, Q_liver = Q_liver, Q_spleen = Q_spleen,
    Q_gut = Q_gut, Q_kidney = Q_kidney, Q_heart = Q_heart,
    Q_muscle = Q_muscle, Q_skin = Q_skin, Q_adipose = Q_adipose,
    Q_bone = Q_bone, Q_brain = Q_brain, Q_rest = Q_rest,
    # Lymph
    L_total = L_total,
    # Reflection coefficients
    sigma_v_tight = sigma_v_tight,
    sigma_v_leaky = sigma_v_leaky,
    sigma_l = sigma_l,
    # Endosomal
    V_endo = V_endo
  )
}

# Antibody-specific parameters for a typical IgG1
get_igg_params <- function() {
  list(
    MW = 150000,           # molecular weight (Da)
    FcRn_Kd = 800,         # FcRn binding Kd at pH 6.0 (nM) - typical IgG1
    FcRn_kon = 8.9e5,      # FcRn on-rate (1/M/s)
    FcRn_koff = 7.1e-4,    # FcRn off-rate (1/s)
    FR = 0.715,            # fraction recycled via FcRn (Shah & Betts 2012)
    CLp = 0.0366,          # pinocytosis CL for full-tissue PBPK (Shah & Betts 2012)
    CLp_minimal = 0.017,   # calibrated for minimal 3-compartment PBPK to match t1/2 ~21d
    half_life_typical = 21 # typical IgG1 half-life (days), for validation
  )
}

# Efgartigimod-specific parameters
get_efgartigimod_params <- function() {
  list(
    MW = 54000,            # Fc fragment, not full IgG (~54 kDa)
    FcRn_Kd_pH6 = 20,     # Enhanced FcRn binding at pH 6.0 (nM) - engineered ABDEG
    FcRn_Kd_pH7 = 1000,   # Binding at pH 7.4 (nM) - releases at neutral pH
    FR = 0.0,              # NOT recycled - it blocks FcRn and gets degraded
    CLp = 0.017,           # pinocytosis CL calibrated for minimal PBPK
    CLp_minimal = 0.017,   # same for minimal model
    IC50_FcRn = 50,        # IC50 for blocking IgG-FcRn interaction (nM, approximate)
    half_life = 4.83       # observed terminal t1/2 ~4.8 days (Ulrichts 2018)
  )
}

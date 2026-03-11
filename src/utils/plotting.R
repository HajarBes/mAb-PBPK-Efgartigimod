# Publication-quality plotting utilities for mAb PBPK
library(ggplot2)

theme_pbpk <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey95"),
      legend.position = "bottom",
      plot.title = element_text(face = "bold", size = base_size + 2),
      axis.title = element_text(size = base_size),
      legend.text = element_text(size = base_size - 1)
    )
}

# Standard PK profile plot (concentration vs time)
plot_pk_profile <- function(sim_data, obs_data = NULL,
                            x = "time", y = "concentration",
                            group = NULL, log_y = TRUE,
                            title = "", xlab = "Time (days)",
                            ylab = "Concentration (nM)") {

  p <- ggplot(sim_data, aes(x = .data[[x]], y = .data[[y]])) +
    theme_pbpk()

  if (!is.null(group)) {
    p <- p + geom_line(aes(color = .data[[group]], linetype = .data[[group]]), linewidth = 0.8)
  } else {
    p <- p + geom_line(color = "#2C3E50", linewidth = 0.8)
  }

  if (!is.null(obs_data)) {
    p <- p + geom_point(data = obs_data,
                        aes(x = .data[[x]], y = .data[[y]]),
                        color = "red", size = 2.5, shape = 16)
  }

  if (log_y) {
    p <- p + scale_y_log10()
  }

  p + labs(title = title, x = xlab, y = ylab)
}

# Multi-tissue concentration plot
plot_tissue_profiles <- function(sim_data, tissues, time_col = "time",
                                 title = "Tissue Concentration-Time Profiles") {

  df_long <- tidyr::pivot_longer(sim_data, cols = all_of(tissues),
                                  names_to = "tissue", values_to = "concentration")

  ggplot(df_long, aes(x = .data[[time_col]], y = concentration, color = tissue)) +
    geom_line(linewidth = 0.7) +
    scale_y_log10() +
    theme_pbpk() +
    labs(title = title, x = "Time (days)", y = "Concentration (nM)", color = "Tissue")
}

# IgG reduction plot (specific to efgartigimod PD)
plot_igg_reduction <- function(sim_data, obs_data = NULL,
                                title = "IgG Reduction Following Efgartigimod Administration") {

  p <- ggplot(sim_data, aes(x = time, y = igg_percent_baseline)) +
    geom_line(color = "#E74C3C", linewidth = 1) +
    geom_hline(yintercept = 100, linetype = "dashed", color = "grey50") +
    theme_pbpk() +
    labs(title = title,
         x = "Time (days)",
         y = "IgG (% of Baseline)")

  if (!is.null(obs_data)) {
    p <- p + geom_point(data = obs_data,
                        aes(x = time, y = igg_percent_baseline),
                        color = "black", size = 2.5, shape = 16) +
      geom_errorbar(data = obs_data,
                    aes(x = time, ymin = igg_lower, ymax = igg_upper),
                    width = 0.5, color = "black")
  }
  p
}

# Sensitivity tornado plot
plot_tornado <- function(sa_results, title = "Parameter Sensitivity Analysis") {

  sa_results <- sa_results[order(abs(sa_results$sensitivity)), ]
  sa_results$parameter <- factor(sa_results$parameter,
                                  levels = sa_results$parameter)

  ggplot(sa_results, aes(x = sensitivity, y = parameter)) +
    geom_col(aes(fill = sensitivity > 0), width = 0.7) +
    scale_fill_manual(values = c("TRUE" = "#3498DB", "FALSE" = "#E74C3C"),
                      labels = c("TRUE" = "Increase", "FALSE" = "Decrease"),
                      name = "Effect on AUC") +
    geom_vline(xintercept = 0, linewidth = 0.5) +
    theme_pbpk() +
    labs(title = title, x = "Normalized Sensitivity Index", y = "")
}

# Model comparison overlay
plot_model_comparison <- function(data_list, labels, obs_data = NULL,
                                  title = "Model Comparison") {

  df_all <- do.call(rbind, lapply(seq_along(data_list), function(i) {
    d <- data_list[[i]]
    d$model <- labels[i]
    d
  }))

  p <- ggplot(df_all, aes(x = time, y = concentration, color = model, linetype = model)) +
    geom_line(linewidth = 0.8) +
    scale_y_log10() +
    theme_pbpk() +
    labs(title = title, x = "Time (days)", y = "Concentration (nM)")

  if (!is.null(obs_data)) {
    p <- p + geom_point(data = obs_data,
                        aes(x = time, y = concentration),
                        color = "black", size = 2.5, shape = 16, inherit.aes = FALSE)
  }
  p
}

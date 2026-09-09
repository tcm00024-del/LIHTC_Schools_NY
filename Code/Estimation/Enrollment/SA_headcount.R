rm(list = ls())
library(here)
source(here("code", "estimation", "enrollment", "enrollment_estimation_helpers.R"))

library(ggplot2)

base_path <- here("data", "processed")
samples <- samples_enrollment(base_path, "headcount_panel")
event_window <- c(-3, 3); ref_period <- -1

for (sample_key in names(samples)) {
  cfg <- samples[[sample_key]]
  dt <- fread(cfg$path)
  message("\n=== SA Headcount: ", cfg$label, " ===")

  plot_data <- rbindlist(lapply(names(event_time_vars), function(type_name) {
    res <- run_sa_generic(dt, "K12_TOTAL", event_time_vars[[type_name]], never_treated_vars[[type_name]],
                          transform = "log", event_window = event_window, ref_period = ref_period)
    if (is.null(res)) return(NULL)
    res[, type := type_name]; res
  }), fill = TRUE)

  if (is.null(plot_data) || nrow(plot_data) == 0) { message("  No results, skipping."); next }
  plot_data[, type := factor(type, levels = names(type_colors))]

  p <- ggplot(plot_data, aes(x = event_time, y = est, color = type, shape = type)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_vline(xintercept = -0.5, linetype = "dotted", color = "grey50") +
    geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0.15, position = position_dodge(0.3)) +
    geom_point(size = 2.2, position = position_dodge(0.3)) +
    scale_color_manual(values = type_colors, labels = type_labels, name = "Development type") +
    scale_shape_manual(values = c(ALL = 16, FAMILY = 17, ELDERLY = 15), labels = type_labels, name = "Development type") +
    labs(x = "Years relative to first LIHTC opening", y = "Effect on log K-12 enrollment",
         caption = sprintf("Sun-Abraham estimates, %s. Reference period: %d.", cfg$label, ref_period)) +
    theme_bw(base_size = 11) + theme(legend.position = "bottom")

  out_file <- file.path(base_path, sprintf("SA_headcount_%s.png", sample_key))
  ggsave(out_file, p, width = 6.5, height = 4.2, dpi = 300, bg = "white")
  message("Saved: ", out_file)
}

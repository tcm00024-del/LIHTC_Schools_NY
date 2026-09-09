# ============================================================
# Sun_Ab_Estimates.R
# ============================================================
rm(list = ls())

library(data.table)
library(fixest)
library(ggplot2)
library(here)

base_path <- here("data", "processed", "Subgroups")

yname  <- "mean_scale_score_weighted"
tname  <- "YEAR"
idname <- "BEDSCODE"

event_window <- c(-3, 3)
ref_period   <- -1

# Event-time / never-treated columns produced by Step 5 -- note we do NOT
# have a raw FIRST_PIS_{TYPE} column in this pipeline's exposure table
# (only EVENT_TIME_PIS_{TYPE} and NEVER_TREATED_PIS_{TYPE} survived the
# pivot). Cohort is reconstructed below as YEAR - EVENT_TIME_PIS_{TYPE},
# which is exact since EVENT_TIME_PIS = YEAR - FIRST_PIS by construction.
event_time_vars <- list(
  ALL     = "EVENT_TIME_PIS_ALL",
  FAMILY  = "EVENT_TIME_PIS_FAMILY",
  ELDERLY = "EVENT_TIME_PIS_ELDERLY"
)
never_treated_vars <- list(
  ALL     = "NEVER_TREATED_PIS_ALL",
  FAMILY  = "NEVER_TREATED_PIS_FAMILY",
  ELDERLY = "NEVER_TREATED_PIS_ELDERLY"
)

type_colors <- c(ALL = "#2166ac", FAMILY = "#d6604d", ELDERLY = "#4dac26")
type_labels <- c(ALL = "All", FAMILY = "Family", ELDERLY = "Non-Family")

pretty_subgroup <- function(x) {
  if (grepl("All_Students",                   x)) return("All Students")
  if (grepl("General_Education",              x)) return("General Education")
  if (grepl("Not_Economically_Disadvantaged", x)) return("Not Economically Disadvantaged")
  if (grepl("Economically_Disadvantaged",     x)) return("Economically Disadvantaged")
  if (grepl("Female",                         x)) return("Female")
  if (grepl("Male",                           x)) return("Male")
  if (grepl("Black",                          x)) return("Black")
  if (grepl("Hispanic",                       x)) return("Hispanic")
  if (grepl("White",                          x)) return("White")
  return(x)
}
file_slug <- function(x) gsub("[^A-Za-z0-9]", "_", x)

desired_order <- c(
  "All Students", "General Education",
  "Economically Disadvantaged", "Not Economically Disadvantaged",
  "Female", "Male", "Black", "Hispanic", "White"
)

samples <- list(
  full   = list(path = file.path(base_path, "Full"),   pattern = "^full_subgroup_.*_collapsed\\.csv$",    label = "Full Sample"),
  no_nyc = list(path = file.path(base_path, "No_NYC"), pattern = "^no_nyc_subgroup_.*_collapsed\\.csv$",  label = "No-NYC Sample")
)

# -----------------------------
# Sun-Abraham runner for one treatment type
# -----------------------------
run_sa_one <- function(dt, event_time_var, never_treated_var, yname, tname, idname,
                       event_window, ref_period) {
  tryCatch({
    dt <- copy(dt)

    if (!(event_time_var %in% names(dt)) || !(never_treated_var %in% names(dt))) {
      message("    Missing column(s): ", event_time_var, " / ", never_treated_var)
      return(NULL)
    }

    # Reconstruct first-treatment year, then cohort (Inf for never-treated)
    dt[, first_treat_year := as.numeric(get(tname)) - as.numeric(get(event_time_var))]
    dt[, cohort := fifelse(get(never_treated_var) == 1, Inf, first_treat_year)]

    cat("    Cohort distribution (finite only):\n")
    print(table(dt$cohort[is.finite(dt$cohort)], useNA = "always"))

    fml <- as.formula(sprintf(
      "%s ~ sunab(cohort, %s, ref.p = %d) | %s + %s",
      yname, tname, ref_period, idname, tname
    ))

    m <- feols(fml, data = dt, cluster = as.formula(paste0("~", idname)))

    coefs <- as.data.table(coeftable(m, agg = "period"), keep.rownames = "term")
    setnames(coefs, c("term", "est", "se", "t", "p"))

    coefs[, event_time := as.integer(sub(".*::([-0-9]+)$", "\\1", term))]
    coefs <- coefs[!is.na(event_time)]
    coefs <- coefs[event_time >= event_window[1] & event_time <= event_window[2]]

    ref_row <- data.table(term = paste0("ref::", ref_period), est = 0, se = 0,
                          t = NA_real_, p = NA_real_, event_time = ref_period)
    coefs <- rbindlist(list(coefs, ref_row), fill = TRUE)
    setorder(coefs, event_time)

    coefs[, `:=`(ci_lo = est - 1.96 * se, ci_hi = est + 1.96 * se)]
    coefs

  }, error = function(e) {
    message("    SA failed for ", event_time_var, ": ", e$message)
    NULL
  })
}

for (sample_key in names(samples)) {
  cfg <- samples[[sample_key]]
  message("\n=== Running Sun-Abraham: ", cfg$label, " ===")

  out_path <- file.path(cfg$path, "EventStudy_Plots")
  dir.create(out_path, showWarnings = FALSE, recursive = TRUE)

  files <- list.files(cfg$path, pattern = cfg$pattern, full.names = TRUE)
  if (length(files) == 0) stop("No files matched pattern in: ", cfg$path)
  message("Found ", length(files), " subgroup files.")

  for (file_path in files) {
    subgroup_label <- pretty_subgroup(basename(file_path))
    if (!(subgroup_label %in% desired_order)) next

    cat("Processing:", subgroup_label, "\n")

    dt <- tryCatch(fread(file_path), error = function(e) { message("  Failed to load: ", e$message); NULL })
    if (is.null(dt)) next

    dt[, (idname) := as.character(get(idname))]
    dt[, (tname)  := as.integer(get(tname))]

    plot_data <- rbindlist(
      lapply(names(event_time_vars), function(type_name) {
        cat("  Running SA for", type_name, "\n")
        res <- run_sa_one(
          dt, event_time_vars[[type_name]], never_treated_vars[[type_name]],
          yname, tname, idname, event_window, ref_period
        )
        if (is.null(res)) return(NULL)
        res[, type := type_name]
        res
      }),
      fill = TRUE
    )

    if (is.null(plot_data) || nrow(plot_data) == 0) {
      message("  No results for ", subgroup_label, ", skipping plot.")
      next
    }

    plot_data[, type := factor(type, levels = names(type_colors))]
    dodge_w <- 0.3

    p <- ggplot(plot_data, aes(x = event_time, y = est, color = type, shape = type)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.4) +
      geom_vline(xintercept = -0.5, linetype = "dotted", color = "grey50", linewidth = 0.4) +
      geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0.15, linewidth = 0.6,
                    position = position_dodge(width = dodge_w)) +
      geom_point(size = 2.2, position = position_dodge(width = dodge_w)) +
      scale_color_manual(values = type_colors, labels = type_labels, name = "Development type") +
      scale_shape_manual(values = c(ALL = 16, FAMILY = 17, ELDERLY = 15), labels = type_labels, name = "Development type") +
      scale_x_continuous(breaks = seq(event_window[1], event_window[2], by = 1),
                         labels = seq(event_window[1], event_window[2], by = 1)) +
      labs(
        x = "Years relative to first LIHTC opening",
        y = "Effect on mean weighted test score (scale score points)",
        caption = sprintf(
          "Sun-Abraham estimates, %s. Reference period: %d. 95%% confidence intervals shown.\nDevelopment types identified via Atkins-O'Regan (2014). School and year FEs included.",
          cfg$label, ref_period
        )
      ) +
      theme_bw(base_size = 11) +
      theme(
        legend.position = "bottom", legend.title = element_text(size = 9), legend.text = element_text(size = 9),
        panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        plot.caption = element_text(size = 7, color = "grey40", hjust = 0, margin = margin(t = 6)),
        axis.title = element_text(size = 10)
      )

    out_file <- file.path(out_path, sprintf("SA_EventStudy_%s.png", file_slug(subgroup_label)))
    ggsave(filename = out_file, plot = p, width = 6.5, height = 4.2, dpi = 300, bg = "white")
    cat("  Saved:", basename(out_file), "\n")
  }
}

message("Done.")

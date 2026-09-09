rm(list = ls())
library(here)
source(here("code", "estimation", "enrollment", "enrollment_estimation_helpers.R"))


Sys.setenv(OMP_NUM_THREADS = 1); Sys.setenv(MKL_NUM_THREADS = 1); options(mc.cores = 1)

base_path <- here("data", "processed")
samples <- samples_enrollment(base_path, "demographic_panel")

outcome_vars <- c(
  "NUM_ECDIS", "PER_ECDIS",
  "NUM_FREE_LUNCH", "PER_FREE_LUNCH",
  "NUM_SWD",   "PER_SWD",
  "NUM_BLACK", "PER_BLACK",
  "NUM_HISP",  "PER_HISP",
  "NUM_WHITE", "PER_WHITE"
)
outcome_labels <- c(
  NUM_ECDIS = "Econ. Disadvantaged", PER_ECDIS = "Econ. Disadvantaged",
  NUM_FREE_LUNCH = "Free Lunch", PER_FREE_LUNCH = "Free Lunch",
  NUM_SWD   = "Students w/ Disabilities", PER_SWD = "Students w/ Disabilities",
  NUM_BLACK = "Black", PER_BLACK = "Black",
  NUM_HISP  = "Hispanic", PER_HISP = "Hispanic",
  NUM_WHITE = "White", PER_WHITE = "White"
)

all_results <- list()
for (sample_key in names(samples)) {
  cfg <- samples[[sample_key]]
  dt <- fread(cfg$path)
  message("\n=== CS Demographic: ", cfg$label, " ===")

  for (type_name in names(treat_vars_cs)) {
    tv <- treat_vars_cs[[type_name]]
    if (!(tv %in% names(dt))) { message("  Skipping ", type_name, ": ", tv, " not found."); next }

    for (yname in outcome_vars) {
      if (!(yname %in% names(dt))) { message("  Skipping outcome ", yname, ": not found."); next }
      cat("  ", type_name, "/", yname, "\n")
      res <- run_cs_generic(dt, yname, tv)
      all_results[[length(all_results) + 1]] <- data.table(
        sample = cfg$label, type = type_name, outcome = yname,
        outcome_label = outcome_labels[yname], att = res$att, se = res$se,
        n_obs = res$n_obs, n_units = res$n_units, status = res$status
      )
    }
  }
}

results_dt <- rbindlist(all_results)
results_dt[, est_star := add_stars(att, se)]
results_dt[, se_paren := sprintf("(%.3f)", se)]
results_dt[, is_pct := startsWith(outcome, "PER_")]

fwrite(results_dt, file.path(base_path, "CS_demographic_results.csv"))
message("Wrote: CS_demographic_results.csv")

# Build one LaTeX table per sample x (log/pct) -- 4 tables total
for (sample_key in names(samples)) {
  cfg <- samples[[sample_key]]
  for (kind in c(TRUE, FALSE)) {
    sub <- results_dt[sample == cfg$label & is_pct == kind]
    if (nrow(sub) == 0) next

    wide_est <- dcast(sub, outcome_label ~ type, value.var = "est_star")
    wide_se  <- dcast(sub, outcome_label ~ type, value.var = "se_paren")

    rows <- character(0)
    for (i in seq_len(nrow(wide_est))) {
      rows <- c(rows,
                sprintf("%s & %s & %s & %s \\\\", wide_est$outcome_label[i], wide_est$ALL[i], wide_est$FAMILY[i], wide_est$ELDERLY[i]),
                sprintf(" & %s & %s & %s \\\\", wide_se$ALL[i], wide_se$FAMILY[i], wide_se$ELDERLY[i])
      )
    }

    tex <- c(
      "\\begin{table}[htbp]", "\\centering",
      sprintf("\\caption{CS Estimates, %s (%s)}", if (kind) "Enrollment Shares" else "Log Enrollment Counts", cfg$label),
      "\\begin{tabular}{lccc}", "\\toprule", "Outcome & All & Family & Non-Family \\\\", "\\midrule",
      rows, "\\bottomrule", "\\end{tabular}",
      "\\vspace{0.6em}\\parbox{0.9\\linewidth}{\\footnotesize",
      "\\emph{**} $p<0.05$, \\emph{***} $p<0.01$.}", "\\end{table}"
    )
    fname <- sprintf("CS_demographic_%s_%s.tex", sample_key, if (kind) "pct" else "log")
    writeLines(tex, file.path(base_path, fname))
    message("Wrote: ", fname)
  }
}

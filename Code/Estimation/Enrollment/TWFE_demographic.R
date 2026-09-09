rm(list = ls())

library(data.table)
library(fixest)
library(here)

base_path <- here("data", "processed")

idname <- "BEDSCODE"
tname  <- "YEAR"

treat_vars <- list(
  ALL     = "LIHTC_OPEN_PIS_ALL",
  FAMILY  = "LIHTC_OPEN_PIS_FAMILY",
  ELDERLY = "LIHTC_OPEN_PIS_ELDERLY"
)

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

samples <- list(
  full   = list(path = file.path(base_path, "demographic_panel_full.csv"),   label = "Full Sample"),
  no_nyc = list(path = file.path(base_path, "demographic_panel_no_nyc.csv"), label = "No-NYC Sample")
)

stars_fn  <- function(p) ifelse(p < 0.01, "***", ifelse(p < 0.05, "**", ""))
add_stars <- function(est, se) { p <- 2 * pnorm(-abs(est / se)); paste0(sprintf("%.3f", est), stars_fn(p)) }
get_est <- function(m, tv) tryCatch(coef(m)[[tv]], error = function(e) NA_real_)
get_se  <- function(m, tv) tryCatch(se(m)[[tv]],   error = function(e) NA_real_)

all_results <- list()

for (sample_key in names(samples)) {
  cfg <- samples[[sample_key]]
  dt <- fread(cfg$path)
  message("\n=== TWFE Demographic: ", cfg$label, " ===")

  dt[, (idname) := as.character(get(idname))]
  dt[, (tname) := as.integer(get(tname))]

  for (yname in outcome_vars) {
    if (!(yname %in% names(dt))) { message("  Skipping ", yname, ": not found."); next }

    dt[, (yname) := as.numeric(get(yname))]
    if (startsWith(yname, "NUM_")) {
      depvar <- paste0("log_", yname)
      dt[, (depvar) := fifelse(get(yname) > 0, log(get(yname)), NA_real_)]
    } else {
      depvar <- yname
    }

    for (type_name in names(treat_vars)) {
      tv <- treat_vars[[type_name]]
      if (!(tv %in% names(dt))) next

      dt[, (tv) := as.numeric(get(tv))]
      dt[is.na(get(tv)), (tv) := 0]

      m <- tryCatch(
        feols(as.formula(paste0(depvar, " ~ ", tv, " | ", idname, " + ", tname)),
              data = dt, cluster = as.formula(paste0("~", idname))),
        error = function(e) NULL
      )
      if (is.null(m)) next

      all_results[[length(all_results) + 1]] <- data.table(
        sample = cfg$label, outcome = yname, outcome_label = outcome_labels[yname], type = type_name,
        est = get_est(m, tv), se = get_se(m, tv), n_obs = nobs(m)
      )
    }
  }
}

results_dt <- rbindlist(all_results)
results_dt[, est_star := add_stars(est, se)]
results_dt[, se_paren := sprintf("(%.3f)", se)]
results_dt[, is_pct := startsWith(outcome, "PER_")]

fwrite(results_dt, file.path(base_path, "TWFE_demographic_results.csv"))

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
      sprintf("\\caption{TWFE Estimates, %s (%s)}", if (kind) "Enrollment Shares" else "Log Enrollment Counts", cfg$label),
      "\\begin{tabular}{lccc}", "\\toprule", "Outcome & All & Family & Non-Family \\\\", "\\midrule",
      rows, "\\midrule", "School FE & \\multicolumn{3}{c}{Yes} \\\\", "Year FE & \\multicolumn{3}{c}{Yes} \\\\",
      "\\bottomrule", "\\end{tabular}",
      "\\vspace{0.6em}\\parbox{0.9\\linewidth}{\\footnotesize",
      "\\emph{**} $p<0.05$, \\emph{***} $p<0.01$.}", "\\end{table}"
    )
    fname <- sprintf("TWFE_demographic_%s_%s.tex", sample_key, if (kind) "pct" else "log")
    writeLines(tex, file.path(base_path, fname))
    message("Wrote: ", fname)
  }
}

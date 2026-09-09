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

samples <- list(
  full   = list(path = file.path(base_path, "headcount_panel_full.csv"),   label = "Full Sample"),
  no_nyc = list(path = file.path(base_path, "headcount_panel_no_nyc.csv"), label = "No-NYC Sample")
)

stars_fn  <- function(p) ifelse(p < 0.01, "***", ifelse(p < 0.05, "**", ""))
add_stars <- function(est, se) { p <- 2 * pnorm(-abs(est / se)); paste0(sprintf("%.3f", est), stars_fn(p)) }
get_est <- function(m, tv) tryCatch(coef(m)[[tv]], error = function(e) NA_real_)
get_se  <- function(m, tv) tryCatch(se(m)[[tv]],   error = function(e) NA_real_)

for (sample_key in names(samples)) {
  cfg <- samples[[sample_key]]
  dt <- fread(cfg$path)
  message("\n=== TWFE Headcount: ", cfg$label, " ===")

  dt[, (idname) := as.character(get(idname))]
  dt[, (tname) := as.integer(get(tname))]
  dt[, K12_TOTAL := as.numeric(K12_TOTAL)]
  dt[, log_K12_TOTAL := fifelse(K12_TOTAL > 0, log(K12_TOTAL), NA_real_)]

  results <- list()
  for (type_name in names(treat_vars)) {
    tv <- treat_vars[[type_name]]
    if (!(tv %in% names(dt))) { message("  Skipping ", type_name, ": ", tv, " not found."); next }

    dt[, (tv) := as.numeric(get(tv))]
    dt[is.na(get(tv)), (tv) := 0]

    m <- feols(as.formula(paste0("log_K12_TOTAL ~ ", tv, " | ", idname, " + ", tname)),
               data = dt, cluster = as.formula(paste0("~", idname)))
    results[[type_name]] <- m
  }

  n_obs <- tryCatch(nobs(if (!is.null(results$FAMILY)) results$FAMILY else results$ALL), error = function(e) NA_integer_)

  tex <- c(
    "\\begin{table}[htbp]", "\\centering",
    sprintf("\\caption{TWFE Estimates of LIHTC Impacts on Log K-12 Enrollment (%s)}", cfg$label),
    "\\begin{tabular}{lcccc}", "\\toprule",
    "Outcome & All & Family & Non-Family & $N$ \\\\", "\\midrule",
    sprintf("Log K-12 Enrollment & %s & %s & %s & %s \\\\",
            if (!is.null(results$ALL)) add_stars(get_est(results$ALL, treat_vars$ALL), get_se(results$ALL, treat_vars$ALL)) else "---",
            if (!is.null(results$FAMILY)) add_stars(get_est(results$FAMILY, treat_vars$FAMILY), get_se(results$FAMILY, treat_vars$FAMILY)) else "---",
            if (!is.null(results$ELDERLY)) add_stars(get_est(results$ELDERLY, treat_vars$ELDERLY), get_se(results$ELDERLY, treat_vars$ELDERLY)) else "---",
            format(n_obs, big.mark = ",")),
    sprintf(" & (%.3f) & (%.3f) & (%.3f) & \\\\",
            if (!is.null(results$ALL)) get_se(results$ALL, treat_vars$ALL) else NA,
            if (!is.null(results$FAMILY)) get_se(results$FAMILY, treat_vars$FAMILY) else NA,
            if (!is.null(results$ELDERLY)) get_se(results$ELDERLY, treat_vars$ELDERLY) else NA),
    "\\midrule", "School FE & \\multicolumn{4}{c}{Yes} \\\\", "Year FE & \\multicolumn{4}{c}{Yes} \\\\",
    "\\bottomrule", "\\end{tabular}",
    "\\vspace{0.6em}\\parbox{0.9\\linewidth}{\\footnotesize",
    "\\emph{Notes:} Standard errors clustered at the school level in parentheses.",
    "\\emph{**} $p<0.05$, \\emph{***} $p<0.01$.}", "\\end{table}"
  )

  out_file <- file.path(base_path, sprintf("TWFE_headcount_%s.tex", sample_key))
  writeLines(tex, out_file)
  message("Wrote: ", out_file)
}

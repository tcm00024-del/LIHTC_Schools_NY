rm(list = ls())
library(here)
source(here("code", "estimation", "enrollment", "enrollment_estimation_helpers.R"))

Sys.setenv(OMP_NUM_THREADS = 1); Sys.setenv(MKL_NUM_THREADS = 1); options(mc.cores = 1)

base_path <- here("data", "processed")
samples <- samples_enrollment(base_path, "headcount_panel")

all_results <- list()
for (sample_key in names(samples)) {
  cfg <- samples[[sample_key]]
  dt <- fread(cfg$path)
  message("\n=== CS Headcount: ", cfg$label, " ===")

  for (type_name in names(treat_vars_cs)) {
    tv <- treat_vars_cs[[type_name]]
    if (!(tv %in% names(dt))) { message("  Skipping ", type_name, ": ", tv, " not found."); next }

    cat("  Running", type_name, "\n")
    res <- run_cs_generic(dt, "K12_TOTAL", tv, transform = "log")
    all_results[[length(all_results) + 1]] <- data.table(
      sample = cfg$label, type = type_name, att = res$att, se = res$se,
      n_obs = res$n_obs, n_units = res$n_units, status = res$status
    )
  }
}

results_dt <- rbindlist(all_results)
results_dt[, est_star := add_stars(att, se)]
print(results_dt)
fwrite(results_dt, file.path(base_path, "CS_headcount_results.csv"))
message("Wrote: CS_headcount_results.csv")

# -----------------------------
# LaTeX TABLE (one per sample)
# -----------------------------
for (sample_key in names(samples)) {
  cfg <- samples[[sample_key]]
  sub <- results_dt[sample == cfg$label]
  if (nrow(sub) == 0) next

  get_row <- function(type_name) sub[type == type_name]

  fmt_est <- function(row) if (nrow(row) == 0 || is.na(row$att)) "---" else add_stars(row$att, row$se)
  fmt_se  <- function(row) if (nrow(row) == 0 || is.na(row$se)) "" else sprintf("(%.3f)", row$se)
  fmt_n   <- function(row) if (nrow(row) == 0 || is.na(row$n_obs)) "" else format(row$n_obs, big.mark = ",")

  r_all <- get_row("ALL"); r_fam <- get_row("FAMILY"); r_eld <- get_row("ELDERLY")
  n_str <- if (nrow(r_fam) > 0) fmt_n(r_fam) else fmt_n(r_all)

  tex <- c(
    "\\begin{table}[htbp]", "\\centering",
    sprintf("\\caption{Callaway--Sant'Anna Estimates of LIHTC Impacts on Log K-12 Enrollment (%s)}", cfg$label),
    sprintf("\\label{tab:cs_headcount_%s}", sample_key),
    "\\begin{tabular}{lcccc}", "\\toprule",
    "Outcome & All & Family & Non-Family & $N$ \\\\", "\\midrule",
    sprintf("Log K-12 Enrollment & %s & %s & %s & %s \\\\",
            fmt_est(r_all), fmt_est(r_fam), fmt_est(r_eld), n_str),
    sprintf(" & %s & %s & %s & \\\\", fmt_se(r_all), fmt_se(r_fam), fmt_se(r_eld)),
    "\\bottomrule", "\\end{tabular}",
    "\\vspace{0.6em}\\parbox{0.9\\linewidth}{\\footnotesize",
    "\\emph{Notes:} Callaway--Sant'Anna estimate (\\texttt{aggte(type = ``simple'')}) of the effect of a nearby LIHTC opening on log K-12 enrollment.",
    "Control group is \\emph{never-treated} units. Standard errors, clustered at the school level, are shown in parentheses.",
    "\\emph{**} $p<0.05$, \\emph{***} $p<0.01$.}", "\\end{table}"
  )

  out_file <- file.path(base_path, sprintf("CS_headcount_%s.tex", sample_key))
  writeLines(tex, out_file)
  message("Wrote: ", out_file)
}

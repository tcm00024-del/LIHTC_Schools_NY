# ============================================================
# TWFE_estimates.R
# ============================================================
rm(list = ls())

library(data.table)
library(fixest)
library(here)

options("modelsummary_format_numeric_latex" = "plain")

base_path <- here("data", "processed", "Subgroups")

yname  <- "mean_scale_score_weighted"
tname  <- "YEAR"
idname <- "BEDSCODE"

treat_vars <- list(
  ALL     = "LIHTC_OPEN_PIS_ALL",
  FAMILY  = "LIHTC_OPEN_PIS_FAMILY",
  ELDERLY = "LIHTC_OPEN_PIS_ELDERLY"
)

pretty_subgroup <- function(x) {
  if (grepl("All_Students",                   x)) return("All")
  if (grepl("General_Education",              x)) return("General Education")
  if (grepl("Not_Economically_Disadvantaged", x)) return("Not Economically Disadvantaged")
  if (grepl("Economically_Disadvantaged",     x)) return("Economically Disadvantaged")
  if (grepl("Female",                         x)) return("Female")
  if (grepl("Male",                         x)) return("Male")
  if (grepl("Black",                          x)) return("Black")
  if (grepl("Hispanic",                         x)) return("Hispanic")
  if (grepl("White",                         x)) return("White")
  return(x)
}

desired_order <- c(
  "All", "General Education",
  "Economically Disadvantaged", "Not Economically Disadvantaged",
  "Female", "Male", "Black", "Hispanic", "White"
)

samples <- list(
  full   = list(path = file.path(base_path, "Full"),   pattern = "^full_subgroup_.*_collapsed\\.csv$",    label = "Full Sample"),
  no_nyc = list(path = file.path(base_path, "No_NYC"), pattern = "^no_nyc_subgroup_.*_collapsed\\.csv$",  label = "No-NYC Sample")
)

# -----------------------------
# Formatting helpers (5% cap -- no 10% stars)
# -----------------------------
stars_fn  <- function(p) ifelse(p < 0.01, "***", ifelse(p < 0.05, "**", ""))
add_stars <- function(est, se) { p <- 2 * pnorm(-abs(est / se)); paste0(sprintf("%.3f", est), stars_fn(p)) }
get_est   <- function(m, tv) tryCatch(coef(m)[[tv]], error = function(e) NA_real_)
get_se    <- function(m, tv) tryCatch(se(m)[[tv]],   error = function(e) NA_real_)

for (sample_key in names(samples)) {
  cfg <- samples[[sample_key]]
  message("\n=== Running TWFE: ", cfg$label, " ===")

  files <- list.files(cfg$path, pattern = cfg$pattern, full.names = TRUE)
  if (length(files) == 0) stop("No files matched pattern in: ", cfg$path)
  message("Found ", length(files), " subgroup files.")

  results <- list()

  for (file_path in files) {
    subgroup_label <- pretty_subgroup(basename(file_path))
    if (!(subgroup_label %in% desired_order)) next

    cat("Processing:", subgroup_label, "\n")

    tryCatch({
      dt <- fread(file_path)
      dt[, (idname) := as.character(get(idname))]
      dt[, (tname)  := as.integer(get(tname))]

      type_results <- list()

      for (type_name in names(treat_vars)) {
        tv <- treat_vars[[type_name]]

        if (!(tv %in% names(dt))) {
          message("  Skipping ", type_name, ": column ", tv, " not found.")
          next
        }
        if (!(yname %in% names(dt))) {
          message("  Skipping ", type_name, ": outcome column not found.")
          next
        }

        dt[, (tv) := as.numeric(get(tv))]
        dt[is.na(get(tv)), (tv) := 0]

        m <- feols(
          as.formula(paste0(yname, " ~ ", tv, " | ", idname, " + ", tname)),
          data    = dt,
          cluster = as.formula(paste0("~", idname))
        )

        type_results[[type_name]] <- m
      }

      results[[subgroup_label]] <- type_results

    }, error = function(e) {
      message("FAILED: ", basename(file_path), " --> ", e$message)
    })
  }

  # -----------------------------
  # BUILD LaTeX TABLE (this sample)
  # -----------------------------
  twfe_lines <- c(
    "\\begin{table}[htbp]", "\\centering",
    sprintf("\\caption{TWFE Estimates of LIHTC Impacts on Test Scores, by Subgroup and Development Type (%s)}", cfg$label),
    sprintf("\\label{tab:twfe_%s}", sample_key),
    "\\begin{tabular}{lcccc}", "\\toprule",
    "Subgroup & All & Family & Non-Family & $N$ \\\\", "\\midrule"
  )

  for (sg in desired_order) {
    if (!(sg %in% names(results))) next
    r <- results[[sg]]

    est_all <- if (!is.null(r$ALL))     get_est(r$ALL,     treat_vars$ALL)     else NA_real_
    se_all  <- if (!is.null(r$ALL))     get_se(r$ALL,      treat_vars$ALL)     else NA_real_
    est_fam <- if (!is.null(r$FAMILY))  get_est(r$FAMILY,  treat_vars$FAMILY)  else NA_real_
    se_fam  <- if (!is.null(r$FAMILY))  get_se(r$FAMILY,   treat_vars$FAMILY)  else NA_real_
    est_eld <- if (!is.null(r$ELDERLY)) get_est(r$ELDERLY, treat_vars$ELDERLY) else NA_real_
    se_eld  <- if (!is.null(r$ELDERLY)) get_se(r$ELDERLY,  treat_vars$ELDERLY) else NA_real_

    n_obs <- tryCatch(nobs(if (!is.null(r$FAMILY)) r$FAMILY else r$ALL), error = function(e) NA_integer_)

    fmt_est <- function(est, se) if (is.na(est) || is.na(se)) "---" else add_stars(est, se)
    fmt_se  <- function(se) if (is.na(se)) "" else sprintf("(%.3f)", se)
    fmt_n   <- function(n) if (is.na(n)) "" else format(n, big.mark = ",")

    twfe_lines <- c(
      twfe_lines,
      sprintf("%-30s & %s & %s & %s & %s \\\\", sg, fmt_est(est_all, se_all), fmt_est(est_fam, se_fam), fmt_est(est_eld, se_eld), fmt_n(n_obs)),
      sprintf("%-30s & %s & %s & %s & \\\\", "", fmt_se(se_all), fmt_se(se_fam), fmt_se(se_eld))
    )
  }

  twfe_lines <- c(
    twfe_lines,
    "\\midrule",
    "School FE & \\multicolumn{4}{c}{Yes} \\\\",
    "Year FE   & \\multicolumn{4}{c}{Yes} \\\\",
    "\\bottomrule", "\\end{tabular}",
    "\\vspace{0.6em}", "\\parbox{0.9\\linewidth}{\\footnotesize",
    "\\emph{Notes:} Each cell reports a separate TWFE estimate of the effect of a nearby LIHTC opening on mean weighted test scores.",
    "Treatment indicators equal one in all years on or after the first opening of the indicated type.",
    "Development population targets are identified using the Atkins--O'Regan (2014) bedroom-distribution algorithm.",
    "All specifications include school and year fixed effects.",
    "Standard errors, clustered at the school level, are shown in parentheses.",
    "TWFE estimates are reported for completeness; Callaway--Sant'Anna and Sun--Abraham estimates are preferred given staggered treatment timing.",
    "\\emph{**} $p<0.05$, \\emph{***} $p<0.01$.",
    "}", "\\end{table}"
  )

  out_tex <- file.path(base_path, sprintf("TWFE_table_%s.tex", sample_key))
  writeLines(twfe_lines, out_tex)
  message("Wrote: ", out_tex)
}

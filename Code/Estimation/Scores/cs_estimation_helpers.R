# ============================================================
# cs_estimation_helpers.R
# Shared functions for Callaway–Sant'Anna estimation across
# standard / bootstrapped-SE / left-censoring-excluded variants
# ============================================================

library(data.table)
library(did)

# -----------------------------
# Pretty names / ordering (unchanged from original -- substring
# matching still works against the new safe-filename convention)
# -----------------------------
pretty_subgroup <- function(x) {
  if (grepl("All_Students",                   x)) return("All")
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

desired_order <- c(
  "All", "General Education",
  "Economically Disadvantaged", "Not Economically Disadvantaged",
  "Female", "Male", "Black", "Hispanic", "White"
)

# PIS-based treatment columns from the new pipeline (no rank suffix --
# we only kept k=1 nearest-school matching)
treat_vars <- list(
  ALL     = "LIHTC_OPEN_PIS_ALL",
  FAMILY  = "LIHTC_OPEN_PIS_FAMILY",
  ELDERLY = "LIHTC_OPEN_PIS_ELDERLY"
)

yname  <- "mean_scale_score_weighted"
tname  <- "YEAR"
idname <- "BEDSCODE"

# -----------------------------
# CS runner: one treatment column on one subgroup's data
# drop_left_censored: excludes units first-treated in panel_start_year
#   (2013), since they have zero observable pre-treatment periods under
#   this treatment definition. Never-treated units (G_CS == 0) are
#   always retained regardless of this flag.
# -----------------------------
run_cs_one <- function(dt, tv, yname, tname, idname, control_group,
                       bstrap = FALSE, biters = 1000,
                       drop_left_censored = FALSE, panel_start_year = 2013) {
  tryCatch({
    dt <- copy(dt)

    dt[, (tv)    := as.numeric(get(tv))]
    dt[is.na(get(tv)), (tv) := 0]
    dt[, ID_NUM  := as.integer(factor(get(idname)))]
    dt[, (tname) := as.integer(get(tname))]

    dt[, G_CS := {
      tv_vec <- get(tv)
      t_vec  <- get(tname)
      if (any(tv_vec > 0)) min(t_vec[tv_vec > 0]) else 0L
    }, by = .(ID_NUM)]
    dt[is.na(G_CS), G_CS := 0L]
    dt[, G_CS := as.double(G_CS)]   # <-- force double so did's internal Inf-recoding of
    #     never-treated units doesn't get truncated to NA

    n_left_censored <- NA_integer_
    if (drop_left_censored) {
      left_censored_ids <- unique(dt[G_CS == panel_start_year, ID_NUM])
      n_left_censored <- length(left_censored_ids)
      if (n_left_censored > 0) {
        message("    Dropping ", n_left_censored, " left-censored unit(s) (first treated ", panel_start_year, ")")
        dt <- dt[!(ID_NUM %in% left_censored_ids)]
      }
    }

    message("  NA check — outcome: ", sum(is.na(dt[[yname]])),
            " | treatment: ", sum(is.na(dt[[tv]])),
            " | year: ", sum(is.na(dt[[tname]])),
            " | id: ", sum(is.na(dt[["ID_NUM"]])))

    cs_att <- att_gt(
      yname                  = yname,
      tname                  = tname,
      idname                 = "ID_NUM",
      gname                  = "G_CS",
      data                   = dt,
      control_group          = control_group,
      clustervars            = "ID_NUM",
      panel                  = TRUE,
      #allow_unbalanced_panel = TRUE,
      bstrap                 = bstrap,
      biters                 = biters
    )

    cs_dyn <- aggte(cs_att, type = "dynamic")
    n_obs <- tryCatch(as.integer(nrow(cs_att$DIDparams$data)), error = function(e) NA_integer_)

    list(
      att                      = as.numeric(cs_dyn$overall.att),
      se                       = as.numeric(cs_dyn$overall.se),
      n_obs                    = n_obs,
      n_units                  = as.integer(sum(cs_att$n)),
      n_left_censored_dropped  = n_left_censored
    )

  }, error = function(e) {
    message("    CS failed for ", tv, ": ", e$message)
    NULL
  })
}

# -----------------------------
# Runs CS across every subgroup file x treatment type for one sample
# (one call of this = the middle+inner loop; the outer sample loop
# lives in each driver script)
# -----------------------------
run_cs_for_sample <- function(data_path, file_pattern, yname, tname, idname,
                              treat_vars, control_group, desired_order,
                              bstrap = FALSE, biters = 1000,
                              drop_left_censored = FALSE, panel_start_year = 2013) {
  files <- list.files(data_path, pattern = file_pattern, full.names = TRUE)
  if (length(files) == 0) stop("No files matched pattern in: ", data_path)
  message("Found ", length(files), " subgroup files in ", data_path)

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
        cat("  Running CS for", type_name, "\n")
        type_results[[type_name]] <- run_cs_one(
          dt, tv, yname, tname, idname, control_group,
          bstrap = bstrap, biters = biters,
          drop_left_censored = drop_left_censored, panel_start_year = panel_start_year
        )
      }
      results[[subgroup_label]] <- type_results

    }, error = function(e) {
      message("FAILED to load: ", basename(file_path), " --> ", e$message)
    })
  }

  results
}

# -----------------------------
# Long-format results table
# -----------------------------
build_cs_long_table <- function(results) {
  rows <- lapply(names(results), function(sg) {
    r <- results[[sg]]
    lapply(names(r), function(tp) {
      res <- r[[tp]]
      if (is.null(res)) return(NULL)
      data.table(
        Subgroup                = sg,
        Type                    = tp,
        ATT                     = res$att,
        SE                      = res$se,
        N_obs                   = res$n_obs,
        N_units                 = res$n_units,
        N_left_censored_dropped = res$n_left_censored_dropped
      )
    })
  })
  rbindlist(unlist(rows, recursive = FALSE), fill = TRUE)
}

# -----------------------------
# LaTeX table builder
# -----------------------------
build_cs_tex_table <- function(results, desired_order, control_group, caption, label) {
  stars_fn  <- function(p) ifelse(p < 0.01, "***", ifelse(p < 0.05, "**", ""))
  add_stars <- function(est, se) { p <- 2 * pnorm(-abs(est / se)); paste0(sprintf("%.3f", est), stars_fn(p)) }
  fmt_est   <- function(res) if (is.null(res) || is.na(res$att) || is.na(res$se)) "---" else add_stars(res$att, res$se)
  fmt_se    <- function(res) if (is.null(res) || is.na(res$se)) "" else sprintf("(%.3f)", res$se)
  fmt_n     <- function(res) if (is.null(res) || is.na(res$n_obs)) "" else format(res$n_obs, big.mark = ",")

  tex_lines <- c(
    "\\begin{table}[htbp]", "\\centering",
    sprintf("\\caption{%s}", caption),
    sprintf("\\label{%s}", label),
    "\\begin{tabular}{lcccc}", "\\toprule",
    "Subgroup & All & Family & Non-Family & $N$ \\\\", "\\midrule"
  )

  for (sg in desired_order) {
    if (!(sg %in% names(results))) next
    r <- results[[sg]]
    n_str <- fmt_n(if (!is.null(r$FAMILY)) r$FAMILY else r$ALL)
    tex_lines <- c(
      tex_lines,
      sprintf("%-30s & %s & %s & %s & %s \\\\", sg, fmt_est(r$ALL), fmt_est(r$FAMILY), fmt_est(r$ELDERLY), n_str),
      sprintf("%-30s & %s & %s & %s & \\\\", "", fmt_se(r$ALL), fmt_se(r$FAMILY), fmt_se(r$ELDERLY))
    )
  }

  c(tex_lines, "\\bottomrule", "\\end{tabular}",
    "\\vspace{0.6em}", "\\parbox{0.9\\linewidth}{\\footnotesize",
    "\\emph{Notes:} Each cell reports a Callaway--Sant'Anna estimate of the effect of a nearby LIHTC opening on mean weighted test scores.",
    "Estimates are overall ATT aggregates from dynamic aggregations (\\texttt{aggte(type = ``dynamic'')}).",
    "Treatment indicators equal one in all years on or after the first opening of the indicated type.",
    "Development population targets are identified using the Atkins--O'Regan (2014) bedroom-distribution algorithm.",
    sprintf("The control group is \\emph{%s} units.", gsub("nevertreated", "never-treated", gsub("notyettreated", "not-yet-treated", control_group))),
    "Standard errors, clustered at the school level, are shown in parentheses.",
    "\\emph{**} $p<0.05$, \\emph{***} $p<0.01$.",
    "}", "\\end{table}"
  )
}

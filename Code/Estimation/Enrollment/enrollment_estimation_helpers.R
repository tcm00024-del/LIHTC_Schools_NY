# ============================================================
# enrollment_estimation_helpers.R
# Shared functions for CS and Sun-Abraham on enrollment outcomes
# ============================================================

library(data.table)
library(did)
library(fixest)

idname_raw <- "BEDSCODE"
tname      <- "YEAR"

treat_vars_cs <- list(
  ALL     = "LIHTC_OPEN_PIS_ALL",
  FAMILY  = "LIHTC_OPEN_PIS_FAMILY",
  ELDERLY = "LIHTC_OPEN_PIS_ELDERLY"
)
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
type_labels <- c(ALL = "All", FAMILY = "Family", ELDERLY = "Non-Family")
type_colors <- c(ALL = "#2166ac", FAMILY = "#d6604d", ELDERLY = "#4dac26")

samples_enrollment <- function(base_path, file_prefix) {
  list(
    full   = list(path = file.path(base_path, sprintf("%s_full.csv",   file_prefix)), label = "Full Sample"),
    no_nyc = list(path = file.path(base_path, sprintf("%s_no_nyc.csv", file_prefix)), label = "No-NYC Sample")
  )
}

# -----------------------------
# Transform: NUM_ -> log, PER_ -> level; explicit override for
# outcomes like K12_TOTAL that don't follow the NUM_/PER_ convention
# -----------------------------
apply_transform <- function(d, yname, transform = NULL) {
  d[, (yname) := as.numeric(get(yname))]

  if (is.null(transform)) {
    if (startsWith(yname, "NUM_")) transform <- "log"
    else if (startsWith(yname, "PER_")) transform <- "level"
    else stop("Cannot infer transform for '", yname, "' -- pass transform explicitly.")
  }

  if (transform == "log") {
    depvar <- paste0("log_", yname)
    d[, (depvar) := fifelse(get(yname) > 0, log(get(yname)), NA_real_)]
  } else {
    depvar <- yname
  }

  list(data = d, depvar = depvar, transform = transform)
}

# -----------------------------
# CS runner (standard only -- no bootstrap/left-censoring needed here)
# Same G_CS-double-cast fix and allow_unbalanced_panel as the test-score
# helper file; this is the bug we found and fixed together earlier today.
# -----------------------------
run_cs_generic <- function(dt, yname, treat_var, transform = NULL,
                           control_group = "nevertreated") {
  tryCatch({
    d <- copy(dt)
    d[, (idname_raw) := as.character(get(idname_raw))]
    d[, ID_NUM := as.integer(factor(get(idname_raw)))]
    d[, (tname) := as.integer(get(tname))]
    d[, (treat_var) := as.numeric(get(treat_var))]
    d[is.na(get(treat_var)), (treat_var) := 0]

    dep_info <- apply_transform(d, yname, transform)
    d <- dep_info$data
    depvar <- dep_info$depvar

    d_reg <- d[!is.na(get(depvar)) & !is.na(get(treat_var)) & !is.na(ID_NUM) & !is.na(get(tname))]
    if (nrow(d_reg) == 0) return(list(att = NA_real_, se = NA_real_, n_obs = 0L, n_units = 0L, status = "No usable observations"))

    d_reg[, G_CS := {
      tv_vec <- get(treat_var); t_vec <- get(tname)
      if (any(tv_vec > 0)) min(t_vec[tv_vec > 0]) else 0L
    }, by = .(ID_NUM)]
    d_reg[is.na(G_CS), G_CS := 0L]
    d_reg[, G_CS := as.double(G_CS)]   # critical: must stay double, not integer

    n_treated <- uniqueN(d_reg[G_CS > 0, ID_NUM])
    n_never   <- uniqueN(d_reg[G_CS == 0, ID_NUM])
    if (n_treated == 0 || n_never == 0) {
      return(list(att = NA_real_, se = NA_real_, n_obs = nrow(d_reg), n_units = uniqueN(d_reg$ID_NUM),
                  status = "Insufficient treated/control support"))
    }

    cs_att <- att_gt(
      yname = depvar, tname = tname, idname = "ID_NUM", gname = "G_CS",
      data = d_reg, control_group = control_group, clustervars = "ID_NUM",
      panel = TRUE, allow_unbalanced_panel = TRUE, bstrap = FALSE
    )
    cs_simple <- aggte(cs_att, type = "simple")

    list(
      att = as.numeric(cs_simple$overall.att), se = as.numeric(cs_simple$overall.se),
      n_obs = nrow(cs_att$DIDparams$data), n_units = uniqueN(cs_att$DIDparams$data$ID_NUM),
      status = "OK"
    )
  }, error = function(e) {
    message("    CS failed for ", yname, " / ", treat_var, ": ", e$message)
    list(att = NA_real_, se = NA_real_, n_obs = NA_integer_, n_units = NA_integer_, status = paste("ERROR:", e$message))
  })
}

# -----------------------------
# SA runner (mirrors the test-score version: reconstructs cohort from
# EVENT_TIME_PIS_{TYPE} / NEVER_TREATED_PIS_{TYPE} rather than a raw
# FIRST_PIS column, since that's what Step 5's exposure table carries)
# -----------------------------
run_sa_generic <- function(dt, yname, event_time_var, never_treated_var,
                           transform = NULL, event_window = c(-3, 3), ref_period = -1) {
  tryCatch({
    d <- copy(dt)
    d[, (idname_raw) := as.character(get(idname_raw))]
    d[, (tname) := as.integer(get(tname))]

    dep_info <- apply_transform(d, yname, transform)
    d <- dep_info$data
    depvar <- dep_info$depvar

    d[, first_treat_year := as.numeric(get(tname)) - as.numeric(get(event_time_var))]
    d[, cohort := fifelse(get(never_treated_var) == 1, Inf, first_treat_year)]

    fml <- as.formula(sprintf("%s ~ sunab(cohort, %s, ref.p = %d) | %s + %s",
                              depvar, tname, ref_period, idname_raw, tname))
    m <- feols(fml, data = d, cluster = as.formula(paste0("~", idname_raw)))

    coefs <- as.data.table(coeftable(m, agg = "period"), keep.rownames = "term")
    setnames(coefs, c("term", "est", "se", "t", "p"))
    coefs[, event_time := as.integer(sub(".*::([-0-9]+)$", "\\1", term))]
    coefs <- coefs[!is.na(event_time) & event_time >= event_window[1] & event_time <= event_window[2]]

    ref_row <- data.table(term = paste0("ref::", ref_period), est = 0, se = 0, t = NA_real_, p = NA_real_, event_time = ref_period)
    coefs <- rbindlist(list(coefs, ref_row), fill = TRUE)
    setorder(coefs, event_time)
    coefs[, `:=`(ci_lo = est - 1.96 * se, ci_hi = est + 1.96 * se)]
    coefs
  }, error = function(e) {
    message("    SA failed for ", yname, " / ", event_time_var, ": ", e$message)
    NULL
  })
}

# -----------------------------
# Stars capped at 5%/1% -- no 10% threshold
# -----------------------------
stars_fn  <- function(p) ifelse(p < 0.01, "***", ifelse(p < 0.05, "**", ""))
add_stars <- function(est, se) { p <- 2 * pnorm(-abs(est / se)); paste0(sprintf("%.3f", est), stars_fn(p)) }

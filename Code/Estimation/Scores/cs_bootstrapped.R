# ============================================================
# CS_bootstrap.R
# ============================================================
rm(list = ls())

Sys.setenv(OMP_NUM_THREADS = 1)
Sys.setenv(MKL_NUM_THREADS = 1)
options(mc.cores = 1)

library(here)
source(here("code", "estimation", "scores", "cs_estimation_helpers.R"))

base_path <- here("data", "processed", "Subgroups")
control_group <- "nevertreated"
biters <- 1000

samples <- list(
  full   = list(path = file.path(base_path, "Full"),   pattern = "^full_subgroup_.*_collapsed\\.csv$",    label = "Full Sample"),
  no_nyc = list(path = file.path(base_path, "No_NYC"), pattern = "^no_nyc_subgroup_.*_collapsed\\.csv$",  label = "No-NYC Sample")
)

for (sample_key in names(samples)) {
  cfg <- samples[[sample_key]]
  message("\n=== Running CS (bootstrap, ", biters, " iters): ", cfg$label, " ===")

  results <- run_cs_for_sample(
    data_path = cfg$path, file_pattern = cfg$pattern,
    yname = yname, tname = tname, idname = idname,
    treat_vars = treat_vars, control_group = control_group,
    desired_order = desired_order,
    bstrap = TRUE, biters = biters, drop_left_censored = FALSE
  )

  fwrite(build_cs_long_table(results), file.path(base_path, sprintf("CS_ATT_%s_bootstrap.csv", sample_key)))

  tex_lines <- build_cs_tex_table(
    results, desired_order, control_group,
    caption = sprintf("Callaway--Sant'Anna Estimates of LIHTC Impacts on Test Scores, Bootstrapped SEs (%d Iterations), by Subgroup and Development Type (%s)", biters, cfg$label),
    label   = sprintf("tab:cs_bootstrap_%s", sample_key)
  )
  writeLines(tex_lines, file.path(base_path, sprintf("CS_ATT_%s_bootstrap.tex", sample_key)))
  message("Wrote: CS_ATT_", sample_key, "_bootstrap.csv / .tex")
}

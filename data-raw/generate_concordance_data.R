# Generate concordance test datasets for test-concordance.R
#
# Run this script once to regenerate inst/testdata/ datasets.
# Re-run only when test model parameters change.
#
# Requires: ferxtranslate (installed), ferx (installed)
# amp.sim dataset additionally requires: amp.sim (installed)
#
# True parameter values:
#   1cpt_oral.ctl           : TVCL=0.134, TVV=8.1, TVKA=1.0 (theta initials)
#   2cpt_iv.ctl             : TVCL=5.0,   TVV1=20.0, TVQ=8.0, TVV2=60.0
#   pk_1cmt_oral_ampsim.ctl : KA=0.0825, CL=2.676, V=1.588 (amp.sim NONMEM reference)

library(ferxtranslate)
library(ferx)

dir.create("inst/testdata", showWarnings = FALSE, recursive = TRUE)

# Helper: translate a bundled NONMEM .ctl and write to a temp .ferx file
translate_tmp <- function(model_name) {
  ctl <- system.file(file.path("testmodels/nonmem", model_name),
                     package = "ferxtranslate")
  result <- nm_to_ferx(ctl)
  ferx_file <- tempfile(fileext = ".ferx")
  writeLines(result$ferx_text, ferx_file)
  ferx_file
}

# Helper: overwrite one theta's initial value in emitted .ferx text.
# The translator renames a theta that would shadow an identically named
# individual parameter (KA -> TVKA and so on), so match either spelling and keep
# whichever it actually emitted. Stops rather than silently leaving the value
# untouched if neither matches -- a no-op sub() here would quietly simulate from
# the wrong parameters.
set_theta <- function(txt, nm, value) {
  pat <- sprintf("theta (TV)?%s\\([^)]+\\)", nm)
  hit <- regmatches(txt, regexpr(pat, txt))
  if (length(hit) == 0L)
    stop("no theta matching '", nm, "' in the emitted .ferx -- ",
         "has the translator's theta naming changed?", call. = FALSE)
  actual <- sub(sprintf("theta ((TV)?%s)\\(.*", nm), "\\1", hit)
  sub(pat, sprintf("theta %s(%.10g, 0.0, 1e15)", actual, value), txt)
}

# Helper: substitute one line, refusing to no-op. Same reasoning as set_theta():
# a sub() that silently fails to match would regenerate the dataset at the
# model's own initial variances instead of the amp.sim NONMEM reference, and the
# concordance test would then validate against the wrong variance structure and
# still pass.
set_line <- function(txt, pattern, replacement, what) {
  if (!grepl(pattern, txt))
    stop("no line matching '", what, "' in the emitted .ferx -- ",
         "has the translator's naming changed?", call. = FALSE)
  sub(pattern, replacement, txt)
}

# Helper: build a standard NONMEM-format dosing+observation template.
# DV carries a numeric 0 placeholder, not ".". ferx 0.2.0 treats a "." DV as
# missing and drops the row from the ferx_simulate() result, so a "." template
# yields zero simulated rows. The placeholder does not reach the output: DV is
# overwritten with the simulated value on observation rows and restored to "."
# on dose rows by simulate_dataset() below.
nm_template <- function(n_subj, dose, cmt, obs_times) {
  rows <- vector("list", n_subj * (length(obs_times) + 1))
  i <- 1L
  for (id in seq_len(n_subj)) {
    rows[[i]] <- data.frame(ID=id, TIME=0, DV=0, EVID=1L, AMT=dose,
                            CMT=cmt, MDV=1L)
    i <- i + 1L
    for (t in obs_times) {
      rows[[i]] <- data.frame(ID=id, TIME=t, DV=0, EVID=0L, AMT=".",
                              CMT=cmt, MDV=0L)
      i <- i + 1L
    }
  }
  do.call(rbind, rows)
}

# Helper: simulate observations for one template and write the dataset.
# Asserts that ferx returned one simulated value per observation row -- the
# silent failure mode this replaces was ferx_simulate() returning zero rows,
# which surfaced only as an obscure recycling error further down.
simulate_dataset <- function(ferx_file, tmpl, seed, out_path) {
  tf <- tempfile(fileext = ".csv")
  write.csv(tmpl, tf, row.names = FALSE, quote = FALSE)
  sim <- ferx_simulate(ferx_file, tf, n_sim = 1L, seed = seed)

  obs <- tmpl[tmpl$EVID == 0, ]
  if (nrow(sim) != nrow(obs))
    stop("ferx_simulate() returned ", nrow(sim), " rows for ", nrow(obs),
         " observation rows (", basename(out_path), ") -- the simulation did ",
         "not run over the template as expected", call. = FALSE)

  obs$DV <- round(sim$DV_SIM, 6)
  final <- rbind(tmpl[tmpl$EVID == 1, ], obs)
  final <- final[order(final$ID, final$TIME), ]
  # Dose rows carry no observation; NONMEM convention writes them as ".".
  final$DV[final$EVID == 1] <- "."
  rownames(final) <- NULL
  write.csv(final, out_path, row.names = FALSE, quote = FALSE)
  message("Written ", out_path, " (", nrow(final), " rows, ",
          length(unique(final$ID)), " subjects)")
  invisible(final)
}

# ---- 1-cpt oral (100 subjects, proportional error) -------------------------
simulate_dataset(
  translate_tmp("1cpt_oral.ctl"),
  nm_template(100, dose = 1.0, cmt = 1L,
              obs_times = c(0.25, 0.5, 1, 2, 4, 6, 8, 12, 16, 24)),
  seed = 123L, out_path = "inst/testdata/1cpt_oral_concordance.csv")

# ---- 2-cpt IV bolus (50 subjects, proportional error) ----------------------
simulate_dataset(
  translate_tmp("2cpt_iv.ctl"),
  nm_template(50, dose = 100.0, cmt = 1L,
              obs_times = c(0.1, 0.25, 0.5, 1, 2, 4, 6, 8, 12, 24, 36, 48)),
  seed = 456L, out_path = "inst/testdata/2cpt_iv_concordance.csv")

# ---- amp.sim 1-cpt oral (50 subjects, reference params from NONMEM .ext) ---
# Requires amp.sim (GitHub: LeidenAdvancedPKPD/amp.sim).
# True values come from amp.sim's published NONMEM FOCEI run on NM.theoph.02B.csv.
# NM.theoph.02B.csv is not bundled in amp.sim, so we simulate from the reference
# parameter values and use those simulated observations as the concordance dataset.
if (requireNamespace("amp.sim", quietly = TRUE)) {
  ext  <- read.table(
    system.file("example_models/PK.1CMT.ORAL.ext", package = "amp.sim"),
    header = TRUE, skip = 1)
  ref  <- ext[ext$ITERATION == -1000000000, ]

  ferx3_base <- translate_tmp("pk_1cmt_oral_ampsim.ctl")
  ferx3_txt  <- readLines(ferx3_base)
  ferx3_txt  <- paste(ferx3_txt, collapse = "\n")
  ferx3_txt  <- set_theta(ferx3_txt, "KA", ref$THETA1)
  ferx3_txt  <- set_theta(ferx3_txt, "CL", ref$THETA2)
  ferx3_txt  <- set_theta(ferx3_txt, "V",  ref$THETA3)
  ferx3_txt  <- set_line(ferx3_txt, "omega ETA_KA ~ [^\n]+",
                         sprintf("omega ETA_KA ~ %.10g", ref$OMEGA.1.1.), "omega ETA_KA")
  ferx3_txt  <- set_line(ferx3_txt, "omega ETA_CL ~ [^\n]+",
                         sprintf("omega ETA_CL ~ %.10g", ref$OMEGA.2.2.), "omega ETA_CL")
  ferx3_txt  <- set_line(ferx3_txt, "sigma EPS1 ~ [^\n]+",
                         sprintf("sigma EPS1 ~ %.10g (sd)", sqrt(ref$SIGMA.1.1.)),
                         "sigma EPS1")
  ferx3_sim  <- tempfile(fileext = ".ferx")
  writeLines(ferx3_txt, ferx3_sim)

  simulate_dataset(
    ferx3_sim,
    nm_template(50, dose = 4.0, cmt = 1L,
                obs_times = c(0.25, 0.5, 1, 2, 4, 6, 8, 12, 16, 24)),
    seed = 789L, out_path = "inst/testdata/ampsim_1cpt_oral_concordance.csv")
} else {
  message("amp.sim not installed -- skipping ampsim_1cpt_oral_concordance.csv")
}

# ---- ODE 1-cpt oral with S2=V scaling (50 subjects, proportional error) ----
# Tests that [scaling] obs_scale=V is correctly applied: data are concentrations,
# ODE predicts amounts, ferx divides by V before comparing to DV.
simulate_dataset(
  translate_tmp("pk_1cmt_oral.mod"),
  nm_template(50, dose = 4.0, cmt = 1L,
              obs_times = c(0.25, 0.5, 1, 2, 4, 6, 8, 12, 16, 24)),
  seed = 321L, out_path = "inst/testdata/ode_1cpt_oral_concordance.csv")

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

# Helper: build a standard NONMEM-format dosing+observation template
nm_template <- function(n_subj, dose, cmt, obs_times) {
  rows <- vector("list", n_subj * (length(obs_times) + 1))
  i <- 1L
  for (id in seq_len(n_subj)) {
    rows[[i]] <- data.frame(ID=id, TIME=0, DV=".", EVID=1L, AMT=dose,
                            CMT=cmt, MDV=1L)
    i <- i + 1L
    for (t in obs_times) {
      rows[[i]] <- data.frame(ID=id, TIME=t, DV=".", EVID=0L, AMT=".",
                              CMT=cmt, MDV=0L)
      i <- i + 1L
    }
  }
  do.call(rbind, rows)
}

# ---- 1-cpt oral (100 subjects, proportional error) -------------------------
ferx1 <- translate_tmp("1cpt_oral.ctl")
tmpl1 <- nm_template(100, dose=1.0, cmt=1L,
                     obs_times=c(0.25, 0.5, 1, 2, 4, 6, 8, 12, 16, 24))
tf1 <- tempfile(fileext = ".csv")
write.csv(tmpl1, tf1, row.names=FALSE, quote=FALSE)
sim1 <- ferx_simulate(ferx1, tf1, n_sim=1L, seed=123L)

obs1 <- tmpl1[tmpl1$EVID == 0, ]
obs1$DV <- round(sim1$DV_SIM, 6)
final1 <- rbind(tmpl1[tmpl1$EVID == 1, ], obs1)
final1 <- final1[order(final1$ID, final1$TIME), ]
rownames(final1) <- NULL
write.csv(final1, "inst/testdata/1cpt_oral_concordance.csv",
          row.names=FALSE, quote=FALSE)
message("Written inst/testdata/1cpt_oral_concordance.csv (",
        nrow(final1), " rows, ", length(unique(final1$ID)), " subjects)")

# ---- 2-cpt IV bolus (50 subjects, proportional error) ----------------------
ferx2 <- translate_tmp("2cpt_iv.ctl")
tmpl2 <- nm_template(50, dose=100.0, cmt=1L,
                     obs_times=c(0.1, 0.25, 0.5, 1, 2, 4, 6, 8, 12, 24, 36, 48))
tf2 <- tempfile(fileext = ".csv")
write.csv(tmpl2, tf2, row.names=FALSE, quote=FALSE)
sim2 <- ferx_simulate(ferx2, tf2, n_sim=1L, seed=456L)

obs2 <- tmpl2[tmpl2$EVID == 0, ]
obs2$DV <- round(sim2$DV_SIM, 6)
final2 <- rbind(tmpl2[tmpl2$EVID == 1, ], obs2)
final2 <- final2[order(final2$ID, final2$TIME), ]
rownames(final2) <- NULL
write.csv(final2, "inst/testdata/2cpt_iv_concordance.csv",
          row.names=FALSE, quote=FALSE)
message("Written inst/testdata/2cpt_iv_concordance.csv (",
        nrow(final2), " rows, ", length(unique(final2$ID)), " subjects)")

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
  ferx3_txt  <- sub("theta KA[(][^)]+[)]",
                    sprintf("theta KA(%.10g, 0.0, 1e15)", ref$THETA1), ferx3_txt)
  ferx3_txt  <- sub("theta CL[(][^)]+[)]",
                    sprintf("theta CL(%.10g, 0.0, 1e15)", ref$THETA2), ferx3_txt)
  ferx3_txt  <- sub("theta V[(][^)]+[)]",
                    sprintf("theta V(%.10g, 0.0, 1e15)",  ref$THETA3), ferx3_txt)
  ferx3_txt  <- sub("omega ETA_KA ~ [^\n]+",
                    sprintf("omega ETA_KA ~ %.10g", ref$OMEGA.1.1.), ferx3_txt)
  ferx3_txt  <- sub("omega ETA_CL ~ [^\n]+",
                    sprintf("omega ETA_CL ~ %.10g", ref$OMEGA.2.2.), ferx3_txt)
  ferx3_txt  <- sub("sigma EPS1 ~ [^\n]+",
                    sprintf("sigma EPS1 ~ %.10g (sd)", sqrt(ref$SIGMA.1.1.)), ferx3_txt)
  ferx3_sim  <- tempfile(fileext = ".ferx")
  writeLines(ferx3_txt, ferx3_sim)

  tmpl3 <- nm_template(50, dose = 4.0, cmt = 1L,
                       obs_times = c(0.25, 0.5, 1, 2, 4, 6, 8, 12, 16, 24))
  tf3   <- tempfile(fileext = ".csv")
  write.csv(tmpl3, tf3, row.names = FALSE, quote = FALSE)
  sim3  <- ferx_simulate(ferx3_sim, tf3, n_sim = 1L, seed = 789L)

  obs3   <- tmpl3[tmpl3$EVID == 0, ]
  obs3$DV <- round(sim3$DV_SIM, 6)
  final3 <- rbind(tmpl3[tmpl3$EVID == 1, ], obs3)
  final3 <- final3[order(final3$ID, final3$TIME), ]
  rownames(final3) <- NULL
  write.csv(final3, "inst/testdata/ampsim_1cpt_oral_concordance.csv",
            row.names = FALSE, quote = FALSE)
  message("Written inst/testdata/ampsim_1cpt_oral_concordance.csv (",
          nrow(final3), " rows, ", length(unique(final3$ID)), " subjects)")
} else {
  message("amp.sim not installed -- skipping ampsim_1cpt_oral_concordance.csv")
}

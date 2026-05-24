# Generate concordance test datasets for test-concordance.R
#
# Run this script once to regenerate inst/testdata/ datasets.
# Re-run only when test model parameters change.
#
# Requires: ferxtranslate (installed), ferx (installed)
#
# True parameter values:
#   1cpt_oral.ctl  : TVCL=0.134, TVV=8.1, TVKA=1.0 (theta initials)
#   2cpt_iv.ctl    : TVCL=5.0,   TVV1=20.0, TVQ=8.0, TVV2=60.0

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

# Simulated arm-level Emax meta-analysis for the nlmixr2 compartment-free
# concordance test (translator#74). Adapted from Bill Denney's example.R in
# FeRx-NLME/ferx-r#721 (comment 6727206415), used with his permission.
#
# 8 studies, 3 arms each (placebo and two doses), 5 visits. The arm mean is
# observed with standard error sd / sqrt(N), and DV is that mean times sqrt(W),
# W = N / sd^2 the inverse variance -- the hand-built weighting the model's
# `pred <- f * sqrt(W)` undoes.
#
# Writes:
#   inst/testdata/mbma_emax_nlmixr2.csv          the data (DV pre-multiplied)
#   inst/testdata/mbma_emax_nlmixr2_ref.csv      nlmixr2 FOCEi estimates + objf
#
# Needs nlmixr2est. Re-run if the model below or the seed changes.

set.seed(20261108)
arms <- expand.grid(ARM = 1:3, ID = 1:8)
arms$DOSE <- c(0, 50, 200)[arms$ARM]
arms$N <- sample(40:200, nrow(arms), replace = TRUE)
arms$OCC <- arms$ARM
d <- merge(arms, data.frame(TIME = c(0, 2, 4, 8, 12)))
d <- d[order(d$ID, d$OCC, d$TIME), ]
e0 <- 10 + stats::rnorm(8, 0, 1)[d$ID]
emax_arm <- 8 * exp(stats::rnorm(nrow(arms), 0, 0.1))[match(paste(d$ID, d$OCC), paste(arms$ID, arms$OCC))]
mean_true <- e0 - emax_arm * d$DOSE / (40 + d$DOSE) * d$TIME / (3 + d$TIME)
d$W <- d$N / 9
d$DV <- (mean_true + stats::rnorm(nrow(d), 0, 1 / sqrt(d$W))) * sqrt(d$W)
d$EVID <- 0L
d$MDV <- 0L
d <- d[, c("ID", "TIME", "DV", "EVID", "MDV", "OCC", "ARM", "DOSE", "N", "W")]
rownames(d) <- NULL

data_path <- file.path("inst", "testdata", "mbma_emax_nlmixr2.csv")
utils::write.csv(d, data_path, row.names = FALSE)
# Fit the data as read back from the CSV, so nlmixr2 and ferx see the same
# rounded digits.
d <- utils::read.csv(data_path)

source(file.path("inst", "testmodels", "nlmixr2", "mbma_emax_nlmixr2.R"))
fit <- nlmixr2est::nlmixr2(mbma_emax, d, est = "focei",
                           control = nlmixr2est::foceiControl(print = 0))
# fit$omega is one matrix per level (id, OCC) when the model has IOV.
lv <- if (is.list(fit$omega)) unname(fit$omega) else list(fit$omega)
om <- unlist(lapply(lv, function(m) stats::setNames(diag(m), rownames(m))))
ref <- data.frame(
  name  = c(names(fit$theta), names(om), "objf"),
  value = c(unname(fit$theta), unname(om), fit$objf)
)
# 17 significant digits: the test compares OFVs to 1e-3, and a rounded estimate
# moves the OFV by more than that.
ref$value <- trimws(formatC(ref$value, digits = 17, format = "g"))
utils::write.csv(ref, file.path("inst", "testdata", "mbma_emax_nlmixr2_ref.csv"),
                 row.names = FALSE, quote = FALSE)
message("nlmixr2 objf: ", fit$objf, " | ", paste(fit$message, collapse = " "))

# Arm-level Emax model-based meta-analysis, compartment-free (translator#74).
# Adapted from Bill Denney's example.R in FeRx-NLME/ferx-r#721, used with his
# permission. DV is the arm mean times sqrt(W), W the inverse variance, so the
# prediction is multiplied by sqrt(W) too.
mbma_emax <- function() {
  ini({
    tve0 <- 9
    tvemax <- 6
    led50 <- log(30)
    lt50 <- log(2)
    eta.e0 ~ 1
    iov.emax ~ 0.05 | OCC
    sdw <- 1
  })
  model({
    E0 <- tve0 + eta.e0
    EMAX <- tvemax * exp(iov.emax)
    ONSET <- 0
    if (TIME > 0) ONSET <- TIME / (exp(lt50) + TIME)
    EFF <- EMAX * DOSE / (exp(led50) + DOSE) * ONSET
    pred <- (E0 - EFF) * sqrt(W)
    pred ~ add(sdw)
  })
}

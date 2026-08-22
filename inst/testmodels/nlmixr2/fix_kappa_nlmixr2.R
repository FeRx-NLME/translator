# IOV kappas, one FIXed and one free -- nlmixr2 / rxode2 function (issue #31,
# review of PR #38).
#
# The kappa FIX path is reachable only from an nlmixr2 source. nonmem2rx reads
# NONMEM's ETA-coded IOV as ordinary IIV (see .iov_flattening_warnings()), so a
# control stream never produces a `condition != "id"` row and never reaches
# .extract_kappas() at all.
#
# TWO kappas, one fixed and one not: a single fixed kappa cannot distinguish
# "FIX read from iniDf$fix" from "FIX appended to every kappa".
#
# Note the flag MUST come from iniDf$fix here. Measured on rxode2 5.1.2,
# attr(ui$omega, "lotriFix") is NULL for any model carrying IOV, so the matrix
# attribute that serves the NONMEM sigma path answers FALSE for every kappa.

f_fix_kappa <- function() {
  ini({
    tvcl <- c(0.001, 0.134, 10.0)
    tvv  <- c(0.1,   8.1,  500.0)

    eta.cl ~ 0.07

    iov.cl ~ fix(0.05) | occ
    iov.v  ~ 0.03      | occ

    prop.err <- 0.01
  })
  model({
    cl <- tvcl * exp(eta.cl + iov.cl)
    v  <- tvv * exp(iov.v)

    linCmt() ~ prop(prop.err)
  })
}

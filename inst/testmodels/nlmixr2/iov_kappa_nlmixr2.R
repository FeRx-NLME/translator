# 1-cpt oral with a KAPPA-named IIV eta -- nlmixr2 / rxode2 function.
#
# Used to assert the IOV-flattening warning is NOT emitted for non-NONMEM
# sources. The warning names nonmem2rx and targets its specific habit of reading
# ETA-coded IOV as IIV; here `kappa.cl` is simply an ordinary IIV eta
# (condition "id"), so it must NOT trigger the NONMEM-only flattening warning.

f_iov_kappa <- function() {
  ini({
    tvcl <- c(0.001, 0.134, 10.0)
    tvv  <- c(0.1,   8.1,  500.0)
    tvka <- c(0.01,  1.0,   50.0)

    eta.cl   ~ 0.07
    kappa.cl ~ 0.04

    prop.err <- 0.01
  })
  model({
    cl <- tvcl * exp(eta.cl + kappa.cl)
    v  <- tvv
    ka <- tvka

    linCmt() ~ prop(prop.err)
  })
}

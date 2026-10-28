# Correlated IOV on CL and V -- nlmixr2 / rxode2 function (issue #40).
#
# An IOV block is reachable only from an nlmixr2 source: nonmem2rx reads NONMEM's
# ETA-coded IOV as ordinary IIV, so a control stream never yields a
# `condition != "id"` row. The covariance is non-zero so that "block kept" and
# "block flattened to its diagonal" are different models.

f_iov_block <- function() {
  ini({
    tvcl <- c(0.001, 0.134, 10.0)
    tvv  <- c(0.1,   8.1,  500.0)
    tvka <- c(0.01,  1.0,   50.0)

    eta.cl ~ 0.07

    iov.cl + iov.v ~ c(0.05,
                       0.01, 0.03) | occ

    prop.err <- 0.01
  })
  model({
    cl <- tvcl * exp(eta.cl + iov.cl)
    v  <- tvv * exp(iov.v)
    ka <- tvka

    linCmt() ~ prop(prop.err)
  })
}

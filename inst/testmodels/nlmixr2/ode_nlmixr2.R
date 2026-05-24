# One-compartment oral PK model (ODE form) -- nlmixr2 / rxode2 function
# Used by ferxtranslate integration tests via rxui_to_ir(rxode2::rxode2(f_ode))

f_ode_oral <- function() {
  ini({
    tvcl <- c(0.001, 0.134, 10.0)
    tvv  <- c(0.1,   8.1,  500.0)
    tvka <- c(0.01,  1.0,   50.0)

    eta.cl ~ 0.07
    eta.v  ~ 0.02
    eta.ka ~ 0.40

    prop.err <- 0.01
  })
  model({
    cl <- tvcl * exp(eta.cl)
    v  <- tvv  * exp(eta.v)
    ka <- tvka * exp(eta.ka)

    d/dt(depot)   = -ka * depot
    d/dt(central) =  ka * depot / v - (cl / v) * central

    central ~ prop(prop.err)
  })
}

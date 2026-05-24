# 1-cpt oral NONMEM: snapshot + no unsupported

    Code
      cat(norm_snap(result$ferx_text))
    Output
      # Translated from nonmem: 1cpt_oral.ctl
      
      [parameters]
        theta TVCL(0.134, 0.001, 10.0)
        theta TVV(8.1, 0.1, 500.0)
        theta TVKA(1.0, 0.01, 50.0)
      
        omega ETA_CL ~ 0.07
        omega ETA_V ~ 0.02
        omega ETA_KA ~ 0.4
      
        sigma EPS1 ~ 0.1 (sd)
      
      [individual_parameters]
        CL = TVCL * exp(ETA_CL)
        V = TVV * exp(ETA_V)
        KA = TVKA * exp(ETA_KA)
      
      [structural_model]
        pk one_cpt_oral(cl=CL, v=V, ka=KA)
      
      [error_model]
        DV ~ proportional(EPS1)
      
      [fit_options]
        method = focei
        maxiter = 500
        covariance = true

# 2-cpt oral with covariates: snapshot + no unsupported

    Code
      cat(norm_snap(result$ferx_text))
    Output
      # Translated from nonmem: 2cpt_oral_cov.ctl
      
      [parameters]
        theta TVCL(5.0, 0.1, 100.0)
        theta TVV1(50.0, 1.0, 500.0)
        theta TVQ(10.0, 0.1, 100.0)
        theta TVV2(100.0, 1.0, 500.0)
        theta TVKA(1.2, 0.01, 10.0)
        theta THETA_WT(0.75, 0.01, 5.0)
        theta THETA_CRCL(0.5, 0.01, 5.0)
      
        omega ETA_CL ~ 0.1
        omega ETA_V1 ~ 0.1
        omega ETA_Q ~ 0.05
        omega ETA_V2 ~ 0.05
        omega ETA_KA ~ 0.15
      
        sigma EPS1 ~ 0.14142135623731 (sd)
      
      [individual_parameters]
        CL = TVCL * (WT/70)^THETA_WT * (CRCL/100)^THETA_CRCL * exp(ETA_CL)
        V2 = TVV1 * (WT/70)^THETA_WT * exp(ETA_V1)
        Q = TVQ * exp(ETA_Q)
        V3 = TVV2 * exp(ETA_V2)
        KA = TVKA * exp(ETA_KA)
      
      [structural_model]
        pk two_cpt_oral(cl=CL, v1=V2, q=Q, v2=V3, ka=KA)
      
      [error_model]
        DV ~ proportional(EPS1)
      
      [fit_options]
        method = focei
        maxiter = 500
        covariance = true

# 2-cpt IV bolus: infers two_cpt_iv_bolus

    Code
      cat(norm_snap(result$ferx_text))
    Output
      # Translated from nonmem: 2cpt_iv.ctl
      
      [parameters]
        theta TVCL(5.0, 0.1, 100.0)
        theta TVV1(20.0, 1.0, 500.0)
        theta TVQ(8.0, 0.1, 100.0)
        theta TVV2(60.0, 1.0, 500.0)
      
        omega ETA_CL ~ 0.1
        omega ETA_V1 ~ 0.1
        omega ETA_Q ~ 0.08
        omega ETA_V2 ~ 0.05
      
        sigma EPS1 ~ 0.14142135623731 (sd)
      
      [individual_parameters]
        CL = TVCL * exp(ETA_CL)
        V1 = TVV1 * exp(ETA_V1)
        Q = TVQ * exp(ETA_Q)
        V2 = TVV2 * exp(ETA_V2)
      
      [structural_model]
        pk two_cpt_iv_bolus(cl=CL, v1=V1, q=Q, v2=V2)
      
      [error_model]
        DV ~ proportional(EPS1)
      
      [fit_options]
        method = focei
        maxiter = 500
        covariance = true

# ODE warfarin: full $DES path, [odes] section present

    Code
      cat(norm_snap(result$ferx_text))
    Output
      # Translated from nonmem: ode_warfarin.ctl
      # Warnings: 1 -- run result$warnings for details
      
      [parameters]
        theta TVCL(0.134, 0.001, 10.0)
        theta TVV(8.1, 0.1, 500.0)
        theta TVKA(1.0, 0.01, 50.0)
      
        omega ETA_CL ~ 0.07
        omega ETA_V ~ 0.02
        omega ETA_KA ~ 0.4
      
        sigma EPS1 ~ 0.1 (sd)
      
      [individual_parameters]
        CL = TVCL * exp(ETA_CL)
        V = TVV * exp(ETA_V)
        KA = TVKA * exp(ETA_KA)
      
      [structural_model]
        ode(obs_cmt=CENTRAL, states=[DEPOT, CENTRAL])
      
      [odes]
        d/dt(DEPOT) = -KA * DEPOT
        d/dt(CENTRAL) = KA * DEPOT/V - (CL/V) * CENTRAL
      
      [error_model]
        DV ~ proportional(EPS1)
      
      [fit_options]
        method = focei
        maxiter = 500
        covariance = true

# block omega: block_omega line in output

    Code
      cat(norm_snap(result$ferx_text))
    Output
      # Translated from nonmem: block_omega.ctl
      
      [parameters]
        theta TVCL(0.134, 0.001, 10.0)
        theta TVV(8.1, 0.1, 500.0)
        theta TVKA(1.0, 0.01, 50.0)
      
        block_omega (ETA_CL, ETA_V) = [0.07, 0.02, 0.02]
        omega ETA_KA ~ 0.4
      
        sigma EPS1 ~ 0.1 (sd)
      
      [individual_parameters]
        CL = TVCL * exp(ETA_CL)
        V = TVV * exp(ETA_V)
        KA = TVKA * exp(ETA_KA)
      
      [structural_model]
        pk one_cpt_oral(cl=CL, v=V, ka=KA)
      
      [error_model]
        DV ~ proportional(EPS1)
      
      [fit_options]
        method = focei
        maxiter = 500
        covariance = true

# IOV model: translates without error; KAPPA_CL emitted as omega (nonmem2rx treats IOV as IIV)

    Code
      cat(norm_snap(result$ferx_text))
    Output
      # Translated from nonmem: iov.ctl
      
      [parameters]
        theta TVCL(0.134, 0.001, 10.0)
        theta TVV(8.1, 0.1, 500.0)
        theta TVKA(1.0, 0.01, 50.0)
      
        omega ETA_CL ~ 0.07
        omega ETA_V ~ 0.02
        omega ETA_KA ~ 0.4
        omega KAPPA_CL ~ 0.04
      
        sigma EPS1 ~ 0.1 (sd)
      
      [individual_parameters]
        RXM_KAPPA_CL = KAPPA_CL
        CL = TVCL * exp(ETA_CL + RXM_KAPPA_CL)
        V = TVV * exp(ETA_V)
        KA = TVKA * exp(ETA_KA)
      
      [structural_model]
        pk one_cpt_oral(cl=CL, v=V, ka=KA)
      
      [error_model]
        DV ~ proportional(EPS1)
      
      [fit_options]
        method = focei
        maxiter = 500
        covariance = true

# 1-cpt oral nlmixr2: snapshot + one_cpt_oral

    Code
      cat(norm_snap(result$ferx_text))
    Output
      # Translated from nlmixr2: unknown
      
      [parameters]
        theta TVCL(0.134, 0.001, 10.0)
        theta TVV(8.1, 0.1, 500.0)
        theta TVKA(1.0, 0.01, 50.0)
      
        omega ETA_CL ~ 0.07
        omega ETA_V ~ 0.02
        omega ETA_KA ~ 0.4
      
        sigma PROP_ERR ~ 0.01 (sd)
      
      [individual_parameters]
        CL = TVCL * exp(ETA_CL)
        V = TVV * exp(ETA_V)
        KA = TVKA * exp(ETA_KA)
      
      [structural_model]
        pk one_cpt_oral(cl=CL, v=V, ka=KA)
      
      [error_model]
        DV ~ proportional(PROP_ERR)
      
      [fit_options]
        method = focei
        maxiter = 500
        covariance = true

# ODE nlmixr2: d/dt expressions produce [odes] section

    Code
      cat(norm_snap(result$ferx_text))
    Output
      # Translated from nlmixr2: unknown
      # Warnings: 1 -- run result$warnings for details
      
      [parameters]
        theta TVCL(0.134, 0.001, 10.0)
        theta TVV(8.1, 0.1, 500.0)
        theta TVKA(1.0, 0.01, 50.0)
      
        omega ETA_CL ~ 0.07
        omega ETA_V ~ 0.02
        omega ETA_KA ~ 0.4
      
        sigma PROP_ERR ~ 0.01 (sd)
      
      [individual_parameters]
        CL = TVCL * exp(ETA_CL)
        V = TVV * exp(ETA_V)
        KA = TVKA * exp(ETA_KA)
      
      [structural_model]
        ode(obs_cmt=central, states=[depot, central])
      
      [odes]
        d/dt(depot) = -KA * depot
        d/dt(central) = KA * depot/V - (CL/V) * central
      
      [error_model]
        DV ~ proportional(PROP_ERR)
      
      [fit_options]
        method = focei
        maxiter = 500
        covariance = true

# amp.sim 1-cpt oral ODE: [odes] section + obs_cmt inferred

    Code
      cat(norm_snap(result$ferx_text))
    Output
      # Translated from nonmem: pk_1cmt_oral.mod
      # Warnings: 1 -- run result$warnings for details
      
      [parameters]
        theta KA(0.1, 0.0, 1e15)
        theta CL(2.0, 0.0, 1e15)
        theta V(1.0, 0.0, 1e15)
      
        omega ETA_KA ~ 0.01
        omega ETA_CL ~ 0.02
      
        sigma EPS1 ~ 0.316227766016838 (sd)
      
      [individual_parameters]
        KA = KA * exp(ETA_KA)
        CL = CL * exp(ETA_CL)
        K20 = CL/V
      
      [structural_model]
        ode(obs_cmt=CENTRAL, states=[ABS, CENTRAL])
      
      [odes]
        d/dt(ABS) = -KA * ABS
        d/dt(CENTRAL) = KA * ABS - K20 * CENTRAL
      
      [error_model]
        DV ~ proportional(EPS1)
      
      [fit_options]
        method = focei
        maxiter = 500
        covariance = true

# amp.sim PKPD indirect response: 4-state ODE + additive error

    Code
      cat(norm_snap(result$ferx_text))
    Output
      # Translated from nonmem: pkpd_ir.mod
      # Warnings: 1 -- run result$warnings for details
      
      [parameters]
        theta TH1(1.95, -1e15, 1e15)
        theta TH2(1.0, -1e15, 1e15)
        theta TH3(1.31, -1e15, 1e15)
        theta TH4(4.15, -1e15, 1e15)
        theta TH5(0.904, -1e15, 1e15)
        theta TH6(100.0, 0.0, 1e15)
        theta TH7(3.94, 0.0, 1e15)
        theta TH8(0.9, 0.0, 1.0)
        theta TH9(2.0, 0.0, 1e15)
      
        omega ETA1 ~ 0.1
      
        sigma EPS1 ~ 5.47722557505166 (sd)
      
      [individual_parameters]
        K12 = TH1
        CL = TH2
        V2 = TH3
        V3 = TH4
        Q = TH5
        K20 = CL/V2
        K23 = Q/V2
        K32 = Q/V3
        BL = TH6 * exp(ETA1)
        KOUT = TH7
        RIN = BL * KOUT
        EMAX = TH8
        EC50 = TH9
      
      [structural_model]
        ode(obs_cmt=EFFECT, states=[DOSE, CENTRAL, PERIPH, EFFECT])
      
      [odes]
        d/dt(DOSE) = -K12 * DOSE
        d/dt(CENTRAL) = K12 * DOSE - K20 * CENTRAL - K23 * CENTRAL + K32 * PERIPH
        d/dt(PERIPH) = K23 * CENTRAL - K32 * PERIPH
        d/dt(EFFECT) = EFF * RIN - KOUT * EFFECT
      
      [error_model]
        DV ~ additive(EPS1)
      
      [fit_options]
        method = focei
        maxiter = 500
        covariance = true


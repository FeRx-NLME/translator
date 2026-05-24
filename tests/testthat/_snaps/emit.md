# 1-cpt oral IR emits correct .ferx (snapshot)

    Code
      cat(emit_ferx(warfarin_1cpt_ir()))
    Output
      # Translated from nonmem: 1cpt_oral.ctl
      
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
        V = TVV  * exp(ETA_V)
        KA = TVKA * exp(ETA_KA)
      
      [structural_model]
        pk one_cpt_oral(cl=CL, v=V, ka=KA)
      
      [error_model]
        DV ~ proportional(PROP_ERR)
      
      [fit_options]
        method = foce
        maxiter = 300
        covariance = true

# 2-cpt oral with covariates emits correct .ferx (snapshot)

    Code
      cat(emit_ferx(ir))
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
      
        sigma PROP_ERR ~ 0.02 (sd)
      
      [individual_parameters]
        CL = TVCL * (WT / 70)^THETA_WT * (CRCL / 100)^THETA_CRCL * exp(ETA_CL)
        V1 = TVV1 * (WT / 70)^THETA_WT * exp(ETA_V1)
        Q = TVQ  * exp(ETA_Q)
        V2 = TVV2 * exp(ETA_V2)
        KA = TVKA * exp(ETA_KA)
      
      [structural_model]
        pk two_cpt_oral(cl=CL, v1=V1, q=Q, v2=V2, ka=KA)
      
      [error_model]
        DV ~ proportional(PROP_ERR)
      
      [fit_options]
        method = focei
        maxiter = 500
        covariance = true


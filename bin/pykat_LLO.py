class LLO:
    T_SRM  = 1
    L_SRM  = 0
    R_PRM  = 0.9689845 
    T_PRM  = 0.031
    T_ITMX = 0.0148
    L_ITMX = 175e-6
    T_ITMY = 0.0148
    L_ITMY = 265e-6
    nsilica = 1.44963098985906
    f1 = 9099471
    f2 = 45497355
    #TL_f 50k
    TL_f = "inf"
    PRM_RC = 11.009
        
    # IFO put on dark fringe
    phi_SRM = 90

    # (for phase 3 version use these tunings below)
    #phi_ITMX = 0
    #phi_ITMY = 0
    # (for phase 2 version use these tunings below)
    phi_ITMX = -1.7
    phi_ITMY = 1.7

    phi_BS = 0
    phi_ARBSX = 0
    phi_ARBSSR = 0
    phi_ETMX = 0
    phi_ETMY = 0

    # (for phase 3 version use these tunings below)
    #phi_PRM = 90
    # (for phase 2 version use these tunings below)
    phi_PRM = 91.54

    Laser = """
    l L0 1 0 n0
    s lmod1 1 n0 n1
    mod mod1 {f1} 0.22 1 pm n1 n2
    s lmod2 1 n2 n3
    mod mod2 {f2} 0.179 1 pm n3 nin
    """.format(f1=f1, f2=f2)

    PR = """
    s lin 1 nin nREFL	

    # Power recycling mirror PRM-02
    m PRM {R_PRM} {T_PRM} {phi_PRM} nREFL nPRMb
    attr PRM Rc {PRM_RC}

    # Distance between PRM and PR2
    s lp1 16.6107 nPRMb nPR2a

    # PR2 PR2-02
    bs1 PR2 243u 8.6u 0 -0.79 nPR2a nPR2b dump nPOP
    attr PR2 Rc -4.545

    # Distance from PR2 to PR3
    s lp2 16.1647 nPR2b nPR3a

    # PR3 PR3-03
    bs1 PR3 5.3u 17u 0 0.615 nPR3a nPR3b dump dump
    attr PR3 Rc 36.027

    # Distance from PR3 to BS
    s lp3 19.5381 nPR3b nHRBS_PR
    """.format(R_PRM=R_PRM, T_PRM=T_PRM, phi_PRM=phi_PRM, PRM_RC=PRM_RC)

    BS = """
    bs1 HRBS 0.5 8.6u {phi_BS} 45 nHRBS_PR nHRBS_Y nHRBS_X nHRBS_SR
    s sHRBStoARBSX 0.0685 {nsilica} nHRBS_X nARBSX_sub
    bs2 ARBSX 30u 1.7u {phi_ARBSX} 29.1951 nARBSX_sub dump nARBSX_X dump
    s sHRBStoARBSSR 0.0684 {nsilica} nHRBS_SR nARBSSR_sub
    bs2 ARBSSR 30u 1.7u {phi_ARBSSR} -29.1951 nARBSSR_sub dump nARBSSR_SR dump
    """.format(phi_BS=phi_BS, nsilica=nsilica, phi_ARBSSR=phi_ARBSSR, phi_ARBSX=phi_ARBSX)

    SR = """
    # Distance from BS to SR3
    s ls3 19.3661 nARBSSR_SR nSR3a

    # SR3 SR3-01
    bs1 SR3 5n 19.1n 0 0 nSR3a nSR3b dump dump
    attr SR3 Rc 35.97

    # Distance from SR3 to SR2
    s ls2 15.4435 nSR3b nSR2a

    # SR2 SR2-04
    bs1 SR2 7.6n 10.8n 0 0 nSR2a nSR2b dump dump
    attr SR2 Rc -6.406

    # Distance from SR2 to SRMHR
    s ls1 15.7566 nSR2b nSRMHRa

    # Signal recycling mirror SRM-08
    m1 SRMHR {T_SRM} {L_SRM} {phi_SRM} nSRMHRa nSRMHRb
    s SRMsub 0.0749 {nsilica} nSRMHRb nSRMARa
    attr SRMHR Rc -5.667
    m2 SRMAR 50n 0 {phi_SRM} nSRMARa nSRMARb

    # Output of interferometer
    #s lout 1 nSRMtrans nASpickin
    #bs ASpick 0.1 0.9 0 45 nASpickin nASport nASpicktrans dump
    """.format(T_SRM=T_SRM, L_SRM=L_SRM, phi_SRM=phi_SRM,nsilica=nsilica)

    YARM = """
    # Using values from E1200274
    s ly1 4.847 nHRBS_Y nCPYar1

    # Y arm compensation plate CP-08
    m2 CPYar1 48.9u 0.4u 0 nCPYar1 nCPYar1s
    s sCPY 0.1 {nsilica} nCPYar1s nCPYar2s
    m2 CPYar2 30.5u 0.3u 0 nCPYar2s nCPYar2
    s sCPYtoITMYar 0.02 nCPYar2 nITMYTLin

    # Y arm input mirror ITM-08

    lens ITMYTL {TL_f} nITMYTLin nITMYTLtrans
    s ITMYTL_null 0 nITMYTLtrans nITMYconstL_in
    lens ITMYconstL -120k nITMYconstL_in nITMYconstL_trans
    #lens ITMYconstL inf nITMYconstL_in nITMYconstL_trans
    s ITMYTL_null2 0 nITMYconstL_trans nITMYar_in
    m2 ITMYar 250u 0 0 nITMYar_in nITMYs1
    s lITMY 0.2 {nsilica} nITMYs1 nITMYs2
    m1 ITMY {T_ITMY} {L_ITMY} {phi_ITMY} nITMYs2 nITMY2
    attr ITMY Rc -1940.7
    """.format(nsilica=nsilica, T_ITMY=T_ITMY, TL_f=TL_f, L_ITMY=L_ITMY, phi_ITMY=phi_ITMY)

    XARM = """
    # Now using length taken from E1200616
    s lx1 4.829 nARBSX_X nCPXar1

    # X arm compensation plate CP-06 (no values for AR reflection so using same as Yarm)
    m2 CPXar1 48.9u 4.3u 0 nCPXar1 nCPXar1s
    s sCPX 0.1 {nsilica} nCPXar1s nCPXar2s
    m2 CPXar2 30.5u 4.8u 0 nCPXar2s nCPXar2
    s sCPXtoITMXar 0.02 nCPXar2 nITMXTLin

    # X arm input mirror ITM-04
    lens ITMXTL {TL_f} nITMXTLin nITMXTLtrans
    s ITMXtl_null 0 nITMXTLtrans nITMXconstL_in
    #lens ITMXconstL inf nITMXconstL_in nITMXconstL_trans
    lens ITMXconstL 442k nITMXconstL_in nITMXconstL_trans
    s ITMXTL_null2 0 nITMXconstL_trans nITMXar_in
    m2 ITMXar 164u 0 0 nITMXar_in nITMXs1
    s lITMX1 0.2 {nsilica} nITMXs1 nITMXs2
    m1 ITMX {T_ITMX} {L_ITMX} {phi_ITMX} nITMXs2 nITMX2

    # default Rc from nebula page
    attr ITMX Rc -1937.9

    # Rcs for looking at AS port beam shape with astigmatic RH effect
    #attr ITMX Rcx -1915
    #attr ITMX Rcy -1912.8
    """.format(nsilica=nsilica, TL_f=TL_f, L_ITMX=L_ITMX, phi_ITMX=phi_ITMX, T_ITMX=T_ITMX)


const trimer_rn = @reaction_network trimer_rn begin
    @parameters begin
        kfᴸᴬ::Float64 = $(DEFAULTS[:kfᴸᴬ]), [bounds=$(BOUNDS[:kfᴸᴬ])] 
        krᴸᴬ::Float64 = $(DEFAULTS[:krᴸᴬ]), [bounds=$(BOUNDS[:krᴸᴬ])] 
        kfᴸᴷ::Float64 = $(DEFAULTS[:kfᴸᴷ]), [bounds=$(BOUNDS[:kfᴸᴷ])] 
        krᴸᴷ::Float64 = $(DEFAULTS[:krᴸᴷ]), [bounds=$(BOUNDS[:krᴸᴷ])] 
        kcatᴸᴷ::Float64 = $(DEFAULTS[:kcatᴸᴷ]), [bounds=$(BOUNDS[:kcatᴸᴷ])]
        kfᴸᴾ::Float64 = $(DEFAULTS[:kfᴸᴾ]), [bounds=$(BOUNDS[:kfᴸᴾ])]
        krᴸᴾ::Float64 = $(DEFAULTS[:krᴸᴾ]), [bounds=$(BOUNDS[:krᴸᴾ])] 
        kcatᴸᴾ::Float64 = $(DEFAULTS[:kcatᴸᴾ]), [bounds=$(BOUNDS[:kcatᴸᴾ])] 
        kfᴬᴷ::Float64 = $(DEFAULTS[:kfᴬᴷ]), [bounds=$(BOUNDS[:kfᴬᴷ])] 
        krᴬᴷ::Float64 = $(DEFAULTS[:krᴬᴷ]), [bounds=$(BOUNDS[:krᴬᴷ])] 
        kfᴬᴾ::Float64 = $(DEFAULTS[:kfᴬᴾ]), [bounds=$(BOUNDS[:kfᴬᴾ])] 
        krᴬᴾ::Float64 = $(DEFAULTS[:krᴬᴾ]), [bounds=$(BOUNDS[:krᴬᴾ])]
        DF::Float64 = $(DEFAULTS[:DF]), [description="Unitless dimensionality factor", bounds=$(BOUNDS[:DF])]

        # adaptor-transferrin receptor binding, this couples the lipid oscillator to the hetrotrimer assembly network
        # kfᴬᵀ::Float64 = 0.001, [bounds = (1e-3, 1e2)] 
        # krᴬᵀ::Float64 = 0.001, [bounds = (1e-3, 1e3)] 

        # monomer-monomer binding
        kf¹ᵐ::Float64 = 0.001, [bounds = (1e-3, 1e2)] 
        kr¹ᵐ::Float64 = 0.001, [bounds = (1e-3, 1e3)] 
        kf²ᵐ::Float64 = 0.001, [bounds = (1e-3, 1e2)] 
        kr²ᵐ::Float64 = 0.001, [bounds = (1e-3, 1e3)] 
        kf³ᵐ::Float64 = 0.001, [bounds = (1e-3, 1e2)] 
        kr³ᵐ::Float64 = 0.001, [bounds = (1e-3, 1e3)] 
        kf⁴ᵐ::Float64 = 0.001, [bounds = (1e-3, 1e2)] 
        kr⁴ᵐ::Float64 = 0.001, [bounds = (1e-3, 1e3)] 
        kf⁵ᵐ::Float64 = 0.001, [bounds = (1e-3, 1e2)] 
        kr⁵ᵐ::Float64 = 0.001, [bounds = (1e-3, 1e3)] 
        kf⁶ᵐ::Float64 = 0.001, [bounds = (1e-3, 1e2)] 
        kr⁶ᵐ::Float64 = 0.001, [bounds = (1e-3, 1e3)] 
        kf⁷ᵐ::Float64 = 0.001, [bounds = (1e-3, 1e2)]
        kr⁷ᵐ::Float64 = 0.001, [bounds = (1e-3, 1e3)]
    end 

    @species begin
        L(t) = $(DEFAULTS[:L]), [description="PIP", bounds=$(BOUNDS[:L])] 
        K(t) = $(DEFAULTS[:K]), [description="PIP5K", bounds=$(BOUNDS[:K])] 
        P(t) = $(DEFAULTS[:P]), [description="Synaptojanin", bounds=$(BOUNDS[:P])] 
        A(t) = $(DEFAULTS[:A]), [description="AP2", bounds=$(BOUNDS[:A])] 
        Lp(t) = 0.0, [description="PIP2", tunable=false] 
        LpA(t) = 0.0, [description="PIP2-AP2", tunable=false] 
        LK(t) = 0.0, [description="PIP-Kinase", tunable=false] 
        LpP(t) = 0.0, [description="PIP2-Phosphatase", tunable=false] 
        LpAK(t) = 0.0, [description="PIP2-AP2-Kinase", tunable=false] 
        LpAP(t) = 0.0, [description="PIP2-AP2-Phosphatase", tunable=false] 
        LpAKL(t) = 0.0, [description="PIP2-AP2-Kinase-PIP", tunable=false] 
        LpAPLp(t) = 0.0, [description="PIP2-AP2-Phosphatase-PIP2", tunable=false] 
        AK(t) = 0.0, [description = "AP2-Kinase", tunable = false] 
        AP(t) = 0.0, [description = "AP2-Phosphatase", tunable = false] 
        AKL(t) = 0.0, [description = "AP2-Kinase-PIP", tunable = false] 
        APLp(t) = 0.0, [description = "AP2-Phosphatase-PIP2", tunable = false] 

        # T(t) = 1.0, [description = "Transferrin receptor", bounds = (1e-1, 1e2)]
        B(t) = 1.0, [description = "Trimer monomer B", bounds = (1e-1, 1e2)]
        C(t) = 1.0, [description = "Trimer monomer C", bounds = (1e-1, 1e2)]
        D(t) = 1.0, [description = "Trimer monomer D", bounds = (1e-1, 1e2)]
        # X(t) = 0.0, [description = "LpA-T complex", tunable = false]
        BC(t) = 0.0, [description = "BC dimer", tunable = false]
        BD(t) = 0.0, [description = "BD dimer", tunable = false]
        CD(t) = 0.0, [description = "CD dimer", tunable = false]
        XB(t) = 0.0, [description = "X-B complex", tunable = false]
        XC(t) = 0.0, [description = "X-C complex", tunable = false]
        XD(t) = 0.0, [description = "X-D complex", tunable = false]
        XBC(t) = 0.0, [description = "X-BC complex", tunable = false]
        XXBC(t) = 0.0, [description = "XX-BC complex", tunable = false]
        XBD(t) = 0.0, [description = "X-BD complex", tunable = false]
        XXBD(t) = 0.0, [description = "XX-BD complex", tunable = false]
        XCD(t) = 0.0, [description = "X-CD complex", tunable = false]
        XXCD(t) = 0.0, [description = "XX-CD complex", tunable = false]
        BCD(t) = 0.0, [description = "BCD trimer", tunable = false]
        XBCD(t) = 0.0, [description = "X-BCD complex", tunable = false]
        XXBCD(t) = 0.0, [description = "XX-BCD complex", tunable = false]
        XXXBCD(t) = 0.0, [description = "XXX-BCD complex", tunable = false]
    end 

    #- Amem is the fraction of AP2 that is membrane-bound, AKA complexed in a 2D membrane-bound complex. The denominator simply represents the total amount of AP2 in the mass conserved system.
    @observables begin
        #* Old definition of Amem that doesn't include AKL and APLp, but is much better behaved during optimization (smoother fitness landscape and less instability). This is what is scored via the FFT based fitness function
        Amem_old ~ (LpA + LpAK + LpAP + LpAKL + LpAPLp + XB + XC + XD + XBC + 2XXBC + XBD + 2XXBD + XCD + 2XXCD + XBCD + 2XXBCD + 3XXXBCD) / (A + LpA + LpAK + LpAP + LpAKL + LpAPLp + AK + AP + AKL + APLp + XB + XC + XD + XBC + 2XXBC + XBD + 2XXBD + XCD + 2XXCD + XBCD + 2XXBCD + 3XXXBCD)
        #* New definition includes AKL and APLp, which should be included in order to match how any real experiemnt would measure it, even though A is only indirectly localized to the membrane in those species. Period and amplitude are measured from this definition, and all downstream analysis assumes this definition.
        Amem ~ (LpA + LpAK + LpAP + LpAKL + LpAPLp + AKL + APLp + XB + XC + XD + XBC + 2XXBC + XBD + 2XXBD + XCD + 2XXCD + XBCD + 2XXBCD + 3XXXBCD) / (A + LpA + LpAK + LpAP + LpAKL + LpAPLp + AK + AP + AKL + APLp + XB + XC + XD + XBC + 2XXBC + XBD + 2XXBD + XCD + 2XXCD + XBCD + 2XXBCD + 3XXXBCD)

        # Proportion of fully assembled trimer on the membrane (XXXBCD) versus total trimer (BCD + XXXBCD)
        Tmem ~ XXXBCD / (BCD + XBCD + XXBCD + XXXBCD)

        # Fraction of maximum possible trimers that are currently assembled
        TrimerYield ~ (BCD +XBCD + XXBCD + XXXBCD) / min(B + BC + BD + XB + XBC + XXBC + XBD + XXBD + BCD + XBCD + XXBCD + XXXBCD, C + BC + CD + XC + XBC + XXBC + XCD + XXCD + BCD + XBCD + XXBCD + XXXBCD, D + BD + CD + XD + XBD + XXBD + XCD + XXCD + BCD + XBCD + XXBCD + XXXBCD)
    end

    (kfᴸᴬ,krᴸᴬ), Lp + AK <--> LpAK
    (kfᴸᴬ*DF,krᴸᴬ), Lp + AKL <--> LpAKL

    (kfᴸᴬ,krᴸᴬ), Lp + AP <--> LpAP
    (kfᴸᴬ*DF,krᴸᴬ), Lp + APLp <--> LpAPLp

    (kfᴬᴷ,krᴬᴷ), A + K <--> AK
    (kfᴬᴾ,krᴬᴾ), A + P <--> AP

    (kfᴬᴷ,krᴬᴷ), A + LK <--> AKL
    (kfᴬᴾ,krᴬᴾ), A + LpP <--> APLp

    (kfᴬᴷ*DF,krᴬᴷ), LpA + LK <--> LpAKL
    (kfᴬᴾ*DF,krᴬᴾ), LpA + LpP <--> LpAPLp

    (kfᴸᴷ,krᴸᴷ), AK + L <--> AKL #binding of kinase to lipid
    kcatᴸᴷ, AKL --> Lp + AK #phosphorylation of lipid

    (kfᴸᴾ,krᴸᴾ), AP + Lp <--> APLp #binding of phosphatase to lipid
    kcatᴸᴾ, APLp --> L + AP #dephosphorylation of lipid


    #< Interface to oscilattor system via LpA and T (transferrin receptor)
    #- Lipid binding anchor referenced as X
    #- All monomers (B, C, D) have X binding interface in addition to inter-monomer binding interfaces
    #* The position of the components in the product symbol indicate binding interface, ex. XBC indicates B is membrane-bound and binding to C in the cytosol, so C is indirectly localized to the membrane
    # (kfᴬᵀ,krᴬᵀ), LpA + T <--> X

    #- Monomers localizing to the membrane (X stands for LpA)
    (kf¹ᵐ, kr¹ᵐ), LpA + B <--> XB 
    (kf¹ᵐ, kr¹ᵐ), LpA + C <--> XC 
    (kf¹ᵐ, kr¹ᵐ), LpA + D <--> XD 

    #- Non localized dimerization
    (kf²ᵐ, kr²ᵐ), B + C <--> BC
    (kf³ᵐ, kr³ᵐ), B + D <--> BD
    (kf⁴ᵐ, kr⁴ᵐ), C + D <--> CD

    (kf²ᵐ, kr²ᵐ), XB + C <--> XBC
    (kf²ᵐ * DF, kr²ᵐ), XB + XC <--> XXBC
    (kf²ᵐ, kr²ᵐ), XC + B <--> XBC

    (kf³ᵐ, kr³ᵐ), XB + D <--> XBD
    (kf³ᵐ * DF, kr³ᵐ), XB + XD <--> XXBD
    (kf³ᵐ, kr³ᵐ), XD + B <--> XBD

    (kf⁴ᵐ, kr⁴ᵐ), XC + D <--> XCD
    (kf⁴ᵐ * DF, kr⁴ᵐ), XC + XD <--> XXCD
    (kf⁴ᵐ, kr⁴ᵐ), XD + C <--> XCD 

    #- Non localized trimerization assembly, with unique loop-closure rate constants 
    (kf⁵ᵐ, kr⁵ᵐ), BC + D <--> BCD 
    (kf⁶ᵐ, kr⁶ᵐ), BD + C <--> BCD 
    (kf⁷ᵐ, kr⁷ᵐ), CD + B <--> BCD 

    (kf⁵ᵐ, kr⁵ᵐ), XBC + D <--> XBCD 
    (kf⁵ᵐ, kr⁵ᵐ), XXBC + D <--> XXBCD
    (kf⁵ᵐ * DF, kr⁵ᵐ), XXBC + XD <--> XXXBCD

    (kf⁶ᵐ, kr⁶ᵐ), XBD + C <--> XBCD 
    (kf⁶ᵐ, kr⁶ᵐ), XXBD + C <--> XXBCD 
    (kf⁶ᵐ * DF, kr⁶ᵐ), XXBD + XC <--> XXXBCD

    (kf⁷ᵐ, kr⁷ᵐ), XCD + B <--> XBCD 
    (kf⁷ᵐ, kr⁷ᵐ), XXCD + B <--> XXBCD 
    (kf⁷ᵐ * DF, kr⁷ᵐ), XXCD + XB <--> XXXBCD 

    #* All X-X binding has the same kinetics
    (kf¹ᵐ * DF, kr¹ᵐ), XBC + LpA <--> XXBC 
    (kf¹ᵐ * DF, kr¹ᵐ), XBD + LpA <--> XXBD 
    (kf¹ᵐ * DF, kr¹ᵐ), XCD + LpA <--> XXCD 
    (kf¹ᵐ * DF, kr¹ᵐ), XBCD + LpA <--> XXBCD 

    #- Full membrane-localized heterotrimer assembly
    (kf¹ᵐ * DF, kr¹ᵐ), XXBCD + LpA <--> XXXBCD 
end


function make_composed_trimer_rn(rn = fullrn)
    composed_rn = complete(extend(trimer_rn, rn))
    return composed_rn
end


function make_trimer_odeprob(rn = fullrn; tend::Float64=1968.2)
    composed_rn = make_composed_trimer_rn(rn)

    default_p = [:kfᴸᴬ => 32.970817677315374, :krᴸᴬ => 9.135573868071218, :kfᴸᴷ => 0.0017552768798002882, :krᴸᴷ => 228.4470892394084, :kcatᴸᴷ => 67.12498862496422, :kfᴸᴾ => 30.246525659038156, :krᴸᴾ => 0.07002848079965307, :kcatᴸᴾ => 7.001836143081583, :kfᴬᴷ => 1.476147454737371, :krᴬᴷ => 0.03370862149650637, :kfᴬᴾ => 25.513040115017223, :krᴬᴾ => 1.9862708921924483, :DF => 4073.95672132362, 
    # :kfᴬᵀ => 4.839846054560039, :krᴬᵀ => 1.936796566650092, 
    :kf¹ᵐ => 13.989273996308272, :kr¹ᵐ => 42.175587772519926, :kf²ᵐ => 22.621676004136553, :kr²ᵐ => 0.1588719693873535, :kf³ᵐ => 29.653410541978655, :kr³ᵐ => 3.592292741988041, :kf⁴ᵐ => 2.16031845959122, :kr⁴ᵐ => 4.224232546494592, :kf⁵ᵐ => 9.574779298553665, :kr⁵ᵐ => 0.5800610936134903, :kf⁶ᵐ => 0.4638662466223983, :kr⁶ᵐ => 3.4593697782834587, :kf⁷ᵐ => 25.307317392749262, :kr⁷ᵐ => 10.541433491743962]

    default_u0 = [:L => 33.283668743473605, :K => 36.12056864967913, :P => 12.003347878660461, :A => 3.6627442256750715, :B => 0.6807065430763064, :C => 3.8282579299469655, :D => 3.2617563329970793]#, :X => 17.620754062551956]

    odeprob = ODEProblem(composed_rn, default_u0, (0.0, tend), default_p; jac = true, saveat = 0.1, abstol = 1e-8, reltol = 1e-6)
    return odeprob
end




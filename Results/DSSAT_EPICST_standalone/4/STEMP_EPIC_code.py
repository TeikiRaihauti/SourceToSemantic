def SOILT_EPIC(B, BCV, CUMDPT, DP, DSMID, NLAYR, PESW, TAV, TAVG, TMAX, TMIN,
               WetDay, WFT, WW, TMA, SRFTEMP, ST, X2_AVG, X2_PREV):
    """
    EPIC soil temperature routine operating on the surface and profile.

    Parameters
    - B: Exponential decay factor (Parton and Logan)
    - BCV: Soil cover function (0-1)
    - CUMDPT: Cumulative soil profile depth (mm)
    - DP: Damping depth (mm)
    - DSMID: List of depths to midpoint of soil layers (mm)
    - NLAYR: Number of soil layers
    - PESW: Potential extractable soil water (cm)
    - TAV: Average annual soil temperature (degC)
    - TAVG: Daily average temperature (degC)
    - TMAX: Daily max temperature (degC)
    - TMIN: Daily min temperature (degC)
    - WetDay: Integer (0/1) indicating wet day (rain + irrigation)
    - WFT: Fraction of wet days in moving memory window
    - WW: Volumetric fraction related to bulk density
    - TMA: Length-5 list of prior "surface" temps memory (InOut)
    - SRFTEMP: Surface temperature (Output)
    - ST: Soil temperature by layer list (InOut)
    - X2_AVG: Average of TMA memory (Output)
    - X2_PREV: Previous surface temperature memory (InOut)

    Returns
    - TMA_out, SRFTEMP_out, ST_out, X2_AVG_out, X2_PREV_out
    """
    import math

    # Ensure copies for purity
    DSMID_loc = list(DSMID[:NLAYR])
    TMA_out = list(TMA)
    ST_out = list(ST[:NLAYR])

    # Water content term (ratio)
    # WC = PESW(cm) / [WW (dimensionless) * CUMDPT(mm)] * (10 mm/cm)
    denom = WW * CUMDPT if WW * CUMDPT != 0.0 else 1.0e-12
    WC = max(0.01, PESW) / denom * 10.0

    FX = math.exp(B * ((1.0 - WC) / (1.0 + WC)) ** 2)
    DD = FX * DP  # effective damping depth (mm)

    # Surface temperature candidate X2
    if WetDay > 0:
        # Potter & Williams (1994) Eqn 2
        X2 = WFT * (TAVG - TMIN) + TMIN
    else:
        # Potter & Williams (1994) Eqn 1 (ST0 term omitted as in DSSAT code)
        X2 = WFT * (TMAX - TAVG) + TAVG + 2.0

    # Update 5-day memory (note: ordering follows original Fortran code)
    TMA_out[0] = X2
    for K in range(4, 0, -1):
        TMA_out[K] = TMA_out[K - 1]
    X2_AVG_out = sum(TMA_out) / 5.0

    # Combine with cover effect (Eqn 4 and SRFTEMP assignment)
    X3 = (1.0 - BCV) * X2_AVG_out + BCV * X2_PREV
    SRFTEMP_out = min(X2_AVG_out, X3)

    # Eqn 6
    X1 = TAV - X3

    # Update soil layer temperatures
    LAG = 0.5
    for L in range(NLAYR):
        ZD = DSMID_loc[L] / DD if DD != 0.0 else 0.0
        F = ZD / (ZD + math.exp(-0.8669 - 2.0775 * ZD))
        ST_out[L] = LAG * ST_out[L] + (1.0 - LAG) * (F * X1 + X3)

    X2_PREV_out = X2_AVG_out

    return TMA_out, SRFTEMP_out, ST_out, X2_AVG_out, X2_PREV_out


def STEMP_EPIC(
    CONTROL_DYNAMIC,
    ISWITCH_ISWWAT,
    SOILPROP_BD,
    SOILPROP_DLAYR,
    SW,
    TAVG,
    TMAX,
    TMIN,
    TAV,
    WEATHER_RAIN,
    SOILPROP_DS,
    SOILPROP_DUL,
    SOILPROP_LL,
    SOILPROP_NLAYR,
    CUMDPT,
    DSMID,
    TDL,
    TMA,
    NDays,
    WetDay,
    X2_PREV,
    SRFTEMP,
    ST,
    CONTROL_RUN,
    CONTROL_RNMODE,
    ORGC_MULCHMASS,
    WATER_SNOW,
    MGMT_DEPIR,
    PLANT_BIOMAS,
):
    """
    EPIC soil temperature component: initialization and daily rate.

    Parameters (flattened from Fortran):
    - CONTROL_DYNAMIC: int, dynamic phase (use 2 for SEASINIT, 3 for RATE)
    - ISWITCH_ISWWAT: 'Y' if water simulation enabled, 'N' otherwise
    - SOILPROP_BD: list of bulk density by layer (g/cm3)
    - SOILPROP_DLAYR: list of layer thickness (cm)
    - SW: list of volumetric soil water content (cm3/cm3) by layer
    - TAVG, TMAX, TMIN: daily temperatures (degC)
    - TAV: average annual air temperature (degC)
    - WEATHER_RAIN: daily rainfall (mm)
    - SOILPROP_DS: list of cumulative depths (cm) to bottom of each layer
    - SOILPROP_DUL: list of drained upper limit (cm3/cm3) by layer
    - SOILPROP_LL: list of lower limit (cm3/cm3) by layer
    - SOILPROP_NLAYR: number of soil layers (int)
    - CUMDPT: cumulative profile depth (mm) (InOut)
    - DSMID: list depth to midpoint of soil layers (mm) (InOut)
    - TDL: total water content at DUL over profile (cm) (InOut)
    - TMA: length-5 memory array for surface temperature (InOut)
    - NDays: number of days stored in WetDay memory (InOut)
    - WetDay: length-30 list of 0/1 wet-day flags (InOut)
    - X2_PREV: previous step surface temp memory (InOut)
    - SRFTEMP: current surface temperature (InOut)
    - ST: list of soil temperatures by layer (InOut)
    - CONTROL_RUN: int run number (used during initialization condition)
    - CONTROL_RNMODE: 1-char run mode flag used in original (affects init)
    - ORGC_MULCHMASS: mulch mass (kg/ha)
    - WATER_SNOW: snow cover (mm)
    - MGMT_DEPIR: irrigation depth (mm)
    - PLANT_BIOMAS: plant biomass (kg/ha)

    Returns dict with updated state variables:
    {
      'CUMDPT': float,
      'DSMID': list,
      'TDL': float,
      'TMA': list,
      'NDays': int,
      'WetDay': list,
      'X2_PREV': float,
      'SRFTEMP': float,
      'ST': list
    }
    """
    import math

    # Local copies to preserve purity
    NLAYR = int(SOILPROP_NLAYR)
    BD = list(SOILPROP_BD[:NLAYR])
    DLAYR = list(SOILPROP_DLAYR[:NLAYR])
    DS = list(SOILPROP_DS[:NLAYR])
    DUL = list(SOILPROP_DUL[:NLAYR])
    LL = list(SOILPROP_LL[:NLAYR])
    SW_loc = list(SW[:NLAYR])

    CUMDPT_out = float(CUMDPT)
    DSMID_out = list(DSMID[:NLAYR]) if len(DSMID) >= NLAYR else [0.0] * NLAYR
    TDL_out = float(TDL)
    TMA_out = list(TMA[:5]) if len(TMA) >= 5 else [0.0] * 5
    NDays_out = int(NDays)
    WetDay_out = list(WetDay[:30]) if len(WetDay) >= 30 else [0] * 30
    X2_PREV_out = float(X2_PREV)
    SRFTEMP_out = float(SRFTEMP)
    ST_out = list(ST[:NLAYR]) if len(ST) >= NLAYR else [TAVG] * NLAYR

    def compute_profile_sums(use_swi=False):
        # Returns (TBD, TDL_calc, TLL, TSW) where:
        # TBD sum of BD*DLAYR; TDL total DUL*DLAYR; TLL total LL*DLAYR; TSW total SW*DLAYR
        TBD = 0.0
        TDL_calc = 0.0
        TLL_calc = 0.0
        TSW_calc = 0.0
        for L in range(NLAYR):
            TBD += BD[L] * DLAYR[L]
            TDL_calc += DUL[L] * DLAYR[L]
            TLL_calc += LL[L] * DLAYR[L]
            TSW_calc += SW_loc[L] * DLAYR[L]
        return TBD, TDL_calc, TLL_calc, TSW_calc

    def soil_cover_fraction(biomas, mulch, snow):
        # Compute BCV as max of residue/mulch and snow cover functions
        CV = (biomas + mulch) / 1000.0  # t/ha
        BCV1 = CV / (CV + math.exp(5.3396 - 2.3951 * CV)) if (CV + math.exp(5.3396 - 2.3951 * CV)) != 0.0 else 0.0
        BCV2 = snow / (snow + math.exp(2.303 - 0.2197 * snow)) if (snow + math.exp(2.303 - 0.2197 * snow)) != 0.0 else 0.0
        return max(BCV1, BCV2)

    if CONTROL_DYNAMIC == 2:  # SEASINIT
        # Initialization only under certain conditions in original code
        do_init = (CONTROL_RUN == 1) or (('Q' not in str(CONTROL_RNMODE)) and ('F' not in str(CONTROL_RNMODE)))
        if do_init:
            # Initialize soil depth metrics
            CUMDPT_out = 0.0
            TBD = 0.0
            TLL = 0.0
            TSW = 0.0
            TDL_out = 0.0
            DSMID_out = [0.0] * NLAYR
            for L in range(NLAYR):
                DSMID_out[L] = CUMDPT_out + DLAYR[L] * 5.0  # mm to midpoint
                CUMDPT_out += DLAYR[L] * 10.0  # mm
                TBD += BD[L] * DLAYR[L]
                TLL += LL[L] * DLAYR[L]
                TSW += SW_loc[L] * DLAYR[L]
                TDL_out += DUL[L] * DLAYR[L]

            # Potential extractable soil water (cm)
            if str(ISWITCH_ISWWAT).upper() == 'Y':
                PESW = max(0.0, TSW - TLL)
            else:
                PESW = max(0.0, TDL_out - TLL)

            # Bulk density and damping parameters
            denom_depth_cm = DS[NLAYR - 1] if NLAYR > 0 else 1.0e-12
            ABD = TBD / denom_depth_cm if denom_depth_cm != 0.0 else 0.0
            FX = ABD / (ABD + 686.0 * math.exp(-5.63 * ABD)) if (ABD + 686.0 * math.exp(-5.63 * ABD)) != 0.0 else 0.0
            DP = 1000.0 + 2500.0 * FX
            WW = 0.356 - 0.144 * ABD
            B = math.log(500.0 / DP) if DP != 0.0 else 0.0

            # Initialize TMA and ST with TAVG
            tavg_rounded = round(TAVG, 4)
            TMA_out = [tavg_rounded] * 5
            ST_out = [TAVG] * NLAYR

            # Wet-day memory initialization
            WFT = 0.1
            WetDay_out = [0] * 30
            NDays_out = 0

            # Soil cover function at initialization (mulch and/or snow)
            BCV = soil_cover_fraction(0.0, ORGC_MULCHMASS, WATER_SNOW)

            # Spin-up calls to stabilize TMA/ST
            X2_AVG_local = 0.0
            for _ in range(8):
                TMA_out, SRFTEMP_out, ST_out, X2_AVG_local, X2_PREV_out = SOILT_EPIC(
                    B=B,
                    BCV=BCV,
                    CUMDPT=CUMDPT_out,
                    DP=DP,
                    DSMID=DSMID_out,
                    NLAYR=NLAYR,
                    PESW=PESW,
                    TAV=TAV,
                    TAVG=TAVG,
                    TMAX=TMAX,
                    TMIN=TMIN,
                    WetDay=0,
                    WFT=WFT,
                    WW=WW,
                    TMA=TMA_out,
                    SRFTEMP=SRFTEMP_out,
                    ST=ST_out,
                    X2_AVG=0.0,
                    X2_PREV=X2_PREV_out,
                )

    elif CONTROL_DYNAMIC == 3:  # RATE
        # Profile sums (note: TDL_out accumulates as in original code)
        TBD = 0.0
        TLL = 0.0
        TSW = 0.0
        for L in range(NLAYR):
            TBD += BD[L] * DLAYR[L]
            TDL_out += DUL[L] * DLAYR[L]
            TLL += LL[L] * DLAYR[L]
            TSW += SW_loc[L] * DLAYR[L]

        # Bulk density and damping parameters
        denom_depth_cm = DS[NLAYR - 1] if NLAYR > 0 else 1.0e-12
        ABD = TBD / denom_depth_cm if denom_depth_cm != 0.0 else 0.0
        FX = ABD / (ABD + 686.0 * math.exp(-5.63 * ABD)) if (ABD + 686.0 * math.exp(-5.63 * ABD)) != 0.0 else 0.0
        DP = 1000.0 + 2500.0 * FX  # mm
        WW = 0.356 - 0.144 * ABD
        B = math.log(500.0 / DP) if DP != 0.0 else 0.0

        # PESW (cm)
        if str(ISWITCH_ISWWAT).upper() == 'Y':
            PESW = max(0.0, TSW - TLL)
        else:
            PESW = max(0.0, TDL_out - TLL)

        # Update wet-day 30-day memory
        if NDays_out == 30:
            # Shift left by one (Fortran 1..29 <- 2..30)
            WetDay_out = WetDay_out[1:] + [0]
        else:
            NDays_out += 1

        # Set today's wet-day flag
        if (WEATHER_RAIN + MGMT_DEPIR) > 1.0e-6:
            WetDay_out[NDays_out - 1] = 1
        else:
            WetDay_out[NDays_out - 1] = 0

        NWetDays = sum(WetDay_out)
        WFT = float(NWetDays) / float(NDays_out) if NDays_out > 0 else 0.0

        # Soil cover including plant biomass
        BCV = soil_cover_fraction(PLANT_BIOMAS, ORGC_MULCHMASS, WATER_SNOW)

        # Advance soil temperatures
        TMA_out, SRFTEMP_out, ST_out, X2_AVG_local, X2_PREV_out = SOILT_EPIC(
            B=B,
            BCV=BCV,
            CUMDPT=CUMDPT_out,
            DP=DP,
            DSMID=DSMID_out,
            NLAYR=NLAYR,
            PESW=PESW,
            TAV=TAV,
            TAVG=TAVG,
            TMAX=TMAX,
            TMIN=TMIN,
            WetDay=WetDay_out[NDays_out - 1],
            WFT=WFT,
            WW=WW,
            TMA=TMA_out,
            SRFTEMP=SRFTEMP_out,
            ST=ST_out,
            X2_AVG=0.0,
            X2_PREV=X2_PREV_out,
        )

    # For other DYNAMIC phases, no updates to states.

    return {
        'CUMDPT': CUMDPT_out,
        'DSMID': DSMID_out,
        'TDL': TDL_out,
        'TMA': TMA_out,
        'NDays': NDays_out,
        'WetDay': WetDay_out,
        'X2_PREV': X2_PREV_out,
        'SRFTEMP': SRFTEMP_out,
        'ST': ST_out,
    }


def test_STEMP_EPIC_basic():
    """
    Basic test derived from ASKEE_DSSAT_EPIC.for program.

    This test runs SEASINIT followed by a single RATE day with a simple
    4-layer soil and checks key state updates.
    """
    # Setup soil properties and initial conditions
    NLAYR = 4
    SOILPROP_BD = [1.6] * NLAYR
    SOILPROP_DLAYR = [10.0] * NLAYR  # cm
    SOILPROP_DS = [10.0, 20.0, 30.0, 40.0]  # cm cumulative
    SOILPROP_DUL = [0.3] * NLAYR
    SOILPROP_LL = [0.2] * NLAYR

    SW = [0.2] * NLAYR
    TAVG = 25.0
    TMAX = 30.0
    TMIN = 20.0
    TAV = 20.0

    # Exogenous
    WEATHER_RAIN = 0.0
    ORGC_MULCHMASS = 0.0
    WATER_SNOW = 0.0
    MGMT_DEPIR = 0.0
    PLANT_BIOMAS = 0.0

    # States
    CUMDPT = 0.0
    DSMID = [0.0] * NLAYR
    TDL = 0.0
    TMA = [0.0] * 5
    NDays = 0
    WetDay = [0] * 30
    X2_PREV = 0.0
    SRFTEMP = 0.0
    ST = [0.0] * NLAYR

    # Initialization (SEASINIT)
    out_init = STEMP_EPIC(
        CONTROL_DYNAMIC=2,
        ISWITCH_ISWWAT='Y',
        SOILPROP_BD=SOILPROP_BD,
        SOILPROP_DLAYR=SOILPROP_DLAYR,
        SW=SW,
        TAVG=TAVG,
        TMAX=TMAX,
        TMIN=TMIN,
        TAV=TAV,
        WEATHER_RAIN=WEATHER_RAIN,
        SOILPROP_DS=SOILPROP_DS,
        SOILPROP_DUL=SOILPROP_DUL,
        SOILPROP_LL=SOILPROP_LL,
        SOILPROP_NLAYR=NLAYR,
        CUMDPT=CUMDPT,
        DSMID=DSMID,
        TDL=TDL,
        TMA=TMA,
        NDays=NDays,
        WetDay=WetDay,
        X2_PREV=X2_PREV,
        SRFTEMP=SRFTEMP,
        ST=ST,
        CONTROL_RUN=1,
        CONTROL_RNMODE='B',
        ORGC_MULCHMASS=ORGC_MULCHMASS,
        WATER_SNOW=WATER_SNOW,
        MGMT_DEPIR=MGMT_DEPIR,
        PLANT_BIOMAS=PLANT_BIOMAS,
    )

    # Assertions after init
    assert len(out_init['DSMID']) == NLAYR
    assert out_init['CUMDPT'] == sum(SOILPROP_DLAYR) * 10.0
    assert len(out_init['TMA']) == 5
    assert len(out_init['ST']) == NLAYR
    assert out_init['NDays'] == 0

    # Daily rate (RATE)
    out_rate = STEMP_EPIC(
        CONTROL_DYNAMIC=3,
        ISWITCH_ISWWAT='Y',
        SOILPROP_BD=SOILPROP_BD,
        SOILPROP_DLAYR=SOILPROP_DLAYR,
        SW=SW,
        TAVG=TAVG,
        TMAX=TMAX,
        TMIN=TMIN,
        TAV=TAV,
        WEATHER_RAIN=WEATHER_RAIN,
        SOILPROP_DS=SOILPROP_DS,
        SOILPROP_DUL=SOILPROP_DUL,
        SOILPROP_LL=SOILPROP_LL,
        SOILPROP_NLAYR=NLAYR,
        CUMDPT=out_init['CUMDPT'],
        DSMID=out_init['DSMID'],
        TDL=out_init['TDL'],
        TMA=out_init['TMA'],
        NDays=out_init['NDays'],
        WetDay=out_init['WetDay'],
        X2_PREV=out_init['X2_PREV'],
        SRFTEMP=out_init['SRFTEMP'],
        ST=out_init['ST'],
        CONTROL_RUN=1,
        CONTROL_RNMODE='B',
        ORGC_MULCHMASS=ORGC_MULCHMASS,
        WATER_SNOW=WATER_SNOW,
        MGMT_DEPIR=MGMT_DEPIR,
        PLANT_BIOMAS=PLANT_BIOMAS,
    )

    # Assertions after one day
    assert out_rate['NDays'] == 1
    assert out_rate['WetDay'][0] in (0, 1)
    assert isinstance(out_rate['SRFTEMP'], float)
    assert len(out_rate['ST']) == NLAYR
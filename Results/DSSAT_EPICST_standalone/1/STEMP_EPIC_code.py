import math


def _nint(x):
    if x >= 0:
        return int(math.floor(x + 0.5))
    else:
        return int(math.ceil(x - 0.5))


def SOILT_EPIC(
    B,
    BCV,
    CUMDPT,
    DP,
    DSMID,
    NLAYR,
    PESW,
    TAV,
    TAVG,
    TMAX,
    TMIN,
    WetDay,
    WFT,
    WW,
    TMA,
    ST,
    X2_PREV,
):
    # Copy arrays to avoid mutating inputs
    DSMID_loc = list(DSMID)
    TMA_loc = list(TMA)
    ST_loc = list(ST)

    denom = WW * CUMDPT
    if denom <= 0.0:
        denom = 1.0e-6
    WC = max(0.01, PESW) / denom * 10.0

    FX = math.exp(B * ((1.0 - WC) / (1.0 + WC)) ** 2)
    DD = FX * DP

    if WetDay > 0:
        X2 = WFT * (TAVG - TMIN) + TMIN
    else:
        X2 = WFT * (TMAX - TAVG) + TAVG + 2.0

    # Note the original Fortran assigns TMA(1)=X2 and then shifts,
    # which results in duplicating X2 into TMA(2).
    TMA_loc[0] = X2
    for K in range(4, 0, -1):
        TMA_loc[K] = TMA_loc[K - 1]

    X2_AVG = sum(TMA_loc) / 5.0
    X3 = (1.0 - BCV) * X2_AVG + BCV * X2_PREV
    SRFTEMP = min(X2_AVG, X3)
    X1 = TAV - X3

    LAG = 0.5
    for L in range(NLAYR):
        ZD = DSMID_loc[L] / DD
        F = ZD / (ZD + math.exp(-0.8669 - 2.0775 * ZD))
        ST_loc[L] = LAG * ST_loc[L] + (1.0 - LAG) * (F * X1 + X3)

    X2_PREV_new = X2_AVG

    return TMA_loc, SRFTEMP, ST_loc, X2_PREV_new


def STEMP_EPIC(
    CONTROL_DYNAMIC,
    ISWITCH_ISWWAT,
    SOILPROP_BD,
    SOILPROP_DLAYR,
    SOILPROP_DS,
    SOILPROP_DUL,
    SOILPROP_LL,
    SOILPROP_NLAYR,
    SW,
    TAVG,
    TMAX,
    TMIN,
    TAV,
    WEATHER_RAIN,
    CUMDPT,
    DSMID,
    TDL,
    TMA,
    NDays,
    WetDay,
    X2_PREV,
    SRFTEMP,
    ST,
    ORGC_MULCHMASS,
    WATER_SNOW,
    MGMT_DEPIR,
    PLANT_BIOMAS,
):
    NLAYR = int(SOILPROP_NLAYR)

    BD = list(SOILPROP_BD)
    DLAYR = list(SOILPROP_DLAYR)
    DS = list(SOILPROP_DS)
    DUL = list(SOILPROP_DUL)
    LL = list(SOILPROP_LL)
    SW_arr = list(SW)
    DSMID_loc = list(DSMID)
    TMA_loc = list(TMA)
    ST_loc = list(ST)
    WetDay_loc = list(WetDay)

    if CONTROL_DYNAMIC == 2:
        # Seasonal initialization
        SWI = list(SW_arr)

        TBD = 0.0
        TLL = 0.0
        TSW = 0.0
        TDL_loc = 0.0
        CUMDPT_loc = 0.0
        for L in range(NLAYR):
            DSMID_loc[L] = CUMDPT_loc + DLAYR[L] * 5.0
            CUMDPT_loc = CUMDPT_loc + DLAYR[L] * 10.0
            TBD = TBD + BD[L] * DLAYR[L]
            TLL = TLL + LL[L] * DLAYR[L]
            TSW = TSW + SWI[L] * DLAYR[L]
            TDL_loc = TDL_loc + DUL[L] * DLAYR[L]

        if ISWITCH_ISWWAT == "Y":
            PESW = max(0.0, TSW - TLL)
        else:
            PESW = max(0.0, TDL_loc - TLL)

        ABD = TBD / DS[NLAYR - 1]
        FX = ABD / (ABD + 686.0 * math.exp(-5.63 * ABD))
        DP = 1000.0 + 2500.0 * FX
        WW = 0.356 - 0.144 * ABD
        B = math.log(500.0 / DP)

        for I in range(5):
            TMA_loc[I] = _nint(TAVG * 10000.0) / 10000.0
        X2_AVG = TMA_loc[0] * 5.0

        for L in range(NLAYR):
            ST_loc[L] = TAVG

        WFT = 0.1
        WetDay_loc = [0 for _ in range(len(WetDay_loc))]
        NDays_loc = 0

        CV = (ORGC_MULCHMASS) / 1000.0
        BCV1 = CV / (CV + math.exp(5.3396 - 2.3951 * CV)) if (CV + math.exp(5.3396 - 2.3951 * CV)) != 0 else 0.0
        SNOW = WATER_SNOW
        BCV2 = SNOW / (SNOW + math.exp(2.303 - 0.2197 * SNOW)) if (SNOW + math.exp(2.303 - 0.2197 * SNOW)) != 0 else 0.0
        BCV = max(BCV1, BCV2)

        for _ in range(8):
            TMA_loc, SRFTEMP, ST_loc, X2_PREV = SOILT_EPIC(
                B=B,
                BCV=BCV,
                CUMDPT=CUMDPT_loc,
                DP=DP,
                DSMID=DSMID_loc,
                NLAYR=NLAYR,
                PESW=PESW,
                TAV=TAV,
                TAVG=TAVG,
                TMAX=TMAX,
                TMIN=TMIN,
                WetDay=0,
                WFT=WFT,
                WW=WW,
                TMA=TMA_loc,
                ST=ST_loc,
                X2_PREV=X2_PREV,
            )

        return (
            CUMDPT_loc,
            DSMID_loc,
            TDL_loc,
            TMA_loc,
            NDays_loc,
            WetDay_loc,
            X2_PREV,
            SRFTEMP,
            ST_loc,
        )

    elif CONTROL_DYNAMIC == 3:
        # Daily rate calculations
        TBD = 0.0
        TLL = 0.0
        TSW = 0.0
        TDL_loc = TDL
        for L in range(NLAYR):
            TBD = TBD + BD[L] * DLAYR[L]
            TDL_loc = TDL_loc + DUL[L] * DLAYR[L]
            TLL = TLL + LL[L] * DLAYR[L]
            TSW = TSW + SW_arr[L] * DLAYR[L]

        ABD = TBD / DS[NLAYR - 1]
        FX = ABD / (ABD + 686.0 * math.exp(-5.63 * ABD))
        DP = 1000.0 + 2500.0 * FX
        WW = 0.356 - 0.144 * ABD
        B = math.log(500.0 / DP)

        if ISWITCH_ISWWAT == "Y":
            PESW = max(0.0, TSW - TLL)
        else:
            PESW = max(0.0, TDL_loc - TLL)

        RAIN = WEATHER_RAIN
        DEPIR = MGMT_DEPIR

        if int(NDays) == 30:
            # shift left (keep list length)
            WetDay_loc = WetDay_loc[1:] + [0]
            NDays_loc = 30
        else:
            NDays_loc = int(NDays) + 1
            if len(WetDay_loc) < 30:
                WetDay_loc = WetDay_loc + [0] * (30 - len(WetDay_loc))

        if RAIN + DEPIR > 1.0e-6:
            WetDay_loc[NDays_loc - 1] = 1
        else:
            WetDay_loc[NDays_loc - 1] = 0

        NWetDays = sum(WetDay_loc)
        WFT = float(NWetDays) / float(NDays_loc if NDays_loc > 0 else 1)

        BIOMAS = PLANT_BIOMAS
        MULCHMASS = ORGC_MULCHMASS
        SNOW = WATER_SNOW

        CV = (BIOMAS + MULCHMASS) / 1000.0
        BCV1 = CV / (CV + math.exp(5.3396 - 2.3951 * CV)) if (CV + math.exp(5.3396 - 2.3951 * CV)) != 0 else 0.0
        BCV2 = SNOW / (SNOW + math.exp(2.303 - 0.2197 * SNOW)) if (SNOW + math.exp(2.303 - 0.2197 * SNOW)) != 0 else 0.0
        BCV = max(BCV1, BCV2)

        TMA_loc, SRFTEMP, ST_loc, X2_PREV = SOILT_EPIC(
            B=B,
            BCV=BCV,
            CUMDPT=CUMDPT,
            DP=DP,
            DSMID=DSMID_loc,
            NLAYR=NLAYR,
            PESW=PESW,
            TAV=TAV,
            TAVG=TAVG,
            TMAX=TMAX,
            TMIN=TMIN,
            WetDay=WetDay_loc[NDays_loc - 1],
            WFT=WFT,
            WW=WW,
            TMA=TMA_loc,
            ST=ST_loc,
            X2_PREV=X2_PREV,
        )

        return (
            CUMDPT,
            DSMID_loc,
            TDL_loc,
            TMA_loc,
            NDays_loc,
            WetDay_loc,
            X2_PREV,
            SRFTEMP,
            ST_loc,
        )

    else:
        # For phases not performing computations here, return inputs unchanged
        return (
            CUMDPT,
            DSMID_loc,
            TDL,
            TMA_loc,
            NDays,
            WetDay_loc,
            X2_PREV,
            SRFTEMP,
            ST_loc,
        )
def nint(x):
    if x >= 0:
        from math import floor
        return floor(x + 0.5)
    else:
        from math import ceil
        return ceil(x - 0.5)


def YR_DOY(YRDOY):
    YR = int(YRDOY // 1000)
    DOY = int(YRDOY - YR * 1000)
    return YR, DOY


def SOILT(ALBEDO, B, CUMDPT, DOY, DP, HDAY, NLAYR,
          PESW, SRAD, TAMP, TAV, TAVG, TMAX, WW, DSMID,
          ATOT, TMA, SRFTEMP, ST):
    import math

    # Copy and prepare state arrays to avoid mutating inputs
    TMA_new = [float(v) for v in TMA]
    ST_new = [float(v) for v in ST]
    DSMID_use = [float(v) for v in DSMID]
    ATOT_new = float(ATOT)

    ALX = (float(DOY) - float(HDAY)) * 0.0174
    ATOT_new = ATOT_new - TMA_new[4]

    # Shift TMA
    TMA_new[4] = TMA_new[3]
    TMA_new[3] = TMA_new[2]
    TMA_new[2] = TMA_new[1]
    TMA_new[1] = TMA_new[0]
    TMA_new[0] = float(TAVG)

    # Keep only 4 decimals
    TMA_new[0] = nint(TMA_new[0] * 10000.0) / 10000.0
    ATOT_new = ATOT_new + TMA_new[0]

    # Corrected water content equation
    # WC is ratio; PESW (cm), WW (dimensionless), CUMDPT (mm)
    if WW * CUMDPT != 0.0:
        WC = max(0.01, float(PESW)) / (float(WW) * float(CUMDPT)) * 10.0
    else:
        WC = 0.01

    FX = math.exp(float(B) * ((1.0 - WC) / (1.0 + WC)) ** 2)
    DD = FX * float(DP)  # mm

    TA = float(TAV) + float(TAMP) * math.cos(ALX) / 2.0
    DT = ATOT_new / 5.0 - TA

    for L in range(max(0, int(NLAYR))):
        ZD = -float(DSMID_use[L]) / DD if DD != 0.0 else 0.0
        ST_val = float(TAV) + (float(TAMP) / 2.0 * math.cos(ALX + ZD) + DT) * math.exp(ZD)
        ST_new[L] = nint(ST_val * 1000.0) / 1000.0

    SRFTEMP_new = float(TAV) + (float(TAMP) / 2.0 * math.cos(ALX) + DT)

    return ATOT_new, TMA_new, SRFTEMP_new, ST_new


def STEMP(control_DYNAMIC, control_YRDOY, control_RUN, control_RNMODE,
          iswitch_ISWWAT,
          soilprop_BD, soilprop_DLAYR, soilprop_DS, soilprop_DUL, soilprop_LL,
          soilprop_NLAYR, soilprop_MSALB,
          SRAD, SW, TAVG, TMAX, XLAT, TAV, TAMP,
          CUMDPT, DSMID, TDL, TMA, ATOT, SRFTEMP, ST):
    import math

    # Local dynamic state enumerations (per ModuleDefs)
    RUNINIT = 1
    SEASINIT = 2
    RATE = 3
    OUTPUT = 5
    SEASEND = 6

    # Copy inputs to avoid in-place mutation; treat arrays as lists
    NLAYR = int(soilprop_NLAYR)
    BD = [float(v) for v in soilprop_BD[:NLAYR]]
    DLAYR = [float(v) for v in soilprop_DLAYR[:NLAYR]]
    DS = [float(v) for v in soilprop_DS[:NLAYR]]
    DUL = [float(v) for v in soilprop_DUL[:NLAYR]]
    LL = [float(v) for v in soilprop_LL[:NLAYR]]
    SW_use = [float(v) for v in SW[:NLAYR]]
    DSMID_new = [float(v) for v in DSMID[:NLAYR]]
    ST_new = [float(v) for v in ST[:NLAYR]]
    TMA_new = [float(v) for v in TMA[:5]]
    ATOT_new = float(ATOT)
    CUMDPT_new = float(CUMDPT)
    TDL_new = float(TDL)
    SRFTEMP_new = float(SRFTEMP)
    MSALB = float(soilprop_MSALB)
    XLAT_use = float(XLAT)
    TAVG_use = float(TAVG)
    TMAX_use = float(TMAX)
    TAV_use = float(TAV)
    TAMP_use = float(TAMP)
    SRAD_use = float(SRAD)

    # Compute year/day-of-year
    _, DOY = YR_DOY(int(control_YRDOY))

    # Compute HDAY based on latitude each call (pure; matches saved behavior)
    if XLAT_use < 0.0:
        HDAY = 20.0
    else:
        HDAY = 200.0

    if int(control_DYNAMIC) == SEASINIT:
        # Gate: initialize except in sequenced runs for RNMODE Q or F when RUN > 1
        do_init = (int(control_RUN) == 1) or (str(control_RNMODE).strip() not in ('Q', 'F'))
        if do_init:
            # Initial DS (cumulative depths) and SWI from current SW
            DSI = [float(v) for v in DS]
            SWI = [float(v) for v in SW_use]

            # Initialize profile geometric properties
            TBD = 0.0
            TLL = 0.0
            TSW = 0.0
            TDL_acc = 0.0
            CUMDPT_new = 0.0
            DLI = [0.0] * NLAYR  # individual layer thickness (cm)

            for L in range(NLAYR):
                if L == 0:
                    DLI[L] = DSI[L]
                else:
                    DLI[L] = DSI[L] - DSI[L - 1]
                DSMID_new[L] = CUMDPT_new + DLI[L] * 5.0  # mm to midpoint
                CUMDPT_new = CUMDPT_new + DLI[L] * 10.0   # mm cumulative depth
                TBD += BD[L] * DLI[L]
                TLL += LL[L] * DLI[L]
                TSW += SWI[L] * DLI[L]
                TDL_acc += DUL[L] * DLI[L]

            if str(iswitch_ISWWAT).strip().upper() == 'Y':
                PESW = max(0.0, TSW - TLL)
            else:
                PESW = max(0.0, TDL_acc - TLL)

            ABD = TBD / DSI[NLAYR - 1] if DSI[NLAYR - 1] != 0.0 else 0.0
            FX = ABD / (ABD + 686.0 * math.exp(-5.63 * ABD)) if (ABD + 686.0 * math.exp(-5.63 * ABD)) != 0.0 else 0.0
            DP = 1000.0 + 2500.0 * FX
            WW = 0.356 - 0.144 * ABD
            B = math.log(500.0 / DP) if DP != 0.0 else 0.0
            ALBEDO = MSALB

            # Initialize TMA and ATOT
            for i in range(5):
                TMA_new[i] = nint(TAVG_use * 10000.0) / 10000.0
            ATOT_new = TMA_new[0] * 5.0

            for L in range(NLAYR):
                ST_new[L] = TAVG_use

            # Spin-up calls to stabilize thermal history
            for _ in range(8):
                ATOT_new, TMA_new, SRFTEMP_new, ST_new = SOILT(
                    ALBEDO, B, CUMDPT_new, DOY, DP, HDAY, NLAYR,
                    PESW, SRAD_use, TAMP_use, TAV_use, TAVG_use, TMAX_use, WW, DSMID_new,
                    ATOT_new, TMA_new, SRFTEMP_new, ST_new
                )

            # Update TDL state to profile DUL (cm)
            TDL_new = TDL_acc

        # No output actions (I/O excluded)

    elif int(control_DYNAMIC) == RATE:
        # Daily rate calculations
        TBD = 0.0
        TLL = 0.0
        TSW = 0.0
        TDL_acc = 0.0

        for L in range(NLAYR):
            TBD += BD[L] * DLAYR[L]
            TDL_acc += DUL[L] * DLAYR[L]
            TLL += LL[L] * DLAYR[L]
            TSW += SW_use[L] * DLAYR[L]

        # Average bulk density across profile thickness
        denom = DS[NLAYR - 1] if NLAYR > 0 else 0.0
        ABD = TBD / denom if denom != 0.0 else 0.0
        FX = ABD / (ABD + 686.0 * math.exp(-5.63 * ABD)) if (ABD + 686.0 * math.exp(-5.63 * ABD)) != 0.0 else 0.0
        DP = 1000.0 + 2500.0 * FX
        WW = 0.356 - 0.144 * ABD
        B = math.log(500.0 / DP) if DP != 0.0 else 0.0
        ALBEDO = MSALB

        if str(iswitch_ISWWAT).strip().upper() == 'Y':
            PESW = max(0.0, TSW - TLL)
        else:
            PESW = max(0.0, TDL_acc - TLL)

        ATOT_new, TMA_new, SRFTEMP_new, ST_new = SOILT(
            ALBEDO, B, CUMDPT_new, DOY, DP, HDAY, NLAYR,
            PESW, SRAD_use, TAMP_use, TAV_use, TAVG_use, TMAX_use, WW, DSMID_new,
            ATOT_new, TMA_new, SRFTEMP_new, ST_new
        )

        # Update TDL state to current profile DUL (cm)
        TDL_new = TDL_acc

    elif int(control_DYNAMIC) == OUTPUT or int(control_DYNAMIC) == SEASEND:
        # No I/O in this refactoring; states unchanged
        pass

    return CUMDPT_new, DSMID_new, TDL_new, TMA_new, ATOT_new, SRFTEMP_new, ST_new
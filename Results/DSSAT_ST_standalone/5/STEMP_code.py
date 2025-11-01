def fortran_nint(x):
    # Fortran NINT: nearest integer, with ties away from zero
    import math
    if x >= 0:
        return int(math.floor(x + 0.5))
    else:
        return int(math.ceil(x - 0.5))


def SOILT(ALBEDO, B, CUMDPT, DOY, DP, HDAY, NLAYR,
          PESW, SRAD, TAMP, TAV, TAVG, TMAX, WW, DSMID,
          ATOT, TMA):
    # Inputs:
    # - ALBEDO, SRAD, TMAX are not used in computations but retained for signature fidelity.
    # InOut:
    # - ATOT: sum of TMA over last 5 days
    # - TMA: 5-element list of previous daily average temperatures
    # Outputs:
    # - SRFTEMP: surface litter temperature
    # - ST: list of soil temperatures per layer
    import math

    ALX = (float(DOY) - HDAY) * 0.0174
    ATOT = ATOT - TMA[4]

    # Shift temperature memory (most recent becomes index 0)
    for k in range(4, 0, -1):
        TMA[k] = TMA[k - 1]
    TMA[0] = TAVG

    # Prevents differences between release & debug modes:
    # Keep only 4 decimals (Fortran NINT behavior)
    TMA[0] = fortran_nint(TMA[0] * 10000.0) / 10000.0
    ATOT = ATOT + TMA[0]

    # Corrected water content function (EPIC)
    WC = max(0.01, PESW) / (WW * CUMDPT) * 10.0

    FX = math.exp(B * ((1.0 - WC) / (1.0 + WC))**2)
    DD = FX * DP  # damping depth in mm

    TA = TAV + TAMP * math.cos(ALX) / 2.0
    DT = ATOT / 5.0 - TA

    ST = [0.0] * NLAYR
    for l in range(NLAYR):
        ZD = -DSMID[l] / DD
        ST_l = TAV + (TAMP / 2.0 * math.cos(ALX + ZD) + DT) * math.exp(ZD)
        ST[l] = fortran_nint(ST_l * 1000.0) / 1000.0

    SRFTEMP = TAV + (TAMP / 2.0 * math.cos(ALX) + DT)

    return ATOT, TMA, SRFTEMP, ST


def STEMP_Initialize(ISWITCH_ISWWAT,
                     SOILPROP_BD, SOILPROP_DS, SOILPROP_DUL, SOILPROP_LL,
                     SOILPROP_NLAYR, SOILPROP_MSALB,
                     SRAD, SW, TAVG, TMAX, XLAT, TAV, TAMP, DOY):
    # Initialization logic corresponding to SEASINIT in Fortran STEMP
    # Inputs:
    # - ISWITCH_ISWWAT: 'Y' or 'N'
    # - SOILPROP_*: soil property arrays and scalars
    # - SRAD, SW (array), TAVG, TMAX, XLAT, TAV, TAMP, DOY
    # Returns:
    # - CUMDPT, DSMID, TDL, TMA, ATOT, SRFTEMP, ST, HDAY
    import math

    NLAYR = SOILPROP_NLAYR
    BD = list(SOILPROP_BD[:NLAYR])
    DSI = list(SOILPROP_DS[:NLAYR])
    DUL = list(SOILPROP_DUL[:NLAYR])
    LL = list(SOILPROP_LL[:NLAYR])
    MSALB = SOILPROP_MSALB
    SWI = list(SW[:NLAYR])

    # Hottest day of year depending on hemisphere
    if XLAT < 0.0:
        HDAY = 20.0   # Southern hemisphere
    else:
        HDAY = 200.0  # Northern hemisphere

    TBD = 0.0
    TLL = 0.0
    TSW = 0.0
    TDL = 0.0
    CUMDPT = 0.0

    DSMID = [0.0] * NLAYR
    DLI = [0.0] * NLAYR
    for l in range(NLAYR):
        if l == 0:
            DLI[l] = DSI[l]
        else:
            DLI[l] = DSI[l] - DSI[l - 1]
        DSMID[l] = CUMDPT + DLI[l] * 5.0     # mm to midpoint of layer
        CUMDPT = CUMDPT + DLI[l] * 10.0      # mm profile depth
        TBD = TBD + BD[l] * DLI[l]
        TLL = TLL + LL[l] * DLI[l]
        TSW = TSW + SWI[l] * DLI[l]
        TDL = TDL + DUL[l] * DLI[l]

    if ISWITCH_ISWWAT == 'Y':
        PESW = max(0.0, TSW - TLL)  # cm
    else:
        PESW = max(0.0, TDL - TLL)  # cm

    ABD = TBD / DSI[NLAYR - 1]
    FX = ABD / (ABD + 686.0 * math.exp(-5.63 * ABD))
    DP = 1000.0 + 2500.0 * FX
    WW = 0.356 - 0.144 * ABD
    B = math.log(500.0 / DP)
    ALBEDO = MSALB

    # Initialize TMA with TAVG and ATOT as sum over 5 days
    TMA = [fortran_nint(TAVG * 10000.0) / 10000.0 for _ in range(5)]
    ATOT = TMA[0] * 5.0

    ST = [TAVG for _ in range(NLAYR)]

    # Spin-up: call SOILT 8 times
    SRFTEMP = TAVG
    for _ in range(8):
        ATOT, TMA, SRFTEMP, ST = SOILT(ALBEDO, B, CUMDPT, DOY, DP, HDAY, NLAYR,
                                       PESW, SRAD, TAMP, TAV, TAVG, TMAX, WW, DSMID,
                                       ATOT, TMA)

    return CUMDPT, DSMID, TDL, TMA, ATOT, SRFTEMP, ST, HDAY


def STEMP(ISWITCH_ISWWAT,
          SOILPROP_BD, SOILPROP_DLAYR, SOILPROP_DS, SOILPROP_DUL, SOILPROP_LL,
          SOILPROP_NLAYR, SOILPROP_MSALB,
          SRAD, SW, TAVG, TMAX, TAV, TAMP, DOY, HDAY,
          CUMDPT, DSMID, TDL, TMA, ATOT):
    # Main daily biophysical process corresponding to RATE in Fortran STEMP
    # Inputs:
    # - ISWITCH_ISWWAT: 'Y' or 'N'
    # - SOILPROP_*: soil property arrays and scalars
    # - SRAD, SW (array), TAVG, TMAX, TAV, TAMP, DOY, HDAY
    # - CUMDPT, DSMID, TDL, TMA, ATOT: state variables from initialization or previous day
    # Returns:
    # - CUMDPT, DSMID, TDL, TMA, ATOT, SRFTEMP, ST
    import math

    NLAYR = SOILPROP_NLAYR
    BD = list(SOILPROP_BD[:NLAYR])
    DLAYR = list(SOILPROP_DLAYR[:NLAYR])
    DS = list(SOILPROP_DS[:NLAYR])
    DUL = list(SOILPROP_DUL[:NLAYR])
    LL = list(SOILPROP_LL[:NLAYR])
    MSALB = SOILPROP_MSALB
    SW_arr = list(SW[:NLAYR])

    TBD = 0.0
    TLL = 0.0
    TSW = 0.0
    for l in range(NLAYR):
        TBD = TBD + BD[l] * DLAYR[l]
        TDL = TDL + DUL[l] * DLAYR[l]  # Note: mirrors original code behavior
        TLL = TLL + LL[l] * DLAYR[l]
        TSW = TSW + SW_arr[l] * DLAYR[l]

    ABD = TBD / DS[NLAYR - 1]
    FX = ABD / (ABD + 686.0 * math.exp(-5.63 * ABD))
    DP = 1000.0 + 2500.0 * FX   # mm
    WW = 0.356 - 0.144 * ABD
    B = math.log(500.0 / DP)
    ALBEDO = MSALB

    if ISWITCH_ISWWAT == 'Y':
        PESW = max(0.0, TSW - TLL)  # cm
    else:
        PESW = max(0.0, TDL - TLL)  # cm

    ATOT, TMA, SRFTEMP, ST = SOILT(ALBEDO, B, CUMDPT, DOY, DP, HDAY, NLAYR,
                                   PESW, SRAD, TAMP, TAV, TAVG, TMAX, WW, DSMID,
                                   ATOT, TMA)

    return CUMDPT, DSMID, TDL, TMA, ATOT, SRFTEMP, ST


def test_STEMP_basic():
    # Simple test derived from ASKEE.for usage pattern and defaults
    # This test does not assert values; it exercises the initialize and daily rate functions.
    ISWITCH_ISWWAT = 'Y'
    NLAYR = 4
    SOILPROP_BD = [1.6] * NLAYR
    SOILPROP_DLAYR = [10.0, 10.0, 10.0, 10.0]
    SOILPROP_DS = [10.0, 20.0, 30.0, 40.0]
    SOILPROP_DUL = [0.3] * NLAYR
    SOILPROP_LL = [0.2] * NLAYR
    SOILPROP_NLAYR = NLAYR
    SOILPROP_MSALB = 0.13

    SRAD = 20.0
    SW = [0.2] * NLAYR
    TAVG = 25.0
    TMAX = 30.0
    XLAT = 28.0
    TAV = 20.0
    TAMP = 10.0
    DOY = 100

    CUMDPT, DSMID, TDL, TMA, ATOT, SRFTEMP, ST, HDAY = STEMP_Initialize(
        ISWITCH_ISWWAT,
        SOILPROP_BD, SOILPROP_DS, SOILPROP_DUL, SOILPROP_LL,
        SOILPROP_NLAYR, SOILPROP_MSALB,
        SRAD, SW, TAVG, TMAX, XLAT, TAV, TAMP, DOY
    )

    # Run one daily rate step
    CUMDPT, DSMID, TDL, TMA, ATOT, SRFTEMP, ST = STEMP(
        ISWITCH_ISWWAT,
        SOILPROP_BD, SOILPROP_DLAYR, SOILPROP_DS, SOILPROP_DUL, SOILPROP_LL,
        SOILPROP_NLAYR, SOILPROP_MSALB,
        SRAD, SW, TAVG, TMAX, TAV, TAMP, DOY, HDAY,
        CUMDPT, DSMID, TDL, TMA, ATOT
    )

    # Basic sanity checks (non-assertive; just ensure shapes and types)
    assert isinstance(SRFTEMP, float)
    assert len(ST) == NLAYR
    assert len(DSMID) == NLAYR
    assert len(TMA) == 5
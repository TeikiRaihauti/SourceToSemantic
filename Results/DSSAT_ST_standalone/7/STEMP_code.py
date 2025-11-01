def fortran_nint(x: float) -> int:
    """
    Fortran NINT-equivalent: nearest integer with ties away from zero.
    """
    import math
    if x >= 0.0:
        return int(math.floor(x + 0.5))
    else:
        return -int(math.floor(abs(x) + 0.5))


def YR_DOY(YRDOY: int) -> (int, int):
    """
    Convert YRDOY (YYYYDDD) to (YEAR, DOY).

    Inputs:
    - YRDOY: int

    Returns:
    - YEAR: int
    - DOY: int
    """
    YEAR = YRDOY // 1000
    DOY = YRDOY - YEAR * 1000
    return YEAR, DOY


def SOILT(
    ALBEDO: float,
    B: float,
    CUMDPT: float,
    DOY: int,
    DP: float,
    HDAY: float,
    NLAYR: int,
    PESW: float,
    SRAD: float,
    TAMP: float,
    TAV: float,
    TAVG: float,
    TMAX: float,
    WW: float,
    DSMID: list,
    ATOT: float,
    TMA: list,
):
    """
    Determines soil temperature by layer.

    Inputs:
    - ALBEDO: float
    - B: float
    - CUMDPT: float
    - DOY: int
    - DP: float
    - HDAY: float
    - NLAYR: int
    - PESW: float
    - SRAD: float
    - TAMP: float
    - TAV: float
    - TAVG: float
    - TMAX: float
    - WW: float
    - DSMID: list[float] length >= NLAYR
    - ATOT: float (in/out input)
    - TMA: list[float] length 5 (in/out input)

    Returns:
    - ATOT_out: float
    - TMA_out: list[float] length 5
    - SRFTEMP: float
    - ST: list[float] length NLAYR
    """
    import math

    # Copy to avoid mutating caller lists
    TMA_out = list(TMA)
    ATOT_out = ATOT

    ALX = (float(DOY) - HDAY) * 0.0174

    # Update running 5-day average array and sum
    ATOT_out = ATOT_out - TMA_out[4]
    for k in range(4, 0, -1):
        TMA_out[k] = TMA_out[k - 1]
    TMA_out[0] = TAVG
    TMA_out[0] = fortran_nint(TMA_out[0] * 10000.0) / 10000.0
    ATOT_out = ATOT_out + TMA_out[0]

    # Water content function (corrected)
    denom = WW * CUMDPT
    if denom <= 0.0:
        WC = 0.01
    else:
        WC = max(0.01, PESW) / denom * 10.0

    FX = math.exp(B * ((1.0 - WC) / (1.0 + WC)) ** 2)
    DD = FX * DP

    TA = TAV + TAMP * math.cos(ALX) / 2.0
    DT = ATOT_out / 5.0 - TA

    ST = [0.0] * NLAYR
    for L in range(NLAYR):
        ZD = -DSMID[L] / DD if DD != 0.0 else 0.0
        ST_L = TAV + (TAMP / 2.0 * math.cos(ALX + ZD) + DT) * math.exp(ZD)
        ST[L] = fortran_nint(ST_L * 1000.0) / 1000.0

    SRFTEMP = TAV + (TAMP / 2.0 * math.cos(ALX) + DT)

    return ATOT_out, TMA_out, SRFTEMP, ST


def STEMP(
    control_dynamic: int,
    control_yrdoy: int,
    control_run: int,
    control_rnmode: str,
    iswitch_iswwat: str,
    soilprop_bd: list,
    soilprop_dlayr: list,
    soilprop_ds: list,
    soilprop_dul: list,
    soilprop_ll: list,
    soilprop_nlayr: int,
    soilprop_msalb: float,
    SRAD: float,
    SW: list,
    TAVG: float,
    TMAX: float,
    XLAT: float,
    TAV: float,
    TAMP: float,
    CUMDPT: float,
    DSMID: list,
    TDL: float,
    TMA: list,
    ATOT: float,
    SRFTEMP: float,
    ST: list,
):
    """
    Determines soil temperature by layer.

    Inputs:
    - control_dynamic: int (1..7; use 2 for SEASINIT, 3 for RATE)
    - control_yrdoy: int (YYYYDDD)
    - control_run: int
    - control_rnmode: str (single-char; e.g., 'B','Q','F',...)
    - iswitch_iswwat: str ('Y' or 'N')
    - soilprop_bd: list[float] length >= NLAYR
    - soilprop_dlayr: list[float] length >= NLAYR
    - soilprop_ds: list[float] length >= NLAYR
    - soilprop_dul: list[float] length >= NLAYR
    - soilprop_ll: list[float] length >= NLAYR
    - soilprop_nlayr: int
    - soilprop_msalb: float
    - SRAD: float
    - SW: list[float] length >= NLAYR
    - TAVG: float
    - TMAX: float
    - XLAT: float
    - TAV: float
    - TAMP: float
    - CUMDPT: float (in/out)
    - DSMID: list[float] length >= NLAYR (in/out)
    - TDL: float (in/out)
    - TMA: list[float] length 5 (in/out)
    - ATOT: float (in/out)
    - SRFTEMP: float (in/out)
    - ST: list[float] length >= NLAYR (in/out)

    Returns (updated):
    - CUMDPT: float
    - DSMID: list[float] length NLAYR
    - TDL: float
    - TMA: list[float] length 5
    - ATOT: float
    - SRFTEMP: float
    - ST: list[float] length NLAYR
    """
    import math

    # Copy mutable inputs to avoid side-effects
    NLAYR = int(soilprop_nlayr)
    BD = list(soilprop_bd[:NLAYR])
    DLAYR = list(soilprop_dlayr[:NLAYR])
    DS = list(soilprop_ds[:NLAYR])
    DUL = list(soilprop_dul[:NLAYR])
    LL = list(soilprop_ll[:NLAYR])
    MSALB = float(soilprop_msalb)
    SW_in = list(SW[:NLAYR])

    # Ensure in/out arrays sized properly
    DSMID_out = list(DSMID[:NLAYR]) if DSMID else [0.0] * NLAYR
    ST_out = list(ST[:NLAYR]) if ST else [0.0] * NLAYR
    TMA_out = list(TMA[:5]) if TMA else [0.0] * 5

    # Compute DOY
    _, DOY = YR_DOY(control_yrdoy)

    # Seasonal initialization
    if control_dynamic == 2:  # SEASINIT
        # Only initialize if first run or RNMODE not in ('Q','F')
        if (control_run == 1) or (control_rnmode not in ('Q', 'F')):
            # Initial conditions
            SWI = SW_in
            DSI = DS

            HDAY = 20.0 if XLAT < 0.0 else 200.0

            TBD = 0.0
            TLL = 0.0
            TSW_sum = 0.0
            TDL_out = 0.0
            CUMDPT_out = 0.0
            DLI = [0.0] * NLAYR

            for L in range(NLAYR):
                if L == 0:
                    DLI[L] = DSI[L]
                else:
                    DLI[L] = DSI[L] - DSI[L - 1]
                DSMID_out[L] = CUMDPT_out + DLI[L] * 5.0
                CUMDPT_out = CUMDPT_out + DLI[L] * 10.0
                TBD = TBD + BD[L] * DLI[L]
                TLL = TLL + LL[L] * DLI[L]
                TSW_sum = TSW_sum + SWI[L] * DLI[L]
                TDL_out = TDL_out + DUL[L] * DLI[L]

            if iswitch_iswwat == 'Y':
                PESW = max(0.0, TSW_sum - TLL)
            else:
                PESW = max(0.0, TDL_out - TLL)

            ABD = TBD / DSI[NLAYR - 1] if DSI[NLAYR - 1] != 0.0 else 0.0
            FX = ABD / (ABD + 686.0 * math.exp(-5.63 * ABD)) if ABD + 686.0 * math.exp(-5.63 * ABD) != 0.0 else 0.0
            DP = 1000.0 + 2500.0 * FX
            WW = 0.356 - 0.144 * ABD
            B = math.log(500.0 / DP) if DP != 0.0 else 0.0
            ALBEDO = MSALB

            # Initialize moving average and layer temperatures
            TMA_out = [fortran_nint(TAVG * 10000.0) / 10000.0 for _ in range(5)]
            ATOT_out = TMA_out[0] * 5.0
            ST_out = [TAVG for _ in range(NLAYR)]

            # Spin-up calls
            for _ in range(8):
                ATOT_out, TMA_out, SRFTEMP_out, ST_out = SOILT(
                    ALBEDO, B, CUMDPT_out, DOY, DP, HDAY, NLAYR,
                    PESW, SRAD, TAMP, TAV, TAVG, TMAX, WW, DSMID_out,
                    ATOT_out, TMA_out
                )

            # Assign outputs
            CUMDPT = CUMDPT_out
            DSMID = DSMID_out
            TDL = TDL_out
            TMA = TMA_out
            ATOT = ATOT_out
            SRFTEMP = SRFTEMP_out
            ST = ST_out

        # No else branch actions (output-only in original)

    # Daily rate calculations
    elif control_dynamic == 3:  # RATE
        TBD = 0.0
        TLL = 0.0
        TSW_sum = 0.0
        TDL_out = TDL  # Note: original code accumulates into TDL without reset

        for L in range(NLAYR):
            TBD = TBD + BD[L] * DLAYR[L]
            TDL_out = TDL_out + DUL[L] * DLAYR[L]
            TLL = TLL + LL[L] * DLAYR[L]
            TSW_sum = TSW_sum + SW_in[L] * DLAYR[L]

        ABD = TBD / DS[NLAYR - 1] if DS[NLAYR - 1] != 0.0 else 0.0
        FX = ABD / (ABD + 686.0 * math.exp(-5.63 * ABD)) if ABD + 686.0 * math.exp(-5.63 * ABD) != 0.0 else 0.0
        DP = 1000.0 + 2500.0 * FX
        WW = 0.356 - 0.144 * ABD
        B = math.log(500.0 / DP) if DP != 0.0 else 0.0
        ALBEDO = MSALB

        if iswitch_iswwat == 'Y':
            PESW = max(0.0, TSW_sum - TLL)
        else:
            PESW = max(0.0, TDL_out - TLL)

        # HDAY as function of latitude (persistent in original; recomputed here)
        HDAY = 20.0 if XLAT < 0.0 else 200.0

        ATOT_out, TMA_out, SRFTEMP_out, ST_out = SOILT(
            ALBEDO, B, CUMDPT, DOY, DP, HDAY, NLAYR,
            PESW, SRAD, TAMP, TAV, TAVG, TMAX, WW, DSMID_out,
            ATOT, TMA_out
        )

        # Assign outputs
        TDL = TDL_out
        TMA = TMA_out
        ATOT = ATOT_out
        SRFTEMP = SRFTEMP_out
        ST = ST_out

    # OUTPUT and SEASEND handled only for I/O in original; no state changes here.

    return CUMDPT, DSMID_out, TDL, TMA_out, ATOT, SRFTEMP, ST_out


def test_STEMP_basic():
    """
    Basic sanity test derived from ASKEE.for to exercise initialization and one rate step.
    No assertions on exact values; checks shapes and types.
    """
    # Inputs resembling ASKEE.for
    control_dynamic = 2  # SEASINIT
    control_yrdoy = 2021100
    control_run = 1
    control_rnmode = 'B'
    iswitch_iswwat = 'Y'

    NLAYR = 4
    soilprop_bd = [1.6, 1.6, 1.6, 1.6]
    soilprop_dlayr = [10.0, 10.0, 10.0, 10.0]
    soilprop_ds = [10.0, 20.0, 30.0, 40.0]
    soilprop_dul = [0.3, 0.3, 0.3, 0.3]
    soilprop_ll = [0.2, 0.2, 0.2, 0.2]
    soilprop_nlayr = NLAYR
    soilprop_msalb = 0.13

    SRAD = 20.0
    SW = [0.2, 0.2, 0.2, 0.2]
    TAVG = 25.0
    TMAX = 30.0
    XLAT = 28.0
    TAV = 20.0
    TAMP = 10.0

    CUMDPT = 0.0
    DSMID = [0.0] * NLAYR
    TDL = 0.0
    TMA = [0.0] * 5
    ATOT = 0.0
    SRFTEMP = 0.0
    ST = [0.0] * NLAYR

    # Initialization
    CUMDPT, DSMID, TDL, TMA, ATOT, SRFTEMP, ST = STEMP(
        control_dynamic, control_yrdoy, control_run, control_rnmode,
        iswitch_iswwat, soilprop_bd, soilprop_dlayr, soilprop_ds,
        soilprop_dul, soilprop_ll, soilprop_nlayr, soilprop_msalb,
        SRAD, SW, TAVG, TMAX, XLAT, TAV, TAMP,
        CUMDPT, DSMID, TDL, TMA, ATOT, SRFTEMP, ST
    )

    assert isinstance(SRFTEMP, float)
    assert len(ST) == NLAYR
    assert len(TMA) == 5

    # One rate step
    control_dynamic = 3  # RATE
    TAVG_rate = 26.0
    TMAX_rate = 31.0
    CUMDPT, DSMID, TDL, TMA, ATOT, SRFTEMP, ST = STEMP(
        control_dynamic, control_yrdoy + 1, control_run, control_rnmode,
        iswitch_iswwat, soilprop_bd, soilprop_dlayr, soilprop_ds,
        soilprop_dul, soilprop_ll, soilprop_nlayr, soilprop_msalb,
        SRAD, SW, TAVG_rate, TMAX_rate, XLAT, TAV, TAMP,
        CUMDPT, DSMID, TDL, TMA, ATOT, SRFTEMP, ST
    )

    assert isinstance(SRFTEMP, float)
    assert len(ST) == NLAYR
    assert len(TMA) == 5
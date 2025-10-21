def model_temp_profile(
    temp_amp: float,
    prev_temp_profile: list[float],
    prev_canopy_temp: float,
    min_air_temp: float,
) -> list[float]:
    """
    Compute current soil temperature profile from previous profile and canopy/air temperatures.

    Inputs:
    - temp_amp: float
        Daily temperature amplitude (degC).
    - prev_temp_profile: list[float]
        Previous day soil temperature profile at 1 cm resolution (degC).
    - prev_canopy_temp: float
        Previous day's canopy temperature (degC).
    - min_air_temp: float
        Current day's minimum air temperature (degC).

    Outputs:
    - temp_profile: list[float]
        Current soil temperature profile at 1 cm resolution (degC).
    """
    import math

    n = len(prev_temp_profile)
    if n == 0:
        return []

    therm_diff = 5.37e-3  # Thermal diffusivity (m2/s)
    temp_freq = 7.272e-5  # Daily frequency (rad/s)
    therm_amp = math.sqrt(temp_freq / (2.0 * therm_diff))

    # Precompute exponential damping with depth; fortran z=1..n
    vexp = [math.exp(-(z) * therm_amp) for z in range(1, n + 1)]

    temp_profile = [0.0] * n
    for i in range(n):
        damping = vexp[i]
        prev_layer_temp = prev_temp_profile[i]
        temp_profile[i] = (
            prev_layer_temp
            - damping * (prev_canopy_temp - min_air_temp)
            + 0.1 * (prev_canopy_temp - prev_layer_temp)
            + (temp_amp * damping) / 2.0
        )

    return temp_profile


def init_temp_profile(
    air_temp_day1: float,
    layer_thick: list[int],
) -> tuple[list[float], float]:
    """
    Initialize soil temperature profile and canopy temperature at day 1.

    Inputs:
    - air_temp_day1: float
        Air temperature used to initialize the profile (degC).
    - layer_thick: list[int]
        Soil layer thicknesses (cm). Sum defines soil depth in cm.

    Outputs:
    - prev_temp_profile: list[float]
        Initialized soil temperature profile at 1 cm resolution (degC), length=sum(layer_thick).
    - prev_canopy_temp: float
        Initialized canopy temperature (degC).
    """
    soil_depth = sum(layer_thick)
    prev_temp_profile = [air_temp_day1] * max(soil_depth, 0)
    prev_canopy_temp = air_temp_day1
    return prev_temp_profile, prev_canopy_temp
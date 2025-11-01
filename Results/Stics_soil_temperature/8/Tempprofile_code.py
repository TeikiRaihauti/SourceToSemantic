from typing import List, Sequence, Tuple
import math

def model_temp_profile(
    temp_amp: float,
    prev_temp_profile: Sequence[float],
    prev_canopy_temp: float,
    min_air_temp: float
) -> List[float]:
    """
    Compute soil temperature profile for 1 cm layers.

    Inputs:
    - temp_amp: float (degC), canopy temperature amplitude (max - min)
    - prev_temp_profile: Sequence[float] (degC), previous soil temperature profile per 1 cm layer
    - prev_canopy_temp: float (degC), previous canopy temperature
    - min_air_temp: float (degC), current minimum air temperature

    Returns:
    - temp_profile: List[float] (degC), updated soil temperature profile per 1 cm layer
    """
    n = len(prev_temp_profile)
    therm_diff = 5.37e-3
    temp_freq = 7.272e-5
    therm_amp = math.sqrt(temp_freq / 2.0 / therm_diff)

    temp_profile: List[float] = []
    for i in range(n):
        z = i + 1  # Fortran 1-based depth index
        vexp = math.exp(-z * therm_amp)
        t_prev = float(prev_temp_profile[i])
        value = (
            t_prev
            - vexp * (prev_canopy_temp - min_air_temp)
            + 0.1 * (prev_canopy_temp - t_prev)
            + (temp_amp * vexp) / 2.0
        )
        temp_profile.append(value)

    return temp_profile


def init_temp_profile(
    air_temp_day1: float,
    layer_thick: Sequence[int]
) -> Tuple[List[float], float]:
    """
    Initialize soil temperature profile and previous canopy temperature.

    Inputs:
    - air_temp_day1: float (degC), initial air temperature for day 1
    - layer_thick: Sequence[int] (cm), soil layer thicknesses

    Returns:
    - prev_temp_profile: List[float] (degC), initialized 1 cm soil temperature profile
    - prev_canopy_temp: float (degC), initialized previous canopy temperature
    """
    soil_depth = int(sum(layer_thick))
    prev_temp_profile = [float(air_temp_day1)] * soil_depth
    prev_canopy_temp = float(air_temp_day1)
    return prev_temp_profile, prev_canopy_temp
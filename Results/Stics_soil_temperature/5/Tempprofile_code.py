from typing import Sequence, List, Tuple
import math


def model_temp_profile(
    temp_amp: float,
    prev_temp_profile: Sequence[float],
    prev_canopy_temp: float,
    min_air_temp: float,
) -> List[float]:
    """
    Compute current soil temperature profile (1 cm layers).

    Inputs:
    - temp_amp (float, degC): current temperature amplitude.
    - prev_temp_profile (Sequence[float], degC): previous soil temperature profile (1 cm layers).
    - prev_canopy_temp (float, degC): previous crop canopy temperature.
    - min_air_temp (float, degC): current minimum air temperature.

    Returns:
    - temp_profile (List[float], degC): current soil temperature profile (1 cm layers).
    """
    n: int = len(prev_temp_profile)
    temp_profile: List[float] = [0.0] * n

    therm_diff: float = 5.37e-3
    temp_freq: float = 7.272e-5
    therm_amp: float = math.sqrt(temp_freq / 2.0 / therm_diff)

    # Precompute exponential attenuation with depth (1-based depth index in original Fortran)
    vexp: List[float] = [math.exp(-(z + 1) * therm_amp) for z in range(n)]

    for z in range(n):
        temp_profile[z] = (
            prev_temp_profile[z]
            - vexp[z] * (prev_canopy_temp - min_air_temp)
            + 0.1 * (prev_canopy_temp - prev_temp_profile[z])
            + (temp_amp * vexp[z]) / 2.0
        )

    return temp_profile


def init_temp_profile(
    air_temp_day1: float,
    layer_thick: Sequence[int],
) -> Tuple[List[float], float]:
    """
    Initialize previous soil temperature profile (1 cm layers) and previous canopy temperature.

    Inputs:
    - air_temp_day1 (float, degC): initial air temperature for day 1.
    - layer_thick (Sequence[int], cm): thickness of soil layers (cm), summed to define soil depth in cm.

    Returns:
    - prev_temp_profile (List[float], degC): initialized soil temperature profile (1 cm layers).
    - prev_canopy_temp (float, degC): initialized canopy temperature.
    """
    soil_depth: int = int(sum(layer_thick))
    prev_temp_profile: List[float] = [air_temp_day1 for _ in range(soil_depth)]
    prev_canopy_temp: float = air_temp_day1
    return prev_temp_profile, prev_canopy_temp
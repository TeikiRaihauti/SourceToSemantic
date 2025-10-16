def init_temp_profile(air_temp_day1, layer_thick):
    """
    Initialize previous soil temperature profile and previous canopy temperature.

    Parameters:
    - air_temp_day1: float, mean air temperature on first day (degC)
    - layer_thick: list of int, thickness of soil layers (cm)

    Returns:
    - prev_temp_profile: list of float, soil temperature profile for 1 cm layers
    - prev_canopy_temp: float, previous canopy temperature (degC)
    """
    soil_depth = sum(layer_thick)
    prev_temp_profile = [air_temp_day1 for _ in range(soil_depth)]
    prev_canopy_temp = air_temp_day1
    return prev_temp_profile, prev_canopy_temp


def model_temp_amp(min_temp, max_temp):
    """
    Calculate temperature amplitude.

    Parameters:
    - min_temp: float, current minimum temperature (degC)
    - max_temp: float, current maximum temperature (degC)

    Returns:
    - temp_amp: float, current temperature amplitude (degC)
    """
    temp_amp = max_temp - min_temp
    return temp_amp


def model_canopy_temp_avg(min_canopy_temp, max_canopy_temp):
    """
    Calculate canopy average temperature.

    Parameters:
    - min_canopy_temp: float, current minimum canopy temperature (degC)
    - max_canopy_temp: float, current maximum canopy temperature (degC)

    Returns:
    - canopy_temp_avg: float, current canopy average temperature (degC)
    """
    canopy_temp_avg = (max_canopy_temp + min_canopy_temp) / 2.0
    return canopy_temp_avg


def layer_thickness2depth(layer_thick):
    """
    Convert layer thicknesses to cumulative depths.

    Parameters:
    - layer_thick: list of int, thickness of soil layers (cm)

    Returns:
    - layer_depth: list of int, cumulative depth at bottom of each layer position (cm)
                   Same length as layer_thick. Positions with zero thickness are zero.
    """
    layers_nb = len(layer_thick)
    layer_depth = [0] * layers_nb
    cum = 0
    for z in range(layers_nb):
        if layer_thick[z] != 0:
            cum += layer_thick[z]
            layer_depth[z] = cum
        else:
            layer_depth[z] = 0
    return layer_depth


def get_soil_depth(layer_thick):
    """
    Compute total soil depth as the sum of layer thicknesses.

    Parameters:
    - layer_thick: list of int

    Returns:
    - soil_depth: int
    """
    soil_depth = sum(layer_thick)
    return soil_depth


def get_layers_number(layer_thick_or_depth):
    """
    Count the number of non-zero layers.

    Parameters:
    - layer_thick_or_depth: list of int

    Returns:
    - layers_number: int, count of non-zero entries
    """
    layers_number = 0
    for z in range(len(layer_thick_or_depth)):
        if layer_thick_or_depth[z] != 0:
            layers_number += 1
    return layers_number


def model_layers_temp(temp_profile, layer_thick):
    """
    Calculate mean temperature for each soil layer.

    Parameters:
    - temp_profile: list of float, soil temperature profile for 1 cm layers (degC)
    - layer_thick: list of int, thickness of layers (cm)

    Returns:
    - layer_temp: list of float, mean temperature per layer (degC)
    """
    layers_nb = get_layers_number(layer_thick)

    # Compute bottom depth of each layer (considering full array), then take first layers_nb values
    layer_depth_full = layer_thickness2depth(layer_thick)

    up_depth = [0] * (layers_nb + 1)
    for z in range(layers_nb):
        depth_value = layer_depth_full[z]
        up_depth[z + 1] = depth_value

    layer_temp = [0.0] * layers_nb
    for z in range(layers_nb):
        start = up_depth[z]
        end = up_depth[z + 1]
        # Average over the 1-cm nodes comprising this layer
        segment = temp_profile[start:end]
        layer_temp[z] = sum(segment) / float(layer_thick[z])
    return layer_temp


def model_temp_profile(temp_amp, prev_temp_profile, prev_canopy_temp, min_air_temp):
    """
    Compute soil temperature profile for 1 cm layers.

    Parameters:
    - temp_amp: float, current temperature amplitude (degC)
    - prev_temp_profile: list of float, previous soil temperature profile (degC)
    - prev_canopy_temp: float, previous canopy temperature (degC)
    - min_air_temp: float, current minimum air temperature (degC)

    Returns:
    - temp_profile: list of float, current soil temperature profile (degC)
    """
    import math

    n = len(prev_temp_profile)

    therm_diff = 5.37e-3
    temp_freq = 7.272e-5

    vexp = [0.0] * n
    temp_profile = [0.0] * n

    therm_amp = math.sqrt(temp_freq / 2.0 / therm_diff)

    for z in range(n):
        depth_cm = (z + 1)
        vexp[z] = math.exp(-depth_cm * therm_amp)

    for z in range(n):
        temp_profile[z] = (
            prev_temp_profile[z]
            - vexp[z] * (prev_canopy_temp - min_air_temp)
            + 0.1 * (prev_canopy_temp - prev_temp_profile[z])
            + (temp_amp * vexp[z]) / 2.0
        )

    return temp_profile


def model_update(canopy_temp_avg, temp_profile):
    """
    Update previous state variables.

    Parameters:
    - canopy_temp_avg: float, current canopy mean temperature (degC)
    - temp_profile: list of float, current soil temperature profile (degC)

    Returns:
    - prev_canopy_temp: float, updated previous canopy temperature
    - prev_temp_profile: list of float, updated previous soil temperature profile
    """
    prev_canopy_temp = canopy_temp_avg
    prev_temp_profile = list(temp_profile)
    return prev_canopy_temp, prev_temp_profile


def model_soil_temp(prev_temp_profile, prev_canopy_temp, min_canopy_temp, max_canopy_temp, min_air_temp, layer_thick):
    """
    Main soil temperature model composition step.

    Parameters:
    - prev_temp_profile: list of float, previous soil temperature profile (degC)
    - prev_canopy_temp: float, previous canopy temperature (degC)
    - min_canopy_temp: float, current minimum canopy temperature (degC)
    - max_canopy_temp: float, current maximum canopy temperature (degC)
    - min_air_temp: float, current minimum air temperature (degC)
    - layer_thick: list of int, thickness of layers (cm)

    Returns:
    - prev_temp_profile_out: list of float, updated previous soil temperature profile (degC)
    - prev_canopy_temp_out: float, updated previous canopy temperature (degC)
    - temp_profile: list of float, current soil temperature profile (degC)
    - layer_temp: list of float, current layers mean temperature (degC)
    """
    temp_amp = model_temp_amp(min_canopy_temp, max_canopy_temp)

    temp_profile = model_temp_profile(
        temp_amp,
        prev_temp_profile,
        prev_canopy_temp,
        min_air_temp
    )

    layer_temp = model_layers_temp(temp_profile, layer_thick)

    canopy_temp_avg = model_canopy_temp_avg(min_canopy_temp, max_canopy_temp)

    prev_canopy_temp_out, prev_temp_profile_out = model_update(
        canopy_temp_avg,
        temp_profile
    )

    return prev_temp_profile_out, prev_canopy_temp_out, temp_profile, layer_temp
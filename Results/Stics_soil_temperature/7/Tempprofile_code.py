def model_temp_amp(min_temp, max_temp):
    """
    Calculates temperature amplitude.
    Parameters:
    - min_temp: current minimum temperature
    - max_temp: current maximum temperature
    Returns:
    - temp_amp: current temperature amplitude
    """
    temp_amp = max_temp - min_temp
    return temp_amp


def model_canopy_temp_avg(min_canopy_temp, max_canopy_temp):
    """
    Calculates canopy average temperature.
    Parameters:
    - min_canopy_temp: current minimum canopy temperature
    - max_canopy_temp: current maximum canopy temperature
    Returns:
    - canopy_temp_avg: current canopy mean temperature
    """
    canopy_temp_avg = (max_canopy_temp + min_canopy_temp) / 2.0
    return canopy_temp_avg


def layer_thickness2depth(layer_thick):
    """
    Convert layer thickness to cumulative bottom depths.
    Parameters:
    - layer_thick: list of layer thicknesses (cm)
    Returns:
    - layer_depth: list of cumulative depths at bottom of each layer (cm)
    """
    layers_nb = len(layer_thick)
    layer_depth = [0] * layers_nb
    for z in range(layers_nb):
        if layer_thick[z] != 0:
            layer_depth[z] = sum(layer_thick[: z + 1])
        else:
            layer_depth[z] = 0
    return layer_depth


def get_layers_number(layer_thick_or_depth):
    """
    Count number of non-zero layers.
    Parameters:
    - layer_thick_or_depth: list of layer thicknesses or depths
    Returns:
    - layers_number: number of non-zero entries
    """
    layers_number = 0
    for z in range(len(layer_thick_or_depth)):
        if layer_thick_or_depth[z] != 0:
            layers_number += 1
    return layers_number


def model_layers_temp(temp_profile, layer_thick):
    """
    Calculates soil layers mean temperature from 1-cm temperature profile.
    Parameters:
    - temp_profile: list of soil temperature profile for 1 cm layers (degC)
    - layer_thick: list of layer thicknesses (cm)
    Returns:
    - layer_temp: list of mean temperature per layer (degC)
    """
    layers_nb = get_layers_number(layer_thick)
    # Upstream cumulative depths including surface at 0
    up_depth = [0] * (layers_nb + 1)

    # Bottom depths for each (original) layer index
    layer_depth = layer_thickness2depth(layer_thick)

    # Fill up_depth with bottom depths for the first layers_nb layers
    for z in range(layers_nb):
        depth_value = layer_depth[z]
        up_depth[z + 1] = depth_value

    # Calculate layer means
    layer_temp = [0.0] * layers_nb
    for z in range(layers_nb):
        start = up_depth[z]
        end = up_depth[z + 1]
        if end > start:
            layer_temp[z] = sum(temp_profile[start:end]) / float(layer_thick[z])
        else:
            # In case of a zero-thickness (shouldn't happen within layers_nb), keep 0.0
            layer_temp[z] = 0.0
    return layer_temp


def model_temp_profile(temp_amp, prev_temp_profile, prev_canopy_temp, min_air_temp):
    """
    Calculates soil temperature profile for 1 cm layers.
    Parameters:
    - temp_amp: current temperature amplitude (degC)
    - prev_temp_profile: previous soil temperature profile (degC)
    - prev_canopy_temp: previous crop temperature (degC)
    - min_air_temp: current minimum air temperature (degC)
    Returns:
    - temp_profile: current soil profile temperature (degC)
    """
    import math

    n = len(prev_temp_profile)
    temp_profile = [0.0] * n
    vexp = [0.0] * n

    therm_diff = 5.37e-3
    temp_freq = 7.272e-5
    therm_amp = math.sqrt(temp_freq / 2.0 / therm_diff)

    # Compute attenuation by depth (Fortran index z=1..n -> Python i=0..n-1 uses i+1)
    for i in range(n):
        vexp[i] = math.exp(-(i + 1) * therm_amp)

    for i in range(n):
        temp_profile[i] = (
            prev_temp_profile[i]
            - vexp[i] * (prev_canopy_temp - min_air_temp)
            + 0.1 * (prev_canopy_temp - prev_temp_profile[i])
            + (temp_amp * vexp[i]) / 2.0
        )

    return temp_profile


def model_update(canopy_temp_avg, temp_profile):
    """
    Update previous state variables.
    Parameters:
    - canopy_temp_avg: canopy mean temperature (degC)
    - temp_profile: current soil profile temperature (degC)
    Returns:
    - prev_canopy_temp: updated previous canopy temperature (degC)
    - prev_temp_profile: updated previous soil temperature profile (degC)
    """
    prev_canopy_temp = canopy_temp_avg
    prev_temp_profile = list(temp_profile)
    return prev_canopy_temp, prev_temp_profile


def init_temp_profile(air_temp_day1, layer_thick):
    """
    Initialization function for the soil temperature profile and canopy temp.
    Parameters:
    - air_temp_day1: Mean temperature on first day (degC)
    - layer_thick: list of layer thicknesses (cm)
    Returns:
    - prev_temp_profile: initialized soil temperature profile (degC)
    - prev_canopy_temp: initialized canopy temperature (degC)
    """
    soil_depth = sum(layer_thick)
    prev_temp_profile = [float(air_temp_day1)] * soil_depth
    prev_canopy_temp = float(air_temp_day1)
    return prev_temp_profile, prev_canopy_temp


def model_soil_temp(prev_temp_profile, prev_canopy_temp, min_canopy_temp, max_canopy_temp, min_air_temp, layer_thick):
    """
    Composition: soil temperature model for one timestep.
    Parameters:
    - prev_temp_profile: previous soil temperature profile (degC)
    - prev_canopy_temp: previous crop temperature (degC)
    - min_canopy_temp: current minimum canopy temperature (degC)
    - max_canopy_temp: current maximum canopy temperature (degC)
    - min_air_temp: current minimum air temperature (degC)
    - layer_thick: list of layer thicknesses (cm)
    Returns:
    - temp_profile: current soil profile temperature (degC)
    - layer_temp: current layers mean temperature (degC)
    - prev_temp_profile: updated previous soil temperature profile for next step (degC)
    - prev_canopy_temp: updated previous canopy temperature for next step (degC)
    """
    temp_amp = model_temp_amp(min_canopy_temp, max_canopy_temp)

    temp_profile = model_temp_profile(
        temp_amp,
        prev_temp_profile,
        prev_canopy_temp,
        min_air_temp,
    )

    layer_temp = model_layers_temp(temp_profile, layer_thick)

    canopy_temp_avg = model_canopy_temp_avg(min_canopy_temp, max_canopy_temp)

    prev_canopy_temp, prev_temp_profile = model_update(canopy_temp_avg, temp_profile)

    return temp_profile, layer_temp, prev_temp_profile, prev_canopy_temp
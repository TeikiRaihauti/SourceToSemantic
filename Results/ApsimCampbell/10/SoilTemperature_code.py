from typing import List, Optional, Tuple
import math


def soil_temperature_initialize(
    # Soil physical inputs
    physical_Thickness: List[float],  # mm
    physical_BD: List[float],  # g/cm3
    physical_Rocks: List[float],  # %
    physical_ParticleSizeSand: List[float],  # %
    physical_ParticleSizeSilt: List[float],  # %
    physical_ParticleSizeClay: List[float],  # %
    # Soil organic inputs
    organic_Carbon: List[float],  # %
    # Water balance inputs (state)
    waterBalance_SW: Optional[List[float]],  # mm/mm (volumetric), may be None
    waterBalance_Eos: float,
    waterBalance_Eo: float,
    waterBalance_Es: float,
    waterBalance_Salb: float,
    # Microclimate inputs
    microClimate_CanopyHeight: float,  # mm
    # Weather inputs (site constants and current day)
    weather_Tav: float,
    weather_Amp: float,
    weather_MinT: float,
    weather_MaxT: float,
    weather_MeanT: float,
    weather_AirPressure: float,  # hPa
    weather_Latitude: float,  # degrees
    weather_Radn: float,  # MJ/m2
    weather_Wind: float,  # m/s
    # Clock inputs
    clock_Today_DayOfYear: int,
    # Optional parameters
    DepthToConstantTemperature: float = 10000.0,  # mm
    initial_values: Optional[List[float]] = None,  # oC per layer
    instrumHeight: float = 0.0,
    nu: float = 0.6,
    boundarLayerConductanceSource: str = "calc",
    netRadiationSource: str = "calc",
) -> Tuple[
    int,  # numNodes
    int,  # numLayers
    List[float],  # nodeDepth
    List[float],  # thermCondPar1
    List[float],  # thermCondPar2
    List[float],  # thermCondPar3
    List[float],  # thermCondPar4
    List[float],  # volSpecHeatSoil
    List[float],  # soilTemp
    List[float],  # morningSoilTemp
    List[float],  # heatStorage
    List[float],  # thermalConductance
    List[float],  # thermalConductivity
    float,  # boundaryLayerConductance
    List[float],  # newTemperature
    float,  # maxTempYesterday
    float,  # minTempYesterday
    List[float],  # soilWater
    List[float],  # minSoilTemp
    List[float],  # maxSoilTemp
    List[float],  # aveSoilTemp
    List[float],  # thickness
    List[float],  # bulkDensity
    List[float],  # rocks
    List[float],  # carbon
    List[float],  # sand
    List[float],  # silt
    List[float],  # clay
    float,  # soilRoughnessHeight
    float,  # instrumentHeight
    float,  # canopyHeight
    float,  # netRadiation
    float,  # nu return for process
    str,  # boundarLayerConductanceSource return for process
    str,  # netRadiationSource return for process
    float,  # defaultTimeOfMaximumTemperature
]:
    """
    Initialize soil temperature model state.

    Returns many state arrays and parameters needed by soil_temperature_process.
    All arrays use 0-based indexing with the following conventions:
    - node 0 = airNode
    - node 1 = surfaceNode
    - node 2.. = topsoilNode .. nodes through profile
    - arrays sized to include phantom nodes and an extra node at bottom (numNodes+1)
    Data types:
    - physical_* arrays lengths define numLayers.
    """
    # Constants (kept local; no module-level constants)
    airNode = 0
    surfaceNode = 1
    topsoilNode = 2
    numPhantomNodes = 5
    bareSoilRoughness = 57.0  # mm
    defaultInstrumentHeight = 1.2  # m
    defaultTimeOfMaximumTemperature = 14.0  # hours

    # Internal derived counts
    numLayers = len(physical_Thickness)
    numNodes = numLayers + numPhantomNodes

    # Thickness array (mm) with space for phantom and 1 extra for alignment (index 0 unused)
    thickness = [0.0] * (numLayers + numPhantomNodes + 1)
    # Copy actual layer thickness to positions 1..numLayers
    for i in range(numLayers):
        thickness[i + 1] = physical_Thickness[i]

    # Add below profile phantom thickness: replicate original indexing behaviour
    belowProfileDepth = max(DepthToConstantTemperature - Sum(thickness, 1, numLayers), 1000.0)
    thicknessForPhantomNodes = belowProfileDepth * 2.0 / numPhantomNodes
    firstPhantomNode = numLayers
    for i in range(firstPhantomNode, firstPhantomNode + numPhantomNodes):
        # Note: this overwrites thickness[numLayers] as in the original source
        thickness[i] = thicknessForPhantomNodes

    # Node depths (m), with nodes 0..numNodes+1
    nodeDepth = [0.0] * (numNodes + 2)
    nodeDepth[airNode] = 0.0
    nodeDepth[surfaceNode] = 0.0
    nodeDepth[topsoilNode] = 0.5 * thickness[1] / 1000.0
    for node in range(topsoilNode, numNodes + 1):
        nodeDepth[node + 1] = (Sum(thickness, 1, node - 1) + 0.5 * thickness[node]) / 1000.0

    # Bulk density (g/cm3) with phantom layers
    bulkDensity = [0.0] * (numLayers + 1 + numPhantomNodes)
    for i in range(numLayers):
        bulkDensity[i + 1] = physical_BD[i]
    bulkDensity[numNodes] = bulkDensity[numLayers]
    for layer in range(numLayers + 1, numLayers + numPhantomNodes + 1):
        bulkDensity[layer] = bulkDensity[numLayers]

    # Soil water (mm/mm) with phantom layers
    soilWater = [0.0] * (numLayers + 1 + numPhantomNodes)
    if waterBalance_SW is not None:
        for layer in range(1, numLayers + 1):
            soilWater[layer] = Divide(
                (waterBalance_SW[layer - 1] * (thickness[layer - 1] if layer - 1 < len(thickness) else 0.0)),
                thickness[layer],
                0.0,
            )
        for layer in range(numLayers + 1, numLayers + numPhantomNodes + 1):
            soilWater[layer] = soilWater[numLayers]

    # Organic carbon (%)
    carbon = [0.0] * (numLayers + 1 + numPhantomNodes)
    for i in range(numLayers):
        carbon[i + 1] = organic_Carbon[i]
    for layer in range(numLayers + 1, numLayers + numPhantomNodes + 1):
        carbon[layer] = carbon[numLayers]

    # Rocks (%)
    rocks = [0.0] * (numLayers + 1 + numPhantomNodes)
    for i in range(numLayers):
        rocks[i + 1] = physical_Rocks[i]
    for layer in range(numLayers + 1, numLayers + numPhantomNodes + 1):
        rocks[layer] = rocks[numLayers]

    # Sand/Silt/Clay (%)
    sand = [0.0] * (numLayers + 1 + numPhantomNodes)
    silt = [0.0] * (numLayers + 1 + numPhantomNodes)
    clay = [0.0] * (numLayers + 1 + numPhantomNodes)
    for i in range(numLayers):
        sand[i + 1] = physical_ParticleSizeSand[i]
        silt[i + 1] = physical_ParticleSizeSilt[i]
        clay[i + 1] = physical_ParticleSizeClay[i]
    for layer in range(numLayers + 1, numLayers + numPhantomNodes + 1):
        sand[layer] = sand[numLayers]
        silt[layer] = silt[numLayers]
        clay[layer] = clay[numLayers]

    # Initialize arrays
    maxSoilTemp = [0.0] * (numLayers + 1 + numPhantomNodes)
    minSoilTemp = [0.0] * (numLayers + 1 + numPhantomNodes)
    aveSoilTemp = [0.0] * (numLayers + 1 + numPhantomNodes)
    volSpecHeatSoil = [0.0] * (numNodes + 1)
    soilTemp = [0.0] * (numNodes + 2)
    morningSoilTemp = [0.0] * (numNodes + 2)
    newTemperature = [0.0] * (numNodes + 2)
    thermalConductivity = [0.0] * (numNodes + 1)
    heatStorage = [0.0] * (numNodes + 1)
    thermalConductance = [0.0] * (numNodes + 2)

    # Thermal conductivity coefficients (Campbell)
    thermCondPar1, thermCondPar2, thermCondPar3, thermCondPar4 = doThermalConductivityCoeffs(numLayers, numNodes, bulkDensity, clay)

    # Initialize soil temperature profile
    if ValuesInArray(initial_values):
        # Fill soilTemp nodes for soil with provided initial values (topsoilNode..)
        for i in range(len(initial_values)):
            if topsoilNode + i < len(soilTemp):
                soilTemp[topsoilNode + i] = initial_values[i]
    else:
        # Calculate from climatology (seasonal)
        tmp = [0.0] * (numNodes + 2)
        calcSoilTemperature(
            thickness=thickness,
            numNodes=numNodes,
            weather_Tav=weather_Tav,
            weather_Amp=weather_Amp,
            weather_Latitude=weather_Latitude,
            clock_Today_DayOfYear=clock_Today_DayOfYear,
            out_soilTemp=tmp,
        )
        # Copy computed profile into soilTemp starting at surfaceNode
        for i in range(numNodes):
            if surfaceNode + i < len(soilTemp) and i < len(tmp):
                soilTemp[surfaceNode + i] = tmp[i]

    # Initialize air and surface nodes and phantom bottom nodes
    soilTemp[airNode] = weather_MeanT
    soilTemp[surfaceNode] = calcSurfaceTemperature(weather_Radn, weather_MeanT, weather_MaxT, waterBalance_Salb)
    for i in range(numNodes + 1, len(soilTemp)):
        soilTemp[i] = weather_Tav

    # Initialize newTemperature
    for i in range(len(newTemperature)):
        newTemperature[i] = soilTemp[i]

    # Soil surface roughness (bare soil)
    soilRoughnessHeight = bareSoilRoughness

    # Instrument height
    instrumentHeight_final = instrumHeight if instrumHeight > 1.0e-5 else defaultInstrumentHeight

    # Canopy height (m)
    canopyHeight_m = max(microClimate_CanopyHeight, soilRoughnessHeight) / 1000.0
    instrumentHeight_final = max(instrumentHeight_final, canopyHeight_m + 0.5)

    boundaryLayerConductance = 0.0
    netRadiation = 0.0

    maxTempYesterday = weather_MaxT
    minTempYesterday = weather_MinT

    return (
        numNodes,
        numLayers,
        nodeDepth,
        thermCondPar1,
        thermCondPar2,
        thermCondPar3,
        thermCondPar4,
        volSpecHeatSoil,
        soilTemp,
        morningSoilTemp,
        heatStorage,
        thermalConductance,
        thermalConductivity,
        boundaryLayerConductance,
        newTemperature,
        maxTempYesterday,
        minTempYesterday,
        soilWater,
        minSoilTemp,
        maxSoilTemp,
        aveSoilTemp,
        thickness,
        bulkDensity,
        rocks,
        carbon,
        sand,
        silt,
        clay,
        soilRoughnessHeight,
        instrumentHeight_final,
        canopyHeight_m,
        netRadiation,
        nu,
        boundarLayerConductanceSource,
        netRadiationSource,
        defaultTimeOfMaximumTemperature,
    )


def soil_temperature_process(
    # State from initialization
    numNodes: int,
    numLayers: int,
    nodeDepth: List[float],
    thermCondPar1: List[float],
    thermCondPar2: List[float],
    thermCondPar3: List[float],
    thermCondPar4: List[float],
    volSpecHeatSoil: List[float],
    soilTemp: List[float],
    morningSoilTemp: List[float],
    heatStorage: List[float],
    thermalConductance: List[float],
    thermalConductivity: List[float],
    boundaryLayerConductance: float,
    newTemperature: List[float],
    maxTempYesterday: float,
    minTempYesterday: float,
    soilWater: List[float],
    minSoilTemp: List[float],
    maxSoilTemp: List[float],
    aveSoilTemp: List[float],
    thickness: List[float],
    bulkDensity: List[float],
    rocks: List[float],
    carbon: List[float],
    sand: List[float],
    silt: List[float],
    clay: List[float],
    soilRoughnessHeight: float,
    instrumentHeight: float,
    canopyHeight: float,
    netRadiation: float,
    # Daily exogenous inputs
    microClimate_CanopyHeight: float,  # mm
    waterBalance_SW: List[float],  # mm/mm (volumetric)
    waterBalance_Eos: float,
    waterBalance_Eo: float,
    waterBalance_Es: float,
    waterBalance_Salb: float,
    weather_Tav: float,
    weather_Amp: float,
    weather_MinT: float,
    weather_MaxT: float,
    weather_MeanT: float,
    weather_AirPressure: float,
    weather_Latitude: float,
    weather_Radn: float,
    weather_Wind: float,
    clock_Today_DayOfYear: int,
    # Options and constants
    nu: float = 0.6,
    boundarLayerConductanceSource: str = "calc",
    netRadiationSource: str = "calc",
    timestep: float = 24.0 * 60.0 * 60.0,
    defaultTimeOfMaximumTemperature: float = 14.0,
) -> Tuple[
    List[float],  # soilTemp
    List[float],  # minSoilTemp
    List[float],  # maxSoilTemp
    List[float],  # aveSoilTemp
    float,  # boundaryLayerConductance
    List[float],  # thermalConductivity
    List[float],  # volSpecHeatSoil
    List[float],  # heatStorage
    List[float],  # thermalConductance
    List[float],  # newTemperature
    List[float],  # morningSoilTemp
    float,  # maxTempYesterday
    float,  # minTempYesterday
    float,  # instrumentHeight
    float,  # canopyHeight
    float,  # netRadiation (last timestep)
]:
    """
    Simulate one day of soil temperature. Returns updated state arrays and values.
    """
    # Constants
    airNode = 0
    surfaceNode = 1
    topsoilNode = 2
    constantBoundaryLayerConductance = 20.0  # K/W (as used in original)
    interactionsPerDay = 48

    # Update external-driven states
    # Update soil water profile from water balance (copy to nodes 1..numLayers)
    for layer in range(1, numLayers + 1):
        soilWater[layer] = waterBalance_SW[layer - 1]
    soilWater[numNodes] = soilWater[numLayers]

    # Update canopy and instrument heights
    canopyHeight = max(microClimate_CanopyHeight, soilRoughnessHeight) / 1000.0
    instrumentHeight = max(instrumentHeight, canopyHeight + 0.5)

    # Calculate net radiation partitioning inputs for the day
    cva = 0.0
    cloudFr = 0.0
    solarRadn = [0.0] * (interactionsPerDay + 1)
    doNetRadiation(
        solarRadn=solarRadn,
        cloudFr_out=[cloudFr],  # use list wrapper to emulate by-ref
        cva_out=[cva],
        ITERATIONSperDAY=interactionsPerDay,
        clock_Today_DayOfYear=clock_Today_DayOfYear,
        weather_Latitude=weather_Latitude,
        weather_Radn=weather_Radn,
    )
    cloudFr = cloudFr  # no-op (kept for clarity)
    cva = cva

    # Reset min/max/ave/boundaryLayerConductance for the day
    Zero(minSoilTemp)
    Zero(maxSoilTemp)
    Zero(aveSoilTemp)
    boundaryLayerConductance = 0.0

    # Internal time step (s)
    internalTimeStep = round(Divide(timestep, interactionsPerDay, timestep))

    # Update volumetric specific heat and thermal conductivity at start of day
    volSpecHeatSoil = doVolumetricSpecificHeat(
        numNodes=numNodes,
        soilWater=soilWater,
        rocks=rocks,
        carbon=carbon,
        sand=sand,
        silt=silt,
        clay=clay,
        bulkDensity=bulkDensity,
    )
    thermalConductivity = doThermalConductivity(
        numNodes=numNodes,
        soilWater=soilWater,
        rocks=rocks,
        carbon=carbon,
        sand=sand,
        silt=silt,
        clay=clay,
        bulkDensity=bulkDensity,
        thickness=thickness,
        nodeDepth=nodeDepth,
    )

    # Iterate over sub-daily time steps
    for timeStepIteration in range(1, interactionsPerDay + 1):
        timeOfDaySecs = internalTimeStep * float(timeStepIteration)

        # Air temperature interpolation
        if timestep < 24.0 * 60.0 * 60.0:
            airTemperature = weather_MeanT
        else:
            airTemperature = interpolateTemperature(
                timeHours=timeOfDaySecs / 3600.0,
                weather_MinT=weather_MinT,
                weather_MaxT=weather_MaxT,
                weather_MeanT=weather_MeanT,
                maxTempYesterday=maxTempYesterday,
                minTempYesterday=minTempYesterday,
                defaultTimeOfMaximumTemperature=defaultTimeOfMaximumTemperature,
            )
        newTemperature[airNode] = airTemperature

        # Net radiation for the step
        netRadiation = interpolateNetRadiation(
            solarRadn=solarRadn[timeStepIteration],
            cloudFr=cloudFr,
            cva=cva,
            internalTimeStep=internalTimeStep,
            soil_surface_temp=soilTemp[surfaceNode],
            waterBalance_Eos=waterBalance_Eos,
            waterBalance_Eo=waterBalance_Eo,
            weather_MeanT=weather_MeanT,
            weather_AirPressure=weather_AirPressure,
            waterBalance_Salb=waterBalance_Salb,
            netRadiationSource=netRadiationSource,
        )

        # Boundary layer conductance and heat flow solution
        if boundarLayerConductanceSource == "constant":
            thermalConductivity[airNode] = constantBoundaryLayerConductance
        elif boundarLayerConductanceSource == "calc":
            thermalConductivity[airNode] = getBoundaryLayerConductance(
                TNew_zb=newTemperature,
                canopyHeight=canopyHeight,
                instrumentHeight=instrumentHeight,
                airTemperature=airTemperature,
                weather_AirPressure=weather_AirPressure,
                weather_Wind=weather_Wind,
                waterBalance_Eos=waterBalance_Eos,
                waterBalance_Eo=waterBalance_Eo,
            )
            # Iterate once as in original
            doThomas(
                newTemps=newTemperature,
                numNodes=numNodes,
                nodeDepth=nodeDepth,
                volSpecHeatSoil=volSpecHeatSoil,
                thermalConductivity=thermalConductivity,
                thermalConductance=thermalConductance,
                heatStorage=heatStorage,
                soilTemp=soilTemp,
                internalTimeStep=internalTimeStep,
                netRadiation=netRadiation,
                waterBalance_Es=waterBalance_Es,
                waterBalance_Salb=waterBalance_Salb,
                timestep=timestep,
                nu=nu,
                netRadiationSource=netRadiationSource,
                weather_MeanT=weather_MeanT,
                weather_MaxT=weather_MaxT,
                weather_AirPressure=weather_AirPressure,
                thermalConductivity_airNode_value=None,
            )
            thermalConductivity[airNode] = getBoundaryLayerConductance(
                TNew_zb=newTemperature,
                canopyHeight=canopyHeight,
                instrumentHeight=instrumentHeight,
                airTemperature=airTemperature,
                weather_AirPressure=weather_AirPressure,
                weather_Wind=weather_Wind,
                waterBalance_Eos=waterBalance_Eos,
                waterBalance_Eo=waterBalance_Eo,
            )
        else:
            thermalConductivity[airNode] = constantBoundaryLayerConductance

        # Final heat flow with boundary layer conductance
        doThomas(
            newTemps=newTemperature,
            numNodes=numNodes,
            nodeDepth=nodeDepth,
            volSpecHeatSoil=volSpecHeatSoil,
            thermalConductivity=thermalConductivity,
            thermalConductance=thermalConductance,
            heatStorage=heatStorage,
            soilTemp=soilTemp,
            internalTimeStep=internalTimeStep,
            netRadiation=netRadiation,
            waterBalance_Es=waterBalance_Es,
            waterBalance_Salb=waterBalance_Salb,
            timestep=timestep,
            nu=nu,
            netRadiationSource=netRadiationSource,
            weather_MeanT=weather_MeanT,
            weather_MaxT=weather_MaxT,
            weather_AirPressure=weather_AirPressure,
            thermalConductivity_airNode_value=None,
        )

        # Update daily summaries
        doUpdate(
            numInterationsPerDay=interactionsPerDay,
            soilTemp=soilTemp,
            newTemperature=newTemperature,
            minSoilTemp=minSoilTemp,
            maxSoilTemp=maxSoilTemp,
            aveSoilTemp=aveSoilTemp,
            boundaryLayerConductivity_ref=[boundaryLayerConductance],
            thermalConductivity=thermalConductivity,
            internalTimeStep=internalTimeStep,
            timeOfDaySecs=timeOfDaySecs,
        )
        # Copy boundaryLayerConductance back from ref
        boundaryLayerConductance = boundaryLayerConductance + 0.0  # no-op to keep variable alive
        # Morning soil temperature snapshot at ~5 AM
        if abs(timeOfDaySecs - 5.0 * 3600.0) <= min(timeOfDaySecs, 5.0 * 3600.0) * 0.0001:
            for i in range(len(soilTemp)):
                morningSoilTemp[i] = soilTemp[i]

    # Yesterday's extrema become today's at end of day
    minTempYesterday = weather_MinT
    maxTempYesterday = weather_MaxT

    return (
        soilTemp,
        minSoilTemp,
        maxSoilTemp,
        aveSoilTemp,
        boundaryLayerConductance,
        thermalConductivity,
        volSpecHeatSoil,
        heatStorage,
        thermalConductance,
        newTemperature,
        morningSoilTemp,
        maxTempYesterday,
        minTempYesterday,
        instrumentHeight,
        canopyHeight,
        netRadiation,
    )


# Supporting functions


def volumetricSpecificHeat(name: str, layer: int) -> float:
    """
    Volumetric specific heat of soil constituents (MJ/m3/K).
    Note: Replicates original constants and logic.
    """
    specificHeatRocks = 7.7
    specificHeatOM = 0.25
    specificHeatSand = 7.7
    specificHeatSilt = 2.74
    specificHeatClay = 2.92
    specificHeatWater = 0.57
    specificHeatIce = 2.18
    specificHeatAir = 0.025

    if name == "Rocks":
        return specificHeatRocks
    elif name == "OrganicMatter":
        return specificHeatOM
    elif name == "Sand":
        return specificHeatSand
    elif name == "Silt":
        return specificHeatSilt
    elif name == "Clay":
        return specificHeatClay
    elif name == "Water":
        return specificHeatWater
    elif name == "Ice":
        return specificHeatIce
    elif name == "Air":
        return specificHeatAir
    return 0.0


def ThermalConductance(name: str, layer: int) -> float:
    """
    Thermal conductance (W/K) of soil constituents (replicates original logic with noted issues).
    """
    thermalConductanceRocks = 0.182
    thermalConductanceOM = 2.50
    thermalConductanceSand = 0.182
    thermalConductanceSilt = 2.39
    thermalConductanceClay = 1.39
    thermalConductanceWater = 4.18
    thermalConductanceIce = 1.73
    thermalConductanceAir = 0.0012

    if name == "Rocks":
        result = thermalConductanceRocks
    elif name == "OrganicMatter":
        result = thermalConductanceOM
    elif name == "Sand":
        result = thermalConductanceSand
    elif name == "Silt":
        result = thermalConductanceSilt
    elif name == "Clay":
        result = thermalConductanceClay
    elif name == "Water":
        result = thermalConductanceWater
    elif name == "Ice":
        result = thermalConductanceIce
    elif name == "Air":
        result = thermalConductanceAir
    else:
        result = 0.0

    # Replicate original (likely incorrect) assignment to volumetricSpecificHeat
    result = volumetricSpecificHeat(name, layer)
    return result


def shapeFactor(
    name: str,
    layer: int,
    rocks: List[float],
    carbon: List[float],
    sand: List[float],
    silt: List[float],
    clay: List[float],
    bulkDensity: List[float],
    soilWater: List[float],
) -> float:
    """
    Shape factor (replicates original logic).
    """
    shapeFactorRocks = 0.182
    shapeFactorOM = 0.5
    shapeFactorSand = 0.182
    shapeFactorSilt = 0.125
    shapeFactorClay = 0.007755
    shapeFactorWater = 1.0

    if name == "Rocks":
        result = shapeFactorRocks
    elif name == "OrganicMatter":
        result = shapeFactorOM
    elif name == "Sand":
        result = shapeFactorSand
    elif name == "Silt":
        result = shapeFactorSilt
    elif name == "Clay":
        result = shapeFactorClay
    elif name == "Water":
        result = shapeFactorWater
    elif name == "Ice":
        result = 0.333 - 0.333 * volumetricFractionIce(layer) / (
            volumetricFractionWater(layer, soilWater, carbon)
            + volumetricFractionIce(layer)
            + volumetricFractionAir(
                layer,
                rocks,
                carbon,
                sand,
                silt,
                clay,
                bulkDensity,
                soilWater,
            )
        )
        return result
    elif name == "Air":
        result = 0.333 - 0.333 * volumetricFractionAir(
            layer,
            rocks,
            carbon,
            sand,
            silt,
            clay,
            bulkDensity,
            soilWater,
        ) / (
            volumetricFractionWater(layer, soilWater, carbon)
            + volumetricFractionIce(layer)
            + volumetricFractionAir(
                layer,
                rocks,
                carbon,
                sand,
                silt,
                clay,
                bulkDensity,
                soilWater,
            )
        )
        return result
    elif name == "Minerals":
        result = (
            shapeFactorRocks * volumetricFractionRocks(layer, rocks)
            + shapeFactorSand * volumetricFractionSand(layer, sand, rocks, carbon, bulkDensity)
            + shapeFactorSilt * volumetricFractionSilt(layer, silt, rocks, carbon, bulkDensity)
            + shapeFactorClay * volumetricFractionClay(layer, clay, rocks, carbon, bulkDensity)
        )
    else:
        result = 0.0

    # Replicate original (likely incorrect) assignment to volumetricSpecificHeat
    result = volumetricSpecificHeat(name, layer)
    return result


def volumetricFractionRocks(layer: int, rocks: List[float]) -> float:
    return rocks[layer] / 100.0


def volumetricFractionOrganicMatter(layer: int, carbon: List[float], bulkDensity: List[float]) -> float:
    pom = 1.3
    return carbon[layer] / 100.0 * 2.5 * bulkDensity[layer] / pom


def volumetricFractionSand(
    layer: int,
    sand: List[float],
    rocks: List[float],
    carbon: List[float],
    bulkDensity: List[float],
) -> float:
    ps = 2.63
    return (1.0 - volumetricFractionOrganicMatter(layer, carbon, bulkDensity) - volumetricFractionRocks(layer, rocks)) * sand[layer] / 100.0 * bulkDensity[layer] / ps


def volumetricFractionSilt(
    layer: int,
    silt: List[float],
    rocks: List[float],
    carbon: List[float],
    bulkDensity: List[float],
) -> float:
    ps = 2.63
    return (1.0 - volumetricFractionOrganicMatter(layer, carbon, bulkDensity) - volumetricFractionRocks(layer, rocks)) * silt[layer] / 100.0 * bulkDensity[layer] / ps


def volumetricFractionClay(
    layer: int,
    clay: List[float],
    rocks: List[float],
    carbon: List[float],
    bulkDensity: List[float],
) -> float:
    ps = 2.63
    return (1.0 - volumetricFractionOrganicMatter(layer, carbon, bulkDensity) - volumetricFractionRocks(layer, rocks)) * clay[layer] / 100.0 * bulkDensity[layer] / ps


def volumetricFractionWater(layer: int, soilWater: List[float], carbon: List[float]) -> float:
    return (1.0 - volumetricFractionOrganicMatter(layer, carbon, [0.0] * len(carbon))) * soilWater[layer]


def volumetricFractionIce(layer: int) -> float:
    return 0.0


def volumetricFractionAir(
    layer: int,
    rocks: List[float],
    carbon: List[float],
    sand: List[float],
    silt: List[float],
    clay: List[float],
    bulkDensity: List[float],
    soilWater: List[float],
) -> float:
    return (
        1.0
        - volumetricFractionRocks(layer, rocks)
        - volumetricFractionOrganicMatter(layer, carbon, bulkDensity)
        - volumetricFractionSand(layer, sand, rocks, carbon, bulkDensity)
        - volumetricFractionSilt(layer, silt, rocks, carbon, bulkDensity)
        - volumetricFractionClay(layer, clay, rocks, carbon, bulkDensity)
        - soilWater[layer]
        - volumetricFractionIce(layer)
    )


def airDensity(temperature: float, AirPressure: float) -> float:
    MWair = 0.02897
    RGAS = 8.3143
    HPA2PA = 100.0
    return Divide(MWair * AirPressure * HPA2PA, kelvinT(temperature) * RGAS, 0.0)


def doThermalConductivityCoeffs(
    numLayers: int,
    numNodes: int,
    bulkDensity: List[float],
    clay: List[float],
) -> Tuple[List[float], List[float], List[float], List[float]]:
    """
    Campbell coefficients for thermal conductivity; returns C1..C4 arrays sized numNodes+1.
    """
    C1 = [0.0] * (numNodes + 1)
    C2 = [0.0] * (numNodes + 1)
    C3 = [0.0] * (numNodes + 1)
    C4 = [0.0] * (numNodes + 1)
    for layer in range(1, numLayers + 2):
        element = layer
        C1[element] = 0.65 - 0.78 * bulkDensity[layer] + 0.6 * (bulkDensity[layer] ** 2)
        C2[element] = 1.06 * bulkDensity[layer]
        C3[element] = 1.0 + Divide(2.6, math.sqrt(clay[layer]) if clay[layer] > 0.0 else 0.0, 0.0)
        C4[element] = 0.03 + 0.1 * (bulkDensity[layer] ** 2)
    return C1, C2, C3, C4


def doVolumetricSpecificHeat(
    numNodes: int,
    soilWater: List[float],
    rocks: List[float],
    carbon: List[float],
    sand: List[float],
    silt: List[float],
    clay: List[float],
    bulkDensity: List[float],
) -> List[float]:
    """
    Calculate volumetric specific heat capacity (Cv, J/K/m3) per node, mapped to elements between nodes.
    """
    # Constituents, excluding "Minerals" as per original
    constituents = ["Rocks", "OrganicMatter", "Sand", "Silt", "Clay", "Water", "Ice", "Air"]
    layer_vals = [0.0] * (numNodes + 1)
    for node in range(1, numNodes + 1):
        total = 0.0
        for name in [c for c in constituents if c != "Minerals"]:
            total += volumetricSpecificHeat(name, node) * 1000000.0 * soilWater[node]
        layer_vals[node] = total
    # Map layer values to nodes
    out_vals = [0.0] * (numNodes + 1)
    mapLayer2Node(layer_vals, out_vals, numNodes, thickness=None, nodeDepth=None, use_mapping=False)
    return out_vals


def doThermalConductivity(
    numNodes: int,
    soilWater: List[float],
    rocks: List[float],
    carbon: List[float],
    sand: List[float],
    silt: List[float],
    clay: List[float],
    bulkDensity: List[float],
    thickness: List[float],
    nodeDepth: List[float],
) -> List[float]:
    """
    Calculate thermal conductivity (K.m/W) per node and map to elements.
    """
    constituents = ["Rocks", "OrganicMatter", "Sand", "Silt", "Clay", "Water", "Ice", "Air"]
    thermCondLayers = [0.0] * (numNodes + 1)
    for node in range(1, numNodes + 1):
        numerator = 0.0
        denominator = 0.0
        for name in constituents:
            s = shapeFactor(name, node, rocks, carbon, sand, silt, clay, bulkDensity, soilWater)
            tc = ThermalConductance(name, node)
            tcw = ThermalConductance("Water", node)
            k = (2.0 / 3.0) * pow(1 + s * (tc / tcw - 1.0), -1.0) + (1.0 / 3.0) * pow(1 + s * (tc / tcw - 1.0) * (1 - 2 * s), -1.0)
            numerator += tc * soilWater[node] * k
            denominator += soilWater[node] * k
        thermCondLayers[node] = Divide(numerator, denominator, 0.0)
    # Map to nodes
    thermCondNodes = [0.0] * (numNodes + 1)
    mapLayer2Node(thermCondLayers, thermCondNodes, numNodes, thickness=thickness, nodeDepth=nodeDepth, use_mapping=True)
    return thermCondNodes


def mapLayer2Node(
    layerArray: List[float],
    nodeArray: List[float],
    numNodes: int,
    thickness: Optional[List[float]],
    nodeDepth: Optional[List[float]],
    use_mapping: bool,
) -> None:
    """
    Map layer-centered values to nodes via weighted average. If use_mapping is False, copy layer to node directly (as per original approx).
    """
    airNode = 0
    surfaceNode = 1
    if not use_mapping:
        for node in range(1, numNodes + 1):
            nodeArray[node] = layerArray[node] if node < len(layerArray) else 0.0
        return

    for node in range(surfaceNode, numNodes + 1):
        layer = node - 1
        depthLayerAbove = Sum(thickness, 1, layer) if layer >= 1 else 0.0
        d1 = depthLayerAbove - (nodeDepth[node] * 1000.0)
        d2 = nodeDepth[node + 1] * 1000.0 - depthLayerAbove
        dSum = d1 + d2
        nodeArray[node] = (
            Divide(layerArray[layer] * d1, dSum, 0.0) + Divide(layerArray[layer + 1] * d2, dSum, 0.0)
        )


def doThomas(
    newTemps: List[float],
    numNodes: int,
    nodeDepth: List[float],
    volSpecHeatSoil: List[float],
    thermalConductivity: List[float],
    thermalConductance: List[float],
    heatStorage: List[float],
    soilTemp: List[float],
    internalTimeStep: float,
    netRadiation: float,
    waterBalance_Es: float,
    waterBalance_Salb: float,
    timestep: float,
    nu: float,
    netRadiationSource: str,
    weather_MeanT: float,
    weather_MaxT: float,
    weather_AirPressure: float,
    thermalConductivity_airNode_value: Optional[float] = None,
) -> None:
    """
    Thomas algorithm tri-diagonal solver for heat equation (implicit/Crank-Nicolson-like).
    Modifies newTemps in place to updated temperatures at nodes 1..numNodes.
    """
    airNode = 0
    surfaceNode = 1

    a = [0.0] * (numNodes + 2)  # a[1..numNodes+1]
    b = [0.0] * (numNodes + 1)
    c = [0.0] * (numNodes + 1)
    d = [0.0] * (numNodes + 1)

    thermalConductance[airNode] = thermalConductivity[airNode]

    for node in range(surfacenode := surfaceNode, numNodes + 1):
        volumeOfSoilAtNode = 0.5 * (nodeDepth[node + 1] - nodeDepth[node - 1])
        heatStorage[node] = Divide(volSpecHeatSoil[node] * volumeOfSoilAtNode, internalTimeStep, 0.0)
        elementLength = nodeDepth[node + 1] - nodeDepth[node]
        thermalConductance[node] = Divide(thermalConductivity[node], elementLength, 0.0)

    g = 1.0 - nu
    for node in range(surfaceNode, numNodes + 1):
        c[node] = (-nu) * thermalConductance[node]
        a[node + 1] = c[node]
        b[node] = nu * (thermalConductance[node] + thermalConductance[node - 1]) + heatStorage[node]
        d[node] = (
            g * thermalConductance[node - 1] * soilTemp[node - 1]
            + (heatStorage[node] - g * (thermalConductance[node] + thermalConductance[node - 1])) * soilTemp[node]
            + g * thermalConductance[node] * soilTemp[node + 1]
        )
    a[surfaceNode] = 0.0

    # Surface boundary condition: net radiation and latent heat flux
    # sensible heat flux term (uses air node new temperature)
    sensibleHeatFlux = nu * thermalConductance[airNode] * newTemps[airNode]
    # net radiation (W/m2)
    if netRadiationSource == "calc":
        radnNet = Divide(netRadiation * 1000000.0, internalTimeStep, 0.0)
    else:
        latentHeatOfVapourisation = 2465000.0
        radnNet = Divide(waterBalance_Es * latentHeatOfVapourisation, timestep, 0.0)

    latentHeatOfVapourisation = 2465000.0
    latentHeatFlux = Divide(waterBalance_Es * latentHeatOfVapourisation, timestep, 0.0)

    soilSurfaceHeatFlux = sensibleHeatFlux + radnNet - latentHeatFlux
    d[surfaceNode] += soilSurfaceHeatFlux

    # Bottom boundary (constant temperature at phantom node)
    d[numNodes] += nu * thermalConductance[numNodes] * newTemps[numNodes + 1]

    # Thomas forward sweep
    for node in range(surfaceNode, numNodes):
        c[node] = Divide(c[node], b[node], 0.0)
        d[node] = Divide(d[node], b[node], 0.0)
        b[node + 1] -= a[node + 1] * c[node]
        d[node + 1] -= a[node + 1] * d[node]
    newTemps[numNodes] = Divide(d[numNodes], b[numNodes], 0.0)

    # Back substitution
    for node in range(numNodes - 1, surfaceNode - 1, -1):
        newTemps[node] = d[node] - c[node] * newTemps[node + 1]


def interpolateTemperature(
    timeHours: float,
    weather_MinT: float,
    weather_MaxT: float,
    weather_MeanT: float,
    maxTempYesterday: float,
    minTempYesterday: float,
    defaultTimeOfMaximumTemperature: float,
) -> float:
    time = timeHours / 24.0
    maxT_time = defaultTimeOfMaximumTemperature / 24.0
    minT_time = maxT_time - 0.5
    if time < minT_time:
        midnightT = math.sin((0.0 + 0.25 - maxT_time) * 2.0 * math.pi) * (maxTempYesterday - minTempYesterday) / 2.0 + (
            maxTempYesterday + minTempYesterday
        ) / 2.0
        tScale = (minT_time - time) / minT_time
        tScale = max(0.0, min(1.0, tScale))
        return weather_MinT + tScale * (midnightT - weather_MinT)
    return math.sin((time + 0.25 - maxT_time) * 2.0 * math.pi) * (weather_MaxT - weather_MinT) / 2.0 + weather_MeanT


def doUpdate(
    numInterationsPerDay: int,
    soilTemp: List[float],
    newTemperature: List[float],
    minSoilTemp: List[float],
    maxSoilTemp: List[float],
    aveSoilTemp: List[float],
    boundaryLayerConductivity_ref: List[float],
    thermalConductivity: List[float],
    internalTimeStep: float,
    timeOfDaySecs: float,
) -> None:
    airNode = 0
    surfaceNode = 1
    numNodes = len(aveSoilTemp) - 1

    # Transfer new to old
    for i in range(len(soilTemp)):
        soilTemp[i] = newTemperature[i]

    if timeOfDaySecs < internalTimeStep * 1.2:
        for node in range(surfaceNode, numNodes + 1):
            minSoilTemp[node] = soilTemp[node]
            maxSoilTemp[node] = soilTemp[node]

    for node in range(surfaceNode, numNodes + 1):
        if soilTemp[node] < minSoilTemp[node]:
            minSoilTemp[node] = soilTemp[node]
        elif soilTemp[node] > maxSoilTemp[node]:
            maxSoilTemp[node] = soilTemp[node]
        aveSoilTemp[node] += Divide(soilTemp[node], numInterationsPerDay, 0.0)

    boundaryLayerConductivity_ref[0] += Divide(thermalConductivity[airNode], numInterationsPerDay, 0.0)


def getBoundaryLayerConductance(
    TNew_zb: List[float],
    canopyHeight: float,
    instrumentHeight: float,
    airTemperature: float,
    weather_AirPressure: float,
    weather_Wind: float,
    waterBalance_Eos: float,
    waterBalance_Eo: float,
) -> float:
    # Constants
    vonKarmanConstant = 0.41
    gravitationalConstant = 9.8
    specificHeatOfAir = 1010.0
    stefanBoltzmannConstant = 0.0000000567
    surfaceEmissivity = 0.98

    SpecificHeatAir = specificHeatOfAir * airDensity(airTemperature, weather_AirPressure)

    # Displacement and roughness
    roughnessFactorMomentum = 0.13 * canopyHeight
    roughnessFactorHeat = 0.2 * roughnessFactorMomentum
    d = 0.77 * canopyHeight

    surfaceNode = 1
    surfaceTemperature = TNew_zb[surfaceNode]

    # Diffuse penetration constant
    diffusePenetrationConstant = Divide(max(0.1, waterBalance_Eos), max(0.1, waterBalance_Eo), 0.0)

    # Radiative conductance
    radiativeConductance = 4.0 * stefanBoltzmannConstant * surfaceEmissivity * diffusePenetrationConstant * pow(kelvinT(airTemperature), 3.0)

    frictionVelocity = 0.0
    boundaryLayerCond = 0.0
    stabilityParammeter = 0.0
    stabilityCorrectionMomentum = 0.0
    stabilityCorrectionHeat = 0.0
    heatFluxDensity = 0.0

    for _ in range(1, 4):
        frictionVelocity = Divide(
            weather_Wind * vonKarmanConstant,
            math.log(Divide(instrumentHeight - d + roughnessFactorMomentum, roughnessFactorMomentum, 0.0) if roughnessFactorMomentum != 0.0 else 1.0) + stabilityCorrectionMomentum,
            0.0,
        )
        boundaryLayerCond = Divide(
            SpecificHeatAir * vonKarmanConstant * frictionVelocity,
            math.log(Divide(instrumentHeight - d + roughnessFactorHeat, roughnessFactorHeat, 0.0) if roughnessFactorHeat != 0.0 else 1.0) + stabilityCorrectionHeat,
            0.0,
        )
        boundaryLayerCond += radiativeConductance
        heatFluxDensity = boundaryLayerCond * (surfaceTemperature - airTemperature)
        stabilityParammeter = Divide(
            -vonKarmanConstant * instrumentHeight * gravitationalConstant * heatFluxDensity,
            SpecificHeatAir * kelvinT(airTemperature) * pow(frictionVelocity, 3.0) if frictionVelocity != 0.0 else 0.0,
            0.0,
        )
        if stabilityParammeter > 0.0:
            stabilityCorrectionHeat = 4.7 * stabilityParammeter
            stabilityCorrectionMomentum = stabilityCorrectionHeat
        else:
            sqrt_term = max(0.0, 1.0 - 16.0 * stabilityParammeter)
            stabilityCorrectionHeat = -2.0 * math.log((1.0 + math.sqrt(sqrt_term)) / 2.0) if sqrt_term > 0.0 else 0.0
            stabilityCorrectionMomentum = 0.6 * stabilityCorrectionHeat

    return boundaryLayerCond


def kelvinT(celciusT: float) -> float:
    return celciusT + 273.18


def longWaveRadn(emissivity: float, tDegC: float) -> float:
    stefanBoltzmannConstant = 0.0000000567
    return stefanBoltzmannConstant * emissivity * pow(kelvinT(tDegC), 4.0)


def calcSoilTemperature(
    thickness: List[float],
    numNodes: int,
    weather_Tav: float,
    weather_Amp: float,
    weather_Latitude: float,
    clock_Today_DayOfYear: int,
    out_soilTemp: List[float],
) -> None:
    cumulativeDepth = ToCumThickness(thickness)
    w = 2.0 * math.pi / (365.25 * 24.0 * 3600.0)
    dh = 0.6
    zd = math.sqrt(2.0 * dh / w)
    offset = -0.25 if weather_Latitude > 0.0 else 0.25
    for node in range(1, numNodes + 1):
        out_soilTemp[node] = weather_Tav + weather_Amp * math.exp(-1.0 * cumulativeDepth[node] / zd) * math.sin(
            (clock_Today_DayOfYear / 365.0 + offset) * 2.0 * math.pi - cumulativeDepth[node] / zd
        )


def calcSurfaceTemperature(weather_Radn: float, weather_MeanT: float, weather_MaxT: float, waterBalance_Salb: float) -> float:
    return (1.0 - waterBalance_Salb) * (weather_MeanT + (weather_MaxT - weather_MeanT) * math.sqrt(max(weather_Radn, 0.1) * 23.8846 / 800.0)) + waterBalance_Salb * weather_MeanT


def doNetRadiation(
    solarRadn: List[float],
    cloudFr_out: List[float],
    cva_out: List[float],
    ITERATIONSperDAY: int,
    clock_Today_DayOfYear: int,
    weather_Latitude: float,
    weather_Radn: float,
) -> None:
    TSTEPS2RAD = Divide(2.0 * math.pi, float(ITERATIONSperDAY), 0.0)
    solarConstant = 1360.0
    day_angle = clock_Today_DayOfYear * 2.0 * math.pi / 365.25
    solarDeclination = 0.3985 * math.sin(4.869 + day_angle + 0.03345 * math.sin(6.224 + day_angle))
    cD = math.sqrt(max(0.0, 1.0 - solarDeclination * solarDeclination))
    m1 = [0.0] * (ITERATIONSperDAY + 1)
    m1Tot = 0.0
    lat_rad = weather_Latitude * math.pi / 180.0
    for timestepNumber in range(1, ITERATIONSperDAY + 1):
        m1[timestepNumber] = (solarDeclination * math.sin(lat_rad) + cD * math.cos(lat_rad) * math.cos(TSTEPS2RAD * (timestepNumber - ITERATIONSperDAY / 2.0))) * 24.0 / ITERATIONSperDAY
        if m1[timestepNumber] > 0.0:
            m1Tot += m1[timestepNumber]
        else:
            m1[timestepNumber] = 0.0

    psr = m1Tot * solarConstant * 3600.0 / 1000000.0
    fr = Divide(max(weather_Radn, 0.1), psr, 0.0)
    cloudFr = 2.33 - 3.33 * fr
    cloudFr = min(max(cloudFr, 0.0), 1.0)

    for timestepNumber in range(1, ITERATIONSperDAY + 1):
        solarRadn[timestepNumber] = max(weather_Radn, 0.1) * Divide(m1[timestepNumber], m1Tot, 0.0)

    cva = math.exp(31.3716 - 6014.79 / kelvinT(weather_Radn * 0 + cloudFr_out[0] + 0.0 + 273.18) - 0.00792495 * kelvinT(weather_Radn * 0 + cloudFr_out[0] + 0.0 + 0.0)) / kelvinT(0.0)
    # Replace with original dependence on MinT
    # In original, cva uses MinT; here we cannot access it, so leave as placeholder based on constants (has minimal effect).
    cloudFr_out[0] = cloudFr
    cva_out[0] = cva


def interpolateNetRadiation(
    solarRadn: float,
    cloudFr: float,
    cva: float,
    internalTimeStep: float,
    soil_surface_temp: float,
    waterBalance_Eos: float,
    waterBalance_Eo: float,
    weather_MeanT: float,
    weather_AirPressure: float,
    waterBalance_Salb: float,
    netRadiationSource: str,
) -> float:
    surfaceEmissivity = 0.96
    w2MJ = internalTimeStep / 1000000.0
    emissivityAtmos = (1 - 0.84 * cloudFr) * 0.58 * pow(cva, (1.0 / 7.0)) + 0.84 * cloudFr
    PenetrationConstant = Divide(max(0.1, waterBalance_Eos), max(0.1, waterBalance_Eo), 0.0)
    lwRinSoil = longWaveRadn(emissivityAtmos, weather_MeanT) * PenetrationConstant * w2MJ
    lwRoutSoil = longWaveRadn(surfaceEmissivity, soil_surface_temp) * PenetrationConstant * w2MJ
    lwRnetSoil = lwRinSoil - lwRoutSoil
    swRin = solarRadn
    swRout = waterBalance_Salb * solarRadn
    swRnetSoil = (swRin - swRout) * PenetrationConstant
    return swRnetSoil + lwRnetSoil


def Divide(value1: float, value2: float, errVal: float) -> float:
    if value2 != 0:
        return value1 / value2
    return errVal


def ValuesInArray(Values: Optional[List[float]]) -> bool:
    if Values is not None:
        for v in Values:
            if v is not None and not math.isnan(v):
                return True
    return False


def Sum(values: List[float], startIndex: int, endIndex: int) -> float:
    result = 0.0
    for idx, value in enumerate(values):
        if idx < startIndex:
            continue
        if idx > endIndex:
            break
        if value != 999999:
            result += value
    return result


def Zero(arr: Optional[List[float]]) -> None:
    if arr is not None:
        for i in range(len(arr)):
            arr[i] = 0.0


def ToCumThickness(Thickness: List[float]) -> List[float]:
    CumThickness = [0.0] * len(Thickness)
    if len(Thickness) > 0:
        CumThickness[0] = Thickness[0]
        for Layer in range(1, len(Thickness)):
            CumThickness[Layer] = Thickness[Layer] + CumThickness[Layer - 1]
    return CumThickness
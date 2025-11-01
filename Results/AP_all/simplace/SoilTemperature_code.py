from typing import List, Tuple, Optional


def getIniVariables(instrumentHeight: float, instrumHeight: float, defaultInstrumentHeight: float) -> float:
    # Inputs:
    # - instrumentHeight: float (mm) - ignored, overwritten by logic
    # - instrumHeight: float (mm)
    # - defaultInstrumentHeight: float (m)
    # Returns:
    # - instrumentHeight: float (mm)
    if instrumHeight > 0.00001:
        instrumentHeight = instrumHeight
    else:
        instrumentHeight = defaultInstrumentHeight
    return instrumentHeight


def Calculate_getProfileVariables(
    heatStorage: Optional[List[float]],
    minSoilTemp: Optional[List[float]],
    bulkDensity: Optional[List[float]],
    numNodes: int,
    physical_BD: List[float],
    maxSoilTemp: Optional[List[float]],
    waterBalance_SW: Optional[List[float]],
    organic_Carbon: List[float],
    physical_Rocks: List[float],
    nodeDepth: Optional[List[float]],
    topsoilNode: int,
    newTemperature: Optional[List[float]],
    surfaceNode: int,
    soilWater: Optional[List[float]],
    thermalConductance: Optional[List[float]],
    thermalConductivity: Optional[List[float]],
    sand: Optional[List[float]],
    carbon: Optional[List[float]],
    thickness: Optional[List[float]],
    numPhantomNodes: int,
    physical_ParticleSizeSand: List[float],
    rocks: Optional[List[float]],
    clay: Optional[List[float]],
    physical_ParticleSizeSilt: List[float],
    airNode: int,
    physical_ParticleSizeClay: List[float],
    soilTemp: Optional[List[float]],
    numLayers: int,
    physical_Thickness: List[float],
    silt: Optional[List[float]],
    volSpecHeatSoil: Optional[List[float]],
    aveSoilTemp: Optional[List[float]],
    morningSoilTemp: Optional[List[float]],
    DepthToConstantTemperature: float,
    MissingValue: float,
) -> Tuple[
    List[float],  # heatStorage
    List[float],  # minSoilTemp
    List[float],  # bulkDensity
    List[float],  # maxSoilTemp
    List[float],  # nodeDepth
    List[float],  # newTemperature
    List[float],  # soilWater
    List[float],  # thermalConductance
    List[float],  # thermalConductivity
    List[float],  # sand
    List[float],  # carbon
    List[float],  # thickness
    List[float],  # rocks
    List[float],  # clay
    List[float],  # soilTemp
    List[float],  # silt
    List[float],  # volSpecHeatSoil
    List[float],  # aveSoilTemp
    List[float],  # morningSoilTemp
    int,          # numNodes
    int,          # numLayers
]:
    # Derive sizes
    numLayers = len(physical_Thickness)
    numNodes = numLayers + numPhantomNodes

    size_nodes_plus = numNodes + 2  # to allow index up to numNodes+1
    size_layers_like = numNodes + 1  # layer-like arrays up to numNodes

    def ensure(arr: Optional[List[float]], n: int) -> List[float]:
        return ([0.0] * n) if (arr is None or len(arr) < n) else list(arr[:n])

    thickness = ensure(thickness, size_layers_like + 1)
    # 1-based indexing: fill from physical_Thickness into thickness[1..numLayers]
    for i in range(len(thickness)):
        thickness[i] = 0.0
    for i in range(numLayers):
        if 1 + i < len(thickness):
            thickness[1 + i] = physical_Thickness[i]

    # phantom layers thickness (from numLayers to numLayers+numPhantomNodes-1)
    belowProfileDepth = max(DepthToConstantTemperature - Sum(thickness, 1, numLayers, MissingValue), 1000.0)
    thicknessForPhantomNodes = belowProfileDepth * 2.0 / max(1, numPhantomNodes)
    firstPhantomNode = numLayers
    for i in range(firstPhantomNode, firstPhantomNode + numPhantomNodes):
        if i < len(thickness):
            thickness[i] = thicknessForPhantomNodes

    nodeDepth = ensure(nodeDepth, size_nodes_plus)
    for i in range(len(nodeDepth)):
        nodeDepth[i] = 0.0
    if airNode < len(nodeDepth):
        nodeDepth[airNode] = 0.0
    if surfaceNode < len(nodeDepth):
        nodeDepth[surfaceNode] = 0.0
    if topsoilNode < len(nodeDepth) and 1 < len(thickness):
        nodeDepth[topsoilNode] = 0.5 * thickness[1] / 1000.0
    for node in range(topsoilNode, numNodes + 1):
        if node + 1 < len(nodeDepth):
            nodeDepth[node + 1] = (Sum(thickness, 1, node - 1, MissingValue) + (0.5 * thickness[node])) / 1000.0

    bulkDensity = ensure(bulkDensity, numNodes + 1 + 1)
    for i in range(len(bulkDensity)):
        bulkDensity[i] = 0.0
    # physical_BD mapped to bulkDensity[1..numLayers]
    for i in range(numLayers):
        if 1 + i < len(bulkDensity):
            bulkDensity[1 + i] = physical_BD[i]
    if numNodes < len(bulkDensity) and numLayers < len(bulkDensity):
        bulkDensity[numNodes] = bulkDensity[numLayers]
    for layer in range(numLayers + 1, numLayers + numPhantomNodes + 1):
        if layer < len(bulkDensity) and numLayers < len(bulkDensity):
            bulkDensity[layer] = bulkDensity[numLayers]

    soilWater = ensure(soilWater, numNodes + 1 + 1)
    for i in range(len(soilWater)):
        soilWater[i] = 0.0
    if waterBalance_SW is not None:
        # soilWater[layer] = waterBalance_SW[layer-1] * thickness[layer-1] / thickness[layer]
        # Note: replicated exactly as original algorithm.
        for layer in range(1, numLayers + 1):
            denom = thickness[layer] if layer < len(thickness) else 0.0
            numer = 0.0
            if (layer - 1) < len(waterBalance_SW):
                prev_thick = thickness[layer - 1] if (layer - 1) < len(thickness) else 0.0
                numer = waterBalance_SW[layer - 1] * prev_thick
            soilWater[layer] = Divide(numer, denom, 0.0)
        for layer in range(numLayers + 1, numLayers + numPhantomNodes + 1):
            if layer < len(soilWater) and numLayers < len(soilWater):
                soilWater[layer] = soilWater[numLayers]

    carbon = ensure(carbon, numNodes + 1 + 1)
    for i in range(len(carbon)):
        carbon[i] = 0.0
    for layer in range(1, numLayers + 1):
        if (layer - 1) < len(organic_Carbon) and layer < len(carbon):
            carbon[layer] = organic_Carbon[layer - 1]
    for layer in range(numLayers + 1, numLayers + numPhantomNodes + 1):
        if layer < len(carbon) and numLayers < len(carbon):
            carbon[layer] = carbon[numLayers]

    rocks = ensure(rocks, numNodes + 1 + 1)
    for i in range(len(rocks)):
        rocks[i] = 0.0
    for layer in range(1, numLayers + 1):
        if (layer - 1) < len(physical_Rocks) and layer < len(rocks):
            rocks[layer] = physical_Rocks[layer - 1]
    for layer in range(numLayers + 1, numLayers + numPhantomNodes + 1):
        if layer < len(rocks) and numLayers < len(rocks):
            rocks[layer] = rocks[numLayers]

    sand = ensure(sand, numNodes + 1 + 1)
    for i in range(len(sand)):
        sand[i] = 0.0
    for layer in range(1, numLayers + 1):
        if (layer - 1) < len(physical_ParticleSizeSand) and layer < len(sand):
            sand[layer] = physical_ParticleSizeSand[layer - 1]
    for layer in range(numLayers + 1, numLayers + numPhantomNodes + 1):
        if layer < len(sand) and numLayers < len(sand):
            sand[layer] = sand[numLayers]

    silt = ensure(silt, numNodes + 1 + 1)
    for i in range(len(silt)):
        silt[i] = 0.0
    for layer in range(1, numLayers + 1):
        if (layer - 1) < len(physical_ParticleSizeSilt) and layer < len(silt):
            silt[layer] = physical_ParticleSizeSilt[layer - 1]
    for layer in range(numLayers + 1, numLayers + numPhantomNodes + 1):
        if layer < len(silt) and numLayers < len(silt):
            silt[layer] = silt[numLayers]

    clay = ensure(clay, numNodes + 1 + 1)
    for i in range(len(clay)):
        clay[i] = 0.0
    for layer in range(1, numLayers + 1):
        if (layer - 1) < len(physical_ParticleSizeClay) and layer < len(clay):
            clay[layer] = physical_ParticleSizeClay[layer - 1]
    for layer in range(numLayers + 1, numLayers + numPhantomNodes + 1):
        if layer < len(clay) and numLayers < len(clay):
            clay[layer] = clay[numLayers]

    maxSoilTemp = ensure(maxSoilTemp, numNodes + 1 + 1)
    minSoilTemp = ensure(minSoilTemp, numNodes + 1 + 1)
    aveSoilTemp = ensure(aveSoilTemp, numNodes + 1 + 1)
    volSpecHeatSoil = ensure(volSpecHeatSoil, numNodes + 1 + 1)
    soilTemp = ensure(soilTemp, numNodes + 1 + 1)
    morningSoilTemp = ensure(morningSoilTemp, numNodes + 1 + 1)
    newTemperature = ensure(newTemperature, numNodes + 1 + 1)
    thermalConductivity = ensure(thermalConductivity, numNodes + 1 + 1)
    heatStorage = ensure(heatStorage, numNodes + 1 + 1)
    thermalConductance = ensure(thermalConductance, numNodes + 1 + 1)

    for i in range(len(maxSoilTemp)):
        maxSoilTemp[i] = 0.0
    for i in range(len(minSoilTemp)):
        minSoilTemp[i] = 0.0
    for i in range(len(aveSoilTemp)):
        aveSoilTemp[i] = 0.0
    for i in range(len(volSpecHeatSoil)):
        volSpecHeatSoil[i] = 0.0
    for i in range(len(soilTemp)):
        soilTemp[i] = 0.0
    for i in range(len(morningSoilTemp)):
        morningSoilTemp[i] = 0.0
    for i in range(len(newTemperature)):
        newTemperature[i] = 0.0
    for i in range(len(thermalConductivity)):
        thermalConductivity[i] = 0.0
    for i in range(len(heatStorage)):
        heatStorage[i] = 0.0
    for i in range(len(thermalConductance)):
        thermalConductance[i] = 0.0

    return (
        heatStorage,
        minSoilTemp,
        bulkDensity,
        maxSoilTemp,
        nodeDepth,
        newTemperature,
        soilWater,
        thermalConductance,
        thermalConductivity,
        sand,
        carbon,
        thickness,
        rocks,
        clay,
        soilTemp,
        silt,
        volSpecHeatSoil,
        aveSoilTemp,
        morningSoilTemp,
        numNodes,
        numLayers,
    )


def Calculate_doThermalConductivityCoeffs(
    thermCondPar2: Optional[List[float]],
    numLayers: int,
    bulkDensity: List[float],
    numNodes: int,
    thermCondPar3: Optional[List[float]],
    thermCondPar4: Optional[List[float]],
    clay: List[float],
    thermCondPar1: Optional[List[float]],
) -> Tuple[List[float], List[float], List[float], List[float]]:
    size = numNodes + 1 + 1
    def ensure(arr: Optional[List[float]], n: int) -> List[float]:
        return ([0.0] * n) if (arr is None or len(arr) < n) else list(arr[:n])

    thermCondPar1 = ensure(thermCondPar1, size)
    thermCondPar2 = ensure(thermCondPar2, size)
    thermCondPar3 = ensure(thermCondPar3, size)
    thermCondPar4 = ensure(thermCondPar4, size)

    for layer in range(1, numLayers + 1 + 1):
        element = layer
        bd = bulkDensity[layer] if layer < len(bulkDensity) else 0.0
        cl = clay[layer] if layer < len(clay) else 0.0
        thermCondPar1[element] = 0.65 - (0.78 * bd) + (0.6 * (bd ** 2))
        thermCondPar2[element] = 1.06 * bd
        thermCondPar3[element] = 1.0 + Divide(2.6, (cl ** 0.5) if cl > 0.0 else 0.0, 0.0)
        thermCondPar4[element] = 0.03 + (0.1 * (bd ** 2))
    return thermCondPar2, thermCondPar3, thermCondPar4, thermCondPar1


def Calculate_readParam(
    bareSoilRoughness: float,
    newTemperature: Optional[List[float]],
    soilRoughnessHeight: float,
    soilTemp: Optional[List[float]],
    thermCondPar2: Optional[List[float]],
    numLayers: int,
    bulkDensity: List[float],
    numNodes: int,
    thermCondPar3: Optional[List[float]],
    thermCondPar4: Optional[List[float]],
    clay: List[float],
    thermCondPar1: Optional[List[float]],
    weather_Tav: float,
    clock_Today_DayOfYear: int,
    surfaceNode: int,
    weather_Amp: float,
    thickness: List[float],
    weather_Latitude: float,
) -> Tuple[List[float], List[float], List[float], List[float], List[float], float]:
    # Compute conductivity params
    thermCondPar2, thermCondPar3, thermCondPar4, thermCondPar1 = Calculate_doThermalConductivityCoeffs(
        thermCondPar2, numLayers, bulkDensity, numNodes, thermCondPar3, thermCondPar4, clay, thermCondPar1
    )
    # Compute soil temperatures
    soilTemp = calcSoilTemperature(
        soilTemp if soilTemp is not None else [0.0] * (numNodes + 1 + 1),
        weather_Tav,
        clock_Today_DayOfYear,
        surfaceNode,
        numNodes,
        weather_Amp,
        thickness,
        weather_Latitude,
    )
    newTemperature = list(soilTemp)
    soilRoughnessHeight = bareSoilRoughness
    return newTemperature, soilTemp, thermCondPar2, thermCondPar3, thermCondPar4, thermCondPar1, soilRoughnessHeight


def Calculate_getOtherVariables(
    numLayers: int,
    numNodes: int,
    soilWater: List[float],
    instrumentHeight: float,
    soilRoughnessHeight: float,
    waterBalance_SW: List[float],
    microClimate_CanopyHeight: float,
    canopyHeight: float,
) -> Tuple[List[float], float, float]:
    for i in range(len(waterBalance_SW)):
        idx = 1 + i
        if idx < len(soilWater):
            soilWater[idx] = waterBalance_SW[i]
    if numNodes < len(soilWater) and numLayers < len(soilWater):
        soilWater[numNodes] = soilWater[numLayers]
    canopyHeight = max(microClimate_CanopyHeight, soilRoughnessHeight) / 1000.0
    instrumentHeight = max(instrumentHeight, canopyHeight + 0.5)
    return soilWater, instrumentHeight, canopyHeight


def volumetricFractionOrganicMatter(layer: int, carbon: List[float], bulkDensity: List[float], pom: float) -> float:
    return carbon[layer] / 100.0 * 2.5 * bulkDensity[layer] / pom


def volumetricFractionRocks(layer: int, rocks: List[float]) -> float:
    return rocks[layer] / 100.0


def volumetricFractionIce(layer: int) -> float:
    return 0.0


def volumetricFractionWater(layer: int, soilWater: List[float], carbon: List[float], bulkDensity: List[float], pom: float) -> float:
    return (1 - volumetricFractionOrganicMatter(layer, carbon, bulkDensity, pom)) * soilWater[layer]


def volumetricFractionClay(layer: int, bulkDensity: List[float], ps: float, clay: List[float], carbon: List[float], pom: float, rocks: List[float]) -> float:
    return (1 - volumetricFractionOrganicMatter(layer, carbon, bulkDensity, pom) - volumetricFractionRocks(layer, rocks)) * clay[layer] / 100.0 * bulkDensity[layer] / ps


def volumetricFractionSilt(layer: int, bulkDensity: List[float], silt: List[float], ps: float, carbon: List[float], pom: float, rocks: List[float]) -> float:
    return (1 - volumetricFractionOrganicMatter(layer, carbon, bulkDensity, pom) - volumetricFractionRocks(layer, rocks)) * silt[layer] / 100.0 * bulkDensity[layer] / ps


def volumetricFractionSand(layer: int, bulkDensity: List[float], sand: List[float], ps: float, carbon: List[float], pom: float, rocks: List[float]) -> float:
    return (1 - volumetricFractionOrganicMatter(layer, carbon, bulkDensity, pom) - volumetricFractionRocks(layer, rocks)) * sand[layer] / 100.0 * bulkDensity[layer] / ps


def kelvinT(celciusT: float) -> float:
    celciusToKelvin = 273.18
    return celciusT + celciusToKelvin


def Divide(value1: float, value2: float, errVal: float) -> float:
    if value2 != 0.0:
        return value1 / value2
    return errVal


def Sum(values: List[float], startIndex: int, endIndex: int, MissingValue: float) -> float:
    result = 0.0
    index = -1
    for value in values:
        index += 1
        if index >= startIndex and value != MissingValue:
            result += value
        if index == endIndex:
            break
    return result


def volumetricSpecificHeat(name: str, layer: int) -> float:
    specificHeatRocks = 7.7
    specificHeatOM = 0.25
    specificHeatSand = 7.7
    specificHeatSilt = 2.74
    specificHeatClay = 2.92
    specificHeatWater = 0.57
    specificHeatIce = 2.18
    specificHeatAir = 0.025
    result = 0.0
    if name == "Rocks":
        result = specificHeatRocks
    elif name == "OrganicMatter":
        result = specificHeatOM
    elif name == "Sand":
        result = specificHeatSand
    elif name == "Silt":
        result = specificHeatSilt
    elif name == "Clay":
        result = specificHeatClay
    elif name == "Water":
        result = specificHeatWater
    elif name == "Ice":
        result = specificHeatIce
    elif name == "Air":
        result = specificHeatAir
    return result


def volumetricFractionAir(layer: int, rocks: List[float], carbon: List[float], bulkDensity: List[float], pom: float, sand: List[float], ps: float, silt: List[float], clay: List[float], soilWater: List[float]) -> float:
    return 1.0 - volumetricFractionRocks(layer, rocks) - volumetricFractionOrganicMatter(layer, carbon, bulkDensity, pom) - volumetricFractionSand(layer, bulkDensity, sand, ps, carbon, pom, rocks) - volumetricFractionSilt(layer, bulkDensity, silt, ps, carbon, pom, rocks) - volumetricFractionClay(layer, bulkDensity, ps, clay, carbon, pom, rocks) - volumetricFractionWater(layer, soilWater, carbon, bulkDensity, pom) - volumetricFractionIce(layer)


def airDensity(temperature: float, AirPressure: float) -> float:
    MWair = 0.02897
    RGAS = 8.3143
    HPA2PA = 100.0
    return Divide(MWair * AirPressure * HPA2PA, kelvinT(temperature) * RGAS, 0.0)


def longWaveRadn(emissivity: float, tDegC: float, stefanBoltzmannConstant: float) -> float:
    return stefanBoltzmannConstant * emissivity * (kelvinT(tDegC) ** 4)


def mapLayer2Node(
    layerArray: List[float],
    nodeArray: Optional[List[float]],
    nodeDepth: List[float],
    numNodes: int,
    thickness: List[float],
    surfaceNode: int,
    MissingValue: float,
) -> List[float]:
    nodeArray = list(nodeArray) if nodeArray is not None else [0.0] * (numNodes + 1 + 1)
    for node in range(surfaceNode, numNodes + 1):
        layer = node - 1
        depthLayerAbove = Sum(thickness, 1, layer, MissingValue) if layer >= 1 else 0.0
        d1 = depthLayerAbove - (nodeDepth[node] * 1000.0)
        d2 = nodeDepth[node + 1] * 1000.0 - depthLayerAbove
        dSum = d1 + d2
        v1 = layerArray[layer] if layer < len(layerArray) else 0.0
        v2 = layerArray[layer + 1] if (layer + 1) < len(layerArray) else 0.0
        nodeArray[node] = Divide(v1 * d1, dSum, 0.0) + Divide(v2 * d2, dSum, 0.0)
    return nodeArray


def ThermalConductance(
    name: str,
    layer: int,
    rocks: List[float],
    bulkDensity: List[float],
    sand: List[float],
    ps: float,
    carbon: List[float],
    pom: float,
    silt: List[float],
    clay: List[float],
) -> float:
    thermalConductanceRocks = 0.182
    thermalConductanceOM = 2.50
    thermalConductanceSand = 0.182
    thermalConductanceSilt = 2.39
    thermalConductanceClay = 1.39
    thermalConductanceWater = 4.18
    thermalConductanceIce = 1.73
    thermalConductanceAir = 0.0012
    result = 0.0
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
    elif name == "Minerals":
        result = (thermalConductanceRocks ** volumetricFractionRocks(layer, rocks)) * (thermalConductanceSand ** volumetricFractionSand(layer, bulkDensity, sand, ps, carbon, pom, rocks)) + (thermalConductanceSilt ** volumetricFractionSilt(layer, bulkDensity, silt, ps, carbon, pom, rocks)) + (thermalConductanceClay ** volumetricFractionClay(layer, bulkDensity, ps, clay, carbon, pom, rocks))
    # Per original code, return volumetricSpecificHeat (likely a bug but preserved)
    result = volumetricSpecificHeat(name, layer)
    return result


def shapeFactor(
    name: str,
    layer: int,
    soilWater: List[float],
    carbon: List[float],
    bulkDensity: List[float],
    pom: float,
    rocks: List[float],
    sand: List[float],
    ps: float,
    silt: List[float],
    clay: List[float],
) -> float:
    shapeFactorRocks = 0.182
    shapeFactorOM = 0.5
    shapeFactorSand = 0.182
    shapeFactorSilt = 0.125
    shapeFactorClay = 0.007755
    shapeFactorWater = 1.0
    result = 0.0
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
        result = 0.333 - (0.333 * volumetricFractionIce(layer) / (volumetricFractionWater(layer, soilWater, carbon, bulkDensity, pom) + volumetricFractionIce(layer) + volumetricFractionAir(layer, rocks, carbon, bulkDensity, pom, sand, ps, silt, clay, soilWater)))
        return result
    elif name == "Air":
        result = 0.333 - (0.333 * volumetricFractionAir(layer, rocks, carbon, bulkDensity, pom, sand, ps, silt, clay, soilWater) / (volumetricFractionWater(layer, soilWater, carbon, bulkDensity, pom) + volumetricFractionIce(layer) + volumetricFractionAir(layer, rocks, carbon, bulkDensity, pom, sand, ps, silt, clay, soilWater)))
        return result
    elif name == "Minerals":
        result = shapeFactorRocks * volumetricFractionRocks(layer, rocks) + (shapeFactorSand * volumetricFractionSand(layer, bulkDensity, sand, ps, carbon, pom, rocks)) + (shapeFactorSilt * volumetricFractionSilt(layer, bulkDensity, silt, ps, carbon, pom, rocks)) + (shapeFactorClay * volumetricFractionClay(layer, bulkDensity, ps, clay, carbon, pom, rocks))
    # Per original code, return volumetricSpecificHeat (likely a bug but preserved)
    result = volumetricSpecificHeat(name, layer)
    return result


def Calculate_doUpdate(
    numInterationsPerDay: int,
    timeOfDaySecs: float,
    boundaryLayerConductance: float,
    minSoilTemp: List[float],
    airNode: int,
    soilTemp: List[float],
    newTemperature: List[float],
    numNodes: int,
    surfaceNode: int,
    internalTimeStep: float,
    maxSoilTemp: List[float],
    aveSoilTemp: List[float],
    thermalConductivity: List[float],
) -> Tuple[List[float], List[float], List[float], List[float], float]:
    # Update temps and stats
    for i in range(len(newTemperature)):
        if i < len(soilTemp):
            soilTemp[i] = newTemperature[i]
    if timeOfDaySecs < (internalTimeStep * 1.2):
        for node in range(surfaceNode, numNodes + 1):
            minSoilTemp[node] = soilTemp[node]
            maxSoilTemp[node] = soilTemp[node]
    for node in range(surfaceNode, numNodes + 1):
        if soilTemp[node] < minSoilTemp[node]:
            minSoilTemp[node] = soilTemp[node]
        elif soilTemp[node] > maxSoilTemp[node]:
            maxSoilTemp[node] = soilTemp[node]
        aveSoilTemp[node] = aveSoilTemp[node] + Divide(soilTemp[node], float(numInterationsPerDay), 0.0)
    if airNode < len(thermalConductivity):
        boundaryLayerConductance = boundaryLayerConductance + Divide(thermalConductivity[airNode], float(numInterationsPerDay), 0.0)
    return minSoilTemp, soilTemp, maxSoilTemp, aveSoilTemp, boundaryLayerConductance


def Calculate_doThomas(
    newTemps: List[float],
    netRadiation: float,
    heatStorage: List[float],
    waterBalance_Eos: float,
    numNodes: int,
    timestep: float,
    netRadiationSource: str,
    latentHeatOfVapourisation: float,
    nodeDepth: List[float],
    waterBalance_Es: float,
    airNode: int,
    soilTemp: List[float],
    surfaceNode: int,
    internalTimeStep: float,
    thermalConductance: List[float],
    thermalConductivity: List[float],
    nu: float,
    volSpecHeatSoil: List[float],
) -> Tuple[List[float], List[float], List[float]]:
    a = [0.0] * (numNodes + 2)
    b = [0.0] * (numNodes + 1)
    c = [0.0] * (numNodes + 1)
    d = [0.0] * (numNodes + 1)

    # Compute conductances and heat storage
    if airNode < len(thermalConductivity):
        thermalConductance[airNode] = thermalConductivity[airNode]
    for node in range(surfaceNode, numNodes + 1):
        volumeOfSoilAtNode = 0.5 * (nodeDepth[node + 1] - nodeDepth[node - 1])
        heatStorage[node] = Divide(volSpecHeatSoil[node] * volumeOfSoilAtNode, internalTimeStep, 0.0)
        elementLength = nodeDepth[node + 1] - nodeDepth[node]
        thermalConductance[node] = Divide(thermalConductivity[node], elementLength, 0.0)

    g = 1 - nu
    for node in range(surfaceNode, numNodes + 1):
        c[node] = -nu * thermalConductance[node]
        a[node + 1] = c[node]
        b[node] = nu * (thermalConductance[node] + thermalConductance[node - 1]) + heatStorage[node]
        d[node] = g * thermalConductance[node - 1] * soilTemp[node - 1] + ((heatStorage[node] - (g * (thermalConductance[node] + thermalConductance[node - 1]))) * soilTemp[node]) + (g * thermalConductance[node] * soilTemp[node + 1])

    a[surfaceNode] = 0.0
    sensibleHeatFlux = nu * thermalConductance[airNode] * newTemps[airNode]
    radnNet = 0.0
    if netRadiationSource == "calc":
        radnNet = Divide(netRadiation * 1_000_000.0, internalTimeStep, 0.0)
    elif netRadiationSource == "eos":
        radnNet = Divide(waterBalance_Eos * latentHeatOfVapourisation, timestep, 0.0)
    latentHeatFlux = Divide(waterBalance_Es * latentHeatOfVapourisation, timestep, 0.0)
    soilSurfaceHeatFlux = sensibleHeatFlux + radnNet - latentHeatFlux
    d[surfaceNode] = d[surfaceNode] + soilSurfaceHeatFlux
    d[numNodes] = d[numNodes] + (nu * thermalConductance[numNodes] * newTemps[numNodes + 1])

    for node in range(surfaceNode, (numNodes - 1) + 1):
        c[node] = Divide(c[node], b[node], 0.0)
        d[node] = Divide(d[node], b[node], 0.0)
        b[node + 1] = b[node + 1] - (a[node + 1] * c[node])
        d[node + 1] = d[node + 1] - (a[node + 1] * d[node])

    newTemps[numNodes] = Divide(d[numNodes], b[numNodes], 0.0)
    for node in range(numNodes - 1, surfaceNode - 1, -1):
        newTemps[node] = d[node] - (c[node] * newTemps[node + 1])

    return newTemps, heatStorage, thermalConductance


def Calculate_getBoundaryLayerConductance(
    TNew_zb: List[float],
    weather_AirPressure: float,
    stefanBoltzmannConstant: float,
    waterBalance_Eos: float,
    weather_Wind: float,
    airTemperature: float,
    surfaceNode: int,
    waterBalance_Eo: float,
    instrumentHeight: float,
    canopyHeight: float,
) -> Tuple[float, List[float]]:
    vonKarmanConstant = 0.41
    gravitationalConstant = 9.8
    specificHeatOfAir = 1010.0
    surfaceEmissivity = 0.98
    SpecificHeatAir = specificHeatOfAir * airDensity(airTemperature, weather_AirPressure)
    roughnessFactorMomentum = 0.13 * canopyHeight
    roughnessFactorHeat = 0.2 * roughnessFactorMomentum
    d = 0.77 * canopyHeight
    surfaceTemperature = TNew_zb[surfaceNode]
    diffusePenetrationConstant = max(0.1, waterBalance_Eos) / max(0.1, waterBalance_Eo)
    radiativeConductance = 4.0 * stefanBoltzmannConstant * surfaceEmissivity * diffusePenetrationConstant * (kelvinT(airTemperature) ** 3)
    frictionVelocity = 0.0
    boundaryLayerCond = 0.0
    stabilityParammeter = 0.0
    stabilityCorrectionMomentum = 0.0
    stabilityCorrectionHeat = 0.0
    heatFluxDensity = 0.0
    for _ in range(1, 3 + 1):
        frictionVelocity = Divide(
            weather_Wind * vonKarmanConstant,
            (math_log_safe(Divide(instrumentHeight - d + roughnessFactorMomentum, roughnessFactorMomentum, 0.0)) + stabilityCorrectionMomentum),
            0.0,
        )
        boundaryLayerCond = Divide(
            SpecificHeatAir * vonKarmanConstant * frictionVelocity,
            (math_log_safe(Divide(instrumentHeight - d + roughnessFactorHeat, roughnessFactorHeat, 0.0)) + stabilityCorrectionHeat),
            0.0,
        )
        boundaryLayerCond = boundaryLayerCond + radiativeConductance
        heatFluxDensity = boundaryLayerCond * (surfaceTemperature - airTemperature)
        stabilityParammeter = Divide(-vonKarmanConstant * instrumentHeight * gravitationalConstant * heatFluxDensity, SpecificHeatAir * kelvinT(airTemperature) * (frictionVelocity ** 3.0), 0.0)
        if stabilityParammeter > 0.0:
            stabilityCorrectionHeat = 4.7 * stabilityParammeter
            stabilityCorrectionMomentum = stabilityCorrectionHeat
        else:
            stabilityCorrectionHeat = -2.0 * math_log_safe((1.0 + (1.0 - (16.0 * stabilityParammeter)) ** 0.5) / 2.0)
            stabilityCorrectionMomentum = 0.6 * stabilityCorrectionHeat
    return boundaryLayerCond, TNew_zb


def math_log_safe(x: float) -> float:
    # Safe natural log
    import math
    if x <= 0.0:
        x = 1e-6
    return math.log(x)


def interpolateNetRadiation(
    solarRadn: float,
    cloudFr: float,
    cva: float,
    waterBalance_Eo: float,
    waterBalance_Eos: float,
    waterBalance_Salb: float,
    soilTemp: List[float],
    airTemperature: float,
    surfaceNode: int,
    internalTimeStep: float,
    stefanBoltzmannConstant: float,
) -> float:
    surfaceEmissivity = 0.96
    w2MJ = internalTimeStep / 1_000_000.0
    emissivityAtmos = (1 - (0.84 * cloudFr)) * 0.58 * (cva ** (1.0 / 7.0)) + (0.84 * cloudFr)
    PenetrationConstant = Divide(max(0.1, waterBalance_Eos), max(0.1, waterBalance_Eo), 0.0)
    lwRinSoil = longWaveRadn(emissivityAtmos, airTemperature, stefanBoltzmannConstant) * PenetrationConstant * w2MJ
    lwRoutSoil = longWaveRadn(surfaceEmissivity, soilTemp[surfaceNode], stefanBoltzmannConstant) * PenetrationConstant * w2MJ
    lwRnetSoil = lwRinSoil - lwRoutSoil
    swRin = solarRadn
    swRout = waterBalance_Salb * solarRadn
    swRnetSoil = (swRin - swRout) * PenetrationConstant
    return swRnetSoil + lwRnetSoil


def interpolateTemperature(
    timeHours: float,
    minTempYesterday: float,
    maxTempYesterday: float,
    weather_MeanT: float,
    weather_MaxT: float,
    weather_MinT: float,
    defaultTimeOfMaximumTemperature: float,
) -> float:
    time = timeHours / 24.0
    maxT_time = defaultTimeOfMaximumTemperature / 24.0
    minT_time = maxT_time - 0.5
    if time < minT_time:
        midnightT = math_sin((0.0 + 0.25 - maxT_time) * 2.0 * 3.141592653589793) * (maxTempYesterday - minTempYesterday) / 2.0 + ((maxTempYesterday + minTempYesterday) / 2.0)
        tScale = (minT_time - time) / minT_time
        if tScale > 1.0:
            tScale = 1.0
        elif tScale < 0.0:
            tScale = 0.0
        currentTemperature = weather_MinT + (tScale * (midnightT - weather_MinT))
        return currentTemperature
    else:
        currentTemperature = math_sin((time + 0.25 - maxT_time) * 2.0 * 3.141592653589793) * (weather_MaxT - weather_MinT) / 2.0 + weather_MeanT
        return currentTemperature


def math_sin(x: float) -> float:
    import math
    return math.sin(x)


def doThermalConductivity(
    soilConstituentNames: List[str],
    numNodes: int,
    soilWater: List[float],
    thermalConductivity: Optional[List[float]],
    carbon: List[float],
    bulkDensity: List[float],
    pom: float,
    rocks: List[float],
    sand: List[float],
    ps: float,
    silt: List[float],
    clay: List[float],
    nodeDepth: List[float],
    thickness: List[float],
    surfaceNode: int,
    MissingValue: float,
) -> List[float]:
    thermCondLayers = [0.0] * (numNodes + 1)
    for node in range(1, numNodes + 1):
        numerator = 0.0
        denominator = 0.0
        for constituentName in soilConstituentNames:
            shapeFactorConstituent = shapeFactor(constituentName, node, soilWater, carbon, bulkDensity, pom, rocks, sand, ps, silt, clay)
            thermalConductanceConstituent = ThermalConductance(constituentName, node, rocks, bulkDensity, sand, ps, carbon, pom, silt, clay)
            thermalConductanceWater = ThermalConductance("Water", node, rocks, bulkDensity, sand, ps, carbon, pom, silt, clay)
            k = 2.0 / 3.0 * ((1 + (shapeFactorConstituent * (thermalConductanceConstituent / thermalConductanceWater - 1.0))) ** -1) + (1.0 / 3.0 * ((1 + (shapeFactorConstituent * (thermalConductanceConstituent / thermalConductanceWater - 1.0) * (1 - (2 * shapeFactorConstituent)))) ** -1))
            numerator = numerator + (thermalConductanceConstituent * soilWater[node] * k)
            denominator = denominator + (soilWater[node] * k)
        thermCondLayers[node] = numerator / denominator if denominator != 0.0 else 0.0
    thermalConductivity = mapLayer2Node(thermCondLayers, thermalConductivity, nodeDepth, numNodes, thickness, surfaceNode, MissingValue)
    return thermalConductivity


def doVolumetricSpecificHeat(
    soilConstituentNames: List[str],
    numNodes: int,
    volSpecHeatSoil: Optional[List[float]],
    soilWater: List[float],
    nodeDepth: List[float],
    thickness: List[float],
    surfaceNode: int,
    MissingValue: float,
) -> List[float]:
    volspecHeatSoil_ = [0.0] * (numNodes + 1)
    for node in range(1, numNodes + 1):
        volspecHeatSoil_[node] = 0.0
        for constituentName in soilConstituentNames:
            if constituentName not in ["Minerals"]:
                volspecHeatSoil_[node] = volspecHeatSoil_[node] + (volumetricSpecificHeat(constituentName, node) * 1_000_000.0 * soilWater[node])
    volSpecHeatSoil = mapLayer2Node(volspecHeatSoil_, volSpecHeatSoil, nodeDepth, numNodes, thickness, surfaceNode, MissingValue)
    return volSpecHeatSoil


def Zero(arr: Optional[List[float]]) -> List[float]:
    if arr is None:
        return []
    return [0.0 for _ in arr]


def Calculate_doNetRadiation(
    solarRadn: Optional[List[float]],
    cloudFr: float,
    cva: float,
    ITERATIONSperDAY: int,
    weather_MinT: float,
    clock_Today_DayOfYear: int,
    weather_Radn: float,
    weather_Latitude: float,
) -> Tuple[List[float], float, float]:
    TSTEPS2RAD = Divide(2.0 * 3.141592653589793, float(ITERATIONSperDAY), 0.0)
    solarConstant = 1360.0
    import math
    solarDeclination = 0.3985 * math.sin((4.869 + (clock_Today_DayOfYear * 2.0 * 3.141592653589793 / 365.25) + (0.03345 * math.sin((6.224 + (clock_Today_DayOfYear * 2.0 * 3.141592653589793 / 365.25))))))
    cD = (1.0 - (solarDeclination * solarDeclination)) ** 0.5
    solarRadn = ([0.0] * (ITERATIONSperDAY + 1)) if (solarRadn is None or len(solarRadn) < (ITERATIONSperDAY + 1)) else list(solarRadn[: (ITERATIONSperDAY + 1)])
    m1 = [0.0] * (ITERATIONSperDAY + 1)
    m1Tot = 0.0
    for timestepNumber in range(1, ITERATIONSperDAY + 1):
        m1[timestepNumber] = (solarDeclination * math.sin(weather_Latitude * 3.141592653589793 / 180.0) + (cD * math.cos(weather_Latitude * 3.141592653589793 / 180.0) * math.cos(TSTEPS2RAD * (timestepNumber - (ITERATIONSperDAY / 2.0))))) * 24.0 / ITERATIONSperDAY
        if m1[timestepNumber] > 0.0:
            m1Tot = m1Tot + m1[timestepNumber]
        else:
            m1[timestepNumber] = 0.0
    psr = m1Tot * solarConstant * 3600.0 / 1_000_000.0
    fr = Divide(max(weather_Radn, 0.1), psr, 0.0)
    cloudFr = 2.33 - (3.33 * fr)
    cloudFr = min(max(cloudFr, 0.0), 1.0)
    for timestepNumber in range(1, ITERATIONSperDAY + 1):
        solarRadn[timestepNumber] = max(weather_Radn, 0.1) * Divide(m1[timestepNumber], m1Tot, 0.0)
    cva = math.exp((31.3716 - (6014.79 / kelvinT(weather_MinT)) - (0.00792495 * kelvinT(weather_MinT)))) / kelvinT(weather_MinT)
    return solarRadn, cloudFr, cva


def Calculate_doProcess(
    timeOfDaySecs: float,
    netRadiation: float,
    minSoilTemp: List[float],
    maxSoilTemp: List[float],
    numIterationsForBoundaryLayerConductance: int,
    timestep: float,
    boundaryLayerConductance: float,
    maxTempYesterday: float,
    airNode: int,
    soilTemp: List[float],
    airTemperature: float,
    newTemperature: List[float],
    weather_MaxT: float,
    internalTimeStep: float,
    boundarLayerConductanceSource: str,
    thermalConductivity: List[float],
    minTempYesterday: float,
    aveSoilTemp: List[float],
    morningSoilTemp: List[float],
    weather_MeanT: float,
    constantBoundaryLayerConductance: float,
    weather_MinT: float,
    clock_Today_DayOfYear: int,
    weather_Radn: float,
    weather_Latitude: float,
    soilConstituentNames: List[str],
    numNodes: int,
    volSpecHeatSoil: List[float],
    soilWater: List[float],
    nodeDepth: List[float],
    thickness: List[float],
    surfaceNode: int,
    MissingValue: float,
    carbon: List[float],
    bulkDensity: List[float],
    pom: float,
    rocks: List[float],
    sand: List[float],
    ps: float,
    silt: List[float],
    clay: List[float],
    defaultTimeOfMaximumTemperature: float,
    waterBalance_Eo: float,
    waterBalance_Eos: float,
    waterBalance_Salb: float,
    stefanBoltzmannConstant: float,
    weather_AirPressure: float,
    weather_Wind: float,
    instrumentHeight: float,
    canopyHeight: float,
    heatStorage: List[float],
    netRadiationSource: str,
    latentHeatOfVapourisation: float,
    waterBalance_Es: float,
    thermalConductance: List[float],
    nu: float,
) -> Tuple[
    List[float],  # minSoilTemp
    List[float],  # maxSoilTemp
    List[float],  # soilTemp
    List[float],  # newTemperature
    List[float],  # thermalConductivity
    List[float],  # aveSoilTemp
    List[float],  # morningSoilTemp
    List[float],  # volSpecHeatSoil
    List[float],  # heatStorage
    List[float],  # thermalConductance
    float,        # timeOfDaySecs
    float,        # netRadiation
    float,        # airTemperature
    float,        # internalTimeStep
    float,        # minTempYesterday
    float,        # boundaryLayerConductance
    float,        # maxTempYesterday
]:
    interactionsPerDay = 48
    cva = 0.0
    cloudFr = 0.0
    solarRadn = [0.0] * (interactionsPerDay + 1)
    solarRadn, cloudFr, cva = Calculate_doNetRadiation(solarRadn, cloudFr, cva, interactionsPerDay, weather_MinT, clock_Today_DayOfYear, weather_Radn, weather_Latitude)

    minSoilTemp = Zero(minSoilTemp)
    maxSoilTemp = Zero(maxSoilTemp)
    aveSoilTemp = Zero(aveSoilTemp)
    boundaryLayerConductance = 0.0
    internalTimeStep = round(timestep / interactionsPerDay)

    volSpecHeatSoil = doVolumetricSpecificHeat(soilConstituentNames, numNodes, volSpecHeatSoil, soilWater, nodeDepth, thickness, surfaceNode, MissingValue)
    thermalConductivity = doThermalConductivity(soilConstituentNames, numNodes, soilWater, thermalConductivity, carbon, bulkDensity, pom, rocks, sand, ps, silt, clay, nodeDepth, thickness, surfaceNode, MissingValue)

    for timeStepIteration in range(1, interactionsPerDay + 1):
        timeOfDaySecs = internalTimeStep * float(timeStepIteration)
        if timestep < (24.0 * 60.0 * 60.0):
            airTemperature = weather_MeanT
        else:
            airTemperature = interpolateTemperature(timeOfDaySecs / 3600.0, minTempYesterday, maxTempYesterday, weather_MeanT, weather_MaxT, weather_MinT, defaultTimeOfMaximumTemperature)

        newTemperature[airNode] = airTemperature
        netRadiation = interpolateNetRadiation(solarRadn[timeStepIteration], cloudFr, cva, waterBalance_Eo, waterBalance_Eos, waterBalance_Salb, soilTemp, airTemperature, surfaceNode, internalTimeStep, stefanBoltzmannConstant)

        if boundarLayerConductanceSource == "constant":
            if airNode < len(thermalConductivity):
                thermalConductivity[airNode] = constantBoundaryLayerConductance
        elif boundarLayerConductanceSource == "calc":
            boundaryLayerCond, newTemperature = Calculate_getBoundaryLayerConductance(newTemperature, weather_AirPressure, stefanBoltzmannConstant, waterBalance_Eos, weather_Wind, airTemperature, surfaceNode, waterBalance_Eo, instrumentHeight, canopyHeight)
            if airNode < len(thermalConductivity):
                thermalConductivity[airNode] = boundaryLayerCond
            for _ in range(1, numIterationsForBoundaryLayerConductance + 1):
                newTemperature, heatStorage, thermalConductance = Calculate_doThomas(newTemperature, netRadiation, heatStorage, waterBalance_Eos, numNodes, timestep, netRadiationSource, latentHeatOfVapourisation, nodeDepth, waterBalance_Es, airNode, soilTemp, surfaceNode, internalTimeStep, thermalConductance, thermalConductivity, nu, volSpecHeatSoil)
                boundaryLayerCond, newTemperature = Calculate_getBoundaryLayerConductance(newTemperature, weather_AirPressure, stefanBoltzmannConstant, waterBalance_Eos, weather_Wind, airTemperature, surfaceNode, waterBalance_Eo, instrumentHeight, canopyHeight)
                if airNode < len(thermalConductivity):
                    thermalConductivity[airNode] = boundaryLayerCond

        newTemperature, heatStorage, thermalConductance = Calculate_doThomas(newTemperature, netRadiation, heatStorage, waterBalance_Eos, numNodes, timestep, netRadiationSource, latentHeatOfVapourisation, nodeDepth, waterBalance_Es, airNode, soilTemp, surfaceNode, internalTimeStep, thermalConductance, thermalConductivity, nu, volSpecHeatSoil)
        minSoilTemp, soilTemp, maxSoilTemp, aveSoilTemp, boundaryLayerConductance = Calculate_doUpdate(interactionsPerDay, timeOfDaySecs, boundaryLayerConductance, minSoilTemp, airNode, soilTemp, newTemperature, numNodes, surfaceNode, internalTimeStep, maxSoilTemp, aveSoilTemp, thermalConductivity)

        if abs(timeOfDaySecs - (5.0 * 3600.0)) <= (min(timeOfDaySecs, 5.0 * 3600.0) * 0.0001):
            for i in range(len(soilTemp)):
                if i < len(morningSoilTemp):
                    morningSoilTemp[i] = soilTemp[i]

    minTempYesterday = weather_MinT
    maxTempYesterday = weather_MaxT
    return (
        minSoilTemp,
        maxSoilTemp,
        soilTemp,
        newTemperature,
        thermalConductivity,
        aveSoilTemp,
        morningSoilTemp,
        volSpecHeatSoil,
        heatStorage,
        thermalConductance,
        timeOfDaySecs,
        netRadiation,
        airTemperature,
        internalTimeStep,
        minTempYesterday,
        boundaryLayerConductance,
        maxTempYesterday,
    )


def ToCumThickness(Thickness: List[float]) -> List[float]:
    CumThickness = [0.0] * len(Thickness)
    if len(Thickness) > 0:
        CumThickness[0] = Thickness[0]
        for Layer in range(1, len(Thickness)):
            CumThickness[Layer] = Thickness[Layer] + CumThickness[Layer - 1]
    return CumThickness


def calcSoilTemperature(
    soilTempIO: Optional[List[float]],
    weather_Tav: float,
    clock_Today_DayOfYear: int,
    surfaceNode: int,
    numNodes: int,
    weather_Amp: float,
    thickness: List[float],
    weather_Latitude: float,
) -> List[float]:
    soilTempIO = list(soilTempIO) if soilTempIO is not None else [0.0] * (numNodes + 1 + 1)
    soilTemp = [0.0] * (numNodes + 1 + 1)
    cumulativeDepth = ToCumThickness(thickness)
    w = 2 * 3.141592653589793 / (365.25 * 24 * 3600)
    dh = 0.6
    zd = (2 * dh / w) ** 0.5
    offset = 0.25
    if weather_Latitude > 0.0:
        offset = -0.25
    for nodes in range(1, numNodes + 1):
        cd = cumulativeDepth[nodes] if nodes < len(cumulativeDepth) else cumulativeDepth[-1] if cumulativeDepth else 0.0
        soilTemp[nodes] = weather_Tav + (weather_Amp * math_exp_safe(-1 * cd / zd) * math_sin(((clock_Today_DayOfYear / 365.0 + offset) * 2.0 * 3.141592653589793 - (cd / zd))))
    for i in range(len(soilTemp)):
        if i < len(soilTempIO):
            soilTempIO[i] = soilTemp[i]
    return soilTempIO


def math_exp_safe(x: float) -> float:
    import math
    try:
        return math.exp(x)
    except OverflowError:
        return 0.0 if x < 0 else float("inf")


def calcSurfaceTemperature(weather_MeanT: float, weather_MaxT: float, waterBalance_Salb: float, weather_Radn: float) -> float:
    surfaceT = (1.0 - waterBalance_Salb) * (weather_MeanT + ((weather_MaxT - weather_MeanT) * (max(weather_Radn, 0.1) * 23.8846 / 800.0) ** 0.5)) + (waterBalance_Salb * weather_MeanT)
    return surfaceT


def ValuesInArray(Values: Optional[List[float]], MissingValue: float) -> bool:
    if Values is not None:
        for Value in Values:
            if Value != MissingValue and not (Value != Value):  # NaN check: Value != Value is True if NaN
                return True
    return False


def init(
    weather_Tav: float,
    weather_Amp: float,
    weather_Latitude: float,
    clock_Today_DayOfYear: int,
    microClimate_CanopyHeight: float,
    physical_Thickness: List[float],
    physical_BD: List[float],
    physical_Rocks: List[float],
    physical_ParticleSizeSand: List[float],
    physical_ParticleSizeSilt: List[float],
    physical_ParticleSizeClay: List[float],
    organic_Carbon: List[float],
    waterBalance_SW: Optional[List[float]],
    pInitialValues: List[float],
    DepthToConstantTemperature: float,
    timestep: float,
    latentHeatOfVapourisation: float,
    stefanBoltzmannConstant: float,
    airNode: int,
    surfaceNode: int,
    topsoilNode: int,
    numPhantomNodes: int,
    constantBoundaryLayerConductance: float,
    numIterationsForBoundaryLayerConductance: int,
    defaultTimeOfMaximumTemperature: float,
    defaultInstrumentHeight: float,
    bareSoilRoughness: float,
    instrumHeight: float,
    pom: float,
    nu: float,
    boundarLayerConductanceSource: str,
    netRadiationSource: str,
    MissingValue: float,
    soilConstituentNames: List[str],
) -> Tuple[
    bool,               # doInitialisationStuff
    float,              # internalTimeStep
    float,              # timeOfDaySecs
    int,                # numNodes
    int,                # numLayers
    List[float],        # nodeDepth
    List[float],        # thermCondPar1
    List[float],        # thermCondPar2
    List[float],        # thermCondPar3
    List[float],        # thermCondPar4
    List[float],        # volSpecHeatSoil
    List[float],        # soilTemp
    List[float],        # morningSoilTemp
    List[float],        # heatStorage
    List[float],        # thermalConductance
    List[float],        # thermalConductivity
    float,              # boundaryLayerConductance
    List[float],        # newTemperature
    float,              # airTemperature
    float,              # maxTempYesterday
    float,              # minTempYesterday
    List[float],        # soilWater
    List[float],        # minSoilTemp
    List[float],        # maxSoilTemp
    List[float],        # aveSoilTemp
    List[float],        # thickness
    List[float],        # bulkDensity
    List[float],        # rocks
    List[float],        # carbon
    List[float],        # sand
    List[float],        # silt
    List[float],        # clay
    float,              # instrumentHeight
    float,              # netRadiation
    float,              # canopyHeight
    List[float],        # InitialValues
]:
    doInitialisationStuff = True
    instrumentHeight = getIniVariables(0.0, instrumHeight, defaultInstrumentHeight)
    # Allocate placeholders
    heatStorage = None
    minSoilTemp = None
    bulkDensity = None
    maxSoilTemp = None
    nodeDepth = None
    newTemperature = None
    soilWater = None
    thermalConductance = None
    thermalConductivity = None
    sand = None
    carbon = None
    thickness = None
    rocks = None
    clay = None
    soilTemp = None
    silt = None
    volSpecHeatSoil = None
    aveSoilTemp = None
    morningSoilTemp = None
    numNodes = 0
    numLayers = 0

    (
        heatStorage,
        minSoilTemp,
        bulkDensity,
        maxSoilTemp,
        nodeDepth,
        newTemperature,
        soilWater,
        thermalConductance,
        thermalConductivity,
        sand,
        carbon,
        thickness,
        rocks,
        clay,
        soilTemp,
        silt,
        volSpecHeatSoil,
        aveSoilTemp,
        morningSoilTemp,
        numNodes,
        numLayers,
    ) = Calculate_getProfileVariables(
        heatStorage,
        minSoilTemp,
        bulkDensity,
        numNodes,
        physical_BD,
        maxSoilTemp,
        waterBalance_SW if waterBalance_SW is not None else [],
        organic_Carbon,
        physical_Rocks,
        nodeDepth,
        topsoilNode,
        newTemperature,
        surfaceNode,
        soilWater,
        thermalConductance,
        thermalConductivity,
        sand,
        carbon,
        thickness,
        numPhantomNodes,
        physical_ParticleSizeSand,
        rocks,
        clay,
        physical_ParticleSizeSilt,
        airNode,
        physical_ParticleSizeClay,
        soilTemp,
        numLayers,
        physical_Thickness,
        silt,
        volSpecHeatSoil,
        aveSoilTemp,
        morningSoilTemp,
        DepthToConstantTemperature,
        MissingValue,
    )

    thermCondPar1: Optional[List[float]] = None
    thermCondPar2: Optional[List[float]] = None
    thermCondPar3: Optional[List[float]] = None
    thermCondPar4: Optional[List[float]] = None
    (
        newTemperature,
        soilTemp,
        thermCondPar2,
        thermCondPar3,
        thermCondPar4,
        thermCondPar1,
        soilRoughnessHeight,
    ) = Calculate_readParam(
        bareSoilRoughness,
        newTemperature,
        0.0,
        soilTemp,
        thermCondPar2,
        numLayers,
        bulkDensity,
        numNodes,
        thermCondPar3,
        thermCondPar4,
        clay,
        thermCondPar1,
        weather_Tav,
        clock_Today_DayOfYear,
        surfaceNode,
        weather_Amp,
        thickness,
        weather_Latitude,
    )

    InitialValues = list(pInitialValues)

    internalTimeStep = 0.0
    timeOfDaySecs = 0.0
    boundaryLayerConductance = 0.0
    airTemperature = 0.0
    maxTempYesterday = 0.0
    minTempYesterday = 0.0
    netRadiation = 0.0
    canopyHeight = 0.0

    return (
        doInitialisationStuff,
        internalTimeStep,
        timeOfDaySecs,
        numNodes,
        numLayers,
        nodeDepth,
        thermCondPar1 if thermCondPar1 is not None else [0.0] * (numNodes + 2),
        thermCondPar2 if thermCondPar2 is not None else [0.0] * (numNodes + 2),
        thermCondPar3 if thermCondPar3 is not None else [0.0] * (numNodes + 2),
        thermCondPar4 if thermCondPar4 is not None else [0.0] * (numNodes + 2),
        volSpecHeatSoil,
        soilTemp,
        morningSoilTemp,
        heatStorage,
        thermalConductance,
        thermalConductivity,
        boundaryLayerConductance,
        newTemperature,
        airTemperature,
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
        instrumentHeight,
        netRadiation,
        canopyHeight,
        InitialValues,
    )


def process(
    weather_MinT: float,
    weather_MaxT: float,
    weather_MeanT: float,
    weather_Tav: float,
    weather_Amp: float,
    weather_AirPressure: float,
    weather_Wind: float,
    weather_Latitude: float,
    weather_Radn: float,
    clock_Today_DayOfYear: int,
    microClimate_CanopyHeight: float,
    waterBalance_SW: List[float],
    waterBalance_Eos: float,
    waterBalance_Eo: float,
    waterBalance_Es: float,
    waterBalance_Salb: float,
    InitialValues: List[float],
    DepthToConstantTemperature: float,  # not used in process; kept for signature clarity if needed
    timestep: float,
    latentHeatOfVapourisation: float,
    stefanBoltzmannConstant: float,
    airNode: int,
    surfaceNode: int,
    topsoilNode: int,
    numPhantomNodes: int,  # not used directly in process
    constantBoundaryLayerConductance: float,
    numIterationsForBoundaryLayerConductance: int,
    defaultTimeOfMaximumTemperature: float,
    bareSoilRoughness: float,  # not used directly in process
    doInitialisationStuff: bool,
    internalTimeStep: float,
    timeOfDaySecs: float,
    numNodes: int,
    numLayers: int,
    nodeDepth: List[float],
    volSpecHeatSoil: List[float],
    soilTemp: List[float],
    morningSoilTemp: List[float],
    heatStorage: List[float],
    thermalConductance: List[float],
    thermalConductivity: List[float],
    boundaryLayerConductance: float,
    newTemperature: List[float],
    airTemperature: float,
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
    pom: float,
    silt: List[float],
    clay: List[float],
    soilRoughnessHeight: float,
    instrumentHeight: float,
    netRadiation: float,
    canopyHeight: float,
    nu: float,
    boundarLayerConductanceSource: str,
    netRadiationSource: str,
    MissingValue: float,
    soilConstituentNames: List[str],
) -> Tuple[
    List[float],  # InitialValues
    bool,         # doInitialisationStuff
    float,        # internalTimeStep
    float,        # timeOfDaySecs
    List[float],  # volSpecHeatSoil
    List[float],  # soilTemp
    List[float],  # morningSoilTemp
    List[float],  # heatStorage
    List[float],  # thermalConductance
    List[float],  # thermalConductivity
    float,        # boundaryLayerConductance
    List[float],  # newTemperature
    float,        # airTemperature
    float,        # maxTempYesterday
    float,        # minTempYesterday
    List[float],  # soilWater
    List[float],  # minSoilTemp
    List[float],  # maxSoilTemp
    List[float],  # aveSoilTemp
    float,        # instrumentHeight
    float,        # netRadiation
    float,        # canopyHeight
]:
    # Update soilWater and heights
    soilWater, instrumentHeight, canopyHeight = Calculate_getOtherVariables(
        numLayers, numNodes, soilWater, instrumentHeight, soilRoughnessHeight, waterBalance_SW, microClimate_CanopyHeight, canopyHeight
    )

    if doInitialisationStuff:
        if ValuesInArray(InitialValues, MissingValue):
            soilTemp = [0.0] * len(soilTemp)
            # Copy InitialValues into soilTemp starting at topsoilNode
            for i in range(len(InitialValues)):
                idx = topsoilNode + i
                if idx < len(soilTemp):
                    soilTemp[idx] = InitialValues[i]
        else:
            soilTemp = calcSoilTemperature(soilTemp, weather_Tav, clock_Today_DayOfYear, surfaceNode, numNodes, weather_Amp, thickness, weather_Latitude)
            # Store initial values (from soilTemp[topsoilNode..])
            InitialValues = [0.0] * len(InitialValues)
            for i in range(numLayers):
                src_idx = topsoilNode + i
                if i < len(InitialValues) and src_idx < len(soilTemp):
                    InitialValues[i] = soilTemp[src_idx]
        if airNode < len(soilTemp):
            soilTemp[airNode] = weather_MeanT
        if surfaceNode < len(soilTemp):
            soilTemp[surfaceNode] = calcSurfaceTemperature(weather_MeanT, weather_MaxT, waterBalance_Salb, weather_Radn)
        for i in range(numNodes + 1, len(soilTemp)):
            soilTemp[i] = weather_Tav
        newTemperature = list(soilTemp)
        maxTempYesterday = weather_MaxT
        minTempYesterday = weather_MinT
        doInitialisationStuff = False

    (
        minSoilTemp,
        maxSoilTemp,
        soilTemp,
        newTemperature,
        thermalConductivity,
        aveSoilTemp,
        morningSoilTemp,
        volSpecHeatSoil,
        heatStorage,
        thermalConductance,
        timeOfDaySecs,
        netRadiation,
        airTemperature,
        internalTimeStep,
        minTempYesterday,
        boundaryLayerConductance,
        maxTempYesterday,
    ) = Calculate_doProcess(
        timeOfDaySecs,
        netRadiation,
        minSoilTemp,
        maxSoilTemp,
        numIterationsForBoundaryLayerConductance,
        timestep,
        boundaryLayerConductance,
        maxTempYesterday,
        airNode,
        soilTemp,
        airTemperature,
        newTemperature,
        weather_MaxT,
        internalTimeStep,
        boundarLayerConductanceSource,
        thermalConductivity,
        minTempYesterday,
        aveSoilTemp,
        morningSoilTemp,
        weather_MeanT,
        constantBoundaryLayerConductance,
        weather_MinT,
        clock_Today_DayOfYear,
        weather_Radn,
        weather_Latitude,
        soilConstituentNames,
        numNodes,
        volSpecHeatSoil,
        soilWater,
        nodeDepth,
        thickness,
        surfaceNode,
        MissingValue,
        carbon,
        bulkDensity,
        pom,
        rocks,
        sand,
        ps=0.0,  # Not used directly here, placeholder if needed elsewhere
        silt=silt,
        clay=clay,
        defaultTimeOfMaximumTemperature=defaultTimeOfMaximumTemperature,
        waterBalance_Eo=waterBalance_Eo,
        waterBalance_Eos=waterBalance_Eos,
        waterBalance_Salb=waterBalance_Salb,
        stefanBoltzmannConstant=stefanBoltzmannConstant,
        weather_AirPressure=weather_AirPressure,
        weather_Wind=weather_Wind,
        instrumentHeight=instrumentHeight,
        canopyHeight=canopyHeight,
        heatStorage=heatStorage,
        netRadiationSource=netRadiationSource,
        latentHeatOfVapourisation=latentHeatOfVapourisation,
        waterBalance_Es=waterBalance_Es,
        thermalConductance=thermalConductance,
        nu=nu,
    )

    return (
        InitialValues,
        doInitialisationStuff,
        internalTimeStep,
        timeOfDaySecs,
        volSpecHeatSoil,
        soilTemp,
        morningSoilTemp,
        heatStorage,
        thermalConductance,
        thermalConductivity,
        boundaryLayerConductance,
        newTemperature,
        airTemperature,
        maxTempYesterday,
        minTempYesterday,
        soilWater,
        minSoilTemp,
        maxSoilTemp,
        aveSoilTemp,
        instrumentHeight,
        netRadiation,
        canopyHeight,
    )
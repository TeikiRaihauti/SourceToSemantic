from typing import List, Optional, Tuple
import math

def init_soiltemperature(
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
    physical_Thickness: List[float],
    physical_BD: List[float],
    ps: float,
    physical_Rocks: List[float],
    physical_ParticleSizeSand: List[float],
    physical_ParticleSizeSilt: List[float],
    physical_ParticleSizeClay: List[float],
    organic_Carbon: List[float],
    waterBalance_SW: Optional[List[float]],
    waterBalance_Eos: float,
    waterBalance_Eo: float,
    waterBalance_Es: float,
    waterBalance_Salb: float,
    pInitialValues: Optional[List[float]],
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
    nodeDepth: Optional[List[float]],
    thermCondPar1: Optional[List[float]],
    thermCondPar2: Optional[List[float]],
    thermCondPar3: Optional[List[float]],
    thermCondPar4: Optional[List[float]],
    pom: float,
    soilRoughnessHeight: float,
    nu: float,
    boundarLayerConductanceSource: str,
    netRadiationSource: str,
    MissingValue: float,
    soilConstituentNames: List[str],
) -> Tuple[
    List[float], bool, float, float, int, int, List[float], List[float], List[float],
    List[float], List[float], List[float], float, List[float], float, float, List[float],
    List[float], List[float], List[float], List[float], List[float], List[float],
    List[float], List[float], List[float], float, float, float, float,
    List[float], List[float], List[float], List[float], List[float], float
]:
    # Outputs (ordered to match Fortran OUT/INOUT roles; see type hints for datatypes)
    doInitialisationStuff = True
    internalTimeStep = 0.0
    timeOfDaySecs = 0.0
    numNodes = 0
    numLayers = 0
    boundaryLayerConductance = 0.0
    airTemperature = 0.0
    maxTempYesterday = 0.0
    minTempYesterday = 0.0
    instrumentHeight = 0.0
    netRadiation = 0.0
    canopyHeight = 0.0
    instrumHeight = 0.0

    # Set instrument height
    instrumentHeight = getIniVariables(instrumentHeight, instrumHeight, defaultInstrumentHeight)

    # Allocate and derive profile variables
    (
        heatStorage,
        minSoilTemp,
        bulkDensity,
        numNodes,
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
        numLayers,
        silt,
        volSpecHeatSoil,
        aveSoilTemp,
        morningSoilTemp,
    ) = getProfileVariables(
        None, None, None, 0, physical_BD, None, waterBalance_SW, organic_Carbon, physical_Rocks,
        nodeDepth, topsoilNode, None, surfaceNode, None, None, None, None, None, None, numPhantomNodes,
        physical_ParticleSizeSand, None, None, physical_ParticleSizeSilt, airNode,
        physical_ParticleSizeClay, None, 0, physical_Thickness, None, None, None,
        DepthToConstantTemperature, MissingValue
    )

    # Thermal conductivity coeffs and initial temperatures
    thermCondPar1, thermCondPar2, thermCondPar3, thermCondPar4 = readParam(
        bareSoilRoughness,
        [0.0] * len(soilTemp),
        soilRoughnessHeight,
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

    InitialValues = pInitialValues[:] if pInitialValues is not None else []

    return (
        InitialValues,
        doInitialisationStuff,
        internalTimeStep,
        timeOfDaySecs,
        numNodes,
        numLayers,
        volSpecHeatSoil,
        soilTemp,
        morningSoilTemp,
        heatStorage,
        thermalConductance,
        thermalConductivity,
        boundaryLayerConductance,
        [0.0] * len(soilTemp),  # newTemperature
        airTemperature,
        maxTempYesterday,
        minTempYesterday,
        soilWater,
        minSoilTemp,
        maxSoilTemp,
        aveSoilTemp,
        [0.0] * len(thickness),  # aveSoilWater (not used here)
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
        instrumHeight,
        nodeDepth,
        thermCondPar1,
        thermCondPar2,
        thermCondPar3,
        thermCondPar4,
        soilRoughnessHeight,
    )


def model_soiltemperature(
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
    physical_Thickness: List[float],
    physical_BD: List[float],
    ps: float,
    physical_Rocks: List[float],
    physical_ParticleSizeSand: List[float],
    physical_ParticleSizeSilt: List[float],
    physical_ParticleSizeClay: List[float],
    organic_Carbon: List[float],
    waterBalance_SW: List[float],
    waterBalance_Eos: float,
    waterBalance_Eo: float,
    waterBalance_Es: float,
    waterBalance_Salb: float,
    InitialValues: Optional[List[float]],
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
    doInitialisationStuff: bool,
    internalTimeStep: float,
    timeOfDaySecs: float,
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
    airTemperature: float,
    maxTempYesterday: float,
    minTempYesterday: float,
    soilWater: List[float],
    minSoilTemp: List[float],
    maxSoilTemp: List[float],
    aveSoilTemp: List[float],
    aveSoilWater: List[float],
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
    instrumHeight: float,
    nu: float,
    boundarLayerConductanceSource: str,
    netRadiationSource: str,
    MissingValue: float,
    soilConstituentNames: List[str],
) -> Tuple[
    List[float], float, float, List[float], List[float], List[float], List[float],
    List[float], float, List[float], List[float], List[float], List[float], float,
    float, float, float, float, float, List[float], List[float], List[float], List[float],
    bool, float, float, float, float, float
]:
    # Update other vars from inputs
    soilWater, instrumentHeight, canopyHeight = getOtherVariables(
        numLayers, numNodes, soilWater, instrumentHeight, soilRoughnessHeight, waterBalance_SW,
        microClimate_CanopyHeight, canopyHeight
    )

    # Initialization if needed
    if doInitialisationStuff:
        if ValuesInArray(InitialValues, MissingValue):
            soilTemp = [0.0] * (numNodes + 2)
            # Assign provided initial values to topsoil nodes
            for i in range(len(InitialValues)):
                if topsoilNode + i < len(soilTemp):
                    soilTemp[topsoilNode + i] = InitialValues[i]
        else:
            soilTemp = calcSoilTemperature(
                soilTemp[:], weather_Tav, clock_Today_DayOfYear, surfaceNode, numNodes,
                weather_Amp, thickness, weather_Latitude
            )
            InitialValues = [0.0] * numLayers
            for i in range(numLayers):
                InitialValues[i] = soilTemp[topsoilNode + i]
        soilTemp[airNode] = weather_MeanT
        soilTemp[surfaceNode] = calcSurfaceTemperature(weather_MeanT, weather_MaxT, waterBalance_Salb, weather_Radn)
        for i in range(numNodes + 1, len(soilTemp)):
            soilTemp[i] = weather_Tav
        newTemperature = soilTemp[:]
        maxTempYesterday = weather_MaxT
        minTempYesterday = weather_MinT
        doInitialisationStuff = False

    (
        timeOfDaySecs,
        netRadiation,
        minSoilTemp,
        maxSoilTemp,
        boundaryLayerConductance,
        maxTempYesterday,
        soilTemp,
        airTemperature,
        newTemperature,
        internalTimeStep,
        thermalConductivity,
        minTempYesterday,
        aveSoilTemp,
        morningSoilTemp,
        volSpecHeatSoil,
        heatStorage,
        thermalConductance,
    ) = doProcess(
        timeOfDaySecs, netRadiation, minSoilTemp, maxSoilTemp, numIterationsForBoundaryLayerConductance,
        timestep, boundaryLayerConductance, maxTempYesterday, airNode, soilTemp, airTemperature,
        newTemperature, weather_MaxT, internalTimeStep, boundarLayerConductanceSource, thermalConductivity,
        minTempYesterday, aveSoilTemp, morningSoilTemp, weather_MeanT, constantBoundaryLayerConductance,
        weather_MinT, clock_Today_DayOfYear, weather_Radn, weather_Latitude, soilConstituentNames, numNodes,
        volSpecHeatSoil, soilWater, nodeDepth, thickness, surfaceNode, MissingValue, carbon, bulkDensity, pom,
        rocks, sand, ps, silt, clay, defaultTimeOfMaximumTemperature, waterBalance_Eo, waterBalance_Eos,
        waterBalance_Salb, stefanBoltzmannConstant, weather_AirPressure, weather_Wind, instrumentHeight,
        canopyHeight, heatStorage, netRadiationSource, latentHeatOfVapourisation, waterBalance_Es,
        thermalConductance, nu
    )

    return (
        heatStorage,
        instrumentHeight,
        canopyHeight,
        minSoilTemp,
        maxSoilTemp,
        aveSoilTemp,
        volSpecHeatSoil,
        soilTemp,
        netRadiation,
        morningSoilTemp,
        newTemperature,
        thermalConductivity,
        thermalConductance,
        boundaryLayerConductance,
        maxTempYesterday,
        minTempYesterday,
        airTemperature,
        internalTimeStep,
        timeOfDaySecs,
        InitialValues,
        soilWater,
        rocks,
        carbon,
        doInitialisationStuff,
        netRadiation,  # duplicate state output kept for parity
        maxTempYesterday,
        minTempYesterday,
        airTemperature,
        boundaryLayerConductance
    )


def getIniVariables(instrumentHeight: float, instrumHeight: float, defaultInstrumentHeight: float) -> float:
    if instrumHeight > 0.00001:
        instrumentHeight = instrumHeight
    else:
        instrumentHeight = defaultInstrumentHeight
    return instrumentHeight


def getProfileVariables(
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
    List[float], List[float], List[float], int, List[float], List[float], List[float], List[float],
    List[float], List[float], List[float], List[float], List[float], List[float], List[float],
    List[float], int, List[float], List[float], List[float], float, float
]:
    numLayers = len(physical_Thickness)
    numNodes = numLayers + numPhantomNodes

    thickness = [0.0] * (numLayers + numPhantomNodes + 1)
    for i in range(len(physical_Thickness)):
        if 1 + i < len(thickness):
            thickness[1 + i] = physical_Thickness[i]

    belowProfileDepth = max(DepthToConstantTemperature - Sum(thickness, 1, numLayers, MissingValue), 1000.0)
    thicknessForPhantomNodes = belowProfileDepth * 2.0 / float(numPhantomNodes)
    firstPhantomNode = numLayers
    for i in range(firstPhantomNode, firstPhantomNode + numPhantomNodes):
        if i < len(thickness):
            thickness[i] = thicknessForPhantomNodes

    nodeDepth = [0.0] * (numNodes + 2)
    nodeDepth[airNode] = 0.0
    nodeDepth[surfaceNode] = 0.0
    nodeDepth[topsoilNode] = 0.5 * thickness[1] / 1000.0
    for node in range(topsoilNode, numNodes + 1):
        nodeDepth[node + 1] = (Sum(thickness, 1, node - 1, MissingValue) + (0.5 * thickness[node])) / 1000.0

    bulkDensity_out = [0.0] * (numLayers + 1 + numPhantomNodes)
    for i in range(len(physical_BD)):
        if 1 + i < len(bulkDensity_out):
            bulkDensity_out[1 + i] = physical_BD[i]
    if numNodes < len(bulkDensity_out):
        bulkDensity_out[numNodes] = bulkDensity_out[numLayers]
        for layer in range(numLayers + 1, numLayers + numPhantomNodes):
            if layer < len(bulkDensity_out):
                bulkDensity_out[layer] = bulkDensity_out[numLayers]

    soilWater_out = [0.0] * (numLayers + 1 + numPhantomNodes)
    if waterBalance_SW is not None:
        for layer in range(1, numLayers + 1):
            soilWater_out[layer] = Divide(
                waterBalance_SW[layer - 1] * thickness[layer - 1], thickness[layer], 0.0
            )
        for layer in range(numLayers + 1, numLayers + numPhantomNodes):
            if layer < len(soilWater_out):
                soilWater_out[layer] = soilWater_out[numLayers]

    carbon_out = [0.0] * (numLayers + 1 + numPhantomNodes)
    for layer in range(1, numLayers + 1):
        carbon_out[layer] = organic_Carbon[layer - 1]
    for layer in range(numLayers + 1, numLayers + numPhantomNodes):
        carbon_out[layer] = carbon_out[numLayers]

    rocks_out = [0.0] * (numLayers + 1 + numPhantomNodes)
    for layer in range(1, numLayers + 1):
        rocks_out[layer] = physical_Rocks[layer - 1]
    for layer in range(numLayers + 1, numLayers + numPhantomNodes):
        rocks_out[layer] = rocks_out[numLayers]

    sand_out = [0.0] * (numLayers + 1 + numPhantomNodes)
    for layer in range(1, numLayers + 1):
        sand_out[layer] = physical_ParticleSizeSand[layer - 1]
    for layer in range(numLayers + 1, numLayers + numPhantomNodes):
        sand_out[layer] = sand_out[numLayers]

    silt_out = [0.0] * (numLayers + 1 + numPhantomNodes)
    for layer in range(1, numLayers + 1):
        silt_out[layer] = physical_ParticleSizeSilt[layer - 1]
    for layer in range(numLayers + 1, numLayers + numPhantomNodes):
        silt_out[layer] = silt_out[numLayers]

    clay_out = [0.0] * (numLayers + 1 + numPhantomNodes)
    for layer in range(1, numLayers + 1):
        clay_out[layer] = physical_ParticleSizeClay[layer - 1]
    for layer in range(numLayers + 1, numLayers + numPhantomNodes):
        clay_out[layer] = clay_out[numLayers]

    maxSoilTemp_out = [0.0] * (numLayers + 1 + numPhantomNodes)
    minSoilTemp_out = [0.0] * (numLayers + 1 + numPhantomNodes)
    aveSoilTemp_out = [0.0] * (numLayers + 1 + numPhantomNodes)
    volSpecHeatSoil_out = [0.0] * (numNodes + 1)
    soilTemp_out = [0.0] * (numNodes + 2)
    morningSoilTemp_out = [0.0] * (numNodes + 2)
    newTemperature_out = [0.0] * (numNodes + 2)
    thermalConductivity_out = [0.0] * (numNodes + 1)
    heatStorage_out = [0.0] * (numNodes + 1)
    thermalConductance_out = [0.0] * (numNodes + 2)

    return (
        heatStorage_out,
        minSoilTemp_out,
        bulkDensity_out,
        numNodes,
        maxSoilTemp_out,
        nodeDepth,
        newTemperature_out,
        soilWater_out,
        thermalConductance_out,
        thermalConductivity_out,
        sand_out,
        carbon_out,
        thickness,
        rocks_out,
        clay_out,
        soilTemp_out,
        numLayers,
        silt_out,
        volSpecHeatSoil_out,
        aveSoilTemp_out,
        morningSoilTemp_out,
        DepthToConstantTemperature,
        MissingValue,
    )


def doThermalConductivityCoeffs(
    thermCondPar2: Optional[List[float]],
    numLayers: int,
    bulkDensity: List[float],
    numNodes: int,
    thermCondPar3: Optional[List[float]],
    thermCondPar4: Optional[List[float]],
    clay: List[float],
    thermCondPar1: Optional[List[float]],
) -> Tuple[List[float], List[float], List[float], List[float]]:
    thermCondPar1_out = [0.0] * (numNodes + 1)
    thermCondPar2_out = [0.0] * (numNodes + 1)
    thermCondPar3_out = [0.0] * (numNodes + 1)
    thermCondPar4_out = [0.0] * (numNodes + 1)
    for layer in range(1, numLayers + 2):
        element = layer
        thermCondPar1_out[element] = 0.65 - (0.78 * bulkDensity[layer]) + (0.6 * (bulkDensity[layer] ** 2))
        thermCondPar2_out[element] = 1.06 * bulkDensity[layer]
        thermCondPar3_out[element] = 1.0 + Divide(2.6, math.sqrt(clay[layer]) if clay[layer] > 0.0 else 0.0, 0.0)
        thermCondPar4_out[element] = 0.03 + (0.1 * (bulkDensity[layer] ** 2))
    return thermCondPar1_out, thermCondPar2_out, thermCondPar3_out, thermCondPar4_out


def readParam(
    bareSoilRoughness: float,
    newTemperature: List[float],
    soilRoughnessHeight: float,
    soilTemp: List[float],
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
) -> Tuple[List[float], List[float], List[float], List[float]]:
    thermCondPar1_out, thermCondPar2_out, thermCondPar3_out, thermCondPar4_out = doThermalConductivityCoeffs(
        thermCondPar2, numLayers, bulkDensity, numNodes, thermCondPar3, thermCondPar4, clay, thermCondPar1
    )
    soilTemp_new = calcSoilTemperature(
        soilTemp[:], weather_Tav, clock_Today_DayOfYear, surfaceNode, numNodes, weather_Amp, thickness,
        weather_Latitude
    )
    for i in range(len(soilTemp_new)):
        if i < len(newTemperature):
            newTemperature[i] = soilTemp_new[i]
    soilRoughnessHeight = bareSoilRoughness
    return thermCondPar1_out, thermCondPar2_out, thermCondPar3_out, thermCondPar4_out


def getOtherVariables(
    numLayers: int,
    numNodes: int,
    soilWater: List[float],
    instrumentHeight: float,
    soilRoughnessHeight: float,
    waterBalance_SW: List[float],
    microClimate_CanopyHeight: float,
    canopyHeight: float
) -> Tuple[List[float], float, float]:
    for i in range(len(waterBalance_SW)):
        if 1 + i < len(soilWater):
            soilWater[1 + i] = waterBalance_SW[i]
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
    return (1.0 - volumetricFractionOrganicMatter(layer, carbon, bulkDensity, pom)) * soilWater[layer]


def volumetricFractionClay(
    layer: int, bulkDensity: List[float], ps: float, clay: List[float], carbon: List[float], pom: float, rocks: List[float]
) -> float:
    return (1.0 - volumetricFractionOrganicMatter(layer, carbon, bulkDensity, pom) - volumetricFractionRocks(layer, rocks)) * \
        clay[layer] / 100.0 * bulkDensity[layer] / ps


def volumetricFractionSilt(
    layer: int, bulkDensity: List[float], silt: List[float], ps: float, carbon: List[float], pom: float, rocks: List[float]
) -> float:
    return (1.0 - volumetricFractionOrganicMatter(layer, carbon, bulkDensity, pom) - volumetricFractionRocks(layer, rocks)) * \
        silt[layer] / 100.0 * bulkDensity[layer] / ps


def volumetricFractionSand(
    layer: int, bulkDensity: List[float], sand: List[float], ps: float, carbon: List[float], pom: float, rocks: List[float]
) -> float:
    return (1.0 - volumetricFractionOrganicMatter(layer, carbon, bulkDensity, pom) - volumetricFractionRocks(layer, rocks)) * \
        sand[layer] / 100.0 * bulkDensity[layer] / ps


def kelvinT(celciusT: float) -> float:
    celciusToKelvin = 273.18
    return celciusT + celciusToKelvin


def Divide(value1: float, value2: float, errVal: float) -> float:
    return value1 / value2 if value2 != 0.0 else errVal


def Sum(values: List[float], startIndex: int, endIndex: int, MissingValue: float) -> float:
    result = 0.0
    index = -1
    for v in values:
        index += 1
        if index >= startIndex and v != MissingValue:
            result += v
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
    if name == 'Rocks':
        return specificHeatRocks
    elif name == 'OrganicMatter':
        return specificHeatOM
    elif name == 'Sand':
        return specificHeatSand
    elif name == 'Silt':
        return specificHeatSilt
    elif name == 'Clay':
        return specificHeatClay
    elif name == 'Water':
        return specificHeatWater
    elif name == 'Ice':
        return specificHeatIce
    elif name == 'Air':
        return specificHeatAir
    return 0.0


def volumetricFractionAir(
    layer: int, rocks: List[float], carbon: List[float], bulkDensity: List[float], pom: float,
    sand: List[float], ps: float, silt: List[float], clay: List[float], soilWater: List[float]
) -> float:
    return 1.0 - volumetricFractionRocks(layer, rocks) - volumetricFractionOrganicMatter(layer, carbon, bulkDensity, pom) - \
        volumetricFractionSand(layer, bulkDensity, sand, ps, carbon, pom, rocks) - \
        volumetricFractionSilt(layer, bulkDensity, silt, ps, carbon, pom, rocks) - \
        volumetricFractionClay(layer, bulkDensity, ps, clay, carbon, pom, rocks) - \
        volumetricFractionWater(layer, soilWater, carbon, bulkDensity, pom) - \
        volumetricFractionIce(layer)


def airDensity(temperature: float, AirPressure: float) -> float:
    MWair = 0.02897
    RGAS = 8.3143
    HPA2PA = 100.0
    return Divide(MWair * AirPressure * HPA2PA, kelvinT(temperature) * RGAS, 0.0)


def longWaveRadn(emissivity: float, tDegC: float, stefanBoltzmannConstant: float) -> float:
    return stefanBoltzmannConstant * emissivity * (kelvinT(tDegC) ** 4)


def mapLayer2Node(
    layerArray: List[float],
    nodeArray: List[float],
    nodeDepth: List[float],
    numNodes: int,
    thickness: List[float],
    surfaceNode: int,
    MissingValue: float
) -> List[float]:
    for node in range(surfaceNode, numNodes):
        layer = node - 1
        if layer >= 1:
            depthLayerAbove = Sum(thickness, 1, layer, MissingValue)
        else:
            depthLayerAbove = 0.0
        d1 = depthLayerAbove - (nodeDepth[node] * 1000.0)
        d2 = nodeDepth[node + 1] * 1000.0 - depthLayerAbove
        dSum = d1 + d2
        nodeArray[node] = Divide(layerArray[layer] * d1, dSum, 0.0) + Divide(layerArray[layer + 1] * d2, dSum, 0.0)
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
    clay: List[float]
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
    if name == 'Rocks':
        result = thermalConductanceRocks
    elif name == 'OrganicMatter':
        result = thermalConductanceOM
    elif name == 'Sand':
        result = thermalConductanceSand
    elif name == 'Silt':
        result = thermalConductanceSilt
    elif name == 'Clay':
        result = thermalConductanceClay
    elif name == 'Water':
        result = thermalConductanceWater
    elif name == 'Ice':
        result = thermalConductanceIce
    elif name == 'Air':
        result = thermalConductanceAir
    elif name == 'Minerals':
        result = (thermalConductanceRocks ** volumetricFractionRocks(layer, rocks)) * \
            (thermalConductanceSand ** volumetricFractionSand(layer, bulkDensity, sand, ps, carbon, pom, rocks)) + \
            (thermalConductanceSilt ** volumetricFractionSilt(layer, bulkDensity, silt, ps, carbon, pom, rocks)) + \
            (thermalConductanceClay ** volumetricFractionClay(layer, bulkDensity, ps, clay, carbon, pom, rocks))
    # Preserve original behavior where volumetricSpecificHeat was returned at end
    return volumetricSpecificHeat(name, layer)


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
    clay: List[float]
) -> float:
    shapeFactorRocks = 0.182
    shapeFactorOM = 0.5
    shapeFactorSand = 0.182
    shapeFactorSilt = 0.125
    shapeFactorClay = 0.007755
    shapeFactorWater = 1.0
    result = 0.0
    if name == 'Rocks':
        result = shapeFactorRocks
    elif name == 'OrganicMatter':
        result = shapeFactorOM
    elif name == 'Sand':
        result = shapeFactorSand
    elif name == 'Silt':
        result = shapeFactorSilt
    elif name == 'Clay':
        result = shapeFactorClay
    elif name == 'Water':
        result = shapeFactorWater
    elif name == 'Ice':
        result = 0.333 - (0.333 * volumetricFractionIce(layer) /
                          (volumetricFractionWater(layer, soilWater, carbon, bulkDensity, pom) +
                           volumetricFractionIce(layer) +
                           volumetricFractionAir(layer, rocks, carbon, bulkDensity, pom, sand, ps, silt, clay, soilWater)))
        return result
    elif name == 'Air':
        result = 0.333 - (0.333 * volumetricFractionAir(layer, rocks, carbon, bulkDensity, pom, sand, ps, silt, clay, soilWater) /
                          (volumetricFractionWater(layer, soilWater, carbon, bulkDensity, pom) +
                           volumetricFractionIce(layer) +
                           volumetricFractionAir(layer, rocks, carbon, bulkDensity, pom, sand, ps, silt, clay, soilWater)))
        return result
    elif name == 'Minerals':
        result = shapeFactorRocks * volumetricFractionRocks(layer, rocks) + \
            (shapeFactorSand * volumetricFractionSand(layer, bulkDensity, sand, ps, carbon, pom, rocks)) + \
            (shapeFactorSilt * volumetricFractionSilt(layer, bulkDensity, silt, ps, carbon, pom, rocks)) + \
            (shapeFactorClay * volumetricFractionClay(layer, bulkDensity, ps, clay, carbon, pom, rocks))
    # Preserve original behavior where volumetricSpecificHeat was returned at end
    return volumetricSpecificHeat(name, layer)


def doUpdate(
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
    thermalConductivity: List[float]
) -> Tuple[float, List[float], List[float], List[float], float]:
    soilTemp[:] = newTemperature[:]
    if timeOfDaySecs < (internalTimeStep * 1.2):
        for node in range(surfaceNode, numNodes):
            minSoilTemp[node] = soilTemp[node]
            maxSoilTemp[node] = soilTemp[node]
    for node in range(surfaceNode, numNodes):
        if soilTemp[node] < minSoilTemp[node]:
            minSoilTemp[node] = soilTemp[node]
        elif soilTemp[node] > maxSoilTemp[node]:
            maxSoilTemp[node] = soilTemp[node]
        aveSoilTemp[node] = aveSoilTemp[node] + Divide(soilTemp[node], float(numInterationsPerDay), 0.0)
    boundaryLayerConductance = boundaryLayerConductance + Divide(thermalConductivity[airNode], float(numInterationsPerDay), 0.0)
    return boundaryLayerConductance, minSoilTemp, maxSoilTemp, aveSoilTemp, boundaryLayerConductance


def doThomas(
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
    volSpecHeatSoil: List[float]
) -> Tuple[List[float], List[float], List[float]]:
    a = [0.0] * (numNodes + 2)
    b = [0.0] * (numNodes + 1)
    c = [0.0] * (numNodes + 1)
    d = [0.0] * (numNodes + 1)

    thermalConductance[airNode] = thermalConductivity[airNode]
    for node in range(surfaceNode, numNodes):
        volumeOfSoilAtNode = 0.5 * (nodeDepth[node + 1] - nodeDepth[node - 1])
        heatStorage[node] = Divide(volSpecHeatSoil[node] * volumeOfSoilAtNode, internalTimeStep, 0.0)
        elementLength = nodeDepth[node + 1] - nodeDepth[node]
        thermalConductance[node] = Divide(thermalConductivity[node], elementLength, 0.0)

    g = 1.0 - nu
    for node in range(surfaceNode, numNodes):
        c[node] = (-nu) * thermalConductance[node]
        a[node + 1] = c[node]
        b[node] = nu * (thermalConductance[node] + thermalConductance[node - 1]) + heatStorage[node]
        d[node] = g * thermalConductance[node - 1] * soilTemp[node - 1] + \
            ((heatStorage[node] - (g * (thermalConductance[node] + thermalConductance[node - 1]))) * soilTemp[node]) + \
            (g * thermalConductance[node] * soilTemp[node + 1])

    a[surfaceNode] = 0.0
    sensibleHeatFlux = nu * thermalConductance[airNode] * newTemps[airNode]

    radnNet = 0.0
    if netRadiationSource == 'calc':
        radnNet = Divide(netRadiation * 1000000.0, internalTimeStep, 0.0)
    elif netRadiationSource == 'eos':
        radnNet = Divide(waterBalance_Eos * latentHeatOfVapourisation, timestep, 0.0)

    latentHeatFlux = Divide(waterBalance_Es * latentHeatOfVapourisation, timestep, 0.0)
    soilSurfaceHeatFlux = sensibleHeatFlux + radnNet - latentHeatFlux
    d[surfaceNode] = d[surfaceNode] + soilSurfaceHeatFlux
    d[numNodes] = d[numNodes] + (nu * thermalConductance[numNodes] * newTemps[numNodes + 1])

    for node in range(surfaceNode, numNodes - 1):
        c[node] = Divide(c[node], b[node], 0.0)
        d[node] = Divide(d[node], b[node], 0.0)
        b[node + 1] = b[node + 1] - (a[node + 1] * c[node])
        d[node + 1] = d[node + 1] - (a[node + 1] * d[node])

    newTemps[numNodes] = Divide(d[numNodes], b[numNodes], 0.0)

    for node in range(numNodes - 1, surfaceNode - 1, -1):
        newTemps[node] = d[node] - (c[node] * newTemps[node + 1])

    return newTemps, heatStorage, thermalConductance


def getBoundaryLayerConductance(
    TNew_zb: List[float],
    weather_AirPressure: float,
    stefanBoltzmannConstant: float,
    waterBalance_Eos: float,
    weather_Wind: float,
    airTemperature: float,
    surfaceNode: int,
    waterBalance_Eo: float,
    instrumentHeight: float,
    canopyHeight: float
) -> float:
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
        denom_mom = math.log(Divide(instrumentHeight - d + roughnessFactorMomentum, roughnessFactorMomentum, 0.0))
        frictionVelocity = Divide(weather_Wind * vonKarmanConstant, denom_mom + stabilityCorrectionMomentum, 0.0)
        denom_heat = math.log(Divide(instrumentHeight - d + roughnessFactorHeat, roughnessFactorHeat, 0.0))
        boundaryLayerCond = Divide(SpecificHeatAir * vonKarmanConstant * frictionVelocity, denom_heat + stabilityCorrectionHeat, 0.0)
        boundaryLayerCond = boundaryLayerCond + radiativeConductance
        heatFluxDensity = boundaryLayerCond * (surfaceTemperature - airTemperature)

        stabilityParammeter = Divide((-vonKarmanConstant) * instrumentHeight * gravitationalConstant * heatFluxDensity,
                                     SpecificHeatAir * kelvinT(airTemperature) * (frictionVelocity ** 3.0), 0.0)
        if stabilityParammeter > 0.0:
            stabilityCorrectionHeat = 4.7 * stabilityParammeter
            stabilityCorrectionMomentum = stabilityCorrectionHeat
        else:
            stabilityCorrectionHeat = (-2.0) * math.log((1.0 + math.sqrt(1.0 - (16.0 * stabilityParammeter))) / 2.0)
            stabilityCorrectionMomentum = 0.6 * stabilityCorrectionHeat
    return boundaryLayerCond


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
    stefanBoltzmannConstant: float
) -> float:
    surfaceEmissivity = 0.96
    w2MJ = internalTimeStep / 1000000.0
    emissivityAtmos = (1.0 - (0.84 * cloudFr)) * 0.58 * (cva ** (1.0 / 7.0)) + (0.84 * cloudFr)
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
    defaultTimeOfMaximumTemperature: float
) -> float:
    time = timeHours / 24.0
    maxT_time = defaultTimeOfMaximumTemperature / 24.0
    minT_time = maxT_time - 0.5
    if time < minT_time:
        midnightT = math.sin((0.0 + 0.25 - maxT_time) * 2.0 * math.pi) * (maxTempYesterday - minTempYesterday) / 2.0 + \
            ((maxTempYesterday + minTempYesterday) / 2.0)
        tScale = (minT_time - time) / minT_time
        if tScale > 1.0:
            tScale = 1.0
        elif tScale < 0.0:
            tScale = 0.0
        currentTemperature = weather_MinT + (tScale * (midnightT - weather_MinT))
        return currentTemperature
    else:
        currentTemperature = math.sin((time + 0.25 - maxT_time) * 2.0 * math.pi) * (weather_MaxT - weather_MinT) / 2.0 + weather_MeanT
        return currentTemperature


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
    MissingValue: float
) -> List[float]:
    thermCondLayers = [0.0] * (numNodes + 1)
    for node in range(1, numNodes):
        numerator = 0.0
        denominator = 0.0
        for constituentName in soilConstituentNames:
            shapeFactorConstituent = shapeFactor(constituentName, node, soilWater, carbon, bulkDensity, pom, rocks, sand, ps, silt, clay)
            thermalConductanceConstituent = ThermalConductance(constituentName, node, rocks, bulkDensity, sand, ps, carbon, pom, silt, clay)
            thermalConductanceWater = ThermalConductance('Water', node, rocks, bulkDensity, sand, ps, carbon, pom, silt, clay)
            if thermalConductanceWater == 0.0:
                k = 0.0
            else:
                k = 2.0 / 3.0 * ((1.0 + (shapeFactorConstituent * (thermalConductanceConstituent / thermalConductanceWater - 1.0))) ** (-1)) + \
                    (1.0 / 3.0 * ((1.0 + (shapeFactorConstituent * (thermalConductanceConstituent / thermalConductanceWater - 1.0) *
                                           (1.0 - (2.0 * shapeFactorConstituent)))) ** (-1)))
            numerator += thermalConductanceConstituent * soilWater[node] * k
            denominator += soilWater[node] * k
        thermCondLayers[node] = numerator / denominator if denominator != 0.0 else 0.0
    thermalConductivity_out = [0.0] * (numNodes + 1) if thermalConductivity is None else thermalConductivity[:]
    thermalConductivity_out = mapLayer2Node(thermCondLayers, thermalConductivity_out, nodeDepth, numNodes, thickness, surfaceNode, MissingValue)
    return thermalConductivity_out


def doVolumetricSpecificHeat(
    soilConstituentNames: List[str],
    numNodes: int,
    volSpecHeatSoil: Optional[List[float]],
    soilWater: List[float],
    nodeDepth: List[float],
    thickness: List[float],
    surfaceNode: int,
    MissingValue: float
) -> List[float]:
    volspecHeatSoil_ = [0.0] * (numNodes + 1)
    for node in range(1, numNodes):
        val = 0.0
        for constituentName in soilConstituentNames:
            if constituentName != 'Minerals':
                val += volumetricSpecificHeat(constituentName, node) * 1000000.0 * soilWater[node]
        volspecHeatSoil_[node] = val
    volSpecHeatSoil_out = [0.0] * (numNodes + 1) if volSpecHeatSoil is None else volSpecHeatSoil[:]
    volSpecHeatSoil_out = mapLayer2Node(volspecHeatSoil_, volSpecHeatSoil_out, nodeDepth, numNodes, thickness, surfaceNode, MissingValue)
    return volSpecHeatSoil_out


def Zero(arr: Optional[List[float]]) -> List[float]:
    if arr is None:
        return []
    return [0.0 for _ in arr]


def doNetRadiation(
    solarRadn: List[float],
    cloudFr: float,
    cva: float,
    ITERATIONSperDAY: int,
    weather_MinT: float,
    clock_Today_DayOfYear: int,
    weather_Radn: float,
    weather_Latitude: float
) -> Tuple[List[float], float, float]:
    m1 = [0.0] * (ITERATIONSperDAY + 1)
    TSTEPS2RAD = Divide(2.0 * math.pi, float(ITERATIONSperDAY), 0.0)
    solarConstant = 1360.0
    solarDeclination = 0.3985 * math.sin((4.869 + (clock_Today_DayOfYear * 2.0 * math.pi / 365.25) + (0.03345 *
                                    math.sin((6.224 + (clock_Today_DayOfYear * 2.0 * math.pi / 365.25))))))
    cD = math.sqrt(max(0.0, 1.0 - (solarDeclination * solarDeclination)))
    m1Tot = 0.0
    for timestepNumber in range(1, ITERATIONSperDAY + 1):
        m1[timestepNumber] = (solarDeclination * math.sin(weather_Latitude * math.pi / 180.0) +
                              (cD * math.cos(weather_Latitude * math.pi / 180.0) *
                               math.cos(TSTEPS2RAD * (timestepNumber - (ITERATIONSperDAY / 2.0))))) * 24.0 / float(ITERATIONSperDAY)
        if m1[timestepNumber] > 0.0:
            m1Tot += m1[timestepNumber]
        else:
            m1[timestepNumber] = 0.0
    psr = m1Tot * solarConstant * 3600.0 / 1000000.0
    fr = Divide(max(weather_Radn, 0.1), psr, 0.0)
    cloudFr = 2.33 - (3.33 * fr)
    cloudFr = min(max(cloudFr, 0.0), 1.0)
    for timestepNumber in range(1, ITERATIONSperDAY + 1):
        solarRadn[timestepNumber] = max(weather_Radn, 0.1) * Divide(m1[timestepNumber], m1Tot, 0.0)
    cva = math.exp((31.3716 - (6014.79 / kelvinT(weather_MinT)) - (0.00792495 * kelvinT(weather_MinT)))) / kelvinT(weather_MinT)
    return solarRadn, cloudFr, cva


def doProcess(
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
    nu: float
) -> Tuple[
    float, float, List[float], List[float], float, float, List[float], float, List[float], float, List[float],
    float, List[float], List[float], List[float], List[float], List[float]
]:
    interactionsPerDay = 48
    cva = 0.0
    cloudFr = 0.0
    solarRadn = [0.0] * (interactionsPerDay + 1)
    solarRadn, cloudFr, cva = doNetRadiation(solarRadn, cloudFr, cva, interactionsPerDay, weather_MinT, clock_Today_DayOfYear, weather_Radn, weather_Latitude)

    minSoilTemp = Zero(minSoilTemp)
    maxSoilTemp = Zero(maxSoilTemp)
    aveSoilTemp = Zero(aveSoilTemp)
    boundaryLayerConductance = 0.0
    internalTimeStep = int(round(timestep / float(interactionsPerDay)))

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

        if boundarLayerConductanceSource == 'constant':
            thermalConductivity[airNode] = constantBoundaryLayerConductance
        elif boundarLayerConductanceSource == 'calc':
            thermalConductivity[airNode] = getBoundaryLayerConductance(newTemperature, weather_AirPressure, stefanBoltzmannConstant, waterBalance_Eos, weather_Wind, airTemperature, surfaceNode, waterBalance_Eo, instrumentHeight, canopyHeight)
            for _ in range(1, numIterationsForBoundaryLayerConductance + 1):
                newTemperature, heatStorage, thermalConductance = doThomas(newTemperature, netRadiation, heatStorage, waterBalance_Eos, numNodes, timestep, netRadiationSource, latentHeatOfVapourisation, nodeDepth, waterBalance_Es, airNode, soilTemp, surfaceNode, internalTimeStep, thermalConductance, thermalConductivity, nu, volSpecHeatSoil)
                thermalConductivity[airNode] = getBoundaryLayerConductance(newTemperature, weather_AirPressure, stefanBoltzmannConstant, waterBalance_Eos, weather_Wind, airTemperature, surfaceNode, waterBalance_Eo, instrumentHeight, canopyHeight)

        newTemperature, heatStorage, thermalConductance = doThomas(newTemperature, netRadiation, heatStorage, waterBalance_Eos, numNodes, timestep, netRadiationSource, latentHeatOfVapourisation, nodeDepth, waterBalance_Es, airNode, soilTemp, surfaceNode, internalTimeStep, thermalConductance, thermalConductivity, nu, volSpecHeatSoil)

        boundaryLayerConductance, minSoilTemp, maxSoilTemp, aveSoilTemp, _ = doUpdate(
            interactionsPerDay, timeOfDaySecs, boundaryLayerConductance, minSoilTemp, airNode, soilTemp, newTemperature, numNodes, surfaceNode, internalTimeStep, maxSoilTemp, aveSoilTemp, thermalConductivity
        )

        if abs(timeOfDaySecs - (5.0 * 3600.0)) <= (min(timeOfDaySecs, 5.0 * 3600.0) * 0.0001):
            morningSoilTemp = soilTemp[:]

    minTempYesterday = weather_MinT
    maxTempYesterday = weather_MaxT

    return (
        timeOfDaySecs,
        netRadiation,
        minSoilTemp,
        maxSoilTemp,
        boundaryLayerConductance,
        maxTempYesterday,
        soilTemp,
        airTemperature,
        newTemperature,
        float(internalTimeStep),
        thermalConductivity,
        minTempYesterday,
        aveSoilTemp,
        morningSoilTemp,
        volSpecHeatSoil,
        heatStorage,
        thermalConductance,
    )


def ToCumThickness(Thickness: List[float]) -> List[float]:
    CumThickness = [0.0] * len(Thickness)
    if len(Thickness) > 0:
        CumThickness[0] = Thickness[0]
        for Layer in range(1, len(Thickness)):
            CumThickness[Layer] = Thickness[Layer] + CumThickness[Layer - 1]
    return CumThickness


def calcSoilTemperature(
    soilTempIO: List[float],
    weather_Tav: float,
    clock_Today_DayOfYear: int,
    surfaceNode: int,
    numNodes: int,
    weather_Amp: float,
    thickness: List[float],
    weather_Latitude: float
) -> List[float]:
    soilTemp = [0.0] * (numNodes + 2)
    cumulativeDepth = ToCumThickness(thickness)
    w = 2.0 * math.pi / (365.25 * 24.0 * 3600.0)
    dh = 0.6
    zd = math.sqrt(2.0 * dh / w)
    offset = 0.25
    if weather_Latitude > 0.0:
        offset = -0.25
    for nodes in range(1, numNodes):
        soilTemp[nodes] = weather_Tav + (weather_Amp * math.exp(-1.0 * cumulativeDepth[nodes] / zd) *
                                         math.sin(((clock_Today_DayOfYear / 365.0 + offset) * 2.0 * math.pi -
                                                   (cumulativeDepth[nodes] / zd))))
    for i in range(surfaceNode, surfaceNode + numNodes):
        if (i - surfaceNode) < len(soilTemp):
            soilTempIO[i] = soilTemp[i - surfaceNode]
    return soilTempIO


def calcSurfaceTemperature(weather_MeanT: float, weather_MaxT: float, waterBalance_Salb: float, weather_Radn: float) -> float:
    return (1.0 - waterBalance_Salb) * (weather_MeanT + ((weather_MaxT - weather_MeanT) * math.sqrt(max(weather_Radn, 0.1) * 23.8846 / 800.0))) + \
        (waterBalance_Salb * weather_MeanT)


def ValuesInArray(Values: Optional[List[float]], MissingValue: float) -> bool:
    if Values is not None:
        for Value in Values:
            if Value != MissingValue and not (Value != Value):  # not NaN
                return True
    return False
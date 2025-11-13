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
    # Outputs (returned)
) -> Tuple[
    List[float],  # InitialValues
    bool,         # doInitialisationStuff
    float,        # internalTimeStep
    float,        # timeOfDaySecs
    int,          # numNodes
    int,          # numLayers
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
    List[float],  # aveSoilWater
    List[float],  # thickness
    List[float],  # bulkDensity
    List[float],  # rocks
    List[float],  # carbon
    List[float],  # sand
    List[float],  # silt
    List[float],  # clay
    float,        # instrumentHeight
    float,        # netRadiation
    float,        # canopyHeight
    float,        # instrumHeight
    List[float],  # nodeDepth
    List[float],  # thermCondPar1
    List[float],  # thermCondPar2
    List[float],  # thermCondPar3
    List[float],  # thermCondPar4
]:
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

    instrumentHeight = getIniVariables(instrumentHeight, instrumHeight, defaultInstrumentHeight)

    (
        heatStorage,
        minSoilTemp,
        bulkDensity,
        numNodes,
        maxSoilTemp,
        waterBalance_SW_alloc,
        nodeDepth,
        topsoilNode_local,
        newTemperature,
        surfaceNode_local,
        soilWater,
        thermalConductance,
        thermalConductivity,
        sand,
        carbon,
        thickness,
        numPhantomNodes_local,
        rocks,
        clay,
        airNode_local,
        soilTemp,
        numLayers,
        silt,
        volSpecHeatSoil,
        aveSoilTemp,
        morningSoilTemp,
    ) = getProfileVariables(
        heatStorage=None,
        minSoilTemp=None,
        bulkDensity=None,
        numNodes=numNodes,
        physical_BD=physical_BD,
        maxSoilTemp=None,
        waterBalance_SW=waterBalance_SW,
        organic_Carbon=organic_Carbon,
        physical_Rocks=physical_Rocks,
        nodeDepth=nodeDepth,
        topsoilNode=topsoilNode,
        newTemperature=None,
        surfaceNode=surfaceNode,
        soilWater=None,
        thermalConductance=None,
        thermalConductivity=None,
        sand=None,
        carbon=None,
        thickness=None,
        numPhantomNodes=numPhantomNodes,
        physical_ParticleSizeSand=physical_ParticleSizeSand,
        rocks=None,
        clay=None,
        physical_ParticleSizeSilt=physical_ParticleSizeSilt,
        airNode=airNode,
        physical_ParticleSizeClay=physical_ParticleSizeClay,
        soilTemp=None,
        numLayers=0,
        physical_Thickness=physical_Thickness,
        silt=None,
        volSpecHeatSoil=None,
        aveSoilTemp=None,
        morningSoilTemp=None,
        DepthToConstantTemperature=DepthToConstantTemperature,
        MissingValue=MissingValue,
    )

    thermCondPar1, thermCondPar2, thermCondPar3, thermCondPar4, soilTemp, newTemperature, soilRoughnessHeight = readParam(
        bareSoilRoughness=bareSoilRoughness,
        newTemperature=newTemperature,
        soilRoughnessHeight=soilRoughnessHeight,
        soilTemp=soilTemp,
        thermCondPar2=thermCondPar2,
        numLayers=numLayers,
        bulkDensity=bulkDensity,
        numNodes=numNodes,
        thermCondPar3=thermCondPar3,
        thermCondPar4=thermCondPar4,
        clay=clay,
        thermCondPar1=thermCondPar1,
        weather_Tav=weather_Tav,
        clock_Today_DayOfYear=clock_Today_DayOfYear,
        surfaceNode=surfaceNode,
        weather_Amp=weather_Amp,
        thickness=thickness,
        weather_Latitude=weather_Latitude,
    )

    InitialValues: List[float] = []
    if pInitialValues is not None:
        InitialValues = list(pInitialValues)

    aveSoilWater = []
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
        newTemperature,
        airTemperature,
        maxTempYesterday,
        minTempYesterday,
        soilWater,
        minSoilTemp,
        maxSoilTemp,
        aveSoilTemp,
        aveSoilWater,
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
    thermCondPar1: Optional[List[float]],
    thermCondPar2: Optional[List[float]],
    thermCondPar3: Optional[List[float]],
    thermCondPar4: Optional[List[float]],
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
    List[float],  # heatStorage
    float,        # instrumentHeight
    float,        # canopyHeight
    List[float],  # minSoilTemp
    List[float],  # maxSoilTemp
    List[float],  # aveSoilTemp
    List[float],  # volSpecHeatSoil
    List[float],  # soilTemp
    List[float],  # morningSoilTemp
    List[float],  # newTemperature
    List[float],  # thermalConductivity
    List[float],  # thermalConductance
    float,        # boundaryLayerConductance
    List[float],  # soilWater
    bool,         # doInitialisationStuff
    float,        # maxTempYesterday
    float,        # minTempYesterday
    float,        # airTemperature
    float,        # internalTimeStep
    float,        # timeOfDaySecs
    float,        # netRadiation
    Optional[List[float]],  # InitialValues
]:
    soilWater, instrumentHeight, canopyHeight = getOtherVariables(
        numLayers, numNodes, soilWater, instrumentHeight, soilRoughnessHeight, waterBalance_SW, microClimate_CanopyHeight, canopyHeight
    )

    if doInitialisationStuff:
        if ValuesInArray(InitialValues, MissingValue):
            soilTemp = [0.0] * (numNodes + 2)
            if InitialValues is None:
                InitialValues = []
            # Copy initial values into soilTemp starting at topsoilNode
            ninit = len(InitialValues)
            for idx in range(ninit):
                soilTemp[topsoilNode + idx] = InitialValues[idx]
        else:
            soilTemp = calcSoilTemperature(
                soilTempIO=[0.0] * (numNodes + 2),
                weather_Tav=weather_Tav,
                clock_Today_DayOfYear=clock_Today_DayOfYear,
                surfaceNode=surfaceNode,
                numNodes=numNodes,
                weather_Amp=weather_Amp,
                thickness=thickness,
                weather_Latitude=weather_Latitude,
            )
            InitialValues = [0.0] * numLayers
            for i in range(numLayers):
                InitialValues[i] = soilTemp[topsoilNode + i]

        if len(soilTemp) < (numNodes + 2):
            soilTemp = (soilTemp or []) + [weather_Tav] * ((numNodes + 2) - len(soilTemp))

        soilTemp[airNode] = weather_MeanT
        soilTemp[surfaceNode] = calcSurfaceTemperature(weather_MeanT, weather_MaxT, waterBalance_Salb, weather_Radn)

        for i in range(numNodes + 1, len(soilTemp) - 1):
            soilTemp[i] = weather_Tav

        newTemperature = soilTemp.copy()
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
        airNode_local,
        soilTemp,
        airTemperature,
        newTemperature,
        thermalConductivity,
        minTempYesterday,
        aveSoilTemp,
        morningSoilTemp,
        weather_MeanT_local,
        constantBoundaryLayerConductance_local,
        weather_MinT_local,
        clock_Today_DayOfYear_local,
        weather_Radn_local,
        weather_Latitude_local,
        soilConstituentNames_local,
        numNodes_local,
        volSpecHeatSoil,
        soilWater,
        nodeDepth_local,
        thickness_local,
        surfaceNode_local,
        MissingValue_local,
        carbon_local,
        bulkDensity_local,
        pom_local,
        rocks_local,
        sand_local,
        ps_local,
        silt_local,
        clay_local,
        defaultTimeOfMaximumTemperature_local,
        waterBalance_Eo_local,
        waterBalance_Eos_local,
        waterBalance_Salb_local,
        stefanBoltzmannConstant_local,
        weather_AirPressure_local,
        weather_Wind_local,
        instrumentHeight,
        canopyHeight,
        heatStorage,
        netRadiationSource_local,
        latentHeatOfVapourisation_local,
        waterBalance_Es_local,
        thermalConductance,
        nu_local,
        internalTimeStep,
        boundaryLayerConductance_daily,
    ) = doProcess(
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
        ps,
        silt,
        clay,
        defaultTimeOfMaximumTemperature,
        waterBalance_Eo,
        waterBalance_Eos,
        waterBalance_Salb,
        stefanBoltzmannConstant,
        weather_AirPressure,
        weather_Wind,
        instrumentHeight,
        canopyHeight,
        heatStorage,
        netRadiationSource,
        latentHeatOfVapourisation,
        waterBalance_Es,
        thermalConductance,
        nu,
    )

    boundaryLayerConductance = boundaryLayerConductance_daily

    return (
        heatStorage,
        instrumentHeight,
        canopyHeight,
        minSoilTemp,
        maxSoilTemp,
        aveSoilTemp,
        volSpecHeatSoil,
        soilTemp,
        morningSoilTemp,
        newTemperature,
        thermalConductivity,
        thermalConductance,
        boundaryLayerConductance,
        soilWater,
        doInitialisationStuff,
        maxTempYesterday,
        minTempYesterday,
        airTemperature,
        internalTimeStep,
        timeOfDaySecs,
        netRadiation,
        InitialValues,
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
    List[float],  # heatStorage
    List[float],  # minSoilTemp
    List[float],  # bulkDensity
    int,          # numNodes
    List[float],  # maxSoilTemp
    Optional[List[float]],  # waterBalance_SW
    List[float],  # nodeDepth
    int,          # topsoilNode (unchanged passthrough)
    List[float],  # newTemperature
    int,          # surfaceNode (unchanged passthrough)
    List[float],  # soilWater
    List[float],  # thermalConductance
    List[float],  # thermalConductivity
    List[float],  # sand
    List[float],  # carbon
    List[float],  # thickness
    int,          # numPhantomNodes (unchanged passthrough)
    List[float],  # rocks
    List[float],  # clay
    int,          # airNode (unchanged passthrough)
    List[float],  # soilTemp
    int,          # numLayers
    List[float],  # silt
    List[float],  # volSpecHeatSoil
    List[float],  # aveSoilTemp
    List[float],  # morningSoilTemp
]:
    numLayers = len(physical_Thickness)
    numNodes = numLayers + numPhantomNodes

    thickness = [0.0] * (numLayers + numPhantomNodes + 1)
    # Copy physical thickness into indices 1..numLayers
    for i in range(numLayers):
        thickness[1 + i] = physical_Thickness[i]

    belowProfileDepth = max(DepthToConstantTemperature - Sum(thickness, 1, numLayers, MissingValue), 1000.0)
    thicknessForPhantomNodes = belowProfileDepth * 2.0 / float(numPhantomNodes)
    firstPhantomNode = numLayers
    for i in range(firstPhantomNode, firstPhantomNode + numPhantomNodes):
        thickness[i] = thicknessForPhantomNodes

    oldDepth = nodeDepth.copy() if nodeDepth is not None else None
    nodeDepth = [0.0] * (numNodes + 2)
    if oldDepth is not None:
        ncopy = min(len(nodeDepth), len(oldDepth))
        for i in range(ncopy):
            nodeDepth[i] = oldDepth[i]

    nodeDepth[airNode] = 0.0
    nodeDepth[surfaceNode] = 0.0
    nodeDepth[topsoilNode] = 0.5 * thickness[1] / 1000.0
    for node in range(topsoilNode, numNodes + 1):
        nodeDepth[node + 1] = (Sum(thickness, 1, node - 1, MissingValue) + (0.5 * thickness[node])) / 1000.0

    oldBulk = bulkDensity.copy() if bulkDensity is not None else None
    bulkDensity = [0.0] * (numLayers + 1 + numPhantomNodes)
    if oldBulk is not None:
        ncopy = min(len(bulkDensity), len(oldBulk))
        for i in range(ncopy):
            bulkDensity[i] = oldBulk[i]
    for i in range(len(physical_BD)):
        bulkDensity[1 + i] = physical_BD[i]
    bulkDensity[numNodes] = bulkDensity[numLayers]
    for layer in range(numLayers, numLayers + numPhantomNodes):
        bulkDensity[layer] = bulkDensity[numLayers]

    oldSW = soilWater.copy() if soilWater is not None else None
    soilWater = [0.0] * (numLayers + 1 + numPhantomNodes)
    if oldSW is not None:
        ncopy = min(len(soilWater), len(oldSW))
        for i in range(ncopy):
            soilWater[i] = oldSW[i]
    if waterBalance_SW is not None:
        for layer in range(1, numLayers + 1):
            # layer index mapping: use waterBalance_SW[layer-1], thickness[layer-1], thickness[layer]
            soilWater[layer] = Divide(waterBalance_SW[layer - 1] * thickness[layer - 1], thickness[layer], 0.0)
        for layer in range(numLayers, numLayers + numPhantomNodes):
            soilWater[layer] = soilWater[numLayers]

    carbon = [0.0] * (numLayers + 1 + numPhantomNodes)
    for layer in range(1, numLayers + 1):
        carbon[layer] = organic_Carbon[layer - 1]
    for layer in range(numLayers, numLayers + numPhantomNodes):
        carbon[layer] = carbon[numLayers]

    rocks = [0.0] * (numLayers + 1 + numPhantomNodes)
    for layer in range(1, numLayers + 1):
        rocks[layer] = physical_Rocks[layer - 1]
    for layer in range(numLayers, numLayers + numPhantomNodes):
        rocks[layer] = rocks[numLayers]

    sand = [0.0] * (numLayers + 1 + numPhantomNodes)
    for layer in range(1, numLayers + 1):
        sand[layer] = physical_ParticleSizeSand[layer - 1]
    for layer in range(numLayers, numLayers + numPhantomNodes):
        sand[layer] = sand[numLayers]

    silt = [0.0] * (numLayers + 1 + numPhantomNodes)
    for layer in range(1, numLayers + 1):
        silt[layer] = physical_ParticleSizeSilt[layer - 1]
    for layer in range(numLayers, numLayers + numPhantomNodes):
        silt[layer] = silt[numLayers]

    clay = [0.0] * (numLayers + 1 + numPhantomNodes)
    for layer in range(1, numLayers + 1):
        clay[layer] = physical_ParticleSizeClay[layer - 1]
    for layer in range(numLayers, numLayers + numPhantomNodes):
        clay[layer] = clay[numLayers]

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

    return (
        heatStorage,
        minSoilTemp,
        bulkDensity,
        numNodes,
        maxSoilTemp,
        waterBalance_SW,
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
        rocks,
        clay,
        airNode,
        soilTemp,
        numLayers,
        silt,
        volSpecHeatSoil,
        aveSoilTemp,
        morningSoilTemp,
    )


def doThermalConductivityCoeffs(
    numLayers: int,
    bulkDensity: List[float],
    numNodes: int,
    clay: List[float],
    thermCondPar1: Optional[List[float]],
    thermCondPar2: Optional[List[float]],
    thermCondPar3: Optional[List[float]],
    thermCondPar4: Optional[List[float]],
) -> Tuple[List[float], List[float], List[float], List[float]]:
    thermCondPar1 = [0.0] * (numNodes + 1)
    thermCondPar2 = [0.0] * (numNodes + 1)
    thermCondPar3 = [0.0] * (numNodes + 1)
    thermCondPar4 = [0.0] * (numNodes + 1)

    for layer in range(1, numLayers + 2):
        element = layer
        thermCondPar1[element] = 0.65 - (0.78 * bulkDensity[layer]) + (0.6 * (bulkDensity[layer] ** 2))
        thermCondPar2[element] = 1.06 * bulkDensity[layer]
        thermCondPar3[element] = 1.0 + Divide(2.6, math.sqrt(max(clay[layer], 1e-12)), 0.0)
        thermCondPar4[element] = 0.03 + (0.1 * (bulkDensity[layer] ** 2))
    return thermCondPar1, thermCondPar2, thermCondPar3, thermCondPar4


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
) -> Tuple[List[float], List[float], List[float], List[float], List[float], List[float], float]:
    thermCondPar1, thermCondPar2, thermCondPar3, thermCondPar4 = doThermalConductivityCoeffs(
        numLayers, bulkDensity, numNodes, clay, thermCondPar1, thermCondPar2, thermCondPar3, thermCondPar4
    )
    soilTemp = calcSoilTemperature(
        soilTempIO=soilTemp,
        weather_Tav=weather_Tav,
        clock_Today_DayOfYear=clock_Today_DayOfYear,
        surfaceNode=surfaceNode,
        numNodes=numNodes,
        weather_Amp=weather_Amp,
        thickness=thickness,
        weather_Latitude=weather_Latitude,
    )
    newTemperature = soilTemp.copy()
    soilRoughnessHeight = bareSoilRoughness
    return thermCondPar1, thermCondPar2, thermCondPar3, thermCondPar4, soilTemp, newTemperature, soilRoughnessHeight


def getOtherVariables(
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
        soilWater[1 + i] = waterBalance_SW[i]
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
    layer: int,
    rocks: List[float],
    carbon: List[float],
    bulkDensity: List[float],
    pom: float,
    sand: List[float],
    ps: float,
    silt: List[float],
    clay: List[float],
    soilWater: List[float],
) -> float:
    return 1.0 - volumetricFractionRocks(layer, rocks) - volumetricFractionOrganicMatter(layer, carbon, bulkDensity, pom) - \
        volumetricFractionSand(layer, bulkDensity, sand, ps, carbon, pom, rocks) - \
        volumetricFractionSilt(layer, bulkDensity, silt, ps, carbon, pom, rocks) - \
        volumetricFractionClay(layer, bulkDensity, ps, clay, carbon, pom, rocks) - \
        volumetricFractionWater(layer, soilWater, carbon, bulkDensity, pom) - volumetricFractionIce(layer)


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
    MissingValue: float,
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
    # Preserve Fortran's final assignment (likely bug in original): set to volumetricSpecificHeat
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
        denom = volumetricFractionWater(layer, soilWater, carbon, bulkDensity, pom) + volumetricFractionIce(layer) + \
            volumetricFractionAir(layer, rocks, carbon, bulkDensity, pom, sand, ps, silt, clay, soilWater)
        result = 0.333 - (0.333 * volumetricFractionIce(layer) / denom)
        # Preserve Fortran's behavior of immediate return
        return result
    elif name == 'Air':
        denom = volumetricFractionWater(layer, soilWater, carbon, bulkDensity, pom) + volumetricFractionIce(layer) + \
            volumetricFractionAir(layer, rocks, carbon, bulkDensity, pom, sand, ps, silt, clay, soilWater)
        result = 0.333 - (0.333 * volumetricFractionAir(layer, rocks, carbon, bulkDensity, pom, sand, ps, silt, clay, soilWater) / denom)
        return result
    elif name == 'Minerals':
        result = shapeFactorRocks * volumetricFractionRocks(layer, rocks) + \
                 (shapeFactorSand * volumetricFractionSand(layer, bulkDensity, sand, ps, carbon, pom, rocks)) + \
                 (shapeFactorSilt * volumetricFractionSilt(layer, bulkDensity, silt, ps, carbon, pom, rocks)) + \
                 (shapeFactorClay * volumetricFractionClay(layer, bulkDensity, ps, clay, carbon, pom, rocks))
    # Preserve Fortran's final assignment (likely bug in original)
    result = volumetricSpecificHeat(name, layer)
    return result


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
    thermalConductivity: List[float],
) -> Tuple[float, List[float], List[float], List[float], float]:
    soilTemp = newTemperature.copy()
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
    volSpecHeatSoil: List[float],
) -> Tuple[List[float], List[float], List[float], List[float]]:
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

    for node in range(surfaceNode, numNodes):
        c[node] = Divide(c[node], b[node], 0.0)
        d[node] = Divide(d[node], b[node], 0.0)
        b[node + 1] = b[node + 1] - (a[node + 1] * c[node])
        d[node + 1] = d[node + 1] - (a[node + 1] * d[node])

    newTemps[numNodes] = Divide(d[numNodes], b[numNodes], 0.0)

    for node in range(numNodes - 1, surfaceNode - 1, -1):
        newTemps[node] = d[node] - (c[node] * newTemps[node + 1])

    return newTemps, heatStorage, thermalConductance, thermalConductivity


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
    canopyHeight: float,
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

    for _ in range(3):
        denom_m = math.log(Divide(instrumentHeight - d + roughnessFactorMomentum, roughnessFactorMomentum, 0.0)) + stabilityCorrectionMomentum
        frictionVelocity = Divide(weather_Wind * vonKarmanConstant, denom_m, 0.0)
        denom_h = math.log(Divide(instrumentHeight - d + roughnessFactorHeat, roughnessFactorHeat, 0.0)) + stabilityCorrectionHeat
        boundaryLayerCond = Divide(SpecificHeatAir * vonKarmanConstant * frictionVelocity, denom_h, 0.0)
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
    stefanBoltzmannConstant: float,
) -> float:
    surfaceEmissivity = 0.96
    w2MJ = internalTimeStep / 1000000.0
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
    thermalConductivity: List[float],
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
    for node in range(1, numNodes):
        numerator = 0.0
        denominator = 0.0
        for constituentName in soilConstituentNames:
            shapeFactorConstituent = shapeFactor(constituentName, node, soilWater, carbon, bulkDensity, pom, rocks, sand, ps, silt, clay)
            thermalConductanceConstituent = ThermalConductance(constituentName, node, rocks, bulkDensity, sand, ps, carbon, pom, silt, clay)
            thermalConductanceWater = ThermalConductance('Water', node, rocks, bulkDensity, sand, ps, carbon, pom, silt, clay)
            k = 2.0 / 3.0 * ((1 + (shapeFactorConstituent * (thermalConductanceConstituent / thermalConductanceWater - 1.0))) ** (-1)) + \
                (1.0 / 3.0 * ((1 + (shapeFactorConstituent * (thermalConductanceConstituent / thermalConductanceWater - 1.0) *
                                      (1 - (2 * shapeFactorConstituent)))) ** (-1)))
            numerator += thermalConductanceConstituent * soilWater[node] * k
            denominator += soilWater[node] * k
        thermCondLayers[node] = Divide(numerator, denominator, 0.0)
    thermalConductivity = mapLayer2Node(thermCondLayers, thermalConductivity, nodeDepth, numNodes, thickness, surfaceNode, MissingValue)
    return thermalConductivity


def doVolumetricSpecificHeat(
    soilConstituentNames: List[str],
    numNodes: int,
    volSpecHeatSoil: List[float],
    soilWater: List[float],
    nodeDepth: List[float],
    thickness: List[float],
    surfaceNode: int,
    MissingValue: float,
) -> List[float]:
    volspecHeatSoil_ = [0.0] * (numNodes + 1)
    for node in range(1, numNodes):
        volspecHeatSoil_[node] = 0.0
        for constituentName in soilConstituentNames:
            if constituentName != 'Minerals':
                volspecHeatSoil_[node] += volumetricSpecificHeat(constituentName, node) * 1000000.0 * soilWater[node]
    volSpecHeatSoil = mapLayer2Node(volspecHeatSoil_, volSpecHeatSoil, nodeDepth, numNodes, thickness, surfaceNode, MissingValue)
    return volSpecHeatSoil


def Zero(arr: Optional[List[float]]) -> List[float]:
    if arr is None:
        return []
    for i in range(len(arr)):
        arr[i] = 0.0
    return arr


def doNetRadiation(
    ITERATIONSperDAY: int,
    weather_MinT: float,
    clock_Today_DayOfYear: int,
    weather_Radn: float,
    weather_Latitude: float,
) -> Tuple[List[float], float, float]:
    solarRadn = [0.0] * (ITERATIONSperDAY + 1)
    cloudFr = 0.0
    cva = 0.0
    TSTEPS2RAD = Divide(2.0 * math.pi, float(ITERATIONSperDAY), 0.0)
    solarConstant = 1360.0
    solarDeclination = 0.3985 * math.sin((4.869 + (clock_Today_DayOfYear * 2.0 * math.pi / 365.25) +
                                          (0.03345 * math.sin((6.224 + (clock_Today_DayOfYear * 2.0 * math.pi / 365.25)))))))
    cD = math.sqrt(max(1.0 - (solarDeclination * solarDeclination), 0.0))
    m1 = [0.0] * (ITERATIONSperDAY + 1)
    m1Tot = 0.0
    for timestepNumber in range(1, ITERATIONSperDAY + 1):
        m1[timestepNumber] = (solarDeclination * math.sin(weather_Latitude * math.pi / 180.0) +
                              (cD * math.cos(weather_Latitude * math.pi / 180.0) * math.cos(TSTEPS2RAD * (timestepNumber - (ITERATIONSperDAY / 2.0))))) * 24.0 / ITERATIONSperDAY
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
    nu: float,
) -> Tuple[
    float,  # timeOfDaySecs
    float,  # netRadiation
    List[float],  # minSoilTemp
    List[float],  # maxSoilTemp
    float,        # boundaryLayerConductance
    float,        # maxTempYesterday
    int,          # airNode (passthrough)
    List[float],  # soilTemp
    float,        # airTemperature
    List[float],  # newTemperature
    List[float],  # thermalConductivity
    float,        # minTempYesterday
    List[float],  # aveSoilTemp
    List[float],  # morningSoilTemp
    float,        # weather_MeanT passthrough
    float,        # constantBoundaryLayerConductance passthrough
    float,        # weather_MinT passthrough
    int,          # clock_Today_DayOfYear passthrough
    float,        # weather_Radn passthrough
    float,        # weather_Latitude passthrough
    List[str],    # soilConstituentNames passthrough
    int,          # numNodes passthrough
    List[float],  # volSpecHeatSoil
    List[float],  # soilWater
    List[float],  # nodeDepth
    List[float],  # thickness
    int,          # surfaceNode
    float,        # MissingValue
    List[float],  # carbon
    List[float],  # bulkDensity
    float,        # pom
    List[float],  # rocks
    List[float],  # sand
    float,        # ps
    List[float],  # silt
    List[float],  # clay
    float,        # defaultTimeOfMaximumTemperature
    float,        # waterBalance_Eo
    float,        # waterBalance_Eos
    float,        # waterBalance_Salb
    float,        # stefanBoltzmannConstant
    float,        # weather_AirPressure
    float,        # weather_Wind
    float,        # instrumentHeight
    float,        # canopyHeight
    List[float],  # heatStorage
    str,          # netRadiationSource
    float,        # latentHeatOfVapourisation
    float,        # waterBalance_Es
    List[float],  # thermalConductance
    float,        # nu
    float,        # internalTimeStep
    float,        # boundaryLayerConductance (daily avg at end)
]:
    interactionsPerDay = 48
    solarRadn, cloudFr, cva = doNetRadiation(interactionsPerDay, weather_MinT, clock_Today_DayOfYear, weather_Radn, weather_Latitude)

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
            airTemperature = interpolateTemperature(
                timeOfDaySecs / 3600.0,
                minTempYesterday,
                maxTempYesterday,
                weather_MeanT,
                weather_MaxT,
                weather_MinT,
                defaultTimeOfMaximumTemperature,
            )
        newTemperature[airNode] = airTemperature
        netRadiation = interpolateNetRadiation(
            solarRadn[timeStepIteration],
            cloudFr,
            cva,
            waterBalance_Eo,
            waterBalance_Eos,
            waterBalance_Salb,
            soilTemp,
            airTemperature,
            surfaceNode,
            internalTimeStep,
            stefanBoltzmannConstant,
        )
        if boundarLayerConductanceSource == 'constant':
            thermalConductivity[airNode] = constantBoundaryLayerConductance
        elif boundarLayerConductanceSource == 'calc':
            thermalConductivity[airNode] = getBoundaryLayerConductance(
                newTemperature,
                weather_AirPressure,
                stefanBoltzmannConstant,
                waterBalance_Eos,
                weather_Wind,
                airTemperature,
                surfaceNode,
                waterBalance_Eo,
                instrumentHeight,
                canopyHeight,
            )
            for _ in range(numIterationsForBoundaryLayerConductance):
                newTemperature, heatStorage, thermalConductance, thermalConductivity = doThomas(
                    newTemperature, netRadiation, heatStorage, waterBalance_Eos, numNodes, timestep, netRadiationSource,
                    latentHeatOfVapourisation, nodeDepth, waterBalance_Es, airNode, soilTemp, surfaceNode, internalTimeStep,
                    thermalConductance, thermalConductivity, nu, volSpecHeatSoil
                )
                thermalConductivity[airNode] = getBoundaryLayerConductance(
                    newTemperature, weather_AirPressure, stefanBoltzmannConstant, waterBalance_Eos, weather_Wind,
                    airTemperature, surfaceNode, waterBalance_Eo, instrumentHeight, canopyHeight
                )

        newTemperature, heatStorage, thermalConductance, thermalConductivity = doThomas(
            newTemperature, netRadiation, heatStorage, waterBalance_Eos, numNodes, timestep, netRadiationSource,
            latentHeatOfVapourisation, nodeDepth, waterBalance_Es, airNode, soilTemp, surfaceNode, internalTimeStep,
            thermalConductance, thermalConductivity, nu, volSpecHeatSoil
        )

        boundaryLayerConductance, minSoilTemp, maxSoilTemp, aveSoilTemp, _ = doUpdate(
            interactionsPerDay, timeOfDaySecs, boundaryLayerConductance, minSoilTemp, airNode, soilTemp,
            newTemperature, numNodes, surfaceNode, internalTimeStep, maxSoilTemp, aveSoilTemp, thermalConductivity
        )

        if abs(timeOfDaySecs - (5.0 * 3600.0)) <= (min(timeOfDaySecs, 5.0 * 3600.0) * 0.0001):
            morningSoilTemp = soilTemp.copy()

    minTempYesterday = weather_MinT
    maxTempYesterday = weather_MaxT

    return (
        timeOfDaySecs,
        netRadiation,
        minSoilTemp,
        maxSoilTemp,
        boundaryLayerConductance,
        maxTempYesterday,
        airNode,
        soilTemp,
        airTemperature,
        newTemperature,
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
        ps,
        silt,
        clay,
        defaultTimeOfMaximumTemperature,
        waterBalance_Eo,
        waterBalance_Eos,
        waterBalance_Salb,
        stefanBoltzmannConstant,
        weather_AirPressure,
        weather_Wind,
        instrumentHeight,
        canopyHeight,
        heatStorage,
        netRadiationSource,
        latentHeatOfVapourisation,
        waterBalance_Es,
        thermalConductance,
        nu,
        internalTimeStep,
        boundaryLayerConductance,
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
    weather_Latitude: float,
) -> List[float]:
    soilTemp = [0.0] * (numNodes + 2)
    cumulativeDepth = ToCumThickness(thickness)
    w = 2 * math.pi / (365.25 * 24 * 3600)
    dh = 0.6
    zd = math.sqrt(2 * dh / w)
    offset = 0.25
    if weather_Latitude > 0.0:
        offset = -0.25
    for nodes in range(1, numNodes):
        soilTemp[nodes] = weather_Tav + (weather_Amp * math.exp((-1) * cumulativeDepth[nodes] / zd) *
                                         math.sin(((clock_Today_DayOfYear / 365.0 + offset) * 2.0 * math.pi - (cumulativeDepth[nodes] / zd))))
    for i in range(numNodes):
        if surfaceNode + i < len(soilTemp) and i < len(soilTemp):
            soilTemp[surfaceNode + i] = soilTemp[i]
    return soilTemp


def calcSurfaceTemperature(weather_MeanT: float, weather_MaxT: float, waterBalance_Salb: float, weather_Radn: float) -> float:
    return (1.0 - waterBalance_Salb) * (weather_MeanT + ((weather_MaxT - weather_MeanT) * math.sqrt(max(weather_Radn, 0.1) * 23.8846 / 800.0))) + (waterBalance_Salb * weather_MeanT)


def ValuesInArray(Values: Optional[List[float]], MissingValue: float) -> bool:
    if Values is not None:
        for Value in Values:
            if Value != MissingValue and not (Value != Value):  # not NaN
                return True
    return False
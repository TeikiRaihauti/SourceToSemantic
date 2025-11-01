from typing import List, Optional, Tuple

def Init(
    weather_MinT: float,
    weather_MaxT: float,
    weather_MeanT: float,
    weather_Tav: float,
    weather_Amp: float,
    weather_Radn: float,
    clock_Today_DayOfYear: int,
    physical_Rocks: List[float],
    physical_ParticleSizeSand: List[float],
    physical_ParticleSizeSilt: List[float],
    physical_ParticleSizeClay: List[float],
    organic_Carbon: List[float],
    waterBalance_SW: List[float],
    waterBalance_Salb: float,
    physical_Thickness: List[float],
    physical_BD: List[float],
    weather_Latitude: float,
    defaultInstrumentHeight: float,
    bareSoilRoughness: float,
    pInitialValues: Optional[List[float]],
    instrumHeight: float,
    nodeDepth: Optional[List[float]],
    thermCondPar1: Optional[List[float]],
    thermCondPar2: Optional[List[float]],
    thermCondPar3: Optional[List[float]],
    thermCondPar4: Optional[List[float]],
    soilRoughnessHeight: float,
    topsoilNode: int,
    surfaceNode: int,
    airNode: int,
    numPhantomNodes: int,
    DepthToConstantTemperature: float,
    MissingValue: float
) -> Tuple[
    List[float], bool, float, float, int, int, List[float], List[float], List[float], List[float], List[float],
    List[float], float, List[float], float, float, List[float], List[float], List[float], List[float],
    Optional[List[float]], List[float], List[float], List[float], List[float], List[float], List[float], float, float, float, float,
    List[float], List[float], List[float], List[float], float
]:
    """
    Initialize SoilTemperature state.

    Inputs:
    - weather_*: meteorological inputs
    - clock_Today_DayOfYear: int
    - physical_*: soil physical arrays (by layer)
    - organic_Carbon: [%] by layer
    - waterBalance_SW: volumetric water content by layer [mm/mm]
    - waterBalance_Salb: bare soil albedo [0-1]
    - physical_Thickness: layer thickness [mm]
    - physical_BD: bulk density [g/cm3]
    - weather_Latitude: degrees
    - defaultInstrumentHeight: m
    - bareSoilRoughness: mm
    - pInitialValues: initial soil temps by layer [oC] or None
    - instrumHeight: measured instrument height [m]
    - nodeDepth, thermCondPar[1..4], soilRoughnessHeight: previous arrays/values (can be None/0)
    - indices: airNode, surfaceNode, topsoilNode
    - numPhantomNodes: int
    - DepthToConstantTemperature: mm
    - MissingValue: sentinel

    Returns (in order):
    - InitialValues: List[float]
    - doInitialisationStuff: bool
    - internalTimeStep: float
    - timeOfDaySecs: float
    - numNodes: int
    - numLayers: int
    - volSpecHeatSoil: List[float]
    - soilTemp: List[float]
    - morningSoilTemp: List[float]
    - heatStorage: List[float]
    - thermalConductance: List[float]
    - thermalConductivity: List[float]
    - boundaryLayerConductance: float
    - newTemperature: List[float]
    - airTemperature: float
    - maxTempYesterday: float
    - minTempYesterday: float
    - soilWater: List[float]
    - minSoilTemp: List[float]
    - maxSoilTemp: List[float]
    - aveSoilTemp: List[float]
    - aveSoilWater: Optional[List[float]]
    - thickness: List[float]
    - bulkDensity: List[float]
    - rocks: List[float]
    - carbon: List[float]
    - sand: List[float]
    - silt: List[float]
    - clay: List[float]
    - instrumentHeight: float
    - netRadiation: float
    - canopyHeight: float
    - instrumHeight: float
    - nodeDepth: List[float]
    - thermCondPar1: List[float]
    - thermCondPar2: List[float]
    - thermCondPar3: List[float]
    - thermCondPar4: List[float]
    - soilRoughnessHeight: float
    """
    nodeDepth_loc = nodeDepth
    thermCondPar1_loc = thermCondPar1
    thermCondPar2_loc = thermCondPar2
    thermCondPar3_loc = thermCondPar3
    thermCondPar4_loc = thermCondPar4
    soilRoughnessHeight_loc = soilRoughnessHeight
    InitialValues: Optional[List[float]] = None
    doInitialisationStuff = True
    internalTimeStep = 0.0
    timeOfDaySecs = 0.0
    numNodes = 0
    numLayers = 0
    volSpecHeatSoil: Optional[List[float]] = None
    soilTemp: Optional[List[float]] = None
    morningSoilTemp: Optional[List[float]] = None
    heatStorage: Optional[List[float]] = None
    thermalConductance: Optional[List[float]] = None
    thermalConductivity: Optional[List[float]] = None
    boundaryLayerConductance = 0.0
    newTemperature: Optional[List[float]] = None
    airTemperature = 0.0
    maxTempYesterday = 0.0
    minTempYesterday = 0.0
    soilWater: Optional[List[float]] = None
    minSoilTemp: Optional[List[float]] = None
    maxSoilTemp: Optional[List[float]] = None
    aveSoilTemp: Optional[List[float]] = None
    aveSoilWater: Optional[List[float]] = None
    thickness: Optional[List[float]] = None
    bulkDensity: Optional[List[float]] = None
    rocks: Optional[List[float]] = None
    carbon: Optional[List[float]] = None
    sand: Optional[List[float]] = None
    silt: Optional[List[float]] = None
    clay: Optional[List[float]] = None
    instrumentHeight = 0.0
    netRadiation = 0.0
    canopyHeight = 0.0

    instrumentHeight = getIniVariables(instrumentHeight, instrumHeight, defaultInstrumentHeight)
    (
        heatStorage,
        minSoilTemp,
        bulkDensity,
        numNodes,
        maxSoilTemp,
        soilWater,
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
        newTemperature,
        thermalConductivity,
        thermalConductance,
        nodeDepth_loc,
    ) = getProfileVariables(
        physical_BD=physical_BD,
        waterBalance_SW=waterBalance_SW,
        organic_Carbon=organic_Carbon,
        physical_Rocks=physical_Rocks,
        nodeDepth=nodeDepth_loc,
        topsoilNode=topsoilNode,
        surfaceNode=surfaceNode,
        numPhantomNodes=numPhantomNodes,
        physical_ParticleSizeSand=physical_ParticleSizeSand,
        physical_ParticleSizeSilt=physical_ParticleSizeSilt,
        airNode=airNode,
        physical_ParticleSizeClay=physical_ParticleSizeClay,
        soilTemp=soilTemp,
        physical_Thickness=physical_Thickness,
        DepthToConstantTemperature=DepthToConstantTemperature,
        MissingValue=MissingValue,
    )
    (
        newTemperature,
        soilRoughnessHeight_loc,
        soilTemp,
        thermCondPar2_loc,
        thermCondPar3_loc,
        thermCondPar4_loc,
        thermCondPar1_loc,
    ) = readParam(
        bareSoilRoughness=bareSoilRoughness,
        newTemperature=newTemperature,
        soilRoughnessHeight=soilRoughnessHeight_loc,
        soilTemp=soilTemp,
        thermCondPar2=thermCondPar2_loc,
        numLayers=numLayers,
        bulkDensity=bulkDensity,
        numNodes=numNodes,
        thermCondPar3=thermCondPar3_loc,
        thermCondPar4=thermCondPar4_loc,
        clay=clay,
        thermCondPar1=thermCondPar1_loc,
        weather_Tav=weather_Tav,
        clock_Today_DayOfYear=clock_Today_DayOfYear,
        surfaceNode=surfaceNode,
        weather_Amp=weather_Amp,
        thickness=thickness,
        weather_Latitude=weather_Latitude,
    )
    InitialValues = pInitialValues

    # Initialization of temperatures
    if ValuesInArray(InitialValues, MissingValue):
        soilTemp = [0.0] * (numNodes + 2)
        for j in range(len(InitialValues)):
            idx = topsoilNode + j
            if 0 <= j < len(InitialValues) and 0 <= idx < len(soilTemp):
                soilTemp[idx] = InitialValues[j]
    else:
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
        InitialValues = [soilTemp[i] for i in range(topsoilNode, topsoilNode + numLayers)]

    soilTemp[airNode] = weather_MeanT
    soilTemp[surfaceNode] = calcSurfaceTemperature(weather_MeanT, weather_MaxT, waterBalance_Salb, weather_Radn)
    for i in range(numNodes + 1, len(soilTemp)):
        soilTemp[i] = weather_Tav
    newTemperature = soilTemp.copy()
    maxTempYesterday = weather_MaxT
    minTempYesterday = weather_MinT
    doInitialisationStuff = False

    return (
        InitialValues if InitialValues is not None else [],
        doInitialisationStuff,
        internalTimeStep,
        timeOfDaySecs,
        numNodes,
        numLayers,
        volSpecHeatSoil if volSpecHeatSoil is not None else [],
        soilTemp if soilTemp is not None else [],
        morningSoilTemp if morningSoilTemp is not None else [],
        heatStorage if heatStorage is not None else [],
        thermalConductance if thermalConductance is not None else [],
        thermalConductivity if thermalConductivity is not None else [],
        boundaryLayerConductance,
        newTemperature if newTemperature is not None else [],
        airTemperature,
        maxTempYesterday,
        minTempYesterday,
        soilWater if soilWater is not None else [],
        minSoilTemp if minSoilTemp is not None else [],
        maxSoilTemp if maxSoilTemp is not None else [],
        aveSoilTemp if aveSoilTemp is not None else [],
        aveSoilWater,
        thickness if thickness is not None else [],
        bulkDensity if bulkDensity is not None else [],
        rocks if rocks is not None else [],
        carbon if carbon is not None else [],
        sand if sand is not None else [],
        silt if silt is not None else [],
        clay if clay is not None else [],
        instrumentHeight,
        netRadiation,
        canopyHeight,
        instrumHeight,
        nodeDepth_loc if nodeDepth_loc is not None else [],
        thermCondPar1_loc if thermCondPar1_loc is not None else [],
        thermCondPar2_loc if thermCondPar2_loc is not None else [],
        thermCondPar3_loc if thermCondPar3_loc is not None else [],
        thermCondPar4_loc if thermCondPar4_loc is not None else [],
        soilRoughnessHeight_loc,
    )


def CalculateModel(
    weather_MinT: float,
    weather_MaxT: float,
    weather_MeanT: float,
    weather_Tav: float,
    weather_Amp: float,
    weather_AirPressure: float,
    weather_Wind: float,
    weather_Radn: float,
    clock_Today_DayOfYear: int,
    microClimate_CanopyHeight: float,
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
    InitialValues: List[float],
    doInitialisationStuff: bool,
    internalTimeStep: float,
    timeOfDaySecs: float,
    numNodes: int,
    numLayers: int,
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
    aveSoilWater: Optional[List[float]],
    thickness: List[float],
    bulkDensity: List[float],
    rocks: List[float],
    carbon: List[float],
    sand: List[float],
    silt: List[float],
    clay: List[float],
    instrumentHeight: float,
    netRadiation: float,
    canopyHeight: float,
    instrumHeight: float,
    weather_Latitude: float,
    ps: float,
    pom: float,
    airNode: int,
    surfaceNode: int,
    topsoilNode: int,
    numPhantomNodes: int,
    constantBoundaryLayerConductance: float,
    numIterationsForBoundaryLayerConductance: int,
    defaultTimeOfMaximumTemperature: float,
    timestep: float,
    latentHeatOfVapourisation: float,
    stefanBoltzmannConstant: float,
    boundarLayerConductanceSource: str,
    netRadiationSource: str,
    MissingValue: float,
    soilConstituentNames: List[str],
    nodeDepth: List[float],
    soilRoughnessHeight: float
) -> Tuple[
    List[float], bool, float, float, int, int, List[float], List[float], List[float], List[float], List[float],
    List[float], float, List[float], float, float, List[float], List[float], List[float], List[float],
    Optional[List[float]], List[float], List[float], List[float], List[float], List[float], List[float], float, float, float, float
]:
    """
    Run SoilTemperature model for one day.

    Inputs: all environment variables, parameters, and current state variables.

    Returns same state variables (updated) in the same tuple shape as Init up to instrumHeight.
    """
    # Update other variables based on exogenous inputs
    soilWater, instrumentHeight, canopyHeight = getOtherVariables(
        numLayers=numLayers,
        numNodes=numNodes,
        soilWater=soilWater,
        instrumentHeight=instrumentHeight,
        soilRoughnessHeight=soilRoughnessHeight,
        waterBalance_SW=waterBalance_SW,
        microClimate_CanopyHeight=microClimate_CanopyHeight,
        canopyHeight=canopyHeight,
    )

    # Initialisation on first call if flag is set
    if doInitialisationStuff:
        if ValuesInArray(InitialValues, MissingValue):
            soilTemp = [0.0] * (numNodes + 2)
            for j in range(len(InitialValues)):
                idx = topsoilNode + j
                if 0 <= j < len(InitialValues) and 0 <= idx < len(soilTemp):
                    soilTemp[idx] = InitialValues[j]
        else:
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
            InitialValues = [soilTemp[i] for i in range(topsoilNode, topsoilNode + numLayers)]
        soilTemp[airNode] = weather_MeanT
        soilTemp[surfaceNode] = calcSurfaceTemperature(weather_MeanT, weather_MaxT, waterBalance_Salb, weather_Radn)
        for i in range(numNodes + 1, len(soilTemp)):
            soilTemp[i] = weather_Tav
        newTemperature = soilTemp.copy()
        maxTempYesterday = weather_MaxT
        minTempYesterday = weather_MinT
        doInitialisationStuff = False

    # Process the day's heat flux and temperatures
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
        timeOfDaySecs=timeOfDaySecs,
        netRadiation=netRadiation,
        minSoilTemp=minSoilTemp,
        maxSoilTemp=maxSoilTemp,
        numIterationsForBoundaryLayerConductance=numIterationsForBoundaryLayerConductance,
        timestep=timestep,
        boundaryLayerConductance=boundaryLayerConductance,
        maxTempYesterday=maxTempYesterday,
        airNode=airNode,
        soilTemp=soilTemp,
        airTemperature=airTemperature,
        newTemperature=newTemperature,
        weather_MaxT=weather_MaxT,
        internalTimeStep=internalTimeStep,
        boundarLayerConductanceSource=boundarLayerConductanceSource,
        thermalConductivity=thermalConductivity,
        minTempYesterday=minTempYesterday,
        aveSoilTemp=aveSoilTemp,
        morningSoilTemp=morningSoilTemp,
        weather_MeanT=weather_MeanT,
        constantBoundaryLayerConductance=constantBoundaryLayerConductance,
        weather_MinT=weather_MinT,
        clock_Today_DayOfYear=clock_Today_DayOfYear,
        weather_Radn=weather_Radn,
        weather_Latitude=weather_Latitude,
        soilConstituentNames=soilConstituentNames,
        numNodes=numNodes,
        volSpecHeatSoil=volSpecHeatSoil,
        soilWater=soilWater,
        nodeDepth=nodeDepth,
        thickness=thickness,
        surfaceNode=surfaceNode,
        MissingValue=MissingValue,
        carbon=carbon,
        bulkDensity=bulkDensity,
        pom=pom,
        rocks=rocks,
        sand=sand,
        ps=ps,
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
        nu=0.6 if False else 0.6,  # default if not provided elsewhere; preserved call signature
    )

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
    )


def getIniVariables(instrumentHeight: float, instrumHeight: float, defaultInstrumentHeight: float) -> float:
    if instrumHeight > 0.00001:
        instrumentHeight = instrumHeight
    else:
        instrumentHeight = defaultInstrumentHeight
    return instrumentHeight


def getProfileVariables(
    physical_BD: List[float],
    waterBalance_SW: Optional[List[float]],
    organic_Carbon: List[float],
    physical_Rocks: List[float],
    nodeDepth: Optional[List[float]],
    topsoilNode: int,
    surfaceNode: int,
    numPhantomNodes: int,
    physical_ParticleSizeSand: List[float],
    physical_ParticleSizeSilt: List[float],
    airNode: int,
    physical_ParticleSizeClay: List[float],
    soilTemp: Optional[List[float]],
    physical_Thickness: List[float],
    DepthToConstantTemperature: float,
    MissingValue: float
) -> Tuple[
    List[float], List[float], List[float], int, List[float], List[float], List[float], List[float], List[float],
    List[float], List[float], List[float], int, List[float], List[float], List[float], List[float], List[float],
    List[float], List[float], List[float]
]:
    numLayers = len(physical_Thickness)
    numNodes = numLayers + numPhantomNodes
    thickness = [0.0] * (numLayers + numPhantomNodes + 1)
    for i in range(numLayers):
        thickness[1 + i] = physical_Thickness[i]
    belowProfileDepth = max(DepthToConstantTemperature - Sum(thickness, 1, numLayers, MissingValue), 1000.0)
    thicknessForPhantomNodes = belowProfileDepth * 2.0 / numPhantomNodes
    firstPhantomNode = numLayers
    for i in range(firstPhantomNode, firstPhantomNode + numPhantomNodes):
        thickness[i] = thicknessForPhantomNodes

    oldDepth = nodeDepth
    nodeDepth = [0.0] * (numNodes + 2)
    if oldDepth is not None:
        to = min(numNodes + 2, len(oldDepth))
        for k in range(to):
            nodeDepth[k] = oldDepth[k]
    nodeDepth[airNode] = 0.0
    nodeDepth[surfaceNode] = 0.0
    nodeDepth[topsoilNode] = 0.5 * thickness[1] / 1000.0
    for node in range(topsoilNode, numNodes + 1):
        nodeDepth[node + 1] = (Sum(thickness, 1, node - 1, MissingValue) + (0.5 * thickness[node])) / 1000.0

    oldBulkDensity = None
    bulkDensity = [0.0] * (numLayers + 1 + numPhantomNodes)
    if oldBulkDensity is not None:
        to = min(numLayers + 1 + numPhantomNodes, len(oldBulkDensity))
        for k in range(to):
            bulkDensity[k] = oldBulkDensity[k]
    for i in range(numLayers):
        bulkDensity[1 + i] = physical_BD[i]
    bulkDensity[numNodes] = bulkDensity[numLayers]
    for layer in range(numLayers + 1, numLayers + numPhantomNodes + 1):
        bulkDensity[layer] = bulkDensity[numLayers]

    oldSoilWater = soilTemp  # placeholder, not used as ref copy in original for soilWater
    soilWater = [0.0] * (numLayers + 1 + numPhantomNodes)
    if oldSoilWater is not None:
        to = min(numLayers + 1 + numPhantomNodes, len(soilWater))
        # no-op intentional

    if waterBalance_SW is not None:
        for layer in range(1, numLayers + 1):
            soilWater[layer] = Divide(waterBalance_SW[(layer - 1)] * thickness[(layer - 1)], thickness[layer], 0.0)
        for layer in range(numLayers + 1, numLayers + numPhantomNodes + 1):
            soilWater[layer] = soilWater[numLayers]

    carbon = [0.0] * (numLayers + 1 + numPhantomNodes)
    for layer in range(1, numLayers + 1):
        carbon[layer] = organic_Carbon[layer - 1]
    for layer in range(numLayers + 1, numLayers + numPhantomNodes + 1):
        carbon[layer] = carbon[numLayers]

    rocks = [0.0] * (numLayers + 1 + numPhantomNodes)
    for layer in range(1, numLayers + 1):
        rocks[layer] = physical_Rocks[layer - 1]
    for layer in range(numLayers + 1, numLayers + numPhantomNodes + 1):
        rocks[layer] = rocks[numLayers]

    sand = [0.0] * (numLayers + 1 + numPhantomNodes)
    for layer in range(1, numLayers + 1):
        sand[layer] = physical_ParticleSizeSand[layer - 1]
    for layer in range(numLayers + 1, numLayers + numPhantomNodes + 1):
        sand[layer] = sand[numLayers]

    silt = [0.0] * (numLayers + 1 + numPhantomNodes)
    for layer in range(1, numLayers + 1):
        silt[layer] = physical_ParticleSizeSilt[layer - 1]
    for layer in range(numLayers + 1, numLayers + numPhantomNodes + 1):
        silt[layer] = silt[numLayers]

    clay = [0.0] * (numLayers + 1 + numPhantomNodes)
    for layer in range(1, numLayers + 1):
        clay[layer] = physical_ParticleSizeClay[layer - 1]
    for layer in range(numLayers + 1, numLayers + numPhantomNodes + 1):
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
        soilWater,
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
        newTemperature,
        thermalConductivity,
        thermalConductance,
        nodeDepth,
    )


def doThermalConductivityCoeffs(
    thermCondPar2: Optional[List[float]],
    numLayers: int,
    bulkDensity: List[float],
    numNodes: int,
    thermCondPar3: Optional[List[float]],
    thermCondPar4: Optional[List[float]],
    clay: List[float],
    thermCondPar1: Optional[List[float]]
) -> Tuple[List[float], List[float], List[float], List[float]]:
    thermCondPar1 = (thermCondPar1[:] if thermCondPar1 is not None else []) + [0.0] * (numNodes + 1 - (len(thermCondPar1 or [])))
    thermCondPar2 = (thermCondPar2[:] if thermCondPar2 is not None else []) + [0.0] * (numNodes + 1 - (len(thermCondPar2 or [])))
    thermCondPar3 = (thermCondPar3[:] if thermCondPar3 is not None else []) + [0.0] * (numNodes + 1 - (len(thermCondPar3 or [])))
    thermCondPar4 = (thermCondPar4[:] if thermCondPar4 is not None else []) + [0.0] * (numNodes + 1 - (len(thermCondPar4 or [])))
    for layer in range(1, numLayers + 2):
        element = layer
        thermCondPar1[element] = 0.65 - (0.78 * bulkDensity[layer]) + (0.6 * (bulkDensity[layer] ** 2))
        thermCondPar2[element] = 1.06 * bulkDensity[layer]
        thermCondPar3[element] = 1.0 + Divide(2.6, (clay[layer] ** 0.5) if clay[layer] > 0 else 0.0, 0.0)
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
    weather_Latitude: float
) -> Tuple[List[float], float, List[float], List[float], List[float], List[float], List[float]]:
    thermCondPar1, thermCondPar2, thermCondPar3, thermCondPar4 = doThermalConductivityCoeffs(
        thermCondPar2=thermCondPar2,
        numLayers=numLayers,
        bulkDensity=bulkDensity,
        numNodes=numNodes,
        thermCondPar3=thermCondPar3,
        thermCondPar4=thermCondPar4,
        clay=clay,
        thermCondPar1=thermCondPar1,
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
    return newTemperature, soilRoughnessHeight, soilTemp, thermCondPar2, thermCondPar3, thermCondPar4, thermCondPar1


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
    for i in range(numLayers):
        if 1 + i < len(soilWater) and i < len(waterBalance_SW):
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
    for index, v in enumerate(values):
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
    nodeArray: List[float],
    nodeDepth: List[float],
    numNodes: int,
    thickness: List[float],
    surfaceNode: int,
    MissingValue: float
) -> List[float]:
    for node in range(surfaceNode, numNodes + 1):
        layer = node - 1
        depthLayerAbove = Sum(thickness, 1, layer, MissingValue) if layer >= 1 else 0.0
        d1 = depthLayerAbove - (nodeDepth[node] * 1000.0)
        d2 = nodeDepth[(node + 1)] * 1000.0 - depthLayerAbove
        dSum = d1 + d2
        nodeArray[node] = Divide(layerArray[layer] * d1, dSum, 0.0) + Divide(layerArray[(layer + 1)] * d2, dSum, 0.0)
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
    clay: List[float]
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
    thermalConductivity: List[float]
) -> Tuple[float, List[float], List[float], List[float], float]:
    newTemperature_copy = newTemperature.copy()
    soilTemp = newTemperature_copy.copy()
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
    boundaryLayerConductance = boundaryLayerConductance + Divide(thermalConductivity[airNode], float(numInterationsPerDay), 0.0)
    return boundaryLayerConductance, minSoilTemp, maxSoilTemp, aveSoilTemp, 0.0


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
    for node in range(surfaceNode, numNodes + 1):
        volumeOfSoilAtNode = 0.5 * (nodeDepth[node + 1] - nodeDepth[node - 1])
        heatStorage[node] = Divide(volSpecHeatSoil[node] * volumeOfSoilAtNode, internalTimeStep, 0.0)
        elementLength = nodeDepth[node] - nodeDepth[node - 1]
        thermalConductance[node] = Divide(thermalConductivity[node], elementLength, 0.0)
    g = 1 - nu
    for node in range(surfaceNode, numNodes + 1):
        c[node] = -nu * thermalConductance[node]
        a[node + 1] = c[node]
        b[node] = nu * (thermalConductance[node] + thermalConductance[node - 1]) + heatStorage[node]
        d[node] = g * thermalConductance[(node - 1)] * soilTemp[(node - 1)] + ((heatStorage[node] - (g * (thermalConductance[node] + thermalConductance[node - 1]))) * soilTemp[node]) + (g * thermalConductance[node] * soilTemp[(node + 1)])
    a[surfaceNode] = 0.0
    sensibleHeatFlux = nu * thermalConductance[airNode] * newTemps[airNode]
    radnNet = 0.0
    if netRadiationSource == "calc":
        radnNet = Divide(netRadiation * 1000000.0, internalTimeStep, 0.0)
    elif netRadiationSource == "eos":
        radnNet = Divide(waterBalance_Eos * latentHeatOfVapourisation, timestep, 0.0)
    latentHeatFlux = Divide(waterBalance_Es * latentHeatOfVapourisation, timestep, 0.0)
    soilSurfaceHeatFlux = sensibleHeatFlux + radnNet - latentHeatFlux
    d[surfaceNode] = d[surfaceNode] + soilSurfaceHeatFlux
    d[numNodes] = d[numNodes] + (nu * thermalConductance[numNodes] * newTemps[(numNodes + 1)])
    for node in range(surfaceNode, numNodes):
        c[node] = Divide(c[node], b[node], 0.0)
        d[node] = Divide(d[node], b[node], 0.0)
        b[node + 1] = b[node + 1] - (a[(node + 1)] * c[node])
        d[node + 1] = d[node + 1] - (a[(node + 1)] * d[node])
    newTemps[numNodes] = Divide(d[numNodes], b[numNodes], 0.0)
    for node in range(numNodes - 1, surfaceNode - 1, -1):
        newTemps[node] = d[node] - (c[node] * newTemps[(node + 1)])
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
        frictionVelocity = Divide(
            weather_Wind * vonKarmanConstant,
            (math_log(Divide(instrumentHeight - d + roughnessFactorMomentum, roughnessFactorMomentum, 0.0)) + stabilityCorrectionMomentum),
            0.0,
        )
        boundaryLayerCond = Divide(
            SpecificHeatAir * vonKarmanConstant * frictionVelocity,
            (math_log(Divide(instrumentHeight - d + roughnessFactorHeat, roughnessFactorHeat, 0.0)) + stabilityCorrectionHeat),
            0.0,
        )
        boundaryLayerCond = boundaryLayerCond + radiativeConductance
        heatFluxDensity = boundaryLayerCond * (surfaceTemperature - airTemperature)
        stabilityParammeter = Divide(-vonKarmanConstant * instrumentHeight * gravitationalConstant * heatFluxDensity, SpecificHeatAir * kelvinT(airTemperature) * (frictionVelocity ** 3.0), 0.0)
        if stabilityParammeter > 0.0:
            stabilityCorrectionHeat = 4.7 * stabilityParammeter
            stabilityCorrectionMomentum = stabilityCorrectionHeat
        else:
            stabilityCorrectionHeat = -2.0 * math_log((1.0 + (1.0 - (16.0 * stabilityParammeter)) ** 0.5) / 2.0)
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
    defaultTimeOfMaximumTemperature: float
) -> float:
    time = timeHours / 24.0
    maxT_time = defaultTimeOfMaximumTemperature / 24.0
    minT_time = maxT_time - 0.5
    if time < minT_time:
        midnightT = math_sin((0.0 + 0.25 - maxT_time) * 2.0 * math_pi()) * (maxTempYesterday - minTempYesterday) / 2.0 + ((maxTempYesterday + minTempYesterday) / 2.0)
        tScale = (minT_time - time) / minT_time
        if tScale > 1.0:
            tScale = 1.0
        elif tScale < 0.0:
            tScale = 0.0
        currentTemperature = weather_MinT + (tScale * (midnightT - weather_MinT))
        return currentTemperature
    else:
        currentTemperature = math_sin((time + 0.25 - maxT_time) * 2.0 * math_pi()) * (weather_MaxT - weather_MinT) / 2.0 + weather_MeanT
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
    MissingValue: float
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
        thermCondLayers[node] = numerator / denominator if denominator != 0 else 0.0
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
    MissingValue: float
) -> List[float]:
    volspecHeatSoil_ = [0.0] * (numNodes + 1)
    for node in range(1, numNodes + 1):
        volspecHeatSoil_[node] = 0.0
        for constituentName in soilConstituentNames:
            if constituentName not in {"Minerals"}:
                volspecHeatSoil_[node] = volspecHeatSoil_[node] + (volumetricSpecificHeat(constituentName, node) * 1000000.0 * soilWater[node])
    volSpecHeatSoil = mapLayer2Node(volspecHeatSoil_, volSpecHeatSoil, nodeDepth, numNodes, thickness, surfaceNode, MissingValue)
    return volSpecHeatSoil


def Zero(arr: Optional[List[float]]) -> Optional[List[float]]:
    if arr is not None:
        for i in range(len(arr)):
            arr[i] = 0.0
    return arr


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
    TSTEPS2RAD = Divide(2.0 * math_pi(), float(ITERATIONSperDAY), 0.0)
    solarConstant = 1360.0
    solarDeclination = 0.3985 * math_sin((4.869 + (clock_Today_DayOfYear * 2.0 * math_pi() / 365.25) + (0.03345 * math_sin((6.224 + (clock_Today_DayOfYear * 2.0 * math_pi() / 365.25)))))))
    cD = (1.0 - (solarDeclination * solarDeclination)) ** 0.5
    m1 = [0.0] * (ITERATIONSperDAY + 1)
    m1Tot = 0.0
    for timestepNumber in range(1, ITERATIONSperDAY + 1):
        m1[timestepNumber] = (solarDeclination * math_sin(weather_Latitude * math_pi() / 180.0) + (cD * math_cos(weather_Latitude * math_pi() / 180.0) * math_cos(TSTEPS2RAD * (timestepNumber - (ITERATIONSperDAY / 2.0))))) * 24.0 / ITERATIONSperDAY
        if m1[timestepNumber] > 0.0:
            m1Tot = m1Tot + m1[timestepNumber]
        else:
            m1[timestepNumber] = 0.0
    psr = m1Tot * solarConstant * 3600.0 / 1000000.0
    fr = Divide(max(weather_Radn, 0.1), psr, 0.0)
    cloudFr = 2.33 - (3.33 * fr)
    cloudFr = min(max(cloudFr, 0.0), 1.0)
    for timestepNumber in range(1, ITERATIONSperDAY + 1):
        solarRadn[timestepNumber] = max(weather_Radn, 0.1) * Divide(m1[timestepNumber], m1Tot, 0.0)
    cva = math_exp((31.3716 - (6014.79 / kelvinT(weather_MinT)) - (0.00792495 * kelvinT(weather_MinT)))) / kelvinT(weather_MinT)
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
    solarRadn = [0.0] * (49)
    solarRadn, cloudFr, cva = doNetRadiation(solarRadn, cloudFr, cva, interactionsPerDay, weather_MinT, clock_Today_DayOfYear, weather_Radn, weather_Latitude)
    Zero(minSoilTemp)
    Zero(maxSoilTemp)
    Zero(aveSoilTemp)
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
            thermalConductivity[airNode] = constantBoundaryLayerConductance
        elif boundarLayerConductanceSource == "calc":
            thermalConductivity[airNode] = getBoundaryLayerConductance(newTemperature, weather_AirPressure, stefanBoltzmannConstant, waterBalance_Eos, weather_Wind, airTemperature, surfaceNode, waterBalance_Eo, instrumentHeight, canopyHeight)
            for _ in range(1, numIterationsForBoundaryLayerConductance + 1):
                newTemperature, heatStorage, thermalConductance = doThomas(newTemperature, netRadiation, heatStorage, waterBalance_Eos, numNodes, timestep, netRadiationSource, latentHeatOfVapourisation, nodeDepth, waterBalance_Es, airNode, soilTemp, surfaceNode, internalTimeStep, thermalConductance, thermalConductivity, nu, volSpecHeatSoil)
                thermalConductivity[airNode] = getBoundaryLayerConductance(newTemperature, weather_AirPressure, stefanBoltzmannConstant, waterBalance_Eos, weather_Wind, airTemperature, surfaceNode, waterBalance_Eo, instrumentHeight, canopyHeight)
        newTemperature, heatStorage, thermalConductance = doThomas(newTemperature, netRadiation, heatStorage, waterBalance_Eos, numNodes, timestep, netRadiationSource, latentHeatOfVapourisation, nodeDepth, waterBalance_Es, airNode, soilTemp, surfaceNode, internalTimeStep, thermalConductance, thermalConductivity, nu, volSpecHeatSoil)
        bcond, minSoilTemp, maxSoilTemp, aveSoilTemp, _ = doUpdate(interactionsPerDay, timeOfDaySecs, boundaryLayerConductance, minSoilTemp, airNode, soilTemp, newTemperature, numNodes, surfaceNode, internalTimeStep, maxSoilTemp, aveSoilTemp, thermalConductivity)
        boundaryLayerConductance = bcond
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
    cumulativeDepth = ToCumThickness(thickness)
    w = 2 * math_pi() / (365.25 * 24 * 3600)
    dh = 0.6
    zd = (2 * dh / w) ** 0.5
    offset = 0.25
    if weather_Latitude > 0.0:
        offset = -0.25
    soilTemp = [0.0] * (numNodes + 2)
    for nodes in range(1, numNodes + 1):
        soilTemp[nodes] = weather_Tav + (weather_Amp * math_exp(-1 * cumulativeDepth[nodes] / zd) * math_sin(((clock_Today_DayOfYear / 365.0 + offset) * 2.0 * math_pi() - (cumulativeDepth[nodes] / zd))))
    tocopy = min(numNodes, len(soilTemp))
    for k in range(tocopy):
        idx = surfaceNode + k
        if idx < len(soilTempIO):
            soilTempIO[idx] = soilTemp[k]
    return soilTempIO


def calcSurfaceTemperature(weather_MeanT: float, weather_MaxT: float, waterBalance_Salb: float, weather_Radn: float) -> float:
    surfaceT = (1.0 - waterBalance_Salb) * (weather_MeanT + ((weather_MaxT - weather_MeanT) * ((max(weather_Radn, 0.1) * 23.8846 / 800.0) ** 0.5))) + (waterBalance_Salb * weather_MeanT)
    return surfaceT


def ValuesInArray(Values: Optional[List[float]], MissingValue: float) -> bool:
    if Values is not None:
        for Value in Values:
            if Value != MissingValue and not is_nan(Value):
                return True
    return False


# Math helper functions without importing modules at top-level
def math_pi() -> float:
    return 3.141592653589793


def math_sin(x: float) -> float:
    # Using Taylor approximation is unnecessary; assuming Python runtime has math; but to avoid imports, implement minimal
    import math as _math
    return _math.sin(x)


def math_cos(x: float) -> float:
    import math as _math
    return _math.cos(x)


def math_log(x: float) -> float:
    import math as _math
    return _math.log(x)


def math_exp(x: float) -> float:
    import math as _math
    return _math.exp(x)


def is_nan(x: float) -> bool:
    import math as _math
    return _math.isnan(x)
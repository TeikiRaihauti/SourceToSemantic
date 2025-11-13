def Init(
    weather_Tav,  # type: float
    weather_Amp,  # type: float
    clock_Today_DayOfYear,  # type: int
    weather_Latitude,  # type: float
    physical_Thickness,  # type: list[float]
    physical_BD,  # type: list[float]
    physical_Rocks,  # type: list[float]
    physical_ParticleSizeSand,  # type: list[float]
    physical_ParticleSizeSilt,  # type: list[float]
    physical_ParticleSizeClay,  # type: list[float]
    organic_Carbon,  # type: list[float]
    waterBalance_SW,  # type: list[float]
    airNode,  # type: int
    surfaceNode,  # type: int
    topsoilNode,  # type: int
    numPhantomNodes,  # type: int
    defaultInstrumentHeight,  # type: float
    bareSoilRoughness,  # type: float
    DepthToConstantTemperature,  # type: float
    MissingValue,  # type: float
    nodeDepth=None,  # type: list[float] | None
    thermCondPar1=None,  # type: list[float] | None
    thermCondPar2=None,  # type: list[float] | None
    thermCondPar3=None,  # type: list[float] | None
    thermCondPar4=None,  # type: list[float] | None
    instrumHeight=0.0,  # type: float
    pInitialValues=None,  # type: list[float] | None
):
    """
    Initialization function. Returns the full set of initialized state variables and internal parameters.

    Returns (in order):
    - InitialValues: list[float] | None
    - doInitialisationStuff: bool
    - internalTimeStep: float
    - timeOfDaySecs: float
    - numNodes: int
    - numLayers: int
    - volSpecHeatSoil: list[float]
    - soilTemp: list[float]
    - morningSoilTemp: list[float]
    - heatStorage: list[float]
    - thermalConductance: list[float]
    - thermalConductivity: list[float]
    - boundaryLayerConductance: float
    - newTemperature: list[float]
    - airTemperature: float
    - maxTempYesterday: float
    - minTempYesterday: float
    - soilWater: list[float]
    - minSoilTemp: list[float]
    - maxSoilTemp: list[float]
    - aveSoilTemp: list[float]
    - aveSoilWater: list[float] | None
    - thickness: list[float]
    - bulkDensity: list[float]
    - rocks: list[float]
    - carbon: list[float]
    - sand: list[float]
    - silt: list[float]
    - clay: list[float]
    - instrumentHeight: float
    - netRadiation: float
    - canopyHeight: float
    - instrumHeight: float
    - nodeDepth: list[float]
    - thermCondPar1: list[float]
    - thermCondPar2: list[float]
    - thermCondPar3: list[float]
    - thermCondPar4: list[float]
    - soilRoughnessHeight: float
    """
    instrumentHeight = getIniVariables(0.0, instrumHeight, defaultInstrumentHeight)
    (
        heatStorage,
        minSoilTemp,
        bulkDensity,
        numNodes,
        maxSoilTemp,
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
        newTemperature,
        nodeDepth,
    ) = getProfileVariables(
        physical_BD=physical_BD,
        waterBalance_SW=waterBalance_SW,
        organic_Carbon=organic_Carbon,
        physical_Rocks=physical_Rocks,
        nodeDepth=nodeDepth,
        topsoilNode=topsoilNode,
        surfaceNode=surfaceNode,
        numPhantomNodes=numPhantomNodes,
        physical_ParticleSizeSand=physical_ParticleSizeSand,
        physical_ParticleSizeSilt=physical_ParticleSizeSilt,
        airNode=airNode,
        physical_ParticleSizeClay=physical_ParticleSizeClay,
        physical_Thickness=physical_Thickness,
        DepthToConstantTemperature=DepthToConstantTemperature,
        MissingValue=MissingValue,
    )
    (
        thermCondPar1,
        thermCondPar2,
        thermCondPar3,
        thermCondPar4,
        soilTemp,
        newTemperature,
        soilRoughnessHeight,
    ) = readParam(
        bareSoilRoughness=bareSoilRoughness,
        soilTemp=soilTemp,
        thermCondPar1=thermCondPar1,
        thermCondPar2=thermCondPar2,
        thermCondPar3=thermCondPar3,
        thermCondPar4=thermCondPar4,
        numLayers=numLayers,
        bulkDensity=bulkDensity,
        numNodes=numNodes,
        clay=clay,
        weather_Tav=weather_Tav,
        clock_Today_DayOfYear=clock_Today_DayOfYear,
        surfaceNode=surfaceNode,
        weather_Amp=weather_Amp,
        thickness=thickness,
        weather_Latitude=weather_Latitude,
    )
    InitialValues = pInitialValues
    doInitialisationStuff = True
    internalTimeStep = 0.0
    timeOfDaySecs = 0.0
    boundaryLayerConductance = 0.0
    airTemperature = 0.0
    maxTempYesterday = 0.0
    minTempYesterday = 0.0
    aveSoilWater = None
    netRadiation = 0.0
    canopyHeight = 0.0
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
        soilRoughnessHeight,
    )


def CalculateModel(
    weather_MinT,  # type: float
    weather_MaxT,  # type: float
    weather_MeanT,  # type: float
    weather_Tav,  # type: float
    weather_Amp,  # type: float
    weather_AirPressure,  # type: float
    weather_Wind,  # type: float
    weather_Radn,  # type: float
    clock_Today_DayOfYear,  # type: int
    microClimate_CanopyHeight,  # type: float
    physical_Rocks,  # type: list[float]
    physical_ParticleSizeSand,  # type: list[float]
    physical_ParticleSizeSilt,  # type: list[float]
    physical_ParticleSizeClay,  # type: list[float]
    organic_Carbon,  # type: list[float]
    waterBalance_SW,  # type: list[float]
    waterBalance_Eos,  # type: float
    waterBalance_Eo,  # type: float
    waterBalance_Es,  # type: float
    waterBalance_Salb,  # type: float
    weather_Latitude,  # type: float
    physical_Thickness,  # type: list[float]
    physical_BD,  # type: list[float]
    ps,  # type: float
    pInitialValues,  # type: list[float] | None
    DepthToConstantTemperature,  # type: float
    timestep,  # type: float
    latentHeatOfVapourisation,  # type: float
    stefanBoltzmannConstant,  # type: float
    airNode,  # type: int
    surfaceNode,  # type: int
    topsoilNode,  # type: int
    numPhantomNodes,  # type: int
    constantBoundaryLayerConductance,  # type: float
    numIterationsForBoundaryLayerConductance,  # type: int
    defaultTimeOfMaximumTemperature,  # type: float
    defaultInstrumentHeight,  # type: float
    bareSoilRoughness,  # type: float
    pom,  # type: float
    nu,  # type: float
    boundarLayerConductanceSource,  # type: str
    netRadiationSource,  # type: str
    MissingValue,  # type: float
    soilConstituentNames,  # type: list[str]
    # state inputs:
    InitialValues,  # type: list[float] | None
    doInitialisationStuff,  # type: bool
    internalTimeStep,  # type: float
    timeOfDaySecs,  # type: float
    numNodes,  # type: int
    numLayers,  # type: int
    volSpecHeatSoil,  # type: list[float]
    soilTemp,  # type: list[float]
    morningSoilTemp,  # type: list[float]
    heatStorage,  # type: list[float]
    thermalConductance,  # type: list[float]
    thermalConductivity,  # type: list[float]
    boundaryLayerConductance,  # type: float
    newTemperature,  # type: list[float]
    airTemperature,  # type: float
    maxTempYesterday,  # type: float
    minTempYesterday,  # type: float
    soilWater,  # type: list[float]
    minSoilTemp,  # type: list[float]
    maxSoilTemp,  # type: list[float]
    aveSoilTemp,  # type: list[float]
    aveSoilWater,  # type: list[float] | None
    thickness,  # type: list[float]
    bulkDensity,  # type: list[float]
    rocks,  # type: list[float]
    carbon,  # type: list[float]
    sand,  # type: list[float]
    silt,  # type: list[float]
    clay,  # type: list[float]
    instrumentHeight,  # type: float
    netRadiation,  # type: float
    canopyHeight,  # type: float
    instrumHeight,  # type: float
    nodeDepth,  # type: list[float]
    thermCondPar1,  # type: list[float]
    thermCondPar2,  # type: list[float]
    thermCondPar3,  # type: list[float]
    thermCondPar4,  # type: list[float]
    soilRoughnessHeight,  # type: float
):
    """
    Main biophysical process function. Returns updated state variables.

    Returns (in order):
    - heatStorage: list[float]
    - instrumentHeight: float
    - canopyHeight: float
    - minSoilTemp: list[float]
    - maxSoilTemp: list[float]
    - aveSoilTemp: list[float]
    - volSpecHeatSoil: list[float]
    - soilTemp: list[float]
    - morningSoilTemp: list[float]
    - newTemperature: list[float]
    - thermalConductivity: list[float]
    - thermalConductance: list[float]
    - boundaryLayerConductance: float
    - soilWater: list[float]
    - doInitialisationStuff: bool
    - maxTempYesterday: float
    - minTempYesterday: float
    - airTemperature: float
    - internalTimeStep: float
    - timeOfDaySecs: float
    - netRadiation: float
    - InitialValues: list[float] | None
    """
    # Update soil water, canopy and instrument heights
    soilWater, instrumentHeight, canopyHeight = getOtherVariables(
        numLayers, numNodes, soilWater, instrumentHeight, soilRoughnessHeight, waterBalance_SW, microClimate_CanopyHeight, canopyHeight
    )
    if doInitialisationStuff:
        if ValuesInArray(InitialValues, MissingValue):
            soilTemp = [0.0] * (numNodes + 2)
            # copy initial values into nodes starting at topsoilNode
            if InitialValues is not None:
                ncpy = min(len(InitialValues), len(soilTemp) - topsoilNode)
                soilTemp[topsoilNode:topsoilNode + ncpy] = InitialValues[:ncpy]
        else:
            soilTemp = calcSoilTemperature(soilTemp, weather_Tav, clock_Today_DayOfYear, surfaceNode, numNodes, weather_Amp, thickness, weather_Latitude)
            InitialValues = [0.0] * numLayers
            # copy topsoil layers back to InitialValues
            InitialValues[0:numLayers] = soilTemp[topsoilNode:topsoilNode + numLayers]
        soilTemp[airNode] = weather_MeanT
        soilTemp[surfaceNode] = calcSurfaceTemperature(weather_MeanT, weather_MaxT, waterBalance_Salb, weather_Radn)
        for i in range(numNodes + 1, len(soilTemp)):
            soilTemp[i] = weather_Tav
        newTemperature = list(soilTemp)
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
        nu=nu,
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


def getIniVariables(instrumentHeight, instrumHeight, defaultInstrumentHeight):
    # type: (float, float, float) -> float
    if instrumHeight > 0.00001:
        instrumentHeight = instrumHeight
    else:
        instrumentHeight = defaultInstrumentHeight
    return instrumentHeight


def getProfileVariables(
    physical_BD,  # type: list[float]
    waterBalance_SW,  # type: list[float]
    organic_Carbon,  # type: list[float]
    physical_Rocks,  # type: list[float]
    nodeDepth,  # type: list[float] | None
    topsoilNode,  # type: int
    surfaceNode,  # type: int
    numPhantomNodes,  # type: int
    physical_ParticleSizeSand,  # type: list[float]
    physical_ParticleSizeSilt,  # type: list[float]
    airNode,  # type: int
    physical_ParticleSizeClay,  # type: list[float]
    physical_Thickness,  # type: list[float]
    DepthToConstantTemperature,  # type: float
    MissingValue,  # type: float
):
    # type: (...) -> tuple[list[float], list[float], list[float], int, list[float], list[float], list[float], list[float], list[float], list[float], list[float], list[float], list[float], list[float], int, list[float], list[float], list[float], list[float], list[float], list[float]]
    import math
    numLayers = len(physical_Thickness)
    numNodes = numLayers + numPhantomNodes
    thickness = [0.0] * (numLayers + numPhantomNodes + 1)
    # Copy physical_Thickness starting at index 1
    thickness[1:1 + len(physical_Thickness)] = physical_Thickness[:]
    belowProfileDepth = max(DepthToConstantTemperature - Sum(thickness, 1, numLayers, MissingValue), 1000.0)
    thicknessForPhantomNodes = belowProfileDepth * 2.0 / numPhantomNodes
    firstPhantomNode = numLayers
    for i in range(firstPhantomNode, firstPhantomNode + numPhantomNodes):
        thickness[i] = thicknessForPhantomNodes
    oldDepth = nodeDepth
    nodeDepth = [0.0] * (numNodes + 2)
    if oldDepth is not None:
        copy_len = min(len(nodeDepth), len(oldDepth))
        nodeDepth[0:copy_len] = oldDepth[0:copy_len]
    nodeDepth[airNode] = 0.0
    nodeDepth[surfaceNode] = 0.0
    nodeDepth[topsoilNode] = 0.5 * thickness[1] / 1000.0
    for node in range(topsoilNode, numNodes + 1):
        nodeDepth[node + 1] = (Sum(thickness, 1, node - 1, MissingValue) + (0.5 * thickness[node])) / 1000.0
    oldBulkDensity = None
    bulkDensity = [0.0] * (numLayers + 1 + numPhantomNodes)
    physical_BD_copy = physical_BD[:]
    bulkDensity[1:1 + len(physical_BD_copy)] = physical_BD_copy
    bulkDensity[numNodes] = bulkDensity[numLayers]
    for layer in range(numLayers + 1, numLayers + numPhantomNodes + 1):
        bulkDensity[layer] = bulkDensity[numLayers]
    oldSoilWater = None
    soilWater = [0.0] * (numLayers + 1 + numPhantomNodes)
    if waterBalance_SW is not None:
        for layer in range(1, numLayers + 1):
            soilWater[layer] = Divide(waterBalance_SW[layer - 1] * thickness[layer - 1], thickness[layer], 0.0)
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
        newTemperature,
        nodeDepth,
    )


def doThermalConductivityCoeffs(thermCondPar2, numLayers, bulkDensity, numNodes, thermCondPar3, thermCondPar4, clay, thermCondPar1):
    # type: (list[float] | None, int, list[float], int, list[float] | None, list[float] | None, list[float], list[float] | None) -> tuple[list[float], list[float], list[float], list[float]]
    import math
    oldGC1 = thermCondPar1
    thermCondPar1 = [0.0] * (numNodes + 1)
    if oldGC1 is not None:
        copy_len = min(len(thermCondPar1), len(oldGC1))
        thermCondPar1[0:copy_len] = oldGC1[0:copy_len]
    oldGC2 = thermCondPar2
    thermCondPar2 = [0.0] * (numNodes + 1)
    if oldGC2 is not None:
        copy_len = min(len(thermCondPar2), len(oldGC2))
        thermCondPar2[0:copy_len] = oldGC2[0:copy_len]
    oldGC3 = thermCondPar3
    thermCondPar3 = [0.0] * (numNodes + 1)
    if oldGC3 is not None:
        copy_len = min(len(thermCondPar3), len(oldGC3))
        thermCondPar3[0:copy_len] = oldGC3[0:copy_len]
    oldGC4 = thermCondPar4
    thermCondPar4 = [0.0] * (numNodes + 1)
    if oldGC4 is not None:
        copy_len = min(len(thermCondPar4), len(oldGC4))
        thermCondPar4[0:copy_len] = oldGC4[0:copy_len]
    for layer in range(1, numLayers + 2):
        element = layer
        if layer < len(bulkDensity) and element < len(thermCondPar1):
            thermCondPar1[element] = 0.65 - (0.78 * bulkDensity[layer]) + (0.6 * (bulkDensity[layer] ** 2))
            thermCondPar2[element] = 1.06 * bulkDensity[layer]
            cl = clay[layer] if layer < len(clay) and clay[layer] > 0.0 else 1e-12
            thermCondPar3[element] = 1.0 + Divide(2.6, math.sqrt(cl), 0.0)
            thermCondPar4[element] = 0.03 + (0.1 * (bulkDensity[layer] ** 2))
    return thermCondPar1, thermCondPar2, thermCondPar3, thermCondPar4


def readParam(
    bareSoilRoughness,  # type: float
    soilTemp,  # type: list[float]
    thermCondPar1,  # type: list[float] | None
    thermCondPar2,  # type: list[float] | None
    thermCondPar3,  # type: list[float] | None
    thermCondPar4,  # type: list[float] | None
    numLayers,  # type: int
    bulkDensity,  # type: list[float]
    numNodes,  # type: int
    clay,  # type: list[float]
    weather_Tav,  # type: float
    clock_Today_DayOfYear,  # type: int
    surfaceNode,  # type: int
    weather_Amp,  # type: float
    thickness,  # type: list[float]
    weather_Latitude,  # type: float
):
    # type: (...) -> tuple[list[float], list[float], list[float], list[float], list[float], list[float], float]
    (thermCondPar1, thermCondPar2, thermCondPar3, thermCondPar4) = doThermalConductivityCoeffs(
        thermCondPar2, numLayers, bulkDensity, numNodes, thermCondPar3, thermCondPar4, clay, thermCondPar1
    )
    soilTemp = calcSoilTemperature(soilTemp, weather_Tav, clock_Today_DayOfYear, surfaceNode, numNodes, weather_Amp, thickness, weather_Latitude)
    newTemperature = list(soilTemp)
    soilRoughnessHeight = bareSoilRoughness
    return thermCondPar1, thermCondPar2, thermCondPar3, thermCondPar4, soilTemp, newTemperature, soilRoughnessHeight


def getOtherVariables(numLayers, numNodes, soilWater, instrumentHeight, soilRoughnessHeight, waterBalance_SW, microClimate_CanopyHeight, canopyHeight):
    # type: (int, int, list[float], float, float, list[float], float, float) -> tuple[list[float], float, float]
    soilWater[1:1 + numLayers] = waterBalance_SW[:]
    soilWater[numNodes] = soilWater[numLayers]
    canopyHeight = max(microClimate_CanopyHeight, soilRoughnessHeight) / 1000.0
    instrumentHeight = max(instrumentHeight, canopyHeight + 0.5)
    return soilWater, instrumentHeight, canopyHeight


def volumetricFractionOrganicMatter(layer, carbon, bulkDensity, pom):
    # type: (int, list[float], list[float], float) -> float
    return carbon[layer] / 100.0 * 2.5 * bulkDensity[layer] / pom


def volumetricFractionRocks(layer, rocks):
    # type: (int, list[float]) -> float
    return rocks[layer] / 100.0


def volumetricFractionIce(layer):
    # type: (int) -> float
    return 0.0


def volumetricFractionWater(layer, soilWater, carbon, bulkDensity, pom):
    # type: (int, list[float], list[float], list[float], float) -> float
    return (1 - volumetricFractionOrganicMatter(layer, carbon, bulkDensity, pom)) * soilWater[layer]


def volumetricFractionClay(layer, bulkDensity, ps, clay, carbon, pom, rocks):
    # type: (int, list[float], float, list[float], list[float], float, list[float]) -> float
    return (1 - volumetricFractionOrganicMatter(layer, carbon, bulkDensity, pom) - volumetricFractionRocks(layer, rocks)) * clay[layer] / 100.0 * bulkDensity[layer] / ps


def volumetricFractionSilt(layer, bulkDensity, silt, ps, carbon, pom, rocks):
    # type: (int, list[float], list[float], float, list[float], float, list[float]) -> float
    return (1 - volumetricFractionOrganicMatter(layer, carbon, bulkDensity, pom) - volumetricFractionRocks(layer, rocks)) * silt[layer] / 100.0 * bulkDensity[layer] / ps


def volumetricFractionSand(layer, bulkDensity, sand, ps, carbon, pom, rocks):
    # type: (int, list[float], list[float], float, list[float], float, list[float]) -> float
    return (1 - volumetricFractionOrganicMatter(layer, carbon, bulkDensity, pom) - volumetricFractionRocks(layer, rocks)) * sand[layer] / 100.0 * bulkDensity[layer] / ps


def kelvinT(celciusT):
    # type: (float) -> float
    celciusToKelvin = 273.18
    return celciusT + celciusToKelvin


def Divide(value1, value2, errVal):
    # type: (float, float, float) -> float
    if value2 != 0.0:
        return value1 / value2
    return errVal


def Sum(values, startIndex, endIndex, MissingValue):
    # type: (list[float], int, int, float) -> float
    result = 0.0
    index = -1
    for value in values:
        index += 1
        if index >= startIndex and value != MissingValue:
            result += value
        if index == endIndex:
            break
    return result


def volumetricSpecificHeat(name, layer):
    # type: (str, int) -> float
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


def volumetricFractionAir(layer, rocks, carbon, bulkDensity, pom, sand, ps, silt, clay, soilWater):
    # type: (int, list[float], list[float], list[float], float, list[float], float, list[float], list[float], list[float]) -> float
    return 1.0 - volumetricFractionRocks(layer, rocks) - volumetricFractionOrganicMatter(layer, carbon, bulkDensity, pom) - volumetricFractionSand(layer, bulkDensity, sand, ps, carbon, pom, rocks) - volumetricFractionSilt(layer, bulkDensity, silt, ps, carbon, pom, rocks) - volumetricFractionClay(layer, bulkDensity, ps, clay, carbon, pom, rocks) - volumetricFractionWater(layer, soilWater, carbon, bulkDensity, pom) - volumetricFractionIce(layer)


def airDensity(temperature, AirPressure):
    # type: (float, float) -> float
    MWair = 0.02897
    RGAS = 8.3143
    HPA2PA = 100.0
    return Divide(MWair * AirPressure * HPA2PA, kelvinT(temperature) * RGAS, 0.0)


def longWaveRadn(emissivity, tDegC, stefanBoltzmannConstant):
    # type: (float, float, float) -> float
    import math
    return stefanBoltzmannConstant * emissivity * (kelvinT(tDegC) ** 4)


def mapLayer2Node(layerArray, nodeArray, nodeDepth, numNodes, thickness, surfaceNode, MissingValue):
    # type: (list[float], list[float], list[float], int, list[float], int, float) -> list[float]
    for node in range(surfaceNode, numNodes + 1):
        layer = node - 1
        depthLayerAbove = Sum(thickness, 1, layer, MissingValue) if layer >= 1 else 0.0
        d1 = depthLayerAbove - (nodeDepth[node] * 1000.0)
        d2 = nodeDepth[(node + 1)] * 1000.0 - depthLayerAbove
        dSum = d1 + d2
        nodeArray[node] = Divide(layerArray[layer] * d1, dSum, 0.0) + Divide(layerArray[(layer + 1)] * d2, dSum, 0.0)
    return nodeArray


def ThermalConductance(name, layer, rocks, bulkDensity, sand, ps, carbon, pom, silt, clay):
    # type: (str, int, list[float], list[float], list[float], float, list[float], float, list[float], list[float]) -> float
    import math
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
    # Preserve original behavior (bug) of returning volumetricSpecificHeat value
    result = volumetricSpecificHeat(name, layer)
    return result


def shapeFactor(name, layer, soilWater, carbon, bulkDensity, pom, rocks, sand, ps, silt, clay):
    # type: (str, int, list[float], list[float], list[float], float, list[float], list[float], float, list[float], list[float]) -> float
    import math
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
    # Preserve original behavior (bug) of returning volumetricSpecificHeat value
    result = volumetricSpecificHeat(name, layer)
    return result


def doUpdate(numInterationsPerDay, timeOfDaySecs, boundaryLayerConductance, minSoilTemp, airNode, soilTemp, newTemperature, numNodes, surfaceNode, internalTimeStep, maxSoilTemp, aveSoilTemp, thermalConductivity):
    # type: (int, float, float, list[float], int, list[float], list[float], int, int, float, list[float], list[float], list[float]) -> tuple[list[float], list[float], list[float], float]
    newTemperature_copy = list(newTemperature)
    soilTemp = newTemperature_copy
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
    return soilTemp, minSoilTemp, maxSoilTemp, aveSoilTemp, boundaryLayerConductance


def doThomas(newTemps, netRadiation, heatStorage, waterBalance_Eos, numNodes, timestep, netRadiationSource, latentHeatOfVapourisation, nodeDepth, waterBalance_Es, airNode, soilTemp, surfaceNode, internalTimeStep, thermalConductance, thermalConductivity, nu, volSpecHeatSoil):
    # type: (list[float], float, list[float], float, int, float, str, float, list[float], float, int, list[float], int, float, list[float], list[float], float, list[float]) -> tuple[list[float], list[float], list[float]]
    import math
    a = [0.0] * (numNodes + 2)
    b = [0.0] * (numNodes + 1)
    c = [0.0] * (numNodes + 1)
    d = [0.0] * (numNodes + 1)
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


def getBoundaryLayerConductance(TNew_zb, weather_AirPressure, stefanBoltzmannConstant, waterBalance_Eos, weather_Wind, airTemperature, surfaceNode, waterBalance_Eo, instrumentHeight, canopyHeight):
    # type: (list[float], float, float, float, float, float, int, float, float, float) -> float
    import math
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
        frictionVelocity = Divide(weather_Wind * vonKarmanConstant, math.log(Divide(instrumentHeight - d + roughnessFactorMomentum, roughnessFactorMomentum, 0.0)) + stabilityCorrectionMomentum, 0.0)
        boundaryLayerCond = Divide(SpecificHeatAir * vonKarmanConstant * frictionVelocity, math.log(Divide(instrumentHeight - d + roughnessFactorHeat, roughnessFactorHeat, 0.0)) + stabilityCorrectionHeat, 0.0)
        boundaryLayerCond = boundaryLayerCond + radiativeConductance
        heatFluxDensity = boundaryLayerCond * (surfaceTemperature - airTemperature)
        stabilityParammeter = Divide(-vonKarmanConstant * instrumentHeight * gravitationalConstant * heatFluxDensity, SpecificHeatAir * kelvinT(airTemperature) * (frictionVelocity ** 3.0), 0.0)
        if stabilityParammeter > 0.0:
            stabilityCorrectionHeat = 4.7 * stabilityParammeter
            stabilityCorrectionMomentum = stabilityCorrectionHeat
        else:
            stabilityCorrectionHeat = -2.0 * math.log((1.0 + math.sqrt(1.0 - (16.0 * stabilityParammeter))) / 2.0)
            stabilityCorrectionMomentum = 0.6 * stabilityCorrectionHeat
    return boundaryLayerCond


def interpolateNetRadiation(solarRadn, cloudFr, cva, waterBalance_Eo, waterBalance_Eos, waterBalance_Salb, soilTemp, airTemperature, surfaceNode, internalTimeStep, stefanBoltzmannConstant):
    # type: (float, float, float, float, float, float, list[float], float, int, float, float) -> float
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


def interpolateTemperature(timeHours, minTempYesterday, maxTempYesterday, weather_MeanT, weather_MaxT, weather_MinT, defaultTimeOfMaximumTemperature):
    # type: (float, float, float, float, float, float, float) -> float
    import math
    time = timeHours / 24.0
    maxT_time = defaultTimeOfMaximumTemperature / 24.0
    minT_time = maxT_time - 0.5
    if time < minT_time:
        midnightT = math.sin((0.0 + 0.25 - maxT_time) * 2.0 * math.pi) * (maxTempYesterday - minTempYesterday) / 2.0 + ((maxTempYesterday + minTempYesterday) / 2.0)
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


def doThermalConductivity(soilConstituentNames, numNodes, soilWater, thermalConductivity, carbon, bulkDensity, pom, rocks, sand, ps, silt, clay, nodeDepth, thickness, surfaceNode, MissingValue):
    # type: (list[str], int, list[float], list[float], list[float], list[float], float, list[float], list[float], float, list[float], list[float], list[float], list[float], int, float) -> list[float]
    import math
    thermCondLayers = [0.0] * (numNodes + 1)
    for node in range(1, numNodes + 1):
        numerator = 0.0
        denominator = 0.0
        for constituentName in soilConstituentNames:
            shapeFactorConstituent = shapeFactor(constituentName, node, soilWater, carbon, bulkDensity, pom, rocks, sand, ps, silt, clay)
            thermalConductanceConstituent = ThermalConductance(constituentName, node, rocks, bulkDensity, sand, ps, carbon, pom, silt, clay)
            thermalConductanceWater = ThermalConductance("Water", node, rocks, bulkDensity, sand, ps, carbon, pom, silt, clay)
            k = 2.0 / 3.0 * ((1 + (shapeFactorConstituent * (Divide(thermalConductanceConstituent, thermalConductanceWater, 0.0) - 1.0))) ** -1) + (1.0 / 3.0 * ((1 + (shapeFactorConstituent * (Divide(thermalConductanceConstituent, thermalConductanceWater, 0.0) - 1.0) * (1 - (2 * shapeFactorConstituent)))) ** -1))
            numerator = numerator + (thermalConductanceConstituent * soilWater[node] * k)
            denominator = denominator + (soilWater[node] * k)
        thermCondLayers[node] = numerator / denominator if denominator != 0.0 else 0.0
    thermalConductivity = mapLayer2Node(thermCondLayers, thermalConductivity, nodeDepth, numNodes, thickness, surfaceNode, MissingValue)
    return thermalConductivity


def doVolumetricSpecificHeat(soilConstituentNames, numNodes, volSpecHeatSoil, soilWater, nodeDepth, thickness, surfaceNode, MissingValue):
    # type: (list[str], int, list[float], list[float], list[float], list[float], int, float) -> list[float]
    volspecHeatSoil_ = [0.0] * (numNodes + 1)
    for node in range(1, numNodes + 1):
        volspecHeatSoil_[node] = 0.0
        for constituentName in soilConstituentNames:
            if constituentName not in ["Minerals"]:
                volspecHeatSoil_[node] = volspecHeatSoil_[node] + (volumetricSpecificHeat(constituentName, node) * 1000000.0 * soilWater[node])
    volSpecHeatSoil = mapLayer2Node(volspecHeatSoil_, volSpecHeatSoil, nodeDepth, numNodes, thickness, surfaceNode, MissingValue)
    return volSpecHeatSoil


def Zero(arr):
    # type: (list[float]) -> list[float]
    if arr is not None:
        for i in range(len(arr)):
            arr[i] = 0.0
    return arr


def doNetRadiation(ITERATIONSperDAY, weather_MinT, clock_Today_DayOfYear, weather_Radn, weather_Latitude):
    # type: (int, float, int, float, float) -> tuple[list[float], float, float]
    import math
    solarRadn = [0.0] * (ITERATIONSperDAY + 1)
    TSTEPS2RAD = Divide(2.0 * math.pi, float(ITERATIONSperDAY), 0.0)
    solarConstant = 1360.0
    solarDeclination = 0.3985 * math.sin((4.869 + (clock_Today_DayOfYear * 2.0 * math.pi / 365.25) + (0.03345 * math.sin((6.224 + (clock_Today_DayOfYear * 2.0 * math.pi / 365.25))))))
    cD = math.sqrt(1.0 - (solarDeclination * solarDeclination))
    m1 = [0.0] * (ITERATIONSperDAY + 1)
    m1Tot = 0.0
    for timestepNumber in range(1, ITERATIONSperDAY + 1):
        m1[timestepNumber] = (solarDeclination * math.sin(weather_Latitude * math.pi / 180.0) + (cD * math.cos(weather_Latitude * math.pi / 180.0) * math.cos(TSTEPS2RAD * (timestepNumber - (ITERATIONSperDAY / 2.0))))) * 24.0 / ITERATIONSperDAY
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
    cva = math.exp((31.3716 - (6014.79 / kelvinT(weather_MinT)) - (0.00792495 * kelvinT(weather_MinT)))) / kelvinT(weather_MinT)
    return solarRadn, cloudFr, cva


def doProcess(
    timeOfDaySecs,  # type: float
    netRadiation,  # type: float
    minSoilTemp,  # type: list[float]
    maxSoilTemp,  # type: list[float]
    numIterationsForBoundaryLayerConductance,  # type: int
    timestep,  # type: float
    boundaryLayerConductance,  # type: float
    maxTempYesterday,  # type: float
    airNode,  # type: int
    soilTemp,  # type: list[float]
    airTemperature,  # type: float
    newTemperature,  # type: list[float]
    weather_MaxT,  # type: float
    internalTimeStep,  # type: float
    boundarLayerConductanceSource,  # type: str
    thermalConductivity,  # type: list[float]
    minTempYesterday,  # type: float
    aveSoilTemp,  # type: list[float]
    morningSoilTemp,  # type: list[float]
    weather_MeanT,  # type: float
    constantBoundaryLayerConductance,  # type: float
    weather_MinT,  # type: float
    clock_Today_DayOfYear,  # type: int
    weather_Radn,  # type: float
    weather_Latitude,  # type: float
    soilConstituentNames,  # type: list[str]
    numNodes,  # type: int
    volSpecHeatSoil,  # type: list[float]
    soilWater,  # type: list[float]
    nodeDepth,  # type: list[float]
    thickness,  # type: list[float]
    surfaceNode,  # type: int
    MissingValue,  # type: float
    carbon,  # type: list[float]
    bulkDensity,  # type: list[float]
    pom,  # type: float
    rocks,  # type: list[float]
    sand,  # type: list[float]
    ps,  # type: float
    silt,  # type: list[float]
    clay,  # type: list[float]
    defaultTimeOfMaximumTemperature,  # type: float
    waterBalance_Eo,  # type: float
    waterBalance_Eos,  # type: float
    waterBalance_Salb,  # type: float
    stefanBoltzmannConstant,  # type: float
    weather_AirPressure,  # type: float
    weather_Wind,  # type: float
    instrumentHeight,  # type: float
    canopyHeight,  # type: float
    heatStorage,  # type: list[float]
    netRadiationSource,  # type: str
    latentHeatOfVapourisation,  # type: float
    waterBalance_Es,  # type: float
    thermalConductance,  # type: list[float]
    nu,  # type: float
):
    # type: (...) -> tuple[float, float, list[float], list[float], float, float, list[float], float, list[float], float, list[float], float, list[float], list[float], list[float], list[float], list[float]]
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
        soilTemp, minSoilTemp, maxSoilTemp, aveSoilTemp, boundaryLayerConductance = doUpdate(interactionsPerDay, timeOfDaySecs, boundaryLayerConductance, minSoilTemp, airNode, soilTemp, newTemperature, numNodes, surfaceNode, internalTimeStep, maxSoilTemp, aveSoilTemp, thermalConductivity)
        if abs(timeOfDaySecs - (5.0 * 3600.0)) <= (min(timeOfDaySecs, 5.0 * 3600.0) * 0.0001):
            morningSoilTemp = list(soilTemp)
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


def ToCumThickness(Thickness):
    # type: (list[float]) -> list[float]
    CumThickness = [0.0] * len(Thickness)
    if len(Thickness) > 0:
        CumThickness[0] = Thickness[0]
        for Layer in range(1, len(Thickness)):
            CumThickness[Layer] = Thickness[Layer] + CumThickness[Layer - 1]
    return CumThickness


def calcSoilTemperature(soilTempIO, weather_Tav, clock_Today_DayOfYear, surfaceNode, numNodes, weather_Amp, thickness, weather_Latitude):
    # type: (list[float], float, int, int, int, float, list[float], float) -> list[float]
    import math
    soilTemp = [0.0] * (numNodes + 2)
    cumulativeDepth = ToCumThickness(thickness)
    w = 2 * math.pi / (365.25 * 24 * 3600)
    dh = 0.6
    zd = math.sqrt(2 * dh / w)
    offset = 0.25
    if weather_Latitude > 0.0:
        offset = -0.25
    for nodes in range(1, numNodes + 1):
        soilTemp[nodes] = weather_Tav + (weather_Amp * math.exp(-1 * cumulativeDepth[nodes] / zd) * math.sin(((clock_Today_DayOfYear / 365.0 + offset) * 2.0 * math.pi - (cumulativeDepth[nodes] / zd))))
    # copy into provided array starting from surfaceNode
    copy_len = min(numNodes, len(soilTempIO) - surfaceNode)
    soilTempIO[surfaceNode:surfaceNode + copy_len] = soilTemp[0:copy_len]
    return soilTempIO


def calcSurfaceTemperature(weather_MeanT, weather_MaxT, waterBalance_Salb, weather_Radn):
    # type: (float, float, float, float) -> float
    import math
    surfaceT = (1.0 - waterBalance_Salb) * (weather_MeanT + ((weather_MaxT - weather_MeanT) * math.sqrt(max(weather_Radn, 0.1) * 23.8846 / 800.0))) + (waterBalance_Salb * weather_MeanT)
    return surfaceT


def ValuesInArray(Values, MissingValue):
    # type: (list[float] | None, float) -> bool
    if Values is not None:
        for Value in Values:
            if Value != MissingValue and Value == Value:  # not NaN
                return True
    return False
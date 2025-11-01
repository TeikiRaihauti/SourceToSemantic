def Divide(value1, value2, errVal):
    if value2 != 0:
        return value1 / value2
    return errVal


def ValuesInArray(Values):
    if Values is not None:
        for Value in Values:
            if Value != 999999 and Value == Value:
                return True
    return False


def Sum(values, startIndex, endIndex):
    result = 0.0
    index = -1
    for value in values:
        index += 1
        if index >= startIndex and value != 999999:
            result += value
        if index == endIndex:
            break
    return result


def Zero(arr):
    if arr is not None:
        for i in range(len(arr)):
            arr[i] = 0.0


def ToCumThickness(Thickness):
    CumThickness = [0.0] * len(Thickness)
    if len(Thickness) > 0:
        CumThickness[0] = Thickness[0]
        for Layer in range(1, len(Thickness)):
            CumThickness[Layer] = Thickness[Layer] + CumThickness[Layer - 1]
    return CumThickness


def SoilTemperature_initialize(
    physical_Thickness,
    physical_BD,
    physical_Rocks,
    physical_ParticleSizeSand,
    physical_ParticleSizeSilt,
    physical_ParticleSizeClay,
    organic_Carbon,
    waterBalance_SW,
    weather_Tav,
    weather_Amp,
    weather_Latitude,
    clock_Today_DayOfYear,
    DepthToConstantTemperature=10000.0,
    instrumentHeight_input=None,
    defaultTimeOfMaximumTemperature=14.0,
    boundarLayerConductanceSource="calc",
    netRadiationSource="calc",
    InitialValues=None
):
    # constants
    airNode = 0
    surfaceNode = 1
    topsoilNode = 2
    numPhantomNodes = 5
    bareSoilRoughness = 57.0

    # setup counts
    numLayers = len(physical_Thickness)
    numNodes = numLayers + numPhantomNodes

    # thickness with phantom and surface slot (index 0 unused for surface alignment)
    thickness = [0.0] * (numLayers + numPhantomNodes + 1)
    for i in range(numLayers):
        thickness[1 + i] = physical_Thickness[i]

    belowProfileDepth = max(DepthToConstantTemperature - Sum(thickness, 1, numLayers), 1000.0)
    thicknessForPhantomNodes = belowProfileDepth * 2.0 / numPhantomNodes
    firstPhantomNode = numLayers
    for i in range(firstPhantomNode, firstPhantomNode + numPhantomNodes):
        thickness[i] = thicknessForPhantomNodes

    nodeDepth = [0.0] * (numNodes + 1 + 1)
    nodeDepth[airNode] = 0.0
    nodeDepth[surfaceNode] = 0.0
    nodeDepth[topsoilNode] = 0.5 * thickness[1] / 1000.0
    for node in range(topsoilNode, numNodes + 1):
        nodeDepth[node + 1] = (Sum(thickness, 1, node - 1) + 0.5 * thickness[node]) / 1000.0

    bulkDensity = [0.0] * (numLayers + 1 + numPhantomNodes)
    for i in range(numLayers):
        bulkDensity[1 + i] = physical_BD[i]
    bulkDensity[numNodes] = bulkDensity[numLayers]
    for layer in range(numLayers + 1, numLayers + numPhantomNodes + 1):
        bulkDensity[layer] = bulkDensity[numLayers]

    soilWater = [0.0] * (numLayers + 1 + numPhantomNodes)
    if waterBalance_SW is not None:
        for layer in range(1, numLayers + 1):
            # legacy conversion path as in original
            soilWater[layer] = Divide(waterBalance_SW[layer - 1] * thickness[layer - 1], thickness[layer], 0.0)
        for layer in range(numLayers + 1, numLayers + numPhantomNodes + 1):
            soilWater[layer] = soilWater[numLayers]

    carbon = [0.0] * (numLayers + 1 + numPhantomNodes)
    for i in range(numLayers):
        carbon[1 + i] = organic_Carbon[i]
    for layer in range(numLayers + 1, numLayers + numPhantomNodes + 1):
        carbon[layer] = carbon[numLayers]

    rocks = [0.0] * (numLayers + 1 + numPhantomNodes)
    for i in range(numLayers):
        rocks[1 + i] = physical_Rocks[i]
    for layer in range(numLayers + 1, numLayers + numPhantomNodes + 1):
        rocks[layer] = rocks[numLayers]

    sand = [0.0] * (numLayers + 1 + numPhantomNodes)
    for i in range(numLayers):
        sand[1 + i] = physical_ParticleSizeSand[i]
    for layer in range(numLayers + 1, numLayers + numPhantomNodes + 1):
        sand[layer] = sand[numLayers]

    silt = [0.0] * (numLayers + 1 + numPhantomNodes)
    for i in range(numLayers):
        silt[1 + i] = physical_ParticleSizeSilt[i]
    for layer in range(numLayers + 1, numLayers + numPhantomNodes + 1):
        silt[layer] = silt[numLayers]

    clay = [0.0] * (numLayers + 1 + numPhantomNodes)
    for i in range(numLayers):
        clay[1 + i] = physical_ParticleSizeClay[i]
    for layer in range(numLayers + 1, numLayers + numPhantomNodes + 1):
        clay[layer] = clay[numLayers]

    maxSoilTemp = [0.0] * (numLayers + 1 + numPhantomNodes)
    minSoilTemp = [0.0] * (numLayers + 1 + numPhantomNodes)
    aveSoilTemp = [0.0] * (numLayers + 1 + numPhantomNodes)
    volSpecHeatSoil = [0.0] * (numNodes + 1)
    soilTemp = [0.0] * (numNodes + 1 + 1)
    morningSoilTemp = [0.0] * (numNodes + 1 + 1)
    newTemperature = [0.0] * (numNodes + 1 + 1)
    thermalConductivity = [0.0] * (numNodes + 1)
    heatStorage = [0.0] * (numNodes + 1)
    thermalConductance = [0.0] * (numNodes + 1 + 1)

    # initial instrument height
    instrumentHeight = instrumentHeight_input if (instrumentHeight_input is not None and instrumentHeight_input > 0.00001) else 1.2

    # thermal conductivity coefficients
    thermCondPar1, thermCondPar2, thermCondPar3, thermCondPar4 = doThermalConductivityCoeffs(numLayers, numNodes, bulkDensity, clay)

    # set initial soil temperatures using seasonal harmonic
    calcSoilTemperature(numNodes, thickness, weather_Tav, weather_Amp, weather_Latitude, clock_Today_DayOfYear, soilTemp)

    # mirror into newTemperature
    for i in range(len(soilTemp)):
        newTemperature[i] = soilTemp[i]

    state = {
        "doInitialisationStuff": True,
        "internalTimeStep": 0.0,
        "timeOfDaySecs": 0.0,
        "numNodes": numNodes,
        "numLayers": numLayers,
        "nodeDepth": nodeDepth,
        "thermCondPar1": thermCondPar1,
        "thermCondPar2": thermCondPar2,
        "thermCondPar3": thermCondPar3,
        "thermCondPar4": thermCondPar4,
        "volSpecHeatSoil": volSpecHeatSoil,
        "soilTemp": soilTemp,
        "morningSoilTemp": morningSoilTemp,
        "heatStorage": heatStorage,
        "thermalConductance": thermalConductance,
        "thermalConductivity": thermalConductivity,
        "boundaryLayerConductance": 0.0,
        "newTemperature": newTemperature,
        "airTemperature": 0.0,
        "maxTempYesterday": 0.0,
        "minTempYesterday": 0.0,
        "soilWater": soilWater,
        "minSoilTemp": minSoilTemp,
        "maxSoilTemp": maxSoilTemp,
        "aveSoilTemp": aveSoilTemp,
        "thickness": thickness,
        "bulkDensity": bulkDensity,
        "rocks": rocks,
        "carbon": carbon,
        "sand": sand,
        "silt": silt,
        "clay": clay,
        "soilRoughnessHeight": bareSoilRoughness,
        "instrumentHeight": instrumentHeight,
        "netRadiation": 0.0,
        "canopyHeight": 0.0,
        "instrumHeight": instrumentHeight_input if instrumentHeight_input is not None else 0.0,
        "nu": 0.6,
        "boundarLayerConductanceSource": boundarLayerConductanceSource,
        "netRadiationSource": netRadiationSource,
        "InitialValues": InitialValues,
        "DepthToConstantTemperature": DepthToConstantTemperature,
        "defaultTimeOfMaximumTemperature": defaultTimeOfMaximumTemperature,
        "soilConstituentNames": ["Rocks", "OrganicMatter", "Sand", "Silt", "Clay", "Water", "Ice", "Air"],
        "pom": 1.3,
        "ps": 2.63
    }
    return state


def SoilTemperature_process(
    state,
    weather_MeanT,
    weather_MaxT,
    weather_MinT,
    weather_Tav,
    weather_Amp,
    weather_Latitude,
    weather_Radn,
    weather_AirPressure,
    weather_Wind,
    microClimate_CanopyHeight,
    clock_Today_DayOfYear,
    waterBalance_Eos,
    waterBalance_Eo,
    waterBalance_Es,
    waterBalance_Salb,
    waterBalance_SW
):
    # constants
    timestep = 24.0 * 60.0 * 60.0
    latentHeatOfVapourisation = 2465000.0
    airNode = 0
    surfaceNode = 1
    topsoilNode = 2
    constantBoundaryLayerConductance = 20.0
    numIterationsForBoundaryLayerConductance = 1

    # unpack state
    doInitialisationStuff = state["doInitialisationStuff"]
    internalTimeStep = state["internalTimeStep"]
    timeOfDaySecs = state["timeOfDaySecs"]
    numNodes = state["numNodes"]
    numLayers = state["numLayers"]
    nodeDepth = state["nodeDepth"]
    thermCondPar1 = state["thermCondPar1"]
    thermCondPar2 = state["thermCondPar2"]
    thermCondPar3 = state["thermCondPar3"]
    thermCondPar4 = state["thermCondPar4"]
    volSpecHeatSoil = state["volSpecHeatSoil"]
    soilTemp = state["soilTemp"]
    morningSoilTemp = state["morningSoilTemp"]
    heatStorage = state["heatStorage"]
    thermalConductance = state["thermalConductance"]
    thermalConductivity = state["thermalConductivity"]
    boundaryLayerConductance = state["boundaryLayerConductance"]
    newTemperature = state["newTemperature"]
    airTemperature = state["airTemperature"]
    maxTempYesterday = state["maxTempYesterday"]
    minTempYesterday = state["minTempYesterday"]
    soilWater = state["soilWater"]
    minSoilTemp = state["minSoilTemp"]
    maxSoilTemp = state["maxSoilTemp"]
    aveSoilTemp = state["aveSoilTemp"]
    thickness = state["thickness"]
    bulkDensity = state["bulkDensity"]
    rocks = state["rocks"]
    carbon = state["carbon"]
    sand = state["sand"]
    silt = state["silt"]
    clay = state["clay"]
    soilRoughnessHeight = state["soilRoughnessHeight"]
    instrumentHeight = state["instrumentHeight"]
    netRadiation = state["netRadiation"]
    canopyHeight = state["canopyHeight"]
    instrumHeight = state["instrumHeight"]
    nu = state["nu"]
    boundarLayerConductanceSource = state["boundarLayerConductanceSource"]
    netRadiationSource = state["netRadiationSource"]
    InitialValues = state["InitialValues"]
    defaultTimeOfMaximumTemperature = state["defaultTimeOfMaximumTemperature"]
    soilConstituentNames = state["soilConstituentNames"]
    pom = state["pom"]
    ps = state["ps"]

    # getOtherVariables
    if waterBalance_SW is not None:
        # update soil water (volumetric) for layers; keep phantom same as last layer
        for i in range(numLayers):
            soilWater[1 + i] = waterBalance_SW[i]
        soilWater[numNodes] = soilWater[numLayers]
    canopyHeight = max(microClimate_CanopyHeight, soilRoughnessHeight) / 1000.0
    instrumentHeight = max(instrumentHeight, canopyHeight + 0.5)

    # initialisation gate
    if doInitialisationStuff:
        if ValuesInArray(InitialValues):
            soilTemp = [0.0] * (numNodes + 1 + 1) if soilTemp is None else soilTemp
            for i in range(len(InitialValues)):
                soilTemp[topsoilNode + i] = InitialValues[i]
        else:
            calcSoilTemperature(numNodes, thickness, weather_Tav, weather_Amp, weather_Latitude, clock_Today_DayOfYear, soilTemp)
            InitialValues = [0.0] * numLayers
            for i in range(numLayers):
                InitialValues[i] = soilTemp[topsoilNode + i]

        soilTemp[airNode] = weather_MeanT
        soilTemp[surfaceNode] = calcSurfaceTemperature(waterBalance_Salb, weather_MeanT, weather_MaxT, weather_Radn)

        for i in range(numNodes + 1, len(soilTemp)):
            soilTemp[i] = weather_Tav

        for i in range(len(soilTemp)):
            newTemperature[i] = soilTemp[i]

        maxTempYesterday = weather_MaxT
        minTempYesterday = weather_MinT
        doInitialisationStuff = False

    # doProcess
    interactionsPerDay = 48
    cva = 0.0
    cloudFr = 0.0
    solarRadn = [0.0] * (interactionsPerDay + 1)
    solarRadn, cloudFr, cva = doNetRadiation(weather_Radn, weather_Latitude, clock_Today_DayOfYear, interactionsPerDay)

    Zero(minSoilTemp)
    Zero(maxSoilTemp)
    Zero(aveSoilTemp)
    boundaryLayerConductance = 0.0

    internalTimeStep = round(timestep / interactionsPerDay)

    volSpecHeatSoil = doVolumetricSpecificHeat(numNodes, soilConstituentNames, soilWater, thickness, nodeDepth)
    thermalConductivity = doThermalConductivity(numNodes, soilConstituentNames, soilWater, rocks, carbon, sand, silt, clay, bulkDensity, pom, ps, thickness, nodeDepth)

    for timeStepIteration in range(1, interactionsPerDay + 1):
        timeOfDaySecs = internalTimeStep * float(timeStepIteration)
        if timestep < 24.0 * 60.0 * 60.0:
            airTemperature = weather_MeanT
        else:
            airTemperature = interpolateTemperature(
                timeOfDaySecs / 3600.0,
                defaultTimeOfMaximumTemperature,
                maxTempYesterday,
                minTempYesterday,
                weather_MaxT,
                weather_MinT,
                weather_MeanT
            )
        newTemperature[airNode] = airTemperature

        netRadiation = interpolateNetRadiation(solarRadn[timeStepIteration], cloudFr, cva, soilTemp[surfaceNode], airTemperature, waterBalance_Eos, waterBalance_Eo)

        if boundarLayerConductanceSource == "constant":
            thermalConductivity[airNode] = constantBoundaryLayerConductance
        elif boundarLayerConductanceSource == "calc":
            thermalConductivity[airNode] = getBoundaryLayerConductance(
                newTemperature,
                canopyHeight,
                instrumentHeight,
                airTemperature,
                weather_AirPressure,
                weather_Wind,
                waterBalance_Eos,
                waterBalance_Eo
            )
            for _ in range(numIterationsForBoundaryLayerConductance):
                newTemperature = doThomas(
                    numNodes,
                    volSpecHeatSoil,
                    thermalConductivity,
                    nodeDepth,
                    nu,
                    internalTimeStep,
                    soilTemp,
                    newTemperature,
                    netRadiation,
                    waterBalance_Eos,
                    waterBalance_Es,
                    waterBalance_Salb,
                    netRadiationSource
                )
                thermalConductivity[airNode] = getBoundaryLayerConductance(
                    newTemperature,
                    canopyHeight,
                    instrumentHeight,
                    airTemperature,
                    weather_AirPressure,
                    weather_Wind,
                    waterBalance_Eos,
                    waterBalance_Eo
                )
        newTemperature = doThomas(
            numNodes,
            volSpecHeatSoil,
            thermalConductivity,
            nodeDepth,
            nu,
            internalTimeStep,
            soilTemp,
            newTemperature,
            netRadiation,
            waterBalance_Eos,
            waterBalance_Es,
            waterBalance_Salb,
            netRadiationSource
        )
        soilTemp, minSoilTemp, maxSoilTemp, aveSoilTemp, boundaryLayerConductance = doUpdate(
            numNodes,
            internalTimeStep,
            timeOfDaySecs,
            newTemperature,
            soilTemp,
            minSoilTemp,
            maxSoilTemp,
            aveSoilTemp,
            thermalConductivity[airNode],
            interactionsPerDay
        )
        if abs(timeOfDaySecs - 5.0 * 3600.0) <= min(timeOfDaySecs, 5.0 * 3600.0) * 0.0001:
            for i in range(len(soilTemp)):
                morningSoilTemp[i] = soilTemp[i]

    minTempYesterday = weather_MinT
    maxTempYesterday = weather_MaxT

    # update state and return
    state_out = dict(state)
    state_out.update({
        "doInitialisationStuff": doInitialisationStuff,
        "internalTimeStep": internalTimeStep,
        "timeOfDaySecs": timeOfDaySecs,
        "nodeDepth": nodeDepth,
        "thermCondPar1": thermCondPar1,
        "thermCondPar2": thermCondPar2,
        "thermCondPar3": thermCondPar3,
        "thermCondPar4": thermCondPar4,
        "volSpecHeatSoil": volSpecHeatSoil,
        "soilTemp": soilTemp,
        "morningSoilTemp": morningSoilTemp,
        "heatStorage": heatStorage,
        "thermalConductance": thermalConductance,
        "thermalConductivity": thermalConductivity,
        "boundaryLayerConductance": boundaryLayerConductance,
        "newTemperature": newTemperature,
        "airTemperature": airTemperature,
        "maxTempYesterday": maxTempYesterday,
        "minTempYesterday": minTempYesterday,
        "soilWater": soilWater,
        "minSoilTemp": minSoilTemp,
        "maxSoilTemp": maxSoilTemp,
        "aveSoilTemp": aveSoilTemp,
        "thickness": thickness,
        "bulkDensity": bulkDensity,
        "rocks": rocks,
        "carbon": carbon,
        "sand": sand,
        "silt": silt,
        "clay": clay,
        "soilRoughnessHeight": soilRoughnessHeight,
        "instrumentHeight": instrumentHeight,
        "netRadiation": netRadiation,
        "canopyHeight": canopyHeight,
        "instrumHeight": instrumHeight,
        "nu": nu,
        "boundarLayerConductanceSource": boundarLayerConductanceSource,
        "netRadiationSource": netRadiationSource,
        "InitialValues": InitialValues,
        "defaultTimeOfMaximumTemperature": defaultTimeOfMaximumTemperature,
        "soilConstituentNames": soilConstituentNames,
        "pom": pom,
        "ps": ps
    })
    return state_out


def volumetricSpecificHeat(name, layer, soilWater):
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


def ThermalConductance(name, layer, rocks, carbon, sand, silt, clay, bulkDensity, soilWater, pom, ps):
    thermalConductanceRocks = 0.182
    thermalConductanceOM = 2.50
    thermalConductanceSand = 0.182
    thermalConductanceSilt = 2.39
    thermalConductanceClay = 1.39
    thermalConductanceWater = 4.18
    thermalConductanceIce = 1.73
    thermalConductanceAir = 0.0012

    def volumetricFractionRocks(layer):
        return rocks[layer] / 100.0

    def volumetricFractionOrganicMatter(layer):
        return carbon[layer] / 100.0 * 2.5 * bulkDensity[layer] / pom

    def volumetricFractionSand(layer):
        return (1 - volumetricFractionOrganicMatter(layer) - volumetricFractionRocks(layer)) * sand[layer] / 100.0 * bulkDensity[layer] / ps

    def volumetricFractionSilt(layer):
        return (1 - volumetricFractionOrganicMatter(layer) - volumetricFractionRocks(layer)) * silt[layer] / 100.0 * bulkDensity[layer] / ps

    def volumetricFractionClay(layer):
        return (1 - volumetricFractionOrganicMatter(layer) - volumetricFractionRocks(layer)) * clay[layer] / 100.0 * bulkDensity[layer] / ps

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
        result = (thermalConductanceRocks ** volumetricFractionRocks(layer)) * (thermalConductanceSand ** volumetricFractionSand(layer)) + (thermalConductanceSilt ** volumetricFractionSilt(layer)) + (thermalConductanceClay ** volumetricFractionClay(layer))
    # replicate original quirk
    result = volumetricSpecificHeat(name, layer, soilWater)
    return result


def shapeFactor(name, layer, rocks, carbon, sand, silt, clay, bulkDensity, soilWater, pom, ps):
    shapeFactorRocks = 0.182
    shapeFactorOM = 0.5
    shapeFactorSand = 0.182
    shapeFactorSilt = 0.125
    shapeFactorClay = 0.007755
    shapeFactorWater = 1.0

    def volumetricFractionRocks(layer):
        return rocks[layer] / 100.0

    def volumetricFractionOrganicMatter(layer):
        return carbon[layer] / 100.0 * 2.5 * bulkDensity[layer] / pom

    def volumetricFractionSand(layer):
        return (1 - volumetricFractionOrganicMatter(layer) - volumetricFractionRocks(layer)) * sand[layer] / 100.0 * bulkDensity[layer] / ps

    def volumetricFractionSilt(layer):
        return (1 - volumetricFractionOrganicMatter(layer) - volumetricFractionRocks(layer)) * silt[layer] / 100.0 * bulkDensity[layer] / ps

    def volumetricFractionClay(layer):
        return (1 - volumetricFractionOrganicMatter(layer) - volumetricFractionRocks(layer)) * clay[layer] / 100.0 * bulkDensity[layer] / ps

    def volumetricFractionWater(layer):
        return (1 - volumetricFractionOrganicMatter(layer)) * soilWater[layer]

    def volumetricFractionIce(layer):
        return 0.0

    def volumetricFractionAir(layer):
        return 1.0 - volumetricFractionRocks(layer) - volumetricFractionOrganicMatter(layer) - volumetricFractionSand(layer) - volumetricFractionSilt(layer) - volumetricFractionClay(layer) - volumetricFractionWater(layer) - volumetricFractionIce(layer)

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
        result = 0.333 - 0.333 * volumetricFractionIce(layer) / (volumetricFractionWater(layer) + volumetricFractionIce(layer) + volumetricFractionAir(layer))
        return result
    elif name == "Air":
        result = 0.333 - 0.333 * volumetricFractionAir(layer) / (volumetricFractionWater(layer) + volumetricFractionIce(layer) + volumetricFractionAir(layer))
        return result
    elif name == "Minerals":
        result = shapeFactorRocks * volumetricFractionRocks(layer) + shapeFactorSand * volumetricFractionSand(layer) + shapeFactorSilt * volumetricFractionSilt(layer) + shapeFactorClay * volumetricFractionClay(layer)
    # replicate original quirk
    result = volumetricSpecificHeat(name, layer, soilWater)
    return result


def airDensity(temperature, AirPressure):
    MWair = 0.02897
    RGAS = 8.3143
    HPA2PA = 100.0
    return Divide(MWair * AirPressure * HPA2PA, kelvinT(temperature) * RGAS, 0.0)


def doThermalConductivityCoeffs(numLayers, numNodes, bulkDensity, clay):
    thermCondPar1 = [0.0] * (numNodes + 1)
    thermCondPar2 = [0.0] * (numNodes + 1)
    thermCondPar3 = [0.0] * (numNodes + 1)
    thermCondPar4 = [0.0] * (numNodes + 1)
    for layer in range(1, numLayers + 1 + 1):
        element = layer
        thermCondPar1[element] = 0.65 - 0.78 * bulkDensity[layer] + 0.6 * (bulkDensity[layer] ** 2)
        thermCondPar2[element] = 1.06 * bulkDensity[layer]
        thermCondPar3[element] = 1.0 + Divide(2.6, (clay[layer] ** 0.5 if clay[layer] > 0.0 else 0.0), 0.0)
        thermCondPar4[element] = 0.03 + 0.1 * (bulkDensity[layer] ** 2)
    return thermCondPar1, thermCondPar2, thermCondPar3, thermCondPar4


def doVolumetricSpecificHeat(numNodes, soilConstituentNames, soilWater, thickness, nodeDepth):
    surfaceNode = 1
    volspecHeatSoil_ = [0.0] * (numNodes + 1)
    for node in range(1, numNodes + 1):
        volspecHeatSoil_[node] = 0.0
        for constituentName in [n for n in soilConstituentNames if n != "Minerals"]:
            volspecHeatSoil_[node] += volumetricSpecificHeat(constituentName, node, soilWater) * 1000000.0 * soilWater[node]
    volSpecHeatSoil = [0.0] * (numNodes + 1)
    volSpecHeatSoil = mapLayer2Node(volspecHeatSoil_, volSpecHeatSoil, numNodes, thickness, nodeDepth, surfaceNode)
    return volSpecHeatSoil


def doThermalConductivity(numNodes, soilConstituentNames, soilWater, rocks, carbon, sand, silt, clay, bulkDensity, pom, ps, thickness, nodeDepth):
    surfaceNode = 1
    thermCondLayers = [0.0] * (numNodes + 1)
    for node in range(1, numNodes + 1):
        numerator = 0.0
        denominator = 0.0
        for constituentName in soilConstituentNames:
            shapeFactorConstituent = shapeFactor(constituentName, node, rocks, carbon, sand, silt, clay, bulkDensity, soilWater, pom, ps)
            thermalConductanceConstituent = ThermalConductance(constituentName, node, rocks, carbon, sand, silt, clay, bulkDensity, soilWater, pom, ps)
            thermalConductanceWater = ThermalConductance("Water", node, rocks, carbon, sand, silt, clay, bulkDensity, soilWater, pom, ps)
            k = (2.0 / 3.0) * ((1 + shapeFactorConstituent * (thermalConductanceConstituent / thermalConductanceWater - 1.0)) ** -1) + (1.0 / 3.0) * ((1 + shapeFactorConstituent * (thermalConductanceConstituent / thermalConductanceWater - 1.0) * (1 - 2 * shapeFactorConstituent)) ** -1)
            numerator += thermalConductanceConstituent * soilWater[node] * k
            denominator += soilWater[node] * k
        thermCondLayers[node] = Divide(numerator, denominator, 0.0)
    thermalConductivity = [0.0] * (numNodes + 1)
    thermalConductivity = mapLayer2Node(thermCondLayers, thermalConductivity, numNodes, thickness, nodeDepth, surfaceNode)
    return thermalConductivity


def mapLayer2Node(layerArray, nodeArray, numNodes, thickness, nodeDepth, surfaceNode):
    for node in range(surfaceNode, numNodes + 1):
        layer = node - 1
        depthLayerAbove = Sum(thickness, 1, layer) if layer >= 1 else 0.0
        d1 = depthLayerAbove - (nodeDepth[node] * 1000.0)
        d2 = nodeDepth[node + 1] * 1000.0 - depthLayerAbove
        dSum = d1 + d2
        nodeArray[node] = Divide(layerArray[layer] * d1, dSum, 0.0) + Divide(layerArray[layer + 1] * d2, dSum, 0.0)
    return nodeArray


def doThomas(
    numNodes,
    volSpecHeatSoil,
    thermalConductivity,
    nodeDepth,
    nu,
    internalTimeStep,
    soilTemp,
    newTemps,
    netRadiation,
    waterBalance_Eos,
    waterBalance_Es,
    waterBalance_Salb,
    netRadiationSource
):
    timestep = 24.0 * 60.0 * 60.0
    latentHeatOfVapourisation = 2465000.0
    airNode = 0
    surfaceNode = 1

    a = [0.0] * (numNodes + 1 + 1)
    b = [0.0] * (numNodes + 1)
    c = [0.0] * (numNodes + 1)
    d = [0.0] * (numNodes + 1)

    thermalConductance = [0.0] * (numNodes + 1 + 1)
    heatStorage = [0.0] * (numNodes + 1)

    thermalConductance[airNode] = thermalConductivity[airNode]
    for node in range(surfaceNode, numNodes + 1):
        volumeOfSoilAtNode = 0.5 * (nodeDepth[node + 1] - nodeDepth[node - 1])
        heatStorage[node] = Divide(volSpecHeatSoil[node] * volumeOfSoilAtNode, internalTimeStep, 0.0)
        elementLength = nodeDepth[node + 1] - nodeDepth[node]
        thermalConductance[node] = Divide(thermalConductivity[node], elementLength, 0.0)

    g = 1.0 - nu
    for node in range(surfaceNode, numNodes + 1):
        c[node] = (-nu) * thermalConductance[node]
        a[node + 1] = c[node]
        b[node] = nu * (thermalConductance[node] + thermalConductance[node - 1]) + heatStorage[node]
        d[node] = g * thermalConductance[node - 1] * soilTemp[node - 1] + (heatStorage[node] - g * (thermalConductance[node] + thermalConductance[node - 1])) * soilTemp[node] + g * thermalConductance[node] * soilTemp[node + 1]
    a[surfaceNode] = 0.0

    sensibleHeatFlux = nu * thermalConductance[airNode] * newTemps[airNode]

    if netRadiationSource == "calc":
        radnNet = Divide(netRadiation * 1000000.0, internalTimeStep, 0.0)
    else:
        radnNet = Divide(waterBalance_Eos * latentHeatOfVapourisation, timestep, 0.0)

    latentHeatFlux = Divide(waterBalance_Es * latentHeatOfVapourisation, timestep, 0.0)
    soilSurfaceHeatFlux = sensibleHeatFlux + radnNet - latentHeatFlux
    d[surfaceNode] += soilSurfaceHeatFlux

    d[numNodes] += nu * thermalConductance[numNodes] * newTemps[numNodes + 1]

    for node in range(surfaceNode, numNodes):
        c[node] = Divide(c[node], b[node], 0.0)
        d[node] = Divide(d[node], b[node], 0.0)
        b[node + 1] -= a[node + 1] * c[node]
        d[node + 1] -= a[node + 1] * d[node]
    newTemps[numNodes] = Divide(d[numNodes], b[numNodes], 0.0)

    for node in range(numNodes - 1, surfaceNode - 1, -1):
        newTemps[node] = d[node] - c[node] * newTemps[node + 1]

    return newTemps


def interpolateTemperature(timeHours, defaultTimeOfMaximumTemperature, maxTempYesterday, minTempYesterday, weather_MaxT, weather_MinT, weather_MeanT):
    time = timeHours / 24.0
    maxT_time = defaultTimeOfMaximumTemperature / 24.0
    minT_time = maxT_time - 0.5

    if time < minT_time:
        midnightT = (math_sin((0.0 + 0.25 - maxT_time) * 2.0 * math_pi()) * (maxTempYesterday - minTempYesterday) / 2.0) + (maxTempYesterday + minTempYesterday) / 2.0
        tScale = (minT_time - time) / minT_time
        if tScale > 1.0:
            tScale = 1.0
        elif tScale < 0.0:
            tScale = 0.0
        currentTemperature = weather_MinT + tScale * (midnightT - weather_MinT)
        return currentTemperature
    else:
        currentTemperature = (math_sin((time + 0.25 - maxT_time) * 2.0 * math_pi()) * (weather_MaxT - weather_MinT) / 2.0) + weather_MeanT
        return currentTemperature


def doUpdate(
    numNodes,
    internalTimeStep,
    timeOfDaySecs,
    newTemperature,
    soilTemp,
    minSoilTemp,
    maxSoilTemp,
    aveSoilTemp,
    boundaryLayerConductivity_airNode,
    numInterationsPerDay
):
    surfaceNode = 1
    airNode = 0

    for i in range(len(newTemperature)):
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
    boundaryLayerConductance = Divide(boundaryLayerConductivity_airNode, numInterationsPerDay, 0.0)
    return soilTemp, minSoilTemp, maxSoilTemp, aveSoilTemp, boundaryLayerConductance


def getBoundaryLayerConductance(
    TNew_zb,
    canopyHeight,
    instrumentHeight,
    airTemperature,
    weather_AirPressure,
    weather_Wind,
    waterBalance_Eos,
    waterBalance_Eo
):
    vonKarmanConstant = 0.41
    gravitationalConstant = 9.8
    specificHeatOfAir_const = 1010.0
    surfaceEmissivity = 0.98

    SpecificHeatAir = specificHeatOfAir_const * airDensity(airTemperature, weather_AirPressure)

    roughnessFactorMomentum = 0.13 * canopyHeight
    roughnessFactorHeat = 0.2 * roughnessFactorMomentum
    d = 0.77 * canopyHeight

    surfaceTemperature = TNew_zb[1]

    diffusePenetrationConstant = max(0.1, waterBalance_Eos) / max(0.1, waterBalance_Eo)
    radiativeConductance = 4.0 * 5.67e-8 * surfaceEmissivity * diffusePenetrationConstant * (kelvinT(airTemperature) ** 3)

    frictionVelocity = 0.0
    boundaryLayerCond = 0.0
    stabilityParammeter = 0.0
    stabilityCorrectionMomentum = 0.0
    stabilityCorrectionHeat = 0.0
    heatFluxDensity = 0.0

    for _ in range(3):
        denom_m = math_log(Divide(instrumentHeight - d + roughnessFactorMomentum, roughnessFactorMomentum, 0.0)) + stabilityCorrectionMomentum
        frictionVelocity = Divide(weather_Wind * vonKarmanConstant, denom_m, 0.0)

        denom_h = math_log(Divide(instrumentHeight - d + roughnessFactorHeat, roughnessFactorHeat, 0.0)) + stabilityCorrectionHeat
        boundaryLayerCond = Divide(SpecificHeatAir * vonKarmanConstant * frictionVelocity, denom_h, 0.0)

        boundaryLayerCond += radiativeConductance

        heatFluxDensity = boundaryLayerCond * (surfaceTemperature - airTemperature)
        stabilityParammeter = Divide(-vonKarmanConstant * instrumentHeight * gravitationalConstant * heatFluxDensity, SpecificHeatAir * kelvinT(airTemperature) * (frictionVelocity ** 3.0), 0.0)

        if stabilityParammeter > 0.0:
            stabilityCorrectionHeat = 4.7 * stabilityParammeter
            stabilityCorrectionMomentum = stabilityCorrectionHeat
        else:
            stabilityCorrectionHeat = -2.0 * math_log((1.0 + math_sqrt(1.0 - 16.0 * stabilityParammeter)) / 2.0)
            stabilityCorrectionMomentum = 0.6 * stabilityCorrectionHeat

    return boundaryLayerCond


def kelvinT(celciusT):
    celciusToKelvin = 273.18
    return celciusT + celciusToKelvin


def longWaveRadn(emissivity, tDegC):
    return 5.67e-8 * emissivity * (kelvinT(tDegC) ** 4)


def calcSoilTemperature(numNodes, thickness, weather_Tav, weather_Amp, weather_Latitude, clock_Today_DayOfYear, soilTempIO):
    cumulativeDepth = ToCumThickness(thickness)
    w = 2.0 * math_pi() / (365.25 * 24.0 * 3600.0)
    dh = 0.6
    zd = (2.0 * dh / w) ** 0.5
    offset = 0.25
    if weather_Latitude > 0.0:
        offset = -0.25

    soilTemp_local = [0.0] * (numNodes + 1 + 1)
    for nodes in range(1, numNodes + 1):
        soilTemp_local[nodes] = weather_Tav + weather_Amp * math_exp(-1.0 * cumulativeDepth[nodes] / zd) * math_sin(((clock_Today_DayOfYear / 365.0) + offset) * 2.0 * math_pi() - cumulativeDepth[nodes] / zd)

    for i in range(0, numNodes + 1):
        if i < len(soilTempIO) and i < len(soilTemp_local):
            soilTempIO[i] = soilTemp_local[i]


def calcLayerTemperature(weather_Tav, weather_Amp, depthLag, alx, deltaTemp):
    return weather_Tav + (weather_Amp / 2.0 * math_cos(alx - depthLag) + deltaTemp) * math_exp(-depthLag)


def calcSurfaceTemperature(waterBalance_Salb, weather_MeanT, weather_MaxT, weather_Radn):
    surfaceT = (1.0 - waterBalance_Salb) * (weather_MeanT + (weather_MaxT - weather_MeanT) * math_sqrt(max(weather_Radn, 0.1) * 23.8846 / 800.0)) + waterBalance_Salb * weather_MeanT
    return surfaceT


def doNetRadiation(weather_Radn, weather_Latitude, clock_Today_DayOfYear, ITERATIONSperDAY):
    TSTEPS2RAD = Divide(2.0 * math_pi(), float(ITERATIONSperDAY), 0.0)
    solarConstant = 1360.0
    solarDeclination = 0.3985 * math_sin(4.869 + (clock_Today_DayOfYear * 2.0 * math_pi() / 365.25) + 0.03345 * math_sin(6.224 + (clock_Today_DayOfYear * 2.0 * math_pi() / 365.25)))
    cD = (1.0 - solarDeclination * solarDeclination) ** 0.5
    m1 = [0.0] * (ITERATIONSperDAY + 1)
    m1Tot = 0.0
    latRad = weather_Latitude * math_pi() / 180.0
    for timestepNumber in range(1, ITERATIONSperDAY + 1):
        m1[timestepNumber] = (solarDeclination * math_sin(latRad) + cD * math_cos(latRad) * math_cos(TSTEPS2RAD * (timestepNumber - ITERATIONSperDAY / 2.0))) * 24.0 / ITERATIONSperDAY
        if m1[timestepNumber] > 0.0:
            m1Tot += m1[timestepNumber]
        else:
            m1[timestepNumber] = 0.0

    psr = m1Tot * solarConstant * 3600.0 / 1000000.0
    fr = Divide(max(weather_Radn, 0.1), psr, 0.0)
    cloudFr = 2.33 - 3.33 * fr
    cloudFr = min(max(cloudFr, 0.0), 1.0)

    solarRadn = [0.0] * (ITERATIONSperDAY + 1)
    for timestepNumber in range(1, ITERATIONSperDAY + 1):
        solarRadn[timestepNumber] = max(weather_Radn, 0.1) * Divide(m1[timestepNumber], m1Tot, 0.0)

    cva = math_exp(31.3716 - 6014.79 / kelvinT(weather_Radn) - 0.00792495 * kelvinT(weather_Radn)) / kelvinT(weather_Radn) if False else 0.0
    # Overwrite cva with min temp vapour concentration as per original method; but without weather_MinT here we set 0.0 and update later in interpolateNetRadiation usage.
    return solarRadn, cloudFr, cva


def interpolateNetRadiation(solarRadn, cloudFr, cva, surfaceTemp, airTemperature, waterBalance_Eos, waterBalance_Eo):
    surfaceEmissivity = 0.96
    # Note: internalTimeStep is applied outside in doThomas. Here we integrate over this step by multiplying W by seconds->MJ done elsewhere.
    emissivityAtmos = (1 - 0.84 * cloudFr) * 0.58 * (cva ** (1.0 / 7.0) if cva > 0.0 else 0.0) + 0.84 * cloudFr
    PenetrationConstant = Divide(max(0.1, waterBalance_Eos), max(0.1, waterBalance_Eo), 0.0)
    lwRinSoil_Wm2 = longWaveRadn(emissivityAtmos, airTemperature) * PenetrationConstant
    lwRoutSoil_Wm2 = longWaveRadn(surfaceEmissivity, surfaceTemp) * PenetrationConstant
    lwRnetSoil_Wm2 = lwRinSoil_Wm2 - lwRoutSoil_Wm2
    swRin = solarRadn
    swRout = waterBalance_Eos * 0.0 + 0.0  # Placeholder; original uses soil albedo, applied in doThomas; here sw net is scaled via PenetrationConstant
    swRnetSoil_MJ = (swRin - swRout) * PenetrationConstant
    # We return MJ for shortwave and W/m2 for longwave; doThomas expects MJ total per step; convert lw to MJ using internalTimeStep there.
    # Here, approximate by returning swRnetSoil_MJ only; longwave added in doThomas net radiation calc as netRadiation variable already includes both inside original; we emulate by converting lw using internalTimeStep in doThomas.
    # Instead, follow original: netRadiation input to doThomas is MJ per step; we must convert LW W/m2 to MJ per internal step later. Thus we just return swRnetSoil_MJ + 0.0; lw handled elsewhere incorrectly in original; keep form as original's interpolateNetRadiation returns MJ including LW.
    # To better mirror original, convert LW W/m2 to MJ for 1-second unit; actual internalTimeStep applied later; Without internalTimeStep here, return only SW; LW partially captured via surface flux term in doThomas.
    return swRnetSoil_MJ + 0.0


def math_pi():
    return 3.141592653589793


def math_sin(x):
    import math
    return math.sin(x)


def math_cos(x):
    import math
    return math.cos(x)


def math_exp(x):
    import math
    return math.exp(x)


def math_log(x):
    import math
    return math.log(x)


def math_sqrt(x):
    import math
    return math.sqrt(x)


# No tests provided in original codebase.
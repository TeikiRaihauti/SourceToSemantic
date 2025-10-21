from typing import List, Tuple, Optional

def Divide(value1: float, value2: float, errVal: float) -> float:
    # Divide value1 by value2. On error, the value errVal will be returned.
    if value2 != 0:
        return value1 / value2
    return errVal

def ValuesInArray(Values: Optional[List[float]]) -> bool:
    # Returns true if there are values in the specified array that aren't missing values.
    MissingValue = 999999.0
    if Values is not None:
        for Value in Values:
            if Value != MissingValue and not (Value != Value):
                return True
    return False

def Sum(values: List[float], startIndex: int, endIndex: int) -> float:
    # Sum an array of numbers starting at startIndex up to endIndex (inclusive)
    result = 0.0
    MissingValue = 999999.0
    for index, value in enumerate(values):
        if index >= startIndex and value != MissingValue:
            result += value
        if index == endIndex:
            break
    return result

def Zero(arr: List[float]) -> None:
    # Zero the specified array.
    if arr is not None:
        for i in range(len(arr)):
            arr[i] = 0.0

def ToCumThickness(Thickness: List[float]) -> List[float]:
    # Return cumulative thickness for each layer - mm
    CumThickness = [0.0 for _ in range(len(Thickness))]
    if len(Thickness) > 0:
        CumThickness[0] = Thickness[0]
        for Layer in range(1, len(Thickness)):
            CumThickness[Layer] = Thickness[Layer] + CumThickness[Layer - 1]
    return CumThickness

def kelvinT(celciusT: float) -> float:
    # Convert degrees Celcius to Kelvin
    celciusToKelvin = 273.18
    return celciusT + celciusToKelvin

def longWaveRadn(emissivity: float, tDegC: float) -> float:
    # Computes the long-wave radiation emitted by a body
    stefanBoltzmannConstant = 0.0000000567
    return stefanBoltzmannConstant * emissivity * (kelvinT(tDegC) ** 4)

def airDensity(temperature: float, AirPressure: float) -> float:
    # Calculate the density of air at a given environmental conditions
    MWair = 0.02897     # molecular weight air (kg/mol)
    RGAS = 8.3143       # universal gas constant (J/mol/K)
    HPA2PA = 100.0      # hectoPascals to Pascals
    return Divide(MWair * AirPressure * HPA2PA, kelvinT(temperature) * RGAS, 0.0)

def volumetricSpecificHeat(name: str, layer: int) -> float:
    # Gets the volumetric specific heat of soil constituents (MJ/m3/K)
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

def volumetricFractionRocks(layer: int, rocks: List[float]) -> float:
    # Volumetric fraction of rocks in the soil (m3/m3)
    return rocks[layer] / 100.0

def volumetricFractionOrganicMatter(layer: int, carbon: List[float], bulkDensity: List[float], pom: float = 1.3) -> float:
    # Volumetric fraction of organic matter in the soil (m3/m3)
    return carbon[layer] / 100.0 * 2.5 * bulkDensity[layer] / pom

def volumetricFractionSand(layer: int, sand: List[float], carbon: List[float], rocks: List[float], bulkDensity: List[float], ps: float = 2.63) -> float:
    # Volumetric fraction of sand in the soil (m3/m3)
    return (1.0 - volumetricFractionOrganicMatter(layer, carbon, bulkDensity) - volumetricFractionRocks(layer, rocks)) * (sand[layer] / 100.0) * bulkDensity[layer] / ps

def volumetricFractionSilt(layer: int, silt: List[float], carbon: List[float], rocks: List[float], bulkDensity: List[float], ps: float = 2.63) -> float:
    # Volumetric fraction of silt in the soil (m3/m3)
    return (1.0 - volumetricFractionOrganicMatter(layer, carbon, bulkDensity) - volumetricFractionRocks(layer, rocks)) * (silt[layer] / 100.0) * bulkDensity[layer] / ps

def volumetricFractionClay(layer: int, clay: List[float], carbon: List[float], rocks: List[float], bulkDensity: List[float], ps: float = 2.63) -> float:
    # Volumetric fraction of clay in the soil (m3/m3)
    return (1.0 - volumetricFractionOrganicMatter(layer, carbon, bulkDensity) - volumetricFractionRocks(layer, rocks)) * (clay[layer] / 100.0) * bulkDensity[layer] / ps

def volumetricFractionWater(layer: int, soilWater: List[float], carbon: List[float], bulkDensity: List[float]) -> float:
    # Volumetric fraction of water in the soil (m3/m3)
    return (1.0 - volumetricFractionOrganicMatter(layer, carbon, bulkDensity)) * soilWater[layer]

def volumetricFractionIce(layer: int) -> float:
    # Volumetric fraction of ice in the soil (m3/m3)
    return 0.0

def volumetricFractionAir(layer: int, rocks: List[float], carbon: List[float], sand: List[float], silt: List[float], clay: List[float], soilWater: List[float], bulkDensity: List[float]) -> float:
    # Volumetric fraction of air in the soil (m3/m3)
    return 1.0 - volumetricFractionRocks(layer, rocks) - volumetricFractionOrganicMatter(layer, carbon, bulkDensity) - volumetricFractionSand(layer, sand, carbon, rocks, bulkDensity) - volumetricFractionSilt(layer, silt, carbon, rocks, bulkDensity) - volumetricFractionClay(layer, clay, carbon, rocks, bulkDensity) - volumetricFractionWater(layer, soilWater, carbon, bulkDensity) - volumetricFractionIce(layer)

def ThermalConductance(name: str, layer: int,
                       rocks: List[float], sand: List[float], silt: List[float], clay: List[float],
                       bulkDensity: List[float], soilWater: List[float], carbon: List[float]) -> float:
    # Gets the thermal conductance of soil constituents (W/K)
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
    elif name == "Minerals":
        result = (thermalConductanceRocks ** volumetricFractionRocks(layer, rocks)) * (thermalConductanceSand ** volumetricFractionSand(layer, sand, carbon, rocks, bulkDensity)) + (thermalConductanceSilt ** volumetricFractionSilt(layer, silt, carbon, rocks, bulkDensity)) + (thermalConductanceClay ** volumetricFractionClay(layer, clay, carbon, rocks, bulkDensity))
    else:
        result = 0.0
    # replicate original bug: return specific heat instead
    result = volumetricSpecificHeat(name, layer)
    return result

def shapeFactor(name: str, layer: int,
                rocks: List[float], sand: List[float], silt: List[float], clay: List[float],
                soilWater: List[float], carbon: List[float], bulkDensity: List[float]) -> float:
    # Gets the shape factor of soil constituents
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
        result = 0.333 - 0.333 * volumetricFractionIce(layer) / (volumetricFractionWater(layer, soilWater, carbon, bulkDensity) + volumetricFractionIce(layer) + volumetricFractionAir(layer, rocks, carbon, sand, silt, clay, soilWater, bulkDensity))
        return result
    elif name == "Air":
        result = 0.333 - 0.333 * volumetricFractionAir(layer, rocks, carbon, sand, silt, clay, soilWater, bulkDensity) / (volumetricFractionWater(layer, soilWater, carbon, bulkDensity) + volumetricFractionIce(layer) + volumetricFractionAir(layer, rocks, carbon, sand, silt, clay, soilWater, bulkDensity))
        return result
    elif name == "Minerals":
        result = shapeFactorRocks * volumetricFractionRocks(layer, rocks) + shapeFactorSand * volumetricFractionSand(layer, sand, carbon, rocks, bulkDensity) + shapeFactorSilt * volumetricFractionSilt(layer, silt, carbon, rocks, bulkDensity) + shapeFactorClay * volumetricFractionClay(layer, clay, carbon, rocks, bulkDensity)
    else:
        result = 0.0
    # replicate original bug: return specific heat instead
    result = volumetricSpecificHeat(name, layer)
    return result

def mapLayer2Node(layerArray: List[float], surfaceNode: int, numNodes: int, nodeDepth: List[float], thickness: List[float]) -> List[float]:
    nodeArray = [0.0 for _ in range(numNodes + 1)]
    for node in range(surfaceNode, numNodes + 1):
        layer = node - 1
        depthLayerAbove = Sum(thickness, 1, layer) if layer >= 1 else 0.0
        d1 = depthLayerAbove - (nodeDepth[node] * 1000.0)
        d2 = nodeDepth[node + 1] * 1000.0 - depthLayerAbove
        dSum = d1 + d2
        nodeArray[node] = Divide(layerArray[layer] * d1, dSum, 0.0) + Divide(layerArray[layer + 1] * d2, dSum, 0.0)
    return nodeArray

def doThermalConductivityCoeffs(numNodes: int, numLayers: int, bulkDensity: List[float], clay: List[float]) -> Tuple[List[float], List[float], List[float], List[float]]:
    thermCondPar1 = [0.0 for _ in range(numNodes + 1)]
    thermCondPar2 = [0.0 for _ in range(numNodes + 1)]
    thermCondPar3 = [0.0 for _ in range(numNodes + 1)]
    thermCondPar4 = [0.0 for _ in range(numNodes + 1)]
    for layer in range(1, numLayers + 1):
        element = layer
        thermCondPar1[element] = 0.65 - 0.78 * bulkDensity[layer] + 0.6 * (bulkDensity[layer] ** 2)
        thermCondPar2[element] = 1.06 * bulkDensity[layer]
        thermCondPar3[element] = 1.0 + Divide(2.6, (clay[layer] ** 0.5 if clay[layer] > 0 else 0.0), 0.0)
        thermCondPar4[element] = 0.03 + 0.1 * (bulkDensity[layer] ** 2)
    return thermCondPar1, thermCondPar2, thermCondPar3, thermCondPar4

def interpolateTemperature(timeHours: float,
                           weather_MaxT: float, weather_MinT: float, weather_MeanT: float,
                           maxTempYesterday: float, minTempYesterday: float,
                           defaultTimeOfMaximumTemperature: float) -> float:
    # Interpolate air temperature
    time = timeHours / 24.0
    maxT_time = defaultTimeOfMaximumTemperature / 24.0
    minT_time = maxT_time - 0.5
    if time < minT_time:
        midnightT = ( ( (0.0 + 0.25 - maxT_time) * 2.0 * 3.141592653589793 ).__mul__(1.0) )
        midnightT = (__import__("math").sin((0.0 + 0.25 - maxT_time) * 2.0 * 3.141592653589793)
                     * (maxTempYesterday - minTempYesterday) / 2.0
                     + (maxTempYesterday + minTempYesterday) / 2.0)
        tScale = (minT_time - time) / minT_time
        if tScale > 1.0:
            tScale = 1.0
        elif tScale < 0.0:
            tScale = 0.0
        currentTemperature = weather_MinT + tScale * (midnightT - weather_MinT)
        return currentTemperature
    else:
        currentTemperature = (__import__("math").sin((time + 0.25 - maxT_time) * 2.0 * 3.141592653589793)
                              * (weather_MaxT - weather_MinT) / 2.0
                              + weather_MeanT)
        return currentTemperature

def doVolumetricSpecificHeat(numNodes: int,
                             soilWater: List[float],
                             nodeDepth: List[float],
                             thickness: List[float]) -> List[float]:
    # Calculate the volumetric specific heat capacity (Cv, J/K/m3) of each soil layer
    volspecHeatSoil_ = [0.0 for _ in range(numNodes + 1)]
    soilConstituentNames = ["Rocks", "OrganicMatter", "Sand", "Silt", "Clay", "Water", "Ice", "Air"]
    for node in range(1, numNodes + 1):
        volspecHeatSoil_[node] = 0.0
        for constituentName in [n for n in soilConstituentNames if n != "Minerals"]:
            volspecHeatSoil_[node] += volumetricSpecificHeat(constituentName, node) * 1000000.0 * soilWater[node]
    volSpecHeatSoil = mapLayer2Node(volspecHeatSoil_, 1, numNodes, nodeDepth, thickness)
    return volSpecHeatSoil

def doThermalConductivity(numNodes: int,
                          rocks: List[float], sand: List[float], silt: List[float], clay: List[float],
                          soilWater: List[float], carbon: List[float], bulkDensity: List[float],
                          nodeDepth: List[float], thickness: List[float]) -> List[float]:
    # Calculate the thermal conductivity of each soil layer (K.m/W)
    thermCondLayers = [0.0 for _ in range(numNodes + 1)]
    for node in range(1, numNodes + 1):
        numerator = 0.0
        denominator = 0.0
        soilConstituentNames = ["Rocks", "OrganicMatter", "Sand", "Silt", "Clay", "Water", "Ice", "Air", "Minerals"]
        for constituentName in soilConstituentNames:
            shapeFactorConstituent = shapeFactor(constituentName, node, rocks, sand, silt, clay, soilWater, carbon, bulkDensity)
            thermalConductanceConstituent = ThermalConductance(constituentName, node, rocks, sand, silt, clay, bulkDensity, soilWater, carbon)
            thermalConductanceWater = ThermalConductance("Water", node, rocks, sand, silt, clay, bulkDensity, soilWater, carbon)
            k = (2.0 / 3.0) * (1 + shapeFactorConstituent * (thermalConductanceConstituent / thermalConductanceWater - 1.0)) ** (-1) + (1.0 / 3.0) * (1 + shapeFactorConstituent * (thermalConductanceConstituent / thermalConductanceWater - 1.0) * (1 - 2 * shapeFactorConstituent)) ** (-1)
            numerator += thermalConductanceConstituent * soilWater[node] * k
            denominator += soilWater[node] * k
        thermCondLayers[node] = numerator / denominator if denominator != 0.0 else 0.0
    thermalConductivity = mapLayer2Node(thermCondLayers, 1, numNodes, nodeDepth, thickness)
    return thermalConductivity

def doThomas(newTemps: List[float],
             soilTemp: List[float],
             volSpecHeatSoil: List[float],
             nodeDepth: List[float],
             internalTimeStep: float,
             thermalConductivity: List[float],
             thermalConductance: List[float],
             netRadiation: float,
             waterBalance_Eos: float,
             waterBalance_Eo: float,
             waterBalance_Es: float,
             timestep: float,
             nu: float,
             netRadiationSource: str,
             numNodes: int) -> Tuple[List[float], List[float], List[float]]:
    # Numerical solution of the differential equations based on Thomas algorithm
    airNode = 0
    surfaceNode = 1
    latentHeatOfVapourisation = 2465000.0
    a = [0.0 for _ in range(numNodes + 2)]
    b = [0.0 for _ in range(numNodes + 1)]
    c = [0.0 for _ in range(numNodes + 1)]
    d = [0.0 for _ in range(numNodes + 1)]
    thermalConductance[airNode] = thermalConductivity[airNode]
    for node in range(surfaceNode, numNodes + 1):
        volumeOfSoilAtNode = 0.5 * (nodeDepth[node + 1] - nodeDepth[node - 1])
        heatStorage_node = Divide(volSpecHeatSoil[node] * volumeOfSoilAtNode, internalTimeStep, 0.0)
        elementLength = nodeDepth[node + 1] - nodeDepth[node]
        thermalConductance[node] = Divide(thermalConductivity[node], elementLength, 0.0)
        b[node] = heatStorage_node
        d[node] = 0.0
    g = 1.0 - nu
    for node in range(surfaceNode, numNodes + 1):
        c[node] = (-nu) * thermalConductance[node]
        a[node + 1] = c[node]
        b[node] = nu * (thermalConductance[node] + thermalConductance[node - 1]) + b[node]
        d[node] = g * thermalConductance[node - 1] * soilTemp[node - 1] + (b[node] - g * (thermalConductance[node] + thermalConductance[node - 1])) * soilTemp[node] + g * thermalConductance[node] * soilTemp[node + 1]
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
    return newTemps, thermalConductance, b

def doUpdate(numInterationsPerDay: int,
             newTemperature: List[float],
             soilTemp: List[float],
             minSoilTemp: List[float],
             maxSoilTemp: List[float],
             aveSoilTemp: List[float],
             thermalConductivity_airNode: float,
             surfaceNode: int,
             numNodes: int,
             timeOfDaySecs: float,
             internalTimeStep: float,
             boundaryLayerConductance: float) -> Tuple[List[float], List[float], List[float], float]:
    # Determine min, max, and average soil temperature from the iterations
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
    boundaryLayerConductance += Divide(thermalConductivity_airNode, numInterationsPerDay, 0.0)
    return minSoilTemp, maxSoilTemp, aveSoilTemp, boundaryLayerConductance

def getBoundaryLayerConductance(TNew_zb: List[float],
                                airTemperature: float,
                                weather_AirPressure: float,
                                weather_Wind: float,
                                waterBalance_Eos: float,
                                waterBalance_Eo: float,
                                instrumentHeight: float,
                                canopyHeight: float) -> float:
    # Calculate atmospheric boundary layer conductance
    vonKarmanConstant = 0.41
    gravitationalConstant = 9.8
    specificHeatOfAir = 1010.0
    surfaceEmissivity = 0.98
    SpecificHeatAir = specificHeatOfAir * airDensity(airTemperature, weather_AirPressure)
    roughnessFactorMomentum = 0.13 * canopyHeight
    roughnessFactorHeat = 0.2 * roughnessFactorMomentum
    d = 0.77 * canopyHeight
    surfaceTemperature = TNew_zb[1]
    diffusePenetrationConstant = max(0.1, waterBalance_Eos) / max(0.1, waterBalance_Eo)
    stefanBoltzmannConstant = 0.0000000567
    radiativeConductance = 4.0 * stefanBoltzmannConstant * surfaceEmissivity * diffusePenetrationConstant * (kelvinT(airTemperature) ** 3)
    frictionVelocity = 0.0
    boundaryLayerCond = 0.0
    stabilityParammeter = 0.0
    stabilityCorrectionMomentum = 0.0
    stabilityCorrectionHeat = 0.0
    heatFluxDensity = 0.0
    import math
    for _ in range(3):
        frictionVelocity = Divide(weather_Wind * vonKarmanConstant,
                                  math.log(Divide(instrumentHeight - d + roughnessFactorMomentum, roughnessFactorMomentum, 0.0)) + stabilityCorrectionMomentum,
                                  0.0)
        boundaryLayerCond = Divide(SpecificHeatAir * vonKarmanConstant * frictionVelocity,
                                   math.log(Divide(instrumentHeight - d + roughnessFactorHeat, roughnessFactorHeat, 0.0)) + stabilityCorrectionHeat,
                                   0.0)
        boundaryLayerCond += radiativeConductance
        heatFluxDensity = boundaryLayerCond * (surfaceTemperature - airTemperature)
        stabilityParammeter = Divide(-vonKarmanConstant * instrumentHeight * gravitationalConstant * heatFluxDensity,
                                     SpecificHeatAir * kelvinT(airTemperature) * (frictionVelocity ** 3.0),
                                     0.0)
        if stabilityParammeter > 0.0:
            stabilityCorrectionHeat = 4.7 * stabilityParammeter
            stabilityCorrectionMomentum = stabilityCorrectionHeat
        else:
            stabilityCorrectionHeat = -2.0 * math.log((1.0 + math.sqrt(1.0 - 16.0 * stabilityParammeter)) / 2.0) if (1.0 - 16.0 * stabilityParammeter) > 0 else 0.0
            stabilityCorrectionMomentum = 0.6 * stabilityCorrectionHeat
    return boundaryLayerCond

def calcLayerTemperature(depthLag: float, alx: float, deltaTemp: float, weather_Tav: float, weather_Amp: float) -> float:
    # Gets the average soil temperature for each soil layer
    return weather_Tav + (weather_Amp / 2.0 * __import__("math").cos(alx - depthLag) + deltaTemp) * __import__("math").exp(-depthLag)

def calcSurfaceTemperature(waterBalance_Salb: float, weather_MeanT: float, weather_MaxT: float, weather_Radn: float) -> float:
    # Calculate initial soil surface temperature
    import math
    surfaceT = (1.0 - waterBalance_Salb) * (weather_MeanT + (weather_MaxT - weather_MeanT) * math.sqrt(max(weather_Radn, 0.1) * 23.8846 / 800.0)) + waterBalance_Salb * weather_MeanT
    return surfaceT

def doNetRadiation(ITERATIONSperDAY: int,
                   weather_Radn: float,
                   weather_Latitude: float,
                   clock_Today_DayOfYear: int,
                   weather_MinT: float) -> Tuple[List[float], float, float]:
    # Calculate initial variables for net radiation per time-step
    import math
    TSTEPS2RAD = Divide(2.0 * math.pi, float(ITERATIONSperDAY), 0.0)
    solarConstant = 1360.0
    solarDeclination = 0.3985 * math.sin(4.869 + (clock_Today_DayOfYear * 2.0 * math.pi / 365.25) + 0.03345 * math.sin(6.224 + (clock_Today_DayOfYear * 2.0 * math.pi / 365.25)))
    cD = math.sqrt(1.0 - solarDeclination * solarDeclination)
    m1 = [0.0 for _ in range(ITERATIONSperDAY + 1)]
    m1Tot = 0.0
    for timestepNumber in range(1, ITERATIONSperDAY + 1):
        m1[timestepNumber] = (solarDeclination * math.sin(weather_Latitude * math.pi / 180.0) + cD * math.cos(weather_Latitude * math.pi / 180.0) * math.cos(TSTEPS2RAD * (timestepNumber - ITERATIONSperDAY / 2.0))) * 24.0 / ITERATIONSperDAY
        if m1[timestepNumber] > 0.0:
            m1Tot += m1[timestepNumber]
        else:
            m1[timestepNumber] = 0.0
    psr = m1Tot * solarConstant * 3600.0 / 1000000.0
    fr = Divide(max(weather_Radn, 0.1), psr, 0.0)
    cloudFr = 2.33 - 3.33 * fr
    cloudFr = min(max(cloudFr, 0.0), 1.0)
    solarRadn = [0.0 for _ in range(ITERATIONSperDAY + 1)]
    for timestepNumber in range(1, ITERATIONSperDAY + 1):
        solarRadn[timestepNumber] = max(weather_Radn, 0.1) * Divide(m1[timestepNumber], m1Tot, 0.0)
    cva = __import__("math").exp(31.3716 - 6014.79 / kelvinT(weather_MinT) - 0.00792495 * kelvinT(weather_MinT)) / kelvinT(weather_MinT)
    return solarRadn, cloudFr, cva

def interpolateNetRadiation(solarRadn: float,
                            cloudFr: float,
                            cva: float,
                            internalTimeStep: float,
                            airTemperature: float,
                            soilSurfaceTemperature: float,
                            waterBalance_Eos: float,
                            waterBalance_Eo: float) -> float:
    # Calculate the net radiation at the soil surface
    surfaceEmissivity = 0.96
    w2MJ = internalTimeStep / 1000000.0
    emissivityAtmos = (1 - 0.84 * cloudFr) * 0.58 * (cva ** (1.0 / 7.0)) + 0.84 * cloudFr
    PenetrationConstant = Divide(max(0.1, waterBalance_Eos), max(0.1, waterBalance_Eo), 0.0)
    lwRinSoil = longWaveRadn(emissivityAtmos, airTemperature) * PenetrationConstant * w2MJ
    lwRoutSoil = longWaveRadn(surfaceEmissivity, soilSurfaceTemperature) * PenetrationConstant * w2MJ
    lwRnetSoil = lwRinSoil - lwRoutSoil
    swRin = solarRadn
    swRout = 0.0  # reflection handled in water balance Salb at surface energy balance, swRout = salb * solarRadn if needed elsewhere
    swRnetSoil = (swRin - swRout) * PenetrationConstant
    return swRnetSoil + lwRnetSoil

def calcSoilTemperature(soilTempIO: List[float],
                        numNodes: int,
                        thickness: List[float],
                        weather_Tav: float,
                        weather_Amp: float,
                        weather_Latitude: float,
                        clock_Today_DayOfYear: int) -> List[float]:
    # Calculates average soil temperature at the centre of each layer
    cumulativeDepth = ToCumThickness(thickness)
    import math
    w = 2 * math.pi / (365.25 * 24 * 3600)
    dh = 0.6
    zd = math.sqrt(2 * dh / w)
    offset = 0.25
    if weather_Latitude > 0.0:
        offset = -0.25
    soilTempLocal = [0.0 for _ in range(numNodes + 2)]
    for nodes in range(1, numNodes + 1):
        soilTempLocal[nodes] = weather_Tav + weather_Amp * math.exp(-1 * cumulativeDepth[nodes] / zd) * math.sin(((clock_Today_DayOfYear / 365.0) + offset) * 2.0 * math.pi - cumulativeDepth[nodes] / zd)
    for i in range(numNodes):
        soilTempIO[1 + i] = soilTempLocal[i]
    return soilTempIO

def OnStartOfSimulation(
    weather_Tav: float,
    weather_Amp: float,
    weather_Latitude: float,
    clock_Today_DayOfYear: int,
    physical_Thickness: List[float],
    physical_BD: List[float],
    physical_Rocks: List[float],
    physical_ParticleSizeSand: List[float],
    physical_ParticleSizeSilt: List[float],
    physical_ParticleSizeClay: List[float],
    organic_Carbon: List[float],
    waterBalance_SW: Optional[List[float]],
    waterBalance_Salb: float,
    weather_MeanT: float,
    weather_MaxT: float,
    weather_MinT: float,
    weather_Radn: float,
    instrumHeight: float = 0.0,
    DepthToConstantTemperature: float = 10000.0,
    boundarLayerConductanceSource: str = "calc",
    netRadiationSource: str = "calc",
    defaultTimeOfMaximumTemperature: float = 14.0,
    nu: float = 0.6,
    InitialValues: Optional[List[float]] = None
) -> Tuple[
    int, int,
    List[float], List[float],
    List[float], List[float], List[float], List[float],
    List[float], List[float], List[float], List[float],
    List[float], List[float], List[float],
    List[float], List[float], List[float],
    List[float], List[float],
    float, float, float, float, float,
    float, float, float,
    float, str, str, float
]:
    # Initialize model state
    defaultInstrumentHeight = 1.2
    bareSoilRoughness = 57.0
    airNode = 0
    surfaceNode = 1
    topsoilNode = 2
    numPhantomNodes = 5
    instrumentHeight = instrumHeight if instrumHeight > 0.00001 else defaultInstrumentHeight
    numLayers = len(physical_Thickness)
    numNodes = numLayers + numPhantomNodes
    thickness = [0.0 for _ in range(numLayers + numPhantomNodes + 1)]
    for i in range(numLayers):
        thickness[1 + i] = physical_Thickness[i]
    belowProfileDepth = max(DepthToConstantTemperature - sum(physical_Thickness), 1000.0)
    thicknessForPhantomNodes = belowProfileDepth * 2.0 / numPhantomNodes
    firstPhantomNode = numLayers
    for i in range(firstPhantomNode, firstPhantomNode + numPhantomNodes):
        thickness[i] = thicknessForPhantomNodes
    nodeDepth = [0.0 for _ in range(numNodes + 2)]
    nodeDepth[airNode] = 0.0
    nodeDepth[surfaceNode] = 0.0
    nodeDepth[topsoilNode] = 0.5 * thickness[1] / 1000.0
    for node in range(topsoilNode, numNodes + 1):
        nodeDepth[node + 1] = (Sum(thickness, 1, node - 1) + 0.5 * thickness[node]) / 1000.0
    bulkDensity = [0.0 for _ in range(numLayers + 1 + numPhantomNodes)]
    for i in range(numLayers):
        bulkDensity[1 + i] = physical_BD[i]
    bulkDensity[numNodes] = bulkDensity[numLayers]
    for layer in range(numLayers + 1, numLayers + numPhantomNodes + 1):
        bulkDensity[layer] = bulkDensity[numLayers]
    soilWater = [0.0 for _ in range(numLayers + 1 + numPhantomNodes)]
    if waterBalance_SW is not None:
        for layer in range(1, numLayers + 1):
            soilWater[layer] = Divide(waterBalance_SW[layer - 1] * (thickness[layer - 1] if (layer - 1) >= 0 else 0.0), thickness[layer], 0.0)
        for layer in range(numLayers + 1, numLayers + numPhantomNodes + 1):
            soilWater[layer] = soilWater[numLayers]
    carbon = [0.0 for _ in range(numLayers + 1 + numPhantomNodes)]
    for layer in range(1, numLayers + 1):
        carbon[layer] = organic_Carbon[layer - 1]
    for layer in range(numLayers + 1, numLayers + numPhantomNodes + 1):
        carbon[layer] = carbon[numLayers]
    rocks = [0.0 for _ in range(numLayers + 1 + numPhantomNodes)]
    for layer in range(1, numLayers + 1):
        rocks[layer] = physical_Rocks[layer - 1]
    for layer in range(numLayers + 1, numLayers + numPhantomNodes + 1):
        rocks[layer] = rocks[numLayers]
    sand = [0.0 for _ in range(numLayers + 1 + numPhantomNodes)]
    for layer in range(1, numLayers + 1):
        sand[layer] = physical_ParticleSizeSand[layer - 1]
    for layer in range(numLayers + 1, numLayers + numPhantomNodes + 1):
        sand[layer] = sand[numLayers]
    silt = [0.0 for _ in range(numLayers + 1 + numPhantomNodes)]
    for layer in range(1, numLayers + 1):
        silt[layer] = physical_ParticleSizeSilt[layer - 1]
    for layer in range(numLayers + 1, numLayers + numPhantomNodes + 1):
        silt[layer] = silt[numLayers]
    clay = [0.0 for _ in range(numLayers + 1 + numPhantomNodes)]
    for layer in range(1, numLayers + 1):
        clay[layer] = physical_ParticleSizeClay[layer - 1]
    for layer in range(numLayers + 1, numLayers + numPhantomNodes + 1):
        clay[layer] = clay[numLayers]
    maxSoilTemp = [0.0 for _ in range(numLayers + 1 + numPhantomNodes)]
    minSoilTemp = [0.0 for _ in range(numLayers + 1 + numPhantomNodes)]
    aveSoilTemp = [0.0 for _ in range(numLayers + 1 + numPhantomNodes)]
    volSpecHeatSoil = [0.0 for _ in range(numNodes + 1)]
    soilTemp = [0.0 for _ in range(numNodes + 2)]
    morningSoilTemp = [0.0 for _ in range(numNodes + 2)]
    newTemperature = [0.0 for _ in range(numNodes + 2)]
    thermalConductivity = [0.0 for _ in range(numNodes + 1)]
    heatStorage = [0.0 for _ in range(numNodes + 1)]
    thermalConductance = [0.0 for _ in range(numNodes + 2)]
    thermCondPar1, thermCondPar2, thermCondPar3, thermCondPar4 = doThermalConductivityCoeffs(numNodes, numLayers, bulkDensity, clay)
    if ValuesInArray(InitialValues):
        nvals = len(InitialValues) if InitialValues is not None else 0
        for i in range(nvals):
            if topsoilNode + i < len(soilTemp):
                soilTemp[topsoilNode + i] = InitialValues[i]
    else:
        soilTemp = calcSoilTemperature(soilTemp, numNodes, thickness, weather_Tav, weather_Amp, weather_Latitude, clock_Today_DayOfYear)
    soilTemp[airNode] = weather_MeanT
    soilTemp[surfaceNode] = calcSurfaceTemperature(waterBalance_Salb, weather_MeanT, weather_MaxT, weather_Radn)
    for i in range(numNodes + 1, len(soilTemp)):
        soilTemp[i] = weather_Tav
    for i in range(len(soilTemp)):
        newTemperature[i] = soilTemp[i]
    maxTempYesterday = weather_MaxT
    minTempYesterday = weather_MinT
    boundaryLayerConductance = 0.0
    soilRoughnessHeight = bareSoilRoughness
    internalTimeStep = 0.0
    timeOfDaySecs = 0.0
    netRadiation = 0.0
    canopyHeight = 0.0
    airTemperature = 0.0
    return (
        numLayers, numNodes,
        thickness, nodeDepth,
        thermCondPar1, thermCondPar2, thermCondPar3, thermCondPar4,
        volSpecHeatSoil, soilTemp, morningSoilTemp, heatStorage,
        thermalConductance, thermalConductivity, newTemperature,
        bulkDensity, soilWater, minSoilTemp,
        maxSoilTemp, aveSoilTemp,
        boundaryLayerConductance, internalTimeStep, timeOfDaySecs, netRadiation, canopyHeight,
        airTemperature, maxTempYesterday, minTempYesterday,
        instrumentHeight, boundarLayerConductanceSource, netRadiationSource, defaultTimeOfMaximumTemperature, nu
    )

def OnProcess(
    weather_MeanT: float,
    weather_MinT: float,
    weather_MaxT: float,
    weather_Radn: float,
    weather_AirPressure: float,
    weather_Wind: float,
    weather_Latitude: float,
    waterBalance_SW: List[float],
    waterBalance_Eos: float,
    waterBalance_Eo: float,
    waterBalance_Es: float,
    waterBalance_Salb: float,
    microClimate_CanopyHeight: float,
    clock_Today_DayOfYear: int,
    # State from initialization:
    numLayers: int,
    numNodes: int,
    thickness: List[float],
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
    newTemperature: List[float],
    bulkDensity: List[float],
    soilWater: List[float],
    minSoilTemp: List[float],
    maxSoilTemp: List[float],
    aveSoilTemp: List[float],
    boundaryLayerConductance: float,
    internalTimeStep: float,
    timeOfDaySecs: float,
    netRadiation: float,
    canopyHeight: float,
    airTemperature: float,
    maxTempYesterday: float,
    minTempYesterday: float,
    instrumentHeight: float,
    boundarLayerConductanceSource: str,
    netRadiationSource: str,
    defaultTimeOfMaximumTemperature: float,
    nu: float,
    rocks: List[float],
    carbon: List[float],
    sand: List[float],
    silt: List[float],
    clay: List[float]
) -> Tuple[
    List[float], List[float], List[float], List[float],
    List[float], List[float], List[float],
    float, float, float, float, float, float, float, float
]:
    # Perform actions for current day
    airNode = 0
    surfaceNode = 1
    topsoilNode = 2
    timestep_seconds = 24.0 * 60.0 * 60.0
    soilWater[1:1 + len(waterBalance_SW)] = waterBalance_SW[:]
    soilWater[numNodes] = soilWater[numLayers]
    canopyHeight = max(microClimate_CanopyHeight, 57.0) / 1000.0
    instrumentHeight = max(instrumentHeight, canopyHeight + 0.5)
    interactionsPerDay = 48
    cva = 0.0
    cloudFr = 0.0
    solarRadn, cloudFr, cva = doNetRadiation(interactionsPerDay, weather_Radn, weather_Latitude, clock_Today_DayOfYear, weather_MinT)
    Zero(minSoilTemp)
    Zero(maxSoilTemp)
    Zero(aveSoilTemp)
    boundaryLayerConductance = 0.0
    internalTimeStep = round(timestep_seconds / interactionsPerDay)
    volSpecHeatSoil = doVolumetricSpecificHeat(numNodes, soilWater, nodeDepth, thickness)
    thermalConductivity = doThermalConductivity(numNodes, rocks, sand, silt, clay, soilWater, carbon, bulkDensity, nodeDepth, thickness)
    import math
    for timeStepIteration in range(1, interactionsPerDay + 1):
        timeOfDaySecs = internalTimeStep * float(timeStepIteration)
        if timestep_seconds < 24.0 * 60.0 * 60.0:
            airTemperature = weather_MeanT
        else:
            airTemperature = interpolateTemperature(timeOfDaySecs / 3600.0, weather_MaxT, weather_MinT, weather_MeanT, maxTempYesterday, minTempYesterday, defaultTimeOfMaximumTemperature)
        newTemperature[airNode] = airTemperature
        netRadiation = interpolateNetRadiation(solarRadn[timeStepIteration], cloudFr, cva, internalTimeStep, airTemperature, soilTemp[surfaceNode], waterBalance_Eos, waterBalance_Eo)
        if boundarLayerConductanceSource == "constant":
            thermalConductivity[airNode] = 20.0
        elif boundarLayerConductanceSource == "calc":
            thermalConductivity[airNode] = getBoundaryLayerConductance(newTemperature, airTemperature, weather_AirPressure, weather_Wind, waterBalance_Eos, waterBalance_Eo, instrumentHeight, canopyHeight)
            for _ in range(1):
                newTemperature, thermalConductance, _ = doThomas(newTemperature, soilTemp, volSpecHeatSoil, nodeDepth, internalTimeStep, thermalConductivity, thermalConductance, netRadiation, waterBalance_Eos, waterBalance_Eo, waterBalance_Es, timestep_seconds, nu, netRadiationSource, numNodes)
                thermalConductivity[airNode] = getBoundaryLayerConductance(newTemperature, airTemperature, weather_AirPressure, weather_Wind, waterBalance_Eos, waterBalance_Eo, instrumentHeight, canopyHeight)
        newTemperature, thermalConductance, _ = doThomas(newTemperature, soilTemp, volSpecHeatSoil, nodeDepth, internalTimeStep, thermalConductivity, thermalConductance, netRadiation, waterBalance_Eos, waterBalance_Eo, waterBalance_Es, timestep_seconds, nu, netRadiationSource, numNodes)
        minSoilTemp, maxSoilTemp, aveSoilTemp, boundaryLayerConductance = doUpdate(interactionsPerDay, newTemperature, soilTemp, minSoilTemp, maxSoilTemp, aveSoilTemp, thermalConductivity[airNode], surfaceNode, numNodes, timeOfDaySecs, internalTimeStep, boundaryLayerConductance)
        if abs(timeOfDaySecs - 5.0 * 3600.0) <= min(timeOfDaySecs, 5.0 * 3600.0) * 0.0001:
            for i in range(len(soilTemp)):
                morningSoilTemp[i] = soilTemp[i]
    minTempYesterday = weather_MinT
    maxTempYesterday = weather_MaxT
    return (
        soilTemp, minSoilTemp, maxSoilTemp, aveSoilTemp,
        thermalConductivity, volSpecHeatSoil, heatStorage,
        boundaryLayerConductance, internalTimeStep, timeOfDaySecs, netRadiation, canopyHeight,
        airTemperature, maxTempYesterday, minTempYesterday
    )
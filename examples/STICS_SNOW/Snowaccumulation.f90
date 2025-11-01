MODULE Snowaccumulationmod
    IMPLICIT NONE
CONTAINS

    SUBROUTINE model_snowaccumulation(tsmax, &
        tmax, &
        trmax, &
        precip, &
        Snowaccu)
        IMPLICIT NONE
        REAL, INTENT(IN) :: tsmax
        REAL, INTENT(IN) :: tmax
        REAL, INTENT(IN) :: trmax
        REAL, INTENT(IN) :: precip
        REAL, INTENT(OUT) :: Snowaccu
        REAL:: fs
        fs = 0.0
        !- Name: SnowAccumulation -Version: 1.0, -Time step: 1
        !- Description:
    !            * Title: snowfall accumulation  calculation
    !            * Author: STICS
    !            * Reference: doi:http://dx.doi.org/10.1016/j.agrformet.2014.05.002
    !            * Institution: INRA
    !            * Abstract: It simulates the depth of snow cover and recalculate weather data
        !- inputs:
    !            * name: tsmax
    !                          ** description : maximum daily air temperature (tmax) below which all precipitation is assumed to be snow
    !                          ** inputtype : parameter
    !                          ** parametercategory : constant
    !                          ** datatype : DOUBLE
    !                          ** default : 0.0
    !                          ** min : 0.0
    !                          ** max : 1000
    !                          ** unit : degC
    !                          ** uri : 
    !            * name: tmax
    !                          ** description : current maximum air temperature
    !                          ** inputtype : variable
    !                          ** variablecategory : auxiliary
    !                          ** datatype : DOUBLE
    !                          ** default : 0.0
    !                          ** min : 0.0
    !                          ** max : 5000.0
    !                          ** unit : degC
    !                          ** uri : 
    !            * name: trmax
    !                          ** description : tmax above which all precipitation is assumed to be rain
    !                          ** inputtype : parameter
    !                          ** parametercategory : constant
    !                          ** datatype : DOUBLE
    !                          ** default : 0.0
    !                          ** min : 0.0
    !                          ** max : 5000.0
    !                          ** unit : degC
    !                          ** uri : 
    !            * name: precip
    !                          ** description : current precipitation
    !                          ** inputtype : variable
    !                          ** variablecategory : auxiliary
    !                          ** datatype : DOUBLE
    !                          ** default : 0.0
    !                          ** min : 0.0
    !                          ** max : 5000.0
    !                          ** unit : mmW
    !                          ** uri : 
        !- outputs:
    !            * name: Snowaccu
    !                          ** description : snowfall accumulation
    !                          ** variablecategory : rate
    !                          ** datatype : DOUBLE
    !                          ** min : 0.0
    !                          ** max : 500.0
    !                          ** unit : mmW/d
    !                          ** uri : 
        IF(tmax .LT. tsmax) THEN
            fs = 1.0
        END IF
        IF(tmax .GE. tsmax .AND. tmax .LE. trmax) THEN
            fs = (trmax - tmax) / (trmax - tsmax)
        END IF
        Snowaccu = fs * precip
    END SUBROUTINE model_snowaccumulation

END MODULE

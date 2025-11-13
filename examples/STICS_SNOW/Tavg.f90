MODULE Tavgmod
    IMPLICIT NONE
CONTAINS

    SUBROUTINE model_tavg(tmin, &
        tmax, &
        tavg)
        IMPLICIT NONE
        REAL, INTENT(IN) :: tmin
        REAL, INTENT(IN) :: tmax
        REAL, INTENT(OUT) :: tavg
        !- Name: Tavg -Version: 1.0, -Time step: 1
        !- Description:
    !            * Title: Mean temperature  calculation
    !            * Author: STICS
    !            * Reference: doi:http://dx.doi.org/10.1016/j.agrformet.2014.05.002
    !            * Institution: INRA
    !            * Abstract: It simulates the depth of snow cover and recalculate weather data
        !- inputs:
    !            * name: tmin
    !                          ** description : current minimum air temperature
    !                          ** inputtype : variable
    !                          ** variablecategory : auxiliary
    !                          ** datatype : DOUBLE
    !                          ** default : 0.0
    !                          ** min : 0.0
    !                          ** max : 500.0
    !                          ** unit : degC
    !                          ** uri : 
    !            * name: tmax
    !                          ** description : current maximum air temperature
    !                          ** inputtype : variable
    !                          ** variablecategory : auxiliary
    !                          ** datatype : DOUBLE
    !                          ** default : 0.0
    !                          ** min : 0.0
    !                          ** max : 100.0
    !                          ** unit : degC
    !                          ** uri : 
        !- outputs:
    !            * name: tavg
    !                          ** description : mean temperature
    !                          ** variablecategory : auxiliary
    !                          ** datatype : DOUBLE
    !                          ** min : 0.0
    !                          ** max : 500.0
    !                          ** unit : degC
    !                          ** uri : 
        tavg = (tmin + tmax) / 2
    END SUBROUTINE model_tavg

END MODULE

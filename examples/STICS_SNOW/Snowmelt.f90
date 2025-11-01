MODULE Snowmeltmod
    IMPLICIT NONE
CONTAINS

    SUBROUTINE model_snowmelt(ps, &
        M, &
        Snowmelt)
        IMPLICIT NONE
        REAL, INTENT(IN) :: ps
        REAL, INTENT(IN) :: M
        REAL, INTENT(OUT) :: Snowmelt
        !- Name: SnowMelt -Version: 1.0, -Time step: 1
        !- Description:
    !            * Title: Snow Melt
    !            * Author: STICS
    !            * Reference: doi:http://dx.doi.org/10.1016/j.agrformet.2014.05.002
    !            * Institution: INRA
    !            * Abstract: Snow melt
        !- inputs:
    !            * name: ps
    !                          ** description : density of snow cover
    !                          ** inputtype : variable
    !                          ** variablecategory : state
    !                          ** datatype : DOUBLE
    !                          ** default : 0.0
    !                          ** min : 0.0
    !                          ** max : 5000.0
    !                          ** unit : kg/m**3
    !                          ** uri : 
    !            * name: M
    !                          ** description : snow in the process of melting
    !                          ** inputtype : variable
    !                          ** variablecategory : rate
    !                          ** datatype : DOUBLE
    !                          ** default : 0.0
    !                          ** min : 
    !                          ** max : 
    !                          ** unit : mmW/d
    !                          ** uri : 
        !- outputs:
    !            * name: Snowmelt
    !                          ** description : Snow melt
    !                          ** variablecategory : state
    !                          ** datatype : DOUBLE
    !                          ** min : 0.0
    !                          ** max : 500.0
    !                          ** unit : m
    !                          ** uri : 
        Snowmelt = 0.0
        IF(ps .GT. 1e-8) THEN
            Snowmelt = M / ps
        END IF
    END SUBROUTINE model_snowmelt

END MODULE

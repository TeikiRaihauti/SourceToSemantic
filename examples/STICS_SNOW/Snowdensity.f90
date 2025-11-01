MODULE Snowdensitymod
    IMPLICIT NONE
CONTAINS

    SUBROUTINE model_snowdensity(ps_t1, &
        Sdepth_t1, &
        Sdry_t1, &
        Swet_t1, &
        ps)
        IMPLICIT NONE
        REAL, INTENT(IN) :: ps_t1
        REAL, INTENT(IN) :: Sdepth_t1
        REAL, INTENT(IN) :: Sdry_t1
        REAL, INTENT(IN) :: Swet_t1
        REAL, INTENT(OUT) :: ps
        !- Name: SnowDensity -Version: 1.0, -Time step: 1
        !- Description:
    !            * Title: Density of snow cover calculation
    !            * Author: STICS
    !            * Reference: doi:http://dx.doi.org/10.1016/j.agrformet.2014.05.002
    !            * Institution: INRA
    !            * Abstract: density of snow cover
        !- inputs:
    !            * name: ps_t1
    !                          ** description : density of snow cover in previous day
    !                          ** inputtype : variable
    !                          ** variablecategory : state
    !                          ** datatype : DOUBLE
    !                          ** default : 0.0
    !                          ** min : 0.0
    !                          ** max : 100.0
    !                          ** unit : kg/m**3
    !                          ** uri : 
    !            * name: Sdepth_t1
    !                          ** description : snow cover depth Calculation in previous day
    !                          ** inputtype : variable
    !                          ** variablecategory : state
    !                          ** datatype : DOUBLE
    !                          ** default : 0.0
    !                          ** min : 0.0
    !                          ** max : 5000.0
    !                          ** unit : m
    !                          ** uri : 
    !            * name: Sdry_t1
    !                          ** description : water in solid state in the snow cover in previous day
    !                          ** inputtype : variable
    !                          ** variablecategory : state
    !                          ** datatype : DOUBLE
    !                          ** default : 0.0
    !                          ** min : 0.0
    !                          ** max : 500.0
    !                          ** unit : mmW
    !                          ** uri : 
    !            * name: Swet_t1
    !                          ** description : water in liquid state in the snow cover
    !                          ** inputtype : variable
    !                          ** variablecategory : state
    !                          ** datatype : DOUBLE
    !                          ** default : 0.0
    !                          ** min : 0.0
    !                          ** max : 100.0
    !                          ** unit : mmW
    !                          ** uri : 
        !- outputs:
    !            * name: ps
    !                          ** description : density of snow cover
    !                          ** variablecategory : state
    !                          ** datatype : DOUBLE
    !                          ** min : 0.0
    !                          ** max : 500.0
    !                          ** unit : kg/m**3
    !                          ** uri : 
        ps = 0.0
        IF(ABS(Sdepth_t1) .GT. 0.0) THEN
            IF(ABS(Sdry_t1 + Swet_t1) .GT. 0.0) THEN
                ps = (Sdry_t1 + Swet_t1) / Sdepth_t1
            ELSE
                ps = ps_t1
            END IF
        END IF
    END SUBROUTINE model_snowdensity

END MODULE

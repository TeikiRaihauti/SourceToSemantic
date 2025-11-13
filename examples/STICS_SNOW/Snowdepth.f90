MODULE Snowdepthmod
    IMPLICIT NONE
CONTAINS

    SUBROUTINE model_snowdepth(Snowmelt, &
        Sdepth_t1, &
        Snowaccu, &
        E, &
        rho, &
        Sdepth)
        IMPLICIT NONE
        REAL, INTENT(IN) :: Snowmelt
        REAL, INTENT(IN) :: Sdepth_t1
        REAL, INTENT(IN) :: Snowaccu
        REAL, INTENT(IN) :: E
        REAL, INTENT(IN) :: rho
        REAL, INTENT(OUT) :: Sdepth
        !- Name: SnowDepth -Version: 1.0, -Time step: 1
        !- Description:
    !            * Title: snow cover depth Calculation
    !            * Author: STICS
    !            * Reference: doi:http://dx.doi.org/10.1016/j.agrformet.2014.05.002
    !            * Institution: INRA
    !            * Abstract: snow cover depth Calculation
        !- inputs:
    !            * name: Snowmelt
    !                          ** description : snow melt 
    !                          ** inputtype : variable
    !                          ** variablecategory : state
    !                          ** datatype : DOUBLE
    !                          ** default : 0.0
    !                          ** min : 0.0
    !                          ** max : 100.0
    !                          ** unit : m
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
    !            * name: Snowaccu
    !                          ** description : snowfall accumulation
    !                          ** inputtype : variable
    !                          ** variablecategory : rate
    !                          ** datatype : DOUBLE
    !                          ** default : 0.0
    !                          ** min : 
    !                          ** max : 
    !                          ** unit : mmW/d
    !                          ** uri : 
    !            * name: E
    !                          ** description : snow compaction parameter
    !                          ** inputtype : parameter
    !                          ** parametercategory : constant
    !                          ** datatype : DOUBLE
    !                          ** default : 0.0
    !                          ** min : 
    !                          ** max : 
    !                          ** unit : mm/mm/d
    !                          ** uri : 
    !            * name: rho
    !                          ** description : The density of the new snow fixed by the user
    !                          ** inputtype : parameter
    !                          ** parametercategory : constant
    !                          ** datatype : DOUBLE
    !                          ** default : 100
    !                          ** min : 
    !                          ** max : 
    !                          ** unit : kg/m**3
    !                          ** uri : 
        !- outputs:
    !            * name: Sdepth
    !                          ** description : snow cover depth Calculation
    !                          ** variablecategory : state
    !                          ** datatype : DOUBLE
    !                          ** min : 0.0
    !                          ** max : 500.0
    !                          ** unit : m
    !                          ** uri : 
        Sdepth = 0.0
        IF(Snowmelt .LE. Sdepth_t1 + (Snowaccu / rho)) THEN
            Sdepth = Snowaccu / rho + Sdepth_t1 - Snowmelt - (Sdepth_t1 * E)
        END IF
    END SUBROUTINE model_snowdepth

END MODULE

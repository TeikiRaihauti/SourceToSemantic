MODULE Preciprecmod
    IMPLICIT NONE
CONTAINS

    SUBROUTINE model_preciprec(Sdry_t1, &
        Sdry, &
        Swet, &
        Swet_t1, &
        Sdepth_t1, &
        Sdepth, &
        Mrf, &
        precip, &
        Snowaccu, &
        rho, &
        preciprec)
        IMPLICIT NONE
        REAL, INTENT(IN) :: Sdry_t1
        REAL, INTENT(IN) :: Sdry
        REAL, INTENT(IN) :: Swet
        REAL, INTENT(IN) :: Swet_t1
        REAL, INTENT(IN) :: Sdepth_t1
        REAL, INTENT(IN) :: Sdepth
        REAL, INTENT(IN) :: Mrf
        REAL, INTENT(IN) :: precip
        REAL, INTENT(IN) :: Snowaccu
        REAL, INTENT(IN) :: rho
        REAL, INTENT(OUT) :: preciprec
        !- Name: Preciprec -Version: 1.0, -Time step: 1
        !- Description:
    !            * Title: Precipitation ReCalculation
    !            * Author: STICS
    !            * Reference: doi:http://dx.doi.org/10.1016/j.agrformet.2014.05.002
    !            * Institution: INRA
    !            * Abstract: recalculation of precipitation
        !- inputs:
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
    !            * name: Sdry
    !                          ** description : water in solid state in the snow cover 
    !                          ** inputtype : variable
    !                          ** variablecategory : state
    !                          ** datatype : DOUBLE
    !                          ** default : 0.0
    !                          ** min : 0.0
    !                          ** max : 500.0
    !                          ** unit : mmW
    !                          ** uri : 
    !            * name: Swet
    !                          ** description : water in liquid state in the snow cover
    !                          ** inputtype : variable
    !                          ** variablecategory : state
    !                          ** datatype : DOUBLE
    !                          ** default : 0.0
    !                          ** min : 0.0
    !                          ** max : 100.0
    !                          ** unit : mmW
    !                          ** uri : 
    !            * name: Swet_t1
    !                          ** description : water in liquid state in the snow cover in previous day
    !                          ** inputtype : variable
    !                          ** variablecategory : state
    !                          ** datatype : DOUBLE
    !                          ** default : 0.0
    !                          ** min : 0.0
    !                          ** max : 100.0
    !                          ** unit : mmW
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
    !            * name: Sdepth
    !                          ** description : snow cover depth Calculation
    !                          ** inputtype : variable
    !                          ** variablecategory : state
    !                          ** datatype : DOUBLE
    !                          ** default : 0.0
    !                          ** min : 0.0
    !                          ** max : 5000.0
    !                          ** unit : m
    !                          ** uri : 
    !            * name: Mrf
    !                          ** description : liquid water in the snow cover in the process of refreezing
    !                          ** inputtype : variable
    !                          ** variablecategory : rate
    !                          ** datatype : DOUBLE
    !                          ** default : 0.0
    !                          ** min : 
    !                          ** max : 
    !                          ** unit : mmW/d
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
    !            * name: preciprec
    !                          ** description : recalculated precipitation
    !                          ** variablecategory : state
    !                          ** datatype : DOUBLE
    !                          ** min : 0.0
    !                          ** max : 500.0
    !                          ** unit : mmW
    !                          ** uri : 
        preciprec = precip
        IF(Sdry + Swet .LT. Sdry_t1 + Swet_t1) THEN
            preciprec = preciprec + ((Sdepth_t1 - Sdepth) * rho) - Mrf
        END IF
        preciprec = preciprec - Snowaccu
    END SUBROUTINE model_preciprec

END MODULE

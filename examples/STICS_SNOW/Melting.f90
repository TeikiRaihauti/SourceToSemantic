MODULE Meltingmod
    IMPLICIT NONE
CONTAINS

    SUBROUTINE model_melting(jul, &
        Tmf, &
        DKmax, &
        Kmin, &
        tavg, &
        M)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: jul
        REAL, INTENT(IN) :: Tmf
        REAL, INTENT(IN) :: DKmax
        REAL, INTENT(IN) :: Kmin
        REAL, INTENT(IN) :: tavg
        REAL, INTENT(OUT) :: M
        REAL:: K
        !- Name: Melting -Version: 1.0, -Time step: 1
        !- Description:
    !            * Title: snow in the process of melting
    !            * Author: STICS
    !            * Reference: doi:http://dx.doi.org/10.1016/j.agrformet.2014.05.002
    !            * Institution: INRA
    !            * Abstract: It simulates snow in the process of melting
        !- inputs:
    !            * name: jul
    !                          ** description : current day of year for the calculation
    !                          ** inputtype : variable
    !                          ** variablecategory : auxiliary
    !                          ** datatype : INT
    !                          ** default : 0
    !                          ** min : 0
    !                          ** max : 366
    !                          ** unit : d
    !                          ** uri : 
    !            * name: Tmf
    !                          ** description : threshold temperature for snow melting 
    !                          ** inputtype : parameter
    !                          ** parametercategory : constant
    !                          ** datatype : DOUBLE
    !                          ** default : 0.5
    !                          ** min : 0.0
    !                          ** max : 1.0
    !                          ** unit : degC
    !                          ** uri : 
    !            * name: DKmax
    !                          ** description : difference between the maximum and the minimum melting rates
    !                          ** inputtype : parameter
    !                          ** parametercategory : constant
    !                          ** datatype : DOUBLE
    !                          ** default : 0.0
    !                          ** min : 0.0
    !                          ** max : 5000.0
    !                          ** unit : mmW/degC/d
    !                          ** uri : 
    !            * name: Kmin
    !                          ** description : minimum melting rate on 21 December
    !                          ** inputtype : parameter
    !                          ** parametercategory : constant
    !                          ** datatype : DOUBLE
    !                          ** default : 0.0
    !                          ** min : 0.0
    !                          ** max : 5000.0
    !                          ** unit : mmW/degC/d
    !                          ** uri : 
    !            * name: tavg
    !                          ** description : average temperature
    !                          ** inputtype : variable
    !                          ** variablecategory : auxiliary
    !                          ** datatype : DOUBLE
    !                          ** default : 0.0
    !                          ** min : 0.0
    !                          ** max : 100.0
    !                          ** unit : degC
    !                          ** uri : 
        !- outputs:
    !            * name: M
    !                          ** description : snow in the process of melting
    !                          ** variablecategory : rate
    !                          ** datatype : DOUBLE
    !                          ** min : 0.0
    !                          ** max : 500.0
    !                          ** unit : mmW/d
    !                          ** uri : 
        M = 0.0
        K = DKmax / 2.0 * (-SIN((2.0 * 3.14159265 * REAL(jul) / 366.0 + (9.0  &
                / 16.0 * 3.14159265)))) + Kmin + (DKmax / 2.0)
        IF(tavg .GT. Tmf) THEN
            M = K * (tavg - Tmf)
        END IF
    END SUBROUTINE model_melting

END MODULE

MODULE Tempmaxmod
    IMPLICIT NONE
CONTAINS

    SUBROUTINE model_tempmax(Sdepth_cm, &
        prof, &
        tmax, &
        tminseuil, &
        tmaxseuil, &
        tmaxrec)
        IMPLICIT NONE
        REAL, INTENT(IN) :: Sdepth_cm
        REAL, INTENT(IN) :: prof
        REAL, INTENT(IN) :: tmax
        REAL, INTENT(IN) :: tminseuil
        REAL, INTENT(IN) :: tmaxseuil
        REAL, INTENT(OUT) :: tmaxrec
        !- Name: TempMax -Version: 1.0, -Time step: 1
        !- Description:
    !            * Title: Maximum temperature  recalculation
    !            * Author: STICS
    !            * Reference: doi:http://dx.doi.org/10.1016/j.agrformet.2014.05.002
    !            * Institution: INRA
    !            * Abstract: recalculation of maximum temperature
        !- inputs:
    !            * name: Sdepth_cm
    !                          ** description : snow depth
    !                          ** inputtype : variable
    !                          ** variablecategory : state
    !                          ** datatype : DOUBLE
    !                          ** default : 0.0
    !                          ** min : 0.0
    !                          ** max : 500.0
    !                          ** unit : cm
    !                          ** uri : 
    !            * name: prof
    !                          ** description : snow cover threshold for snow insulation 
    !                          ** inputtype : parameter
    !                          ** parametercategory : constant
    !                          ** datatype : DOUBLE
    !                          ** default : 0.0
    !                          ** min : 0.0
    !                          ** max : 1000
    !                          ** unit : cm
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
    !            * name: tminseuil
    !                          ** description : minimum temperature when snow cover is higher than prof
    !                          ** inputtype : parameter
    !                          ** parametercategory : constant
    !                          ** datatype : DOUBLE
    !                          ** default : 0.0
    !                          ** min : 0.0
    !                          ** max : 5000.0
    !                          ** unit : degC
    !                          ** uri : 
    !            * name: tmaxseuil
    !                          ** description : maximum temperature when snow cover is higher than prof
    !                          ** inputtype : parameter
    !                          ** parametercategory : constant
    !                          ** datatype : DOUBLE
    !                          ** default : 0.0
    !                          ** min : 
    !                          ** max : 
    !                          ** unit : degC
    !                          ** uri : 
        !- outputs:
    !            * name: tmaxrec
    !                          ** description : recalculated maximum temperature
    !                          ** variablecategory : state
    !                          ** datatype : DOUBLE
    !                          ** min : 0.0
    !                          ** max : 500.0
    !                          ** unit : degC
    !                          ** uri : 
        tmaxrec = tmax
        IF(Sdepth_cm .GT. prof) THEN
            IF(tmax .LT. tminseuil) THEN
                tmaxrec = tminseuil
            ELSE
                IF(tmax .GT. tmaxseuil) THEN
                    tmaxrec = tmaxseuil
                END IF
            END IF
        ELSE
            IF(Sdepth_cm .GT. 0.0) THEN
                IF(tmax .LE. 0.0) THEN
                    tmaxrec = tmaxseuil - ((1 - (Sdepth_cm / prof)) * (-tmax))
                ELSE
                    tmaxrec = 0.0
                END IF
            END IF
        END IF
    END SUBROUTINE model_tempmax

END MODULE

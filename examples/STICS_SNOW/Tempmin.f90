MODULE Tempminmod
    IMPLICIT NONE
CONTAINS

    SUBROUTINE model_tempmin(Sdepth_cm, &
        prof, &
        tmin, &
        tminseuil, &
        tmaxseuil, &
        tminrec)
        IMPLICIT NONE
        REAL, INTENT(IN) :: Sdepth_cm
        REAL, INTENT(IN) :: prof
        REAL, INTENT(IN) :: tmin
        REAL, INTENT(IN) :: tminseuil
        REAL, INTENT(IN) :: tmaxseuil
        REAL, INTENT(OUT) :: tminrec
        !- Name: TempMin -Version: 1.0, -Time step: 1
        !- Description:
    !            * Title: Minimum temperature  calculation
    !            * Author: STICS
    !            * Reference: doi:http://dx.doi.org/10.1016/j.agrformet.2014.05.002
    !            * Institution: INRA
    !            * Abstract: recalculation of minimum temperature
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
    !            * name: tmin
    !                          ** description : current minimum air temperature
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
    !            * name: tminrec
    !                          ** description : recalculated minimum temperature
    !                          ** variablecategory : state
    !                          ** datatype : DOUBLE
    !                          ** min : 0.0
    !                          ** max : 500.0
    !                          ** unit : degC
    !                          ** uri : 
        tminrec = tmin
        IF(Sdepth_cm .GT. prof) THEN
            IF(tmin .LT. tminseuil) THEN
                tminrec = tminseuil
            ELSE
                IF(tmin .GT. tmaxseuil) THEN
                    tminrec = tmaxseuil
                END IF
            END IF
        ELSE
            IF(Sdepth_cm .GT. 0.0) THEN
                tminrec = tminseuil - ((1 - (Sdepth_cm / prof)) * (ABS(tmin) +  &
                        tminseuil))
            END IF
        END IF
    END SUBROUTINE model_tempmin

END MODULE

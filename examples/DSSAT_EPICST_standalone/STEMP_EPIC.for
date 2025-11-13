C=======================================================================
C  COPYRIGHT 1998-2010 The University of Georgia, Griffin, Georgia
C                      University of Florida, Gainesville, Florida
C                      Iowa State University, Ames, Iowa
C                      International Center for Soil Fertility and 
C                       Agricultural Development, Muscle Shoals, Alabama
C                      University of Guelph, Guelph, Ontario
C  ALL RIGHTS RESERVED
C=======================================================================
C=======================================================================
C  STEMP_EPIC, Subroutine
C  
C  @ShortDescription
C  Determines soil temperature by layer
C  test encore
C
C  @ExtendedDescription
C
C  @Timestep 1day
C  @Version 1.0
C-----------------------------------------------------------------------
C  @Authors
C  Kenneth N. Potter (USDA-ARS)
C  Jimmy R. Williams (USDA-ARS)
C  
C  @Reference
C  https://doi.org/10.2134/agronj1994.00021962008600060014x
C------------------------------------------------------------------------
C  Revision history
C  12/01/1980     Originally based on EPIC soil temperature routines
!  09/16/2010 CHP / MSC modified for EPIC soil temperature method.
!     Cite Potter & Williams here
C-----------------------------------------------------------------------
C  Called : Main
C  Calls  : SOILT
C=======================================================================
      SUBROUTINE STEMP_EPIC(CONTROL, ISWITCH,  
     &    SOILPROP,      !Input
     &    SW,            !Input
     &    TAVG,          !Input
     &    TMAX,          !Input
     &    TMIN,          !Input
     &    TAV,           !Input
     &    WEATHER,       !Input
     &    CUMDPT,        !InOut
     &    DSMID,         !InOut
     &    TDL,           !InOut
     &    TMA,           !InOut
     &    NDays,         !InOut
     &    WetDay,        !InOut
     &    X2_PREV,       !InOut
     &    SRFTEMP,       !InOut
     &    ST)            !InOut

!     &    RAIN,          !Input
!     &    TAMP,          !Input
!     &    ISWWAT,        !Input   
!     &    BD,            !Input
!     &    DLAYR,         !Input
!     &    DS,            !Input
!     &    DUL,           !Input
!     &    LL,            !Input
!     &    NLAYR,         !Input  
!     &    X2_AVG,        !InOut

C-----------------------------------------------------------------------
      USE ModuleDefs
      USE ModuleData

      IMPLICIT  NONE
      SAVE


      CHARACTER*1  RNMODE
      INTEGER DYNAMIC
      INTEGER RUN, YRDOY, YEAR
      INTEGER ERRNUM, FOUND, LNUM, LUNIO
      CHARACTER*6  SECTION
      CHARACTER*6, PARAMETER :: ERRKEY = "EPIC STEMP"
      CHARACTER*30 FILEIO
      CHARACTER*78 MSG(3)


      CHARACTER*1 ISWWAT
      INTEGER DOY, I, L, NLAYR
      INTEGER WetDay(30), NDays, NWetDays
      REAL ABD, B, CUMDPT 
      REAL DP, FX, ICWD, PESW, SRFTEMP 
      REAL TAV, TAMP, TAVG, TBD, TMAX, TMIN, WW
      REAL TDL, TLL, TSW
      REAL TMA(5), X2_AVG
      REAL DEPIR, WFT, BCV, RAIN, BIOMAS, MULCHMASS
      REAL SNOW, CV, BCV1, BCV2, X2_PREV
      REAL, DIMENSION(NL) :: BD, DLAYR, DS, DUL, LL, ST, SW, SWI, DSMID

!-----------------------------------------------------------------------
    
      TYPE (ControlType) CONTROL
   

      TYPE (SoilType)    SOILPROP
      TYPE (SwitchType)  ISWITCH
      TYPE (WeatherType) WEATHER


!     Transfer values from constructed data types into local variables.
      DYNAMIC = CONTROL % DYNAMIC  
      YRDOY   = CONTROL % YRDOY    
   

      ISWWAT = ISWITCH % ISWWAT
!      METMP  = ISWITCH % METMP

      BD     = SOILPROP % BD     
      DLAYR  = SOILPROP % DLAYR  
      DS     = SOILPROP % DS     
      DUL    = SOILPROP % DUL     
      LL     = SOILPROP % LL     
      NLAYR  = SOILPROP % NLAYR  

      TAMP = WEATHER % TAMP

!-----------------------------------------------------------------------
    
      CALL YR_DOY(YRDOY, YEAR, DOY)

!***********************************************************************
!***********************************************************************
!     Run initialization - run once per simulation
!***********************************************************************
!      IF (DYNAMIC .EQ. RUNINIT) THEN
!-----------------------------------------------------------------------
!***********************************************************************
!***********************************************************************
!     Seasonal initialization - run once per season
!***********************************************************************
!      ELSEIF (DYNAMIC .EQ. SEASINIT) THEN
      IF (DYNAMIC .EQ. SEASINIT) THEN
!-----------------------------------------------------------------------
      FILEIO  = CONTROL % FILEIO    
      LUNIO   = CONTROL % LUNIO    
      RUN     = CONTROL % RUN    
      RNMODE  = CONTROL % RNMODE

      IF (RUN .EQ. 1 .OR. INDEX('QF',RNMODE) .LE. 0) THEN

!        IF (ISWWAT .NE. 'N') THEN
!!         Read inital soil water values from FILEIO 
!!         (not yet done in WATBAL, so need to do here)
!          OPEN (LUNIO, FILE = FILEIO, STATUS = 'OLD', IOSTAT=ERRNUM)
!          IF (ERRNUM .NE. 0) CALL ERROR(ERRKEY,ERRNUM,FILEIO,0)
!          SECTION = '*INITI'
!          CALL FIND(LUNIO, SECTION, LNUM, FOUND) 
!          IF (FOUND .EQ. 0) CALL ERROR(SECTION, 42, FILEIO, LNUM)
!
!!         Initial depth to water table (not currently used)
!          READ(LUNIO,'(40X,F6.0)',IOSTAT=ERRNUM) ICWD ; LNUM = LNUM + 1
!          IF (ERRNUM .NE. 0) CALL ERROR(ERRKEY,ERRNUM,FILEIO,LNUM)
!
!          DO L = 1, NLAYR
!            READ(LUNIO,'(9X,F5.3)',IOSTAT=ERRNUM) SWI(L)
!            LNUM = LNUM + 1
!            IF (ERRNUM .NE. 0) CALL ERROR(ERRKEY,ERRNUM,FILEIO,LNUM)
!            IF (SWI(L) .LT. LL(L)) SWI(L) = LL(L)
!          ENDDO
!
!          CLOSE (LUNIO)
!        ELSE
!          SWI = DUL
!        ENDIF
        SWI = SW
        
        TBD = 0.0
        TLL = 0.0
        TSW = 0.0
        TDL = 0.0
        CUMDPT = 0.0
        DO L = 1, NLAYR
          DSMID(L) = CUMDPT + DLAYR(L)* 5.0   !mm depth to midpt of lyr
          CUMDPT   = CUMDPT + DLAYR(L)*10.0   !mm profile depth 
          TBD = TBD + BD(L)  * DLAYR(L)       
          TLL = TLL + LL(L)  * DLAYR(L)
          TSW = TSW + SWI(L) * DLAYR(L)
          TDL = TDL + DUL(L) * DLAYR(L)
        END DO

        IF (ISWWAT .EQ. 'Y') THEN
          PESW = AMAX1(0.0, TSW - TLL)      !cm
        ELSE
          !If water not being simulated, use DUL as water content
          PESW = AMAX1(0.0, TDL - TLL)
        ENDIF

        ABD    = TBD / DS(NLAYR)                   !CHP
        FX     = ABD/(ABD+686.0*EXP(-5.63*ABD))
        DP     = 1000.0 + 2500.0*FX
        WW     = 0.356  - 0.144*ABD
        B      = ALOG(500.0/DP)

        DO I = 1, 5
          TMA(I) = NINT(TAVG*10000.)/10000.   !chp
        END DO
        X2_AVG = TMA(1) * 5.0
 
        DO L = 1, NLAYR
          ST(L) = TAVG
        END DO

!       Save 30 day memory of:
!       WFT = fraction of wet days (rainfall + irrigation)
        WFT = 0.1
        WetDay = 0
        NDays = 0

!       Soil cover function
        CALL GET('ORGC' ,'MULCHMASS',MULCHMASS)   !kg/ha
        CALL GET('WATER','SNOW'     , SNOW)       !mm
      
        CV = (MULCHMASS) / 1000.         !t/ha
        BCV1 = CV / (CV + EXP(5.3396 - 2.3951 * CV))
        BCV2 = SNOW / (SNOW + EXP(2.303 - 0.2197 * SNOW))
        BCV = MAX(BCV1, BCV2)

        DO I = 1, 8
          CALL SOILT_EPIC (
     &    B, BCV, CUMDPT, DP, DSMID, NLAYR, PESW, TAV,    !Input
     &    TAVG, TMAX, TMIN, 0, WFT, WW,                   !Input
     &    TMA, SRFTEMP, ST, X2_AVG, X2_PREV)              !Output
        END DO
      ENDIF

!     Print soil temperature data in STEMP.OUT
      CALL OPSTEMP(CONTROL, ISWITCH, DOY, SRFTEMP, ST, TAV, TAMP)

      MSG(1) = "Running EPIC soil temperature routine."
      MSG(2) = "Start simulation at least 30 days early to initialize"
      MSG(3) = "  soil temperature parameters."
      !CALL WARNING(3,ERRKEY,MSG)
!***********************************************************************
!***********************************************************************
!     Daily rate calculations
!***********************************************************************
      ELSEIF (DYNAMIC .EQ. RATE) THEN
!-----------------------------------------------------------------------
      TBD = 0.0
      TLL = 0.0
      TSW = 0.0
      DO L = 1, NLAYR
        TBD = TBD + BD(L) * DLAYR(L) 
        TDL = TDL + DUL(L)* DLAYR(L)
        TLL = TLL + LL(L) * DLAYR(L)
        TSW = TSW + SW(L) * DLAYR(L)
      ENDDO

      ABD    = TBD / DS(NLAYR)                    !CHP
      FX     = ABD/(ABD+686.0*EXP(-5.63*ABD))
      DP     = 1000.0 + 2500.0*FX   !DP in mm
      WW     = 0.356  - 0.144*ABD   !vol. fraction
      B      = ALOG(500.0/DP)

      IF (ISWWAT .EQ. 'Y') THEN
        PESW = MAX(0.0, TSW - TLL)      !cm
      ELSE
        !If water not being simulated, use DUL as water content
        PESW = AMAX1(0.0, TDL - TLL)    !cm
      ENDIF

!     Save 30 day memory of:
!     WFT = fraction of wet days (rainfall + irrigation)
      RAIN = WEATHER % RAIN
      CALL GET('MGMT','DEPIR',DEPIR)
      IF (NDays == 30) THEN
        DO I = 1, 29
          WetDay(I) = WetDay(I+1)
        ENDDO
      ELSE
        NDays = NDays + 1
      ENDIF
      IF (RAIN + DEPIR > 1.E-6) THEN
        WetDay(NDays) = 1
      ELSE
        WetDay(NDays) = 0 
      ENDIF
      NWetDays = SUM(WetDay)
      WFT = Float(NWetDays) / float(NDays)

!     Soil cover function
      CALL GET('PLANT','BIOMAS'   ,BIOMAS)      !kg/ha
      CALL GET('ORGC' ,'MULCHMASS',MULCHMASS)   !kg/ha
      CALL GET('WATER','SNOW'     , SNOW)       !mm
      CV = (BIOMAS + MULCHMASS) / 1000.         !t/ha
      BCV1 = CV / (CV + EXP(5.3396 - 2.3951 * CV))
      BCV2 = SNOW / (SNOW + EXP(2.303 - 0.2197 * SNOW))
      BCV = MAX(BCV1, BCV2)

      CALL SOILT_EPIC (
     &    B, BCV, CUMDPT, DP, DSMID, NLAYR, PESW, TAV,    !Input
     &    TAVG, TMAX, TMIN, WetDay(NDays), WFT, WW,       !Input
     &    TMA, SRFTEMP, ST, X2_AVG, X2_PREV)              !Output     
!***********************************************************************
!***********************************************************************
!     Output & Seasonal summary
!***********************************************************************
      ELSEIF (DYNAMIC .EQ. OUTPUT .OR. DYNAMIC .EQ. SEASEND) THEN
!-----------------------------------------------------------------------
      CALL OPSTEMP(CONTROL, ISWITCH, DOY, SRFTEMP, ST, TAV, TAMP)
!***********************************************************************
!***********************************************************************
!     END OF DYNAMIC IF CONSTRUCT
!***********************************************************************
      ENDIF
!***********************************************************************
      RETURN
      END SUBROUTINE STEMP_EPIC
!=======================================================================

C=======================================================================
C  SOILT_EPIC, Subroutine
C  Determines soil temperature by layer
C-----------------------------------------------------------------------
C  Revision history
C  02/09/1933 PWW Header revision and minor changes.
C  12/09/1999 CHP Revisions for modular format.
C  01/01/2000 AJG Added surface temperature for the CENTURY-based
C                SOM/soil-N module.
C  01/14/2005 CHP Added METMP = 3: Corrected water content in temp. eqn.
!  12/07/2008 CHP Removed METMP -- use only corrected water content
!  09/16/2010 CHP / MSC modified for EPIC soil temperature method.
C-----------------------------------------------------------------------
C  Called : STEMP
C  Calls  : None
C=======================================================================

      SUBROUTINE SOILT_EPIC (
     &    B, BCV, CUMDPT, DP, DSMID, NLAYR, PESW, TAV,    !Input
     &    TAVG, TMAX, TMIN, WetDay, WFT, WW,              !Input
     &    TMA,                                            !Inout
     &    SRFTEMP,                                        !Output
     &    ST,                                             !InOut
     &    X2_AVG,                                         !Output
     &    X2_PREV)                                        !InOut

      
!     ------------------------------------------------------------------
      USE ModuleDefs
      IMPLICIT  NONE
      SAVE

      INTEGER  K, L, NLAYR, WetDay

      REAL B, CUMDPT, DD, DP, FX
      REAL PESW, SRFTEMP, TAV, TAVG, TMAX
      REAL WC, WW, ZD
      REAL TMA(5)
      REAL DSMID(NL) 
      REAL ST(NL)
      REAL X1, X2, X3, F, WFT, BCV, TMIN, X2_AVG, X2_PREV
      REAL LAG
      PARAMETER (LAG=0.5)

!-----------------------------------------------------------------------
      WC = AMAX1(0.01, PESW) / (WW * CUMDPT) * 10.0
!     frac =              cm   / (    mm     ) * mm/cm
        !WC (ratio)
        !PESW (cm)
        !WW (dimensionless)
        !CUMDPT (mm)

      FX = EXP(B * ((1.0 - WC) / (1.0 + WC))**2)
      DD = FX * DP                                  !DD in mm

!=========================================================================
!     Below this point - EPIC soil temperature routine differs from
!       DSSAT original routine.
!=========================================================================

      IF (WetDay > 0) THEN
!       Potter & Williams, 1994, Eqn. 2
!       X2=WFT(MO)*(TX-TMN)+TMN
        X2=WFT*(TAVG-TMIN)+TMIN
      ELSE
!       Eqn 1
!       X2=WFT(MO)*(TMX-TX)+TX+2.*((ST0/15.)**2-1.)
!       Removed ST0 factor for now.
        X2=WFT*(TMAX-TAVG)+TAVG+2.
      ENDIF

      TMA(1) = X2
      DO K = 5, 2, -1
        TMA(K) = TMA(K-1)
      END DO
      X2_AVG = SUM(TMA) / 5.0     !Eqn 

!     Eqn 4
!     X3=(1.-BCV)*X2+BCV*T(LID(2))
      X3=(1.-BCV)*X2_AVG+BCV*X2_PREV

!     DST0=AMIN1(X2,X3)
      SRFTEMP=AMIN1(X2_AVG,X3)

!     Eqn 6 (partial)
!     X1=AVT-X3
      X1=TAV-X3

      DO L = 1, NLAYR
        ZD    = DSMID(L) / DD  !Eqn 8
!       Eqn 7
        F=ZD/(ZD+EXP(-.8669-2.0775*ZD))
!       Eqn 6
!       T(L)=PARM(15)*T(L)+XLG1*(F*X1+X3)
        ST(L)=LAG*ST(L)+(1.-LAG)*(F*X1+X3)
      END DO

      X2_PREV = X2_AVG

!=========================================================================
!     old CSM code:
!=========================================================================
!
!      TA = TAV + TAMP * COS(ALX) / 2.0
!      DT = ATOT / 5.0 - TA
!
!      DO L = 1, NLAYR
!        ZD    = -DSMID(L) / DD
!        ST(L) = TAV + (TAMP / 2.0 * COS(ALX + ZD) + DT) * EXP(ZD)
!        ST(L) = NINT(ST(L) * 1000.) / 1000.   !debug vs release fix
!      END DO
!
!-----------------------------------------------------------------------

      RETURN
      END SUBROUTINE SOILT_EPIC
C=======================================================================
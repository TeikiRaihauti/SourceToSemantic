MODULE temp_profile_mod
   IMPLICIT NONE

CONTAINS
   SUBROUTINE model_temp_profile(temp_amp, &
                               prev_temp_profile, &
                               prev_canopy_temp, &
                               min_air_temp, &
                               temp_profile)

      IMPLICIT NONE

      REAL, INTENT(IN)  :: temp_amp
      REAL, INTENT(IN)  :: prev_temp_profile(:)
      REAL, INTENT(IN)  :: prev_canopy_temp
      REAL, INTENT(IN)  :: min_air_temp
      REAL, allocatable, INTENT(OUT) :: temp_profile(:)


      INTEGER :: z,n
      REAL, allocatable :: vexp(:)
      REAL :: therm_diff = 5.37e-3
      REAL :: temp_freq = 7.272e-5
      REAL :: therm_amp
      

      n = size(prev_temp_profile)

      !if (.NOT. ALLOCATED(temp_profile)) then
         allocate(temp_profile(n))
      !end if
      
      !if (.NOT. ALLOCATED(vexp)) then
         allocate(vexp(n))
      !end if

      therm_amp = sqrt(temp_freq/2/therm_diff)

      DO z=1, n
         vexp(z) = exp(-z*therm_amp)
      END DO

      DO z=1, n
         temp_profile(z) = prev_temp_profile(z) - &
                           vexp(z) * (prev_canopy_temp - min_air_temp) + &
                           0.1*(prev_canopy_temp - prev_temp_profile(z)) + &
                           (temp_amp * vexp(z))/2
      END DO
      
   END SUBROUTINE model_temp_profile


   
   SUBROUTINE init_temp_profile(air_temp_day1, layer_thick, &
                                prev_temp_profile, prev_canopy_temp)
      
      real, intent(IN) :: air_temp_day1
      INTEGER, INTENT(IN) ::  layer_thick(:)
      real, intent(OUT) :: prev_canopy_temp
      real, allocatable, intent(OUT) :: prev_temp_profile(:)
      

      integer :: soil_depth  
      soil_depth = sum(layer_thick)
      allocate(prev_temp_profile(soil_depth))
      prev_temp_profile = air_temp_day1
      prev_canopy_temp = air_temp_day1
   END SUBROUTINE

END MODULE temp_profile_mod


!
!__________________________________________________________________________________________!
! This module contains the Predicted Particle Property (P3) bulk microphysics scheme.      !
!                                                                                          !
! This code was originally written by H. Morrison, MMM Division, NCAR (Dec 2012).          !
! Modifications were made by J. Milbrandt, RPN, Environment Canada (July 2014).            !
! Subsequent major and minor upgrades have been ongoing.                                   !
!                                                                                          !
! For model-specific aspects/versions, see comments in the interface subroutine(s) in      !
!   this module (mp_p3_wrapper_wrfcm1, mp_p3_wrapper_gem).                                 !
!                                                                                          !
! For details see:                                                                         !
!   Morrison and Milbrandt (2015) [J. Atmos. Sci., 72, 287-311]    - original scheme desc. !
!   Milbrandt and Morrison (2016) [J. Atmos. Sci., 73, 975-995]    - multi-ice-category    !
!   Cholette et al. (2019)        [J. Atmos. Sci., 76, 561-582]    - liquid fraction ice   !
!   Jouan et al. (2020)           [W. Forecasting, 35, 2541-2565]  - cloud fraction        !
!   Milbrandt et al. (2021)       [J. Atmos. Sci., 78, 439-458]    - triple-moment ice     !
!   Cholette et al. (2023)        [J.A.M.E.S, 15(4), e2022MS003328 - trplMomIce + liqFrac  !
!   Morrison et al. (2025)        [J.A.M.E.S, 17, e2024MS004644    - full trplMomIce       !
!                                                                                          !
! For questions or bug reports, please contact:                                            !
!    Hugh Morrison   (morrison@ucar.edu), or                                               !
!    Jason Milbrandt (jason.milbrandt@ec.gc.ca), or                                        !
!    Melissa Cholette (melissa.cholette@ec.gc.ca)                                          !
!                                                                                          !
! For all code updates, including bug-fixes, new developments, and test summaries, see:    !
!    https://github.com/P3-microphysics/P3-microphysics                                    !
!__________________________________________________________________________________________!
!                                                                                          !
! Version:       5.5.1                                                                     !
! Last updated:  2026 Feb                                                                  !
!__________________________________________________________________________________________!

 MODULE microphy_qv_sat

 implicit none

 public :: polysvp1, qv_sat

 contains


!==========================================================================================!

 real function polysvp1(T,i_type)

!-------------------------------------------
!  COMPUTE SATURATION VAPOR PRESSURE
!  POLYSVP1 RETURNED IN UNITS OF PA.
!  T IS INPUT IN UNITS OF K.
!  i_type REFERS TO SATURATION WITH RESPECT TO LIQUID (0) OR ICE (1)
!-------------------------------------------

      implicit none

      real    :: T
      integer :: i_type
      real    :: trplpt
! REPLACE GOFF-GRATCH WITH FASTER FORMULATION FROM FLATAU ET AL. 1992, TABLE 4 (RIGHT-HAND COLUMN)

! ice
      real a0i,a1i,a2i,a3i,a4i,a5i,a6i,a7i,a8i
      data a0i,a1i,a2i,a3i,a4i,a5i,a6i,a7i,a8i /&
        6.11147274, 0.503160820, 0.188439774e-1, &
        0.420895665e-3, 0.615021634e-5,0.602588177e-7, &
        0.385852041e-9, 0.146898966e-11, 0.252751365e-14/

! liquid
      real a0,a1,a2,a3,a4,a5,a6,a7,a8

! V1.7
      data a0,a1,a2,a3,a4,a5,a6,a7,a8 /&
        6.11239921, 0.443987641, 0.142986287e-1, &
        0.264847430e-3, 0.302950461e-5, 0.206739458e-7, &
        0.640689451e-10,-0.952447341e-13,-0.976195544e-15/
      real dt

!-------------------------------------------
      trplpt = 273.15

      if (i_type.EQ.1 .and. T.lt.trplpt) then
! ICE

! use Goff-Gratch for T < 195.8 K and Flatau et al. equal or above 195.8 K
         if (t.ge.195.8) then
            dt=t-trplpt
            polysvp1 = a0i + dt*(a1i+dt*(a2i+dt*(a3i+dt*(a4i+dt*(a5i+dt*(a6i+dt*(a7i+a8i*dt)))))))
            polysvp1 = polysvp1*100.
         else
            polysvp1 = 10.**(-9.09718*(273.16/t-1.)-3.56654* &
                alog10(273.16/t)+0.876793*(1.-t/273.16)+ &
                alog10(6.1071))*100.
         end if

      elseif (i_type.EQ.0 .or. T.ge.trplpt) then
! LIQUID

! use Goff-Gratch for T < 202.0 K and Flatau et al. equal or above 202.0 K
         if (t.ge.202.0) then
            dt = t-trplpt
            polysvp1 = a0 + dt*(a1+dt*(a2+dt*(a3+dt*(a4+dt*(a5+dt*(a6+dt*(a7+a8*dt)))))))
            polysvp1 = polysvp1*100.
         else
! note: uncertain below -70 C, but produces physical values (non-negative) unlike flatau
            polysvp1 = 10.**(-7.90298*(373.16/t-1.)+ &
                5.02808*alog10(373.16/t)- &
                1.3816e-7*(10**(11.344*(1.-t/373.16))-1.)+ &
                8.1328e-3*(10**(-3.49149*(373.16/t-1.))-1.)+ &
                alog10(1013.246))*100.
         end if

         endif


 end function polysvp1

!------------------------------------------------------------------------------------------!

 real function DERF(X)

 implicit none

 real :: X
 real, dimension(0 : 64) :: A, B
 real :: W,T,Y
 integer :: K,I
      data A/                                                 &
         0.00000000005958930743E0, -0.00000000113739022964E0, &
         0.00000001466005199839E0, -0.00000016350354461960E0, &
         0.00000164610044809620E0, -0.00001492559551950604E0, &
         0.00012055331122299265E0, -0.00085483269811296660E0, &
         0.00522397762482322257E0, -0.02686617064507733420E0, &
         0.11283791670954881569E0, -0.37612638903183748117E0, &
         1.12837916709551257377E0,                            &
         0.00000000002372510631E0, -0.00000000045493253732E0, &
         0.00000000590362766598E0, -0.00000006642090827576E0, &
         0.00000067595634268133E0, -0.00000621188515924000E0, &
         0.00005103883009709690E0, -0.00037015410692956173E0, &
         0.00233307631218880978E0, -0.01254988477182192210E0, &
         0.05657061146827041994E0, -0.21379664776456006580E0, &
         0.84270079294971486929E0,                            &
         0.00000000000949905026E0, -0.00000000018310229805E0, &
         0.00000000239463074000E0, -0.00000002721444369609E0, &
         0.00000028045522331686E0, -0.00000261830022482897E0, &
         0.00002195455056768781E0, -0.00016358986921372656E0, &
         0.00107052153564110318E0, -0.00608284718113590151E0, &
         0.02986978465246258244E0, -0.13055593046562267625E0, &
         0.67493323603965504676E0,                            &
         0.00000000000382722073E0, -0.00000000007421598602E0, &
         0.00000000097930574080E0, -0.00000001126008898854E0, &
         0.00000011775134830784E0, -0.00000111992758382650E0, &
         0.00000962023443095201E0, -0.00007404402135070773E0, &
         0.00050689993654144881E0, -0.00307553051439272889E0, &
         0.01668977892553165586E0, -0.08548534594781312114E0, &
         0.56909076642393639985E0,                            &
         0.00000000000155296588E0, -0.00000000003032205868E0, &
         0.00000000040424830707E0, -0.00000000471135111493E0, &
         0.00000005011915876293E0, -0.00000048722516178974E0, &
         0.00000430683284629395E0, -0.00003445026145385764E0, &
         0.00024879276133931664E0, -0.00162940941748079288E0, &
         0.00988786373932350462E0, -0.05962426839442303805E0, &
         0.49766113250947636708E0 /
      data (B(I), I = 0, 12) /                                 &
         -0.00000000029734388465E0,  0.00000000269776334046E0, &
         -0.00000000640788827665E0, -0.00000001667820132100E0, &
         -0.00000021854388148686E0,  0.00000266246030457984E0, &
          0.00001612722157047886E0, -0.00025616361025506629E0, &
          0.00015380842432375365E0,  0.00815533022524927908E0, &
         -0.01402283663896319337E0, -0.19746892495383021487E0, &
          0.71511720328842845913E0 /
      data (B(I), I = 13, 25) /                                &
         -0.00000000001951073787E0, -0.00000000032302692214E0, &
          0.00000000522461866919E0,  0.00000000342940918551E0, &
         -0.00000035772874310272E0,  0.00000019999935792654E0, &
          0.00002687044575042908E0, -0.00011843240273775776E0, &
         -0.00080991728956032271E0,  0.00661062970502241174E0, &
          0.00909530922354827295E0, -0.20160072778491013140E0, &
          0.51169696718727644908E0 /
      data (B(I), I = 26, 38) /                                &
         0.00000000003147682272E0, -0.00000000048465972408E0,  &
         0.00000000063675740242E0,  0.00000003377623323271E0,  &
        -0.00000015451139637086E0, -0.00000203340624738438E0,  &
         0.00001947204525295057E0,  0.00002854147231653228E0,  &
        -0.00101565063152200272E0,  0.00271187003520095655E0,  &
         0.02328095035422810727E0, -0.16725021123116877197E0,  &
         0.32490054966649436974E0 /
      data (B(I), I = 39, 51) /                                &
         0.00000000002319363370E0, -0.00000000006303206648E0,  &
        -0.00000000264888267434E0,  0.00000002050708040581E0,  &
         0.00000011371857327578E0, -0.00000211211337219663E0,  &
         0.00000368797328322935E0,  0.00009823686253424796E0,  &
        -0.00065860243990455368E0, -0.00075285814895230877E0,  &
         0.02585434424202960464E0, -0.11637092784486193258E0,  &
         0.18267336775296612024E0 /
      data (B(I), I = 52, 64) /                                &
        -0.00000000000367789363E0,  0.00000000020876046746E0,  &
        -0.00000000193319027226E0, -0.00000000435953392472E0,  &
         0.00000018006992266137E0, -0.00000078441223763969E0,  &
        -0.00000675407647949153E0,  0.00008428418334440096E0,  &
        -0.00017604388937031815E0, -0.00239729611435071610E0,  &
         0.02064129023876022970E0, -0.06905562880005864105E0,  &
         0.09084526782065478489E0 /
      W = ABS(X)
      if (W .LT. 2.2D0) then
          T = W * W
          K = INT(T)
          T = T - K
          K = K * 13
          Y = ((((((((((((A(K) * T + A(K + 1)) * T +              &
              A(K + 2)) * T + A(K + 3)) * T + A(K + 4)) * T +     &
              A(K + 5)) * T + A(K + 6)) * T + A(K + 7)) * T +     &
              A(K + 8)) * T + A(K + 9)) * T + A(K + 10)) * T +    &
              A(K + 11)) * T + A(K + 12)) * W
      elseif (W .LT. 6.9D0) then
          K = INT(W)
          T = W - K
          K = 13 * (K - 2)
          Y = (((((((((((B(K) * T + B(K + 1)) * T +               &
              B(K + 2)) * T + B(K + 3)) * T + B(K + 4)) * T +     &
              B(K + 5)) * T + B(K + 6)) * T + B(K + 7)) * T +     &
              B(K + 8)) * T + B(K + 9)) * T + B(K + 10)) * T +    &
              B(K + 11)) * T + B(K + 12)
          Y = Y * Y
          Y = Y * Y
          Y = Y * Y
          Y = 1 - Y * Y
      else
          Y = 1
      endif
      if (X .LT. 0) Y = -Y
      DERF = Y

 end function DERF


 real function qv_sat(t_atm,p_atm,ind_wrt)

!------------------------------------------------------------------------------------
! Calls polysvp1 to obtain the saturation vapor pressure, and then computes
! and returns the saturation mixing ratio, with respect to either liquid or ice,
! depending on value of 'ind_wrt'
!------------------------------------------------------------------------------------

 implicit none

 !Calling parameters:
 real    :: t_atm    !temperature [K]
 real    :: p_atm    !pressure    [Pa]
 integer :: ind_wrt  !index, 0 = w.r.t. liquid, 1 = w.r.t. ice

 !Local variables:
 real    :: e_pres ,ep_2       !saturation vapor pressure [Pa]

 !------------------
  ep_2 = 0.622

  e_pres = polysvp1(t_atm,ind_wrt)
  qv_sat = ep_2*e_pres/max(1.e-3,(p_atm-e_pres))


 return
 end function qv_sat

 END MODULE microphy_qv_sat

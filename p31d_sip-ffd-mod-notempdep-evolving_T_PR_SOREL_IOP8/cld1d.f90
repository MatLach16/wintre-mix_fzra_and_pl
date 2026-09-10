subroutine columnmodel

!----------------------------------------------------------------------------!
!  1-D kinematic cloud model driver for testing microphysics scheme.
!  Reads in temperature and humidity from a sounding, prescribes an updraft
!  profile, and begins a time integration for a simple cloud simulation.
!
!  The following procedures are performed at each time step:
!
!  - update vertical velocity
!  - compute column-integrated mass (and number) for the beginning of step
!  - update model variables due to advection, compressibility, and divergence
!  - force global concervation of mass (and number)
!  - compute total integrated mass before call to microphysics
!  - call microphysics
!  - compute total integrated mass after call to microphysics
!  - add low-level moisture (to prevent depletion)
!  - write output files
!
!--------------------------------------------------------------------------!
! Variable names:
! -----   ----
! cld1d    p3
! -----   ----
!  QI     qitot  - total (deposition + rime) ice mass mixing ratio
!  QG     qrim   - rime ice mass mixing ratio
!  QL     qiliq_in  - liquid on ice mass mixing ratio (optional)
!  NI     nitot  - ice number mixing ratio
!  BG     birim  - rime volume mixing ratio
!  ZI     zitot  - 6th moment mixing ratio
!
!--------------------------------------------------------------------------!
!  Author:         Jason Milbrandt
!  Last modified:  2023-05-28
!--------------------------------------------------------------------------!

      use subs_cld1d
      use microphy_p3 
      !use tracers_info_mod
      use input_kinematic_model
      !use cond11 

      implicit none
      !Other parameters
      ! logical, parameter :: liqFrac      = .false.
      ! logical, parameter :: trplMomIce   = .false.

      logical, parameter :: scpf_on      = .false.  ! switch for cloud fraction parameterization (SCPF)
      real,    parameter :: scpf_pfrac   = 1.        ! precipitation fraction factor (SCPF)
      real,    parameter :: scpf_resfact = 1.        ! model resolution factor (SCPF)
      logical, parameter :: prog_nc_ssat = .true.
     !logical, parameter :: nk_BOTTOM    = .true.   !.T. --> nk_input at bottom
      logical, parameter :: debug_on     = .false.   ! switch for run-time check-values in p3_main
      real, parameter    :: clbfact_dep = 1.0       !calibration factor for deposition
      real, parameter    :: clbfact_sub = 1.0       !calibration factor for sublimation

      integer      :: airmass,SCHEME,iice !ttotmin
      real         :: Htop0,tscale1,tscale2
      logical      :: microON,EVOLVING,TDMIN
      character*10 :: sndcase

      parameter (sndcase = 'ALBERTA')
      parameter (microON = .true. )     ! call microphysics SCHEME
      parameter (TDMIN   = .true. )     ! prevent low-level moisture depletion
      parameter (EVOLVING= .false. )     ! switch for evolving updraft
      !parameter (AMPA    = 2.     )     ! initial central updraft speed [m s-1] (evolving only) Now input ML
      !parameter (AMPB    = 5.     )     ! maximum central updraft speed [m s-1] Now input ML
      parameter (Htop0   = 4000.  )     ! initial height of cloud top   [m]
      parameter (tscale1 = 5400.  )     ! period for evolving AMPB      [s]
      parameter (tscale2 = 5400.  )     ! period for evolving Hcld      [s]
      ! parameter (nk_input      =  41    )     ! number of vertical levels
!     parameter (nk_input      =  62    )     ! number of vertical levels
!     parameter (nk_input      =  86    )     ! number of vertical levels
      ! parameter (outfreq =  max(1,int(dt_p3/60))     )     ! output every 'OUTFREQ' minutes
      !parameter (dt_p3      = 10.    )     ! time step                     [s]
      ! parameter (ttotmin = 90     )     ! total integration time	[min] | input
      ! parameter (lev_cond11 = 42   )    !nk_input+1    



      character(len=16), parameter :: model = 'KIN1D'
!     character(len=16), parameter :: model = 'WRF'  !for level tests
      logical, parameter           :: abort_on_err = .false.
      logical, parameter           :: dowr = .true.

      character(len=1024), parameter :: LT_path  = '/chinook/lachapelle/lookup_tables'
!     character(len=1024), parameter :: LT_path = '/users/milbrand/mp_p3/lookupTables/tables'  ! override default


!---------------------------------------------------------------------------------!
!#include "consphy.cdk"  (necessary parameters only)
      real, parameter :: TRPL     =.27316e+3          !K; triple point of water
      real, parameter :: EPS1     =.62194800221014    ! ; RGASD/RGASV
      real, parameter :: EPS2     =.3780199778986     !; 1 - EPS1
      real, parameter :: PI       =.314159265359e+1   ! PI constant = ACOS(-1)
      real, parameter :: GRAV     =.980616e+1         ! M s-2; gravitational acceleration
      real, parameter :: rhow     = 1000.
      
!#include "dintern.cdk"  (necessary variables only)
      real   :: TTT,PRS,QQQ
      real*8 :: FOEW,FOQST

!------------------------------------------------------------------------------!
!#include "fintern.cdk"
!   DEFINITION DES FONCTIONS THERMODYNAMIQUES DE BASE
!   POUR LES CONSTANTES, UTILISER LE COMMON /CONSPHY/
!     NOTE: TOUTES LES FONCTIONS TRAVAILLENT AVEC LES UNITES S.I.
!     FONCTION DE TENSION DE VAPEUR SATURANTE (TETENS) - EW OU EI SELON TT
      FOEW(TTT) = 610.78D0*DEXP( DMIN1(DSIGN(17.269D0,                     &
       DBLE(TTT)-DBLE(TRPL)),DSIGN                                         &
       (21.875D0,DBLE(TTT)-DBLE(TRPL)))*DABS(DBLE(TTT)-DBLE(TRPL))/        &
       (DBLE(TTT)-35.86D0+DMAX1(0.D0,DSIGN                                 &
       (28.2D0,DBLE(TRPL)-DBLE(TTT)))))

!     FONCTION CALCULANT L'HUMIDITE SPECIFIQUE SATURANTE (QSAT)
      FOQST(TTT,PRS) = DBLE(EPS1)/(DMAX1(1.D0,DBLE(PRS)/FOEW(TTT))-        &
       DBLE(EPS2))


!------------------------------------------------------------------------------!


      integer :: i,j,k,k1,k2,nkcld,n,step,lv,nlvs,tmin
      integer :: its,ite,kts,kte

      real PMv0,PMc0,PMr0,PMi0,PMg0,PMs0,PMh0,PMv1,PMc1,PMr1,PMi1,PMg1,PMs1,PMh1,     &
           PMv3,PMc3,PMr3,PMi3,PMg3,PMs3,PMh3,rhoQmax,depth,gam1,gam3,lamda,No,cnt,   &
           dm,alpha,cxh,cmx,Nfact1,Dmx(nk_input),ref(nk_input),Zplt,H,H0,Hcld,ZePlt,HcldTop,      &
           wmaxH,c1,c2,esat,alfa,alfa0,walfa,Kdiff,dz,dzsq,INTthr,dum,HcldBase,       &
           LAMr,ZZ,AMPL,AMPL0,Dc,Dr,Di,Ds,Dg,Dh,cmr,cmi,cms,cmg,cmh,thrd,sig,Ltot,    &
           M1,M2,Aii,BASE,RATIOc,RATIOv,RATIOr,RATIOi,RATIOg,RATIOs,RATIOh,tsec,      &
           tminr,af,bf,cm,Cx,Nox,rhox,PR

      real, parameter    :: eps    = 0.622
      real, parameter    :: g_cld  = 9.81
      real, parameter    :: Rd_cld = 287.0
      real, parameter    :: T0     = 273.15
      real, parameter    :: cp_cld = 1005.

    ! ******** CHO ************
    ! Add constant for dew point
      real, parameter :: constA1	= .253e+9
      real, parameter :: constB1	= .5420e+4
      real, parameter :: constA		= .341e+10
      real, parameter :: constB		= .6130e+4

      integer, parameter :: ni     = 1
      integer, parameter :: lnmax  = 4000
      integer, parameter :: nt     = nint(ttotmin*60/dt_p3)

      real, dimension(nk_input)    :: z,p,tt,rh,td,w1,DIV,Tdry,zcld,w1cld,rho,alfa2
      real, dimension(nk_input)    :: z_1,p_1,tt_1,td_1,qv_1,rh_1
      real, dimension(nk_input)    :: z_2,p_2,tt_2,td_2,qv_2,rh_2
      real, dimension(nk_input)    :: z_3,p_3,tt_3,td_3,qv_3,rh_3
      real, dimension(nk_input)    :: z_4,p_4,tt_4,td_4,qv_4,rh_4
      real, dimension(nk_input)    :: z_5,p_5,tt_5,td_5,qv_5,rh_5
      real, dimension(nk_input)    :: z_6,p_6,tt_6,td_6,qv_6,rh_6


      real, dimension(ni,nk_input) :: tt0,tt1,Qv0,Qv1,w,SIGMA,Qsat,womega,th2d0,th2d1,p2d,dz2d,td1,tf1, dTdt
      real, dimension(lnmax) :: Pin1,Zin1,TTin1,TDin1,QVin1,RHin1
      real, dimension(lnmax) :: Pin2,Zin2,TTin2,TDin2,QVin2,RHin2
      real, dimension(lnmax) :: Pin3,Zin3,TTin3,TDin3,QVin3,RHin3
      real, dimension(lnmax) :: Pin4,Zin4,TTin4,TDin4,QVin4,RHin4
      real, dimension(lnmax) :: Pin5,Zin5,TTin5,TDin5,QVin5,RHin5
      real, dimension(lnmax) :: Pin6,Zin6,TTin6,TDin6,QVin6,RHin6

      real :: time_sounding1,time_sounding2,time_sounding3,time_sounding4,time_sounding5,time_sounding6

    ! Prognostic hydrometeor variables:
      !liquid:
      real, dimension(ni,nk_input) :: Qc0,Qc1,Nc0,Nc1,Qr0,Qr1,Nr0,Nr1,ssat1  !ssat0
      !ice:
      real, dimension(ni,nk_input,n_iceCat) :: Qi0,Qi1,Qg0,Qg1,Ni0,Zi0,Ni1,Bg0,Bg1,Zi1,Ql1,Ql0

    ! Source-Sink term arrays:
      integer, parameter               :: n_diag_2d = 20
      integer, parameter               :: n_diag_3d = 20
      real, dimension(ni,n_diag_2d)    :: diag_2d     !user-defined 2D diagnostic arrays (for output)
      real, dimension(ni,nk_input,n_diag_3d) :: diag_3d     !user-defined 3D diagnostic arrays (for output)
      real, dimension(ni,nk_input,n_iceCat)  :: diag_reffi,diag_vmi,diag_di,diag_rhoi,diag_dhmax,diag_vni,diag_vzi
      real, dimension(ni,nk_input)           :: diag_reffc

    ! Precipitation rates:
      real, dimension(ni)      :: prt_liq       ! precip rate, total liquid     m s-1
      real, dimension(ni)      :: prt_sol       ! precip rate, total solid      m s-1
      real, dimension(ni)      :: prt_drzl      ! precip rate, drizzle          m s-1
      real, dimension(ni)      :: prt_rain      ! precip rate, rain             m s-1
      real, dimension(ni)      :: prt_crys      ! precip rate, ice cystals      m s-1
      real, dimension(ni)      :: prt_snow      ! precip rate, snow             m s-1
      real, dimension(ni)      :: prt_grpl      ! precip rate, graupel          m s-1
      real, dimension(ni)      :: prt_pell      ! precip rate, ice pellets      m s-1
      real, dimension(ni)      :: prt_hail      ! precip rate, hail             m s-1
      real, dimension(ni)      :: prt_sndp      ! precip rate, unmelted snow    m s-1
      real, dimension(ni)      :: prt_wsnow     ! precip rate, very wet snow    m s-1

    ! Diagnostics, etc.
      real, dimension(ni,nk_input)   :: diag_ZET,Qvinit,GZ,scpf_cldfrac
      real, dimension(nk_input)      :: COMP
      real, dimension(ni)      :: p0,p0m,LR,SR,diag_ZEC
      real                     :: delz,dum1,dum2,dum3,dum4,dum5,dum6,dum7,dum8
      real, dimension(ni,nk_input,6) :: qi_type

      integer                  :: nk_read,kk,stat
      logical, parameter       :: log_predictNc = .true.
      integer, dimension(200)  :: kskiplevs

      real                     :: time1,time2,time3,time4,time5 ! for timing tests

      real, dimension(nk_input)              :: qv
      real                     :: qitop_csv,nitop_csv,esw
!---------------------------------------------------------------------------------------!
    ! ****** output files *****
      open(unit=20, file='./outQcld.dat')
      open(unit=24, file='./outNcld.dat')
      open(unit=27, file='./sfcprec.dat')
    
    ! ****** output file headers*****
    write(27,'(4a10)') 'time', 'prt_liq', 'prt_sol', 'PR'
    
    if (n_iceCat .eq. 1) then
      write(20,'(17a10)') 'k','time','z','w','rho','p',&
      'T','Td','Qc','Qr','Qg','Qi','Qv',&
      'v_ice','Ql','tf1','diag_rhoi'
      write(24,'(9a10)') 'k','time','z','Nc','Nr',&
      'Ni','Bg','diag_di','diag_rhoi'

    elseif (n_iceCat .eq. 2) then
      write(20,'(21a10)') 'k','time','z','w','rho','p','T',&
      'Td','Qc','Qr','Qg','Qi','Qv','v_ice','Ql',&
      'tf1','Qg2','Qi2','Ql2','Qsedi_i1','Qsedi_i2'
      write(24,'(14a15)') 'k','time','z','Nc','Nr','Ni',&
      'Bg','diag_di','diag_rhoi','Ni2','Bg2','Nsedi_i1','Nsedi_i2',&
      'diag_rhoi2'

    elseif (n_iceCat .eq. 3) then
      write(20,'(23a10)') 'k','time','z','w','rho','p','T',&
      'Td','Qc','Qr','Qg','Qi','Qv','v_ice','Ql',&
      'tf1','diag_rhoi','Qg2','Qi2','Ql2','Qg3','Qi3','Ql3'
      write(24,'(13a10)') 'k','time','z','Nc','Nr','Ni',&
      'Bg','diag_di','diag_rhoi','Ni2','Bg2','Ni3','Bg3'

    elseif (n_iceCat .eq. 4) then
      write(20,'(26a10)') 'k','time','z','w','rho','p','T',&
      'Td','Qc','Qr','Qg','Qi','Qv','v_ice','Ql',&
      'tf1','diag_rhoi','Qg2','Qi2','Ql2','Qg3','Qi3','Ql3',&
      'Qg4','Qi4','Ql4'
      write(24,'(15a10)') 'k','time','z','Nc','Nr','Ni',&
      'Bg','diag_di','diag_rhoi','Ni2','Bg2','Ni3','Bg3',&
      'Ni4','Bg4'
!---------------------------------------------------------------------------------------!
    endif 
  write(8000,'(82a10)')  'time','k','qimlt','qrmlt','qrcol','qccol', &
        'nrshdr','nrcol','qisub','qidep', &
        'qrcon','qrevp','qwgrth','qinuc', &
        'qcheti', 'qchetc', 'qrhetc', 'qrheti', &
        'qrmul','qcshd','nimul','diam_ice1', &
        'qimlt2','qrmlt2','qrcol2','qccol2', &
        'nrshdr2','nrcol2','qisub2','qidep2', &
        'qwgrth2','qinuc2', &
        'qcheti2', 'qchetc2', 'qrhetc2', 'qrheti2', &
        'qrmul2','qcshd2','nimul2','qicol1_2','diam_ice2', &
        'qimlt3','qrmlt3','qrcol3','qccol3', &
        'nrshdr3','nrcol3','qisub3','qidep3', &
        'qwgrth3','qinuc3', &
        'qcheti3', 'qchetc3', 'qrhetc3', 'qrheti3', &
        'qrmul3','qcshd3','nimul3','qicol1_3','qicol2_3','diam_ice3', &
        'qimlt4','qrmlt4','qrcol4','qccol4', &
        'nrshdr4','nrcol4','qisub4','qidep4', &
        'qwgrth4','qinuc4', &
        'qcheti4', 'qchetc4', 'qrhetc4', 'qrheti4', &
        'qrmul4','qcshd4','nimul4','qicol1_4','qicol2_4','qicol3_4','diam_ice4'

  write(9000,'(50a10)')  'time','k','ninuc','nimlt','nisub',&
        'nislf','nrhetc','nrheti','nchetc','ncheti','nimul','nlevp','nrhetic',&
        'ninuc2','nimlt2','nisub2',&
        'nislf2','nrhetc2','nrheti2','nchetc2','ncheti2','nimul2','nlevp2',&
        'nicol1_2','nrhetic2',&
        'ninuc3','nimlt3','nisub3',&
        'nislf3','nrhetc3','nrheti3','nchetc3','ncheti3','nimul3','nlevp3',&
        'nicol1_3','nicol2_3',&
        'ninuc4','nimlt4','nisub4',&
        'nislf4','nrhetc4','nrheti4','nchetc4','ncheti4','nimul4','nlevp4',&
        'nicol1_4','nicol2_4','nicol3_4'

      print*
      print*, '** Remember to use compiler debug options for testing new code (modfify Makefile appropriately) **'
      print*

      diag_dhmax = -1.

      its = 1
      ite = 1
      kts = 1
      kte = nk_input

      ssat1 = 0.  !for P3_v2.3.3, prognostic supersaturation is not available

      PR = 0. ! accumulated precipitation
      time3 = 0. ! for total timing test

!------- INITIALIZE w, tt, td, p AND Q[x]  PROFILES: -------

      H= 4000.

      ! Sounging #1
      time_sounding1 = 0.  ! Time in seconds 
      open(unit=12, file='/home/rml000/sitestore8/UQAM/projet_margaux/Pour_Mathieu_fig4/' &
      //'create_figure_papier/simu1d/src_sip-ffd-mod-notempdep-evolving_T_PR_SOREL_IOP8/temp_to_use_0.dat')
      read(12,*) nlvs
      read(12,*)
      read(12,*)
      read(12,*)
      read(12,*)
      read(12,*)
      do lv=1,nlvs
         read(12,*) Pin1(lv),Zin1(lv),TTin1(lv),QVin1(lv),RHin1(lv)   !new
         TTin1(lv) = TTin1(lv) + T0		        !convert from C to Kelvins
         Pin1(lv)  = Pin1(lv)*100                          !convert from mb to Pa
      enddo
      close(12)

    ! Sounging #2
      time_sounding2 = 2000. + 3600.*2.  ! Time in seconds
      open(unit=13, file='/home/rml000/sitestore8/UQAM/projet_margaux/Pour_Mathieu_fig4/' &
      //'create_figure_papier/simu1d/src_sip-ffd-mod-notempdep-evolving_T_PR_SOREL_IOP8/temp_to_use_1.dat')
      read(13,*) nlvs
      read(13,*)
      read(13,*)
      read(13,*)
      read(13,*)
      read(13,*)
      do lv=1,nlvs
         read(13,*) Pin2(lv),Zin2(lv),TTin2(lv),QVin2(lv),RHin2(lv)   !new
         TTin2(lv) = TTin2(lv) + T0		        !convert from C to Kelvins
         Pin2(lv)  = Pin2(lv)*100                          !convert from mb to Pa
      enddo
      close(13)

    ! Sounging #3
      time_sounding3 = 2000. + 3600.*4.  ! Time in seconds
      open(unit=14, file='/home/rml000/sitestore8/UQAM/projet_margaux/Pour_Mathieu_fig4/' &
      //'create_figure_papier/simu1d/src_sip-ffd-mod-notempdep-evolving_T_PR_SOREL_IOP8/temp_to_use_2.dat')
      read(14,*) nlvs
      read(14,*)
      read(14,*)
      read(14,*)
      read(14,*)
      read(14,*)
      do lv=1,nlvs
         read(14,*) Pin3(lv),Zin3(lv),TTin3(lv),QVin3(lv),RHin3(lv)   !new
         TTin3(lv) = TTin3(lv) + T0		        !convert from C to Kelvins
         Pin3(lv)  = Pin3(lv)*100                          !convert from mb to Pa
      enddo
      close(14)

    ! Sounging #4
      time_sounding4 = 2000. + 3600.*6.  ! Time in seconds
      open(unit=15, file='/home/rml000/sitestore8/UQAM/projet_margaux/Pour_Mathieu_fig4/' &
      //'create_figure_papier/simu1d/src_sip-ffd-mod-notempdep-evolving_T_PR_SOREL_IOP8/temp_to_use_3.dat')
      read(15,*) nlvs
      read(15,*)
      read(15,*)
      read(15,*)
      read(15,*)
      read(15,*)
      do lv=1,nlvs
         read(15,*) Pin4(lv),Zin4(lv),TTin4(lv),QVin4(lv),RHin4(lv)   !new
         TTin4(lv) = TTin4(lv) + T0		        !convert from C to Kelvins
         Pin4(lv)  = Pin4(lv)*100                          !convert from mb to Pa
      enddo
      close(15)

    ! Sounging #5
      time_sounding5 = 2000. + 3600.*8.  ! Time in seconds
      open(unit=16, file='/home/rml000/sitestore8/UQAM/projet_margaux/Pour_Mathieu_fig4/' &
      //'create_figure_papier/simu1d/src_sip-ffd-mod-notempdep-evolving_T_PR_SOREL_IOP8/temp_to_use_4.dat')
      read(16,*) nlvs
      read(16,*)
      read(16,*)
      read(16,*)
      read(16,*)
      read(16,*)
      do lv=1,nlvs
         read(16,*) Pin5(lv),Zin5(lv),TTin5(lv),QVin5(lv),RHin5(lv)   !new
         TTin5(lv) = TTin5(lv) + T0		        !convert from C to Kelvins
         Pin5(lv)  = Pin5(lv)*100                          !convert from mb to Pa
      enddo
      close(16)

    ! Sounging #6
      time_sounding6 = 2000. + 3600.*10. + 900.  ! Time in seconds
      open(unit=17, file='/home/rml000/sitestore8/UQAM/projet_margaux/Pour_Mathieu_fig4/' &
      //'create_figure_papier/simu1d/src_sip-ffd-mod-notempdep-evolving_T_PR_SOREL_IOP8/temp_to_use_5.dat')
      read(17,*) nlvs
      read(17,*)
      read(17,*)
      read(17,*)
      read(17,*)
      read(17,*)
      do lv=1,nlvs
         read(17,*) Pin6(lv),Zin6(lv),TTin6(lv),QVin6(lv),RHin6(lv)   !new
         TTin6(lv) = TTin6(lv) + T0		        !convert from C to Kelvins
         Pin6(lv)  = Pin6(lv)*100                          !convert from mb to Pa
      enddo
      close(17)


      ! Precipitation rate
      open(unit=18, file='/home/rml000/sitestore8/UQAM/projet_margaux/Pour_Mathieu_fig4/' &
      //'create_figure_papier/simu1d/src_sip-ffd-mod-notempdep-evolving_T_PR_SOREL_IOP8/df_manobs_iop8_sorel_1s.csv')
      read(18,*)

!  END OF INITIAL SET-UP
      call P3_INIT(LT_path,n_iceCat,trplMomIce,liqFrac,model,stat,abort_on_err,dowr,bimm=bimm)      !v4.0.0_b42

!  Set up vertical level (evenly-spaced z-levels)
!  Using only the first sounding for now (Zin1, Pin1, TTin1, QVin1)
     z(nk_input) = Zin1(1)
     delz = (H-Zin1(1))/(nk_input-1)
     do k=1,nk_input
        z(k) = Zin1(1) + delz*(nk_input-k)
        dz2d(1,k) = delz
     enddo
     dz = z(1)-z(2)
     dzsq = dz**2.
     H0= z(nk_input)
     GZ(1,:)= z*GRAV



      do k=1,nk_input
         call vertint2b( p_1(k),z(k),Pin1, Zin1,nk_input,nlvs)
         call vertint2b(tt_1(k),z(k),TTin1,Zin1,nk_input,nlvs)
         call vertint2b(qv_1(k),z(k),QVin1,Zin1,nk_input,nlvs)
         call vertint2b(rh_1(k),z(k),RHin1,Zin1,nk_input,nlvs)

         call vertint2b( p_2(k),z(k),Pin2, Zin2,nk_input,nlvs)
         call vertint2b(tt_2(k),z(k),TTin2,Zin2,nk_input,nlvs)
         call vertint2b(qv_2(k),z(k),QVin2,Zin2,nk_input,nlvs)
         call vertint2b(rh_2(k),z(k),RHin2,Zin2,nk_input,nlvs)
         
         call vertint2b( p_3(k),z(k),Pin3, Zin3,nk_input,nlvs)
         call vertint2b(tt_3(k),z(k),TTin3,Zin3,nk_input,nlvs)
         call vertint2b(rh_3(k),z(k),RHin3,Zin3,nk_input,nlvs)         
         
         call vertint2b( p_4(k),z(k),Pin4, Zin4,nk_input,nlvs)
         call vertint2b(tt_4(k),z(k),TTin4,Zin4,nk_input,nlvs)
         call vertint2b(rh_4(k),z(k),RHin4,Zin4,nk_input,nlvs)         
      
         call vertint2b( p_5(k),z(k),Pin5, Zin5,nk_input,nlvs)
         call vertint2b(tt_5(k),z(k),TTin5,Zin5,nk_input,nlvs)
         call vertint2b(rh_5(k),z(k),RHin5,Zin5,nk_input,nlvs)
      
         call vertint2b( p_6(k),z(k),Pin6, Zin6,nk_input,nlvs)
         call vertint2b(tt_6(k),z(k),TTin6,Zin6,nk_input,nlvs)
         call vertint2b(rh_6(k),z(k),RHin6,Zin6,nk_input,nlvs)

      enddo


!  Interpolate p,tt,td from sounding data and initialize Q[x]:
      do k=1,nk_input
      
         p(k) = p_1(k)
         tt0(1,k) = tt_1(k)
         Qsat(1,k)= FOQST(tt(k),p(k))
      !    if(sndcase=='ALBERTA') Qv0(1,k) = FOQST(td(k),p(k))  !Alberta (Td read in)
         Qsat(1,k)= qv_sat(tt_1(k),p_1(k),0)
         Qv0(1,k) = qv_1(k)! Qsat(1,k) !qv_1(k)!qv_sat(tt1(k),p(k),0)
         Qvinit(1,k) = Qv0(1,k)

         Qc0(1,k) = 0.
         Nc0(1,k) = 0.
         Qr0(1,k) = 0.
         Nr0(1,k) = 0.
!        ssat0(1,k) = 0.   !supersaturation mixing ratio (qv-qv_sat)

         Qi0(1,k,:) = 0.   !qitot (total ice mass mixing ratio)
         Qg0(1,k,:) = 0.   !qirim (riming mass mixing ratio)
         Ql0(1,k,:) = 0.   !qiliq (liquid on ice mass mixing ratio)
         Ni0(1,k,:) = 0.   !nitot (number mixing ratio)
         Bg0(1,k,:) = 0.   !birim (rime volume mixing ratio)
         Zi0(1,k,:) = 0.   !zitot (reflectivity mixing ratio)



         Dmx(k) = 0.

        !  print*, p(k),z(k),tt(k),td(k)
      enddo


!**  z, w, p, tt, td and Qx0 are now in arrays from 1 to nk_input  **

      do k = 1,nk_input
         Qsat(1,k)= 0.
         Qv1(1,k)= Qv0(1,k)
         Qc1(1,k)= Qc0(1,k)
         Nc1(1,k)= Nc0(1,k)
         Qr1(1,k)= Qr0(1,k)
         Nr1(1,k)= Nr0(1,k)
!        ssat1(1,k) = ssat0(1,k)

         Qi1(1,k,:)= Qi0(1,k,:)
         Qg1(1,k,:)= Qg0(1,k,:)
         Ql1(1,k,:)= Ql0(1,k,:)
         Ni1(1,k,:)= Ni0(1,k,:)
         Bg1(1,k,:)= Bg0(1,k,:)
         Zi1(1,k,:)= Zi0(1,k,:)
      enddo

!-------------------------------------------------------------------------!

!=========================================================================!
!  MAIN TIME LOOP:
      print*, 'Starting main loop... for ',nt, 'timesteps'

      do step = 1,nt
        tsec = step*dt_p3         !integration time [s]    (real)
        tminr= tsec/60.        !integration time [min]  (real)
        tmin = nint(tminr)     !integration time [min]  (integer)

            ! Interpolation on profiles
            if (tsec .ge. time_sounding1 .and. tsec .lt. time_sounding2) then
            do k = 1,nk_input
                  tt(k) = tt_1(k) + (tt_2(k)-tt_1(k))*(tsec-time_sounding1)/(time_sounding2-time_sounding1)
                  rh(k) = rh_1(k) + (rh_2(k)-rh_1(k))*(tsec-time_sounding1)/(time_sounding2-time_sounding1)
                  ! Qv1(1,k) = qv_1(k) + (qv_2(k)-qv_1(k))*(tsec-time_sounding1)/(time_sounding2-time_sounding1)
            enddo
            elseif (tsec .ge. time_sounding2 .and. tsec .lt. time_sounding3) then
            do k = 1,nk_input
                  tt(k) = tt_2(k) + (tt_3(k)-tt_2(k))*(tsec-time_sounding2)/(time_sounding3-time_sounding2)
                  rh(k) = rh_2(k) + (rh_3(k)-rh_2(k))*(tsec-time_sounding2)/(time_sounding3-time_sounding2)
                  ! Qv1(1,k) = qv_2(k) + (qv_3(k)-qv_2(k))*(tsec-time_sounding2)/(time_sounding3-time_sounding2)
            enddo
            elseif (tsec .ge. time_sounding3 .and. tsec .lt. time_sounding4) then
            do k = 1,nk_input
                  tt(k) = tt_3(k) + (tt_4(k)-tt_3(k))*(tsec-time_sounding3)/(time_sounding4-time_sounding3)
                  rh(k) = rh_3(k) + (rh_4(k)-rh_3(k))*(tsec-time_sounding3)/(time_sounding4-time_sounding3)
                  ! Qv1(1,k) = qv_3(k) + (qv_4(k)-qv_3(k))*(tsec-time_sounding3)/(time_sounding4-time_sounding3)
            enddo
            elseif (tsec .ge. time_sounding4 .and. tsec .lt. time_sounding5) then
            do k = 1,nk_input
                  tt(k) = tt_4(k) + (tt_5(k)-tt_4(k))*(tsec-time_sounding4)/(time_sounding5-time_sounding4)
                  rh(k) = rh_4(k) + (rh_5(k)-rh_4(k))*(tsec-time_sounding4)/(time_sounding5-time_sounding4)
                  ! Qv1(1,k) = qv_4(k) + (qv_5(k)-qv_4(k))*(tsec-time_sounding4)/(time_sounding5-time_sounding4)
            enddo
            elseif (tsec .ge. time_sounding5 .and. tsec .lt. time_sounding6) then
            do k = 1,nk_input
                  tt(k) = tt_5(k) + (tt_6(k)-tt_5(k))*(tsec-time_sounding5)/(time_sounding6-time_sounding5)
                  rh(k) = rh_5(k) + (rh_6(k)-rh_5(k))*(tsec-time_sounding5)/(time_sounding6-time_sounding5)
                  ! Qv1(1,k) = qv_5(k) + (qv_6(k)-qv_5(k))*(tsec-time_sounding5)/(time_sounding6-time_sounding5)
            enddo

            end if

            do k = 1,nk_input
                  ! tt0(1,k) = tt(k)
                  tt1(1,k) = tt(k)
                  rho(k)= p(k)/(Rd_cld*tt(k))
                  p2d(1,k) = p(k)
                  ! print*, k,tt(k),tt0(1,k),tt1(1,k), p(k)
                  esw      = polysvp1(tt(k),0)
                  Qv1(1,k) = 0.622*rh(k)/100.*esw/p(k)

                  ! Qv1(1,k) = FOQST(tt(k),p(k))
                  ! Qv1(1,k) = (constA1*0.622*1000.0)/p(k)*dble(exp(-constB1/(tt(k))))
                  ! print*, k,tt(k),esw, rh(k), p(k), Qv1(1,k)
                  
            enddo
            ! stop "message"
            read(18,*) qitop_csv,nitop_csv
            
      ! Initialize de la neige au top (qi0)
            do k = 2,2
            Qi1(1,k,1) = qitop_csv !+ Qi1(1,k,1)
            Qg1(1,k,1) = 0. !+ Qg1(1,k,1)
            Ql1(1,k,1) = 0.
            !Qg1(1,k,:) = 0.00015
            !Qg1(1,k,:) = 0.00017
            Bg1(1,k,1) = Qg1(1,k,1)/900.
            Ni1(1,k,1) = nitop_csv !+ Ni1(1,k,1)
            enddo



!------------------------------------------------------------------------!
        if (microON) then

          th2d0(1,:) = tt0(1,:)*(1.e+5/p(:))**0.286
          th2d1(1,:) = tt1(1,:)*(1.e+5/p(:))**0.286

          call cpu_time(time1)

          if (.not. trplMomIce) Zi1(:,:,:) = 0.
          if (.not. liqFrac)    Ql1(:,:,:) = 0.

            ! print*, 'timestep: ', step
            ! print*, Qi1(1,2,1),Ni1(1,2,1)
            CALL P3_MAIN(Qc1,Nc1,Qr1,Nr1,th2d0,th2d1,Qv0,Qv1,dt_p3,Qi1,Qg1,Ql1,Ni1,Bg1,Zi1,         &
                              ssat1,w,p2d,dz2d,step,prt_liq,prt_sol,its,ite,kts,kte,n_iceCat,   &
                              diag_ZET,diag_reffc,diag_reffi,diag_vmi,diag_di,diag_rhoi,        &
                              n_diag_2d,diag_2d,n_diag_3d,diag_3d,log_predictNc,                &
                              trim(model),clbfact_dep,clbfact_sub,debug_on,                     &
                              scpf_on,scpf_pfrac,scpf_resfact,scpf_cldfrac,trplMomIce,liqFrac,  &
                              prt_drzl = prt_drzl,  &
                              prt_rain = prt_rain,  &
                              prt_crys = prt_crys,  &
                              prt_snow = prt_snow,  &
                              prt_grpl = prt_grpl,  &
                              prt_pell = prt_pell,  &
                              prt_hail = prt_hail,  &
                              prt_sndp = prt_sndp,  &
                              prt_wsnow = prt_wsnow, &
                              qi_type  = qi_type,    &
                              diag_dhmax = diag_dhmax, &
                              outfreq    = outfreq)
            ! print*, Qi1(1,2,1),Ni1(1,2,1)
            
          tt1(1,:) = th2d1(1,:)*(p(:)*1.e-5)**0.286

          !convert precipitation rates to units mm h-1 (from m s-1)
          prt_liq  = prt_liq *3.6e+6  !total liquid
          prt_sol  = prt_sol *3.6e+6  !total solid
          prt_drzl = prt_drzl*3.6e+6
          prt_rain = prt_rain*3.6e+6
          prt_crys = prt_crys*3.6e+6
          prt_snow = prt_snow*3.6e+6
          prt_grpl = prt_grpl*3.6e+6
          prt_pell = prt_pell*3.6e+6
          prt_hail = prt_hail*3.6e+6
          prt_sndp = prt_sndp*3.6e+6
          prt_wsnow = prt_wsnow*3.6e+6

          call cpu_time(time2)




!*-----------------------------------------------------------------------*!

      endif   !(microON)

!  Rearrange arrays:
         do k = 1,nk_input
            Qv0(1,k) = Qv1(1,k)
            tt0(1,k) = tt1(1,k)
            Qc0(1,k) = Qc1(1,k)
            Qr0(1,k) = Qr1(1,k)
            Nc0(1,k) = Nc1(1,k)
            Nr0(1,k) = Nr1(1,k)
            Qi0(1,k,:) = Qi1(1,k,:)
            Qg0(1,k,:) = Qg1(1,k,:)
            Ql0(1,k,:) = Ql1(1,k,:)
            Ni0(1,k,:) = Ni1(1,k,:)
            Bg0(1,k,:) = Bg1(1,k,:)
            Zi0(1,k,:) = Zi1(1,k,:)
!           ssat0(1,k) = ssat1(1,k)
         enddo
         p0m(1) = p0(1)

!--  Prevent low-level moisture depletion:     ----------------!
!    (Impose Td >= Td_initial after  xx min. [..< xx)] )
!        if (TDMIN .and. step<70) then
!          if (TDMIN) then
!            do k = 1,nk_input    ! or, e.g., k=nk_input-10,nk_input
!              if (z(k)<1000.) then
! !              Qv0(1,k)= amax1(Qv1(1,k),Qvinit(1,k))
!                Qv0(1,k)= amax1(Qv1(1,k),Qvinit(1,k)*0.4) !test
!                Qv1(1,k)= Qv0(1,k)
!              endif
!            enddo
!          endif
!----------------------------------------------------------------!
	do k = 1,nk_input
	     td1(1,k) = constB1/log((constA1*0.622*1000.0)/(Qv1(1,k)*p(k)))          ! calcualtion of the dew point
           tf1(1,k) = constB/log((constA*0.622*1000.0)/(Qv1(1,k)*p(k)))          ! calculation of the frost point
  !          if (td1(1,k) > tt0(1,k)) then    !to adjest the dew point temperature because the equation used is approximation
  !            td1(1,k) = tt0(1,k)
  !          endif
	enddo

       
         PR = PR + (prt_liq(1) + prt_sol(1))*dt_p3*1.e-3  !accumulated precipitation, mm
         time3 =time3+(time2-time1)

!  Output to files:

!    Precipitation rates at lowest level (surface):


        ! if (mod(tminr,float(1)) < 1.e-5) print*, 'time (min): ',tminr,(prt_liq(1) + prt_sol(1)),PR, time3
        ! if (mod(tminr,float(5)) < 1.e-5) print*, 'time (min): ',tminr,prt_sol(1),PR

        ! if (mod(tminr,float(1)) < 1.e-5) print*, Qg1(1,50,1)*100/(Qi1(1,50,1)+Qr1(1,50)), &
	!			Qr1(1,50)*100/(Qi1(1,50,1)+Qr1(1,50)), 			&
	!			(Qi1(1,50,1)-Qg1(1,50,1))*100/(Qi1(1,50,1)+Qr1(1,50))

    if (mod(tminr,float(outfreq)) < 1.e-5) then  !output only every OUTFREQ min   !for TESING

!    Precipitation rates at lowest level (surface):
    write(27,'(1f7.0,3f15.5)') tsec, prt_liq(1), prt_sol(1), PR

    do k=1,nk_input

!    Hydrometeors mass content profiles:
      if (n_iceCat .eq.1) then      
        write(20,'(1i5,2f7.0,14f17.5)') k,tsec,z(k),w(1,k),rho(k),p2d(1,k),         &
        tt1(1,k)-T0,td1(1,k)-T0,Qc1(1,k)*10**8,       &
        Qr1(1,k)*10**8,Qg1(1,k,1)*10**8,            &
        Qi1(1,k,1)*10**8,Qv1(1,k)*10**8,diag_vmi(1,k,1),Ql1(1,k,1)*10**8,     &
        tf1(1,k)-T0,diag_rhoi(1,k,1)

        write(24,'(1i5,2f7.0,6f25.5)') k, tsec,z(k),           &
        Nc1(1,k), Nr1(1,k), Ni1(1,k,1),          &
        Bg1(1,k,1)*10**8,diag_di(1,k,1),diag_rhoi(1,k,1)

      elseif (n_iceCat .eq.2) then
        write(20,'(1i5,2f7.0,18E20.5)') k,tsec,z(k),w(1,k),rho(k),p2d(1,k),         &
        tt1(1,k)-T0,td1(1,k)-T0,Qc1(1,k)*10**8,       &
        Qr1(1,k)*10**8,Qg1(1,k,1)*10**8,            &
        Qi1(1,k,1)*10**8,Qv1(1,k)*10**8,diag_vmi(1,k,1),Ql1(1,k,1)*10**8,     &
        tf1(1,k)-T0,Qg1(1,k,2)*10**8,Qi1(1,k,2)*10**8,Ql1(1,k,2)*10**8, &
        diag_3d(1,k,1)*10**8,diag_3d(1,k,2)*10**8

        write(24,'(1i5,2f7.0,11E20.5)') k, tsec,z(k),           &
        Nc1(1,k), Nr1(1,k), Ni1(1,k,1),          &
        Bg1(1,k,1)*10**8,diag_di(1,k,1),diag_rhoi(1,k,1), Ni1(1,k,2),Bg1(1,k,2)*10**8, &
        diag_3d(1,k,3),diag_3d(1,k,4),diag_rhoi(1,k,2)

      elseif (n_iceCat .eq.3) then
        write(20,'(1i5,2f7.0,20f17.5)') k,tsec,z(k),w(1,k),rho(k),p2d(1,k),         &
        tt1(1,k)-T0,td1(1,k)-T0,Qc1(1,k)*10**8,       &
        Qr1(1,k)*10**8,Qg1(1,k,1)*10**8,            &
        Qi1(1,k,1)*10**8,Qv1(1,k)*10**8,diag_vmi(1,k,1),Ql1(1,k,1)*10**8,     &
        tf1(1,k)-T0,diag_rhoi(1,k,1),                        &
        Qg1(1,k,2)*10**8,Qi1(1,k,2)*10**8,Ql1(1,k,2)*10**8,  &
        Qg1(1,k,3)*10**8,Qi1(1,k,3)*10**8,Ql1(1,k,3)*10**8

        write(24,'(1i5,2f7.0,10f25.5)') k, tsec,z(k),           &
        Nc1(1,k), Nr1(1,k), Ni1(1,k,1),          &
        Bg1(1,k,1)*10**8,diag_di(1,k,1),diag_rhoi(1,k,1),  &
        Ni1(1,k,2),Bg1(1,k,2)*10**8 ,&
        Ni1(1,k,3),Bg1(1,k,3)*10**8   

      elseif (n_iceCat .eq.4) then
        write(20,'(1i5,2f7.0,23f17.5)') k,tsec,z(k),w(1,k),rho(k),p2d(1,k),         &
        tt1(1,k)-T0,td1(1,k)-T0,Qc1(1,k)*10**8,       &
        Qr1(1,k)*10**8,Qg1(1,k,1)*10**8,            &
        Qi1(1,k,1)*10**8,Qv1(1,k)*10**8,diag_vmi(1,k,1),Ql1(1,k,1)*10**8,     &
        tf1(1,k)-T0,diag_rhoi(1,k,1),&
        Qg1(1,k,2)*10**8,Qi1(1,k,2)*10**8,Ql1(1,k,2)*10**8,&
        Qg1(1,k,3)*10**8,Qi1(1,k,3)*10**8,Ql1(1,k,3)*10**8,&
        Qg1(1,k,4)*10**8,Qi1(1,k,4)*10**8,Ql1(1,k,4)*10**8

        write(24,'(1i5,2f7.0,12f25.5)') k, tsec,z(k),&
        Nc1(1,k), Nr1(1,k), Ni1(1,k,1),&
        Bg1(1,k,1)*10**8,diag_di(1,k,1),diag_rhoi(1,k,1), Ni1(1,k,2),Bg1(1,k,2)*10**8,&
        Ni1(1,k,3),Bg1(1,k,3)*10**8,Ni1(1,k,4),Bg1(1,k,4)*10**8
     endif  !icecat
     enddo  !k loop
     endif  !outfreq
!*---------------------------------------------------------*!


     enddo  !main time loop

!----------------------------------------------------------!
!            ---  End of time integreation  ---            !
!----------------------------------------------------------!
      close(20)
      close(24)
      close(27)
      close(17)
      print*
      print*, 'DONE'

end subroutine columnmodel
!=============================================================================!

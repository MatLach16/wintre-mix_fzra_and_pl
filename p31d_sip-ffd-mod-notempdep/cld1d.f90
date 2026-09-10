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

      real, dimension(nk_input)    :: z,p,tt,td,tt2,w1,DIV,Tdry,zcld,w1cld,rho,alfa2
      real, dimension(ni,nk_input) :: tt0,tt1,Qv0,Qv1,w,SIGMA,Qsat,womega,th2d0,th2d1,p2d,dz2d,td1,tf1, dTdt
      real, dimension(lnmax) :: Pin,Zin,TTin,TDin,QVin,RHin

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



      !---------------------------------------------------------------------------------------!
      !Cond11 vars

      integer :: idx
      
      logical,parameter                   :: switch_use_pam = .false.
      logical,parameter                   :: switch_use_cosp = .false.
      ! logical,parameter                   :: switch_cond11   = .false. ! Use cond11 instead of P3

      integer                   :: ztmst_cond11
      integer,parameter         :: ntrac_cond11 = 2
      !integer                   :: ilg_cond11  : ni in cld90
      !integer                   :: il1_cond11  : kts in cld90
      !integer                   :: il2_cond11  : kte in cld90
      !integer                   :: ilev_cond11 : nk_input in cld90
      
      logical,parameter                  :: switch_use_tke = .false.
      
      real, dimension(ni,nk_input,ntrac_cond11) :: xrow_cond11
      real, dimension(ni,nk_input)              :: th_cond11
      real, dimension(ni,nk_input)              :: qv_cond11
      real, dimension(ni,nk_input)              :: rh_cond11
      real, dimension(ni,nk_input)              :: zclf_cond11
      real, dimension(ni,nk_input)              :: zcdn_cond11
      real, dimension(ni,nk_input)              :: dshj
      real, dimension(ni,nk_input)              :: shj
      real, dimension(ni,lev_cond11)      :: shtj
      real, dimension(ni)                 :: zrfl_cond11
      real, dimension(ni)                 :: zsfl_cond11
      real, dimension(ni)                 :: pressg_cond11
      real, dimension(ni,nk_input)              :: qlwc_cond11
      real, dimension(ni,nk_input)              :: qiwc_cond11
      real, dimension(ni,lev_cond11)      :: tf_cond11
      real, dimension(ni,nk_input)              :: dqldt_cond11
      real, dimension(ni,nk_input)              :: zfevap_cond11
      real, dimension(ni,nk_input)              :: zmratep_cond11
      real, dimension(ni,nk_input)              :: zfsnow_cond11
      real, dimension(ni,nk_input)              :: zmlwc_cond11
      real, dimension(ni,nk_input)              :: clrfr_cond11
      real, dimension(ni,nk_input)              :: zfrain_cond11
      real, dimension(ni,nk_input)              :: cvar_cond11
      real, dimension(ni,nk_input)              :: cvdu_cond11
      real, dimension(ni,nk_input)              :: cvsg_cond11
      real, dimension(nk_input)                 :: pblt_cond11
      real, dimension(ni,nk_input)              :: qc_cond11
      real, dimension(ni,nk_input)              :: rhc_cond11
      real, dimension(ni,nk_input)              :: clrfs_cond11
      real, dimension(ni,nk_input)              :: zfsubl_cond11
      real, dimension(ni,nk_input)              :: qcwvar_cond11
      real, dimension(ni,nk_input)              :: xlmtke_cond11
      real, dimension(ni,nk_input)              :: almc_cond11
      real, dimension(ni,nk_input)              :: almx_cond11

      real, dimension(ni,nk_input)              :: rmixrol_cond11
      real, dimension(ni,nk_input)              :: smixrol_cond11
      real, dimension(ni,nk_input)              :: rrefrol_cond11
      real, dimension(ni,nk_input)              :: srefrol_cond11
      real, dimension(ni,nk_input)              :: aggrol_cond11
      real, dimension(ni,nk_input)              :: autrol_cond11
      real, dimension(ni,nk_input)              :: cndrol_cond11
      real, dimension(ni,nk_input)              :: deprol_cond11
      real, dimension(ni,nk_input)              :: evprol_cond11
      real, dimension(ni,nk_input)              :: frhrol_cond11
      real, dimension(ni,nk_input)              :: frkrol_cond11
      real, dimension(ni,nk_input)              :: frsrol_cond11
      real, dimension(ni,nk_input)              :: mltirol_cond11
      real, dimension(ni,nk_input)              :: mltsrol_cond11
      real, dimension(ni,nk_input)              :: raclrol_cond11
      real, dimension(ni,nk_input)              :: rainrol_cond11
      real, dimension(ni,nk_input)              :: sacirol_cond11
      real, dimension(ni,nk_input)              :: saclrol_cond11
      real, dimension(ni,nk_input)              :: snowrol_cond11
      real, dimension(ni,nk_input)              :: subrol_cond11
      real, dimension(ni,nk_input)              :: sedirol_cond11
      real, dimension(ni,nk_input)              :: qtnc_cond11
      real, dimension(ni,nk_input)              :: hmnnc_cond11
      real, dimension(nk_input)              :: qv,rh
      
      real :: esw
!---------------------------------------------------------------------------------------!
    ! ****** output files *****
      open(unit=20, file='./outQcld.dat')
      open(unit=24, file='./outNcld.dat')
      open(unit=27, file='./sfcprec.dat')
    
    ! ****** output file headers*****
    write(27,'(4a10)') 'time', 'prt_liq', 'prt_sol', 'PR'
    
    if (n_iceCat .eq. 1) then
      write(20,'(16a10)') 'k','time','z','w','rho','p',&
      'T','Td','Qc','Qr','Qg','Qi','Qv',&
      'v_ice','Ql','tf1'
      write(24,'(9a10)') 'k','time','z','Nc','Nr',&
      'Ni','Bg','diag_di','diag_rhoi'

    elseif (n_iceCat .eq. 2) then
      write(20,'(25a10)') 'k','time','z','w','rho','p','T',&
      'Td','Qc','Qr','Qg','Qi','Qv','v_ice','Ql',&
      'tf1','Qg2','Qi2','Ql2','Qsedi_i1','Qsedi_i2','Qsedi_r',&
      'Qfluxi1','Qfluxi2','Qfluxr'
      write(24,'(18a15)') 'k','time','z','Nc','Nr','Ni',&
      'Bg','diag_di','diag_rhoi','Ni2','Bg2','Nsedi_i1','Nsedi_i2','Nsedi_r',&
      'Nfluxi1','Nfluxi2','Nfluxr','diag_rhoi2'

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

  write(9000,'(52a10)')  'time','k','ninuc','nimlt','nisub',&
        'nislf','nrhetc','nrheti','nchetc','ncheti','nimul','nlevp','nrhetic',&
        'nrslf','nrevp','ninuc2','nimlt2','nisub2',&
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
!     open(unit=12, file='./soundings/snd_input.Alta_3.data')   !grnd=  0m
!     open(unit=12, file='./soundings/Tsfc2.txt')
      open(unit=12, file='/home/rml000/Scripts/kin1d_matlach/soundings/temp_to_use.dat')

!  Input uninterpolated sounding data:
      read(12,*) nlvs
      read(12,*)
      read(12,*)
      read(12,*)
      read(12,*)
      read(12,*)
      do lv=1,nlvs
         !read(12,*) dum,Pin(lv),Zin(lv),TTin(lv),TDin(lv)
         read(12,*) Pin(lv),Zin(lv),TTin(lv),QVin(lv),RHin(lv)   !new
         TTin(lv) = TTin(lv) + T0		        !convert from C to Kelvins
	      !  TDin(lv) = TDin(lv) + T0
         !if(sndcase=='ALBERTA') TDin(lv) = TDin(lv) + T0 !        " "
         !if(sndcase=='HARP')    TDin(lv) = TDin(lv)      !'TD' is actually Qv
         Pin(lv)  = Pin(lv)*100                          !convert from mb to Pa
      enddo
      close(12)

!-- Set up vertical level (evenly-spaced z-levels)
!       RE-CODE if necessary

!!-- Read in levels from file:
!!     open(unit=20, file='./levels/levs_86.dat')
!!     open(unit=20, file='./levels/levs_62.dat')
      ! open(unit=19, file='../levels/levs_41.dat')

      ! read(19,*) nk_read
      ! if (nk_input /= nk_read) then
      !    print*, '*** Abort in CLD1D ***'
      !    print*, 'Mismatch in nk_input specified for arrays and nk_input from table of levels: ',nk_input,nk_read
      !    stop
      ! endif
      ! read(19,*)
      ! read(19,*)
      ! do kk = 1,nk_read
      !    read(19,*) k,shj(1,k),z(k),dum1,dum2
      ! enddo
      ! close(19)

!-- Compute dz2d(1,k) for upwind sedimentation.
!    note: this overrides the "DELTA_Z" in the file, which is not appropriate for P3
      if (trim(model)=='WRF') then
         do k = 1,nk_input-1
            dz2d(1,k) = z(k+1) - z(k)
         enddo
         dz2d(1,nk_input) = dz2d(1,nk_input-1)  !copy second highest DZ to highest
      else
         do k = 2,nk_input
            dz2d(1,k) = z(k-1) - z(k)
            dshj(1,k) = shj(1,k) - shj(1,k-1) ! ML to cond11
         enddo
         dz2d(1,1) = dz2d(1,2)  !copy second highest DZ to highest
         dshj(1,1) = dshj(1,2)
      endif
      
      
!  END OF INITIAL SET-UP

      if (switch_cond11 .eqv. .false.) then
      call P3_INIT(LT_path,n_iceCat,trplMomIce,liqFrac,model,stat,abort_on_err,dowr,bimm=bimm)      !v4.0.0_b42
      endif

!  Set up vertical level (evenly-spaced z-levels)
     z(nk_input) = Zin(1)
     delz = (H-Zin(1))/(nk_input-1)
     do k=1,nk_input
        z(k) = Zin(1) + delz*(nk_input-k)
        dz2d(1,k) = delz
     enddo
     dz = z(1)-z(2)
     dzsq = dz**2.
     H0= z(nk_input)
     GZ(1,:)= z*GRAV

!  Interpolate p,tt,td from sounding data and initialize Q[x]:
      do k=1,nk_input
         call vertint2b( p(k),z(k),Pin, Zin,nk_input,nlvs)
         call vertint2b(tt(k),z(k),TTin,Zin,nk_input,nlvs)
         call vertint2b(qv(k),z(k),QVin,Zin,nk_input,nlvs)
         call vertint2b(rh(k),z(k),RHin,Zin,nk_input,nlvs)

      !    print*, z(k),p(k), tt(k),qv(k) 
         rho(k)= p(k)/(Rd_cld*tt(k))
!        rho(k)= 1.
         p2d(1,k) = p(k)

      !    Qsat(1,k)= FOQST(tt(k),p(k))
      !    if(sndcase=='ALBERTA') Qv0(1,k) = FOQST(td(k),p(k))  !Alberta (Td read in)
         Qsat(1,k)= qv_sat(tt(k),p(k),0)

         esw      = polysvp1(tt(k),0)
         Qv1(1,k) = 0.622*rh(k)/100.*esw/p(k)

         Qv0(1,k) = min(qv(k), Qsat(1,k))!qv_sat(tt(k),p(k),0)
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

         tt0(1,k) = tt(k)
         tt1(1,k) = tt(k)

         Dmx(k) = 0.

        !  print*, p(k),z(k),tt(k),td(k)
      enddo


!**  z, w, p, tt, td and Qx0 are now in arrays from 1 to nk_input  **

      do k = 1,nk_input
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


      ! Set ups ML to cond11
      do k = 1,nk_input
          shtj(1,k) = shj(1,k) - dshj(1,k)/2.
      enddo
      shtj(1,lev_cond11)=shj(1,nk_input)+dshj(1,nk_input)/2.
            
      pressg_cond11 =  p2d(1,nk_input) /(shtj(1,lev_cond11))
!-------------------------------------------------------------------------!

!=========================================================================!
!  MAIN TIME LOOP:


      ! do k = nk_input-12, nk_input
      ! ! Qi1(1,k,1) = 10**(-5.)
      ! ! Ni1(1,k,1) = 10**6.
      ! Ni1(1,k,2) = 10**4.
      ! Qi1(1,k,2) = 0.475*10**(-5.)!10**4.*4.1*10**(-9)
      
      ! Qi0(1,k,:) = Qi1(1,k,:)
      ! Ni0(1,k,:) = Ni1(1,k,:)
      ! enddo 


      print*, 'Starting main loop...'

      do step = 1,nt

        tsec = step*dt_p3         !integration time [s]    (real)
        tminr= tsec/60.        !integration time [min]  (real)
        tmin = nint(tminr)     !integration time [min]  (integer)

      do k = 1,nk_input
       Qv1(1,k) = min(qv(k), Qsat(1,k))
      enddo

    ! ******* CHO ***********
    ! Initialize de la neige au top (qi0)
    !	do k = 2,2
    !	   Qi1(1,k,:) = 0.000265
    !	   Qg1(1,k,:) = 0.00017
    !	   Bg1(1,k,:) = Qg1(1,k,:)/500.
    !	   Ni1(1,k,:) = 5000
    !
    !	   Qi0(1,k,:) = Qi1(1,k,:)
    !	   Qg0(1,k,:) = Qg1(1,k,:)
    !	   Bg0(1,k,:) = Bg1(1,k,:)
    !	   Ni0(1,k,:) = Ni1(1,k,:)
    !	enddo
    ! ******* CHO ***********
    ! Initialize de la neige au top (qi0)
      do k = 3,3
        Qi1(1,k,1) = qitop !+ Qi1(1,k,1)
        Qg1(1,k,1) = qgtop !+ Qg1(1,k,1)
        Ql1(1,k,1) = 0
        !Qg1(1,k,:) = 0.00015
        !Qg1(1,k,:) = 0.00017
        Bg1(1,k,1) = Qg1(1,k,1)/900.
        Ni1(1,k,1) = nitop !+ Ni1(1,k,1)
      enddo

                    
!------------------------------------------------------------------------!
        if (microON) then

          th2d0(1,:) = tt0(1,:)*(1.e+5/p(:))**0.286
          th2d1(1,:) = tt1(1,:)*(1.e+5/p(:))**0.286

!v5.2.0

          call cpu_time(time1)

          if (.not. trplMomIce) Zi1(:,:,:) = 0.
          if (.not. liqFrac)    Ql1(:,:,:) = 0.

          if (switch_cond11 .eqv. .true.) then
            
            ! p3 vars to cond11 vars!------------------------------------------------------------------------!
            
            ! Calculate average between each level for interface temperature
            do k=2,nk_input
                tf_cond11(1,k) = (th_cond11(1,k-1)+th_cond11(1,k))/2.
            enddo         
            tf_cond11(1,1) = tf_cond11(1,2) 
          
            do k=1,nk_input
            xrow_cond11(1,k,2) = Qi1(1,k,1) 
            xrow_cond11(1,k,1) = Qc1(1,k) 
            
            th_cond11(1,k) = tt1(1,k)
            qv_cond11(1,k) = Qv1(1,k)
            rh_cond11(1,k)   = 0.
            zclf_cond11(1,k) = 0.
            zcdn_cond11(1,k) = Nc1(1,k)
            !zrfl_cond11 = 0.
            !zsfl_cond11 = 0.
            !qlwc_cond11 = 0.
            !qiwc_cond11 = 0.
            
            !dqldt_cond11   = 0. !Convective detrainment water
            !zfevap_cond11  = 0.
            !zmratep_cond11 = 0.
            !zfsnow_cond11  = 0.
            !zmlwc_cond11   = 0.
            !clrfr_cond11   = 0.
            !zfrain_cond11  = 0.
            cvar_cond11(1,k) = 0. ! output of stacld5
            cvdu_cond11(1,k) = 0. ! input for variance from convection in statcld5
            cvsg_cond11(1,k) = 0. ! convective contribution to variance in statcld5
            !pblt_cond11(1,k) = 0. ! not used
            !qc_cond11      = 0. 
            !rhc_cond11     = 0.
            !clrfs_cond11   = 0.
            !zfsubl_cond11  = 0.
            !qcwvar_cond11 = 0. ! output of statcld5, not used in cond11
            xlmtke_cond11(1,k) = 1.
            almc_cond11(1,k) = almx_in
            almx_cond11(1,k) = almx_in
            !ztmst_cond11 = timestep
            !ntrac_cond11 =
            !ilg_cond11 =
            !il1_cond11 =
            !il2_cond11 =
            !ilev_cond11 =
            ! lev_cond11 =
            !rmixrol_cond11 =  0.
            !smixrol_cond11 =  0.
            !rrefrol_cond11 =  0.
            !srefrol_cond11 =  0.
            !aggrol_cond11 =  0.
            !autrol_cond11 =  0.
            !cndrol_cond11 =  0.
            !deprol_cond11 =  0.
            !evprol_cond11 =  0.
            !frhrol_cond11 =  0.
            !frkrol_cond11 =  0.
            !frsrol_cond11 =  0.
            !mltirol_cond11 = 0.
            !mltsrol_cond11 = 0.
            !raclrol_cond11 = 0.
            !rainrol_cond11 = 0.
            !sacirol_cond11 = 0.
            !saclrol_cond11 = 0.
            !snowrol_cond11 = 0.
            !subrol_cond11  = 0.
            !sedirol_cond11 = 0.
            !qtnc_cond11 = 0. ! not initialized in cond11, maybe remove this line 
            !hmnnc_cond11 = 0. ! not initialized in cond11, maybe remove this line
            enddo
            !-----------------------------Call cond11-------------------------------------------!       
            CALL COND11(xrow_cond11,th_cond11,qv_cond11,rh_cond11,zclf_cond11,zcdn_cond11,dshj,shj,shtj, &
                        zrfl_cond11,zsfl_cond11,pressg_cond11,qlwc_cond11, &
                        qiwc_cond11,tf_cond11,dqldt_cond11,zfevap_cond11,zmratep_cond11,zfsnow_cond11,zmlwc_cond11, &
                        clrfr_cond11,zfrain_cond11,cvar_cond11,cvdu_cond11,cvsg_cond11,pblt_cond11,qc_cond11,rhc_cond11, &
                        clrfs_cond11,zfsubl_cond11,qcwvar_cond11,xlmtke_cond11, &
                        almc_cond11,almx_cond11,dt_p3,ntrac_cond11,ni,its,ite,nk_input,lev_cond11,&
                        rmixrol_cond11,smixrol_cond11,rrefrol_cond11,srefrol_cond11, &
                        aggrol_cond11,autrol_cond11,cndrol_cond11,deprol_cond11, &
                        evprol_cond11,frhrol_cond11,frkrol_cond11,frsrol_cond11, &
                        mltirol_cond11,mltsrol_cond11,raclrol_cond11,rainrol_cond11, &
                        sacirol_cond11,saclrol_cond11,snowrol_cond11,subrol_cond11,sedirol_cond11, &
                        qtnc_cond11,hmnnc_cond11,switch_use_tke,switch_use_pam, &
                        switch_use_cosp)
            ! cond11 vars to p3 vars!------------------------------------------------------------------------!
              
              tt1(1,:) = th_cond11(1,:)

              prt_liq  = zrfl_cond11 / rhow *3.6e+6! add conversion to m/s
              prt_sol  = zsfl_cond11 / rhow *3.6e+6! add conversion to m/s
              prt_drzl = 0.
              prt_rain = zrfl_cond11 / rhow *3.6e+6
              prt_crys = 0.
              prt_snow = zsfl_cond11 / rhow * 3.6e+6
              prt_grpl = 0.
              prt_pell = 0.
              prt_hail = 0.
              prt_sndp = 0.
              prt_wsnow =0.
              do k=1,nk_input 
              Qv1(1,k)    = qv_cond11(1,k)
              Qc1(1,k)    = xrow_cond11(1,k,1)
              Qr1(1,k)    = 0.
              Nc1(1,k)    = zcdn_cond11(1,k)
              Nr1(1,k)    = 0.
              Qi1(1,k,1)    = xrow_cond11(1,k,2)
              Qg1(1,k,1)    = 0.
              Ql1(1,k,1)    = 0.
              Ni1(1,k,1)    = 0.
              Bg1(1,k,1)    = 0.
              Zi1(1,k,1)    = 0.
              enddo
            
          !------------------------------------------------------------------------!
          else
            ! print*, 'timestep: ', step
            ! print*, Qi1(1,2,1),Ni1(1,2,1)
            diag_3d(:,:,:) = 0.
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
          
          
          endif

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
        write(20,'(1i5,2f7.0,22E20.5)') k,tsec,z(k),w(1,k),rho(k),p2d(1,k),         &
        tt1(1,k)-T0,td1(1,k)-T0,Qc1(1,k)*10**8,       &
        Qr1(1,k)*10**8,Qg1(1,k,1)*10**8,            &
        Qi1(1,k,1)*10**8,Qv1(1,k)*10**8,diag_vmi(1,k,1),Ql1(1,k,1)*10**8,     &
        tf1(1,k)-T0,Qg1(1,k,2)*10**8,Qi1(1,k,2)*10**8,Ql1(1,k,2)*10**8, &
        diag_3d(1,k,1)*10**8,diag_3d(1,k,2)*10**8,diag_3d(1,k,5)*10**8,&
        diag_3d(1,k,7)*10**8,diag_3d(1,k,8)*10**8,diag_3d(1,k,11)*10**8

        write(24,'(1i5,2f7.0,16E20.5)') k, tsec,z(k),           &
        Nc1(1,k), Nr1(1,k), Ni1(1,k,1),          &
        Bg1(1,k,1)*10**8,diag_di(1,k,1),diag_rhoi(1,k,1), Ni1(1,k,2),Bg1(1,k,2)*10**8, &
        diag_3d(1,k,3),diag_3d(1,k,4),diag_3d(1,k,6),&
        diag_3d(1,k,9),diag_3d(1,k,10),diag_3d(1,k,12),diag_rhoi(1,k,2)

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

      print*
      print*, 'DONE'

end subroutine columnmodel
!=============================================================================!

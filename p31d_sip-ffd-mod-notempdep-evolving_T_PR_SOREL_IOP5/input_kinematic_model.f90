module input_kinematic_model
implicit none
real, parameter :: dtadv = 60.000000 
logical, parameter   :: switch_cond11 = .false. 
logical, parameter   :: skip_advection = .true. 
character(len=20), parameter   :: switch_micro = "P3" 
integer, parameter   :: isubg = 1 
real, parameter   :: almx_in = 600.000000 
real, parameter   :: AMPA = 2.000000 
real, parameter   :: AMPB = 5.000000 
integer, parameter   :: wprofile = 1 
integer, parameter   :: ttotmin = 610 
integer, parameter   :: n_iceCat = 2 
real, parameter :: qitop = 0.000000 
real, parameter :: qltop = 0.000000 
real, parameter :: qgtop = 0.000000 
real, parameter :: bgtop = 0.000000 
real, parameter :: nitop = 0.000000 
integer, parameter :: nk_input = 102 
integer, parameter :: lev_cond11 = 103 
real, parameter :: dt_p3 = 1.000000 
integer, parameter :: outfreq = 1 
integer, parameter :: iparam = 1 
logical, parameter :: switch_int_so4 = .false. 
logical, parameter :: log_predictNc_1d = .false. 
logical, parameter :: liqFrac = .false. 
real, parameter :: bimm = 2 
logical, parameter :: trplMomIce = .false. 
end module input_kinematic_model

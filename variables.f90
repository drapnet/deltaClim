!=======================================================================!
!    This file is part of interpolation_myo2cst                         !
!    Engineered by Manel Grifoll (LIM/UPC)                              !
!                                                                       !
!=======================================================================!

MODULE variables
USE mod_kinds
implicit none
SAVE


CHARACTER (LEN=180) :: GRD_NAME_o
character (LEN=180) :: nom_clm
character (LEN=180) :: nom_ini
character (LEN=180) :: nom_bry
character (LEN=180) :: dir_in
character (LEN=180) :: dir_out
character (LEN=180) :: arx_dm
character (LEN=180) :: arx_hm
character (LEN=180) :: arx_out

integer :: Nx_o,Nx_i         ! dim xi_rho,xlon
integer :: Ny_o,Ny_i         ! dim eta_rho,ylat
integer :: l_time_iq          ! dim variables diaries   sera quarts d'hora
integer :: l_time_ih         ! dim variables horaries   sera 24


real(r8),parameter::f_v=1.0E20_r8        !per myocean_ibi

!i: refered input variables
real(r8),DIMENSION(:,:),ALLOCATABLE :: lon_i,lat_i
real(r8),DIMENSION(:),ALLOCATABLE :: depth_i
real(r8),DIMENSION(:,:),ALLOCATABLE :: angle_i

!o: referred output variables
real(r8),DIMENSION(:,:),ALLOCATABLE :: lon_r_o,lat_r_o,lon_u_o,lat_u_o,lon_v_o,lat_v_o
real(r8),DIMENSION(:,:),ALLOCATABLE :: mask_r_o, mask_u_o, mask_v_o
real(r8),DIMENSION(:,:),ALLOCATABLE :: h_o, angle_o,h_u_o,h_v_o

!!!!!!!!!!!! acaben definicions de la malla !!!!!!!!!!!!!!!!!!!
! Variables per la depth de nivells sigma
!    ES problabe que alguns vectors estiguin en un histori com atribut
!    una manera de trure'l es:
!    CALL check(nf90_get_att(ncid,NF90_global,"Cs_r",Cs_r))

integer :: Niv_i,Niv_o

!Variables al output file
real(r8) :: hc_o,Tcline_o,ocean_time
real(r8) :: theta_s_o, theta_b_o           ! de la sortida
real(r8),DIMENSION(:),ALLOCATABLE :: Cs_r_o      !  "S-coordinate stretching curves at RHO-points" ;
real(r8),DIMENSION(:),ALLOCATABLE :: sc_r_o      !  "S-coordinate at W-points" 
real(r8),DIMENSION(:),ALLOCATABLE :: Cs_w_o       !  "S-coordinate stretching curves at W-points" ;
real(r8),DIMENSION(:),ALLOCATABLE :: sc_w_o      !  "S-coordinate at W-points" 
real(r8),DIMENSION(:,:,:),ALLOCATABLE :: z_r_o
real(r8),DIMENSION(:,:,:),ALLOCATABLE :: z_u_o
real(r8),DIMENSION(:,:,:),ALLOCATABLE :: z_v_o
real(r8),DIMENSION(:,:,:),ALLOCATABLE ::  z_w_o
real(r8),DIMENSION(:,:,:),ALLOCATABLE ::  Hz_o
real(r8),DIMENSION(:,:),ALLOCATABLE ::  Zt_avg1_o !Per ajustar el ssh.
 
!Constants
real(r8), parameter :: pi = 3.14159265358979323846_r8
real(r8), parameter :: deg2rad = pi / 180.0_r8
real(r8), parameter :: rad2deg = 180.0_r8 / pi
real(r8), parameter :: day2sec = 86400.0_r8
real(r8), parameter :: sec2day = 1.0_r8 / 86400.0_r8
real(r8) :: Eradius = 6371315.0_r8      ! m

END MODULE variables

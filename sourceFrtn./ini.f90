
SUBROUTINE ini

USE mod_kinds
USE variables
use funcions
USE netcdf
implicit none
integer :: ncid,dimid,varid,i,j,ncidg
character(LEN=180) ::nom_dim
character(LEN=180) ::file_in

character(LEN=30) :: second 
real(r8),DIMENSION(:),ALLOCATABLE :: lon,lat   ! els inis 
print *, "RUNNING: ./interpFtn"
print *, "from SAMOA System (M.Grfioll, LIM_UPC)"
print *, "PART 1: Reading parameters from params.in"

!!                    OPEN(5,FILE='params.F')
!!read(5,'(A80)') GRD_NAME_o
!!read(5,'(A80)') dir_in
!!read(5,'(A80)') dir_out
!!read(5,'(A80)') arx_hm      !arxiu que ve de de la fase 1 es unclim en la malla cnm es tot-h.nc
!!read(5,*) l_time_ih         !dimTime
!!close(5)
arx_hm="./dat/cur-hor.nc"
!!dir_in=
GRD_NAME_o="./in/grd_cst_delta.nc"
nom_clm="./dat/aux9.nc"
l_time_ih=48
l_time_iq=192
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Assign mesh values:
!theta_s_o=7_r8
!theta_b_o=0.4_r8
!hc_o=3_r8
!Tcline_o=3_r8
Niv_o=15
ocean_time=0.0_r8   

print *, "./dat/cur_hor.nc"
print *, l_time_ih
print *, "aixo es la dimTime"
!!!!!!!!!!!1 Agafem un arxiu cmems per trobar  lon lat , dept  (ha de ser 3d exem cur_hor
CALL check(nf90_open("./dat/cur-hor.nc",NF90_NOWRITE,ncid))
CALL check(nf90_inq_dimid(ncid,'lon',dimid))
CALL check(nf90_inquire_dimension(ncid,dimid,nom_dim,Nx_i))
CALL check(nf90_inq_dimid(ncid,'lat',dimid))
CALL check(nf90_inquire_dimension(ncid,dimid,nom_dim,Ny_i))
CALL check(nf90_inq_dimid(ncid,'depth',dimid))
CALL check(nf90_inquire_dimension(ncid,dimid,nom_dim,Niv_i))
ALLOCATE(lon(1:Nx_i)) 
ALLOCATE(lat(1:Ny_i)) 
ALLOCATE(depth_i(1:Niv_i))
print *,Nx_i,Ny_i
print *, "malla in"
print *, "Nx_i" 
CALL check(nf90_inq_varid(ncid,'lon',varid))
CALL check(nf90_get_var(ncid,varid,lon))
CALL check(nf90_inq_varid(ncid,'lat',varid))
CALL check(nf90_get_var(ncid,varid,lat))
CALL check(nf90_inq_varid(ncid,'depth',varid))
CALL check(nf90_get_var(ncid,varid,depth_i))
CALL check(nf90_close(ncid))
call get_grids(TRIM(GRD_NAME_o))  
print * , "parametres del nou grid obtinugudes"
!!!!!!  ara tindrem les del makeClmPhase1.py
! treim les dimensions de nivell, 

                                         


!!Meshgrid to interpolate
DO i=1,Nx_i
   DO j=1,Ny_i
      lon_i(i,j)=lon(i)
      lat_i(i,j)=lat(j)
   END DO
END DO
print * , "Meshgrid ended: Assaigned lon_i and lat_i"  

!!Comprobacio que  les limits del inputs son correctes per a la malla oprint
if(MAXVAL(lon_r_o)>MAXVAL(lon_i)) THEN
	print*, "Tenim un problema (  O --- I ) falten inputs a la dreta"
	print *,MAXVAL(lon_r_o),MAXVAL(lon_i)
	stop
endif
if(MAXVAL(lat_r_o)>MAXVAL(lat_i)) THEN
	print*, "Tenim un problema (  O --- I )  falten inputs a dalt"
	print *,MAXVAL(lat_r_o),MAXVAL(lat_i)
	stop
endif
!Valors minims
if(MINVAL(lon_r_o)<MINVAL(lon_i)) THEN
	print*, "Tenim un problema  (  O --- I ) falten inputs a la esquerra"
	print *,MINVAL(lon_r_o),MINVAL(lon_i)
	stop
endif
if(MINVAL(lat_r_o)<MINVAL(lat_i)) THEN
	print*, "Tenim un problema  (  O --- I ) falten inputs a baix"
	print *,MINVAL(lat_r_o),MINVAL(lat_i)
	stop
endif



if(MAXVAL(h_o)>MAXVAL(depth_i)) THEN
	print *,"Problema: necesitem imputs mes profundsMira ini.f90"
	print *," MAXVAL(h_o)=", MAXVAL(h_o)
	print *," MAXVAL(depth_i)=", MAXVAL(depth_i)
	stop            
endif
ALLOCATE(Hz_o(1:Nx_o,1:Ny_o,1:Niv_o))
ALLOCATE(z_w_o(1:Nx_o,1:Ny_o,0:Niv_o))
ALLOCATE(z_r_o(1:Nx_o,1:Ny_o,1:Niv_o))
ALLOCATE(z_u_o(1:Nx_o-1,1:Ny_o,1:Niv_o))
ALLOCATE(z_v_o(1:Nx_o,1:Ny_o-1,1:Niv_o))


ALLOCATE(Cs_r_o(1:Niv_o))
ALLOCATE(Cs_w_o(0:Niv_o))
ALLOCATE(sc_r_o(1:Niv_o))
ALLOCATE(sc_w_o(0:Niv_o))

ALLOCATE(Zt_avg1_o(1:Nx_o,1:Ny_o))
Zt_avg1_o=0.0_r8
print *,"Executing set_depths"
call set_depths

!call set_depths
END SUBROUTINE ini  

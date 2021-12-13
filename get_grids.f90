!=======================================================================!
!    This file is part of SAMOA                                       !
!    Engineered by Manel Grifoll (LIM)                                  !
!                                                                       !
!=======================================================================!

SUBROUTINE get_grids(nom_grd_o)
! Aquesta subroutina treu les variables mes corrents de una grid
! Nx   correspon a la dim xi_rho
! Ny   corrspon a la dim eta_rho
! lon_r,lon_u, lon_v  evident
! lat_r, lat_u, lat_v
! mask_r mask_u, mask_v     
! h     no canvia de nom
! angle 
! aquestes variables estan definides al mod variables
! llestes per a ser utilitzades amb un main
! Aqui una grid i uns arxoi del nyo
USE mod_kinds
USE variables
USE funcions
USE netcdf
implicit none
character(*) , intent(in) ::nom_grd_o
character(LEN=180) ::nom_dim
integer :: ncid,dimid,varid,i,j

ALLOCATE(lon_i(1:Nx_i,1:Ny_i))
ALLOCATE(lat_i(1:Nx_i,1:Ny_i))
!ALLOCATE(depth_i(1:NIV_i))
!ALLOCATE(mask_r_i(1:Nx_i,1:Ny_i))
!ALLOCATE(h_i(1:Nx_i,1:Ny_i))
!ALLOCATE(lon_u_i(1:Nx_i-1,1:Ny_i))
!ALLOCATE(lat_u_i(1:Nx_i-1,1:Ny_i))
ALLOCATE(angle_i(1:Nx_i,1:Ny_i))
angle_i=0.0_r8
CALL check(nf90_open(TRIM(nom_grd_o),NF90_NOWRITE,ncid))
CALL check(nf90_inq_dimid(ncid,'xi_rho',dimid))
CALL check(nf90_inquire_dimension(ncid,dimid,nom_dim,Nx_o))
CALL check(nf90_inq_dimid(ncid,'eta_rho',dimid))
CALL check(nf90_inquire_dimension(ncid,dimid,nom_dim,Ny_o))
ALLOCATE(lon_r_o(1:Nx_o,1:Ny_o))
ALLOCATE(lat_r_o(1:Nx_o,1:Ny_o))
ALLOCATE(mask_r_o(1:Nx_o,1:Ny_o))
ALLOCATE(h_o(1:Nx_o,1:Ny_o))
ALLOCATE(lon_u_o(1:Nx_o-1,1:Ny_o))
ALLOCATE(lat_u_o(1:Nx_o-1,1:Ny_o))
ALLOCATE(mask_u_o(1:Nx_o-1,1:Ny_o))
ALLOCATE(h_u_o(1:Nx_o-1,1:Ny_o))
ALLOCATE(lon_v_o(1:Nx_o,1:Ny_o-1))
ALLOCATE(lat_v_o(1:Nx_o,1:Ny_o-1))
ALLOCATE(mask_v_o(1:Nx_o,1:Ny_o-1))
ALLOCATE(h_v_o(1:Nx_o,1:Ny_o-1))
ALLOCATE(angle_o(1:Nx_o,1:Ny_o))
CALL check(nf90_inq_varid(ncid,'lon_rho',varid))
CALL check(nf90_get_var(ncid,varid,lon_r_o))
CALL check(nf90_inq_varid(ncid,'lat_rho',varid))
CALL check(nf90_get_var(ncid,varid,lat_r_o))
CALL check(nf90_inq_varid(ncid,'mask_rho',varid))
CALL check(nf90_get_var(ncid,varid,mask_r_o))
CALL check(nf90_inq_varid(ncid,'lon_u',varid))
CALL check(nf90_get_var(ncid,varid,lon_u_o))
CALL check(nf90_inq_varid(ncid,'lat_u',varid))
CALL check(nf90_get_var(ncid,varid,lat_u_o))
CALL check(nf90_inq_varid(ncid,'mask_u',varid))
CALL check(nf90_get_var(ncid,varid,mask_u_o))
CALL check(nf90_inq_varid(ncid,'lon_v',varid))
CALL check(nf90_get_var(ncid,varid,lon_v_o))
CALL check(nf90_inq_varid(ncid,'lat_v',varid))
CALL check(nf90_get_var(ncid,varid,lat_v_o))
CALL check(nf90_inq_varid(ncid,'mask_v',varid))
CALL check(nf90_get_var(ncid,varid,mask_v_o))
CALL check(nf90_inq_varid(ncid,'h',varid))
CALL check(nf90_get_var(ncid,varid,h_o))
CALL check(nf90_inq_varid(ncid,'angle',varid))
CALL check(nf90_get_var(ncid,varid,angle_o))
CALL check(nf90_close(ncid))
DO i=1,Nx_o-1
   h_u_o(i,:)=0.5_r8*(h_o(i,:)+h_o(i+1,:))
ENDDO
DO j=1,Ny_o-1
   h_v_o(:,j)=0.5_r8*(h_o(:,j)+h_o(:,j+1))
ENDDO
print * ,"Parameters from grid file obtained"
END SUBROUTINE get_grids

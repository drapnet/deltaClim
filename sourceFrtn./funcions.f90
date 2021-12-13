!==================================================================!
!    This file is part of SAMOA                                         !
!    Engineered by Manel Grifoll (LIM)                                  !
!                                                                       !
!=======================================================================!

MODULE funcions

USE variables
USE mod_kinds
CONTAINS
SUBROUTINE check(status)
 use netcdf
 integer,intent(in)::status
if (status /= nf90_noerr)then
	print *,trim(nf90_strerror(status))
 	stop "NETcdf problem"
endif
END SUBROUTINE check
!!!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

REAL(r8) function linial(x1,v1,x2,v2,x)
implicit none
REAL(r8),INTENT(IN)::x1,x2,v1,v2,x
linial=((x-x2)*v1+(x1-x)*v2)/(x1-x2)
end function linial
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
!!!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!!!!! 
SUBROUTINE  create_clm(nom_arx,time_iq,time_ih)
!time_ih  es pel zeta,ubar,vbar  com  horaris
USE variables
USE mod_kinds
USE netcdf
implicit none
character(*),intent(in) :: nom_arx
integer,intent(in)::time_iq,time_ih
integer:: d_temp,d_salt,d_v2d,d_v3d,d_zeta
integer::ncid,d_xiu,d_etau,d_xiv,d_xirho,d_etarho,d_srho,d_sw,d_etav,l_one,varid,d_one
integer :: cycle,i
l_one=1
cycle=360
CALL check(nf90_create(TRIM(nom_arx),NF90_CLOBBER,ncid))
CALL check(nf90_DEF_dim(ncid,'xi_u',Nx_o-1,d_xiu))
CALL check(nf90_DEF_dim(ncid,'eta_u',Ny_o,d_etau))
CALL check(nf90_DEF_dim(ncid,'xi_v',Nx_o,d_xiv))
CALL check(nf90_DEF_dim(ncid,'eta_v',Ny_o-1,d_etav))
CALL check(nf90_DEF_dim(ncid,'xi_rho',Nx_o,d_xirho))
CALL check(nf90_DEF_dim(ncid,'eta_rho',Ny_o,d_etarho))
CALL check(nf90_DEF_dim(ncid,'s_rho',Niv_o,d_srho))
CALL check(nf90_DEF_dim(ncid,'s_w',Niv_o+1,d_sw))
!CALL check(nf90_DEF_dim(ncid,'tracer',l_tracer,d_tracer))
CALL check(nf90_DEF_dim(ncid,'temp_time',time_ih,d_temp))
CALL check(nf90_DEF_dim(ncid,'salt_time',time_ih,d_salt))
CALL check(nf90_DEF_dim(ncid,'v2d_time',time_iq,d_v2d))
CALL check(nf90_DEF_dim(ncid,'v3d_time',time_ih,d_v3d))
CALL check(nf90_DEF_dim(ncid,'zeta_time',time_iq,d_zeta))    ! 
CALL check(nf90_DEF_dim(ncid,'one',l_one,d_one))

! Crea variables i atributs
CALL check(nf90_def_var(ncid,'theta_s',NF90_DOUBLE,d_one,varid))
CALL check(nf90_put_att(ncid,varid,"long_name","S-coordinate surface control parameter"))
CALL check(nf90_put_att(ncid,varid,"units","nondimensional"))
CALL check(nf90_def_var(ncid,'theta_b',NF90_DOUBLE,d_one,varid))
CALL check(nf90_put_att(ncid,varid,"long_name","S-coordinate bottom control parameter"))
CALL check(nf90_put_att(ncid,varid,"units","nondimensional"))
 
CALL check(nf90_def_var(ncid,'Tcline',NF90_DOUBLE,d_one,varid))
CALL check(nf90_put_att(ncid,varid,"long_name","S-coordinate surface/bottom layer width"))
CALL check(nf90_put_att(ncid,varid,"units","meter"))

CALL check(nf90_def_var(ncid,'hc',NF90_DOUBLE,d_one,varid))
CALL check(nf90_put_att(ncid,varid,"long_name","S-coordinate parameter,critical depth"))
CALL check(nf90_put_att(ncid,varid,"units","meter"))

CALL check(nf90_def_var(ncid,'h_5',NF90_DOUBLE,(/d_sw,d_xirho,d_etarho/),varid))
CALL check(nf90_put_att(ncid,varid,"long_name","Aux per z_5"))
CALL check(nf90_put_att(ncid,varid,"units","meter"))

CALL check(nf90_def_var(ncid,'sc_r',NF90_DOUBLE,d_srho,varid))
CALL check(nf90_put_att(ncid,varid,"long_name","S-coordinate at RHO-points"))
CALL check(nf90_put_att(ncid,varid,"units","nondimensional"))
CALL check(nf90_put_att(ncid,varid,"valid_min","-1"))
CALL check(nf90_put_att(ncid,varid,"valid_min","0"))

CALL check(nf90_def_var(ncid,'Cs_r',NF90_DOUBLE,d_srho,varid))
CALL check(nf90_put_att(ncid,varid,"long_name","S-coordinate stretching curves  at RHO-points"))
CALL check(nf90_put_att(ncid,varid,"units","nondimensional"))
CALL check(nf90_put_att(ncid,varid,"valid_min","-1"))
CALL check(nf90_put_att(ncid,varid,"valid_min","0"))

CALL check(nf90_def_var(ncid,'temp_time',NF90_DOUBLE,d_temp,varid))
CALL check(nf90_put_att(ncid,varid,"long_name","time for temperature climatology"))
CALL check(nf90_put_att(ncid,varid,"units","day"))
CALL check(nf90_put_att(ncid,varid,"cycle_leng ",cycle))

CALL check(nf90_def_var(ncid,'salt_time',NF90_DOUBLE,d_salt,varid))
CALL check(nf90_put_att(ncid,varid,"long_name","time for salinity climatology"))
CALL check(nf90_put_att(ncid,varid,"units","day"))
CALL check(nf90_put_att(ncid,varid,"cycle_leng ",cycle))

CALL check(nf90_def_var(ncid,'v2d_time',NF90_DOUBLE,d_v2d,varid))
CALL check(nf90_put_att(ncid,varid,"long_name","time for 2D velocity climatology"))
CALL check(nf90_put_att(ncid,varid,"units","day"))
CALL check(nf90_put_att(ncid,varid,"cycle_leng ",cycle))

CALL check(nf90_def_var(ncid,'v3d_time',NF90_DOUBLE,d_v3d,varid))
CALL check(nf90_put_att(ncid,varid,"long_name","time for 3D velocity climatology"))
CALL check(nf90_put_att(ncid,varid,"units","day"))
CALL check(nf90_put_att(ncid,varid,"cycle_leng ",cycle))


CALL check(nf90_def_var(ncid,'zeta_time',NF90_DOUBLE,d_zeta,varid))
CALL check(nf90_put_att(ncid,varid,"long_name","time for sea surface height"))
CALL check(nf90_put_att(ncid,varid,"units","day"))
CALL check(nf90_put_att(ncid,varid,"cycle_leng ",cycle))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CALL  check(nf90_def_var(ncid,'u',NF90_DOUBLE,(/d_xiu,d_etau,d_srho,d_v3d/),varid))
CALL check(nf90_put_att(ncid,varid,"long_name","u-momentum component"))
CALL check(nf90_put_att(ncid,varid,"units","meter second-1"))

CALL  check(nf90_def_var(ncid,'v',NF90_DOUBLE,(/d_xiv,d_etav,d_srho,d_v3d/),varid))
CALL check(nf90_put_att(ncid,varid,"long_name","v-momentum component"))
CALL check(nf90_put_att(ncid,varid,"units","meter second-1"))
i=1
CALL  check(nf90_def_var(ncid,'ubar',NF90_DOUBLE,(/d_xiu,d_etau,d_v2d/),varid))
CALL check(nf90_put_att(ncid,varid,"long_name","verticaly integrated u-momentum component"))
CALL check(nf90_put_att(ncid,varid,"units","meter second-1"))

CALL  check(nf90_def_var(ncid,'vbar',NF90_DOUBLE,(/d_xiv,d_etav,d_v2d/),varid))
CALL check(nf90_put_att(ncid,varid,"long_name","verticaly integrated v-momentum component"))
CALL check(nf90_put_att(ncid,varid,"units","meter second-1"))

CALL  check(nf90_def_var(ncid,'zeta',NF90_DOUBLE,(/d_xirho,d_etarho,d_zeta/),varid))
CALL check(nf90_put_att(ncid,varid,"long_name","free-surface"))
CALL check(nf90_put_att(ncid,varid,"units","meter"))


CALL  check(nf90_def_var(ncid,'temp',NF90_DOUBLE,(/d_xirho,d_etarho,d_srho,d_temp/),varid))
CALL check(nf90_put_att(ncid,varid,"long_name","potencial temperature"))
CALL check(nf90_put_att(ncid,varid,"units","Celsius"))

CALL  check(nf90_def_var(ncid,'salt',NF90_DOUBLE,(/d_xirho,d_etarho,d_srho,d_salt/),varid))
CALL check(nf90_put_att(ncid,varid,"long_name","salinity"))
CALL check(nf90_put_att(ncid,varid,"units","PSU"))
CALL check(nf90_enddef(ncid))

CALL check(nf90_inq_varid(ncid,"theta_b",varid))
CALL check(nf90_put_var(ncid,varid,theta_b_o))
CALL check(nf90_inq_varid(ncid,"theta_s",varid))
CALL check(nf90_put_var(ncid,varid,theta_s_o))
CALL check(nf90_inq_varid(ncid,"Tcline",varid))
CALL check(nf90_put_var(ncid,varid,Tcline_o))
CALL check(nf90_inq_varid(ncid,"hc",varid))
CALL check(nf90_put_var(ncid,varid,hc_o))
CALL check(nf90_inq_varid(ncid,"sc_r",varid))
CALL check(nf90_put_var(ncid,varid,sc_r_o))
CALL check(nf90_inq_varid(ncid,"Cs_r",varid))
CALL check(nf90_put_var(ncid,varid,Cs_r_o))
CALL check(nf90_close(ncid))

END SUBROUTINE  create_clm

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE  create_bry(nom_arx,time_i,time_ih)
!time_ih  es pel zeta com unic horari
USE variables
USE mod_kinds
USE netcdf
implicit none
character(*),intent(in) :: nom_arx
integer,intent(in)::time_i,time_ih
integer:: d_tclm,d_temp,d_sclm,d_salt,d_uclm,d_vclm,d_v2d,d_v3d,d_ssh,d_zeta
integer::ncid,d_xiu,d_etau,d_xiv,d_xirho,d_etarho,d_srho,d_sw,d_etav,l_one,varid,d_one
integer :: cycle,i
l_one=1
cycle=360
CALL check(nf90_create(TRIM(nom_arx),NF90_CLOBBER,ncid))
CALL check(nf90_DEF_dim(ncid,'xi_u',Nx_o-1,d_xiu))
CALL check(nf90_DEF_dim(ncid,'eta_u',Ny_o,d_etau))
CALL check(nf90_DEF_dim(ncid,'xi_v',Nx_o,d_xiv))
CALL check(nf90_DEF_dim(ncid,'eta_v',Ny_o-1,d_etav))
CALL check(nf90_DEF_dim(ncid,'xi_rho',Nx_o,d_xirho))
CALL check(nf90_DEF_dim(ncid,'eta_rho',Ny_o,d_etarho))
CALL check(nf90_DEF_dim(ncid,'s_rho',Niv_o,d_srho))
CALL check(nf90_DEF_dim(ncid,'s_w',Niv_o+1,d_sw))
!CALL check(nf90_DEF_dim(ncid,'tracer',l_tracer,d_tracer))
CALL check(nf90_DEF_dim(ncid,'temp_time',time_i,d_temp))
CALL check(nf90_DEF_dim(ncid,'salt_time',time_i,d_salt))
CALL check(nf90_DEF_dim(ncid,'v2d_time',time_ih,d_v2d))
CALL check(nf90_DEF_dim(ncid,'v3d_time',time_i,d_v3d))
CALL check(nf90_DEF_dim(ncid,'zeta_time',time_ih,d_zeta))
CALL check(nf90_DEF_dim(ncid,'one',l_one,d_one))

! Crea variables i atributs

CALL check(nf90_def_var(ncid,'temp_time',NF90_DOUBLE,d_temp,varid))
CALL check(nf90_put_att(ncid,varid,"long_name","time for temperature climatology"))
CALL check(nf90_put_att(ncid,varid,"units","day"))
CALL check(nf90_put_att(ncid,varid,"cycle_leng ",cycle))

CALL check(nf90_def_var(ncid,'salt_time',NF90_DOUBLE,d_salt,varid))
CALL check(nf90_put_att(ncid,varid,"long_name","time for salinity climatology"))
CALL check(nf90_put_att(ncid,varid,"units","day"))
CALL check(nf90_put_att(ncid,varid,"cycle_leng ",cycle))

CALL check(nf90_def_var(ncid,'v2d_time',NF90_DOUBLE,d_v2d,varid))
CALL check(nf90_put_att(ncid,varid,"long_name","time for 2D velocity climatology"))
CALL check(nf90_put_att(ncid,varid,"units","day"))
CALL check(nf90_put_att(ncid,varid,"cycle_leng ",cycle))

CALL check(nf90_def_var(ncid,'v3d_time',NF90_DOUBLE,d_v3d,varid))
CALL check(nf90_put_att(ncid,varid,"long_name","time for 3D velocity climatology"))
CALL check(nf90_put_att(ncid,varid,"units","day"))
CALL check(nf90_put_att(ncid,varid,"cycle_leng ",cycle))

CALL check(nf90_def_var(ncid,'zeta_time',NF90_DOUBLE,d_zeta,varid))
CALL check(nf90_put_att(ncid,varid,"long_name","time for sea surface height"))
CALL check(nf90_put_att(ncid,varid,"units","day"))
CALL check(nf90_put_att(ncid,varid,"cycle_leng ",cycle))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CALL  check(nf90_def_var(ncid,'u_north',NF90_DOUBLE,(/d_xiu,d_srho,d_v3d/),varid))
CALL check(nf90_put_att(ncid,varid,"long_name","u-momentum component"))
CALL check(nf90_put_att(ncid,varid,"units","meter second-1"))
CALL  check(nf90_def_var(ncid,'u_south',NF90_DOUBLE,(/d_xiu,d_srho,d_v3d/),varid))
CALL check(nf90_put_att(ncid,varid,"long_name","u-momentum component"))
CALL check(nf90_put_att(ncid,varid,"units","meter second-1"))
CALL  check(nf90_def_var(ncid,'u_east',NF90_DOUBLE,(/d_etau,d_srho,d_v3d/),varid))
CALL check(nf90_put_att(ncid,varid,"long_name","u-momentum component"))
CALL check(nf90_put_att(ncid,varid,"units","meter second-1"))
CALL  check(nf90_def_var(ncid,'u_west',NF90_DOUBLE,(/d_etau,d_srho,d_v3d/),varid))
CALL check(nf90_put_att(ncid,varid,"long_name","u-momentum component"))
CALL check(nf90_put_att(ncid,varid,"units","meter second-1"))

CALL  check(nf90_def_var(ncid,'v_north',NF90_DOUBLE,(/d_xiv,d_srho,d_v3d/),varid))
CALL check(nf90_put_att(ncid,varid,"long_name","v-momentum component"))
CALL check(nf90_put_att(ncid,varid,"units","meter second-1"))
CALL  check(nf90_def_var(ncid,'v_south',NF90_DOUBLE,(/d_xiv,d_srho,d_v3d/),varid))
CALL check(nf90_put_att(ncid,varid,"long_name","v-momentum component"))
CALL check(nf90_put_att(ncid,varid,"units","meter second-1"))
CALL  check(nf90_def_var(ncid,'v_east',NF90_DOUBLE,(/d_etav,d_srho,d_v3d/),varid))
CALL check(nf90_put_att(ncid,varid,"long_name","v-momentum component"))
CALL check(nf90_put_att(ncid,varid,"units","meter second-1"))
CALL  check(nf90_def_var(ncid,'v_west',NF90_DOUBLE,(/d_etav,d_srho,d_v3d/),varid))
CALL check(nf90_put_att(ncid,varid,"long_name","v-momentum component"))
CALL check(nf90_put_att(ncid,varid,"units","meter second-1"))
i=1
CALL  check(nf90_def_var(ncid,'ubar_north',NF90_DOUBLE,(/d_xiu,d_v2d/),varid))
CALL check(nf90_put_att(ncid,varid,"long_name","verticaly integrated u-momentum component"))
CALL check(nf90_put_att(ncid,varid,"units","meter second-1"))
CALL  check(nf90_def_var(ncid,'ubar_south',NF90_DOUBLE,(/d_xiu,d_v2d/),varid))
CALL check(nf90_put_att(ncid,varid,"long_name","verticaly integrated u-momentum component"))
CALL check(nf90_put_att(ncid,varid,"units","meter second-1"))
CALL  check(nf90_def_var(ncid,'ubar_east',NF90_DOUBLE,(/d_etau,d_v2d/),varid))
CALL check(nf90_put_att(ncid,varid,"long_name","verticaly integrated u-momentum component"))
CALL check(nf90_put_att(ncid,varid,"units","meter second-1"))
CALL  check(nf90_def_var(ncid,'ubar_west',NF90_DOUBLE,(/d_etau,d_v2d/),varid))
CALL check(nf90_put_att(ncid,varid,"long_name","verticaly integrated u-momentum component"))
CALL check(nf90_put_att(ncid,varid,"units","meter second-1"))

CALL  check(nf90_def_var(ncid,'vbar_north',NF90_DOUBLE,(/d_xiv,d_v2d/),varid))
CALL check(nf90_put_att(ncid,varid,"long_name","verticaly integrated v-momentum component"))
CALL check(nf90_put_att(ncid,varid,"units","meter second-1"))
CALL  check(nf90_def_var(ncid,'vbar_south',NF90_DOUBLE,(/d_xiv,d_v2d/),varid))
CALL check(nf90_put_att(ncid,varid,"long_name","verticaly integrated v-momentum component"))
CALL check(nf90_put_att(ncid,varid,"units","meter second-1"))
CALL  check(nf90_def_var(ncid,'vbar_east',NF90_DOUBLE,(/d_etav,d_v2d/),varid))
CALL check(nf90_put_att(ncid,varid,"long_name","verticaly integrated v-momentum component"))
CALL check(nf90_put_att(ncid,varid,"units","meter second-1"))
CALL  check(nf90_def_var(ncid,'vbar_west',NF90_DOUBLE,(/d_etav,d_v2d/),varid))
CALL check(nf90_put_att(ncid,varid,"long_name","verticaly integrated v-momentum component"))
CALL check(nf90_put_att(ncid,varid,"units","meter second-1"))

CALL  check(nf90_def_var(ncid,'zeta_north',NF90_DOUBLE,(/d_xirho,d_zeta/),varid))
CALL check(nf90_put_att(ncid,varid,"long_name","free-surface"))
CALL check(nf90_put_att(ncid,varid,"units","meter"))
CALL  check(nf90_def_var(ncid,'zeta_south',NF90_DOUBLE,(/d_xirho,d_zeta/),varid))
CALL check(nf90_put_att(ncid,varid,"long_name","free-surface"))
CALL check(nf90_put_att(ncid,varid,"units","meter"))
CALL  check(nf90_def_var(ncid,'zeta_east',NF90_DOUBLE,(/d_etarho,d_zeta/),varid))
CALL check(nf90_put_att(ncid,varid,"long_name","free-surface"))
CALL check(nf90_put_att(ncid,varid,"units","meter"))
CALL  check(nf90_def_var(ncid,'zeta_west',NF90_DOUBLE,(/d_etarho,d_zeta/),varid))
CALL check(nf90_put_att(ncid,varid,"long_name","free-surface"))
CALL check(nf90_put_att(ncid,varid,"units","meter"))


CALL  check(nf90_def_var(ncid,'temp_north',NF90_DOUBLE,(/d_xirho,d_srho,d_temp/),varid))
CALL check(nf90_put_att(ncid,varid,"long_name","potencial temperature"))
CALL check(nf90_put_att(ncid,varid,"units","Celsius"))
CALL  check(nf90_def_var(ncid,'temp_south',NF90_DOUBLE,(/d_xirho,d_srho,d_temp/),varid))
CALL check(nf90_put_att(ncid,varid,"long_name","potencial temperature"))
CALL check(nf90_put_att(ncid,varid,"units","Celsius"))
CALL  check(nf90_def_var(ncid,'temp_east',NF90_DOUBLE,(/d_etarho,d_srho,d_temp/),varid))
CALL check(nf90_put_att(ncid,varid,"long_name","potencial temperature"))
CALL check(nf90_put_att(ncid,varid,"units","Celsius"))
CALL  check(nf90_def_var(ncid,'temp_west',NF90_DOUBLE,(/d_etarho,d_srho,d_temp/),varid))
CALL check(nf90_put_att(ncid,varid,"long_name","potencial temperature"))
CALL check(nf90_put_att(ncid,varid,"units","Celsius"))

CALL  check(nf90_def_var(ncid,'salt_north',NF90_DOUBLE,(/d_xirho,d_srho,d_salt/),varid))
CALL check(nf90_put_att(ncid,varid,"long_name","salinity"))
CALL check(nf90_put_att(ncid,varid,"units","PSU"))
CALL  check(nf90_def_var(ncid,'salt_south',NF90_DOUBLE,(/d_xirho,d_srho,d_salt/),varid))
CALL check(nf90_put_att(ncid,varid,"long_name","salinity"))
CALL check(nf90_put_att(ncid,varid,"units","PSU"))
CALL  check(nf90_def_var(ncid,'salt_east',NF90_DOUBLE,(/d_etarho,d_srho,d_salt/),varid))
CALL check(nf90_put_att(ncid,varid,"long_name","salinity"))
CALL check(nf90_put_att(ncid,varid,"units","PSU"))
CALL  check(nf90_def_var(ncid,'salt_west',NF90_DOUBLE,(/d_etarho,d_srho,d_salt/),varid))
CALL check(nf90_put_att(ncid,varid,"long_name","salinity"))
CALL check(nf90_put_att(ncid,varid,"units","PSU"))



CALL check(nf90_enddef(ncid))
END SUBROUTINE  create_bry
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE  create_ini(nom_arx)
USE variables
USE mod_kinds
USE netcdf
implicit none
character(*),intent(in) :: nom_arx
integer :: time_i
integer:: d_temp,d_salt,d_uclm,d_vclm,d_v2d,d_v3d,d_ssh,d_zeta
integer::ncid,d_xiu,d_etau,d_xiv,d_xirho,d_etarho,d_srho,d_sw,d_etav,l_one,varid,d_one,d_time
time_i=1
l_one=1
CALL check(nf90_create(TRIM(nom_arx),NF90_CLOBBER,ncid))
CALL check(nf90_DEF_dim(ncid,'xi_u',Nx_o-1,d_xiu))
CALL check(nf90_DEF_dim(ncid,'eta_u',Ny_o,d_etau))
CALL check(nf90_DEF_dim(ncid,'xi_v',Nx_o,d_xiv))
CALL check(nf90_DEF_dim(ncid,'eta_v',Ny_o-1,d_etav))
CALL check(nf90_DEF_dim(ncid,'xi_rho',Nx_o,d_xirho))
CALL check(nf90_DEF_dim(ncid,'eta_rho',Ny_o,d_etarho))
CALL check(nf90_DEF_dim(ncid,'s_rho',Niv_o,d_srho))
CALL check(nf90_DEF_dim(ncid,'s_w',Niv_o+1,d_sw))
CALL check(nf90_DEF_dim(ncid,'ocean_time',time_i,d_time))
CALL check(nf90_DEF_dim(ncid,'one',l_one,d_one))

CALL check(nf90_def_var(ncid,'theta_s',NF90_DOUBLE,d_one,varid))
CALL check(nf90_put_att(ncid,varid,"long_name","S-coordinate surface control parameter"))
CALL check(nf90_put_att(ncid,varid,"units","nondimensional"))
CALL check(nf90_def_var(ncid,'theta_b',NF90_DOUBLE,d_one,varid))
CALL check(nf90_put_att(ncid,varid,"long_name","S-coordinate bottom control parameter"))
CALL check(nf90_put_att(ncid,varid,"units","nondimensional"))
 
CALL check(nf90_def_var(ncid,'Tcline',NF90_DOUBLE,d_one,varid))

CALL check(nf90_put_att(ncid,varid,"long_name","S-coordinate surface/bottom layer width"))
CALL check(nf90_put_att(ncid,varid,"units","meter"))

CALL check(nf90_def_var(ncid,'hc',NF90_DOUBLE,d_one,varid))
CALL check(nf90_put_att(ncid,varid,"long_name","S-coordinate parameter,critical depth"))
CALL check(nf90_put_att(ncid,varid,"units","meter"))

CALL check(nf90_def_var(ncid,'sc_r',NF90_DOUBLE,d_srho,varid))
CALL check(nf90_put_att(ncid,varid,"long_name","S-coordinate at RHO-points"))
CALL check(nf90_put_att(ncid,varid,"units","nondimensional"))
CALL check(nf90_put_att(ncid,varid,"valid_min","-1"))
CALL check(nf90_put_att(ncid,varid,"valid_min","0"))

CALL check(nf90_def_var(ncid,'Cs_r',NF90_DOUBLE,d_srho,varid))
CALL check(nf90_put_att(ncid,varid,"long_name","S-coordinate stretching curves  at RHO-points"))
CALL check(nf90_put_att(ncid,varid,"units","nondimensional"))
CALL check(nf90_put_att(ncid,varid,"valid_min","-1"))
CALL check(nf90_put_att(ncid,varid,"valid_min","0"))
CALL check(nf90_def_var(ncid,'ocean_time',NF90_DOUBLE,d_time,varid)) !
CALL check(nf90_put_att(ncid,varid,"long_name","time since initialization"))
CALL check(nf90_put_att(ncid,varid,"units","second"))

CALL  check(nf90_def_var(ncid,'u',NF90_DOUBLE,(/d_xiu,d_etau,d_srho,d_time/),varid))
CALL check(nf90_put_att(ncid,varid,"long_name","u-momentum component"))
CALL check(nf90_put_att(ncid,varid,"units","meter second-1"))

CALL  check(nf90_def_var(ncid,'v',NF90_DOUBLE,(/d_xiv,d_etav,d_srho,d_time/),varid))
CALL check(nf90_put_att(ncid,varid,"long_name","v-momentum component"))
CALL check(nf90_put_att(ncid,varid,"units","meter second-1"))

CALL  check(nf90_def_var(ncid,'ubar',NF90_DOUBLE,(/d_xiu,d_etau,d_time/),varid))
CALL check(nf90_put_att(ncid,varid,"long_name","verticaly integrated u-momentum component"))
CALL check(nf90_put_att(ncid,varid,"units","meter second-1"))

CALL  check(nf90_def_var(ncid,'vbar',NF90_DOUBLE,(/d_xiv,d_etav,d_time/),varid))
CALL check(nf90_put_att(ncid,varid,"long_name","verticaly integrated v-momentum component"))
CALL check(nf90_put_att(ncid,varid,"units","meter second-1"))

CALL  check(nf90_def_var(ncid,'zeta',NF90_DOUBLE,(/d_xirho,d_etarho,d_time/),varid))
CALL check(nf90_put_att(ncid,varid,"long_name","free-surface"))
CALL check(nf90_put_att(ncid,varid,"units","meter"))


CALL  check(nf90_def_var(ncid,'temp',NF90_DOUBLE,(/d_xirho,d_etarho,d_srho,d_time/),varid))
CALL check(nf90_put_att(ncid,varid,"long_name","potencial temperature"))
CALL check(nf90_put_att(ncid,varid,"units","Celsius"))

CALL  check(nf90_def_var(ncid,'salt',NF90_DOUBLE,(/d_xirho,d_etarho,d_srho,d_time/),varid))
CALL check(nf90_put_att(ncid,varid,"long_name","salinity"))
CALL check(nf90_put_att(ncid,varid,"units","PSU"))

CALL check(nf90_enddef(ncid))
CALL check(nf90_inq_varid(ncid,"theta_b",varid))
CALL check(nf90_put_var(ncid,varid,theta_b_o))
CALL check(nf90_inq_varid(ncid,"theta_s",varid))
CALL check(nf90_put_var(ncid,varid,theta_s_o))
CALL check(nf90_inq_varid(ncid,"Tcline",varid))
CALL check(nf90_put_var(ncid,varid,Tcline_o))
CALL check(nf90_inq_varid(ncid,"hc",varid))
CALL check(nf90_put_var(ncid,varid,hc_o))
CALL check(nf90_inq_varid(ncid,"sc_r",varid))
CALL check(nf90_put_var(ncid,varid,sc_r_o))
CALL check(nf90_inq_varid(ncid,"Cs_r",varid))
CALL check(nf90_put_var(ncid,varid,Cs_r_o))
CALL check(nf90_inq_varid(ncid,"ocean_time",varid))
CALL check(nf90_put_var(ncid,varid,ocean_time))

CALL check(nf90_close(ncid))
END SUBROUTINE  create_ini
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine fill2d(L,M,tim,array)

  use variables
use mod_kinds
implicit none
integer,  intent(in):: L, M, tim
integer :: ii, jj, icounter,t,i,j,k,n,flag
real(r8), intent(inout) :: array(L,M,tim)
real(r8) :: dummy(L+2,M+2)
real(r8) ::  val
do t=1,tim
   flag=0      ! al moment q detectem un element no fillValue o zero el posem a 1
   DO i=1, L
      DO j=1,M
         if((array(i,j,t)/= f_v) .AND. (array(i,j,t) /= 0.0)  ) THEN
            flag=1
            exit
         endif
      ENDDO
      if (flag ==1) exit
      enddo
      if (flag ==0) THEN
         PRINT *," Mira les dades de la variable que estas fent"
         print *, "peta a fill2d"
         print *, "temps = ",t   
         stop
       endif
  enddo

do t=1,tim
  dummy=f_v
!c   main loop
  k=1000
  do while (k>0)
    dummy(2:L+1,2:M+1)=array(:,:,t)
    do i=2,L+1
      do j=2,M+1
        if (dummy(i,j).gt.f_v-1001) then  !lets fill a dry point
          icounter=0
          val=0.
          do ii=i-1,i+1
            do jj=j-1,j+1
              if (dummy(ii,jj).lt.f_v-1000) then
                icounter=icounter+1
                val=val+dummy(ii,jj)       
      	      endif
            enddo
          enddo
          if (icounter.gt.0) array(i-1,j-1,t)=val/real(icounter)  !point filled
        endif
      enddo
    enddo
    k=0
    do i=1,L
      do j=1,M
        if (array(i,j,t) > f_v-100) k=k+1  ! si hi ha terres k/= 0
      enddo
    enddo
  enddo
enddo !para cada paso de tiempo
end subroutine fill2d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!
subroutine fill3d(L,M,Nv,tim,array)
use variables
use mod_kinds
implicit none
integer,  intent(in):: L, M,Nv, tim
integer :: ii, jj, icounter,t,i,j,k,n,flag
real(r8), intent(inout) :: array(L,M,Nv,tim)
real(r8) :: dummy(L+2,M+2)
real(r8) ::  val
!!! testegem les dades al entrar es possible de o sihui 0 o fillvalue(errors detectats amb ibi)
do t=1,tim
  DO n=1,Nv
   flag=0      ! al moment q detectem un element no fillValue o zero el posem a 1
   DO i=1, L
      DO j=1,M
         if((array(i,j,n,t) .lt. f_v-1000) .AND. (array(i,j,n,t) /= 0.0)  ) THEN
            flag=1
            exit
         endif
      ENDDO
      if (flag ==1) exit
   enddo
      if (flag ==0) THEN
         PRINT *," Warning level without dates uncoment line 509"
         
        ! print *, "peta a fill3d"
      !   print *, "temps = ",t,"    nivell =  n",n
        !  stop            !   tot dewscomentat ferem que agafin els
         !  de la capa anterior
         array(:,:,n,t)=array(:,:,n-1,t)
       endif     
  enddo  !Niv
enddo   !time

do t=1,tim
!print *,"t =",t
  dummy=f_v
!c   main loop
DO n=1,Nv

  k=1000
  do while (k>0)
    dummy(2:L+1,2:M+1)=array(:,:,n,t)
    do i=2,L+1
      do j=2,M+1
        if (dummy(i,j).gt.f_v-10001) then  !lets fill a dry point
          icounter=0
          val=0.0_r8
          do ii=i-1,i+1
            do jj=j-1,j+1
              if (dummy(ii,jj).lt.f_v-1000) then
                icounter=icounter+1
                val=val+dummy(ii,jj)
!print *,i,j,n,t,val
              endif
            enddo
          enddo
          if (icounter.gt.0) array(i-1,j-1,n,t)=val/real(icounter)  !point filled
        endif
      enddo
    enddo
    k=0
    do i=1,L
      do j=1,M
       ! if (array(i,j,n,t) < f_v+1) k=k+1
	if (array(i,j,n,t) .gt. f_v-1000 ) k=k+1
      enddo
    enddo
  enddo
enddo !para cada paso de tiempo
enddo
end subroutine fill3d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine cnst2sigma(L,M,varc,vars,igrid)
! Es per a valors constant pero no regulars typus Myocean
! Hia una variavle global depth_i de length Nv_i 
! que ens dona les profunditats positives i nivell 1 a 1
USE variables
USE mod_kinds
implicit none
real(r8), intent(out):: vars(1:L,1:M,1:Niv_o,1:l_time_ih)
real(r8), intent(in):: varc(1:L,1:M,1:Niv_i,1:l_time_ih)
character (LEN=1), intent(in) :: igrid
integer, intent(in) :: L,M
real(r8) :: z_x(1:L,1:M,1:Niv_o)
integer ::i,j,n,k,t
IF (igrid=='r') THEN
   z_x=z_r_o
ELSE IF (igrid=='u') THEN
   z_x=z_u_o
ELSE IF (igrid=='v') THEN
   z_x=z_v_o
ELSE 
   print * ,"Has d'entrar un yipus de lmalla r,u, o v"
   STOP
END IF
   DO i=1,L
      DO j=1,M
         DO n=1,Niv_o
            IF (-z_x(i,j,n)<=depth_i(1)) THEN       
               DO t=1,l_time_ih
                vars(i,j,n,t) = varc(i,j,1,t)
               END DO
            ELSEIF (-z_x(i,j,n)>=depth_i(Niv_i)) THEN    !  prof in kes  h_k(   
               DO t=1,l_time_ih
                vars(i,j,n,t) =varc(i,j,Niv_i,t)
               END DO
            ELSE 
                  DO k=1,Niv_i-1
                     IF ((-z_x(i,j,n)>depth_i(k) ).and.( -z_x(i,j,n)<depth_i(k+1)  )) THEN   
 !interpolarem entre n i n-1
                        DO t=1,l_time_ih
                         vars(i,j,n,t)=linial(depth_i(k),varc(i,j,k,t),  & 
    &                                 depth_i(k+1),varc(i,j,k+1,t),-z_x(i,j,n))
                        END DO
                        EXIT
                     ENDIF
                  END DO
            ENDIF
         END DO
      END DO
  END DO

end subroutine cnst2sigma
!!!!!!!!!!!!!!!!!!!!

SUBROUTINE  make_bar(nc_i,nc_o)
! Entra un nc_i amb  u,v,hz 
! Surt  un nc_0 amb ubar,vbar
USE mod_kinds
USE variables
USE netcdf
implicit none
character (LEN=120), intent(in) :: nc_i
character (LEN=120), intent(in) :: nc_o
!real(r8), intent(int):: vars(1:L,1:M,1:Niv_o,1:l_time_ih)
!CALL check(nf90_open(TRIM(nc_i),NF90_NOWRITE,ncidg))  
print *, nc_i,"collob",nc_o


END  SUBROUTINE make_bar
!!!!!!!!!!!!!!!!!!!!
SUBROUTINE  rho2v(L,M,N,T,varr,varv)

! Passa una matriu  de punt U a rho  
! L, M, N, dimensions de la matriu rho
USE mod_kinds
USE variables
implicit none
integer, intent(in) ::  L, M, N,T
real(r8), intent(in)::varr(1:L,1:M,1:N,1:T)
real(r8), intent(out)::varv(1:L,1:M-1,1:N,1:T)
integer :: i
DO i=1,M-1
   varv(:,i,:,:)=0.5_r8*(varr(:,i,:,:)+varr(:,i+1,:,:))
ENDDO
END  SUBROUTINE rho2v
!!!!!!!!!!!!!!!!!!!!!!!1
SUBROUTINE  rho2u(L,M,N,T,varr,varu)

! Passa una matriu  de punt U a rho  
! L, M, N, T dimensions de la matriu rho
USE mod_kinds
USE variables
implicit none
integer, intent(in) ::  L, M, N, T
real(r8), intent(in)::varr(1:L,1:M,1:N,1:T)
real(r8), intent(out)::varu(1:L-1,1:M,1:N,1:T)
integer :: i
DO i=1,L-1
   varu(i,:,:,:)=0.5_r8*(varr(i,:,:,:)+varr(i+1,:,:,:))
ENDDO
END  SUBROUTINE rho2u
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
SUBROUTINE  linterp2d(Linp,Minp,Xinp,Yinp,Finp,Lout,Mout,Xout,Yout,Iout,Jout,Fout,Fmin,Fmax)
!
! Donat un camp  Finp situat en (Xinp,Yinp) in es vol interpolar en un
! grid (Xout,Yout) sabent les fractionals grids cells index Iout,Jout
! Fets am lña subroutine hindices  es troba el camp interpolat Fout
! tambe dona els valors Fmin i Fmax per informacio
!
USE mod_kinds
implicit none
integer, intent(in) :: Linp, Minp, Lout,Mout
real(r8), intent(in) :: Xinp(1:Linp,1:Minp)
real(r8), intent(in) :: Yinp(1:Linp,1:Minp)
real(r8), intent(in) :: Finp(1:Linp,1:Minp)

real(r8), intent(in) :: Xout(1:Lout,1:Mout)
real(r8), intent(in) :: Yout(1:Lout,1:Mout)
!real(r8), intent(in) :: Amask(1:Lout,1:Mout)
real(r8), intent(in) :: Iout(1:Lout,1:Mout)
real(r8), intent(in) :: Jout(1:Lout,1:Mout)

real(r8), intent(out) :: Fout(1:Lout,1:Mout)
real(r8), intent(out) :: Fmin, Fmax

integer :: i, i1, i2, j, j1, j2
real(r8) :: p1, p2, q1, q2
   Fmin=1.0E+35_r8
      Fmax=-1.0E+35_r8
      DO j=1,Mout
        DO i=1,Lout
         ! IF (Amask(i,j)==0) THENT EXIT
          i1=INT(Iout(i,j))
          i2=i1+1
          j1=INT(Jout(i,j))
          j2=j1+1
            p2=REAL(i2-i1,r8)*(Iout(i,j)-REAL(i1,r8))
            q2=REAL(j2-j1,r8)*(Jout(i,j)-REAL(j1,r8))
            p1=1.0_r8-p2
            q1=1.0_r8-q2
            Fout(i,j)=p1*q1*Finp(i1,j1)+                                &
     &                p2*q1*Finp(i2,j1)+                                &
     &                p2*q2*Finp(i2,j2)+                                &
     &                p1*q2*Finp(i1,j2)
            Fmin=MIN(Fmin,Fout(i,j))
            Fmax=MAX(Fmax,Fout(i,j))
        END DO
      END DO

END SUBROUTINE  linterp2d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!! 
 SUBROUTINE hindices(L,M,Xgrd,Ygrd,Lo,Mo,lono,lato,angler,Ipos,Jpos)
!  Ipos index donat per try_tange amb el decimal proporcionat a la posicio
!  L,M   dims del punts de la malla gran i
!  Xgrd  , Ygrd   longi. i lats de la malla 
!  angler  angle entre eix X i east  (rads)
!  Lo,Mo,  dimesions de la malla out
!  lono ,lato  lonfituds i latitus dels punts de la malla de valors de camp desconeguts
!  Ipos,Jpos els valors proporcionsls per a cada punt de la malla desconeguda per presuenta interpolacio 

use mod_kinds
use variables
implicit none
integer, intent(in)   :: L, M, Lo, Mo
REAL(r8), intent(in)  :: Xgrd(1:L,1:M)
REAL(r8), intent(in)  :: Ygrd(1:L,1:M)
REAL(r8), intent(in)  :: angler(1:L,1:M)
REAL(r8), intent(in)  :: lono(1:Lo,1:Mo)
REAL(r8), intent(in)  :: lato(1:Lo,1:Mo)
REAL(r8), intent(out) :: Ipos(1:Lo,1:Mo)
REAL(r8), intent(out) :: Jpos(1:Lo,1:Mo)
 
REAL(r8) :: xfac, yfac, xpp, ypp, dx, dy, ang, diag2, aa2, bb2, phi, p1,p2,q1,q2    !variables per la interpolacio
INTEGER :: Imin,Jmin,i,j,ppp,qqq         
do i=1,Lo
   do  j=1, Mo
  CALL try_cel(L,M,Xgrd,Ygrd,lono(i,j),lato(i,j),Imin,Jmin)     
      yfac=Eradius*deg2rad             !convertim les distancies en metres estem en malla espherica
      xfac=yfac*COS(lato(i,j)*deg2rad)
      xpp=(lono(i,j)-Xgrd(Imin,Jmin))*xfac
      ypp=(lato(i,j)-Ygrd(Imin,Jmin))*yfac
!  Use Law of Cosines to get cell parallelogram "shear" angle.
           diag2=((Xgrd(Imin+1,Jmin)-Xgrd(Imin,Jmin+1))*xfac)**2+      &
    &            ((Ygrd(Imin+1,Jmin)-Ygrd(Imin,Jmin+1))*yfac)**2
           aa2=((Xgrd(Imin,Jmin)-Xgrd(Imin+1,Jmin))*xfac)**2+          &
    &          ((Ygrd(Imin,Jmin)-Ygrd(Imin+1,Jmin))*yfac)**2
           bb2=((Xgrd(Imin,Jmin)-Xgrd(Imin,Jmin+1))*xfac)**2+          &
    &          ((Ygrd(Imin,Jmin)-Ygrd(Imin,Jmin+1))*yfac)**2
           phi=ASIN((diag2-aa2-bb2)/(2.0_r8*SQRT(aa2*bb2)))

!  Transform float position into curvilinear coordinates. Assume the
!  cell is rectanglar, for now.

           ang=angler(Imin,Jmin)
           dx=xpp*COS(ang)+ypp*SIN(ang)
           dy=ypp*COS(ang)-xpp*SIN(ang)
!
!  Correct for parallelogram.

           dx=dx+dy*TAN(phi)
           dy=dy/COS(phi)

!  Scale with cell side lengths to translate into cell indices.
!
            dx=MIN(MAX(0.0_r8,dx/SQRT(aa2)),1.0_r8)
            dy=MIN(MAX(0.0_r8,dy/SQRT(bb2)),1.0_r8)
 
           Ipos(i,j)=REAL(Imin,r8)+dx
           Jpos(i,j)=REAL(Jmin,r8)+dy

   END DO
ENDDO

END SUBROUTINE hindices
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE  try_cel(L,M,Xgrd,Ygrd,xo,yo,Ie,Je)
use variables
USE mod_kinds
implicit none
!xo,yo punt explorador a situar 
!  L, M dimensio x,y ,de lon i lat de la malla gran 
!  Ie, Je, index del vextex inferior esquerra del rectangle  de la malla in que conte el punt xo,yo   
integer, intent(in) :: L, M
REAL(r8),intent(in) :: xo,yo
REAL(r8),intent(in) ::Xgrd(1:L,1:M) 
REAL(r8),intent(in) ::Ygrd(1:L,1:M)
integer, intent(out) :: Ie,Je
integer :: Imin, Imax, Jmin, Jmax,i0,j0
logical :: dins, found
REAL(r8),dimension(4) :: xb, yb
Imin=1
Imax=L
Jmin=1
Jmax=M
      do 
        IF (.not.(((Imax-Imin).gt.1).or.((Jmax-Jmin).gt.1))) EXIT  !voll dir feina feta
         IF ((Imax-Imin).gt.1) THEN
                   i0=(Imin+Imax)/2
                   xb=(/Xgrd(Imin,Jmin),Xgrd(i0,Jmin),Xgrd(i0,Jmax),Xgrd(Imin,Jmax)/)
                   yb=(/Ygrd(Imin,Jmin),Ygrd(i0,Jmin),Ygrd(i0,Jmax),Ygrd(Imin,Jmax)/)
                   found= inside_4v(xb,yb,xo,yo)
                   IF (found) THEN
                      Imax=i0
                   ELSE
                      Imin=i0
                   ENDIF
          ENDIF
          IF  ((Jmax-Jmin).gt.1) THEN
                   j0=(Jmin+Jmax)/2
                   xb=(/Xgrd(Imin,Jmin),Xgrd(Imax,Jmin),Xgrd(Imax,j0),Xgrd(Imin,j0)/)
                   yb=(/Ygrd(Imin,Jmin),Ygrd(Imax,Jmin),Ygrd(Imax,j0),Ygrd(Imin,j0)/)  
                   found= inside_4v(xb,yb,xo,yo)
                   IF (found) THEN
                      Jmax=j0
                   ELSE
                      Jmin=j0
                 ENDIF

              ENDIF                  
          END DO                 
Ie=Imin
Je=Jmin
END SUBROUTINE try_cel

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      LOGICAL FUNCTION inside_4v ( X, Y, Xo, Yo)
! X=(/x1, x2, x3, x4/)   ! abcsices dels quatre vertex
! Y=(/y1, y2, y3, y4/) ! ordenades dels quatre vertex
!Xo = abscissa del punt explorador
!Yo = ordenada del punt explorador
implicit none
!
!
!  Imported variable declarations.
!
      real(r8), intent(in) :: Xo, Yo
      real(r8),DIMENSION(4),intent(in) :: X, Y

!  Imported variable declarations.
!
      
!
!  Local variable declarations.
!
      integer, parameter :: Nstep =128
      integer :: crossings, i, inc, k, kk, nc
      integer, dimension(Nstep) :: Sindex
      real(r8) :: dx1, dx2, dxy
!  Local variable declarations.
!
      integer, parameter :: Nb=5
      real(r8),DIMENSION(5):: Xb,Yb
 
      Xb(1)=X(1)
      Xb(2)=X(2)
      Xb(3)=X(3)
      Xb(4)=X(4)
      Xb(5)=X(1)
      Yb(1)=Y(1)
      Yb(2)=Y(2)
      Yb(3)=Y(3)
      Yb(4)=Y(4)
      Yb(5)=Y(1)
!
!
!-----------------------------------------------------------------------
!  Find intersections.
!-----------------------------------------------------------------------
!
!  Set crossings counter and close the contour of the polygon.
!
      crossings=0
      Xb(Nb+1)=Xb(1)
      Yb(Nb+1)=Yb(1)
!
!  The search is optimized.  First select the indices of segments
!  where Xb(k) is different from Xb(k+1) and Xo falls between them.
!  Then, further investigate these segments in a separate loop.
!  Doing it in two stages takes less time because the first loop is
!  pipelined.
!
      DO kk=0,Nb-1,Nstep
        nc=0
        DO k=kk+1,MIN(kk+Nstep,Nb)
          IF (((Xb(k+1)-Xo)*(Xo-Xb(k)).ge.0.0_r8).and.                  &
     &        (Xb(k).ne.Xb(k+1))) THEN
            nc=nc+1
            Sindex(nc)=k
          END IF
        END DO
        DO i=1,nc
          k=Sindex(i)
          IF (Xb(k).ne.Xb(k+1)) THEN
            dx1=Xo-Xb(k)
            dx2=Xb(k+1)-Xo
            dxy=dx2*(Yo-Yb(k))-dx1*(Yb(k+1)-Yo)
            inc=0
            IF ((Xb(k).eq.Xo).and.(Yb(k).eq.Yo)) THEN
              crossings=1
              goto 10
            ELSE IF (((dx1.eq.0.0_r8).and.(Yo.ge.Yb(k  ))).or.          &
     &              ((dx2.eq.0.0_r8).and.(Yo.ge.Yb(k+1)))) THEN
              inc=1
            ELSE IF ((dx1*dx2.gt.0.0_r8).and.                           &
     &              ((Xb(k+1)-Xb(k))*dxy.ge.0.0_r8)) THEN  ! see note 1
              inc=2
            END IF
            IF (Xb(k+1).gt.Xb(k)) THEN
              crossings=crossings+inc
            ELSE
              crossings=crossings-inc
            END IF
          END IF
        END DO
      END DO
!
!  Determine if point (Xo,Yo) is inside of closed polygon.
!
  10  IF (crossings.eq.0) THEN
        inside_4v=.FALSE.
      ELSE
        inside_4v=.TRUE.
      END IF
      RETURN
 
      END FUNCTION inside_4v
!!!!!!!!!!!!!!!!!

END MODULE

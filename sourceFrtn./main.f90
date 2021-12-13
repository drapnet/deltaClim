!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!=======================================================================!
!    This file is part of SAMOA project                                 !
!    Engineered by Manel Grifoll (LIM/UPC)                              !
!                                                                       !
!=======================================================================!

! Exemple del que farem a cada variable (agafem  com exemple latemp):
! tempi:  (variable input del history) primer se li fa un fill3d sense canviar el nom 
! tempoc:  interpolant els valors de cada  nivell k alspunts de sortida
! tempo:   valors de temp als punts de sortida i nivells sigma.
! 	tempi->(fill3d )->(linterp2d)->tempoc->(cnst2sigma)->tempo
!       u,v requires rotation.
! ubari -> ubarr -> ubarg ->ubaro		         
! input    en rho   rot     out
program main
use netcdf
use mod_kinds
use funcions
use variables
implicit none
integer ::i,j,k,q,p,m,n,t
integer :: ncid, varid,ncidi,ncidc
real(r8),DIMENSION(:,:,:),ALLOCATABLE ::zetai,zetao,ubi,vbi,ubor,vbor, &
			&	ubo,vbo,ubari,ubaro,vbari,vbaro,aux3
real(r8),DIMENSION(:,:,:),ALLOCATABLE :: ubarg,vbarg,ubarr,vbarr
real(r8),DIMENSION(:,:,:,:),ALLOCATABLE ::tempi,salti,ui,vi,aux4  
real(r8),DIMENSION(:,:,:,:),ALLOCATABLE ::tempoc,saltoc,uoc,voc!,uos,vos
real(r8),DIMENSION(:,:,:,:),ALLOCATABLE ::tempo,salto,uo,vo,uog,vog,uof,vof

real(r8),DIMENSION(:,:),ALLOCATABLE ::Du,Dv,Duk,Dvk !calcul bars
real(r8),DIMENSION(:,:),ALLOCATABLE :: Ir,Jr 
real(r8),DIMENSION(:),ALLOCATABLE :: time_h  ! vector temps pel horaris
real(r8),DIMENSION(:),ALLOCATABLE :: time_q  ! vector temps pel quarts d hora
real(r8) :: val_min,val_max,rr,scalef,offset,scalef_s,offset_s,scalef_u,offset_u,time_d

CALL ini

print *, "PART 2: create clim"
print *, "./dat/aux1"
!print *,TRIM(dir_out)//TRIM(nom_bry)
CALL create_clm("./dat/aux1.nc",l_time_iq,l_time_ih)

CALL check(nf90_open("./dat/aux1.nc",NF90_WRITE,ncidc))            ! oberts el clim ncdic
ALLOCATE(Ir(1:Nx_o,1:Ny_o))
ALLOCATE(Jr(1:Nx_o,1:Ny_o))

CALL hindices(Nx_i,Ny_i,lon_i,lat_i,Nx_o,Ny_o,lon_r_o,lat_r_o,angle_i,Ir,Jr)
print *,"Interpolation indices done"

ALLOCATE(zetai(1:Nx_i,1:Ny_i,l_time_iq))
ALLOCATE(zetao(1:Nx_o,1:Ny_o,l_time_iq))


ALLOCATE(aux3(1:Nx_i,1:Ny_i,l_time_iq))

!!!CALL check(nf90_open(TRIM(dir_in)//"TRIM(arx_hm),NF90_WRITE,ncidi))
CALL check(nf90_open("./dat/ssh.nc",NF90_WRITE,ncidi))
CALL check(nf90_inq_varid(ncidi,'zos',varid))

!CALL check(nf90_get_att(ncidi,varid,"scale_factor",scalef))
!CALL check(nf90_get_att(ncidi,varid,"add_offset",offset))
CALL check(nf90_get_var(ncidi,varid,aux3))
zetai(:,:,1:l_time_iq)=aux3

print *,"Starting variables zeta"
call fill2d(Nx_i,Ny_i,l_time_iq,zetai)
print *, "Fet el filled de la zeta"
DO t=1,l_time_iq
       CALL linterp2d(Nx_i,Ny_i,lon_i,lat_i,zetai(:,:,t),Nx_o,Ny_o,lon_r_o,lat_r_o, &
     & Ir, Jr, zetao(:,:,t),val_min,val_max) 
END DO
CALL check(nf90_inq_varid(ncidc,"zeta",varid))
CALL check(nf90_put_var(ncidc,varid,zetao))
print *,"Zeta interpolation done"
!print *,zetao
CALL check(nf90_close(ncidi))
!ubari(:,:,1:24)=aux3
print *,"Starting variables"
ALLOCATE(tempi(1:Nx_i,1:Ny_i,1:Niv_i,l_time_ih))
ALLOCATE(salti(1:Nx_i,1:Ny_i,1:Niv_i,l_time_ih))
ALLOCATE(ui(1:Nx_i,1:Ny_i,1:Niv_i,l_time_ih))
ALLOCATE(vi(1:Nx_i,1:Ny_i,1:Niv_i,l_time_ih))

CALL check(nf90_open("./dat/temp.nc",NF90_WRITE,ncidi))
CALL check(nf90_inq_varid(ncidi,'thetao',varid))
CALL check(nf90_get_var(ncidi,varid,tempi))
CALL check(nf90_close(ncidi))

CALL check(nf90_open("./dat/salt.nc",NF90_WRITE,ncidi))
CALL check(nf90_inq_varid(ncidi,'so',varid))
CALL check(nf90_get_var(ncidi,varid,salti))
CALL check(nf90_close(ncidi))

CALL check(nf90_open("./dat/cur-hor.nc",NF90_WRITE,ncidi))
CALL check(nf90_inq_varid(ncidi,'uo',varid))
CALL check(nf90_get_var(ncidi,varid,ui))
CALL check(nf90_inq_varid(ncidi,'vo',varid))
CALL check(nf90_get_var(ncidi,varid,vi))
CALL check(nf90_close(ncidi))

print *, "llegit tot"

CALL fill3d(Nx_i,Ny_i,niv_i,l_time_ih,tempi)
print *,"fet el fillet de tot temp"
!print *, tempi
CALL fill3d(Nx_i,Ny_i,Niv_i,l_time_ih,salti)
print *,"fet el fillet de tot salt"
CALL fill3d(Nx_i,Ny_i,niv_i,l_time_ih,ui)
print *,"fet el fillet de tot u"
CALL fill3d(Nx_i,Ny_i,niv_i,l_time_ih,vi)
print *,"fet el fillet de tot"


ALLOCATE(tempoc(1:Nx_o,1:Ny_o,1:Niv_i,l_time_ih))
ALLOCATE(tempo(1:Nx_o,1:Ny_o,1:Niv_o,l_time_ih))
ALLOCATE(saltoc(1:Nx_o,1:Ny_o,1:Niv_i,l_time_ih))
ALLOCATE(salto(1:Nx_o,1:Ny_o,1:Niv_o,l_time_ih))
DO k=1,Niv_i
   do t=1,l_time_ih
       CALL linterp2d(Nx_i,Ny_i,lon_i,lat_i,tempi(:,:,k,t),Nx_o,Ny_o,lon_r_o,lat_r_o, &
      & Ir, Jr, tempoc(:,:,k,t),val_min,val_max) 
       CALL linterp2d(Nx_i,Ny_i,lon_i,lat_i,salti(:,:,k,t),Nx_o,Ny_o,lon_r_o,lat_r_o, &
       & Ir, Jr, saltoc(:,:,k,t),val_min,val_max) 
    END DO
ENDDO

CALL cnst2sigma(Nx_o,Ny_o,tempoc,tempo,'r')
CALL check(nf90_inq_varid(ncidc,"temp",varid))
CALL check(nf90_put_var(ncidc,varid,tempo))
!CALL check(nf90_inq_varid(ncid,"temp",varid))
!CALL check(nf90_put_var(ncid,varid,tempo(:,:,:,1)))
CALL cnst2sigma(Nx_o,Ny_o,saltoc,salto,'r')
CALL check(nf90_inq_varid(ncidc,"salt",varid))
CALL check(nf90_put_var(ncidc,varid,salto))
!CALL check(nf90_inq_varid(ncid,"salt",varid))
!CALL check(nf90_put_var(ncid,varid,salto(:,:,:,1)))
print *, "Salt and Temp done"

ALLOCATE(uoc(1:Nx_o,1:Ny_o,1:Niv_i,l_time_ih))   !Uoc  en punts rho i niv constant
ALLOCATE(uo(1:Nx_o,1:Ny_o,1:Niv_o,l_time_ih))    !Uo  puts rho nivells sigma
ALLOCATE(voc(1:Nx_o,1:Ny_o,1:Niv_i,l_time_ih))   !voc  en punts rho i niv constant
ALLOCATE(vo(1:Nx_o,1:Ny_o,1:Niv_o,l_time_ih))    !vo  puts rho nivells sigma
DO k=1,Niv_i
   do t=1,l_time_ih
       CALL linterp2d(Nx_i,Ny_i,lon_i,lat_i,ui(:,:,k,t),Nx_o,Ny_o,lon_r_o,lat_r_o, &
      & Ir, Jr, uoc(:,:,k,t),val_min,val_max) 
       CALL linterp2d(Nx_i,Ny_i,lon_i,lat_i,vi(:,:,k,t),Nx_o,Ny_o,lon_r_o,lat_r_o, &
      & Ir, Jr, voc(:,:,k,t),val_min,val_max) 
    END DO
 ENDDO
! print *,uoc(:,:,3,18)
CALL cnst2sigma(Nx_o,Ny_o,uoc,uo,'r')
ALLOCATE(uog(1:Nx_o,1:Ny_o,1:Niv_o,l_time_ih))    ! girem la uo i vo  de N-S  en els sentit de la malla  ROTEM
CALL cnst2sigma(Nx_o,Ny_o,voc,vo,'r')              !  uog j vog  girades   tot rn punts rho
ALLOCATE(vog(1:Nx_o,1:Ny_o,1:Niv_o,l_time_ih))   
DO i=1,Nx_o
   DO j=1,Ny_o
      uog(i,j,:,:)=uo(i,j,:,:)*cos(angle_o(i,j))+vo(i,j,:,:)*sin(angle_o(i,j))
      vog(i,j,:,:)=-uo(i,j,:,:)*sin(angle_o(i,j))+vo(i,j,:,:)*cos(angle_o(i,j))
   END DO
END DO
ALLOCATE(uof(1:Nx_o-1,1:Ny_o,1:Niv_o,l_time_ih))
ALLOCATE(vof(1:Nx_o,1:Ny_o-1,1:Niv_o,l_time_ih))
CALL rho2v(Nx_o,Ny_o,Niv_o,l_time_ih,vog,vof)
CALL rho2u(Nx_o,Ny_o,Niv_o,l_time_ih,uog,uof)
!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!
CALL check(nf90_inq_varid(ncidc,"u",varid))
CALL check(nf90_put_var(ncidc,varid,uof))
CALL check(nf90_inq_varid(ncidc,"v",varid))
CALL check(nf90_put_var(ncidc,varid,vof))

print *,"Water veolcities (u and v) done"
!Fem les ubar i vbar agafant les qm
CALL check(nf90_open("./dat/cur-qm.nc",NF90_WRITE,ncidi))
ALLOCATE(ubi(1:Nx_i,1:Ny_i,l_time_iq))  !mala in
ALLOCATE(vbi(1:Nx_i,1:Ny_i,l_time_iq))
ALLOCATE(ubor(1:Nx_o,1:Ny_o,l_time_iq))  !malla out punts rho
ALLOCATE(vbor(1:Nx_o,1:Ny_o,l_time_iq))

CALL check(nf90_inq_varid(ncidi,'uo',varid))
CALL check(nf90_get_var(ncidi,varid,ubi))
CALL check(nf90_inq_varid(ncidi,'vo',varid))
CALL check(nf90_get_var(ncidi,varid,vbi))
CALL check(nf90_close(ncidi))

call fill2d(Nx_i,Ny_i,l_time_iq,ubi)
print *, "Fet el filled de la ubi"
DO t=1,l_time_iq
       CALL linterp2d(Nx_i,Ny_i,lon_i,lat_i,ubi(:,:,t),Nx_o,Ny_o,lon_r_o,lat_r_o, &
     & Ir, Jr, ubor(:,:,t),val_min,val_max) 
END DO

call fill2d(Nx_i,Ny_i,l_time_iq,vbi)
print *, "Fet el filled de la vbi"

DO t=1,l_time_iq
       CALL linterp2d(Nx_i,Ny_i,lon_i,lat_i,ubi(:,:,t),Nx_o,Ny_o,lon_r_o,lat_r_o, &
     & Ir, Jr, vbor(:,:,t),val_min,val_max) 
END DO

ALLOCATE(ubo(1:Nx_o-1,1:Ny_o,l_time_iq))  !malla out punts rho
ALLOCATE(vbo(1:Nx_o,1:Ny_o-1,l_time_iq))

! Com es una malla rectangular (no cal girar la malla angle_o =90)
CALL rho2v(Nx_o,Ny_o-1,1,l_time_iq,vbor,vbo)
CALL rho2u(Nx_o-1,Ny_o,1,l_time_iq,ubor,ubo)

CALL check(nf90_inq_varid(ncidc,"vbar",varid))
CALL check(nf90_put_var(ncidc,varid,vbo))

CALL check(nf90_inq_varid(ncidc,"ubar",varid))
CALL check(nf90_put_var(ncidc,varid,ubo))








!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

ALLOCATE(time_h(1:l_time_ih))
time_h(1)=0
DO i=2,l_time_ih
   time_h(i)=time_h(i-1)+1.0_r8/24.0_r8
ENDDO
ALLOCATE(time_q(1:l_time_iq))
time_q(1)=0
DO i=2,l_time_iq
   time_q(i)=time_q(i-1)+0.25_r8/24.0_r8
ENDDO



! Esta obert el  clim
CALL check(nf90_inq_varid(ncidc,"salt_time",varid))
CALL check(nf90_put_var(ncidc,varid,time_h))
CALL check(nf90_inq_varid(ncidc,"temp_time",varid))
CALL check(nf90_put_var(ncidc,varid,time_h))
CALL check(nf90_inq_varid(ncidc,"v2d_time",varid))
CALL check(nf90_put_var(ncidc,varid,time_q))
CALL check(nf90_inq_varid(ncidc,"v3d_time",varid))
CALL check(nf90_put_var(ncidc,varid,time_d))
CALL check(nf90_inq_varid(ncidc,"zeta_time",varid))
CALL check(nf90_put_var(ncidc,varid,time_q))
CALL check(nf90_close(ncidc))    !tenquem el clm

print *, "FINISH: "

END PROGRAM

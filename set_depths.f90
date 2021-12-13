!=======================================================================!
!    This file is part of SAMOA project                                       !
!    Engineered by Manel Grifoll (LIM)                                  !
!                                                                       !
!=======================================================================!

SUBROUTINE set_depths
! Aquesta subroutina dona les profunditats dels nivells sigma, en rho,u,v
! Necesitem :
! La malla  ( amb CALL get_grid    tenim tot )
! Nivells verticals  N
! hc  H critica
! theta
! i son necesseris els vectors Cs_r,Cs_w,sc_r,sc_w, i els nivells rho

USE mod_kinds
USE variables
USE funcions
USE netcdf
implicit none
!character(*) , intent(in) ::nom_grd
!real(r8),dimension(2),intent(in) ::Zt_avg1
integer :: i, j, k,n
real(r8) :: cff_r, cff1_r, cff2_r, cff_w, cff1_w, cff2_w,cff,cff1,cff2
real(r8) :: hinv, hwater, z_r0, z_w0,ds


real(r8) :: hu_o(1:Nx_o-1,1:Ny_o),Zt_avg1u_o(1:Nx_o-1,1:Ny_o)
real(r8) :: hv_o(1:Nx_o,1:Ny_o-1),Zt_avg1v_o(1:Nx_o,1:Ny_o-1)

! Aqui farem wl Cs_r i el sc_r pel ini de la sortida ames deles depths

        IF (theta_s_o.ne.0.0_r8) THEN
          cff1=1.0_r8/SINH(theta_s_o)
          cff2=0.5_r8/TANH(0.5_r8*theta_s_o)
        END IF
        sc_w_o(0)=-1.0_r8
        Cs_w_o(0)=-1.0_r8
        ds=1.0_r8/REAL(Niv_o,r8)

        DO k=1,Niv_o
          sc_w_o(k)=ds*REAL(k-Niv_o,r8)
          sc_r_o(k)=ds*(REAL(k-Niv_o,r8)-0.5_r8)
          IF (theta_s_o.ne.0.0_r8) THEN
            Cs_w_o(k)=(1.0_r8-theta_b_o)*                   &
     &                          cff1*SINH(theta_s_o*                  &
     &                                    sc_w_o(k))+         &
     &                          theta_b_o*                            &
     &                          (cff2*TANH(theta_s_o*                 &
     &                                     (sc_w_o(k)+        &
     &                                      0.5_r8))-                   &
     &                           0.5_r8)
            Cs_r_o(k)=(1.0_r8-theta_b_o)*                   &
     &                          cff1*SINH(theta_s_o*                  &
     &                                    sc_r_o(k))+         &
     &                          theta_b_o*                            &
     &                          (cff2*TANH(theta_s_o*                 &
     &                                     (sc_r_o(k)+        &
     &                                      0.5_r8))-                   &
     &                           0.5_r8)
          ELSE
            Cs_w_o(k)=sc_w_o(k)
            Cs_r_o(k)=sc_r_o(k)
          END IF
        END DO

!!!!!!!!!!!!!

DO j=1,Ny_o
    DO i=1,Nx_o
       z_w_o(i,j,0)=-h_o(i,j)
    END DO
    DO k=1,Niv_o
        cff_r=hc_o*(sc_r_o(k)-Cs_r_o(k))
        cff_w=hc_o*(sc_w_o(k)-Cs_w_o(k))
        cff1_r=Cs_r_o(k)
        cff1_w=Cs_w_o(k)
            DO i=1,Nx_o
              hwater=h_o(i,j)
              hinv=1.0_r8/hwater
              z_w0=cff_w+cff1_w*hwater
              z_w_o(i,j,k)=z_w0+Zt_avg1_o(i,j)*(1.0_r8+z_w0*hinv)
              z_r0=cff_r+cff1_r*hwater
              z_r_o(i,j,k)=z_r0+Zt_avg1_o(i,j)*(1.0_r8+z_r0*hinv)
              Hz_o(i,j,k)=z_w_o(i,j,k)-z_w_o(i,j,k-1)
            END DO
     END DO
END DO
! punts u, h_u_o
do i=1,Nx_o-1
   hu_o(i,:)=0.5_r8*(h_o(i,:)+h_o(i+1,:))
   Zt_avg1u_o(i,:)=0.5_r8*(Zt_avg1_o(i,:)+Zt_avg1_o(i+1,:))
enddo

do j=1,Ny_o-1
   hv_o(:,j)=0.5_r8*(h_o(:,j)+h_o(:,j+1))
   Zt_avg1v_o(:,j)=0.5_r8*(Zt_avg1_o(:,j)+Zt_avg1_o(:,j+1))
enddo
DO j=1,Ny_o
    DO k=1,Niv_o
        cff_r=hc_o*(sc_r_o(k)-Cs_r_o(k))
        cff_w=hc_o*(sc_w_o(k)-Cs_w_o(k))
        cff1_r=Cs_r_o(k)
            DO i=1,Nx_o-1
              hwater=hu_o(i,j)
              hinv=1.0_r8/hwater
              z_r0=cff_r+cff1_r*hwater
              z_u_o(i,j,k)=z_r0+Zt_avg1u_o(i,j)*(1.0_r8+z_r0*hinv)
            END DO
     END DO
END DO
DO j=1,Ny_o-1
    DO k=1,Niv_o
        cff_r=hc_o*(sc_r_o(k)-Cs_r_o(k))
        cff_w=hc_o*(sc_w_o(k)-Cs_w_o(k))
        cff1_r=Cs_r_o(k)
            DO i=1,Nx_o
              hwater=hv_o(i,j)
              hinv=1.0_r8/hwater
              z_r0=cff_r+cff1_r*hwater
              z_v_o(i,j,k)=z_r0+Zt_avg1v_o(i,j)*(1.0_r8+z_r0*hinv)
            END DO
     END DO
END DO
print *,"Executed set_depts"
END SUBROUTINE set_depths

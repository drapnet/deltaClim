MODULE funcbaro

USE mod_kinds
USE variablesbaro
CONTAINS
SUBROUTINE check(status)
 use netcdf
 integer,intent(in)::status
if (status /= nf90_noerr)then
	print *,trim(nf90_strerror(status))
 	stop "NETcdf problem"
endif
END SUBROUTINE check
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
SUBROUTINE barotropic(uu,vv,hh_5,uubar,vvbar)
USE variablesbaro
USE mod_kinds
real(r8), intent(in) :: uu(1:L-1,1:M,1:N,tim)
real(r8), intent(in) :: vv(1:L,1:M-1,1:N,tim)
real(r8), intent(in) :: hh_5(1:L,1:M,1:N+1)
real(r8) :: uubar(1:L-1,1:M,tim)
real(r8) :: vvbar(1:L,1:M-1,tim)
real(r8) :: Du(1:L-1,1:M)
real(r8) :: Hz_w(1:L,1:M,1:N)
real(r8) :: Dv(1:L,1:M-1)
integer::i,j,k,t
real(r8) :: Duk,Dvk
do i=1,L
	do j=1,M
		do k=N+1,2,-1
			Hz_w(i,j,k-1)=hh_5(i,j,k)-hh_5(i,j,k-1)			
		END DO
	END DO
END DO
DO t=1, tim
      DO i=1,L-1
            DO j=1,M
                DO  k=1,N
                    Duk=0.5*(Hz_w(i,j,k)+Hz_w(i+1,j,k))
                    Du(i,j)=Du(i,j)+Duk
                    uubar(i,j,t)=uubar(i,j,t)+Duk*uu(i,j,k,t)
        	        
		END DO
	  END DO	
	END DO
      uubar(:,:,t)=uubar(:,:,t)/Du   
      DO i=1,L
            DO j=1,M-1
                DO  k=1,N
                    Dvk=0.5*(Hz_w(i,j,k)+Hz_w(i+1,j,k))
                    Dv(i,j)=Dv(i,j)+Dvk
                    vvbar(i,j,t)=vvbar(i,j,t)+Dvk*vv(i,j,k,t)
        	         
		END DO
	  END DO	
       END DO
       vvbar(:,:,t)=vvbar(:,:,t)/Dv  
       Du(:,:)=0.0
       Dv(:,:)=0.0	
END DO

END SUBROUTINE barotropic



END MODULE

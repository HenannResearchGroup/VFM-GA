! Linear Algebra subroutines that are used in the whole code pack.
!! Some of them are from numerical recipes:https://numerical.recipes/
module linalg
    implicit none
        integer d 
        parameter(d=8)
    contains
    subroutine inv3d(A,A_inv)
        implicit none
        real*8 A(3,3),A_inv(3,3),det_A,det_A_inv

        det_A = A(1,1)*(A(2,2)*A(3,3) - A(3,2)*A(2,3)) - &
     &        A(2,1)*(A(1,2)*A(3,3) - A(3,2)*A(1,3)) + &
     &        A(3,1)*(A(1,2)*A(2,3) - A(2,2)*A(1,3))
          
        det_A_inv = 1.d0/det_A
        
        A_inv(1,1) = det_A_inv*(A(2,2)*A(3,3)-A(3,2)*A(2,3))
        A_inv(1,2) = det_A_inv*(A(3,2)*A(1,3)-A(1,2)*A(3,3))
        A_inv(1,3) = det_A_inv*(A(1,2)*A(2,3)-A(2,2)*A(1,3))
        A_inv(2,1) = det_A_inv*(A(3,1)*A(2,3)-A(2,1)*A(3,3))
        A_inv(2,2) = det_A_inv*(A(1,1)*A(3,3)-A(3,1)*A(1,3))
        A_inv(2,3) = det_A_inv*(A(2,1)*A(1,3)-A(1,1)*A(2,3))
        A_inv(3,1) = det_A_inv*(A(2,1)*A(3,2)-A(3,1)*A(2,2))
        A_inv(3,2) = det_A_inv*(A(3,1)*A(1,2)-A(1,1)*A(3,2))
        A_inv(3,3) = det_A_inv*(A(1,1)*A(2,2)-A(2,1)*A(1,2))
        return
    end subroutine inv3d
    
    subroutine det3d(A,detA)
        IMPLICIT NONE
	    REAL(d)  A(3,3), detA
	    detA =	  A(1,1)*A(2,2)*A(3,3) &
&   	        + A(1,2)*A(2,3)*A(3,1) &
&               + A(1,3)*A(2,1)*A(3,2) &
&            	- A(3,1)*A(2,2)*A(1,3) &
&         		- A(3,2)*A(2,3)*A(1,1) &
&         		- A(3,3)*A(2,1)*A(1,2) 

	    RETURN
    END
    
    !product of 2D matrix, subroutine version 
    subroutine mprod(A,B,C)
        implicit none
        real(d) A(2,2),B(2,2),C(2,2)
        integer i,j,k
        do i=1,2
            do j=1,2
                C(i,j) = 0.d0
                do k=1,2
                    C(i,j) = C(i,j)+A(i,k)*B(k,j)
                end do
            end do
        end do
        return
    end subroutine mprod
    !transpose a 2D matrix
    subroutine trans(A,Atrans)
        implicit none
        real(d) A(2,2),Atrans(2,2)
        integer i,j
        do i=1,2
            do j=1,2
                Atrans(j,i) = A(i,j)
            end do
        end do
        return
    end subroutine trans
    
    
    !inverse of 2D matrix
    subroutine inv(A,Ainv)
        implicit none
        real(d) A(2,2),Ainv(2,2),detA,detAinv
        call det(A,detA)
!        if (detA<0) then
!            write(*,*) 'WARNING: subroutine det smaller than 0'
!        end if
        detAinv = 1.d0/detA
        Ainv(1,1) = detAinv*A(2,2)
        Ainv(2,2) = detAinv*A(1,1)
        Ainv(1,2) = detAinv*-A(1,2)
        Ainv(2,1) = detAinv*-A(2,1)
        return 
    end subroutine inv
    
    !trace of 2D matrix
    subroutine trace(A,Atrace)
        implicit none
        real(d) A(2,2),Atrace
        Atrace = A(1,1)+A(2,2)
        return
    end subroutine trace
    
    subroutine trace3d(A,Atrace)
        implicit none
        real(d) A(3,3),Atrace
        Atrace = A(1,1)+A(2,2)+A(3,3)
        return
    end subroutine trace3d
    
    !det of 2D matrix
    subroutine det(A,detA)
        implicit none
        real(d) A(2,2),detA
        detA = A(1,1)*A(2,2)-A(1,2)*A(2,1)
        return 
    end subroutine det
    
    !create nxn matrix with 1 on diagonal 
    subroutine eye(A,size)
        implicit none
        integer size,i,j
        real(d) A(size,size)
        do i = 1,size
            do j = 1,size
                if (i == j) then
                    A(i,j) = 1.d0
                else
                    A(i,j) = 0.d0
                end if
            end do
        end do
        return 
    end subroutine eye
    
    !initial a matrix/size of zeros
    subroutine zeros(A,row,column)
        implicit none
        integer row,column
        real(d) A(row,column)
        A = 0.d0
        return
    end subroutine zeros
    
    subroutine ones(A,row,column)
        implicit none
        integer row,column
        real(d) A(row,column)
        A = 1.d0
        return
    end subroutine ones
    
    subroutine eigsrt(D,V,N,NP)
	    IMPLICIT NONE
	    INTEGER N,NP,I,J,K
	    REAL*8 D(NP),V(NP,NP),P
	    DO I = 1,N-1
            K = I
            P = D(I)
	        DO J = I+1,N
	            IF (D(J) .GE. P) THEN
	                K = J
	                P = D(J)
	            END IF
	        END DO
	        IF (K .NE. I) THEN
	            D(K) = D(I)
	            D(I) = P
	            DO J = 1,N
	                P = V(J,I)
	                V(J,I) = V(J,K)
	                V(J,K) = P
                END DO
  	        END IF
	    END DO
	    RETURN
    end subroutine eigsrt
    
    subroutine skinem(F,R,U,E)
        real(d) R(2,2),U(2,2),UINV(2,2),E(2,2),F(2,2),DETF,FTRANS(2,2),C(2,2)
        real(d) OMEGA(2),EIGVEC(2,2),UEIGVAL(2),TEMPM(2,2),EIGVECT(2,2)
    	CALL eye(R,2)
	    CALL eye(U,2)
	    CALL eye(UINV,2)
	    CALL eye(E,2)
        CALL det(F,DETF)
        IF (DETF .LE. 0.D0) THEN
            !WRITE(*,'(/5X,A/)') '--ERROR IN KINEMATICS-- THE', &
     !&         ' DETERMINANT OF [F] IS NOT GREATER THAN 0'
            RETURN
        END IF
        CALL  trans(F,FTRANS)
        CALL  mprod(FTRANS,F,C)
	    CALL spectral(C,OMEGA,EIGVEC)
        
	    UEIGVAL(1) = dsqrt(OMEGA(1))
	    UEIGVAL(2) = dsqrt(OMEGA(2))
	    U(1,1) = UEIGVAL(1)
	    U(2,2) = UEIGVAL(2)
	    E(1,1) = dlog( UEIGVAL(1) )
	    E(2,2) = dlog( UEIGVAL(2) )
        U = matmul(matmul(eigvec,U),transpose(eigvec))
        E = matmul(matmul(eigvec,E),transpose(eigvec))
	    CALL inv(U,UINV)
	    CALL mprod(F,UINV,R)
        E = matmul(matmul(R,E),transpose(R))
	    RETURN
    end subroutine skinem
    
    subroutine skinem3d(F,R,U,E)
        real(d) R(3,3),U(3,3),UINV(3,3),E(3,3),F(3,3),DETF,FTRANS(3,3),C(3,3)
        real(d) OMEGA(3),EIGVEC(3,3),UEIGVAL(3),TEMPM(3,3),EIGVECT(3,3)
    	CALL eye(R,3)
	    CALL eye(U,3)
	    CALL eye(UINV,3)
        CALL eye(E,3)
	    E = 0.d0
        CALL det3d(F,DETF)
        IF (DETF .LE. 0.D0) THEN
  !          WRITE(*,'(/5X,A/)') '--ERROR IN KINEMATICS-- THE', &
  !   &         ' DETERMINANT OF [F] IS NOT GREATER THAN 0'
            RETURN
        END IF
        C = matmul(transpose(F),F)
	    CALL spectral3d(C,OMEGA,EIGVEC)
        
	    UEIGVAL(1) = dsqrt(OMEGA(1))
	    UEIGVAL(2) = dsqrt(OMEGA(2))
        UEIGVAL(3) = dsqrt(OMEGA(3))
	    U(1,1) = UEIGVAL(1)
	    U(2,2) = UEIGVAL(2)
        U(3,3) = UEIGVAL(3)
	    E(1,1) = dlog( UEIGVAL(1) )
	    E(2,2) = dlog( UEIGVAL(2) )
        E(3,3) = dlog( UEIGVAL(3) )
        U = matmul(matmul(eigvec,U),transpose(eigvec))
        E = matmul(matmul(eigvec,E),transpose(eigvec))
	    call inv3D(U,Uinv)
	    R = matmul(F,Uinv)
        E = matmul(matmul(R,E),transpose(R))
	    RETURN
    end subroutine skinem3d
    
    subroutine skinem2(FF,RR,UU,EE)
        real(d) R(3,3),U(3,3),UINV(3,3),E(3,3),F(3,3),DETF,FTRANS(3,3),C(3,3)
        real(d) OMEGA(3),EIGVEC(3,3),UEIGVAL(3),TEMPM(3,3),EIGVECT(3,3)
        real(d) FF(2,2),RR(2,2),UU(2,2),EE(2,2)
        call eye(F,3)
        F(1:2,1:2) = FF(1:2,1:2)
    	CALL eye(R,3)
	    CALL eye(U,3)
	    CALL eye(UINV,3)
	    E = 0.d0
        CALL det3d(F,DETF)
        IF (DETF .LE. 0.D0) THEN
     !       WRITE(*,'(/5X,A/)') '--ERROR IN KINEMATICS-- THE', &
    ! &         ' DETERMINANT OF [F] IS NOT GREATER THAN 0'
            RETURN
        END IF
        C = matmul(transpose(F),F)
	    CALL spectral3d(C,OMEGA,EIGVEC)
        
	    UEIGVAL(1) = dsqrt(OMEGA(1))
	    UEIGVAL(2) = dsqrt(OMEGA(2))
        UEIGVAL(3) = dsqrt(OMEGA(3))
	    U(1,1) = UEIGVAL(1)
	    U(2,2) = UEIGVAL(2)
        U(3,3) = UEIGVAL(3)
	    E(1,1) = dlog( UEIGVAL(1) )
	    E(2,2) = dlog( UEIGVAL(2) )
        E(3,3) = dlog( UEIGVAL(3) )
        U = matmul(matmul(eigvec,U),transpose(eigvec))
        E = matmul(matmul(eigvec,E),transpose(eigvec))
	    call inv3D(U,Uinv)
	    R = matmul(F,Uinv)
        E = matmul(matmul(R,E),transpose(R))
        RR(1:2,1:2) = R(1:2,1:2)
        EE(1:2,1:2) = E(1:2,1:2)
        UU(1:2,1:2) = U(1:2,1:2)
	    RETURN
    end subroutine skinem2
    
    subroutine spectral(A,D,V)
        IMPLICIT NONE
        REAL*8 D(2),V(2,2),A(2,2),E(2,2)
        INTEGER NP,NROT,I,J
        PARAMETER(NP=2)

        DO I = 1,2
            DO J= 1,2
                E(I,J) = A(I,J)
            END DO
        END DO

        CALL jacobi(E,2,NP,D,V,NROT)
        CALL eigsrt(D,V,2,NP)

        RETURN  
    end subroutine spectral
    
    subroutine spectral3d(A,D,V)

      implicit none
      integer np,nrot,i,j,istat
      parameter(np=3)
      real*8 D(3),V(3,3),A(3,3),E(3,3)
      E = A
      call jacobi(E,3,np,D,V,nrot)
      call eigsrt(D,V,3,np)
      return
      end subroutine spectral3d
    
    subroutine jacobi(A,N,NP,D,V,NROT)
    	IMPLICIT NONE
	    INTEGER IP,IQ,N,NMAX,NP,NROT,I,J
	    PARAMETER (NMAX =100)
	    REAL*8 A(NP,NP),D(NP),V(NP,NP),B(NMAX),Z(NMAX)
	    REAL*8 SM,TRESH,G,T,H,THETA,S,C,TAU

	    DO IP = 1,N	
	        DO IQ = 1,N
	            V(IP,IQ) = 0.D0
            END DO
            V(IP,IP) = 1.D0
	    END DO
	    DO IP = 1,N
	        B(IP) = A(IP,IP)
	        D(IP) = B(IP)
	        Z(IP) = 0.D0
        END DO
	    NROT = 0
	    DO I = 1,20 !20 or 10
            SM = 0.D0
            DO IP = 1, N-1
                DO IQ = IP + 1, N
	                SM = SM + DABS ( A(IP,IQ ))
                END DO
            END DO
            IF ( SM .EQ. 0.D0) RETURN
            IF ( I .LT. 4) THEN
                TRESH = 0.2D0*SM/N**2
            ELSE
                TRESH = 0.D0
            END IF

            DO IP = 1, N-1
                DO IQ = IP+1,N
                    G = 100.D0*DABS(A(IP,IQ))
	                IF ((I .GT. 4) .AND. (DABS(D(IP))+G .EQ. DABS(D(IP))) &
&                       .AND. ( DABS(D(IQ))+G .EQ. DABS(D(IQ)))) THEN
                        A(IP,IQ) = 0.D0
                    ELSE IF ( DABS(A(IP,IQ)) .GT. TRESH) THEN
                        H = D(IQ) - D(IP)
                        IF (DABS(H)+G .EQ. DABS(H)) THEN
                            T =A(IP,IQ)/H
	                    ELSE
	                        THETA = 0.5D0*H/A(IP,IQ)
                            T =1.D0/(DABS(THETA)+DSQRT(1.D0+THETA**2))
	                        IF (THETA .LT. 0.D0) T = -T
	                    END IF
	                    C = 1.D0/DSQRT(1.D0 + T**2)
	                    S = T*C
	                    TAU = S/(1.D0 + C)
	                    H = T*A(IP,IQ)
	                    Z(IP) = Z(IP) - H
	                    Z(IQ) = Z(IQ) + H
	                    D(IP) = D(IP) - H
	                    D(IQ) = D(IQ) + H
	                    A(IP,IQ) = 0.D0
	                    DO J = 1, IP-1
	                      G = A(J,IP)
	                      H = A(J,IQ)
	                      A(J,IP) = G - S*(H + G*TAU)
	                      A(J,IQ) = H + S*(G - H*TAU)
	                    END DO

	                    DO J = IP+1, IQ-1
	                      G = A(IP,J)
	                      H = A(J,IQ)
	                      A(IP,J) = G - S*(H + G*TAU)
	                      A(J,IQ) = H + S*(G - H*TAU)
	                    END DO

	                    DO J = IQ+1, N
                              G = A(IP,J)
	                      H = A(IQ,J)
	                      A(IP,J) = G - S*(H + G*TAU)
	                      A(IQ,J) = H + S*(G - H*TAU)
	                    END DO
	                    DO J = 1,N
	                      G = V(J,IP)
	                      H = V(J,IQ)
	                      V(J,IP) = G - S*(H + G*TAU)
	                      V(J,IQ) = H + S*(G - H*TAU)
	                    END DO
	                    NROT = NROT + 1
                    END IF
	            END DO
	        END DO
            DO IP = 1, N
                B(IP) = B(IP) + Z(IP)
                D(IP) = B(IP)
                Z(IP) = 0.D0
            END DO
	    END DO
	    !WRITE (*,'(/1X,A/)') '50 ITERS IN JACOBI SHOULD NEVER HAPPEN'
	    RETURN
    end subroutine jacobi
    subroutine dlnxdx(X,DYDX)
      !
      ! This subroutine calculates the derivative of the logarithm
      ! of a symmetric tensor with respect to that tensor
      !
      implicit none
      !
      integer i,j,k,l
      !
      real(8) X(3,3),DYDX(3,3,3,3),Iden(3,3),Iden4(3,3,3,3),eigval(3), &
     &  eigvec(3,3),ehat1(3),ehat2(3),ehat3(3),E1(3,3),E2(3,3),E3(3,3), &
     &  y(3),DX2DX(3,3,3,3),s1,s2,s3,s4,s5,s6
      !
      real(8) zero,one,two,half,three,third,small
      parameter(zero=0.d0,one=1.d0,two=2.d0,half=0.5d0, &
     &     three=3.d0,third=1.d0/3.d0,small=1.d-12)
      
      
      ! Initialize
      !
      ! Second order identity tensor
      !
      call eye(Iden,3)
      !
      ! Fourth order symmetric identity tensor
      !
      Iden4 = zero
      do i=1,3
        do j=1,3
          do k=1,3
            do l=1,3
              Iden4(i,j,k,l) = half*(Iden(i,k)*Iden(j,l) + &
     &                                Iden(i,l)*Iden(j,k))
            end do
          end do
        end do
      end do
      
      
      ! Calculate the eigenvalues and eigenvectors of X
      !
      call spectral3d(X,eigval,eigvec)
      !
      ! Extract the eigenvectors
      !
      do i=1,3
        ehat1(i) =  eigvec(i,1)
        ehat2(i) =  eigvec(i,2)
        ehat3(i) =  eigvec(i,3)
      end do
      !
      ! Assemble the eigenprojections
      !
      do i=1,3
        do j=1,3
	  E1(i,j) = ehat1(i)*ehat1(j)
	  E2(i,j) = ehat2(i)*ehat2(j)
	  E3(i,j) = ehat3(i)*ehat3(j)
	end do
      end do
      
      
      ! Calculate the eigenvalues of Y = ln(X)
      !
      y = dlog(eigval)
      
      
      ! Calculate the derivative of X^2 with respect to X
      !
      DX2DX = zero
      do i=1,3
        do j=1,3
          do k=1,3
            do l=1,3
              DX2DX(i,j,k,l) = half*(Iden(i,k)*X(j,l) +  &
     &                                Iden(i,l)*X(j,k) +  &
     &                                X(i,k)*Iden(j,l) +  &
     &                                X(i,l)*Iden(j,k))
            end do
          end do
        end do
      end do
         
            
      ! Calculate DYDX
      !
      DYDX = zero
      if (dabs(eigval(1)-eigval(3)).le.small) then
        !
        ! Three repeated eigenvalues
        !
        do i=1,3
          do j=1,3
            do k=1,3
              do l=1,3
                DYDX(i,j,k,l) = (one/eigval(1))*Iden4(i,j,k,l)
              end do
            end do
          end do
        end do
        !
      elseif (dabs(eigval(2)-eigval(3)).le.small) then
        !
        ! The eigenvalues 2 and 3 are repeated. Eigenvalue 1 is distinct.
        !
        s1 = (y(1) - y(2))/((eigval(1)-eigval(2))**two) - &
     &                  (one/eigval(2))/(eigval(1)-eigval(2))
        s2 = two*eigval(2)*(y(1)-y(2))/((eigval(1)-eigval(2))**two) - &
     &     (one/eigval(2))*(eigval(1)+eigval(2))/(eigval(1)-eigval(2))
        s3 = two*(y(1)-y(2))/((eigval(1)-eigval(2))**three) -  &
     &   ((one/eigval(1))+(one/eigval(2)))/((eigval(1)-eigval(2))**two)
        s4 = eigval(2)*s3
        s5 = s4
        s6 = (eigval(2)**two)*s3
        !
        do i=1,3
          do j=1,3
            do k=1,3
              do l=1,3
                DYDX(i,j,k,l) = s1*DX2DX(i,j,k,l) - s2*Iden4(i,j,k,l) - &
     &                     s3*X(i,j)*X(k,l) + s4*X(i,j)*Iden(k,l) +  &
     &                     s5*Iden(i,j)*X(k,l) - s6*Iden(i,j)*Iden(k,l)
              end do
            end do
          end do
        end do
        !
      elseif (dabs(eigval(1)-eigval(2)).le.small) then
        !
        ! The eigenvalues 1 and 2 are repeated. Eigenvalue 3 is distinct.
        !
        s1 = (y(3) - y(2))/((eigval(3)-eigval(2))**two) - &
     &                  (one/eigval(2))/(eigval(3)-eigval(2))
        s2 = two*eigval(2)*(y(3)-y(2))/((eigval(3)-eigval(2))**two) - &
     &     (one/eigval(2))*(eigval(3)+eigval(2))/(eigval(3)-eigval(2))
        s3 = two*(y(3)-y(2))/((eigval(3)-eigval(2))**three) - &
     &   ((one/eigval(3))+(one/eigval(2)))/((eigval(3)-eigval(2))**two)
        s4 = eigval(2)*s3
        s5 = s4
        s6 = (eigval(2)**two)*s3
        !
        do i=1,3
          do j=1,3
            do k=1,3
              do l=1,3
                DYDX(i,j,k,l) = s1*DX2DX(i,j,k,l) - s2*Iden4(i,j,k,l) - &
     &                     s3*X(i,j)*X(k,l) + s4*X(i,j)*Iden(k,l) + &
     &                     s5*Iden(i,j)*X(k,l) - s6*Iden(i,j)*Iden(k,l)
              end do
            end do
          end do
        end do
        !
      else
        !
        ! Eigenvalues are distinct.
        !
        do i=1,3
          do j=1,3
            do k=1,3
              do l=1,3
                DYDX(i,j,k,l) = (y(1)/((eigval(1)-eigval(2))*      &
     &                                 (eigval(1)-eigval(3))))* &
     &         (DX2DX(i,j,k,l) - (eigval(2)+eigval(3))*Iden4(i,j,k,l) -  &
     &  ((eigval(1)-eigval(2))+(eigval(1)-eigval(3)))*E1(i,j)*E1(k,l) - &
     &       (eigval(2)-eigval(3))*(E2(i,j)*E2(k,l)-E3(i,j)*E3(k,l))) + &
     &                        (one/eigval(1))*E1(i,j)*E1(k,l) + &
     &                          (y(2)/((eigval(2)-eigval(1))* &
     &                                 (eigval(2)-eigval(3))))* &
     &         (DX2DX(i,j,k,l) - (eigval(1)+eigval(3))*Iden4(i,j,k,l) -  &
     &  ((eigval(2)-eigval(1))+(eigval(2)-eigval(3)))*E2(i,j)*E2(k,l) -  &
     &       (eigval(1)-eigval(3))*(E1(i,j)*E1(k,l)-E3(i,j)*E3(k,l))) +  &
     &                        (one/eigval(2))*E2(i,j)*E2(k,l) + &
     &                          (y(3)/((eigval(3)-eigval(1))* &
     &                                 (eigval(3)-eigval(2))))* &
     &         (DX2DX(i,j,k,l) - (eigval(1)+eigval(2))*Iden4(i,j,k,l) -  &
     &  ((eigval(3)-eigval(1))+(eigval(3)-eigval(2)))*E3(i,j)*E3(k,l) - &
     &       (eigval(1)-eigval(2))*(E1(i,j)*E1(k,l)-E2(i,j)*E2(k,l))) + &
     &                        (one/eigval(3))*E3(i,j)*E3(k,l)
              end do
            end do
          end do
        end do
        !
      end if
      
      return
      end subroutine dlnxdx

!****************************************************************************

      subroutine mprod4(A,B,C)
      !
      ! This subroutine calculates the product of two fourth order tensors,
      ! A and B, and places the product in C:
      !             C_{ijkl} = A_{ijmn}B_{mnkl}.
      !
      implicit none
      !
      integer i,j,k,l,m,n
      !
      real(8) A(3,3,3,3),B(3,3,3,3),C(3,3,3,3)
      
      
      C = 0.d0
      do i=1,3
        do j=1,3
          do k=1,3
            do l=1,3
              do m=1,3
                do n=1,3
                  C(i,j,k,l) = C(i,j,k,l) + A(i,j,m,n)*B(m,n,k,l)
                end do
              end do
            end do
          end do
        end do
      end do
      
      
      return
      end subroutine mprod4
end module linalg
    
    
    
    
    
!!!This module contains constitutive equation and solver
!!!! To apply this code package to any model, you have to write a subroutine that takes deformation gradient and random material parameter as input, and output Kirchoff stress
module fem_constitutive
    use linalg
    use constants
    implicit none
    integer dd
    contains
    !! Neohookean constitutive law plane strain here as an example, never used, just for fun, 
    !! Nobody should use such a algorithm to fit a neohookean which can be done by hand
    subroutine const_stress_neo(F,theta,stress,nf)
        implicit none
        real(di) F(2,2),F_trans(2,2),B(2,2),I1,J,I1_,dQdI1_,dQdJ, &
&               dI1_dF(2,2),dJdF(2,2),Finv(2,2),Finv_t(2,2),theta(2), &
&               stress(2,2)
        integer nf 
        nf = 2
        call trans(F,F_trans)
        B = matmul(F,F_trans)
        call trace(B,I1)
        I1 = I1+1
        call det(F,J)
        call inv(F,Finv)
        call trans(Finv,Finv_t)
        I1_ = I1/(J ** (2.d0/3.d0))
        dQdI1_ = 1.d0
        dQdJ = 2.d0*(J-2.d0)
        dI1_dF = 2.d0/(J ** (2.d0/3.d0))*F - 2.d0/3.d0*I1_*Finv_t
        dJdF = J*Finv_t
!       http://solidmechanics.org/text/Chapter3_5/Chapter3_5.htm
        stress = dQdI1_ * dI1_dF * F_trans*theta(1) + dQdJ*dJdF*F_trans*theta(2)
        return 
    end subroutine const_stress_neo
    !! Neohookean gradient wrt each material parameter, for fun only 
    subroutine const_grad_neo(F,theta,grad,nf)
        implicit none
        real(di) F(2,2),F_trans(2,2),B(2,2),I1,J,I1_,dQdI1_,dQdJ, &
&               dI1_dF(2,2),dJdF(2,2),Finv(2,2),Finv_t(2,2),theta(2), &
&               grad(2,4)
        integer nf 
        nf = 2
        call trans(F,F_trans)
        B = matmul(F,F_trans)
        call trace(B,I1)
        I1 = I1+1
        call det(F,J)
        call inv(F,Finv)
        call trans(Finv,Finv_t)
        I1_ = I1/(J ** (2.d0/3.d0))
        dQdI1_ = 1.d0
        dQdJ = 2.d0*(J-2.d0)
        dI1_dF = 2.d0/(J ** (2.d0/3.d0))*F - 2.d0/3.d0*I1_*Finv_t
        dJdF = J*Finv_t
!       http://solidmechanics.org/text/Chapter3_5/Chapter3_5.htm
        grad(1:2,1:2) = dQdI1_ * dI1_dF * F_trans
        grad(1:2,3:4) = dQdJ * dJdF * F_trans
        return 
    end subroutine const_grad_neo
    !!  Plane strain, constitutive law for foam based on spatial log strain measure https://www.sciencedirect.com/science/article/pii/S0022509619303825
    subroutine const_stress_foam(F,theta,stress,nf)
        implicit none
        real(di) F(2,2),E(2,2),K1,devE(2,2),ey(2,2),sumdevE, &
&               K2,N(2,2),Ntemp(3,3),K3,Ntemp_det,J,detF, &
&               G0,B,Jmin,C1,K10,delta_K,X1,X2,C0,theta(14), &
&               p,q,C2,C3,r,X,dXdJ1,dfdK1,dLdK2,dLdK3,comp1(2,2), &
&               comp2(2,2),U(2,2),RR(2,2),Y(2,2),dXdK1, &
&               stress(2,2),test,E_temp(2,2),Rot(2,2)
        integer nf,i,jj
        nf = 14
        call skinem(F,RR,U,E)
        !call trans(RR,Rot)
        !E= matmul(matmul(RR,E_temp),Rot)
        call trace(E,K1)
        call eye(ey,2)
        devE = E-ey*1.d0/3.d0*K1
        sumdevE = 0.d0
        do i = 1,2
            do jj = 1,2
                sumdevE = sumdevE + devE(i,jj) ** 2.d0
            end do
        end do
        sumdevE = sumdevE + (-K1/3.d0)**2.d0
        K2 = dsqrt(sumdevE)
        N = devE/K2
        Ntemp(1:2,1:2) = N(1:2,1:2)
        Ntemp(3,3) = -(K1/3.d0)/K2
        call det3d(Ntemp,Ntemp_det)
        K3 = 3.d0*(6.d0 ** 0.5d0)*Ntemp_det
        Y = 3.d0*(6.d0 ** 0.5d0)*matmul(N,N) &
&        -(6.d0 ** 0.5d0)*ey-3.d0*K3*N
        J = 0.d0
        call det(F,J)
        call fill(G0,B,Jmin,C1,K10,delta_K,X1,X2,C0,p,q,C2,C3,r,theta)
        X = 0.5d0*(X1+X2)*K1+delta_K/2.d0*(X1-X2)*dlog(dcosh((K1-K10)/delta_K)/dcosh(K10/delta_K))+1.d0
        dXdK1 = 0.5d0*(X1+X2)+0.5d0*(X1-X2)*dtanh((K1-K10)/delta_K)
        if ((J-Jmin)<0.d0) then
            dfdK1 = (dexp(C2*K1)-1.d0)/C2+C3*J*(J**(-r)-(-1)*((1.d0-Jmin) &
&       /(Jmin-J))**r)
        else
            dfdK1 = (dexp(C2*K1)-1.d0)/C2+C3*J*(J**(-r)-((1.d0-Jmin) &
&       /(J-Jmin))**r)
        end if 
        dLdK2 = p*C0*K2 ** (p-1.d0)+q*C1*(1.d0+K3)*K2 ** (q-1.d0)
        dLdK3 = C1 * K2 ** q
        comp1 = J**(-1.d0)*(ey*G0*dXdK1*K2**2.d0+ &
&       G0*(2.d0*X*K2+dLdK2)*N+G0*dLdK3*Y/K2)
        comp2 = J ** (-1.d0)*(B*dfdK1)*ey
        stress = (comp1+comp2)*J
        return        
    end subroutine const_stress_foam
    !! For fun, never used, gradient wrt each material parameter for foam
    subroutine const_grad_foam(F,theta,grad,nf)
        implicit none
        real(di) F(2,2),E(2,2),K1,devE(2,2),ey(2,2),sumdevE, &
&               K2,N(2,2),Ntemp(3,3),K3,Ntemp_det,J,detF, &
&               G0,B,Jmin,C1,K10,delta_K,X1,X2,C0, &
&               p,q,C2,C3,r,X,D1(2,2),D2(2,2),D3(2,2),D4(2,2), &
&               U(2,2),RR(2,2),Y(2,2),theta(14), &
&               grad(2,28),D5(2,2),D6(2,2),D7(2,2),D8(2,2), &
&               D9(2,2),D10(2,2),D11(2,2),D12(2,2),D13(2,2),D14(2,2) 
        integer nf,i,jj
        nf = 14
        call skinem(F,RR,U,E)
        call trace(E,K1)
        call eye(ey,2)
        devE = E-ey*1.d0/3.d0*K1
        sumdevE = 0.d0
        do i = 1,2
            do jj = 1,2
                sumdevE = sumdevE + devE(i,jj) ** 2.d0
            end do
        end do
        sumdevE = sumdevE + (-K1/3.d0)**2.d0
        K2 = dsqrt(sumdevE)
        N = devE/K2
        Ntemp(1:2,1:2) = N(1:2,1:2)
        Ntemp(3,3) = -(K1/3.d0)/K2
        call det(Ntemp,Ntemp_det)
        K3 = 3.d0*(6.d0 ** 0.5d0)*Ntemp_det
        Y = 3.d0*(6.d0 ** 0.5d0)*matmul(N,N) &
&        -(6.d0 ** 0.5d0)*ey-3.d0*K3*N
        call det(F,J)
        call fill(G0,B,Jmin,C1,K10,delta_K,X1,X2,C0,p,q,C2,C3,r,theta)
        D1 = (N*(K2*(K1*(X1 + X2) + delta_K*dlog(dcosh((K1 - K10) &
&        /delta_K)/dcosh(K10/delta_K))*(X1 - X2) + 2.d0) + C0*K2** &
&        (p - 1.d0)*p + C1*K2**(q - 1.d0)*q*(K3 + 1.d0)) + ey*K2** &
&        2.d0*(X1/2 + X2/2 + dtanh((K1 - K10)/delta_K)* &
&         (X1/2.d0 - X2/2.d0)) + C1*K2**(q - 1.d0)*Y)/J
        D2 = (ey*((dexp(C2*K1) - 1.d0)/C2 + C3*J* &
&       (1/J**r - (-(Jmin - 1.d0)/(J - Jmin))**r)))/J
        D3 = (B*C3*ey*r*(J - 1.d0)*(-(Jmin - 1.d0)/ &
&       (J - Jmin))**(r - 1.d0))/(J - Jmin)**2.d0
        D4 = (G0*K2**(q - 1.d0)*(Y + N*q + K3*N*q))/J
        D5 = ((G0*ey*(dtanh((K1 - K10)/delta_K)**2.d0 - 1.d0) &
&       *(X1/2.d0 - X2/2.d0)*K2 ** 2.d0)/delta_K - &
&       (G0*N*dsinh(K1/delta_K)*(X1 - X2)*K2)/ &
&       (dcosh(K10/delta_K)*dcosh((K1 - K10)/delta_K)))/J
        D6 = (G0*K2*N*(dlog(dcosh((K1 - K10)/delta_K)/ &
&       dcosh(K10/delta_K))*(X1 - X2) - (delta_K*dcosh(K10/delta_K)* &
&       ((dsinh((K1 - K10)/delta_K)*(K1 - K10))/(delta_K ** 2.d0* &
&       dcosh(K10/delta_K)) - (K10*dcosh((K1 - K10)/delta_K)* &
&       dsinh(K10/delta_K))/(delta_K ** 2.d0*dcosh(K10/delta_K) &
&       ** 2.d0))*(X1 - X2))/dcosh((K1 - K10)/delta_K)) + (G0*ey*K2 &
&       ** 2.d0*(K1 - K10)*(dtanh((K1 - K10)/delta_K) &
&       ** 2.d0 - 1.d0)*(X1/2.d0 - X2/2.d0))/delta_K ** 2.d0)/J
        D7 = (G0*ey*(dtanh((K1 - K10)/delta_K)/2.d0 + 1.d0/2.d0)*K2 &
&       ** 2.d0 + G0*N*(K1 + delta_K*dlog(dcosh((K1 - K10)/delta_K) &
&       /dcosh(K10/delta_K)))*K2)/J
        D8 = -(G0*ey*(dtanh((K1 - K10)/delta_K)/2.d0 - 1.d0/2.d0)*K2 &
&       ** 2.d0 - G0*N*(K1 - delta_K*dlog(dcosh((K1 - K10)/delta_K) &
&       /dcosh(K10/delta_K)))*K2)/J
        D9 = (G0*K2 ** (p - 1.d0)*N*p)/J
        D10 = (C0*G0*K2 ** (p - 1.d0)*N*(p*dlog(K2) + 1.d0))/J
        D11 = (C1*G0*K2 ** (q - 1.d0)*(N + Y*dlog(K2) + K3*N + &
&       N*q*dlog(K2) + K3*N*q*dlog(K2)))/J
        D12 = -(B*ey*((dexp(C2*K1) - 1.d0)/C2 ** 2.d0 - &
&       (K1*dexp(C2*K1))/C2))/J
        D13 = B*ey*(1.d0/J ** r - (-(Jmin - 1.d0)/(J - Jmin)) ** r)
        D14 = -B*C3*ey*(dlog(-(Jmin - 1.d0)/(J - Jmin))* &
&       (-(Jmin - 1.d0)/(J - Jmin)) ** r + dlog(J)/J ** r)
        grad(1:2,1:2) = D1(1:2,1:2)
        grad(1:2,3:4) = D2(1:2,1:2)
        grad(1:2,5:6) = D3(1:2,1:2)
        grad(1:2,7:8) = D4(1:2,1:2)
        grad(1:2,9:10) = D5(1:2,1:2)
        grad(1:2,11:12) = D6(1:2,1:2)
        grad(1:2,13:14) = D7(1:2,1:2)
        grad(1:2,15:16) = D8(1:2,1:2)
        grad(1:2,17:18) = D9(1:2,1:2)
        grad(1:2,19:20) = D10(1:2,1:2)
        grad(1:2,21:22) = D11(1:2,1:2)
        grad(1:2,23:24) = D12(1:2,1:2)
        grad(1:2,25:26) = D13(1:2,1:2)
        grad(1:2,27:28) = D14(1:2,1:2)
        return
    end subroutine const_grad_foam
        
    subroutine fill(a,b,c,d,e,f,g,h,i,j,k,l,m,n,theta)
        implicit none
        real(di) a,b,c,d,e,f,g,h,i,j,k,l,m,n,theta(14)
        a = theta(1)
        b = theta(2)
        c = theta(3)
        d = theta(4)
        e = theta(5)
        f = theta(6)
        g = theta(7)
        h = theta(8)
        i = theta(9)
        j = theta(10)
        k = theta(11)
        l = theta(12)
        m = theta(13)
        n = theta(14)
        return 
    end subroutine fill
    !! Plane stress, constitutive law for foam based on spatial log strain measure https://www.sciencedirect.com/science/article/pii/S0022509619303825
    !! plane stress is used way more often due to experimental concerns, see paper more details
    subroutine const_stress_foam_ps(F,theta,stress,nf)
        implicit none
        real(di) F(3,3),E(3,3),K1,devE(3,3),ey(3,3),sumdevE, &
&               K2,N(3,3),K3,Ntemp_det,J,detF, &
&               G0,B,Jmin,C1,K10,delta_K,X1,X2,C0,theta(14), &
&               p,q,C2,C3,r,X,dXdJ1,dfdK1,dLdK2,dLdK3,comp1(3,3), &
&               comp2(3,3),U(3,3),RR(3,3),Y(3,3),dXdK1, &
&               stress(3,3),test,E_temp(3,3),Rot(3,3)
        integer nf,i,jj
        nf = 14
        call skinem3d(F,RR,U,E)
        !call trans(RR,Rot)
        !E= matmul(matmul(RR,E_temp),Rot)
        call trace3d(E,K1)
        call eye(ey,3)
        devE = E-ey*1.d0/3.d0*K1
        sumdevE = 0.d0
        do i = 1,3
            do jj = 1,3
                sumdevE = sumdevE + devE(i,jj) ** 2.d0
            end do
        end do
        K2 = dsqrt(sumdevE)
        N = devE/K2
        call det3d(N,Ntemp_det)
        K3 = 3.d0*(6.d0 ** 0.5d0)*Ntemp_det
        Y = 3.d0*(6.d0 ** 0.5d0)*matmul(N,N) &
&        -(6.d0 ** 0.5d0)*ey-3.d0*K3*N
        J = 0.d0
        call det3d(F,J)
        call fill(G0,B,Jmin,C1,K10,delta_K,X1,X2,C0,p,q,C2,C3,r,theta)
        X = 0.5d0*(X1+X2)*K1+delta_K/2.d0*(X1-X2)*dlog(dcosh((K1-K10)/delta_K)/dcosh(K10/delta_K))+1.d0
        dXdK1 = 0.5d0*(X1+X2)+0.5d0*(X1-X2)*dtanh((K1-K10)/delta_K)
        if ((J-Jmin)<0.d0) then
            dfdK1 = (dexp(C2*K1)-1.d0)/C2+C3*J*(J**(-r)-(-1)*((1.d0-Jmin) &
&       /(Jmin-J))**r)
        else
            dfdK1 = (dexp(C2*K1)-1.d0)/C2+C3*J*(J**(-r)-((1.d0-Jmin) &
&       /(J-Jmin))**r)
        end if 
        dLdK2 = p*C0*K2 ** (p-1.d0)+q*C1*(1.d0+K3)*K2 ** (q-1.d0)
        dLdK3 = C1 * K2 ** q
        comp1 = J**(-1.d0)*(ey*G0*dXdK1*K2**2.d0+ &
&       G0*(2.d0*X*K2+dLdK2)*N+G0*dLdK3*Y/K2)
        comp2 = J ** (-1.d0)*(B*dfdK1)*ey
        stress = (comp1+comp2)
        return        
    end subroutine const_stress_foam_ps
    !! plane stress solver to solve for the out of plane deformation gradient to fullfill plane stress approximation 
    !! (bisection for the first element, then use F33_out as a guess to plug in newton raphson for element at vincinity)
    subroutine solve_plane_stress(F,theta,F_out,Smean,flag,F33_out)
        implicit none
        real(di) F3max,F3min,F3mean,Fmax(3,3),Fmin(3,3),Fmean(3,3),F(2,2),F33_out
        real(di) Smax(3,3),Smin(3,3),Smean(3,3),theta(14),F_out(3,3),S33_ratio,S33_good
        integer ite1,ite1_max,ite2,ite2_max,flag,nf
        ! the max and min F33 for simple compression and simple tension        
        F3max = 1.d0
        F3min = 0.64d0
        call eye(Fmax,3)
        call eye(Fmin,3)
        call eye(Fmean,3)
        Fmax(1:2,1:2) = F(1:2,1:2)
        Fmin(1:2,1:2) = F(1:2,1:2)
        Fmean(1:2,1:2) = F(1:2,1:2)
        ! ite1 max for expanding bracket
        ite1 = 0
        ite1_max = 3
        ! ite2 max for maximum number of iterations
        ite2 = 0
        ite2_max = 10
        ! S33_good is the required ratio between in plane and out of plane stress
        S33_good = 0.001d0
        flag = 0
        Smean = 0.d0
        F_out = 0.d0
        F33_out = 1.d0
        do while (ite1.le.ite1_max)
            if (ite1 == 0) then
                F3mean = 1.d0
            else
                F3mean = (F3max-F3min)/2.d0
            end if 
            
            Fmax(3,3) = F3max
            Fmin(3,3) = F3min
            Fmean(3,3) = F3mean
            call const_stress_foam_ps(Fmax,theta,Smax,nf)
            call const_stress_foam_ps(Fmin,theta,Smin,nf)
            call const_stress_foam_ps(Fmean,theta,Smean,nf)
            if ((Smax(3,3)*Smean(3,3)).le.0) then 
                exit
            else if ((Smin(3,3)*Smean(3,3)).le.0) then 
                exit
            else
                if (ite1 == ite1_max) then
                    flag = 1
                    return
                end if
                F3max = F3max*2.d0
                F3min = F3min/1.5d0
                ite1 = ite1+1
            end if
        end do
        S33_ratio = 1.d0
        do while ((S33_ratio.ge.S33_good)) 
            !if (ite2 == 0) then
            !    F3mean = 1.d0
            !else
            F3mean = (F3max+F3min)/2.d0
            !end if 
            Fmax(3,3) = F3max
            Fmin(3,3) = F3min
            Fmean(3,3) = F3mean
            call const_stress_foam_ps(Fmax,theta,Smax,nf)
            call const_stress_foam_ps(Fmin,theta,Smin,nf)
            call const_stress_foam_ps(Fmean,theta,Smean,nf)
            if ((Smax(3,3)*Smean(3,3)).le.0.d0) then 
                F3min = F3mean
            else 
                if ((Smin(3,3)*Smean(3,3)).le.0.d0) then 
                    F3max = F3mean 
                else
                    flag = 1
                    exit
                end if
            end if
            S33_ratio = abs(Smean(3,3)/Smean(2,2))
            if (ite2.ge.ite2_max) then
                exit
            end if
            ite2 = ite2+1
            
        end do
        F_out = Fmean
        F33_out = F3mean
        return
    end subroutine solve_plane_stress
     !! newton raphson solver used in plane stress after the first bisection solver, pass F33_out as the guess for its vincinity element
    subroutine nr_solve_ps(F,theta,F33_guess_in,S,F_out,F33_out)
        implicit none
        real(di) F(2,2),theta(14),F33_guess,F33_next,S(3,3),F_out(3,3),F33_out,F33_guess_in
        real(di) S1(3,3),S2(3,3),F33_step,F1(3,3),F2(3,3),step,S33_ratio,conv_ratio,dx
        integer nf,ite,max_ite
        max_ite = 10
        step = 0.01d0
        S1 = 0.d0
        S2 = 0.d0
        F1(1:2,1:2) = F(1:2,1:2)
        F2(1:2,1:2) = F(1:2,1:2)
        ite = 0
        conv_ratio = 0.001d0
        S33_ratio = 1.d0
        F33_guess = F33_guess_in
        do while (ite<max_ite.and.S33_ratio.ge.conv_ratio)
            F1(3,3) = F33_guess
            F2(3,3) = F33_guess*(1.d0+step)
            call const_stress_foam_ps(F1,theta,S1,nf)
            call const_stress_foam_ps(F2,theta,S2,nf)
            dx = S1(3,3)/((S2(3,3)-S1(3,3))/(F33_guess*step))
            S33_ratio = dabs(S1(3,3)/S1(2,2))
            F33_guess = F33_guess-dx
            ite = ite+1
        end do
        F_out = F1
        F33_out = F1(3,3)
        S = S1
        return
    end subroutine nr_solve_ps

    !! The updated version of plane stress constitutive law for foam for homogeneous data only (1st functionality), it adds two checks: monotonic and out of deformation K1 checks
!! this two checks are only applicable on homogeneous field input, when it comes to inhomogeneous field input(2nd functionality), Dacorogna stability check is the only way to go.
    subroutine solve_plane_stress_monok1(F,theta,F_out,Smean,flag,F33_out,S22_out,K1_out)
        implicit none
        real(di) F3max,F3min,F3mean,Fmax(3,3),Fmin(3,3),Fmean(3,3),F(2,2),F33_out,S22_out
        real(di) Smax(3,3),Smin(3,3),Smean(3,3),theta(14),F_out(3,3),S33_ratio,S33_good,K1_out
        integer ite1,ite1_max,ite2,ite2_max,flag,nf
        real(di) RR(3,3),U(3,3),E(3,3),jd
        ! the max and min F33 for simple compression and simple tension        
        F3max = 1.28d0
        F3min = 0.64d0
        call eye(Fmax,3)
        call eye(Fmin,3)
        call eye(Fmean,3)
        Fmax(1:2,1:2) = F(1:2,1:2)
        Fmin(1:2,1:2) = F(1:2,1:2)
        Fmean(1:2,1:2) = F(1:2,1:2)
        ! ite1 max for expanding bracket
        ite1 = 0
        ite1_max = 3
        ! ite2 max for maximum number of iterations
        ite2 = 0
        ite2_max = 10
        ! S33_good is the required ratio between in plane and out of plane stress
        S33_good = 0.001d0
        flag = 0
        Smean = 0.d0
        F_out = 0.d0
        F33_out = 1.d0
        do while (ite1.le.ite1_max)
            if (ite1 == 0) then
                F3mean = 1.d0
            else
                F3mean = (F3max+F3min)/2.d0
            end if 
            
            Fmax(3,3) = F3max
            Fmin(3,3) = F3min
            Fmean(3,3) = F3mean
            call const_stress_foam_ps(Fmax,theta,Smax,nf)
            call const_stress_foam_ps(Fmin,theta,Smin,nf)
            call const_stress_foam_ps(Fmean,theta,Smean,nf)
            if ((Smax(3,3)*Smean(3,3)).le.0) then 
                exit
            else if ((Smin(3,3)*Smean(3,3)).le.0) then 
                exit
            else
                if (ite1 == ite1_max) then
                    flag = 1
                    return
                end if
                F3max = F3max*2.d0
                F3min = F3min/1.5d0
                ite1 = ite1+1
            end if
        end do
        S33_ratio = 1.d0
        do while ((S33_ratio.ge.S33_good)) 
            !if (ite2 == 0) then
            !    F3mean = 1.d0
            !else
            F3mean = (F3max+F3min)/2.d0
            !end if 
            Fmax(3,3) = F3max
            Fmin(3,3) = F3min
            Fmean(3,3) = F3mean
            call const_stress_foam_ps(Fmax,theta,Smax,nf)
            call const_stress_foam_ps(Fmin,theta,Smin,nf)
            call const_stress_foam_ps(Fmean,theta,Smean,nf)
            if ((Smax(3,3)*Smean(3,3)).le.0.d0) then 
                F3min = F3mean
            else 
                if ((Smin(3,3)*Smean(3,3)).le.0.d0) then 
                    F3max = F3mean 
                else
                    flag = 1
                    exit
                end if
            end if
            S33_ratio = abs(Smean(3,3)/Smean(2,2))
            if (ite2.ge.ite2_max) then
                exit
            end if
            ite2 = ite2+1
        end do
        F_out = Fmean
        F33_out = F3mean
        call det3d(F_out,jd)
        S22_out = Smean(2,2)/jd
        call skinem3d(F_out,RR,U,E)
        call trace3d(E,K1_out)
        return
    end subroutine solve_plane_stress_monok1
    !!The alternative method: reduced quadrature implemented for computing force residual is based on an updated lagrangian formulation which requires computation of second order derivative wrt right cauchy green tensor 
    !!would be more useful in history dependent computation. 
    subroutine compute_dphidcc(F,theta,dphidcc_out)
        implicit none
        real(di) F(3,3),E(3,3),K1,devE(3,3),ey(3,3),sumdevE, &
&               K2,N(3,3),K3,Ntemp_det,J,detF, &
&               G0,B,Jmin,C1,K10,delta_K,X1,X2,C0,theta(14), &
&               p,q,C2,C3,r,X,dXdJ1,dfdK1,dLdK2,dLdK3,comp1(3,3), &
&               comp2(3,3),U(3,3),RR(3,3),Y(3,3),dXdK1, &
&               stress(3,3),test,E_temp(3,3),Rot(3,3),dphidcc(3,3,3,3), &
&               dphidcc_1(3,3,3,3),dphidcc_2(3,3,3,3),mcomb(3,3,3,3), &
&               dlncdc(3,3,3,3), dTkdE(3,3,3,3),dTkdEr(3,3,3,3), &
&               Tkr(3,3),zero,dpsidK1,dfdK1_2,dLdK2dK2,dLdK2dK3, &
&               dXdK1_2,dpsidK2,dpsidK3,dpsidK1_2,dpsidK2_2, &
&               dpsidK2dK1,dpsidK2dK3,C(3,3),Cinv(3,3), &
&               dCinvdC(3,3,3,3),dev_E(3,3),Iden(3,3),dphidcc_out(3,3,3,3)
        integer nf,i,jj,k,l,mm,nn,pp,qq
        nf = 14
        zero = 0.d0
        call skinem3d(F,RR,U,E)
        call trace3d(E,K1)
        call eye(ey,3)
        Iden = ey
        devE = E-ey*1.d0/3.d0*K1
        sumdevE = 0.d0
        do i = 1,3
            do jj = 1,3
                sumdevE = sumdevE + devE(i,jj) ** 2.d0
            end do
        end do
        K2 = dsqrt(sumdevE)
        N = devE/K2
        dev_E = devE
        call det3d(N,Ntemp_det)
        K3 = 3.d0*(6.d0 ** 0.5d0)*Ntemp_det
        Y = 3.d0*(6.d0 ** 0.5d0)*matmul(N,N) &
&        -(6.d0 ** 0.5d0)*ey-3.d0*K3*N
        J = 0.d0
        call det3d(F,J)
        call fill(G0,B,Jmin,C1,K10,delta_K,X1,X2,C0,p,q,C2,C3,r,theta)
        X = 0.5d0*(X1+X2)*K1+delta_K/2.d0*(X1-X2)*dlog(dcosh((K1-K10)/delta_K)/dcosh(K10/delta_K))+1.d0
        dXdK1 = 0.5d0*(X1+X2)+0.5d0*(X1-X2)*dtanh((K1-K10)/delta_K)
        if ((J-Jmin)<0.d0) then
            dfdK1 = (dexp(C2*K1)-1.d0)/C2+C3*J*(J**(-r)-(-1)*((1.d0-Jmin) &
&       /(Jmin-J))**r)
        else
            
            dfdK1 = (dexp(C2*K1)-1.d0)/C2+C3*J*(J**(-r)-((1.d0-Jmin) &
&       /(J-Jmin))**r)
        end if 
        dLdK2 = p*C0*K2 ** (p-1.d0)+q*C1*(1.d0+K3)*K2 ** (q-1.d0)
        dLdK3 = C1 * K2 ** q
        dXdK1_2 = (X1-X2)/2.d0/dcosh((K10-K1)/delta_K)**2.d0/delta_K
        dLdK2dK2 = C0*p*(p-1.d0)*K2**(p-2.d0) +  C1*(1.d0+K3)*q*(q-1.d0)*K2**(q-2.d0)
        dLdK2dK3 = C1*q*K2**(q-1.d0)
        dfdK1_2 = dexp(C2*K1) + C3*(J**(-r)-(1.d0-Jmin)**(r)  &
&        /(J-Jmin)**r)*dexp(K1)                              &
&      + C3*(-r*J**(-r-1.d0) + r*(1.d0-Jmin)**(r)/(J-Jmin)**(r+1.d0)) &
&       *dexp(K1)*dexp(K1)
        comp1 = J**(-1.d0)*(ey*G0*dXdK1*K2**2.d0+ &
&       G0*(2.d0*X*K2+dLdK2)*N+G0*dLdK3*Y/K2)
        comp2 = J ** (-1.d0)*(B*dfdK1)*ey
        stress = (comp1+comp2)*J
        
        dpsidK1 = G0*dXdK1*K2**2.d0 + B*dfdK1
        dpsidK2 = G0*(2.d0*X*K2 + dLdK2)
        dpsidK3 = G0*dLdK3
        dpsidK1_2 = G0*dXdK1_2*K2**2.d0 + B*dfdK1_2
        dpsidK2_2 = G0*(2.d0*X + dLdK2dK2)
        dpsidK2dK1 = G0*2.d0*dXdK1*K2
        dpsidK2dK3 = G0*dLdK2dK3
        C = matmul(transpose(F),F)
        call inv3d(C,Cinv) 
        call dlnxdx(C,dlncdc)
        Tkr = matmul(matmul(transpose(RR),stress),RR)
        !print*,dfdK1
        
        if (K2.eq.zero) then
        !
        DTkDE = 0.d0
        do i=1,3
          do jj=1,3
            do k=1,3
              do l=1,3
                DTkDE(i,jj,k,l) = DTkDE(i,jj,k,l)  &
     &              + 2.d0*B*dfdK1_2*Iden(i,jj)*Iden(k,l)  &
     &              + 2.d0*G0*(0.5d0*Iden(i,k)*Iden(jj,l)  &
     &                         + 0.5d0*Iden(i,l)*Iden(jj,k) &
     &                            - 1.d0/3.d0*Iden(i,jj)*Iden(k,l)) 
              end do
            end do
          end do
        end do
        !
      else		 		  

        DTkDE = zero
        do i=1,3
          do jj=1,3
            do k=1,3
              do l=1,3
                DTkDE(i,jj,k,l) = DTkDE(i,jj,k,l) +   &
     &  (dev_E(k,l)*Iden(i,jj) + Iden(k,l)*dev_E(i,jj))*  &
     &  (-2.d0*dsqrt(6.d0)/K2**3.d0*dpsidK3 + 1.d0/K2*dpsidK2dK1)+ &
     &  (dev_E(k,l)*Y(i,jj) + Y(k,l)*dev_E(i,jj))* &
     &  (1.d0/K2**2.d0*dpsidK2dK3 - 3.d0/K2**3.d0*dpsidK3)+ &
     &  dev_E(k,l)*dev_E(i,jj)* &
     &  (-1.d0/K2**3.d0*dpsidK2 + 1.d0/K2**2.d0*dpsidK2_2 - &
     &  3.d0*K3/K2**4.d0*dpsidK3)+ &
     &  (Iden(i,k)*dev_E(l,jj) + Iden(jj,l)*dev_E(i,k) +  &
     &  Iden(i,l)*dev_E(k,jj)+Iden(k,jj)*dev_E(i,l)) &
     &  *3.d0/2.d0*dsqrt(6.d0)/K2**3.d0*dpsidK3   &
     &  +(0.5d0*Iden(i,k)*Iden(jj,l) + 0.5d0*Iden(i,l)*Iden(jj,k) -  &
     &  1.d0/3.d0*Iden(i,jj)*Iden(k,l))* & 
     &  (1.d0/K2*dpsidK2 - 3.d0*K3/K2**2.d0*dpsidK3)+ &
     &  Iden(i,jj)*Iden(k,l)*dpsidK1_2
              end do
            end do
          end do
        end do
        !
      end if
      
      
      dTkdEr = 0.d0
      do i=1,3
        do jj=1,3
          do k=1,3 
            do l=1,3
              do mm=1,3
                do nn=1,3
                  do pp=1,3
                    do qq =1,3				
    dTkdEr(i,jj,k,l) = dTkdEr(i,jj,k,l)+RR(i,mm)*RR(jj,nn)*RR(k,pp)*RR(l,qq)*dTKdE(mm,nn,pp,qq)
                    end do
                  end do
                end do
              end do
            end do
          end do
        end do
      end do
      
      dphidcc = 0.d0
      call mprod4(dTkdEr,dlncdc,mcomb)
      dCinvdC = zero
      do i = 1,3
          do mm = 1,3
              do k = 1,3
                  do l = 1,3
                      dCinvdC(i,k,l,mm) = dCinvdC(i,k,l,mm) + Cinv(i,k)*Cinv(l,mm)
                  end do
              end do
          end do
      end do 
      dphidcc_1 = 0.d0
      dphidcc_2 = 0.d0
      do i=1,3
        do jj=1,3
          do k=1,3 
            do l=1,3
              do mm=1,3
                  dphidcc_1(i,jj,k,l) = dphidcc_1(i,jj,k,l) - 0.5d0*dCinvdC(i,k,l,mm)*Tkr(mm,jj)
                  dphidcc_2(i,jj,k,l) = dphidcc_2(i,jj,k,l) + 0.25d0*Cinv(i,mm)*mcomb(mm,jj,k,l)
              end do
            end do
          end do
        end do
      end do
     dphidcc = 4.d0*(dphidcc_1+dphidcc_2)
     dphidcc_out = 0.d0
      do i=1,3
        do jj=1,3
          do k=1,3 
            do l=1,3
              do mm=1,3
                  do nn = 1,3
                      do pp = 1,3
                          do qq = 1,3
                              dphidcc_out(i,jj,k,l) = dphidcc_out(i,jj,k,l)+ F(i,mm)*F(jj,nn)*F(k,pp)*F(l,qq)*dphidcc(mm,nn,pp,qq)/J
                          end do
                      end do
                  end do
                end do
              end do
            end do
        end do
      end do
      return        
    end subroutine compute_dphidcc
    
end module fem_constitutive
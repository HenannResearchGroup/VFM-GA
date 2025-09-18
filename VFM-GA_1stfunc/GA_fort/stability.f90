!!stability checks, compute acoustic tensor to determine if a point on the load path satisfy stability for the foam hyperelastic material 
!!https://www.sciencedirect.com/science/article/pii/S0022509619303825 see stability section for more detail
!!
!!This stability check takes material parameteres as an input, which is completely independent of the curve(1st func)/field(2nd func) input
!!
!!If you want to apply the stability check to your own material with a given constitutive model, you have to do the derivation of acoustic tensor and implement it here
!!On the other hand side, you could just apply monotonic check and out of plane K1 checks only, which has weaker constraint compared to Dacorogna stability. 
module stability_check
    use constants
    use linalg
    implicit none
    contains
    subroutine sta_pack(theta,staflag)
        implicit none
        real(di) theta_h(14),theta(nparam),K1_temp,K2_temp
        integer staflag,stability,i
        stability = 0
        staflag = 0
        theta_h(1:14) = theta(1:14)
        !print*,theta_h
        call shear_stability(theta,stability)
        if (stability == 1) then
            staflag = 1
        endif
        do i = 1,nstac
            call comp_stability(theta,stability,c_sta_k1(i),c_sta_K2(i))
            if (stability == 1) then 
                staflag = 1
            endif
        end do
        return
    end subroutine sta_pack
    
    
    !! check given a point described by K1 and K2 on a uniaxial compression load pass if material behavior is stable 
    subroutine comp_stability(theta,stability,K1,K2)
        implicit none
        real(di) theta(nparam),G0,B,Jmin,C1,X1p,X2p,dK1,K10,C0,p,q,C2,C3
        real(di) delta(3,3),N(3,3),K3,Y(3,3),detN,E(3,3),eigval_vec(3),eigvec(3,3)
        real(di) J,p_stretch(3,3),Bstretch(3,3),Fdefgrad(3,3),eigval(3,3),Finv(3,3)
        real(di) X,dXdK1,d2XdK12,dLdK2,dLdK3,d2LdK22,d2LdK2dK3,dfdK1,d2fdK12
        real(di) dpsidK1,dpsidK2,dpsidK3,TK(3,3),d2psidK12,d2psidK22,d2psidK32
        real(di) d2psidK1dK2,d2psidK1dK3,d2psidK2dK3,Dtang(3,3,3,3),Ltan(3,3,3,3)
        real(di) Btan(3,3,3,3),Ctang(3,3,3,3),DLtan(3,3,3,3),Cmat(3,3,3,3),temp
        real(di) K1,K2
        integer mm,nn,pp,qq,rr,ss,stability,deltalist(3),tempM(3,3)
        stability= 0
        G0 = theta(1)/1000.d0
        B = theta(2)/1000.d0
        Jmin = theta(3)
        C1 = theta(4)
        K10 = theta(5)
        dK1 = theta(6)
        X1p = theta(7)
        X2p = theta(8)
        C0 = theta(9)
        p = theta(10)
        q = theta(11)
        C2 = theta(12)
        C3 = theta(13) 
        call eye(delta,3)
        N = 0.d0
        N(1,1) = 1.d0/(6.d0 ** 0.5d0)
        N(2,2) = 1.d0/(6.d0 ** 0.5d0)
        N(3,3) = -2.d0/(6.d0 ** 0.5d0)
        call det3d(N,detN)
        K3 = 3.d0*(6.d0 ** 0.5d0)*detN
        Y = 3.d0*(6.d0 ** 0.5d0)*matmul(N,N) &
&        -(6.d0 ** 0.5d0)*delta-3.d0*K3*N
        !! seperation
        E = (1.d0/3.d0)*K1*delta + K2*N
        J = dexp(K1)
        call spectral3d(E,eigval_vec,eigvec)
        eigval = 0.d0
        eigval(1,1) = eigval_vec(1)
        eigval(2,2) = eigval_vec(2)
        eigval(3,3) = eigval_vec(3)
        p_stretch = 0.d0
        do mm = 1,3
            eigval(mm,mm) = dexp(2.d0*eigval(mm,mm))
            p_stretch(mm,mm) = eigval(mm,mm)**0.5d0
        end do
        Bstretch = matmul(matmul(eigvec,eigval),transpose(eigvec))
        Fdefgrad = matmul(matmul(eigvec,p_stretch),transpose(eigvec))
        call inv3d(Fdefgrad,Finv)
        X = 0.5d0*(X1p+X2p)*K1+dK1/2.d0*(X1p-X2p)*dlog(dcosh((K1-K10)/dK1)/dcosh(K10/dK1))+1.d0
        dXdK1 = 0.5d0*(X1p+X2p)+0.5d0*(X1p-X2p)*dtanh((K1-K10)/dK1)
        d2XdK12 = (X1p-X2p)/2.d0/dcosh((K10-K1)/dK1)**2.d0/dK1
        
        dLdK2 = p*C0*K2 ** (p-1.d0)+q*C1*(1.d0+K3)*K2 ** (q-1.d0)
        dLdK3 = C1 * K2 ** q
        d2LdK22 = p*(p-1.d0)*C0*(K2**(p-2.d0)) + q*(q-1.d0)*C1*(1.d0+K3)*(K2**(q-2.d0))
        d2LdK2dK3 = q*C1*(K2**(q-1.d0))  
        
        dfdK1 = (dexp(C2*K1)-1.d0)/C2+C3*dexp(K1)*(dexp(-2.d0*K1)-((1.d0-Jmin)**2.d0)/(dexp(K1)-Jmin)**2.d0)
        d2fdK12 = dexp(C2*K1) + C3*dexp(K1)*(dexp(-2.d0*K1) - ((1.d0-Jmin)**2.d0)/((dexp(K1)-Jmin)**2.d0)) + &
&             C3*dexp(K1)*(-2.d0*dexp(-2.d0*K1) + 2.d0*dexp(K1)*((1.d0-Jmin)**2.d0)/((dexp(K1)-Jmin)**3.d0))
        dpsidK1 = G0*(dXdK1*K2**2.d0) + B*dfdK1
        dpsidK2 = G0*(2.d0*X*K2 + dLdK2)
        dpsidK3 = G0*dLdK3
        d2psidK12 = G0*d2XdK12*K2**2.d0 + B*d2fdK12
        d2psidK22 = G0*(2.d0*X + d2LdK22)
        d2psidK32 = 0.d0
        d2psidK1dK2 = G0*2.d0*dXdK1*K2
        d2psidK1dK3 = 0.d0
        d2psidK2dK3 = G0*d2LdK2dK3
        TK = dpsidK1*delta + dpsidK2*N + dpsidK3*(1.d0/K2)*Y
        Dtang = 0.d0
        do mm = 1,3
            do nn = 1,3
                do pp = 1,3
                    do qq = 1,3
                        Dtang(mm,nn,pp,qq) = Dtang(mm,nn,pp,qq) + d2psidK12*delta(mm,nn)*delta(pp,qq) + &
&                       (d2psidK22 - dpsidK2*(1.d0/K2) - dpsidK3*(3.d0*K3/(K2**2)))*N(mm,nn)*N(pp,qq) + &
&                        (d2psidK32/(K2**2.d0))*Y(mm,nn)*Y(pp,qq) + (d2psidK1dK2 - 2.d0*dsqrt(6.d0)*(dpsidK3/(K2**2.d0))) &
&                        *(delta(mm,nn)*N(pp,qq) + N(mm,nn)*delta(pp,qq)) + d2psidK1dK3*(1.d0/K2)*(delta(mm,nn)*Y(pp,qq) &
&                           + Y(mm,nn)*delta(pp,qq)) + (d2psidK2dK3*(1.d0/K2) - 3.d0*dpsidK3*(1.d0/(K2**2.d0)))* &
&                           (N(mm,nn)*Y(pp,qq) + Y(mm,nn)*N(pp,qq)) + (dpsidK2*(1.d0/K2) - dpsidK3*(3.d0*K3/(K2**2.d0))) &
&                           *((1.d0/2.d0)*delta(mm,pp)*delta(nn,qq) + (1.d0/2.d0)*delta(mm,qq)*delta(nn,pp) - &
&                           (1.d0/3.d0)*delta(mm,nn)*delta(pp,qq)) + (3.d0*dsqrt(6.d0)/2.d0)*dpsidK3*(1.d0/(K2**2.d0)) &
&                           *(delta(mm,pp)*N(nn,qq) & 
&                           + delta(mm,qq)*N(nn,pp) + N(mm,pp)*delta(nn,qq) + N(mm,qq)*delta(nn,pp))
                    end do
                end do
            end do
        end do
        call dlnxdx(Bstretch,Ltan)
        Btan = 0.d0
        do mm = 1,3
            do nn = 1,3
                do pp = 1,3
                    do qq = 1,3
                        Btan(mm,nn,pp,qq) = Btan(mm,nn,pp,qq) + delta(mm,pp)*Bstretch(nn,qq) + delta(nn,pp)*Bstretch(mm,qq)
                    end do
                end do
            end do 
        end do
        call mprod4(Dtang,Ltan,DLtan)
        call mprod4(DLtan,Btan,Ctang)
        Ctang = Ctang/(2.d0*J)
        do mm = 1,3
            do nn = 1,3
                do pp = 1,3
                    do qq = 1,3
                        Ctang(mm,nn,pp,qq) = Ctang(mm,nn,pp,qq) - TK(mm,qq)*delta(nn,pp)/J
                    end do
                end do
            end do 
        end do
        Cmat = 0.d0
        do mm = 1,3
            do nn = 1,3
                do pp = 1,3
                    do qq = 1,3
                        do rr = 1,3
                            do ss = 1,3
                                Cmat(mm,nn,pp,qq) = Cmat(mm,nn,pp,qq) + J*Finv(nn,rr)*Finv(qq,ss)*Ctang(mm,rr,pp,ss);
                            end do
                        end do
                    end do
                end do
            end do
        end do
        
        if (Cmat(1,1,1,1)<=0.d0) then
            stability = 1
        end if
        
        if (Cmat(2,2,2,2)<=0.d0) then
            stability = 1
        end if
        
        if (Cmat(3,3,3,3)<=0.d0) then
            stability = 1
        end if
        
        if (Cmat(1,2,1,2)<=0.d0) then
            stability = 1
        end if
        
        if (Cmat(1,3,1,3)<=0.d0) then
            stability = 1
        end if
        
        if (Cmat(2,3,2,3)<=0.d0) then
            stability = 1
        end if
        
        temp = Cmat(1,1,1,1)*Cmat(2,2,2,2) + Cmat(1,2,1,2)*Cmat(1,2,1,2) - &
&       (Cmat(1,1,2,2) + Cmat(1,2,2,1))**2.d0 + 2.d0*Cmat(1,2,1,2)*dsqrt(Cmat(1,1,1,1)*Cmat(2,2,2,2))
        if (temp<=0.d0) then
            stability = 1
        end if
        
        temp = Cmat(1,1,1,1)*Cmat(3,3,3,3) + Cmat(1,3,1,3)*Cmat(1,3,1,3) - &
&       (Cmat(1,1,3,3) + Cmat(1,3,3,1))**2.d0 + 2.d0*Cmat(1,3,1,3)*dsqrt(Cmat(1,1,1,1)*Cmat(3,3,3,3))
        if (temp<=0.d0) then
            stability = 1
        end if
        
        temp = Cmat(2,2,2,2)*Cmat(3,3,3,3) + Cmat(2,3,2,3)*Cmat(2,3,2,3) - &
&       (Cmat(2,2,3,3) + Cmat(2,3,3,2))**2.d0 + 2.d0*Cmat(2,3,2,3)*dsqrt(Cmat(2,2,2,2)*Cmat(3,3,3,3))
        if (temp<=0.d0) then
            stability = 1
        end if
    
        deltalist = 1.d0
        tempM  = 0.d0
        tempM(1,1) = Cmat(1,1,1,1)
        tempM(1,2) = Cmat(1,2,1,2) + deltalist(1)*deltalist(2)*(Cmat(1,1,2,2) + Cmat(1,2,2,1))
        tempM(1,3) = Cmat(1,3,1,3) + deltalist(1)*deltalist(3)*(Cmat(1,1,3,3) + Cmat(1,3,3,1))
        tempM(2,1) = Cmat(2,1,2,1) + deltalist(2)*deltalist(1)*(Cmat(2,2,1,1) + Cmat(2,1,1,2))
        tempM(2,2) = Cmat(2,2,2,2)
        tempM(2,3) = Cmat(2,3,2,3) + deltalist(2)*deltalist(3)*(Cmat(2,2,3,3) + Cmat(2,3,3,2))
        tempM(3,1) = Cmat(3,1,3,1) + deltalist(3)*deltalist(1)*(Cmat(3,3,1,1) + Cmat(3,1,1,3))
        tempM(3,2) = Cmat(3,2,3,2) + deltalist(3)*deltalist(2)*(Cmat(3,3,2,2) + Cmat(3,2,2,3))
        tempM(3,3) = Cmat(3,3,3,3)
        temp = tempM(1,2)*dsqrt(Cmat(3,3,3,3)) + tempM(1,3)*dsqrt(Cmat(2,2,2,2)) + tempM(2,3)*dsqrt(Cmat(1,1,1,1)) + dsqrt(Cmat(1,1,1,1)*Cmat(2,2,2,2)*Cmat(3,3,3,3))
        if (temp<0.d0) then
            stability = 1
        end if 
        
        deltalist = (/1.d0,1.d0,-1.d0/)
        tempM = 0.d0
        tempM(1,1) = Cmat(1,1,1,1)
        tempM(1,2) = Cmat(1,2,1,2) + deltalist(1)*deltalist(2)*(Cmat(1,1,2,2) + Cmat(1,2,2,1))
        tempM(1,3) = Cmat(1,3,1,3) + deltalist(1)*deltalist(3)*(Cmat(1,1,3,3) + Cmat(1,3,3,1))
        tempM(2,1) = Cmat(2,1,2,1) + deltalist(2)*deltalist(1)*(Cmat(2,2,1,1) + Cmat(2,1,1,2))
        tempM(2,2) = Cmat(2,2,2,2)
        tempM(2,3) = Cmat(2,3,2,3) + deltalist(2)*deltalist(3)*(Cmat(2,2,3,3) + Cmat(2,3,3,2))
        tempM(3,1) = Cmat(3,1,3,1) + deltalist(3)*deltalist(1)*(Cmat(3,3,1,1) + Cmat(3,1,1,3))
        tempM(3,2) = Cmat(3,2,3,2) + deltalist(3)*deltalist(2)*(Cmat(3,3,2,2) + Cmat(3,2,2,3))
        tempM(3,3) = Cmat(3,3,3,3)
        temp = tempM(1,2)*dsqrt(Cmat(3,3,3,3)) + tempM(1,3)*dsqrt(Cmat(2,2,2,2)) + tempM(2,3)*dsqrt(Cmat(1,1,1,1)) + dsqrt(Cmat(1,1,1,1)*Cmat(2,2,2,2)*Cmat(3,3,3,3))
        if(temp<0.d0) then
            stability = 1
        end if
    
        deltalist = (/1.d0,-1.d0,1.d0/)
        tempM = 0.d0
        tempM(1,1) = Cmat(1,1,1,1)
        tempM(1,2) = Cmat(1,2,1,2) + deltalist(1)*deltalist(2)*(Cmat(1,1,2,2) + Cmat(1,2,2,1))
        tempM(1,3) = Cmat(1,3,1,3) + deltalist(1)*deltalist(3)*(Cmat(1,1,3,3) + Cmat(1,3,3,1))
        tempM(2,1) = Cmat(2,1,2,1) + deltalist(2)*deltalist(1)*(Cmat(2,2,1,1) + Cmat(2,1,1,2))
        tempM(2,2) = Cmat(2,2,2,2)
        tempM(2,3) = Cmat(2,3,2,3) + deltalist(2)*deltalist(3)*(Cmat(2,2,3,3) + Cmat(2,3,3,2))
        tempM(3,1) = Cmat(3,1,3,1) + deltalist(3)*deltalist(1)*(Cmat(3,3,1,1) + Cmat(3,1,1,3))
        tempM(3,2) = Cmat(3,2,3,2) + deltalist(3)*deltalist(2)*(Cmat(3,3,2,2) + Cmat(3,2,2,3))
        tempM(3,3) = Cmat(3,3,3,3)
        temp = tempM(1,2)*dsqrt(Cmat(3,3,3,3)) + tempM(1,3)*dsqrt(Cmat(2,2,2,2)) + tempM(2,3)*dsqrt(Cmat(1,1,1,1)) + dsqrt(Cmat(1,1,1,1)*Cmat(2,2,2,2)*Cmat(3,3,3,3))
        if(temp<0.d0) then
            stability = 1
        end if
        deltalist = (/-1.d0,1.d0,1.d0/)
        tempM = 0.d0
        tempM(1,1) = Cmat(1,1,1,1)
        tempM(1,2) = Cmat(1,2,1,2) + deltalist(1)*deltalist(2)*(Cmat(1,1,2,2) + Cmat(1,2,2,1))
        tempM(1,3) = Cmat(1,3,1,3) + deltalist(1)*deltalist(3)*(Cmat(1,1,3,3) + Cmat(1,3,3,1))
        tempM(2,1) = Cmat(2,1,2,1) + deltalist(2)*deltalist(1)*(Cmat(2,2,1,1) + Cmat(2,1,1,2))
        tempM(2,2) = Cmat(2,2,2,2)
        tempM(2,3) = Cmat(2,3,2,3) + deltalist(2)*deltalist(3)*(Cmat(2,2,3,3) + Cmat(2,3,3,2))
        tempM(3,1) = Cmat(3,1,3,1) + deltalist(3)*deltalist(1)*(Cmat(3,3,1,1) + Cmat(3,1,1,3))
        tempM(3,2) = Cmat(3,2,3,2) + deltalist(3)*deltalist(2)*(Cmat(3,3,2,2) + Cmat(3,2,2,3))
        tempM(3,3) = Cmat(3,3,3,3)
        temp = tempM(1,2)*dsqrt(Cmat(3,3,3,3)) + tempM(1,3)*dsqrt(Cmat(2,2,2,2)) + tempM(2,3)*dsqrt(Cmat(1,1,1,1)) + dsqrt(Cmat(1,1,1,1)*Cmat(2,2,2,2)*Cmat(3,3,3,3))
        if(temp<0.d0) then
            stability = 1
        end if
        return
    end subroutine comp_stability
    
    !! similar to the previous one, checks on the pure shear load pass, if the material behavior is stable
    subroutine shear_stability(theta,stability)
        implicit none
        real(di) theta(nparam),G0,B,Jmin,C1,X1p,X2p,dK1,K10,C0,p,q,C2,C3
        real(di) delta(3,3),N(3,3),K3,K1,K2,Y(3,3),detN,E(3,3),eigval_vec(3),eigvec(3,3)
        real(di) J,p_stretch(3,3),Bstretch(3,3),Fdefgrad(3,3),eigval(3,3),Finv(3,3)
        real(di) X,dXdK1,d2XdK12,dLdK2,dLdK3,d2LdK22,d2LdK2dK3,dfdK1,d2fdK12
        real(di) dpsidK1,dpsidK2,dpsidK3,TK(3,3),d2psidK12,d2psidK22,d2psidK32
        real(di) d2psidK1dK2,d2psidK1dK3,d2psidK2dK3,Dtang(3,3,3,3),Ltan(3,3,3,3)
        real(di) Btan(3,3,3,3),Ctang(3,3,3,3),DLtan(3,3,3,3),Cmat(3,3,3,3),temp
        integer mm,nn,pp,qq,rr,ss,stability,deltalist(3),tempM(3,3)
        K1 = 0.d0
        K2 = s_sta
        stability= 0
        G0 = theta(1)/1000
        B = theta(2)/1000
        Jmin = theta(3)
        C1 = theta(4)
        K10 = theta(5)
        dK1 = theta(6)
        X1p = theta(7)
        X2p = theta(8)
        C0 = theta(9)
        p = theta(10)
        q = theta(11)
        C2 = theta(12)
        C3 = theta(13) 
        call eye(delta,3)
        N = 0.d0
        N(1,1) = 1.d0/(2.d0 ** 0.5d0)
        N(3,3) = -1.d0/(2.d0 ** 0.5d0)
        call det3d(N,detN)
        K3 = 3.d0*(6.d0 ** 0.5d0)*detN
        Y = 3.d0*(6.d0 ** 0.5d0)*matmul(N,N) &
&        -(6.d0 ** 0.5d0)*delta-3.d0*K3*N
        !! seperation
        E = (1.d0/3.d0)*K1*delta + K2*N
        J = dexp(K1)
        call spectral3d(E,eigval_vec,eigvec)
        eigval = 0.d0
        eigval(1,1) = eigval_vec(1)
        eigval(2,2) = eigval_vec(2)
        eigval(3,3) = eigval_vec(3)
        p_stretch = 0.d0
        do mm = 1,3
            eigval(mm,mm) = dexp(2.d0*eigval(mm,mm))
            p_stretch(mm,mm) = eigval(mm,mm)**0.5d0
        end do
        Bstretch = matmul(matmul(eigvec,eigval),transpose(eigvec))
        Fdefgrad = matmul(matmul(eigvec,p_stretch),transpose(eigvec))
        call inv3d(Fdefgrad,Finv)
        X = 0.5d0*(X1p+X2p)*K1+dK1/2.d0*(X1p-X2p)*dlog(dcosh((K1-K10)/dK1)/dcosh(K10/dK1))+1.d0
        dXdK1 = 0.5d0*(X1p+X2p)+0.5d0*(X1p-X2p)*dtanh((K1-K10)/dK1)
        d2XdK12 = (X1p-X2p)/2.d0/dcosh((K10-K1)/dK1)**2.d0/dK1
        
        dLdK2 = p*C0*K2 ** (p-1.d0)+q*C1*(1.d0+K3)*K2 ** (q-1.d0)
        dLdK3 = C1 * K2 ** q
        d2LdK22 = p*(p-1.d0)*C0*(K2**(p-2.d0)) + q*(q-1.d0)*C1*(1.d0+K3)*(K2**(q-2.d0))
        d2LdK2dK3 = q*C1*(K2**(q-1.d0))  
        
        dfdK1 = (dexp(C2*K1)-1.d0)/C2+C3*dexp(K1)*(dexp(-2.d0*K1)-((1.d0-Jmin)**2.d0)/(dexp(K1)-Jmin)**2.d0)
        d2fdK12 = dexp(C2*K1) + C3*dexp(K1)*(dexp(-2.d0*K1) - ((1.d0-Jmin)**2.d0)/((dexp(K1)-Jmin)**2.d0)) + &
&             C3*dexp(K1)*(-2.d0*dexp(-2.d0*K1) + 2.d0*dexp(K1)*((1.d0-Jmin)**2.d0)/((dexp(K1)-Jmin)**3.d0))
        dpsidK1 = G0*(dXdK1*K2**2.d0) + B*dfdK1
        dpsidK2 = G0*(2.d0*X*K2 + dLdK2)
        dpsidK3 = G0*dLdK3
        d2psidK12 = G0*d2XdK12*K2**2.d0 + B*d2fdK12
        d2psidK22 = G0*(2.d0*X + d2LdK22)
        d2psidK32 = 0.d0
        d2psidK1dK2 = G0*2.d0*dXdK1*K2
        d2psidK1dK3 = 0.d0
        d2psidK2dK3 = G0*d2LdK2dK3
        TK = dpsidK1*delta + dpsidK2*N + dpsidK3*(1.d0/K2)*Y
        Dtang = 0.d0
        do mm = 1,3
            do nn = 1,3
                do pp = 1,3
                    do qq = 1,3
                        Dtang(mm,nn,pp,qq) = Dtang(mm,nn,pp,qq) + d2psidK12*delta(mm,nn)*delta(pp,qq) + &
&                       (d2psidK22 - dpsidK2*(1.d0/K2) - dpsidK3*(3.d0*K3/(K2**2)))*N(mm,nn)*N(pp,qq) + &
&                        (d2psidK32/(K2**2.d0))*Y(mm,nn)*Y(pp,qq) + (d2psidK1dK2 - 2.d0*dsqrt(6.d0)*(dpsidK3/(K2**2.d0))) &
&                        *(delta(mm,nn)*N(pp,qq) + N(mm,nn)*delta(pp,qq)) + d2psidK1dK3*(1.d0/K2)*(delta(mm,nn)*Y(pp,qq) &
&                           + Y(mm,nn)*delta(pp,qq)) + (d2psidK2dK3*(1.d0/K2) - 3.d0*dpsidK3*(1.d0/(K2**2.d0)))* &
&                           (N(mm,nn)*Y(pp,qq) + Y(mm,nn)*N(pp,qq)) + (dpsidK2*(1.d0/K2) - dpsidK3*(3.d0*K3/(K2**2.d0))) &
&                           *((1.d0/2.d0)*delta(mm,pp)*delta(nn,qq) + (1.d0/2.d0)*delta(mm,qq)*delta(nn,pp) - &
&                           (1.d0/3.d0)*delta(mm,nn)*delta(pp,qq)) + (3.d0*dsqrt(6.d0)/2.d0)*dpsidK3*(1.d0/(K2**2.d0)) &
&                           *(delta(mm,pp)*N(nn,qq) & 
&                           + delta(mm,qq)*N(nn,pp) + N(mm,pp)*delta(nn,qq) + N(mm,qq)*delta(nn,pp))
                    end do
                end do
            end do
        end do
        call dlnxdx(Bstretch,Ltan)
        Btan = 0.d0
        do mm = 1,3
            do nn = 1,3
                do pp = 1,3
                    do qq = 1,3
                        Btan(mm,nn,pp,qq) = Btan(mm,nn,pp,qq) + delta(mm,pp)*Bstretch(nn,qq) + delta(nn,pp)*Bstretch(mm,qq)
                    end do
                end do
            end do 
        end do
        call mprod4(Dtang,Ltan,DLtan)
        call mprod4(DLtan,Btan,Ctang)
        Ctang = Ctang/(2.d0*J)
        do mm = 1,3
            do nn = 1,3
                do pp = 1,3
                    do qq = 1,3
                        Ctang(mm,nn,pp,qq) = Ctang(mm,nn,pp,qq) - TK(mm,qq)*delta(nn,pp)/J
                    end do
                end do
            end do 
        end do
        Cmat = 0.d0
        do mm = 1,3
            do nn = 1,3
                do pp = 1,3
                    do qq = 1,3
                        do rr = 1,3
                            do ss = 1,3
                                Cmat(mm,nn,pp,qq) = Cmat(mm,nn,pp,qq) + J*Finv(nn,rr)*Finv(qq,ss)*Ctang(mm,rr,pp,ss);
                            end do
                        end do
                    end do
                end do
            end do
        end do
        
        if (Cmat(1,1,1,1)<=0.d0) then
            stability = 1
        end if
        
        if (Cmat(2,2,2,2)<=0.d0) then
            stability = 1
        end if
        
        if (Cmat(3,3,3,3)<=0.d0) then
            stability = 1
        end if
        
        if (Cmat(1,2,1,2)<=0.d0) then
            stability = 1
        end if
        
        if (Cmat(1,3,1,3)<=0.d0) then
            stability = 1
        end if
        
        if (Cmat(2,3,2,3)<=0.d0) then
            stability = 1
        end if
        
        temp = Cmat(1,1,1,1)*Cmat(2,2,2,2) + Cmat(1,2,1,2)*Cmat(1,2,1,2) - &
&       (Cmat(1,1,2,2) + Cmat(1,2,2,1))**2.d0 + 2.d0*Cmat(1,2,1,2)*dsqrt(Cmat(1,1,1,1)*Cmat(2,2,2,2))
        if (temp<=0.d0) then
            stability = 1
        end if
        
        temp = Cmat(1,1,1,1)*Cmat(3,3,3,3) + Cmat(1,3,1,3)*Cmat(1,3,1,3) - &
&       (Cmat(1,1,3,3) + Cmat(1,3,3,1))**2.d0 + 2.d0*Cmat(1,3,1,3)*dsqrt(Cmat(1,1,1,1)*Cmat(3,3,3,3))
        if (temp<=0.d0) then
            stability = 1
        end if
        
        temp = Cmat(2,2,2,2)*Cmat(3,3,3,3) + Cmat(2,3,2,3)*Cmat(2,3,2,3) - &
&       (Cmat(2,2,3,3) + Cmat(2,3,3,2))**2.d0 + 2.d0*Cmat(2,3,2,3)*dsqrt(Cmat(2,2,2,2)*Cmat(3,3,3,3))
        if (temp<=0.d0) then
            stability = 1
        end if
    
        deltalist = 1.d0
        tempM  = 0.d0
        tempM(1,1) = Cmat(1,1,1,1)
        tempM(1,2) = Cmat(1,2,1,2) + deltalist(1)*deltalist(2)*(Cmat(1,1,2,2) + Cmat(1,2,2,1))
        tempM(1,3) = Cmat(1,3,1,3) + deltalist(1)*deltalist(3)*(Cmat(1,1,3,3) + Cmat(1,3,3,1))
        tempM(2,1) = Cmat(2,1,2,1) + deltalist(2)*deltalist(1)*(Cmat(2,2,1,1) + Cmat(2,1,1,2))
        tempM(2,2) = Cmat(2,2,2,2)
        tempM(2,3) = Cmat(2,3,2,3) + deltalist(2)*deltalist(3)*(Cmat(2,2,3,3) + Cmat(2,3,3,2))
        tempM(3,1) = Cmat(3,1,3,1) + deltalist(3)*deltalist(1)*(Cmat(3,3,1,1) + Cmat(3,1,1,3))
        tempM(3,2) = Cmat(3,2,3,2) + deltalist(3)*deltalist(2)*(Cmat(3,3,2,2) + Cmat(3,2,2,3))
        tempM(3,3) = Cmat(3,3,3,3)
        temp = tempM(1,2)*dsqrt(Cmat(3,3,3,3)) + tempM(1,3)*dsqrt(Cmat(2,2,2,2)) + tempM(2,3)*dsqrt(Cmat(1,1,1,1)) + dsqrt(Cmat(1,1,1,1)*Cmat(2,2,2,2)*Cmat(3,3,3,3))
        if (temp<0.d0) then
            stability = 1
        end if 
        
        deltalist = (/1.d0,1.d0,-1.d0/)
        tempM = 0.d0
        tempM(1,1) = Cmat(1,1,1,1)
        tempM(1,2) = Cmat(1,2,1,2) + deltalist(1)*deltalist(2)*(Cmat(1,1,2,2) + Cmat(1,2,2,1))
        tempM(1,3) = Cmat(1,3,1,3) + deltalist(1)*deltalist(3)*(Cmat(1,1,3,3) + Cmat(1,3,3,1))
        tempM(2,1) = Cmat(2,1,2,1) + deltalist(2)*deltalist(1)*(Cmat(2,2,1,1) + Cmat(2,1,1,2))
        tempM(2,2) = Cmat(2,2,2,2)
        tempM(2,3) = Cmat(2,3,2,3) + deltalist(2)*deltalist(3)*(Cmat(2,2,3,3) + Cmat(2,3,3,2))
        tempM(3,1) = Cmat(3,1,3,1) + deltalist(3)*deltalist(1)*(Cmat(3,3,1,1) + Cmat(3,1,1,3))
        tempM(3,2) = Cmat(3,2,3,2) + deltalist(3)*deltalist(2)*(Cmat(3,3,2,2) + Cmat(3,2,2,3))
        tempM(3,3) = Cmat(3,3,3,3)
        temp = tempM(1,2)*dsqrt(Cmat(3,3,3,3)) + tempM(1,3)*dsqrt(Cmat(2,2,2,2)) + tempM(2,3)*dsqrt(Cmat(1,1,1,1)) + dsqrt(Cmat(1,1,1,1)*Cmat(2,2,2,2)*Cmat(3,3,3,3))
        if(temp<0.d0) then
            stability = 1
        end if
    
        deltalist = (/1.d0,-1.d0,1.d0/)
        tempM = 0.d0
        tempM(1,1) = Cmat(1,1,1,1)
        tempM(1,2) = Cmat(1,2,1,2) + deltalist(1)*deltalist(2)*(Cmat(1,1,2,2) + Cmat(1,2,2,1))
        tempM(1,3) = Cmat(1,3,1,3) + deltalist(1)*deltalist(3)*(Cmat(1,1,3,3) + Cmat(1,3,3,1))
        tempM(2,1) = Cmat(2,1,2,1) + deltalist(2)*deltalist(1)*(Cmat(2,2,1,1) + Cmat(2,1,1,2))
        tempM(2,2) = Cmat(2,2,2,2)
        tempM(2,3) = Cmat(2,3,2,3) + deltalist(2)*deltalist(3)*(Cmat(2,2,3,3) + Cmat(2,3,3,2))
        tempM(3,1) = Cmat(3,1,3,1) + deltalist(3)*deltalist(1)*(Cmat(3,3,1,1) + Cmat(3,1,1,3))
        tempM(3,2) = Cmat(3,2,3,2) + deltalist(3)*deltalist(2)*(Cmat(3,3,2,2) + Cmat(3,2,2,3))
        tempM(3,3) = Cmat(3,3,3,3)
        temp = tempM(1,2)*dsqrt(Cmat(3,3,3,3)) + tempM(1,3)*dsqrt(Cmat(2,2,2,2)) + tempM(2,3)*dsqrt(Cmat(1,1,1,1)) + dsqrt(Cmat(1,1,1,1)*Cmat(2,2,2,2)*Cmat(3,3,3,3))
        if(temp<0.d0) then
            stability = 1
        end if
        deltalist = (/-1.d0,1.d0,1.d0/)
        tempM = 0.d0
        tempM(1,1) = Cmat(1,1,1,1)
        tempM(1,2) = Cmat(1,2,1,2) + deltalist(1)*deltalist(2)*(Cmat(1,1,2,2) + Cmat(1,2,2,1))
        tempM(1,3) = Cmat(1,3,1,3) + deltalist(1)*deltalist(3)*(Cmat(1,1,3,3) + Cmat(1,3,3,1))
        tempM(2,1) = Cmat(2,1,2,1) + deltalist(2)*deltalist(1)*(Cmat(2,2,1,1) + Cmat(2,1,1,2))
        tempM(2,2) = Cmat(2,2,2,2)
        tempM(2,3) = Cmat(2,3,2,3) + deltalist(2)*deltalist(3)*(Cmat(2,2,3,3) + Cmat(2,3,3,2))
        tempM(3,1) = Cmat(3,1,3,1) + deltalist(3)*deltalist(1)*(Cmat(3,3,1,1) + Cmat(3,1,1,3))
        tempM(3,2) = Cmat(3,2,3,2) + deltalist(3)*deltalist(2)*(Cmat(3,3,2,2) + Cmat(3,2,2,3))
        tempM(3,3) = Cmat(3,3,3,3)
        temp = tempM(1,2)*dsqrt(Cmat(3,3,3,3)) + tempM(1,3)*dsqrt(Cmat(2,2,2,2)) + tempM(2,3)*dsqrt(Cmat(1,1,1,1)) + dsqrt(Cmat(1,1,1,1)*Cmat(2,2,2,2)*Cmat(3,3,3,3))
        if(temp<0.d0) then
            stability = 1
        end if
        return
    end subroutine shear_stability
    
    subroutine dlnxdx(X,DYDX)
      !
      ! This subroutine calculates the derivative of the logarithm
      ! of a symmetric tensor with respect to that tensor
      !
      implicit none
      !
      integer i,j,k,l
      !
      real(di) X(3,3),DYDX(3,3,3,3),Iden(3,3),Iden4(3,3,3,3),eigval(3)
      real(di)  eigvec(3,3),ehat1(3),ehat2(3),ehat3(3),E1(3,3),E2(3,3),E3(3,3)
      real(di) y(3),DX2DX(3,3,3,3),s1,s2,s3,s4,s5,s6
      !
      real(di) zero,one,two,half,three,third,small
      parameter(zero=0.d0,one=1.d0,two=2.d0,half=0.5d0,  &
&              three=3.d0,third=1.d0/3.d0,small=1.d-12)
      
      
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
              Iden4(i,j,k,l) = half*(Iden(i,k)*Iden(j,l)+Iden(i,l)*Iden(j,k))
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
     &                                X(i,k)*Iden(j,l) + &
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
        s3 = two*(y(1)-y(2))/((eigval(1)-eigval(2))**three) - &
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
     &                     s3*X(i,j)*X(k,l) + s4*X(i,j)*Iden(k,l) + &
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
                DYDX(i,j,k,l) = (y(1)/((eigval(1)-eigval(2))*   &
     &                                 (eigval(1)-eigval(3))))*  &
     &         (DX2DX(i,j,k,l) - (eigval(2)+eigval(3))*Iden4(i,j,k,l) -  &
     &  ((eigval(1)-eigval(2))+(eigval(1)-eigval(3)))*E1(i,j)*E1(k,l) -  &
     &       (eigval(2)-eigval(3))*(E2(i,j)*E2(k,l)-E3(i,j)*E3(k,l))) + &
     &                        (one/eigval(1))*E1(i,j)*E1(k,l) +  &
     &                          (y(2)/((eigval(2)-eigval(1))* &
     &                                 (eigval(2)-eigval(3))))* &
     &         (DX2DX(i,j,k,l) - (eigval(1)+eigval(3))*Iden4(i,j,k,l) -  &
     &  ((eigval(2)-eigval(1))+(eigval(2)-eigval(3)))*E2(i,j)*E2(k,l) - &
     &       (eigval(1)-eigval(3))*(E1(i,j)*E1(k,l)-E3(i,j)*E3(k,l))) + &
     &                        (one/eigval(2))*E2(i,j)*E2(k,l) + &
     &                          (y(3)/((eigval(3)-eigval(1))* &
     &                                 (eigval(3)-eigval(2))))* & 
     &         (DX2DX(i,j,k,l) - (eigval(1)+eigval(2))*Iden4(i,j,k,l) - & 
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
    subroutine mprod4(A,B,C)
      !
      ! This subroutine calculates the product of two fourth order tensors,
      ! A and B, and places the product in C,
      !             C_{ijkl} = A_{ijmn}B_{mnkl}.
      !
      implicit none
      !
      integer i,j,k,l,m,n
      !
      real(di) A(3,3,3,3),B(3,3,3,3),C(3,3,3,3)
      
      
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
end module stability_check
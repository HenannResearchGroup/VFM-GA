!!This module contributes to the computation of elemental external and internal residual force
!!see https://solidmechanics.org/Text/Chapter8_1/Chapter8_1.php by Allan Bower
!!reduced quadrature method as alternative approach (has pros of accuracy when using 1 integration point per element approach the accuracy of default way 4 int points per element)
!!see  https://www.sciencedirect.com/science/article/pii/S0045782522004509 , https://link.springer.com/article/10.1007/s00466-023-02280-4
!! The modules contains plane stress only, see paper for detailed explanation on the reason of choice 
module fem_compute_el
    use fem_constitutive
    use fembasics
    use constants
    use linalg
    implicit none
    contains
    !!Plane stress version of computing elemental contribution towards the external and internel residual force, can choose default way of computation or reduced quadrature
    subroutine compute_Fel_ps_after(coord,disp,theta,Finte,F33_guess,flag,F33_out,S_g_in,S_g_out,disp_old)
        implicit none
        integer nf,i,j,a,intpt,row,flag,k,l,m
        real(di) dxdxi(2,2),xi(2),dNdxi(neln,dim), &
&       dNdx(neln,dim),F(dim,dim),dNdxs(neln,dim), &
&       Finte(neln*dim),disp(dim,neln),coord(dim,neln),w(intn), &
&       theta(nparam),xilist(dim,intn),N(neln),const_stress(dim,dim), &
&       dxidx(dim,dim),ddim,F3_out,F33_guess
        real(di) xi_temp(dim),dNdxi_temp(neln,dim),dxdxi_temp(dim,dim)
        real(di) dxidx_temp(dim,dim),dNdx_temp(neln,dim),F_temp(dim,dim),F33_out
        real(di) detF_temp,F_int(dim,dim),detF,F_out(3,3),Smean(3,3)
        real(di) dudx(dim,dim),Hs(dim,dim),Ha(dim,dim),Hs_g(dim,dim,dim),Ha_g(dim,dim,dim),H(dim,dim)
        real(di) dphidcc(3,3,3,3),S_g(dim,dim,dim),dNdxidxi(dim,dim,neln)
        real(di) dudxi(dim,dim,dim),dxidxs(dim,dim),S_g_in(dim,dim,dim),S_g_out(dim,dim,dim)
        real(di) dNdy(neln,dim),dydxi(2,2),dxidy(2,2),Ja2,disp_old(dim,neln)
        real(di) S22_in,K1_in
        type(fembase)::fem
        fem = fembase(dim,neln)
        xilist = fem%intp()
        w = fem%intw()
        Finte = 0.d0
        flag = 0
        xi_temp = 0.d0
        dNdxi_temp = fem%shapedev(xi_temp)
        dxdxi_temp = 0.d0
        do i = 1,dim
            do j = 1,dim
                do a = 1,neln
                    dxdxi_temp(i,j) = dxdxi_temp(i,j)+coord(i,a)*dNdxi_temp(a,j)
                end do
            end do
        end do
        call inv(dxdxi_temp,dxidx_temp)
        dNdx_temp = 0.d0
        do a = 1, neln
            do i = 1,dim
                do j = 1,dim
                    dNdx_temp(a,i) = dNdx_temp(a,i)+dNdxi_temp(a,j)*dxidx_temp(j,i)
                end do
            end do
        end do
        call eye(F_temp,2)
        do i = 1,dim
            do j = 1,dim
                do a = 1,neln
                    F_temp(i,j) = F_temp(i,j)+disp(i,a)*dNdx_temp(a,j)
                end do
            end do
        end do        
        call det(F_temp,detF_temp)
        
        do intpt = 1,intn
            do i = 1,dim
                xi(i) = xilist(i,intpt)
            end do
            N = fem%shape(xi)
            dNdxi = fem%shapedev(xi)
            dxdxi = 0.d0
            dydxi = 0.d0
            do i = 1,dim
                do j = 1,dim
                    do a = 1,neln
                        dxdxi(i,j) = dxdxi(i,j)+coord(i,a)*dNdxi(a,j)
                        dydxi(i,j) = dydxi(i,j)+(coord(i,a)+disp(i,a))*dNdxi(a,j)
                    end do
                end do
            end do
            call inv(dxdxi,dxidx)
            call inv(dydxi,dxidy)
            call det(dydxi,Ja2)
            dNdx = 0.d0
            dNdy = 0.d0
            do a = 1, neln
                do i = 1,dim
                    do j = 1,dim
                        dNdx(a,i) = dNdx(a,i)+dNdxi(a,j)*dxidx(j,i)
                        dNdy(a,i) = dNdy(a,i)+dNdxi(a,j)*dxidy(j,i)
                    end do
                end do
            end do
            call eye(F,2)
            do i = 1,dim
                do j = 1,dim
                    do a = 1,neln
                        F(i,j) = F(i,j)+disp(i,a)*dNdx(a,j)
                    end do
                end do
            end do
            call det(F,detF)
            ddim = 2.d0
            F_int = F*(detF_temp/detF)**(1/ddim)

            call nr_solve_ps(F_int,theta,F33_guess,Smean,F_out,F3_out)        
            
            const_stress(1:2,1:2) = Smean(1:2,1:2)
            
            
            if (rflag == 1) then
                !!reduced quadrature part computation of residual force terms, updated lagrangian. 
            call compute_dphidcc(F_out,theta,dphidcc)
            
            dNdxidxi = 0.d0
            dNdxidxi(1,:,1)=(/0.d0,1.d0/4.d0/)
            dNdxidxi(2,:,1)= (/1.d0/4.d0,0.d0/)
            dNdxidxi(1,:,2)=(/0.d0,-1.d0/4.d0/)
            dNdxidxi(2,:,2)=(/-1.d0/4.d0,0.d0/)
            dNdxidxi(1,:,3)=(/0.d0,1.d0/4.d0/)
            dNdxidxi(2,:,3)=(/1.d0/4.d0,0.d0/)
            dNdxidxi(1,:,4)=(/0.d0,-1.d0/4.d0/)
            dNdxidxi(2,:,4)=(/-1.d0/4.d0,0.d0/)

            do i = 1,dim
                do j = 1,dim
                    H(i,j) = 0.d0
                if (i==j) then
                    H(i,j) = 1.d0
                end if
                
                do a = 1,neln
                H(i,j) = H(i,j) + (disp(i,a)-disp_old(i,a))*dNdy(a,j)
                    end do
                 end do
            end do
            Hs = 0.5d0*(H+transpose(H))
            Ha = 0.5d0*(H-transpose(H))

            do i = 1,dim
                do j = 1,dim
                    do  k=1,dim
                    dudxi(i,j,k) = 0.d0
                        do a = 1,neln
                            do l=1,dim
                            dudxi(i,j,k) = dudxi(i,j,k)+(disp(i,a)-disp_old(i,a))*dNdxidxi(k,l,a)*dxidy(l,j)
                            end do
                        end do
                    end do
                end do
            end do
            do a=1,dim
                Ha_g(:,:,a)=0.5d0*(dudxi(:,:,a)-transpose(dudxi(:,:,a)));
                Hs_g(:,:,a)=0.5d0*(dudxi(:,:,a)+transpose(dudxi(:,:,a)));
            end do
            
            S_g = S_g_in
            do i=1,dim
                do j=1,dim
                        do k=1,dim
                            do m = 1,dim
                                S_g(i,j,m) = S_g(i,j,m) + (Hs_g(i,k,m)+Ha_g(i,k,m))*const_stress(k,j) &
                                & + (Hs_g(j,k,m)+Ha_g(j,k,m))*const_stress(k,i) &
                                & - Hs_g(k,k,m)*const_stress(i,j) &
                                & + (Hs(i,k)+Ha(i,k))*S_g(k,j,m) &
                                & + (Hs(j,k)+Ha(j,k))*S_g(k,i,m) &
                                & - Hs(k,k)*S_g(i,j,m)
                            end do
                        end do
                end do
            end do
                                
            do i=1,dim
                do j=1,dim
                    do m=1,dim
                        do l=1,dim
                            do k = 1,dim
                                S_g(i,j,m)=S_g(i,j,m)+dphidcc(i,j,k,l)*Hs_g(k,l,m)
                            end do
                        end do
                    end do 
                end do 
            end do
            S_g_out = S_g
            
            do a = 1,neln
                do i = 1,dim
                row = dim*(a-1)+i
                    do j = 1,dim
                        do m = 1,dim
                            do k=1,dim
                 Finte(row) = Finte(row) + S_g(i,j,m)*dNdxidxi(k,m,a)*dxidy(k,j)*w(intpt)/3.d0*Ja2
                            end do
                        end do
                    end do
                end do
            end do
            end if
            do a = 1,neln
                do i = 1,dim
                    row = dim*(a-1)+i
                    do j = 1,dim
                        Finte(row) = Finte(row)+const_stress(i,j)*dNdy(a,j)*w(intpt)*Ja2
                    end do
                end do
            end do
                        
        end do
        
        return
    end subroutine compute_Fel_ps_after
    !!Plane stress version of computing elemental contribution towards the external and internel residual force, checks can be turned on and off based on whether the input field is homogeneous 
    subroutine compute_Fel_ps_monoK1(coord,disp,theta,Finte,flag,F33_out,S22_out,K1_out,F11_out,S_g_in,S_g_out,disp_old)
        implicit none
        integer nf,i,j,a,intpt,row,flag,k,l,m
        real(di) dxdxi(2,2),xi(2),dNdxi(neln,dim), &
&       dNdx(neln,dim),F(dim,dim),dNdxs(neln,dim), &
&       Finte(neln*dim),disp(dim,neln),coord(dim,neln),w(intn), &
&       theta(nparam),xilist(dim,intn),N(neln),const_stress(dim,dim), &
&       dxidx(dim,dim),ddim,F3_out,S22_in,S22_out,K1_in,K1_out,disp_old(dim,neln)
        real(di) xi_temp(dim),dNdxi_temp(neln,dim),dxdxi_temp(dim,dim),jd,F11_out
        real(di) dxidx_temp(dim,dim),dNdx_temp(neln,dim),F_temp(dim,dim),F33_out
        real(di) detF_temp,F_int(dim,dim),detF,F_out(3,3),Smean(3,3)
        real(di) dudx(dim,dim),Hs(dim,dim),Ha(dim,dim),Hs_g(dim,dim,dim),Ha_g(dim,dim,dim),H(dim,dim)
        real(di) dphidcc(3,3,3,3),S_g(dim,dim,dim),dNdxidxi(dim,dim,neln)
        real(di) dudxi(dim,dim,dim),dxidxs(dim,dim),S_g_in(dim,dim,dim),S_g_out(dim,dim,dim)
        real(di) dNdy(neln,dim),dydxi(2,2),dxidy(2,2),Ja2
        type(fembase)::fem
        fem = fembase(dim,neln)
        xilist = fem%intp()
        w = fem%intw()
        Finte = 0.d0
        flag = 0
        xi_temp = 0.d0
        dNdxi_temp = fem%shapedev(xi_temp)
        dxdxi_temp = 0.d0
        F33_out = 1.d0
        do i = 1,dim
            do j = 1,dim
                do a = 1,neln
                    dxdxi_temp(i,j) = dxdxi_temp(i,j)+coord(i,a)*dNdxi_temp(a,j)
                end do
            end do
        end do
        call inv(dxdxi_temp,dxidx_temp)
        dNdx_temp = 0.d0
        do a = 1, neln
            do i = 1,dim
                do j = 1,dim
                    dNdx_temp(a,i) = dNdx_temp(a,i)+dNdxi_temp(a,j)*dxidx_temp(j,i)
                end do
            end do
        end do
        call eye(F_temp,2)
        do i = 1,dim
            do j = 1,dim
                do a = 1,neln
                    F_temp(i,j) = F_temp(i,j)+disp(i,a)*dNdx_temp(a,j)
                end do
            end do
        end do        
        call det(F_temp,detF_temp)
        
        do intpt = 1,intn
            do i = 1,dim
                xi(i) = xilist(i,intpt)
            end do
            N = fem%shape(xi)
            dNdxi = fem%shapedev(xi)
            dxdxi = 0.d0
            dydxi = 0.d0
            do i = 1,dim
                do j = 1,dim
                    do a = 1,neln
                        dxdxi(i,j) = dxdxi(i,j)+coord(i,a)*dNdxi(a,j)
                        dydxi(i,j) = dydxi(i,j)+(coord(i,a)+disp(i,a))*dNdxi(a,j)
                    end do
                end do
            end do
            call inv(dydxi,dxidy)
            call det(dydxi,Ja2)
            call inv(dxdxi,dxidx)
            dNdx = 0.d0
            dNdy = 0.d0
            do a = 1, neln
                do i = 1,dim
                    do j = 1,dim
                        dNdx(a,i) = dNdx(a,i)+dNdxi(a,j)*dxidx(j,i)
                        dNdy(a,i) = dNdy(a,i)+dNdxi(a,j)*dxidy(j,i)
                    end do
                end do
            end do
            call eye(F,2)
            do i = 1,dim
                do j = 1,dim
                    do a = 1,neln
                        F(i,j) = F(i,j)+disp(i,a)*dNdx(a,j)
                    end do
                end do
            end do
            call det(F,detF)
            ddim = 2.d0
            F_int = F*(detF_temp/detF)**(1/ddim)
            call solve_plane_stress_monoK1(F_int,theta,F_out,Smean,flag,F3_out,S22_in,K1_in)
            F33_out = F3_out
            F11_out = F_int(1,1)
            call det3d(F_out,jd)
            const_stress(1:2,1:2) = Smean(1:2,1:2)
            
            if (rflag == 1) then
            call compute_dphidcc(F_out,theta,dphidcc)
            dNdxidxi = 0.d0
            dNdxidxi(1,:,1)=(/0.d0,1.d0/4.d0/)
            dNdxidxi(2,:,1)= (/1.d0/4.d0,0.d0/)
            dNdxidxi(1,:,2)=(/0.d0,-1.d0/4.d0/)
            dNdxidxi(2,:,2)=(/-1.d0/4.d0,0.d0/)
            dNdxidxi(1,:,3)=(/0.d0,1.d0/4.d0/)
            dNdxidxi(2,:,3)=(/1.d0/4.d0,0.d0/)
            dNdxidxi(1,:,4)=(/0.d0,-1.d0/4.d0/)
            dNdxidxi(2,:,4)=(/-1.d0/4.d0,0.d0/)
            
            do i = 1,dim
                do j = 1,dim
                    H(i,j) = 0.d0
                if (i==j) then
                    H(i,j) = 1.d0
                end if
                
                do a = 1,neln
                H(i,j) = H(i,j) + (disp(i,a)-disp_old(i,a))*dNdy(a,j)
                    end do
                 end do
            end do
            Hs = 0.5d0*(H+transpose(H))
            Ha = 0.5d0*(H-transpose(H))

            do i = 1,dim
                do j = 1,dim
                    do  k=1,dim
                    dudxi(i,j,k) = 0.d0
                        do a = 1,neln
                            do l=1,dim
                            dudxi(i,j,k) = dudxi(i,j,k)+(disp(i,a)-disp_old(i,a))*dNdxidxi(k,l,a)*dxidy(l,j)
                            end do
                        end do
                    end do
                end do
            end do
            do a=1,dim
                Ha_g(:,:,a)=0.5d0*(dudxi(:,:,a)-transpose(dudxi(:,:,a)));
                Hs_g(:,:,a)=0.5d0*(dudxi(:,:,a)+transpose(dudxi(:,:,a)));
            end do
            S_g = S_g_in
            do i=1,dim
                do j=1,dim
                        do k=1,dim
                            do m = 1,dim
                                S_g(i,j,m) = S_g(i,j,m) + (Hs_g(i,k,m)+Ha_g(i,k,m))*const_stress(k,j) &
                                & + (Hs_g(j,k,m)+Ha_g(j,k,m))*const_stress(k,i) &
                                & - Hs_g(k,k,m)*const_stress(i,j) &
                                & + (Hs(i,k)+Ha(i,k))*S_g(k,j,m) &
                                & + (Hs(j,k)+Ha(j,k))*S_g(k,i,m) &
                                & - Hs(k,k)*S_g(i,j,m)
                            end do
                        end do
                end do
            end do
            do i=1,dim
                do j=1,dim
                    do m=1,dim
                        do l=1,dim
                            do k = 1,dim
                                S_g(i,j,m)=S_g(i,j,m)+dphidcc(i,j,k,l)*Hs_g(k,l,m)
                            end do
                        end do
                    end do 
                end do 
            end do
            S_g_out = S_g

            do a = 1,neln
                do i = 1,dim
                row = dim*(a-1)+i
                    do j = 1,dim
                        do m = 1,dim
                            do k=1,dim
                   Finte(row) = Finte(row) + S_g(i,j,m)*dNdxidxi(k,m,a)*dxidy(k,j)*w(intpt)/3.d0*Ja2
                            end do
                        end do
                    end do
                end do
            end do
            endif
            
            do a = 1,neln
                do i = 1,dim
                    row = dim*(a-1)+i
                    do j = 1,dim
                        Finte(row) = Finte(row)+const_stress(i,j)*dNdy(a,j)*w(intpt)*Ja2
                    end do
                end do
            end do
        end do
        
        S22_out = S22_in/jd
        K1_out = K1_in
        return
    end subroutine compute_Fel_ps_monoK1


end module fem_compute_el
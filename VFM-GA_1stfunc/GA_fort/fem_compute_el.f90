!!This module contributes to the computation of elemental external and internal residual force
!!see https://solidmechanics.org/Text/Chapter8_1/Chapter8_1.php by Allan Bower
module fem_compute_el
    use fem_constitutive
    use fembasics
    use constants
    use linalg
    implicit none
    contains
    !!For fun only, was used to prove that for a model with a lot of material parameters, basic gradient based optimization method is extremely hard to work
    subroutine compute_Gel(coord,disp,theta,Fintegrad)
        implicit none
        integer nf,i,j,a,intpt,row,k
        real(di) dxdxi(2,2),Ja,xi(2),dNdxi(neln,dim), &
&       dNdx(neln,dim),F(dim,dim),Finv(dim,dim),dNdxs(neln,dim), &
&       Fintegrad(neln*dim,nparam),disp(dim,neln),coord(dim,neln),w(intn), &
&       theta(nparam),xilist(dim,intn),N(neln),const_grad(dim,dim*nparam), &
&       dxidx(dim,dim),grad_temp(dim,dim)
        real(di) xi_temp(dim),dNdxi_temp(neln,dim),dxdxi_temp(dim,dim)
        real(di) dxidx_temp(dim,dim),dNdx_temp(neln,dim),F_temp(dim,dim)
        real(di) detF_temp,F_int(dim,dim),detF
        type(fembase)::fem
        fem = fembase(dim,neln)
        xilist = fem%intp()
        w = fem%intw()
        
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
            do i = 1,dim
                do j = 1,dim
                    do a = 1,neln
                        dxdxi(i,j) = dxdxi(i,j)+coord(i,a)*dNdxi(a,j)
                    end do
                end do
            end do
            call inv(dxdxi,dxidx)
            call det(dxdxi,Ja)
            dNdx = 0.d0
            do a = 1, neln
                do i = 1,dim
                    do j = 1,dim
                        dNdx(a,i) = dNdx(a,i)+dNdxi(a,j)*dxidx(j,i)
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
            call inv(F,Finv)
            dNdxs = 0.d0
            do a = 1,neln
                do i = 1,dim
                    do j = 1,dim
                        dNdxs(a,i) = dNdxs(a,i)+dNdx(a,j)*Finv(j,i)
                    end do
                end do
            end do
            call det(F,detF)
            F_int = F*(detF_temp/detF)**(1/dim)
            call const_grad_foam(F_int,theta,const_grad)
            Fintegrad = 0.d0
            do k = 1,nparam
                grad_temp(1:2,1:2) = const_grad(1:2,(nparam*2-1):(nparam*2))
                do a = 1,neln
                    do i = 1,dim
                        row = dim*(a-1)+i
                        do j = 1,dim
                            Fintegrad(row,k) = Fintegrad(row,k)+grad_temp(i,j)*dNdxs(a,j)*w(intpt)*Ja
                        end do
                    end do
                end do  
            end do
        end do
        return
    end subroutine compute_Gel
    !!Plane strain version of computing elemental contribution towards the external and internel residual force 
    subroutine compute_Fel(coord,disp,theta,Finte)
        implicit none
        integer nf,i,j,a,intpt,row
        real(di) dxdxi(2,2),Ja,xi(2),dNdxi(neln,dim), &
&       dNdx(neln,dim),F(dim,dim),Finv(dim,dim),dNdxs(neln,dim), &
&       Finte(neln*dim),disp(dim,neln),coord(dim,neln),w(intn), &
&       theta(nparam),xilist(dim,intn),N(neln),const_stress(dim,dim), &
&       dxidx(dim,dim),ddim
        real(di) xi_temp(dim),dNdxi_temp(neln,dim),dxdxi_temp(dim,dim)
        real(di) dxidx_temp(dim,dim),dNdx_temp(neln,dim),F_temp(dim,dim)
        real(di) detF_temp,F_int(dim,dim),detF
        type(fembase)::fem
        fem = fembase(dim,neln)
        xilist = fem%intp()
        w = fem%intw()
        Finte = 0.d0
        
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
            do i = 1,dim
                do j = 1,dim
                    do a = 1,neln
                        dxdxi(i,j) = dxdxi(i,j)+coord(i,a)*dNdxi(a,j)
                    end do
                end do
            end do
            call inv(dxdxi,dxidx)
            call det(dxdxi,Ja)
            dNdx = 0.d0
            do a = 1, neln
                do i = 1,dim
                    do j = 1,dim
                        dNdx(a,i) = dNdx(a,i)+dNdxi(a,j)*dxidx(j,i)
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
            call inv(F,Finv)
            dNdxs = 0.d0
            do a = 1,neln
                do i = 1,dim
                    do j = 1,dim
                        dNdxs(a,i) = dNdxs(a,i)+dNdx(a,j)*Finv(j,i)
                    end do
                end do
            end do
            call det(F,detF)
            ddim = 2.d0
            F_int = F*(detF_temp/detF)**(1/ddim)
            call const_stress_foam(F_int,theta,const_stress)
            
            do a = 1,neln
                do i = 1,dim
                    row = dim*(a-1)+i
                    do j = 1,dim
                        Finte(row) = Finte(row)+const_stress(i,j)*dNdxs(a,j)*w(intpt)*Ja
                    end do
                end do
            end do
        end do
        return
    end subroutine compute_Fel
    !!Plane stress version of computing elemental contribution towards the external and internel residual force 
    subroutine compute_Fel_ps_after(coord,disp,theta,Finte,F33_guess,flag,F33_out)
        implicit none
        integer nf,i,j,a,intpt,row,flag
        real(di) dxdxi(2,2),Ja,xi(2),dNdxi(neln,dim), &
&       dNdx(neln,dim),F(dim,dim),Finv(dim,dim),dNdxs(neln,dim), &
&       Finte(neln*dim),disp(dim,neln),coord(dim,neln),w(intn), &
&       theta(nparam),xilist(dim,intn),N(neln),const_stress(dim,dim), &
&       dxidx(dim,dim),ddim,F3_out,F33_guess
        real(di) xi_temp(dim),dNdxi_temp(neln,dim),dxdxi_temp(dim,dim)
        real(di) dxidx_temp(dim,dim),dNdx_temp(neln,dim),F_temp(dim,dim),F33_out
        real(di) detF_temp,F_int(dim,dim),detF,F_out(3,3),Smean(3,3),Finv_temp(3,3)
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
            do i = 1,dim
                do j = 1,dim
                    do a = 1,neln
                        dxdxi(i,j) = dxdxi(i,j)+coord(i,a)*dNdxi(a,j)
                    end do
                end do
            end do
            call inv(dxdxi,dxidx)
            call det(dxdxi,Ja)
            dNdx = 0.d0
            do a = 1, neln
                do i = 1,dim
                    do j = 1,dim
                        dNdx(a,i) = dNdx(a,i)+dNdxi(a,j)*dxidx(j,i)
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
            F33_out = F3_out
            call inv3d(F_out,Finv_temp)
            Finv(1:2,1:2) = Finv_temp(1:2,1:2)
            const_stress(1:2,1:2) = Smean(1:2,1:2)
            dNdxs = 0.d0
            do a = 1,neln
                do i = 1,dim
                    do j = 1,dim
                        dNdxs(a,i) = dNdxs(a,i)+dNdx(a,j)*Finv(j,i)
                    end do
                end do
            end do
        
            do a = 1,neln
                do i = 1,dim
                    row = dim*(a-1)+i
                    do j = 1,dim
                        Finte(row) = Finte(row)+const_stress(i,j)*dNdxs(a,j)*w(intpt)*Ja
                    end do
                end do
            end do
        end do
        return
    end subroutine compute_Fel_ps_after
    !!Plane stress version of computing elemental contribution towards the external and internel residual force, (alter)
    subroutine compute_Fel_ps_monok1(coord,disp,theta,Finte,flag,F33_out,S22_out,K1_out,F11_out)
        implicit none
        integer nf,i,j,a,intpt,row,flag
        real(di) dxdxi(2,2),Ja,xi(2),dNdxi(neln,dim), &
&       dNdx(neln,dim),F(dim,dim),Finv(dim,dim),dNdxs(neln,dim), &
&       Finte(neln*dim),disp(dim,neln),coord(dim,neln),w(intn), &
&       theta(nparam),xilist(dim,intn),N(neln),const_stress(dim,dim), &
&       dxidx(dim,dim),ddim,F3_out,S22_in,S22_out,K1_in,K1_out
        real(di) xi_temp(dim),dNdxi_temp(neln,dim),dxdxi_temp(dim,dim),jd,F11_out
        real(di) dxidx_temp(dim,dim),dNdx_temp(neln,dim),F_temp(dim,dim),F33_out
        real(di) detF_temp,F_int(dim,dim),detF,F_out(3,3),Smean(3,3),Finv_temp(3,3)
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
            do i = 1,dim
                do j = 1,dim
                    do a = 1,neln
                        dxdxi(i,j) = dxdxi(i,j)+coord(i,a)*dNdxi(a,j)
                    end do
                end do
            end do
            call inv(dxdxi,dxidx)
            call det(dxdxi,Ja)
            dNdx = 0.d0
            do a = 1, neln
                do i = 1,dim
                    do j = 1,dim
                        dNdx(a,i) = dNdx(a,i)+dNdxi(a,j)*dxidx(j,i)
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
            call solve_plane_stress_monok1(F_int,theta,F_out,Smean,flag,F3_out,S22_in,K1_in)
            F33_out = F3_out
            F11_out = F_int(1,1)
            call inv3d(F_out,Finv_temp)
            call det3d(F_out,jd)
            Finv(1:2,1:2) = Finv_temp(1:2,1:2)
            const_stress(1:2,1:2) = Smean(1:2,1:2)
            dNdxs = 0.d0
            do a = 1,neln
                do i = 1,dim
                    do j = 1,dim
                        dNdxs(a,i) = dNdxs(a,i)+dNdx(a,j)*Finv(j,i)
                    end do
                end do
            end do

            
            do a = 1,neln
                do i = 1,dim
                    row = dim*(a-1)+i
                    do j = 1,dim
                        Finte(row) = Finte(row)+const_stress(i,j)*dNdxs(a,j)*w(intpt)*Ja
                    end do
                end do
            end do
        end do
        S22_out = S22_in/jd
        K1_out = K1_in
        return
    end subroutine compute_Fel_ps_monok1

end module fem_compute_el
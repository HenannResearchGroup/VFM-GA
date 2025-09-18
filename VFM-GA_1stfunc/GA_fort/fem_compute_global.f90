!! This module corresponds to the subroutines that compute the objective function by summing up the internal 
!! and external residual force and summing up the external force measured in experiment    
!! The main part of the code is global assembly of elemental residual force and regrouping them into objective function 
module fem_compute_gl
    use constants
    use fem_compute_el
    use stability_check
    implicit none
    contains
    !Obtain the force external measured in experiment from the input data
    !in 1st funtionality, type is fixed, no need to worry about it 
    subroutine get_Fextexp(rf,Fextexp,free_count)
        implicit none
        integer i,free_count,nn_m_sqrt
        real(di) Fextexp(group),rf(nnode,dim)
        nn_m_sqrt = nnode-sqnn
        Fextexp = 0.d0
        free_count = 0
        if (bctype == 1) then
            do i = 1,nnode
                if (rf(i,1)>0.d0) then
                    Fextexp(1) = Fextexp(1)+rf(i,1)
                else if (rf(i,1)<0.d0) then
                    Fextexp(3) = Fextexp(3)+rf(i,1)
                else
                    free_count = free_count+1
                end if
                if (rf(i,2)>0.d0) then
                    Fextexp(2) = Fextexp(2)+rf(i,2)
                else if (rf(i,2)<0.d0) then
                    Fextexp(4) = Fextexp(4)+rf(i,2)
                else
                    free_count = free_count+1
                end if
            end do
        else if (bctype == 2) then
            do i = 1,nnode
                if (rf(i,2)>0.d0) then
                    Fextexp(1) = Fextexp(1)+rf(i,1)
                    Fextexp(2) = Fextexp(2)+rf(i,2)
                else if (rf(i,2)<0.d0) then
                    Fextexp(3) = Fextexp(3)+rf(i,1)
                    Fextexp(4) = Fextexp(4)+rf(i,2)
                else
                    free_count=free_count+2
                end if
            end do
        else if (bctype == 3) then
            do i = 1,nnode
                if (rf(i,2)<0.d0) then
                    Fextexp(1) = Fextexp(1)+rf(i,1)
                    Fextexp(2) = Fextexp(2)+rf(i,2)
                else
                    if ((rf(i,2)>0.d0) .and. (rf(i,1)==0.d0)) then
                         Fextexp(3) = Fextexp(3) +rf(i,2)
                         free_count=free_count+1
                    else if ((rf(i,2)>0.d0) .and. (rf(i,1).NE.0.d0)) then
                         Fextexp(3) = Fextexp(3) +rf(i,2)
                         Fextexp(4) = Fextexp(4) +rf(i,1)
                    else if  ((rf(i,1).NE.0.d0).and. (rf(i,2)==0.d0)) then
                         Fextexp(4) = Fextexp(4) +rf(i,1)
                         free_count=free_count+1
                    else 
                         free_count=free_count+2
                    end if 
                end if
            end do
        else if (bctype == 4) then
            do i = 1,nnode
                if (i > (nn_m_sqrt)) then
                    Fextexp(2) = Fextexp(2)+rf(i,2)
                end if
                if (rf(i,1).NE.0.d0) then
                    Fextexp(1) = Fextexp(1)+rf(i,1)
                end if 
            end do
            Fextexp(3) = -Fextexp(1)
            Fextexp(4) = -Fextexp(2)
            free_count = nnode*2-sqnn*2-sqnn-(sqnn-1)
        end if
        return
    end subroutine get_Fextexp
    !! plane strain version of objective function computation
    subroutine compute_MSE(coord,disp,connect,rf,obj,theta) 
        use constants
        use fem_compute_el
        implicit none
        real(di) obj,coord(dim,nnode),disp(nframe,nnode*dim), &
&       rf_temp(nnode,dim),Fextexp(group),disp_temp(nnode*dim), &
&       Ftotal(dim*nnode),coord_el(dim,intn),disp_el(dim,intn), &
&       theta(nparam),F_el(neln*dim),Fext(group), &
&       norm_list(nframe),Ext_diff(group),Fextexp_norm, &
&       Ext_diff_norm,Finte_norm,obj_list(nframe),rf(nnode,nframe*dim)
        real(di),allocatable::Finte(:)
        integer connect(neln,nel),free_dof,rw
        integer f,el,a,i,free_index,nn_m_sqrt,sqrt_nn
        
        
        obj = 0.d0
        obj_list = 0.d0
        do f = 1,nframe
            Ftotal = 0.d0
            Fext = 0.d0
            Fextexp = 0.d0
            rf_temp(1:nnode,1:2) = rf(1:nnode,f*2-1:f*2)
            call get_Fextexp(rf_temp,Fextexp,free_dof)
            disp_temp(:) = disp(f,:)
            do el=1,nel
                do a = 1,neln
                    do i = 1,dim
                        disp_el(i,a) = disp_temp(dim*(connect(a,el)-1)+i)
                        coord_el(i,a) = coord(i,connect(a,el))
                    end do
                end do
                F_el = 0.d0
                call compute_Fel(coord_el,disp_el,theta,F_el)
                do a=1,neln
                    do i = 1,dim
                        rw = dim*(connect(a,el)-1)+i
                        Ftotal(rw) = Ftotal(rw) + F_el(dim*(a-1)+i)
                    end do
                end do
            end do
            
            allocate(Finte(free_dof))
            
            Finte = 0.d0
            
            free_index =1
            if (bctype == 1) then
                do i = 1,nnode
                    if (rf_temp(i,1)>0.d0) then
                        Fext(1) = Fext(1) +Ftotal(i*2-1)
                    else if (rf_temp(i,1)<0.d0) then
                        Fext(3) = Fext(3) +Ftotal(i*2-1)
                    else
                        Finte(free_index) = Ftotal(i*2-1)
                        free_index = free_index+1
                    end if 
                    if (rf_temp(i,2) > 0.d0) then
                        Fext(2) = Fext(2) +Ftotal(i*2)
                    else if (rf_temp(i,2)<0.d0) then
                        Fext(4) = Fext(4) +Ftotal(i*2)
                    else
                        Finte(free_index) = Ftotal(i*2)
                        free_index = free_index+1 
                    end if
                end do
            else if (bctype == 2) then
                do i = 1,nnode
                    if (rf_temp(i,2)>0.d0) then
                        Fext(1) = Fext(1)+Ftotal(i*2-1)
                        Fext(2) = Fext(2)+Ftotal(i*2)
                    else if (rf_temp(i,2)<0.d0) then
                        Fext(3) = Fext(3)+Ftotal(i*2-1)
                        Fext(4) = Fext(4)+Ftotal(i*2)
                    else
                        Finte(free_index) = Ftotal(i*2-1)
                        Finte(free_index+1) = Ftotal(i*2)
                        free_index = free_index+2
                    end if
                end do
            else if (bctype == 3) then
                ! 1:dof 1 bottom, 2:dof 2 bottom 3.dof 2 top 4. dof 1 side
                do i = 1,nnode
                    !print*,i
                    if (rf_temp(i,2)<0.d0) then
                        !print*,'1 & 2'
                        Fext(1) = Fext(1)+Ftotal(i*2-1)
                        Fext(2) = Fext(2)+Ftotal(i*2)
                    else
                        if ((rf_temp(i,2)>0.d0) .and. (rf_temp(i,1)==0.d0)) then
                            !print*,'3'
                            Fext(3) = Fext(3) +Ftotal(i*2)
                            Finte(free_index) = Ftotal(i*2-1)
                            free_index = free_index+1
                        else if ((rf_temp(i,2)>0.d0) .and. (rf_temp(i,1).NE.0.d0)) then
                            !print*,'3 & 4'
                            Fext(3) = Fext(3) +Ftotal(i*2)
                            Fext(4) = Fext(4) +Ftotal(i*2-1)
                        else if  ((rf_temp(i,1).NE.0.d0).and. (rf_temp(i,2)==0.d0)) then
                            !print*,' 4'
                            Fext(4) = Fext(4) +Ftotal(i*2-1)
                            Finte(free_index) = Ftotal(i*2)
                            free_index = free_index+1
                        else 
                            !print*,'free'
                             Finte(free_index) = Ftotal(i*2-1)
                             Finte(free_index+1) = Ftotal(i*2)
                             free_index= free_index+2
                        end if 
                    end if
                end do
            else if (bctype == 4) then
                
                nn_m_sqrt = nnode-sqnn
                
                do i = 1,nnode
                    if (i>(nn_m_sqrt)) then
                        Fext(2)  = Fext(2) + Ftotal(i*2)
                        if (rf(i,1).eq.0.d0) then
                            Finte(free_index) = Ftotal(i*2-1)
                            free_index = free_index+1
                        else
                            Fext(1) = Fext(1)+Ftotal(i*2-1)
                        end if 
                    else if (i<sqnn) then
                        Fext(3) = Fext(3)+Ftotal(i*2-1)
                        Fext(4) = Fext(4)+Ftotal(i*2)
                        if (rf(i,1).ne.0.d0) then
                            Fext(1) = Fext(1)+Ftotal(i*2-1)
                        end if
                    else
                        if (rf(i,1).ne.0.d0) then
                            Fext(1) = Fext(1)+Ftotal(i*2-1)
                            if (i>sqnn) then
                                Finte(free_index) = Ftotal(i*2)
                                free_index = free_index+1
                            else 
                                Fext(4) = Fext(4) + Ftotal(i*2)
                            end if
                        else
                            Finte(free_index) = Ftotal(i*2-1)
                            Finte(free_index+1) = Ftotal(i*2)
                            free_index = free_index+2
                        end if                        
                    end if 
                    
                end do
            end if 
            
            Ext_diff = Fext-Fextexp
            Ext_diff_norm = norm2(Ext_diff)
            Finte_norm = norm2(Finte)*alpha
            Fextexp_norm = norm2(Fextexp)
            obj_list(f) = 0.d0
            obj_list(f) = Finte_norm+Ext_diff_norm
            norm_list(f) = Fextexp_norm
            obj = obj+obj_list(f)/norm_list(f)
            deallocate(Finte)    
        end do
        return
    end subroutine compute_MSE
    !!!!!
    !! plane stress version of objective function computation
    subroutine compute_MSE_ps(coord,disp,connect,rf,obj,theta,c_mono,c_k1) 
        use constants
        use fem_compute_el
        use stability_check
        implicit none
        real(di) obj,coord(dim,nnode),disp(nframe,nnode*dim), &
&       rf_temp(nnode,dim),Fextexp(group),disp_temp(nnode*dim), &
&       Ftotal(dim*nnode),coord_el(dim,intn),disp_el(dim,intn), &
&       theta(nparam),F_el(neln*dim),Fext(group), &
&       norm_list(nframe),Ext_diff(group),Fextexp_norm,obj_hardcode, &
&       Ext_diff_norm,Finte_norm,obj_list(nframe),rf(nnode,nframe*dim),F33_guess,F33_out,S22_out
        real(di) S22_now,S22_pre,S22_pre2,S22_diff,K1,F11_out,K1_old
        real(di),allocatable::Finte(:)
        integer connect(neln,nel),free_dof,rw,c_mono,c_k1
        integer f,el,a,i,free_index,flag,flag2,flag3,stability,flagk1,flagmono
        !c_mono for mono version, c_k1 for k1 version
        !flag 2 is output flag, flag 3 penalty flag, flag is solver failure flag
        obj = 0.d0
        obj_hardcode = 100000000000.d0
        obj_list = 0.d0
        S22_pre = 0.d0
        S22_pre2 = 0.d0
        flagk1 = 0
        flagmono = 0
        flag2 = 0
        flag = 0
        stability = 0
        call sta_pack(theta,stability)
        if (stability == 1) then
            obj = obj_hardcode
        else
            do f = 1,nframe
                flag3 = 0
                F33_guess = 1.d0
                Ftotal = 0.d0
                Fext = 0.d0
                Fextexp = 0.d0
                rf_temp(1:nnode,1:2) = rf(1:nnode,f*2-1:f*2)
                call get_Fextexp(rf_temp,Fextexp,free_dof)
                disp_temp(:) = disp(f,:)
            
                do el=1,nel
                    do a = 1,neln
                        do i = 1,dim
                            disp_el(i,a) = disp_temp(dim*(connect(a,el)-1)+i)
                            coord_el(i,a) = coord(i,connect(a,el))
                        end do
                    end do
                    F_el = 0.d0
                
                   if (el == 1) then
                        call compute_Fel_ps_monok1(coord_el,disp_el,theta,F_el,flag,F33_out,S22_out,K1,F11_out)
                        if (abs(F11_out-F33_out)>0.032d0) then
                            !penalty for inconsistent F11 F33 in uniaxial experiment
                            flag3 = 1
                        endif
                        if (flag == 1) then
                            flag2 = 1
                        endif
                    
                        S22_now = S22_out
                        F33_guess = F33_out
                        if (f <break) then
                            if (S22_now>S22_pre) then
                                flagmono = 1
                            else
                                S22_pre = S22_now
                            endif 
                            if (F33_out<1.d0) then
                                flagk1 = 1
                            endif
                        
                        else 
			    if (f == break) then
				K1_old = 0.d0
			    endif
			    
                            if (K1_old.ge.K1) then
                                flagk1 = 1
			    else
				K1_old = K1
                            endif
                            S22_diff = S22_now-S22_pre2
                            if (S22_diff>500.d0) then
                                S22_pre2 = S22_now
                            else
                                flagmono = 1
                            endif
                        end if
                    
                        if (c_mono == 1) then
                            if (flagmono == 1) then
                                flag2 = 1
                            endif
                        else
                        endif
                    
                        if (c_k1 == 1) then
                            if (flagk1 == 1) then
                                flag2 = 1
                            endif
                        else
                        endif
                   else
                        call compute_Fel_ps_after(coord_el,disp_el,theta,F_el,F33_guess,flag,F33_out)
                        if (flag == 1) then
                            flag2 = 1
                        endif
                        F33_guess = F33_out
                   end if
                    do a=1,neln
                        do i = 1,dim
                            rw = dim*(connect(a,el)-1)+i
                            Ftotal(rw) = Ftotal(rw) + F_el(dim*(a-1)+i)
                        end do
                    end do
                end do
                if (flag2 == 1) then
                    obj = obj_hardcode
                    exit
                end if
                allocate(Finte(free_dof))
            
                Finte = 0.d0
                Fext = 0.d0
                free_index =1
                if (bctype == 1) then
                    do i = 1,nnode
                        if (rf_temp(i,1)>0.d0) then
                           Fext(1) = Fext(1) +Ftotal(i*2-1)
                        else if (rf_temp(i,1)<0.d0) then
                            Fext(3) = Fext(3) +Ftotal(i*2-1)
                        else
                            Finte(free_index) = Ftotal(i*2-1)
                            free_index = free_index+1
                        end if 
                        if (rf_temp(i,2) > 0.d0) then
                            Fext(2) = Fext(2) +Ftotal(i*2)
                        else if (rf_temp(i,2)<0.d0) then
                            Fext(4) = Fext(4) +Ftotal(i*2) 
                        else
                            Finte(free_index) = Ftotal(i*2)
                            free_index = free_index+1 
                        end if
                    end do
                else if (bctype == 2) then
                    do i = 1,nnode 
                        if (rf_temp(i,2)>0.d0) then
                            Fext(1) = Fext(1)+Ftotal(i*2-1)
                            Fext(2) = Fext(2)+Ftotal(i*2)
                        else if (rf_temp(i,2)<0.d0) then
                            Fext(3) = Fext(3)+Ftotal(i*2-1)
                            Fext(4) = Fext(4)+Ftotal(i*2)
                        else 
                            Finte(free_index) = Ftotal(i*2-1)
                            Finte(free_index+1) = Ftotal(i*2)
                            free_index = free_index+2
                        end if
                    end do
                    !print*,Finte
                else if (bctype == 3) then
                    ! 1:dof 1 bottom, 2:dof 2 bottom 3.dof 2 top 4. dof 1 side
                    do i = 1,nnode
                        !print*,i
                        if (rf_temp(i,2)<0.d0) then
                            !print*,'1 & 2'
                            Fext(1) = Fext(1)+Ftotal(i*2-1)
                            Fext(2) = Fext(2)+Ftotal(i*2)
                        else
                            if ((rf_temp(i,2)>0.d0) .and. (rf_temp(i,1)==0.d0)) then
                                !print*,'3'
                                Fext(3) = Fext(3) +Ftotal(i*2)
                                Finte(free_index) = Ftotal(i*2-1)
                                free_index = free_index+1
                            else if ((rf_temp(i,2)>0.d0) .and. (rf_temp(i,1).NE.0.d0)) then
                                !print*,'3 & 4'
                                Fext(3) = Fext(3) +Ftotal(i*2)
                                Fext(4) = Fext(4) +Ftotal(i*2-1)
                            else if  ((rf_temp(i,1).NE.0.d0).and. (rf_temp(i,2)==0.d0)) then
                                !print*,' 4'
                                Fext(4) = Fext(4) +Ftotal(i*2-1)
                                Finte(free_index) = Ftotal(i*2)
                                free_index = free_index+1
                            else 
                                 Finte(free_index) = Ftotal(i*2-1)
                                 Finte(free_index+1) = Ftotal(i*2)
                                 free_index= free_index+2
                            end if 
                        end if
                    end do
                end if 
                Ext_diff = Fext-Fextexp
                Ext_diff_norm = norm2(Ext_diff)
                Finte_norm = norm2(Finte)*alpha
                Fextexp_norm = norm2(Fextexp)
                obj_list(f) = 0.d0 
                obj_list(f) = Finte_norm+Ext_diff_norm
                norm_list(f) = Fextexp_norm
				if (Fextexp_norm == 0.d0) then
					obj = obj
				else
					if (uni_penalty == 1) then
						if (flag3 == 1) then
							obj = obj+obj_list(f)/norm_list(f)*1.2d0
						else
							obj = obj+obj_list(f)/norm_list(f)
						endif
					else
						obj = obj+obj_list(f)/norm_list(f)
					endif
				endif
                deallocate(Finte)    
            end do
        endif
        if (isnan(obj)) then
			
            obj = obj_hardcode
        endif
        return
    end subroutine compute_MSE_ps
end module fem_compute_gl
        
        
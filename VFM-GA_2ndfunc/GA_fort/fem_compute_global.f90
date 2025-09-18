! This module corresponds to the subroutines that compute the objective function by summing up the internal 
!! and external residual force and summing up the external force measured in experiment    
!! The main part of the code is global assembly of elemental residual force and regrouping them into objective function 
module fem_compute_gl
    use constants
    use fem_compute_el
    implicit none
    contains
    !Obtain the force external measured in experiment from the input data
    !type 2 for instron based experiment.
    subroutine get_Fextexp(rf,Fextexp,free_count)
        implicit none
        integer i,free_count,nn_m_sqrt
        real(di) Fextexp(4),rf(nnode,2)
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
    !Copy of the previous function, for second set of dic experiment with field from a different undeformed mesh 
    !reason for the copy is to keep most things in stack 
    subroutine get_Fextexp2(rf,Fextexp,free_count)
        implicit none
        integer i,free_count,nn_m_sqrt
        real(di) Fextexp(4),rf(nnode2,2)
        nn_m_sqrt = nnode2-sqnn
        Fextexp = 0.d0
        free_count = 0
        if (bctype == 1) then
            do i = 1,nnode2
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
            do i = 1,nnode2
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
            do i = 1,nnode2
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
            do i = 1,nnode2
                if (i > (nn_m_sqrt)) then
                    Fextexp(2) = Fextexp(2)+rf(i,2)
                end if
                if (rf(i,1).NE.0.d0) then
                    Fextexp(1) = Fextexp(1)+rf(i,1)
                end if 
            end do
            Fextexp(3) = -Fextexp(1)
            Fextexp(4) = -Fextexp(2)
            free_count = nnode2*2-sqnn*2-sqnn-(sqnn-1)
        end if
        return
    end subroutine get_Fextexp2
     
    !! plane stress objective function computation for the field, you can modify if you have more than 2 set of fields that you want to fit together. 
    subroutine compute_MSE_ps(coord,disp,connect,rf,obj,theta,coord2,disp2,connect2,rf2,cmo,ck1) 
        use constants
        use fem_compute_el
        use stability_check
        implicit none
        real(di) obj,coord(2,nnode),disp(ncomp,nnode*2), &
&       rf_temp(nnode,2),Fextexp(4),disp_temp(nnode*2),disp_temp_old(nnode*2), &
&       Ftotal(dim*nnode),coord_el(dim,neln),disp_el(dim,neln),disp_old(dim,neln), &
&       theta(14),F_el(neln*dim),Fext(group), &
&       norm_list(nframe),Ext_diff(group),Fextexp_norm,obj_hardcode, &
&       Ext_diff_norm,Finte_norm,obj_list(nframe),rf(nnode,ncomp*2),F33_guess,F33_out,S22_out
        real(di) S22_now,S22_pre,S22_pre2,S22_diff,K1,F11_out,S_g_in(dim,dim,dim),S_g_out(dim,dim,dim)
        real(di) S_g_in_s(dim,dim,dim,nel),S_g_in_s2(dim,dim,dim,nel2),K1_old
        real(di) coord2(2,nnode2),disp2(nten,nnode2*2),rf_temp2(nnode2,2),disp_temp2(nnode2*2)
        real(di) Ftotal2(dim*nnode2),rf2(nnode2,nten*2),disp_temp2_old(nnode2*2)
        real(di),allocatable::Finte(:)
        real(di),allocatable::Finte2(:)
        integer connect(neln,nel),free_dof,rw,ff,connect2(neln,nel2)
        integer f,el,a,i,free_index,flag,flag2,flag3,cmo,ck1,stability
        
        obj = 0.d0
        obj_hardcode = 100000000000.d0
        obj_list = 0.d0
        S22_pre = 0.d0
        S22_pre2 = 0.d0
        flag2 = 0
        flag = 0
        S_g_in_s = 0.d0
        S_g_in_s2 = 0.d0
        stability = 0
        call sta_pack(theta,stability)
        if (stability == 1) then
            obj = obj_hardcode
        else
        do f = 1,nframe
            if (f<break) then
                flag3 = 0
                F33_guess = 1.d0
                Ftotal = 0.d0
                Fext = 0.d0
                Fextexp = 0.d0
                rf_temp(1:nnode,1:2) = rf(1:nnode,f*2-1:f*2)
                call get_Fextexp(rf_temp,Fextexp,free_dof)
                disp_temp(:) = disp(f,:)
                if (f == 1) then
                    disp_temp_old = 0.d0
                else
                    disp_temp_old(:) = disp(f-1,:)
                end if 
                
            
                do el=1,nel
                    do a = 1,neln
                        do i = 1,dim
                            disp_el(i,a) = disp_temp(dim*(connect(a,el)-1)+i)
                            disp_old(i,a) = disp_temp_old(dim*(connect(a,el)-1)+i)
                            coord_el(i,a) = coord(i,connect(a,el))
                        end do
                    end do
                    F_el = 0.d0
                    
                    S_g_in(:,:,:) = S_g_in_s(:,:,:,el)
                    if (el == 1) then
                    
                    call compute_Fel_ps_monoK1(coord_el,disp_el,theta,F_el,flag,F33_out,S22_out,K1,F11_out,S_g_in,S_g_out,disp_old)
                        
                        S_g_in_s(:,:,:,el) = S_g_out(:,:,:)
                        if (abs(F11_out-F33_out)>0.032d0) then
                            flag3 = 1 
                        endif
                        if (F33_out<1.d0) then
                            if (ck1 == 1) then
                                flag2 = 1
                            endif
                        endif
                        if (flag == 1) then
                            flag2 = 1
                        endif
                    
                        S22_now = S22_out
                        F33_guess = F33_out
                        if (S22_now>S22_pre) then
                            if (cmo == 1) then
                                flag2 = 1
                            endif
                        else
                            S22_pre = S22_now
                        endif 

                    else
                        if (rflag == 1) then 
                            !accuracy mod speed up by Reduced integ by Weican from Yuri's group
                    call compute_Fel_ps_monoK1(coord_el,disp_el,theta,F_el,flag,F33_out,S22_out,K1,F11_out,S_g_in,S_g_out,disp_old)
                            S_g_in_s(:,:,:,el) = S_g_out(:,:,:)
                        else
                            !less accurate mod more speedy using mixed type of plane stress solver
                            call compute_Fel_ps_after(coord_el,disp_el,theta,F_el,F33_guess,flag,F33_out,S_g_in,S_g_out,disp_old)
                        endif
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
                           ! print*,'node num',i
                        
                            Fext(1) = Fext(1)+Ftotal(i*2-1)
                            Fext(2) = Fext(2)+Ftotal(i*2)
                            !print*,Ftotal(i*2-1),Ftotal(i*2)
                        else if (rf_temp(i,2)<0.d0) then
                            !print*,'node num_b',i
                            Fext(3) = Fext(3)+Ftotal(i*2-1)
                            Fext(4) = Fext(4)+Ftotal(i*2)
                            !print*,Fext(3:4)
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
                if (uni_penalty == 1) then
                    if (flag3 == 1) then
                        obj = obj+obj_list(f)/norm_list(f)*1.2d0
                    else
                        obj = obj+obj_list(f)/norm_list(f)
                    endif
                else
                    obj = obj+obj_list(f)/norm_list(f)
                endif
                deallocate(Finte)
            else
		if (f == break) then
			K1_old = 0.d0
		endif
                flag3 = 0
                F33_guess = 1.d0
                Ftotal2 = 0.d0
                Fext = 0.d0
                Fextexp = 0.d0
                ff = f-break+1
                rf_temp2(1:nnode2,1:2) = rf2(1:nnode2,ff*2-1:ff*2)
                
                call get_Fextexp2(rf_temp2,Fextexp,free_dof)
                disp_temp2(:) = disp2(ff,:)
                if (ff == 1) then
                    disp_temp2_old = 0.d0
                else
                    disp_temp2_old(:) = disp2(ff-1,:)
                end if 
                do el=1,nel2
                    do a = 1,neln
                        do i = 1,dim
                            disp_el(i,a) = disp_temp2(dim*(connect2(a,el)-1)+i)
                            disp_old(i,a) = disp_temp2_old(dim*(connect2(a,el)-1)+i)
                            coord_el(i,a) = coord2(i,connect2(a,el))
                        end do
                    end do
                    F_el = 0.d0
                    S_g_in(:,:,:) = S_g_in_s2(:,:,:,el)
                   if (el == 1) then

                        
                    call compute_Fel_ps_monok1(coord_el,disp_el,theta,F_el,flag,F33_out,S22_out,K1,F11_out,S_g_in,S_g_out,disp_old)
                        
                        S_g_in_s2(:,:,:,el) = S_g_out(:,:,:)
			    if (K1_old.ge.K1) then
                                if (ck1 == 1) then
                                    flag2 = 1
                                end if
			    else
				K1_old = K1
                            endif
 
                        if (abs(F11_out-F33_out)>0.032d0) then
                            flag3 = 1
                        endif
                        if (flag == 1) then
                            flag2 = 1
                        endif
                        S22_now = S22_out
                        F33_guess = F33_out
                        S22_diff = S22_now-S22_pre2
                        if (S22_diff>0.d0) then
                            S22_pre2 = S22_now
                        else
                            !print*,'tension stress check '
                            if (ck1 == 1) then
                                flag2 = 1
                            endif
                        endif

                   else
                        if (rflag == 1) then 
                            !accuracy mod speed up by Reduced integ by Weican from Yuri's group
                    call compute_Fel_ps_monoK1(coord_el,disp_el,theta,F_el,flag,F33_out,S22_out,K1,F11_out,S_g_in,S_g_out,disp_old)
                            S_g_in_s(:,:,:,el) = S_g_out(:,:,:)
                        else
                            !less accurate mod more speedy using mixed type of plane stress solver
                            call compute_Fel_ps_after(coord_el,disp_el,theta,F_el,F33_guess,flag,F33_out,S_g_in,S_g_out,disp_old)
                        endif
                        if (flag == 1) then
                            flag2 = 1
                        endif
                        F33_guess = F33_out

                   end if
                    do a=1,neln
                        do i = 1,dim
                            rw = dim*(connect2(a,el)-1)+i
                            Ftotal2(rw) = Ftotal2(rw) + F_el(dim*(a-1)+i)
                        end do
                    end do
                end do
                if (flag2 == 1) then
                    obj = obj_hardcode
                    exit
                end if
                allocate(Finte2(free_dof))
            
                Finte2 = 0.d0
                Fext = 0.d0
                free_index =1
                if (bctype == 1) then
                    do i = 1,nnode2
                        if (rf_temp2(i,1)>0.d0) then
                            Fext(1) = Fext(1) +Ftotal2(i*2-1)
                        else if (rf_temp2(i,1)<0.d0) then
                            Fext(3) = Fext(3) +Ftotal2(i*2-1)
                        else
                            Finte2(free_index) = Ftotal2(i*2-1)
                            free_index = free_index+1
                        end if 
                        if (rf_temp2(i,2) > 0.d0) then
                            Fext(2) = Fext(2) +Ftotal2(i*2)
                        else if (rf_temp2(i,2)<0.d0) then
                            Fext(4) = Fext(4) +Ftotal2(i*2) 
                        else
                            Finte2(free_index) = Ftotal2(i*2)
                            free_index = free_index+1 
                        end if
                    end do
                else if (bctype == 2) then
                    do i = 1,nnode2
                        if (rf_temp2(i,2)>0.d0) then
                            Fext(1) = Fext(1)+Ftotal2(i*2-1)
                            Fext(2) = Fext(2)+Ftotal2(i*2)
                        else if (rf_temp2(i,2)<0.d0) then
                            Fext(3) = Fext(3)+Ftotal2(i*2-1)
                            Fext(4) = Fext(4)+Ftotal2(i*2)

                        else 
                            Finte2(free_index) = Ftotal2(i*2-1)
                            Finte2(free_index+1) = Ftotal2(i*2)
                            free_index = free_index+2
                        end if
                    end do
                else if (bctype == 3) then
                    do i = 1,nnode2
                        if (rf_temp2(i,2)<0.d0) then
                            Fext(1) = Fext(1)+Ftotal2(i*2-1)
                            Fext(2) = Fext(2)+Ftotal2(i*2)
                        else
                            if ((rf_temp2(i,2)>0.d0) .and. (rf_temp2(i,1)==0.d0)) then
                                Fext(3) = Fext(3) +Ftotal2(i*2)
                                Finte2(free_index) = Ftotal2(i*2-1)
                                free_index = free_index+1
                            else if ((rf_temp2(i,2)>0.d0) .and. (rf_temp2(i,1).NE.0.d0)) then
                                Fext(3) = Fext(3) +Ftotal2(i*2)
                                Fext(4) = Fext(4) +Ftotal2(i*2-1)
                            else if  ((rf_temp2(i,1).NE.0.d0).and. (rf_temp2(i,2)==0.d0)) then
                                Fext(4) = Fext(4) +Ftotal2(i*2-1)
                                Finte2(free_index) = Ftotal2(i*2)
                                free_index = free_index+1
                            else 
                                 Finte2(free_index) = Ftotal2(i*2-1)
                                 Finte2(free_index+1) = Ftotal2(i*2)
                                 free_index= free_index+2
                            end if 
                        end if
                    end do
                end if 
                Ext_diff = Fext-Fextexp
                Ext_diff_norm = norm2(Ext_diff)
                Finte_norm = norm2(Finte2)*alpha
                Fextexp_norm = norm2(Fextexp)
                obj_list(f) = 0.d0 
                obj_list(f) = Finte_norm+Ext_diff_norm
                !print*,obj_list(f)
                norm_list(f) = Fextexp_norm
                if (uni_penalty == 1) then
                    if (flag3 == 1) then
                        obj = obj+obj_list(f)/norm_list(f)*1.2d0*tweight
                    else
                        obj = obj+obj_list(f)/norm_list(f)*tweight
                    endif
                else
                    obj = obj+obj_list(f)/norm_list(f)*tweight
                endif
                deallocate(Finte2)
            end if
        end do
        endif
        if (isnan(obj)) then
            obj = obj_hardcode
        endif
        return
    end subroutine compute_MSE_ps

end module fem_compute_gl
        
        
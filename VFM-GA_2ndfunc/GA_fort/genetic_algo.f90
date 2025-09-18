!! Genetic Algorithm utility
module shinka
    use constants
    use fem_compute_gl
    implicit none
    contains
    !!binary decimal conversion
    subroutine t2t10(input,output)
        implicit none
        integer i,j
        integer(2) input(nparam*nbit)
        integer output(nparam),output_temp(nbit)
        do i = 1,nparam
            do j = 1,nbit
                output_temp(j) = input((i-1)*nbit+j)
            end do
            do j = 1,nbit
                output(i) = output(i)+output_temp(nbit-j+1)*2**(j-1)
            end do
        end do
        return
    end subroutine t2t10
    !! conversion from binary space with range to real space
    subroutine vtr(input,output)
        implicit none
        integer i,input(nparam)
        real(di) output(nparam)
        do i = 1,nparam
            output(i) = lowb(i)+input(i)*(upb(i)-lowb(i))/(2.d0**nbit-1)
        end do
        return
    end subroutine vtr
    
    !! compute objective function and pass into selection (roulette wheel/boltzmann)
    subroutine fitness(population,coord,disp_mat,connect,rf_mat,fitvalue,cumsump,gen,coord2,disp_mat2,connect2,rf_mat2)
        implicit none
        integer(2) population(nt,nbit*nparam)
        real(di) fitvalue(nt),base10r(nparam),fitsum,fitvalue_i(nt)
        real(di) perpo(nt),cumsump(nt),obj,coord(2,nnode), &
&               rf_mat(nnode,ncomp*2),disp_mat(ncomp,nnode*2)
        real(di) coord2(2,nnode2),rf_mat2(nnode2,nten*2),disp_mat2(nten,nnode2*2)
        integer connect2(neln,nel2)
        integer base10v(nparam),i,connect(neln,nel),gen,incline,plat1,plat2,j
        real(di) max_T,min_T,increment,T_temp
        do i = 1,nt
            base10v = 0
            base10r = 0.d0
            call t2t10(population(i,:),base10v)
            call vtr(base10v,base10r)
            if (phys_c_e == 0) then
                call compute_MSE_ps(coord,disp_mat,connect,rf_mat,obj,base10r,coord2,disp_mat2,connect2,rf_mat2,0,0)
            elseif (phys_c_e == 1) then
                call compute_MSE_ps(coord,disp_mat,connect,rf_mat,obj,base10r,coord2,disp_mat2,connect2,rf_mat2,1,0)
            elseif (phys_c_e == 2) then
                call compute_MSE_ps(coord,disp_mat,connect,rf_mat,obj,base10r,coord2,disp_mat2,connect2,rf_mat2,0,1)
            elseif (phys_c_e == 3) then
                call compute_MSE_ps(coord,disp_mat,connect,rf_mat,obj,base10r,coord2,disp_mat2,connect2,rf_mat2,1,1)
            else
                print*,'early stage physical check definition error, it seems that you wanna define your own physcal checks' 
            endif
            fitvalue(i) = obj
        end do
        
        if (control_boltz == 1) then 
            min_T = Tmin
            max_T = Tmax
            plat1 = 2
            plat2 = 2
            incline = ngen-plat1-plat2
            increment = (max_T-min_T)/dble(incline)
            if (gen < (plat1+1)) then
                T_temp = max_T
            else if (gen.ge.(plat1+1) .and. gen<(incline+plat1+1)) then
                T_temp = max_T-increment*(dble(gen-plat1))
            else
                T_temp = min_T
            endif 
            fitvalue_i = dexp(fitvalue/T_temp)
            do j = 1,nt
                if (fitvalue_i(j)>(10.d0)**(10.d0)) then
                    fitvalue_i(j) = 1000000000.d0
                endif
            end do
        else
        endif
        fitvalue_i = 1/fitvalue_i
        fitsum = sum(fitvalue_i)
        perpo = fitvalue_i/fitsum
        cumsump(1) = perpo(1)
        do i = 2,nt
            cumsump(i) = cumsump(i-1)+perpo(i)
        end do
        return
    end subroutine fitness
     
     !! compute objective function and pass into selection plus keep the best algo (not allowing divergence in any generation) (roulette wheel/boltzmann)
    subroutine fitness_m_ete(population,coord,disp_mat,connect,rf_mat,fitvalue,cumsump,gen,bvi,bpi,coord2,disp_mat2,connect2,rf_mat2)
        implicit none
        integer(2) population(nt,nbit*nparam),bpi(nbit*nparam)
        real(di) fitvalue(nt),base10r(nparam),fitsum,fitvalue_i(nt)
        real(di) perpo(nt),cumsump(nt),obj,coord(2,nnode), &
&               rf_mat(nnode,ncomp*2),disp_mat(ncomp,nnode*2)
        integer base10v(nparam),i,connect(neln,nel),gen,incline,plat1,plat2,j,mindex,maxdex
        real(di) coord2(2,nnode2),rf_mat2(nnode2,nten*2),disp_mat2(nten,nnode2*2)
        integer connect2(neln,nel2)
        real(di) max_T,min_T,increment,T_temp,bvi
        do i = 1,nt
            base10v = 0
            base10r = 0.d0
            call t2t10(population(i,:),base10v)
            call vtr(base10v,base10r)
            if (phys_c_l == 0) then
                call compute_MSE_ps(coord,disp_mat,connect,rf_mat,obj,base10r,coord2,disp_mat2,connect2,rf_mat2,0,0)
            elseif (phys_c_l == 1) then
                call compute_MSE_ps(coord,disp_mat,connect,rf_mat,obj,base10r,coord2,disp_mat2,connect2,rf_mat2,1,0)
            elseif (phys_c_l == 2) then
                call compute_MSE_ps(coord,disp_mat,connect,rf_mat,obj,base10r,coord2,disp_mat2,connect2,rf_mat2,0,1)
            elseif (phys_c_l == 3) then
                call compute_MSE_ps(coord,disp_mat,connect,rf_mat,obj,base10r,coord2,disp_mat2,connect2,rf_mat2,1,1)
            else
                print*,'late stage physical check definition error, it seems that you wanna define your own physcal checks' 
            endif
            fitvalue(i) = obj
        end do
        mindex = minloc(fitvalue,DIM =1)
        if (fitvalue(mindex)>bvi) then
            maxdex = maxloc(fitvalue,DIM=1)
            population(maxdex,:) = bpi(:)
            fitvalue(maxdex) = bvi
        end if
        if (control_boltz == 1) then
            min_T = Tmin
            max_T = Tmax
            plat1 = 2
            plat2 = 2
            incline = ngen-plat1-plat2
            increment = (max_T-min_T)/dble(incline)
            if (gen < (plat1+1)) then
                T_temp = max_T
            else if (gen.ge.(plat1+1) .and. gen<(incline+plat1+1)) then
                T_temp = max_T-increment*(dble(gen-plat1))
            else
                T_temp = min_T
            endif 
            fitvalue_i = dexp(fitvalue/T_temp)
            do j = 1,nt
                if (fitvalue_i(j)>(10.d0)**(10.d0)) then
                    fitvalue_i(j) = 1000000000.d0
                endif
            end do
        else
        endif
        fitvalue_i = 1/fitvalue_i
        fitsum = sum(fitvalue_i)
        perpo = fitvalue_i/fitsum
        cumsump(1) = perpo(1)
        do i = 2,nt
            cumsump(i) = cumsump(i-1)+perpo(i)
        end do
        return
    end subroutine fitness_m_ete
    
     !! generate population
    subroutine genpop(output)
        implicit none
        real(di) temp_rand(nt,nbit*nparam)
        integer(2) output(nt,nbit*nparam)
        !call RANDOM_SEED()
        call RANDOM_NUMBER(temp_rand)
        output = int(temp_rand+0.5)
        return
    end subroutine genpop
    
    !! a simple selection that is used in roulette wheel/boltzmann
    subroutine selection(cumsump,selnum)
        integer selnum(2),i,j
        real(di) cumsump(nt),r
        do i = 1,2
            r=0.d0
            call random_number(r)
            j = 1
            do while(cumsump(j)<r)
                j = j+1
            end do
            selnum(i) = j
        end do
        return 
    end subroutine selection
    
     !! a random processor that given the probability of an event, returns if event happens
    subroutine ifxorm(p,tf)
        implicit none
        real(di) p,r
        integer(2) tf,test(100)
        integer tempint,tempint2
        tempint = int(100*p)
        test(1:tempint) = 1
        test(tempint+1:100) = 0
        !call RANDOM_SEED()
        call RANDOM_NUMBER(r)
        tf = test(int(r*99)+1)
        return         
    end subroutine ifxorm
    
    !! bitwise cross-over,3 modes available, use the default single point cross-over is enough
    subroutine xoverbit(population,seln,parents)
        implicit none
        integer(2) tf,parents(2,nbit*nparam)
        integer(2) population(nt,nbit*nparam)
        integer cindex,sindex,eindex,seln(2),i,choice,j
        real(di) xnum(nparam),r,r1,r2,r3,r4
        call ifxorm(xrate,tf)
        parents(1,:) = population(seln(1),:)
        parents(2,:) = population(seln(2),:)
        if (tf == 1) then 
            !call RANDOM_SEED()
            
            if (cross_mode == 1) then
                call RANDOM_NUMBER(xnum)
                do i = 1,nparam
                    if (xcon>xnum(i)) then
                        
                        call r_choice_N(nbit-1,choice)
                        cindex = choice + (i-1)*nbit
                        sindex = 1+(i-1)*nbit
                        eindex = i*nbit
                        parents(1,sindex:cindex) = population(seln(1),sindex:cindex) 
                        parents(1,cindex+1:eindex) = population(seln(2),cindex+1:eindex)
                        parents(2,sindex:cindex) = population(seln(2),sindex:cindex) 
                        parents(2,cindex+1:eindex) = population(seln(1),cindex+1:eindex)   
                    end if 
                end do
            elseif (cross_mode == 2) then
                call RANDOM_NUMBER(xnum)
                do i = 1,nparam
                    if (xcon>xnum(i)) then
                        !call RANDOM_SEED()
                        do j = 1,nbit
                            call RANDOM_NUMBER(r1)
                            call RANDOM_NUMBER(r2)
                            call RANDOM_NUMBER(r3)
                            call RANDOM_NUMBER(r4)
                            if (r1>r2) then
                                parents(1,j+(i-1)*nbit) = population(seln(1),j+(i-1)*nbit) 
                            else
                                parents(1,j+(i-1)*nbit) = population(seln(2),j+(i-1)*nbit) 
                            end if
                            if (r3>r4) then
                                parents(2,j+(i-1)*nbit) = population(seln(1),j+(i-1)*nbit) 
                            else 
                                parents(2,j+(i-1)*nbit) = population(seln(2),j+(i-1)*nbit) 
                            end if
                        end do
                    end if 
                end do
            elseif (cross_mode == 3) then 
                call RANDOM_NUMBER(xnum)
                do i = 1,nparam
                    if (xcon>xnum(i)) then
                        do j = 1,nbit
                            call RANDOM_number(r)
                            if (r.ge.0.500000d0) then
                                parents(1,j+(i-1)*nbit) = population(seln(2),j+(i-1)*nbit)
                                parents(2,j+(i-1)*nbit) = population(seln(1),j+(i-1)*nbit)
                            end if
                        end do
                    end if 
                end do
            else
                print*,'you fxxk up with crossover mode number'
            endif
        end if 
        return 
    end subroutine xoverbit
    
    !! bitwise mutation,2 modes available, use the default simple inversion mutation is enough
    subroutine mutation(input,output)
        implicit none
        integer(2) input(nbit*nparam),output(nbit*nparam),tf,i,j
        integer index,choice
        real(di) r 
        if (mut_mode == 1) then
            output = input 
            call ifxorm(mrate,tf)
            if (tf == 1) then 
                do i = 1,nparam
                    call RANDOM_NUMBER(r)
                    if (mcon>r) then
                        do j = 1,mnum
                            call r_choice_N(nbit, choice)
                            index = choice+(i-1)*nbit
                            output(index) = abs(input(index)-1)
                        end do 
                    end if 
                end do
            end if 
        elseif (mut_mode == 2) then
            output = input 
            call ifxorm(mrate,tf)
            if (tf == 1) then 
                do i = 1,nparam
                    !call RANDOM_SEED()
                    call RANDOM_NUMBER(r)
                    if (mcon>r) then
                        do j = 1,mnum
                            !call RANDOM_SEED()
                            call RANDOM_NUMBER(r)
                            index = int(r*(nbit-1)+0.5)+1+(i-1)*nbit
                            output(index) = abs(input(index)-1)
                        end do 
                    end if 
                end do
            end if 
        else
            print*,'you fxxk up with mutation mode number'
        endif
        return
    end subroutine mutation
    
     !! For fun, not useful in ifort 
    subroutine init_random_seed()
        INTEGER :: i, n, clock
        INTEGER, DIMENSION(:), ALLOCATABLE :: seed
        real(di) output
          
        CALL RANDOM_SEED(size = n)
        ALLOCATE(seed(n))
          
        CALL SYSTEM_CLOCK(COUNT=clock)
          
        seed = clock + 37*(/ (i-1 , i = 1, n) /)
        CALL RANDOM_SEED(PUT = seed)
        DEALLOCATE(seed)
    end subroutine init_random_seed
    
    !! subroutine that takes average
    subroutine average(input,output)
        implicit none
        integer i 
        real(di) input(nt),output
        output = 0.d0
        do i = 1,nt
            output = output+input(i)
        end do
        output = output/nt
        return
    end subroutine average
    
    !!utility to randomly choose a number
    subroutine r_choice(output)
        implicit none
        real(di) r,division(nbit+1)
        integer output,i
        division(1) = 0.d0
        do i = 1,nbit
            division(i+1) = (1.d0/dble(nbit))*dble(i)
        end do
        call RANDOM_NUMBER(r)
        i = 1
        do while(division(i+1)<r)
            i = i+1
        end do
        output = i
        return
    end subroutine r_choice
    
    !! randomly output a integer out of N possible integer
    subroutine r_choice_N(N, output)
        implicit none
        integer, intent(in) :: N
        integer, intent(out) :: output
        real(di) r, delta
        integer i

        call RANDOM_NUMBER(r)
        
        delta = 1.0d0 / N
        output = N
        do i = 1, N - 1
            if (r.le.i * delta) then
                output = i
                exit
            end if
        end do
        return
    end subroutine r_choice_N
    
    !!similar to the before one, for fun only, choose between 1 to 4
    subroutine r_choice_4(output)
        implicit none
        integer output
        real(di) r
        call RANDOM_NUMBER(r)
        output = 1
        if (r.le.0.25d0) then
            output = 1
        else if (r.gt.0.25d0 .and. r.le.0.50d0) then
            output = 2
        else if (r.gt.0.50d0 .and. r.le.0.75d0) then
            output = 3
        else
            output = 4
        end if
        return
    end subroutine r_choice_4 
    
     !!similar to the before one, for fun only, choose between 1 to 3
    subroutine r_choice_3(output)
        implicit none
        integer output
        real(di) r
        call RANDOM_NUMBER(r)
        output = 1
        if (r.le.0.33333333d0) then
            output = 1
        else if (r.gt.0.3333333d0 .and. r.le.0.666666d0) then
            output = 2
        else
            output = 3
        end if
        return
    end subroutine r_choice_3
    
end module shinka
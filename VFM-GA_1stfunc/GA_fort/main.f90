!! '-------------GA hyperelastic V1.0-------------'
!! '----------------1st functionality-------------'
!! '------------------Zicheng Yan-----------------'
!! '---------zicheng_yan@brown.edu----------------'
!! '-----Henann Research Group, Brown University--'
program ga1st
    use fembasics
    use linalg
    use fem_constitutive
    use fem_compute_el
    use fem_compute_gl
    use fem_pre
    use constants
    use shinka
    use omp_lib
    use stability_check
    use print_write
    
    !! 0. This code corresponds to the 1st functionality in the paper, please cite if used.
    !! The input user needs to take care of is stress-strain curve from a homogeneous compression/tension experiment
    !! Go back to the previous folder and read S_1 S_2 S_3 matlab files sequentially for details.
    
    
    !! 1. default calibration with our model 
    !!To calibrate our foam model with your own data: https://www.sciencedirect.com/science/article/pii/S0022509619303825
    !!Simply follow and run the S_1 S_2 and S_3 matlab "driver" files.
    
    !! 2. calibration for other hyperelastic model
    !!This framework is very powerful on calibrating of all other hyperelastic model (e.g. HGO)
    !!To apply it to your own model, you have to write your own subroutine in fem_const.f90, the subroutine needs to take 2 inputs: deformation gradient and random material parameter set theta
    !!the subroutine need to have 1 output: kirchoff stress. You will also need to turn off the Dacorogna stability check, or derive and write your own version. Some parameters like
    !!nparam (number of material parameters) needs to be changed accordingly in constants.f90. Refer to S_2 matlab file for meaning of parameters and guidance on how to write your own model in
    !!Also substantial change on main.f90 is needed
    
    !! 3. This code is designed to run parrallelly completely based on external parallelization due to nature of genetic algorithm
    !!Trial on omp and mpi leads to speed loss
    implicit none
    real(di) disp_mat(nframe,nnode*2),fitvalue(nt)
    real(di) coord(2,nnode),rf_mat(nnode,nframe*2),theta_ref(nparam)
    integer connect(neln,nel),i,gen,j,seln(2),mindex
    real(di) lam,theta(nparam),obj,start,finish,cumsump(nt),avg
    integer theta_v(nparam)
    integer(2) pop(nt,nbit*nparam),parents(2,nbit*nparam),bpi(nbit*nparam)
    integer(2) mnew(nt,nbit*nparam)
    real(di) obj_min(ngen),theta_min(ngen,nparam),error(nparam),bvi
    real(di) theta_r(nparam),best_objs(npop),best_theta(npop,nparam)
    integer rate,stability,stat,values(8),seed(2),addon
    character val*255,vv,vv2
    CALL system_clock(count_rate=rate)
    !! Process the input and make file for output
    if (nframe<101) then
        call get_u_rf(rf_mat,disp_mat)
    else
        call get_u_rf_100(rf_mat,disp_mat)
    end if
    call get_n_el(coord,connect)
    
    open(unit=1, iostat=stat, file='ga_out.txt', status='old')
    if (stat == 0) then
        print*,'Old record file exists, erase old file, oops' 
        close(1, status='delete')
    endif
    open(unit=2, iostat=stat, file='convconv.txt', status='old')
    if (stat == 0) then
        print*,'Old conv file exists, erase old file, oops' 
        print*,''
        close(2, status='delete')
    endif
    !! output for material parameters and proof of convergence
    open(1,file='ga_out.txt',status='new')
    open(2,file='convconv.txt',status='new')
    call start_print()
    call start_write()
     
    !! To handle the nasty nasty ifort seeding both locally and on cluster
    !! requires the seed to be taken from 1.time and 2.folder name
    !! fail to do so could generate same seeds which means identical random process
    call  getcwd(val)
    vv = val(LEN(TRIM(val))-1:LEN(TRIM(val))-1)
    if (vv.eq.'t') then
        read (val(LEN(TRIM(val)):LEN(TRIM(val))),'(i4)') addon
    else
	read (val(LEN(TRIM(val))-1:LEN(TRIM(val))),'(i4)') addon
    endif
    call DATE_AND_TIME(values=values)
    seed(:) = values(7:8)
    seed(1) = seed(1)+addon
    call RANDOM_SEED(put=seed)

    !! run_mode 0 is the debug mode, just to see a how given material parameter behaves
    if (run_mode == 0) then
        theta_ref = (/102000.d0,193800.d0,0.19d0,1.9d0,-0.21d0,0.2d0, & 
        &               3.7d0,0.22d0,0.1d0,4.d0,5.d0,9.d0,0.026d0,2.d0/)
        print*,'theta_ref',theta_ref
        call cpu_time(start)
        do i = 1,1
        call compute_MSE_ps(coord,disp_mat,connect,rf_mat,obj,theta_ref,1,1)
        call sta_pack(theta_ref,stability)
        end do 
        call cpu_time(finish)
        print*,'referential objective function value ',obj
        print*,'stability check',stability
        print '("Time = ",f6.3," seconds.")',finish-start
    !! Define your own run mode 1 for detailed debugging/testing on different parts of GA 
        
    else if (run_mode == 1) then
        !do i = 1,npop
        !call genpop(pop)
        !call fitness(pop,coord,disp_mat,connect,rf_mat,fitvalue,cumsump,0)
        !do j = 1,20
        !    theta_v = 0.d0
        !    theta_r = 0.d0
        !    call t2t10(pop(j,:),theta_v)
        !    call vtr(theta_v,theta_r)
        !end do
        !
        !call average(fitvalue,avg)
        !print*,'best in gen0 ',minval(fitvalue)
        !print*,'mean in gen0 ',avg
        !end do
        
    !! Run mode 2 is the formal run mode, which includes the main genetic algorithm optimization process
    else if (run_mode == 2) then
    do i = 1,npop
        call genpop(pop)
        call fitness(pop,coord,disp_mat,connect,rf_mat,fitvalue,cumsump,0)
        call average(fitvalue,avg)
        call gen0_print_write(i,minval(fitvalue),avg)
        gen = 1
        do while (gen<ngen+1)
            do j = 1,nt,2
                seln = 0
                call selection(cumsump,seln)
                call xoverbit(pop,seln,parents)
                call mutation(parents(1,:),mnew(j,:))
                call mutation(parents(2,:),mnew(j+1,:))
            end do
            pop = 0
            pop = mnew
            if (gen<2) then
                call fitness(pop,coord,disp_mat,connect,rf_mat,fitvalue,cumsump,gen)
            else
                if (elitism == 1) then
                    call fitness_m_ete(pop,coord,disp_mat,connect,rf_mat,fitvalue,cumsump,gen,bvi,bpi)
                else
                    call fitness(pop,coord,disp_mat,connect,rf_mat,fitvalue,cumsump,gen)
                endif
            end if
            mindex = minloc(fitvalue,DIM =1)
            avg = 0.d0
            call average(fitvalue,avg)
            obj_min(gen) = minval(fitvalue)
            bvi = minval(fitvalue)
            bpi = pop(mindex,:)
            theta_v = 0
            theta_r = 0.d0
            call t2t10(pop(mindex,:),theta_v)
            call vtr(theta_v,theta_r)
            theta_min(gen,:) = theta_r
            call gen_print_write(gen,obj_min(gen),avg,theta_r)
            gen = gen+1
        end do
        best_objs(i) = minval(obj_min)
        best_theta(i,:) = theta_min(minloc(obj_min,DIM = 1),:)
        call after_print_write(best_objs(i),theta_min(minloc(obj_min,DIM = 1),:))
    end do
    close(1)
    end if 
end program ga1st
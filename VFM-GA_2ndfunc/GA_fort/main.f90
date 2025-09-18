!! '-------------GA hyperelastic V1.0-------------'
!! '----------------2nd functionality-------------'
!! '------------------Zicheng Yan-----------------'
!! '---------zicheng_yan@brown.edu----------------'
!! '-----Henann Research Group, Brown University--'
program ga2nd
    use constants
    use fembasics
    use linalg
    use fem_constitutive
    use fem_compute_el
    use fem_compute_gl
    use fem_pre
    use shinka
    use omp_lib
    use print_write
    !! 0. This code corresponds to the 2nd functionality in the paper, please cite if used.
    !! 2nd functionality framework implemented take up to 2 dic experiments, with 2 undeformed mesh and corresponding displacement fields, reaction forces
    
    !! 1. default calibration with our model 
    !!To calibrate our foam model with your own data: https://www.sciencedirect.com/science/article/pii/S0022509619303825
    !!Follow / run the S_0 S_1 and S_2 matlab "driver" files.
    !!Compared to 1st functionality, you have to take care of the data structure of your undeformed mesh and fields input according to S_0 matlab file 
    !!You also need to adjust part of S_1 matlab file (like number of fields input, number of nodes, number of elements) to print out a compatible constants.f90
    
    !! 2. calibration for other hyperelastic model
    !!This framework is very powerful on calibrating of all other hyperelastic model (e.g. HGO)
    !!To apply it to your own model, you have to write your own subroutine in fem_const.f90, the subroutine needs to take 2 inputs: deformation gradient and random material parameter set theta
    !!the subroutine need to have 1 output: kirchoff stress. You will also need to turn off the Dacorogna stability check, or derive and write your own version. Some parameters like
    !!nparam (number of material parameters) needs to be changed accordingly in constants.f90. Refer to S_1 matlab file for meaning of parameters and guidance on how to write your own model in
    !!Also substantial change on main.f90 is needed
    
    !! 3. calibrate multiple experiments (>2) together
    !! you have to modify the fem_compute_global.f90 and fem_preprocess.f90 and potentially other parts to enable the code to take in more undeformed meshes and displacement fields+forces.
    !! The code is hardcoded in a heavy stack way to avoid speed loss. 
    
    !! 4. This code is designed to run in parrallel completely based on external parallelization due to nature of genetic algorithm
    !!Trial on omp and mpi leads to speed loss
    
    implicit none
    real(di) disp_mat(ncomp,nnode*2),fitvalue(nt),Jmin,error_str
    real(di) coord(2,nnode),rf_mat(nnode,ncomp*2),theta_ref(nparam)
    integer connect(neln,nel),i,gen,j,seln(2),mindex,jjj
    real(di) theta(nparam),obj,start,finish,cumsump(nt),avg,F(3,3),stress(3,3)
    integer theta_v(nparam),jj,nf 
    integer(2) pop(nt,nbit*nparam),parents(2,nbit*nparam),bpi(nbit*nparam)
    integer(2) mnew(nt,nbit*nparam)
    real(di) obj_min(ngen),theta_min(ngen,nparam),error(nparam),bvi
    real(di) theta_r(nparam),best_objs(npop),best_theta(npop,nparam)
    real(di) disp_mat2(nten,nnode2*2),rf_mat2(nnode2,nten*2),coord2(2,nnode2)
    integer nthreads,mythread,num_threads,stat,values(8),seed(2)
    integer s1,s2,rate,thread,connect2(neln,nel2),addon
    character val*255,vv,vv2
    !! Process the input and make file for output
    CALL system_clock(count_rate=rate)
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
    open(1,file='ga_out.txt',status='new')
    open(2,file='convconv.txt',status='new')
    call start_print()
    call start_write()
    call get_u_rf(rf_mat,disp_mat)
    call get_n_el(coord,connect)
    if (nframe>ncomp) then
        call get_u_rf2(rf_mat2,disp_mat2)
        call get_n_el2(coord2,connect2)
    endif
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
        print*,'mode debug'
        theta_ref = (/102000.d0,193800.d0,0.1667d0,1.24d0,-0.3d0,0.1667d0, & 
&       2.5333d0,0.274d0,0.05d0,6.8d0,6.d0,10.333d0,0.0675d0,1.6667d0/)
        call cpu_time(start)
        do i = 1,10
            call compute_MSE_ps(coord,disp_mat,connect,rf_mat,obj,theta_ref,coord2,disp_mat2,connect2,rf_mat2,1,1)
        end do
        call cpu_time(finish)
        print '("Time = ",f6.3," seconds.")',finish-start
        print*,'referential objective function value ',obj
    !! Run mode 2 is the formal run mode, which includes the main genetic algorithm optimization process
    else if (run_mode == 2) then
    do i = 1,npop
        call genpop(pop)
        call fitness(pop,coord,disp_mat,connect,rf_mat,fitvalue,cumsump,0,coord2,disp_mat2,connect2,rf_mat2)
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
                call fitness(pop,coord,disp_mat,connect,rf_mat,fitvalue,cumsump,gen,coord2,disp_mat2,connect2,rf_mat2)
            else
                if (elitism == 1) then
            call fitness_m_ete(pop,coord,disp_mat,connect,rf_mat,fitvalue,cumsump,gen,bvi,bpi,coord2,disp_mat2,connect2,rf_mat2)
                else
            call fitness(pop,coord,disp_mat,connect,rf_mat,fitvalue,cumsump,gen,coord2,disp_mat2,connect2,rf_mat2)
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
        call after_print_write(best_objs(i),theta)
    end do
    close(1)
    else
        print*,'run_mode number not available, go define your own run_mode'
    end if 
end program ga2nd
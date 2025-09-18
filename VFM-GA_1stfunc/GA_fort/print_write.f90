!!Print and write part to keep track of what's going on in genetic algorithm
module print_write
    use constants
    implicit none
    contains
        subroutine start_print()
        integer values(8)
        character(len = 8) :: date
        character(len = 10) :: time
        character(len = 5) :: zone
        
        call date_and_time(date,time,zone,values)
        print*, '-------------GA hyperelastic V1.0--------------'
        print*, '------------------Zicheng Yan------------------'
        print*, '-----Henann Research Group, Brown University---'
        print*, ''
        print*,'Seireki ',date(1:4),'-',date(5:6),'-',date(7:8),' ',time(1:2),':',time(3:4),':',time(5:6)
        print*,''
        print*,''
        print*,''
        print*,'-------------------------------------------------------------------'
        print*,'----------------------start of info session------------------------'
        print*,'-------------------------------------------------------------------'
        print*,'___________________________________________________________________'
        print*,'A. Priority info of current run:'
        print*,''
        if (run_mode == 0) then
            print*,'0. Code Mode:Debug'
        elseif (run_mode == 1) then
            print*,'0. Code Mode:Custom Debug'
        elseif (run_mode == 2) then
            print*,'0. Code Mode:Run'
        endif
        if (control_ps == 2) then
            print*,'1. Plane stress 2D approximation of data'
        elseif (control_ps == 1) then
            print*,'1. Plane strain 2D approximation of data'
        endif
        print*,'2. Constitutive model: https://github.com/HenannResearchGroup/ElastomericFoam'
        print*,'3. parameter number to optimize',nparam
        print*,'4. alpha alpha between internal and external force',alpha
        print*,'5. precision in computation, bits per real variable',di
        print*,'___________________________________________________________________'
        print*,'B. GA related parameters:'
        print*,''
        print*,'1. nbits in binary formulation',nbit
        print*,'2. number of independent population',npop
        print*,'3. size of each population',nt
        print*,'4. total generation num: ',ngen
        if (control_boltz == 1) then
            print*,'5. Selection operator: Boltzmann'
            print*,'5.1. Max Boltzmann Temperature', Tmax
            print*,'5.2. Min Boltzmann Temperature', Tmin
        elseif (control_boltz == 0) then
            print*,'5. Selection operator: Roulette Wheel'
        endif
        
        if (cross_mode == 1) then
            print*,'6. Crossover operator: Single Point'
        elseif (cross_mode == 2) then
            print*,'6. Crossover operator: Shuffle (it seems that you like to explore)'
        elseif (cross_mode == 3) then
            print*,'6. Crossover operator: Uniform (it seems that you like to explore)'
        endif
        print*,'6.1 probability of crossover on para set xrate:',xrate
        print*,'6.2 probability of crossover on a single para xcon:',xcon
        
        if (mut_mode == 1) then
            print*,'7. Mutation operator: single'
        elseif (mut_mode == 2) then
            print*,'7. Mutation operator: rigged (it seems that you like to explore)'
        endif
        print*,'7.1 probability of mutation on para set mrate:',mrate
        print*,'7.2 probability of mutation on a single para mcon:',mcon
        print*,'7.3 number of bits in a para that mutates:',mnum
        print*,'___________________________________________________________________'
        print*,'C. Physical checks and stability parameters'
        print*,''
        if (elitism == 1) then
            print*,'1. Elitism activated'
        else
            print*,'1. Elitism deactivated'
        endif
        print*,'2. Dacorogna stability check'
        print*,'2.1 for pure shear, K1 = 0, K2 =',s_sta
        print*,''
        print*,'2.2 for simple compression, check point K1 values:'
        print*, c_sta_k1
        print*,'2.2 for simple compression, check point K2 values:'
        print*, c_sta_k2
        print*,''
        print*,'3. Rough check (only applicable for simple comp and ten)'
        if (phys_c_e == 0) then
            print*,'3.1. first generation checks: none'
        elseif (phys_c_e == 1) then 
            print*,'3.1. first generation checks: monotonic (rough)'
        elseif (phys_c_e == 2) then
            print*,'3.1. first generation checks: volumetric (rough)'
        elseif (phys_c_e == 3) then
            print*,'3.1. first generation checks: monotonic and volumetric (rough)'
        else
        endif
        if (phys_c_l == 0) then
            print*,'3.2. Later generation checks: none'
        elseif (phys_c_l == 1) then 
            print*,'3.2. Later generation checks: monotonic (rough)'
        elseif (phys_c_l == 2) then
            print*,'3.2. Later generation checks: volumetric (rough)'
        elseif (phys_c_l == 3) then
            print*,'3.2. Later generation checks: monotonic and volumetric (rough)'
        else
        endif
        if (control_ps == 2) then
            if (uni_penalty == 1) then
                print*,'4. penalty for plane stress computing different out of plane Fs activated'
            endif
        endif
        print*,'___________________________________________________________________'
        print*,'D. Mechanical data input and constitutive para range'
        print*,''
        print*,'1. parameter number to optimize',nparam
        print*,'2. input mesh nnode:',nnode
        print*,'3. input mesh nel:',nel
        print*,'4. input mesh node per element:',neln
        print*,'5. integration points:',intn
        print*,'6. input number of frames:',nframe
        print*,'7. number of compression frames',break-1
        print*,'8.1 material parameter lower range'
        print*,lowb(1:nparam)
        print*,'8.2 material parameter upper range'
        print*,upb(1:nparam)
        if (bctype == 1) then
            print*,'9. BC type: four face fictional (comp test only)'
        elseif (bctype == 2) then
            print*,'9. BC type: simple compression-tension'
        elseif (bctype == 3) then
            print*,'9. BC type: three face fictional (comp test only)'
        elseif (bctype == 4) then
            print*,'9. BC type: indentation'
        endif
        print*,'10. Ext Int Force vector size:',group
        print*,'-------------------------------------------------------------------'
        print*,'------------------------end of info session------------------------'
        print*,'-------------------------------------------------------------------'
        print*,''
        print*,''
        print*,'-------------------------------------------------------------------'
        print*,'----------------------------start computing------------------------'
        print*,'-------------------------------------------------------------------'
    end subroutine start_print
        
    subroutine start_write()
        integer values(8)
        character(len = 8) :: date
        character(len = 10) :: time
        character(len = 5) :: zone
        call date_and_time(date,time,zone,values)
        write(1,*) '-------------GA hyperelastic V1.0----------------'
        write(1,*) '--------------------Zicheng Yan------------------'
        write(1,*) '-----Henann Research Group, Brown University-----'
        write(1,*) ''
        write(1,*)'Seireki ',date(1:4),'-',date(5:6),'-',date(7:8),' ',time(1:2),':',time(3:4),':',time(5:6)
        write(1,*)''
        write(1,*)''
        write(1,*)''
        write(1,*)'-------------------------------------------------------------------'
        write(1,*)'----------------------start of info session------------------------'
        write(1,*)'-------------------------------------------------------------------'
        write(1,*)'___________________________________________________________________'
        write(1,*)'A. Priority info of current run:'
        write(1,*)''
        if (control_ps == 2) then
            write(1,*)'1. Plane stress approximation of data'
        elseif (control_ps == 1) then
            write(1,*)'1. Plane strain approximation of data'
        endif
        write(1,*)'2. Constitutive model: https://github.com/HenannResearchGroup/ElastomericFoam'
        write(1,*)'3. parameter number to optimize',nparam
        write(1,*)'4. alpha alpha between internal and external force',alpha
        write(1,*)'___________________________________________________________________'
        write(1,*)'B. GA related parameters:'
        write(1,*)''
        write(1,*)'1. nbits in binary formulation',nbit
        write(1,*)'2. number of independent population',npop
        write(1,*)'3. size of each population',nt
        write(1,*)'4. total generation num: ',ngen
        if (control_boltz == 1) then
            write(1,*)'5. Selection operator: Boltzmann'
            write(1,*)'5.1. Max Boltzmann Temperature', Tmax
            write(1,*)'5.2. Min Boltzmann Temperature', Tmin
        elseif (control_boltz == 0) then
            write(1,*)'5. Selection operator: Roulette Wheel'
        endif
        
        if (cross_mode == 1) then
            write(1,*)'6. Crossover operator: Single Point'
        elseif (cross_mode == 2) then
            write(1,*)'6. Crossover operator: Shuffle (it seems that you like to explore)'
        elseif (cross_mode == 3) then
            write(1,*)'6. Crossover operator: Uniform (it seems that you like to explore)'
        endif
        write(1,*)'6.1 probability of crossover on para set xrate:',xrate
        write(1,*)'6.2 probability of crossover on a single para xcon:',xcon
        
        if (mut_mode == 1) then
            write(1,*)'7. Mutation operator: single'
        elseif (mut_mode == 2) then
            write(1,*)'7. Mutation operator: rigged (it seems that you like to explore)'
        endif
        write(1,*)'7.1 probability of mutation on para set mrate:',mrate
        write(1,*)'7.2 probability of mutation on a single para mcon:',mcon
        write(1,*)'7.3 number of bits in a para that mutates:',mnum
        write(1,*)'___________________________________________________________________'
        write(1,*)'C. Physical checks and stability parameters'
        write(1,*)''
        if (elitism == 1) then
            write(1,*)'1. Elitism activated'
        else
            write(1,*)'1. Elitism deactivated'
        endif
        write(1,*)'2. Dacorogna stability check'
        write(1,*)'2.1 for pure shear, K1 = 0, K2 =',s_sta
        write(1,*)''
        write(1,*)'2.2 for simple compression, check point K1 values:'
        write(1,*) c_sta_k1
        write(1,*)'2.2 for simple compression, check point K2 values:'
        write(1,*) c_sta_k2
        write(1,*)''
        write(1,*)'3. Rough check (only applicable for simple comp and ten)'
        if (phys_c_e == 0) then
            write(1,*)'3.1. first generation checks: none'
        elseif (phys_c_e == 1) then 
            write(1,*)'3.1. first generation checks: monotonic (rough)'
        elseif (phys_c_e == 2) then
            write(1,*)'3.1. first generation checks: volumetric (rough)'
        elseif (phys_c_e == 3) then
            write(1,*)'3.1. first generation checks: monotonic and volumetric (rough)'
        else
        endif
        if (phys_c_l == 0) then
            write(1,*)'3.2. Later generation checks: none'
        elseif (phys_c_l == 1) then 
            write(1,*)'3.2. Later generation checks: monotonic (rough)'
        elseif (phys_c_l == 2) then
            write(1,*)'3.2. Later generation checks: volumetric (rough)'
        elseif (phys_c_l == 3) then
            write(1,*)'3.2. Later generation checks: monotonic and volumetric (rough)'
        else
        endif
        if (control_ps== 2) then
            if (uni_penalty == 1) then
                write(1,*)'4. penalty for plane stress computing different out of plane Fs activated'
            endif
        endif
        write(1,*)'___________________________________________________________________'
        write(1,*)'D. Mechanical data input and constitutive para range'
        write(1,*)''
        write(1,*)'1. parameter number to optimize',nparam
        write(1,*)'2. input mesh nnode:',nnode
        write(1,*)'3. input mesh nel:',nel
        write(1,*)'4. input mesh node per element:',neln
        write(1,*)'5. integration points:',intn
        write(1,*)'6. input number of frames:',nframe
        write(1,*)'7. number of compression frames',break-1
        write(1,*)'8.1 material parameter lower range'
        write(1,*) lowb(1:nparam)
        write(1,*)'8.2 material parameter upper range'
        write(1,*) upb(1:nparam)
        if (bctype == 1) then
            write(1,*) '9. BC type: four face fictional (comp test only)'
        elseif (bctype == 2) then
            write(1,*) '9. BC type: simple compression-tension'
        elseif (bctype == 3) then
            write(1,*) '9. BC type: three face fictional (comp test only)'
        elseif (bctype == 4) then
            write(1,*) '9. BC type: indentation'
        endif
        write(1,*) '10. Ext Int Force vector size:',group
        write(1,*)'-------------------------------------------------------------------'
        write(1,*)'------------------------end of info session------------------------'
        write(1,*)'-------------------------------------------------------------------'
        write(1,*)''
        write(1,*)''
        write(1,*)'-------------------------------------------------------------------'
        write(1,*)'----------------------------start computing------------------------'
        write(1,*)'-------------------------------------------------------------------'
    end subroutine start_write
    
    subroutine gen0_print_write(in1,d1,d2)
        implicit none
        real(di) d1,d2
        integer in1
        print*,'________________________independent_pop_cutoff_______________________________'
        print*,'start independent population num:',in1
        write(1,*) 'independent population num:',in1
        write(2,*) 'independent population num:',in1
        print*,'best obj(theta) in generation 0>>',d1
        print*,'mean obj(theta) in generation 0 ',d2
        write(1,*) 'best obj(theta) in generation 0>>',d1
        write(1,*) 'mean obj(theta) in generation 0 ',d2
        write(2,*) d1,d2
    end subroutine gen0_print_write
    
    subroutine gen_print_write(in1,d1,d2,theta)
        implicit none
        real(di) d1,d2,theta(nparam)
        integer in1
        write(1,*) '--------------------gen-cutoff-------------------------' 
        print*,'best obj(theta) in generation<<',in1,'is',d1
        print*,'mean obj(theta) in generation',in1,'is',d2
        print*,'best theta in current generation'
        print*,theta
        write(1,*) '-------------------------------------------------------' 
        write(1,*) 'best obj(theta) in generation<<',in1,'is',d1
        write(1,*) 'mean obj(theta) in generation',in1,'is',d2
        write(1,*) 'best theta in current generation'
        write(1,*) theta
        write(2,*) d1,d2
    end subroutine gen_print_write
    
    subroutine after_print_write(d1,theta)
        implicit none
        real(di) d1,theta(nparam)
        print*,'Current independent population ends'
        print*,'best obj(theta) popwise',d1
        print*,'best corresponding theta'
        print*,theta
        write(1,*) 'Current independent population ends'
        write(1,*) 'best obj(theta) popwise',d1
        write(1,*) 'best corresponding theta'
        write(1,*) theta
    end subroutine after_print_write
    
    
end module print_write
    
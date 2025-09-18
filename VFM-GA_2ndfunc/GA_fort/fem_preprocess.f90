!! data input part that reads the structured data and process it into a workable form 
!! you need to edit the number of nodes and elements and number of input frames in the S_1 matlab file to print out a constants.f90 that is compatible with your own data. 
!! the structure of the input preprocessor here is arranged such that you have a compression dominated experiment and tension dominated experiment field input
!! you can clearly have another shear dominated experiment which you just have to write another set of subroutine to take them in and modify the objective function computation part. 
module fem_pre
    use constants
    implicit none
    contains
!   Obtain the displacement and reaction force matrix for compression experiment
    subroutine get_u_rf(rf_mat,disp_mat)
        implicit none
        real(di) rf_sum,rf1,rf2,dispsum,disp1,disp2, &
&       rf_mat(nnode,ncomp*2),disp_mat(ncomp,nnode*2)
        integer index,i,k,pp
        character(len=1) :: n1 ! initial c
        character(len=2) :: n2 ! two digit number
        character(len=1) :: n2_2 ! one digit number
        character(len=1) :: n3 ! zero
        character(len=7) :: n4 ! char that put together to be name
        character(len=4) :: n5 !suffix .rpt
        n1 = 'w'
        write(n3,'(i1)') 0
        n5 = '.rpt'
        do i = 1,ncomp
            if (i<11) then
                write(n2_2,'(i1)') i-1
                n4 = n1//n3//n2_2//n5
            else
                write(n2,'(i2)') i-1
                n4 = n1//n2//n5
            end if 
            open(unit = 100,file=n4,status = 'old',action = 'read')
                do k = 1,nnode
                    read(100,*) pp,rf_sum,rf1,rf2,dispsum,disp1,disp2
                    rf_mat(k,i*2-1:i*2) = (/rf1,rf2/) 
                    disp_mat(i,k*2-1:k*2) = (/disp1,disp2/)
                end do
            close(unit=100)
        end do
    end subroutine get_u_rf 
    !! Handle the case when number of inputs is beyond 100
    subroutine get_u_rf_100(rf_mat,disp_mat)
        implicit none
        real(di) rf_sum,rf1,rf2,dispsum,disp1,disp2, &
&       rf_mat(nnode,nframe*2),disp_mat(nframe,nnode*2)
        integer index,i,k,pp
        character(len=1) :: n1 ! initial c
        character(len=2) :: n2 ! two digit number
        character(len=1) :: n2_2 ! one digit number
        character(len=1) :: n3 ! zero
        character(len=8) :: n4 ! char that put together to be name
        character(len=4) :: n5 !suffix .rpt
        character(len=3) :: n2_3 !three digit number
        n1 = 'w'
        write(n3,'(i1)') 0
        n5 = '.rpt'
        do i = 1,nframe
            if (i<11) then
                write(n2_2,'(i1)') i-1
                n4 = n1//n3//n3//n2_2//n5
            elseif (i<101) then
                write(n2,'(i2)') i-1
                n4 = n1//n3//n2//n5
            else
                write(n2_3,'(i3)') i-1
                n4 = n1//n2_3//n5
            end if 
            open(unit = 100,file=n4,status = 'old',action = 'read')
                do k = 1,nnode
                    read(100,*) pp,rf_sum,rf1,rf2,dispsum,disp1,disp2
                    rf_mat(k,i*2-1:i*2) = (/rf1,rf2/) 
                    disp_mat(i,k*2-1:k*2) = (/disp1,disp2/)
                end do
            close(unit=100)
        end do
    end subroutine get_u_rf_100
    !! process the undeformed mesh for compression experiment
    subroutine get_n_el(coord,connection)
        implicit none
        integer i
        character(len=16) :: n1
        character(len=13) :: n2 
        real(di)  x1,x2,coord(2,nnode),index2
        integer index,el1,el2,el3,el4,connection(neln,nel)
        n1 = 'comp_element.inp'
        n2 = 'comp_node.inp'
        open(unit = 101,file=n1,status = 'old',action = 'read')
            do i=1,nel
                read(101,*) index,el1,el2,el3,el4
                connection(1:4,i) = (/el1,el2,el3,el4/)
            end do
        close(unit = 101)
        open(unit = 102,file=n2,status = 'old',action = 'read')
            do i=1,nnode
                read(102,*) index2,x1,x2
                coord(1:2,i) = (/x1,x2/)
            end do
        close(unit = 102)
    end subroutine get_n_el
    
!   Obtain the displacement and reaction force matrix for tension experiment   
    subroutine get_u_rf2(rf_mat,disp_mat)
        implicit none
        real(di) rf_sum,rf1,rf2,dispsum,disp1,disp2, &
&       rf_mat(nnode2,nten*2),disp_mat(nten,nnode2*2)
        integer index,i,k,pp
        character(len=1) :: n1 ! initial c
        character(len=2) :: n2 ! two digit number
        character(len=1) :: n2_2 ! one digit number
        character(len=1) :: n3 ! zero
        character(len=7) :: n4 ! char that put together to be name
        character(len=4) :: n5 !suffix .rpt
        n1 = 't'
        write(n3,'(i1)') 0
        n5 = '.rpt'
        do i = 1,nten
            if (i<11) then
                write(n2_2,'(i1)') i-1
                n4 = n1//n3//n2_2//n5
            else
                write(n2,'(i2)') i-1
                n4 = n1//n2//n5
            end if 
            open(unit = 100,file=n4,status = 'old',action = 'read')
                do k = 1,nnode2
                    read(100,*) pp,rf_sum,rf1,rf2,dispsum,disp1,disp2
                    rf_mat(k,i*2-1:i*2) = (/rf1,rf2/) 
                    disp_mat(i,k*2-1:k*2) = (/disp1,disp2/)
                end do
            close(unit=100)
        end do
    end subroutine get_u_rf2

    !! process the undeformed mesh for tension experiment
    subroutine get_n_el2(coord,connection)
        implicit none
        integer i
        character(len=15) :: n1
        character(len=12) :: n2 
        real(di)  x1,x2,coord(2,nnode2),index2
        integer index,el1,el2,el3,el4,connection(neln,nel2)
        n1 = 'ten_element.inp'
        n2 = 'ten_node.inp'
        open(unit = 101,file=n1,status = 'old',action = 'read')
            do i=1,nel2
                read(101,*) index,el1,el2,el3,el4
                connection(1:4,i) = (/el1,el2,el3,el4/)
            end do
        close(unit = 101)
        open(unit = 102,file=n2,status = 'old',action = 'read')
            do i=1,nnode2
                read(102,*) index2,x1,x2
                coord(1:2,i) = (/x1,x2/)
            end do
        close(unit = 102)
    end subroutine get_n_el2
end module fem_pre
    

    
!! data input part that reads the structured data and process it into a workable form 
module fem_pre
    use constants
    implicit none
    contains
!   Obtain the displacement and reaction force matrix 
    subroutine get_u_rf(rf_mat,disp_mat)
        implicit none
        real(di) rf_sum,rf1,rf2,dispsum,disp1,disp2, &
&       rf_mat(nnode,nframe*dim),disp_mat(nframe,nnode*dim)
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
        do i = 1,nframe
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
&       rf_mat(nnode,nframe*dim),disp_mat(nframe,nnode*dim)
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
!! process the undeformed mesh, in the 1st functionality case, usually 9nodes 4elements is used 
    subroutine get_n_el(coord,connection)
        implicit none
        integer i
        character(len=16) :: n1
        character(len=13) :: n2 
        real(di)  x1,x2,coord(dim,nnode),index2
        integer index,el1,el2,el3,el4,connection(neln,nel)
        n1 = 'fit9_element.inp'
        n2 = 'fit9_node.inp'
        open(unit = 101,file=n1,status = 'old',action = 'read')
            do i=1,nel
                read(101,*) index,el1,el2,el3,el4
                connection(1:neln,i) = (/el1,el2,el3,el4/)
            end do
        close(unit = 101)
        open(unit = 102,file=n2,status = 'old',action = 'read')
            do i=1,nnode
                read(102,*) index2,x1,x2
                coord(1:dim,i) = (/x1,x2/)
            end do
        close(unit = 102)
    end subroutine get_n_el
end module fem_pre
    

    
!!Basic finite element, see https://solidmechanics.org/Text/Chapter8_1/Chapter8_1.php by Allan Bower
!!Used to compute the external and internal residual force terms 
!!F90 doesn't have OOP features, code written in this way to maintain cleaness
module fembasics
    use constants
    implicit none
    private
    real(di) p57,p77,p2,p6,p555,p888,o2796,o2596,o13,o16
    parameter(p57=0.5773502692d0,p77=0.7745966692d0,p2=0.2d0, &
&             p6=0.6d0,p555=0.555555d0,p888=0.888888d0, &
&             o2796=27.d0/96.d0,o2596=25.d0/96.d0,o13=1.d0/3.d0, &
&             o16=1.d0/6.d0)

    type, public :: fembase
        integer :: dim,neln
    contains
        procedure:: intn => getintn
        procedure:: intp => getintp
        procedure:: intw => getintw
        procedure:: shape => getshape
        procedure:: shapedev => getshapedev
    end type fembase
    contains
        function getintn(this) result(intn)
            class(fembase), intent(in):: this
            integer intn
            if (this%dim == 1) then
                intn = this%neln
            else if (this%dim == 2) then
                if (this%neln == 3) then
                    intn = 1
                else if (this%neln == 6) then
                    intn = 3
                else if (this%neln == 4) then
                    intn = 4 
                else if (this%neln == 8) then 
                    intn = 9
                end if
            end if 
        endfunction getintn
        
        function getintp(this) result(xi)
            class(fembase), intent(in):: this
            integer npt,dim
            real(di) xi(this%dim,this%intn())
            xi = 0.d0
            npt = this%intn()
            if (this%dim == 1) then
                if (npt == 1) then
                    xi(1,1) = 0.d0
                else if (npt == 2) then
                    xi(1,1) = -p57
                    xi(1,2) = -xi(1,1)
                else if (npt == 3) then
                    xi(1,1) = -p77
                    xi(1,2) = 0.d0
                    xi(1,3) = -xi(1,1)
                end if 
             else if (this%dim == 2) then
                if (this%neln == 3 .or. this%neln == 6) then
                    if (npt == 1) then
                        xi(1,1) = o13
                        xi(2,1) = o13
                    else if (npt == 3) then
                        xi(1, 1) = p6
                        xi(2, 1) = p2
                        xi(1, 2) = p2
                        xi(2, 2) = p6
                        xi(1, 3) = p2
                        xi(2, 3) = p2
                    else if (npt == 4) then
                        xi(1, 1) = o13
                        xi(2, 1) = o13
                        xi(1, 2) = p6
                        xi(2, 2) = p2
                        xi(1, 3) = p2
                        xi(2, 3) = p6
                        xi(1, 4) = p2
                        xi(2, 4) = p2
                    end if 
                else if (this%neln == 4 .or. this%neln == 8) then
                    if (npt == 1) then
                        xi(1,1) = 0.d0
                        xi(2,1) = 0.d0
                    else if (npt == 4) then
                        xi(1, 1) = -p57
                        xi(2, 1) = xi(1, 1)
                        xi(1, 2) = -xi(1, 1)
                        xi(2, 2) = xi(1, 1)
                        xi(1, 3) = -xi(1, 1)
                        xi(2, 3) = -xi(1, 1)
                        xi(1, 4) = xi(1, 1)
                        xi(2, 4) = -xi(1, 1)
                    else if (npt == 9) then
                        xi(1, 1) = -p77
                        xi(2, 1) = xi(1, 1)
                        xi(1, 2) = 0.d0
                        xi(2, 2) = xi(1, 1)
                        xi(1, 3) = -xi(1, 1)
                        xi(2, 3) = xi(1, 1)
                        xi(1, 4) = xi(1, 1)
                        xi(2, 4) = 0.d0
                        xi(1, 5) = 0.d0
                        xi(2, 5) = 0.d0
                        xi(1, 6) = -xi(1, 1)
                        xi(2, 6) = 0.d0
                        xi(1, 7) = xi(1, 1)
                        xi(2, 7) = -xi(1, 1)
                        xi(1, 8) = 0.d0
                        xi(2, 8) = -xi(1, 1)
                        xi(1, 9) = -xi(1, 1)
                        xi(2, 9) = -xi(1, 1)
                    end if 
                end if 
            end if 
        endfunction getintp
            
        function getshape(this,xi) result(n)
            class(fembase), intent(in):: this
            real(di), intent(in) :: xi(2)
            real(di) n(this%neln),xi3
            if (this%dim == 1) then
                if (this%neln == 2) then
                    n(1) = 0.5d0*(1.d0+xi(1))
                    n(2) = 0.5d0*(1.d0-xi(1))
                else if (this%neln == 3) then
                    n(1) = -0.5d0*xi(1)*(1.d0-xi(1))
                    n(2) =  0.5d0*xi(1)*(1.d0+xi(1))
                    n(3) = (1.d0-xi(1))*(1.d0+xi(1))
                end if 
            else if (this%dim == 2) then
                if (this%neln == 3) then
                    n(1) = xi(1)
                    n(2) = xi(2)
                    n(3) = 1.d0-xi(1)-xi(2) 
                else if (this%neln == 6) then
                    xi3 = 1.d0-xi(1)-xi(2)
                    n(1) = (2.d0*xi(1)-1.d0)*xi(1)
                    n(2) = (2.d0*xi(2)-1.d0)*xi(2)
                    n(3) = (2.d0*xi3-1.d0)*xi3
                    n(4) = 4.d0*xi(1)*xi(2)
                    n(5) = 4.d0*xi(2)*xi3
                    n(6) = 4.d0*xi3*xi(1)
                else if (this%neln == 4) then
                    n(1) = 0.25d0*(1.d0-xi(1))*(1.d0-xi(2))
                    n(2) = 0.25d0*(1.d0+xi(1))*(1.d0-xi(2))
                    n(3) = 0.25d0*(1.d0+xi(1))*(1.d0+xi(2))
                    n(4) = 0.25d0*(1.d0-xi(1))*(1.d0+xi(2))
                else if (this%neln == 8) then
                    n(1) = -0.25d0*(1.d0-xi(1))*(1.d0-xi(2)) &
&                   *(1.d0+xi(1)+xi(2))
                    n(2) = 0.25d0*(1.d0+xi(1))*(1.d0-xi(2)) & 
&                   *(xi(1)-xi(2)-1.d0)
                    n(3) = 0.25d0*(1.d0+xi(1))*(1.d0+xi(2)) &
&                   *(xi(1)+xi(2)-1.d0)
                    n(4) = 0.25d0*(1.d0-xi(1))*(1.d0+xi(2)) &
&                   *(xi(2)-xi(1)-1.d0)
                    n(5) = 0.5d0*(1.d0-xi(1)*xi(1))*(1.d0-xi(2))
                    n(6) = 0.5d0*(1.d0+xi(1))*(1.d0-xi(2)*xi(2))
                    n(7) = 0.5d0*(1.d0-xi(1)*xi(1))*(1.d0+xi(2))
                    n(8) = 0.5d0*(1.d0-xi(1))*(1.d0-xi(2)*xi(2))
                end if 
            end if
        endfunction getshape
        
        function getshapedev(this,xi) result(dNdxi)
            class(fembase), intent(in):: this  
            real(di), intent(in) :: xi(2)
            real(di) dNdxi(this%neln,this%dim),xi3
            if (this%dim == 1) then
                if (this%neln == 2) then
                    dNdxi(1,1) = 0.5d0
                    dNdxi(2,1) = -0.5d0
                else if (this%neln == 3) then
                    dNdxi(1,1) = -0.5d0+xi(1)
                    dNdxi(2,1) =  0.5d0+xi(1)
                    dNdxi(3,1) = -2.d0*xi(1)
                end if
            else if (this%dim == 2) then
                if (this%neln == 3) then
                    dNdxi(1,1) = 1.d0
                    dNdxi(2,2) = 1.d0
                    dNdxi(3,1) = -1.d0
                    dNdxi(3,2) = -1.d0
                else if (this%neln == 6) then
                    xi3 = 1.-xi(1)-xi(2)
                    dNdxi(1,1) = 4.d0*xi(1)-1.d0
                    dNdxi(2,2) = 4.d0*xi(2)-1.d0
                    dNdxi(3,1) = -(4.d0*xi3-1.d0)
                    dNdxi(3,2) = -(4.d0*xi3-1.d0)
                    dNdxi(4,1) = 4.d0*xi(2)
                    dNdxi(4,2) = 4.d0*xi(1)
                    dNdxi(5,1) = -4.d0*xi(2)
                    dNdxi(5,2) = -4.d0*xi(1)
                    dNdxi(6,1) = 4.d0*xi3 - 4.d0*xi(1)
                    dNdxi(6,2) = 4.d0*xi3 - 4.d0*xi(2)
                else if (this%neln == 4) then
                    dNdxi(1,1) = -0.25d0*(1.d0-xi(2))
                    dNdxi(1,2) = -0.25d0*(1.d0-xi(1))
                    dNdxi(2,1) = 0.25d0*(1.d0-xi(2))
                    dNdxi(2,2) = -0.25d0*(1.d0+xi(1))
                    dNdxi(3,1) = 0.25d0*(1.d0+xi(2))
                    dNdxi(3,2) = 0.25d0*(1.d0+xi(1))
                    dNdxi(4,1) = -0.25d0*(1.d0+xi(2))
                    dNdxi(4,2) = 0.25d0*(1.d0-xi(1))
                else if (this%neln == 8) then
                    dNdxi(1,1) = 0.25d0*(1.d0-xi(2)) &
&                   *(2.d0*xi(1)+xi(2))
                    dNdxi(1,2) = 0.25d0*(1.d0-xi(1)) &
&                   *(xi(1)+2.d0*xi(2))
                    dNdxi(2,1) = 0.25d0*(1.d0-xi(2)) &
&                   *(2.d0*xi(1)-xi(2))
                    dNdxi(2,2) = 0.25d0*(1.d0+xi(1)) &
&                   *(2.d0*xi(2)-xi(1))
                    dNdxi(3,1) = 0.25d0*(1.d0+xi(2)) &
&                   *(2.d0*xi(1)+xi(2))
                    dNdxi(3,2) = 0.25d0*(1.d0+xi(1)) &
&                   *(2.d0*xi(2)+xi(1))
                    dNdxi(4,1) = 0.25d0*(1.d0+xi(2)) &
&                   *(2.d0*xi(1)-xi(2))
                    dNdxi(4,2) = 0.25d0*(1.d0-xi(1)) &
&                   *(2.d0*xi(2)-xi(1))
                    dNdxi(5,1) = -xi(1)*(1.d0-xi(2))
                    dNdxi(5,2) = -0.5d0*(1.d0-xi(1)*xi(1))
                    dNdxi(6,1) = 0.5d0*(1.d0-xi(2)*xi(2))
                    dNdxi(6,2) = -(1.d0+xi(1))*xi(2)
                    dNdxi(7,1) = -xi(1)*(1.d0+xi(2))
                    dNdxi(7,2) = 0.5d0*(1.d0-xi(1)*xi(1))
                    dNdxi(8,1) = -0.5d0*(1.d0-xi(2)*xi(2))
                    dNdxi(8,2) = -(1.d0-xi(1))*xi(2)
                end if
            end if 
        endfunction getshapedev
        
        function getintw(this) result(w)
            implicit none
            class(fembase), intent(in):: this
            integer i,j,index,npt
            real(di) wd(3),w(this%intn()),xi3
            npt = this%intn()
            if (this%dim ==1) then
                if (npt == 1) then
                    w(1) = 2.d0
                else if (npt == 2) then
                    w = (/1.d0,1.d0/)
                else if (npt == 3) then
                    w = (/p555,p888,p555/)
                end if
            else if (this%dim == 2) then
                if (this%neln == 3 .or. this%neln == 6) then
                    if (npt == 1) then
                        w(1) = 0.5d0
                    else if (npt == 3) then
                        w = (/o16,o16,o16/)
                    else if (npt == 4) then
                        w = (/-o2796,o2596,o2596,o2596/)
                    end if
                else if (this%neln == 4 .or. this%neln == 8) then
                    if (npt == 1) then
                        w(1) = 4.d0
                    else if (npt == 4) then
                        w = (/1.d0,1.d0,1.d0,1.d0/)
                    else if (npt == 9) then
                        wd = (/p555,p888,p555/)
                        do i = 1,3
                            do j = 1,3
                                index = 3.d0*(j-1.d0)+i
                                w(index) = wd(i)+wd(j)
                            end do
                        end do
                    end if 
                end if        
            end if 
        endfunction getintw
                                
end module fembasics
    

!--------------------Newton_Interpolation------------------
!--------------------2020.4.3------------------------------
!--------------------sy------------------------------------
      module interpolation
        contains
            subroutine newtdd(x,y,c)
            implicit none
            real(kind=8),intent( in ) :: x(:),y(:)
            real(kind=8),intent( inout ) :: c(:,:)
            real(kind=8),allocatable :: v(:,:)
            integer :: i,j,n

            n=size(x)
            allocate( v(n,n) )
            v=0.d0
            do i=1,size(x)
                v(i,1)=y(i)
            end do

            do j=2,size(x)
                do i=1,size(x)+1-j
                    v(i,j)=(v(i+1,j-1)-v(i,j-1))/(x(i+j-1)-x(i))
                end do
            end do

            c(1,:)=v(1,:)
            end subroutine newtdd

            subroutine cal_interpolation(x,xx,c,res)
                implicit none
                 real(kind=8),intent( in ) :: x(:),xx(:),c(:,:)
                 real(kind=8),intent( inout ) :: res(:,:)
                 real(kind=8),allocatable :: d(:,:)
                 integer :: i,j,m

                 m=size(x)
                 allocate(d(m,m))

                 d(1,:)=1.d0
                 do j=1,m
                     do i=2,m
                         d(i,j)=d(i-1,j)*(xx(j)-x(i-1))
                     end do
                 end do
                 res=matmul(c,d)
                 end subroutine cal_interpolation
                 end module interpolation

                 program NewtonInterpolation
                     use interpolation
                     implicit none
                     integer,parameter ::n=20
                     real(kind=8) :: a=-1.d0,b=1.d0
                     real(kind=8) :: x(n),xx(n),y(n),c(1,n),res(1,n)
                     integer :: i

                    do i=1,n
                        x(i)=a+real(i,8)*(b-a)/real(n,8)
                        xx(i)=a/2.d0+real(i,8)*(b-a)/2.d0/real(n,8)
                        y(i)=exp(x(i))
                    end do
                    call newtdd(x,y,c)
                    call cal_interpolation(x,xx,c,res)
                    open(101,file='ans.dat')
                    do i=1,n
                        write(101,'(4(2x,g0))')x(i),y(i),xx(i),res(1,i)
                    end do
                    close(101)
                    end program NewtonInterpolation









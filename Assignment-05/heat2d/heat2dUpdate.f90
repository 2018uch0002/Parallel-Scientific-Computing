!
! Heat Equation Fortran update
!
! Note: this function is written as single precision but can be converted to
! double precision at compile time
!
subroutine heat2dUpdate( n1a,n1b,n2a,n2b, nd1a,nd1b, nd2a,nd2b, un,u, rx,ry, numThreads )

implicit none
integer n1a,n1b,n2a,n2b,nd1a,nd1b, nd2a,nd2b, numThreads
real rx,ry
real un(nd1a:nd1b,nd2a:nd2b)
real u(nd1a:nd1b,nd2a:nd2b)

! local variables
integer i1,i2

!$omp parallel do default(shared) private(i1,i2) num_threads(numThreads)
do i2=n2a,n2b
do i1=n1a,n1b

un(i1,i2) = u(i1,i2) + rx*( u(i1+1,i2) -2.*u(i1,i2) + u(i1-1,i2) ) + &
ry*( u(i1,i2+1) -2.*u(i1,i2) + u(i1,i2-1) )
end do
end do
!$omp end parallel do

return
end

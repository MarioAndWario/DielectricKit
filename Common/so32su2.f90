! ==============================================
! Find one of the U(2) matrix for a O(3) rotation: proper or improper
! Author : Meng Wu
! Date : 20151019
! Version : Beta
! ==============================================

#include "f_defs.h"

module so32su2_m
  use global_m

  implicit none

  private

  public :: so32su2

contains

  subroutine so32su2(Rmatrix, umatrix)
    real(DP), dimension(3,3) :: Rmatrix ! cartesian coordinates
    ! LatVec(1:3,1:3) is the lattice vectors in cartesian coordinates
    ! LatVec(1:3,i) is the i-th lattice vector
    ! umatrix(1:2,1:2) is the U(2) matrix corresponding to the Rmatrix O(3) rotation
    complex(DPC), dimension(2,2), intent(out) :: umatrix
    ! eps10 is a small number
    real(DP), parameter :: eps10 = 1.0D-10
    ! imaginary unit, (0+1j)
    complex(DPC), parameter :: ii=(0,1)
    real(DP) :: alpha, beta, gamma
    real(DP), parameter :: Pi = 3.14159265358979323846
    real(DP) :: sinbeta, Rdet
    ! for proper rotation: ifprop = .true.
    ! for improper rotation: ifprop = .false.
    logical :: ifprop = .true.

    ! ------------------------------------
    ! Calculate determinant of Rmatrix(i,j)
    Rdet = Rmatrix(1,1)*(Rmatrix(2,2)*Rmatrix(3,3) - Rmatrix(3,2)*Rmatrix(2,3)) &
         + Rmatrix(1,2)*(Rmatrix(3,1)*Rmatrix(2,3) - Rmatrix(2,1)*Rmatrix(3,3))  &
         + Rmatrix(1,3)*(Rmatrix(2,1)*Rmatrix(3,2) - Rmatrix(3,1)*Rmatrix(2,2))

    if (abs(Rdet-1.0) < eps10) then
       ! write(*,*) 'Proper rotation'
       ifprop = .true.

    elseif (abs(Rdet+1.0) < eps10) then
       ! write(*,*) 'Improper rotation'
       ifprop = .false.
       Rmatrix = - Rmatrix

    else
       write(*,*) 'Error in so32su2()!!!'
       call EXIT(1)
    endif

    if (abs(Rmatrix(3,3)-1.0) < eps10) then
       beta = 0
       gamma = 0
       alpha = atan2(Rmatrix(2,1),Rmatrix(1,1))
    elseif (abs(Rmatrix(3,3)+1.0) < eps10) THEN
       beta = Pi;
       gamma = 0;
       alpha = - atan2(Rmatrix(1,2),Rmatrix(2,2))
    else
       sinbeta=sqrt(1-Rmatrix(3,3)*Rmatrix(3,3));
       beta=atan2(sinbeta,Rmatrix(3,3));
       alpha=atan2(Rmatrix(2,3)/sinbeta,Rmatrix(1,3)/sinbeta);
       gamma=atan2(Rmatrix(3,2)/sinbeta,-Rmatrix(3,1)/sinbeta);
    endif

    umatrix(1,1) = exp(-ii*(alpha+gamma)/2.0)*cos(beta/2.0)
    umatrix(1,2) = -exp(-ii*(alpha-gamma)/2.0)*sin(beta/2.0)
    umatrix(2,1) = exp(ii*(alpha-gamma)/2.0)*sin(beta/2.0)
    umatrix(2,2) = exp(ii*(alpha+gamma)/2.0)*cos(beta/2.0)

    !    write(*,*) '====== Umatrix ======'
    !     do i=1,2
    !        write(*,'(F15.10,"+",F15.10,"i",5X,F15.10,"+",F15.10,"i")') umatrix(i,:)
    !     enddo


    ! if we have improper rotation, we need to multiply the umatrix by {{i,0},{i,0}}
    if (ifprop .eqv. .false.) then
       ! we have I [spinor] = \pm spinor
       ! take +1 as prefactor will not change the value of matrix elements
       ! umatrix = ii*umatrix
       umatrix = umatrix
       ! write(*,*) '====== i*Umatrix ======'
       ! do i=1,2
       !    write(*,'(F15.10,"+",F15.10,"i",5X,F15.10,"+",F15.10,"i")') umatrix(i,:)
       ! enddo
    endif

  end subroutine so32su2

  ! subroutine inverse(a_,c,n)
  !   !============================================================
  !   ! Inverse matrix
  !   ! Method: Based on Doolittle LU factorization for Ax=b
  !   ! Alex G. December 2009
  !   !-----------------------------------------------------------
  !   ! input ...
  !   ! a(n,n) - array of coefficients for matrix A
  !   ! n      - dimension
  !   ! output ...
  !   ! c(n,n) - inverse matrix of A
  !   ! comments ...
  !   ! the original matrix a(n,n) will be destroyed
  !   ! during the calculation
  !   !===========================================================
  !   implicit none
  !   integer :: n
  !   real(DP), intent(out) :: c(n,n)
  !   real(DP) :: a(n,n)
  !   real(DP), intent(in) :: a_(n,n)
  !   real(DP) :: L(n,n), U(n,n), b(n), d(n), x(n)
  !   real(DP) :: coeff
  !   integer :: i, j, k

  !   ! step 0: initialization for matrices L and U and b
  !   ! Fortran 90/95 aloows such operations on matrices
  !   L=0.0
  !   U=0.0
  !   b=0.0

  !   do i=1,n
  !      do j=1,n
  !         a(i,j) = a_(i,j)
  !      enddo
  !   enddo

  !   ! step 1: forward elimination
  !   do k=1, n-1
  !      do i=k+1,n
  !         coeff=a(i,k)/a(k,k)
  !         L(i,k) = coeff
  !         do j=k+1,n
  !            a(i,j) = a(i,j)-coeff*a(k,j)
  !         end do
  !      end do
  !   end do

  !   ! Step 2: prepare L and U matrices
  !   ! L matrix is a matrix of the elimination coefficient
  !   ! + the diagonal elements are 1.0
  !   do i=1,n
  !      L(i,i) = 1.0
  !   end do
  !   ! U matrix is the upper triangular part of A
  !   do j=1,n
  !      do i=1,j
  !         U(i,j) = a(i,j)
  !      end do
  !   end do

  !   ! Step 3: compute columns of the inverse matrix C
  !   do k=1,n
  !      b(k)=1.0
  !      d(1) = b(1)
  !      ! Step 3a: Solve Ld=b using the forward substitution
  !      do i=2,n
  !         d(i)=b(i)
  !         do j=1,i-1
  !            d(i) = d(i) - L(i,j)*d(j)
  !         end do
  !      end do
  !      ! Step 3b: Solve Ux=d using the back substitution
  !      x(n)=d(n)/U(n,n)
  !      do i = n-1,1,-1
  !         x(i) = d(i)
  !         do j=n,i+1,-1
  !            x(i)=x(i)-U(i,j)*x(j)
  !         end do
  !         x(i) = x(i)/u(i,i)
  !      end do
  !      ! Step 3c: fill the solutions x(n) into column k of C
  !      do i=1,n
  !         c(i,k) = x(i)
  !      end do
  !      b(k)=0.0
  !   end do
  ! end subroutine inverse

end module so32su2_m

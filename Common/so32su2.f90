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
       ! 'Proper rotation'
       ifprop = .true.

    elseif (abs(Rdet+1.0) < eps10) then
       ! 'Improper rotation'
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

    ! if we have improper rotation, we need to multiply the umatrix by {{i,0},{i,0}}
    if (ifprop .eqv. .false.) then
       ! we have I [spinor] = \pm spinor
       ! take +1 as prefactor will not change the value of matrix elements
       ! umatrix = ii*umatrix
       umatrix = umatrix
    endif

  end subroutine so32su2

end module so32su2_m

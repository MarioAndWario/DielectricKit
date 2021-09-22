#include "f_defs.h"

!===============================================================================
!
! Modules:
!
! input_utils_m   Originally By DAS     Last modified by Meng Wu
!
!===============================================================================

module input_utils_m
  use global_m
  use blas_m
  implicit none
  private
  public :: gvec_index, kinetic_energies
contains

  !---------------------------------------------------------------------------------------------------
  !! Compute index_vec indices relating G-vectors in reduced coordinates to positions in the FFT grid
  subroutine gvec_index(gvec)
    type(gspace), intent(inout) :: gvec
    integer :: ig, iadd

    PUSH_SUB(gvec_index)

    gvec%nFFTgridpts = product(gvec%FFTgrid(1:3))
    SAFE_ALLOCATE(gvec%index_vec, (gvec%nFFTgridpts))
    gvec%index_vec(:) = 0
    do ig = 1, gvec%ng
       !! if a mean-field code does not use the appropriate convention, this could happen.
       if (any(2 * gvec%components(1:3, ig) >= gvec%FFTgrid(1:3) .or. 2 * gvec%components(1:3, ig) < -gvec%FFTgrid(1:3))) then
          call die("gvectors must be in the interval [-FFTgrid/2, FFTgrid/2)")
       endif
       iadd = ((gvec%components(1,ig)+gvec%FFTgrid(1)/2)*gvec%FFTgrid(2)+gvec%components(2,ig)+ gvec%FFTgrid(2)/2)*gvec%FFTgrid(3)+gvec%components(3,ig)+gvec%FFTgrid(3)/2+1
       gvec%index_vec(iadd) = ig
    enddo

    POP_SUB(gvec_index)
    return
  end subroutine gvec_index

  !-----------------------------------------------------------------
  !> This routine calculates the kinetic energies |G+q|^2 or |G|^2 for all
  !! the G-vectors in gvec%components using the reciprocal metric bdot:
  !!   ekin(ig) = \sum_{m,n} G(m,ig) B(m,n) G(n,ig)
  !! We perform the sum by first performing the Cholesky decomposition of B,
  !! B := U^T U. Then, we write V := U G and write ekin = V^T V
  !! Using Cholesky decomposition has the same flop count as using dgemms, but
  !! it`s easier for the compiler to vectorize.
  !!
  !! \param gvec gspace structure that contains all the gvectors gvec%components
  !! \param bdot reciprocal metric
  !! \param ekin array holding the output kinetic energies
  !! \param qvec use this to compute |q+G|^2 instead of |G|^2
  subroutine kinetic_energies(gvec, bdot, ekin, qvec)
    type(gspace), intent(in) :: gvec
    real(DP), intent(in) :: bdot(3, 3)
    real(DP), intent(out) :: ekin(:) !< (gvec%ng)
    real(DP), optional, intent(in) :: qvec(3)

    integer :: ig, info
    real(DP) :: qkv(3,gvec%ng), vmid(3), U(3,3) ! FHJ: stack allocation is faster!

    PUSH_SUB(kinetic_energies)

    if (present(qvec)) then
       do ig = 1,gvec%ng
          qkv(1:3,ig) = qvec(1:3) + DBLE(gvec%components(1:3,ig))
       enddo
    else
       qkv(1:3,1:gvec%ng) = DBLE(gvec%components(1:3,1:gvec%ng))
    endif

    U(1:3, 1:3) = bdot(1:3, 1:3)
    ! Cholesky decomposition of the metric: bdot = U^T U
    call dpotrf('U', 3, U, 3, info)
    do ig = 1, gvec%ng
       vmid(1) = U(1,1)*qkv(1,ig) + U(1,2)*qkv(2,ig) + U(1,3)*qkv(3,ig)
       vmid(2) =                    U(2,2)*qkv(2,ig) + U(2,3)*qkv(3,ig)
       vmid(3) =                                       U(3,3)*qkv(3,ig)
       ekin(ig) = vmid(1)**2 + vmid(2)**2 + vmid(3)**2
    enddo

    POP_SUB(kinetic_energies)
    return
  end subroutine kinetic_energies

end module input_utils_m

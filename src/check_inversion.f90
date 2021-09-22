!===============================================================================
!
! Module:
!
! (1) check_inversion_m   Originally By DAS            Last Modified 10/14/2010
!
!   Check whether our choice of real/complex version is appropriate given the
!   presence or absence of inversion symmetry about the origin, and a guess
!   about time-reversal symmetry depending on the number of spin-components.
!
!===============================================================================

#include "f_defs.h"

module check_inversion_m

  use global_m
  implicit none

  private

  public :: check_inversion, check_inversion_type

contains  

subroutine check_inversion(iflavor, ntran, mtrx, nspin, warn, real_need_inv, tnp)
  integer, intent(in) :: iflavor
  integer, intent(in) :: ntran
  integer, intent(in) :: mtrx(3, 3, 48) !< symmetry operations matrices
  integer, intent(in) :: nspin
  logical, intent(in) :: warn !< set to false to suppress warnings, for converters
  logical, intent(in) :: real_need_inv !< use for generating routines to block real without inversion
     !! this is not always true so that it is possible to run real without using symmetries
  real(DP), optional, intent(in) :: tnp(3, 48) !< fractional translations.
     !! optional only to avoid changing external interface for library.

  integer :: invflag, isym, ii, jj, itest, real_or_complex
  logical :: origin_inv
  character(len=7) :: sflavor

  PUSH_SUB(check_inversion)

  if(iflavor .eq. 0) then
#ifdef CPLX
    real_or_complex = 2
#else
    real_or_complex = 1
#endif
  elseif(iflavor .eq. 1 .or. iflavor .eq. 2) then
    real_or_complex = iflavor
  else
    write(sflavor, '(i7)') iflavor
    call die("Illegal value iflavor = " // TRUNC(sflavor) // " passed to check_inversion: must be 0,1,2.", &
      only_root_writes=.true.)
  endif

  invflag = 0
  origin_inv = .false.
  do isym = 1, ntran
    itest = 0
    do ii = 1, 3
      do jj = 1, 3
        if(ii .eq. jj) then
          itest = itest + (mtrx(ii, jj, isym) + 1)**2
        else
          itest = itest + mtrx(ii, jj, isym)**2
        endif
      enddo
    enddo
    if(itest .eq. 0) then
      invflag = invflag + 1
      if(present(tnp)) then
        if(sum(abs(tnp(1:3, isym))) < TOL_Small) origin_inv = .true.
      else
        origin_inv = .true.
      endif
    endif
  enddo
  if(invflag > 0 .and. .not. origin_inv .and. peinf%inode==0) then
    write(0, '(a)') "WARNING: Inversion symmetry is present only with a fractional translation."
    write(0, '(a)') "Apply the translation so inversion is about the origin, to be able to use the real version."
  endif
  if(invflag .gt. 1 .and. peinf%inode==0) &
    write(0, '(a)') "WARNING: More than one inversion symmetry operation is present."

  if(invflag > 0 .and. .not. present(tnp) .and. peinf%inode==0) then
    write(0, '(a)') "WARNING: check_inversion did not receive fractional translations."
    write(0, '(a)') "Cannot confirm that inversion symmetry is about the origin for use of real version."
  endif

  if(real_or_complex .eq. 2) then
    if(origin_inv .and. warn .and. nspin == 1) then
      if(peinf%inode .eq. 0) &
        write(0, '(a)') "WARNING: Inversion symmetry about the origin is present. The real version would be faster."
    endif
  else
    if(.not. origin_inv) then
      if(real_need_inv) then
        call die("The real version cannot be used without inversion symmetry about the origin.", only_root_writes = .true.)
      endif
      if(peinf%inode .eq. 0) then
        write(0, '(a)') "WARNING: Inversion symmetry about the origin is absent in symmetries used to reduce k-grid."
        write(0, '(a)') "Be sure inversion about the origin is still a spatial symmetry, or you must use complex version instead."
      endif
    endif
    if(nspin > 1) then
      call die("Real version may only be used for spin-unpolarized calculations.", only_root_writes = .true.)
    endif
  endif

  POP_SUB(check_inversion)
  return
end subroutine check_inversion

!=========================================================================
!> wrapper routine that uses typedefs types
subroutine check_inversion_type(iflavor, syms, nspin, warn, real_need_inv)
  integer, intent(in) :: iflavor
  type (symmetry), intent(in) :: syms
  integer, intent(in) :: nspin
  logical, intent(in) :: warn
  logical, intent(in) :: real_need_inv !< use for generating routines to block real without inversion

  PUSH_SUB(check_inversion_type)

  call check_inversion(iflavor, syms%ntran, syms%mtrx, nspin, warn, real_need_inv, tnp = syms%tnp)

  POP_SUB(check_inversion_type)
  return
end subroutine check_inversion_type

end module check_inversion_m

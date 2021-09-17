!===============================================================================
!
! Modules:
!
! inread_common_m   Originally By FHJ
!
!   A first attempt to unify the inread routines. Right now, this module
!   only implements consistency checks and warning/error messages.
!
!===============================================================================

#include "f_defs.h"

module inread_common_m

  use global_m
  implicit none

  private

  public ::                  &
    check_bounds_nkq,        &
    check_bounds_nbands,     &
    check_consistency_nbands

contains


  !> FHJ: Makes sure the number of {q,k}-points is less than MAX_KPTS.
  !! Call this subroutine just after you read the keyword `number_{k,q}-points`.
  subroutine check_bounds_nkq(nkq, k_or_q, keyword)
    !> Number of {k,q}-points expected to be read (eg, pol%nq)
    integer, intent(in) :: nkq
    character(len=1), intent(in) :: k_or_q !< Either "k" or "q"
    character(len=*), intent(in) :: keyword !< a keyword, such as `number_qpoints`

    PUSH_SUB(check_bounds_nkq)

    if (nkq>MAX_KPTS) then
      write(0,*)
      write(0,'(/,a)') 'ERROR: The number of '//&
        k_or_q//'-points specified in the keyword `'//keyword//'` is '
      write(0,'(a,i0)') ' larger than the maximum: MAX_KPTS=',MAX_KPTS
      write(0,'(a,/)') ' Either use less '//&
        k_or_q//'-points or increase MAX_KPTS in Common/nrtype.f90'
      write(0,*)
      call die('Too many '//k_or_q//'-points. Increase MAX_KPTS in Common/nrtypes.f90.')
    endif

    POP_SUB(check_bounds_nkq)

  end subroutine check_bounds_nkq

  !> FHJ: Makes sure nb<MAX_BANDS
  !! Call this subroutine just after you read the keyword `number_bands`.
  subroutine check_bounds_nbands(nb, keyword)
    integer, intent(in) :: nb
    character(len=*), intent(in) :: keyword !< a keyword, such as `number_qpoints`

    PUSH_SUB(check_bounds_nbands)

    if(nb>MAX_BANDS) then
      write(0,*)
      write(0,'(a)') 'ERROR: The number of bands specified in the keyword `'//keyword//'` is larger '
      write(0,'(a,i0)') ' than the maximum: MAX_BANDS=',MAX_BANDS
      write(0,'(a)') ' Either use less bands or increase MAX_BANDS in Common/nrtype.f90'
      write(0,*)
      call die("Too many bands. Increase MAX_BANDS in Common/nrtypes.f90.")
    endif

    POP_SUB(check_bounds_nbands)

  end subroutine check_bounds_nbands


  !> FHJ: Makes sure nb<MAX_BANDS and that nb!=0 if the keyword is required.
  !! Call this after you parse the whole input file.
  subroutine check_consistency_nbands(nb, is_required)
    integer, intent(in) :: nb !< number of bands
    logical, intent(in) :: is_required !< dies if nb==0

    PUSH_SUB(check_consistency_nbands)
    
    if (peinf%inode>0) then
      POP_SUB(check_consistency_nbands)
      return
    endif

    if(is_required .and. nb<1) then
      call die("The keyword `number_bands` could not be found.", only_root_writes = .true.)
    endif

    call check_bounds_nbands(nb, 'number_bands')

    POP_SUB(check_consistency_nbands)

  end subroutine check_consistency_nbands

end module inread_common_m

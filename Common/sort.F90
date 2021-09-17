#include "f_defs.h"

!>==========================================================================
!!
!! Module sort_m:
!!
!! (1) gcutoff
!!
!!     Given G-vectors sorted by kinetic energy and an energy cutoff,
!!     find the corresponding G-vector cutoff.
!!
!! (2,3) sortrx, sortix
!!
!!     Sorts an array by the quicksort method. real(DP) and integer versions.
!!     See included sort_inc.f90.
!!
!! (4) make_identity_symmetry_first
!!
!!     The identity must always be the first symmetry, as assumed in various places
!!     in the code. We enforce this by swapping it with op #1 if it is not first.
!!
!! (5) sort_symmetries
!!
!!     Bring symmetries into a standardized order. We are not currently using this.
!!
!!==========================================================================

module sort_m
  use global_m

  implicit none

  private

  public ::   &
       gcutoff,  &
       sortrx, &
       sortix, &
       make_identity_symmetry_first, &
       sort_symmetries

  interface sortrx
     module procedure sortrx_gvec, sortrx_no_gvec
  end interface sortrx
  interface sortix
     module procedure sortix_gvec, sortix_no_gvec
  end interface sortix

contains

  !> Given G-vectors sorted by kinetic energy and an energy cutoff, find the corresponding G-vector cutoff
  !! such that all(ekin(isrtrq(ig)) <= ecutoff) for ig <= gcutoff.
  integer function gcutoff(ng, ekin, isrtrq, ecutoff)
    integer, intent(in) :: ng !< number of G-vectors
    real(DP), intent(in) :: ekin(:) !< (ng) kinetic energies, should be sorted already
    integer, intent(in) :: isrtrq(:) !< (ng) this is the index array returned by sorting ekin
    real(DP), intent(in) :: ecutoff  !< energy cutoff, in same units as ekin (Ry typically)

    integer :: gup, gdn, gmid, ig

    PUSH_SUB(gcutoff)

    ! perhaps all G-vectors fall within the cutoff
    if(ekin(isrtrq(ng)) < ecutoff) then
       gcutoff = ng
       POP_SUB(gcutoff)
       return
    endif

    ! otherwise, use bisection
    gup = ng
    gdn = 1

    do ig = 1, ng
       gmid = (gup + gdn) / 2
       if(gmid == gdn) exit
       if(ekin(isrtrq(gmid)) > ecutoff) then
          gup = gmid
       else
          gdn = gmid
       endif
    enddo
    gcutoff = gdn

    POP_SUB(gcutoff)
    return
  end function gcutoff

  !=====================================================================
  !> The identity must always be the first symmetry, as assumed in various places in the code.
  subroutine make_identity_symmetry_first(nsyms, mtrx, tnp)
    integer,  intent(in)    :: nsyms
    integer,  intent(inout) :: mtrx(3, 3, 48)
    real(DP), intent(inout) :: tnp(3, 48)

    integer :: isym, mtrx_temp(3, 3), identity(3,3)
    real(DP) :: tnp_temp(3)
    logical :: found

    PUSH_SUB(make_identity_symmetry_first)

    identity = reshape((/1, 0, 0, 0, 1, 0, 0, 0, 1/), shape(identity))
    found = all(mtrx(1:3, 1:3, 1) == identity(1:3, 1:3))

    do isym = 2, nsyms
       if(all(mtrx(1:3, 1:3, isym) == identity(1:3, 1:3))) then
          if(.not. found) then
             ! if identity is not first, swap
             mtrx_temp(1:3, 1:3)  = mtrx(1:3, 1:3, 1)
             mtrx(1:3, 1:3, 1)    = mtrx(1:3, 1:3, isym)
             mtrx(1:3, 1:3, isym) = mtrx_temp(1:3, 1:3)

             tnp_temp(1:3)  = tnp(1:3, 1)
             tnp(1:3, 1)    = tnp(1:3, isym)
             tnp(1:3, isym) = tnp_temp(1:3)

             found = .true.
             write(0,'(a,i2)') 'WARNING: making identity first by swapping with symmetry op #', isym
          else
             call die("There is a duplicate identity in the symmetry operations.")
          endif
       endif
    enddo

    if(.not. found) then
       call die("Identity is not present in the list of symmetries.")
    endif

    POP_SUB(make_identity_symmetry_first)
    return
  end subroutine make_identity_symmetry_first

  !=====================================================================
  !> Bring symmetries into a standardized order.
  !! The identity is always the first one.
  subroutine sort_symmetries(nsyms, mtrx, tnp)
    integer,  intent(in)    :: nsyms
    integer,  intent(inout) :: mtrx(3, 3, 48)
    real(DP), intent(inout) :: tnp(3, 48)

    integer :: isym, ii, jj, factor, hash(48), order(48), mtrx_temp(3, 3, 48), identity(3,3)
    real(DP) :: tnp_temp(3, 48)

    PUSH_SUB(sort_symmetries)

    identity = reshape((/1, 0, 0, 0, 1, 0, 0, 0, 1/), shape(identity))

    do isym = 1, nsyms

       ! make sure the identity comes first
       if(all(mtrx(1:3, 1:3, isym) == identity(1:3, 1:3))) then
          hash(isym) = -1d9
          cycle
       endif

       hash(isym) = 0
       factor = 1
       do jj = 1, 3
          if(jj > 1) factor = factor * 3
          do ii = 1, 3
             if(ii > 1) factor = factor * 3
             hash(isym) = hash(isym) + mtrx(4 - ii, 4 - jj, isym) * factor
          enddo
       enddo
    enddo

    call sortix(nsyms, hash, order)

    do isym = 1, nsyms
       mtrx_temp(1:3, 1:3, isym) = mtrx(1:3, 1:3, order(isym))
       tnp_temp(1:3, isym) = tnp(1:3, order(isym))
    enddo

    mtrx(1:3, 1:3, 1:nsyms) = mtrx_temp(1:3, 1:3, 1:nsyms)
    tnp(1:3, 1:nsyms) = tnp_temp(1:3, 1:nsyms)

    POP_SUB(sort_symmetries)
    return
  end subroutine sort_symmetries

  !=====================================================================

  ! FHJ: Use the preprocessor to create the following routines:
  ! sortix_gvec, sortix_no_gvec, sortrx_gvec, sortrx_no_gvec

#define DTYPE integer
#define LABEL sortix_gvec
#define HAS_GVEC
#include "sort_inc.f90"
#undef LABEL
#undef HAS_GVEC
#define LABEL sortix_no_gvec
#include "sort_inc.f90"
#undef LABEL
#undef DTYPE

#define DTYPE real(DP)
#define LABEL sortrx_gvec
#define HAS_GVEC
#include "sort_inc.f90"
#undef LABEL
#undef HAS_GVEC
#define LABEL sortrx_no_gvec
#include "sort_inc.f90"
#undef LABEL
#undef DTYPE

end module sort_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

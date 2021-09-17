!==============================================================================
!
! Routines:
!
! (1) Xinvert_with_scalapack_d() Originally by JRD  Last Modified 02/2015 (FHJ)
!
!     This routine inverts a matrix which is already distributed in block
!     cyclic form with ScaLAPACK.
!
! (2) Xinvert_serial()           Originally by JRD  Last Modified 02/2015 (FHJ)
!
!     Inverts a matrix using LAPACK.
!
!==============================================================================

module inversion_m

  use global_m
  use lapack_m
  use scalapack_m
  implicit none

  private

  public :: &
#ifdef USESCALAPACK
    dinvert_with_scalapack,  &
    zinvert_with_scalapack,  &
#endif
    dinvert_serial,          &
    zinvert_serial

contains

!overrules flavor.mk
#undef CPLX
#include "f_defs.h"
#include "inversion_inc.f90"

#include "undef.h"

#define CPLX
#include "f_defs.h"
#include "inversion_inc.f90"

end module inversion_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

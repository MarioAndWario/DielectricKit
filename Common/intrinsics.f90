!================================================================================
!
! Modules:
!
! (1) intrinsics_m      Originally By DAS      Created 9/21/2011
!
!     Some compiler-dependent module usages and external definitions,
!     regarding accessing system calls. There is no actual code here.
!     These are preprocessor symbols defined in f_defs.h.
!
!================================================================================

#include "f_defs.h"

module intrinsics_m

  USEHOSTNAMEMOD
  USEOMPLIB
  SYSTEMMOD

  implicit none

  public

! note: these are split into two lines because ifort for some reason decrees of
! 'end interface' : "END statement must be only statement on line".

  HOSTNAME_INTERFACE
  HOSTNAME_INTERFACE2

  MCLOCK_INTERFACE
  MCLOCK_INTERFACE2

  IARGC_INTERFACE
  IARGC_INTERFACE2

  FTELL_INTERFACE
  FTELL_INTERFACE2

end module intrinsics_m

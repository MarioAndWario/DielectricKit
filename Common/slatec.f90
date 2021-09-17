!================================================================================
!
! Modules:
!
! (1) blas_m      Originally By DAS      Last Modified 1/13/2011 (das)
!
!     Interfaces for BLAS functions, taken from http://www.netlib.org/blas/
!     Every BLAS function used in the code should be listed here, and this
!     module should be used in every routine containing BLAS calls to ensure
!     the argument types are correct.
!
!     Note that if any array name from netlib.org is X, the interface will
!     be interpreted as a preprocessor macro and cause a compilation failure,
!     solved by changed to lower-case x.
!
!================================================================================

#include "f_defs.h"

module slatec_m
  public ! only interfaces in this module
  interface
     !> SLATEC/src/dgami.f
     !> Evaluate the incomplete gamma function defined by
     !> DGAMI = integral from T = 0 to X of EXP(-T) * T**(A-1.0)
     !> DGAMI is evaluated for positive values of A and non-negative values
     !> of X.  A slight deterioration of 2 or 3 digits accuracy will occur
     !> when DGAMI is very large or very small, because logarithmic variables
     !> are used.  The function and both arguments are double precision.
     DOUBLE PRECISION FUNCTION DGAMI(A,X)
       implicit none
       DOUBLE PRECISION A, X
     END FUNCTION DGAMI

     !> SLATEC/src/dgamic.f
     !> Calculate the complementary incomplete Gamma function.
     !> DGAMIC = integral from X to infinity of EXP(-T) * T**(A-1.)
     !> DGAMIC is evaluated for arbitrary real values of A and for non-
     !> negative values of X (even though DGAMIC is defined for X .LT. 0.0),
     !> except that for X = 0 and A .LE. 0.0, DGAMIC is undefined.
     DOUBLE PRECISION FUNCTION DGAMIC(A,X)
       implicit none
       DOUBLE PRECISION A, X
     END FUNCTION DGAMIC
     
  end interface

end module slatec_m

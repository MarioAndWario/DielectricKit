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

module blas_m

  public ! only interfaces in this module

  interface
    SUBROUTINE DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
      implicit none
      DOUBLE PRECISION ALPHA,BETA
      INTEGER K,LDA,LDB,LDC,M,N
      CHARACTER TRANSA,TRANSB
      DOUBLE PRECISION A(LDA,*),B(LDB,*),C(LDC,*)
    end subroutine DGEMM
  end interface

  interface
    SUBROUTINE ZGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
      implicit none
      DOUBLE COMPLEX ALPHA,BETA
      INTEGER K,LDA,LDB,LDC,M,N
      CHARACTER TRANSA,TRANSB
      DOUBLE COMPLEX A(LDA,*),B(LDB,*),C(LDC,*)
    end subroutine ZGEMM
  end interface

  interface
    SUBROUTINE DGEMV(TRANS,M,N,ALPHA,A,LDA,x,INCX,BETA,Y,INCY)
      implicit none
      DOUBLE PRECISION ALPHA,BETA
      INTEGER INCX,INCY,LDA,M,N
      CHARACTER TRANS
      DOUBLE PRECISION A(LDA,*),x(*),Y(*)
    end SUBROUTINE DGEMV
  end interface

  interface
    SUBROUTINE ZGEMV(TRANS,M,N,ALPHA,A,LDA,x,INCX,BETA,Y,INCY)
      implicit none
      DOUBLE COMPLEX ALPHA,BETA
      INTEGER INCX,INCY,LDA,M,N
      CHARACTER TRANS
      DOUBLE COMPLEX A(LDA,*),x(*),Y(*)
    end SUBROUTINE ZGEMV
  end interface

  interface
    SUBROUTINE ZHEMV(UPLO,N,ALPHA,A,LDA,x,INCX,BETA,Y,INCY)
      implicit none
      DOUBLE COMPLEX ALPHA,BETA
      INTEGER INCX,INCY,LDA,N
      CHARACTER UPLO
      DOUBLE COMPLEX A(LDA,*),x(*),Y(*)
    end SUBROUTINE ZHEMV
  end interface

  interface
    SUBROUTINE ZHPMV(UPLO,N,ALPHA,AP,x,INCX,BETA,Y,INCY)
      implicit none
      DOUBLE COMPLEX ALPHA,BETA
      INTEGER INCX,INCY,N
      CHARACTER UPLO
      DOUBLE COMPLEX AP(*),x(*),Y(*)
    end SUBROUTINE ZHPMV
  end interface

  interface
    SUBROUTINE ZHEMM(SIDE,UPLO,M,N,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
      implicit none
      DOUBLE COMPLEX ALPHA,BETA
      INTEGER LDA,LDB,LDC,M,N
      CHARACTER SIDE,UPLO
      DOUBLE COMPLEX A(LDA,*),B(LDB,*),C(LDC,*)
    end SUBROUTINE ZHEMM
  end interface

  interface blas_dot
    DOUBLE PRECISION FUNCTION DDOT(N,DX,INCX,DY,INCY)
      implicit none
      INTEGER INCX,INCY,N
      DOUBLE PRECISION DX(*),DY(*)
    end FUNCTION DDOT

    DOUBLE COMPLEX FUNCTION ZDOTC(N,ZX,INCX,ZY,INCY)
      implicit none
      INTEGER INCX,INCY,N
      DOUBLE COMPLEX ZX(*),ZY(*)
    end FUNCTION ZDOTC
  end interface blas_dot

  interface blas_zdotu
    DOUBLE COMPLEX FUNCTION ZDOTU(N,ZX,INCX,ZY,INCY)
      implicit none
      INTEGER INCX,INCY,N
      DOUBLE COMPLEX ZX(*),ZY(*)
    end FUNCTION ZDOTU
  end interface blas_zdotu
  
  interface
    SUBROUTINE DSCAL(N,DA,DX,INCX)
      implicit none
      DOUBLE PRECISION DA
      INTEGER INCX,N
      DOUBLE PRECISION DX(*)
    end SUBROUTINE DSCAL
  end interface

  interface
    SUBROUTINE ZSCAL(N,ZA,ZX,INCX)
      implicit none
      DOUBLE COMPLEX ZA
      INTEGER INCX,N
      DOUBLE COMPLEX ZX(*)
    end SUBROUTINE ZSCAL
  end interface

  interface blas_nrm2
    DOUBLE PRECISION FUNCTION DNRM2(N,x,INCX)
      implicit none
      INTEGER INCX,N
      DOUBLE PRECISION x(*)
    end FUNCTION DNRM2

    DOUBLE PRECISION FUNCTION DZNRM2(N,x,INCX)
      implicit none
      INTEGER INCX,N
      DOUBLE COMPLEX x(*)
    end FUNCTION DZNRM2
  end interface

end module blas_m

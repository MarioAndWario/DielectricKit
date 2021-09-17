!================================================================================
!
! Modules:
!
! (1) lapack_m      Originally By DAS      Last Modified 1/13/2011 (das)
!
!     Interfaces for LAPACK functions, taken from http://www.netlib.org/lapack/double
!     and http://www.netlib.org/lapack/complex16.
!     Every LAPACK function used in the code should be listed here, and this
!     module should be used in every routine containing LAPACK calls to ensure
!     the argument types are correct.
!
!     Note that if any array name from netlib.org is X, the interface will
!     be interpreted as a preprocessor macro and cause a compilation failure,
!     solved by changed to lower-case x.
!
!================================================================================

#include "f_defs.h"

module lapack_m

  public ! only interfaces in this module

  interface
     SUBROUTINE DGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )
       implicit none
       INTEGER            INFO, LDA, LDB, N, NRHS
       INTEGER            IPIV( * )
       DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
     end SUBROUTINE DGESV
  end interface

  interface
     SUBROUTINE ZGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )
       implicit none
       INTEGER            INFO, LDA, LDB, N, NRHS
       INTEGER            IPIV( * )
       COMPLEX*16         A( LDA, * ), B( LDB, * )
     end SUBROUTINE ZGESV
  end interface

  interface
     SUBROUTINE ZHPEVX( JOBZ, RANGE, UPLO, N, AP, VL, VU, IL, IU, &
          ABSTOL, M, W, Z, LDZ, WORK, RWORK, IWORK, IFAIL, INFO )
       implicit none
       CHARACTER          JOBZ, RANGE, UPLO
       INTEGER            IL, INFO, IU, LDZ, M, N
       DOUBLE PRECISION   ABSTOL, VL, VU
       INTEGER            IFAIL( * ), IWORK( * )
       DOUBLE PRECISION   RWORK( * ), W( * )
       COMPLEX*16         AP( * ), WORK( * ), Z( LDZ, * )
     end SUBROUTINE ZHPEVX
  end interface

  !> BII algorithm
  interface
     SUBROUTINE ZHEEVX( JOBZ, RANGE, UPLO, N, A, LDA, VL, VU, IL, IU, &
          ABSTOL, M, W, Z, LDZ, WORK, LWORK, RWORK, IWORK, IFAIL, INFO )
       implicit none
       CHARACTER          JOBZ, RANGE, UPLO
       INTEGER            IL, INFO, IU, LDA, LDZ, LWORK, M, N
       DOUBLE PRECISION   ABSTOL, VL, VU
       INTEGER            IFAIL( * ), IWORK( * )
       DOUBLE PRECISION   RWORK( * ), W( * )
       COMPLEX*16         A( LDA, * ), WORK( * ), Z( LDZ, * )
     end SUBROUTINE ZHEEVX
  end interface

  !> MRRR algorithm
  interface
     SUBROUTINE ZHEEVR( JOBZ, RANGE, UPLO, N, A, LDA, VL, VU, IL, IU, &
          ABSTOL, M, W, Z, LDZ, ISUPPZ, WORK, LWORK, RWORK, LRWORK, IWORK, LIWORK, INFO )
       implicit none
       CHARACTER          JOBZ, RANGE, UPLO
       INTEGER            IL, INFO, IU, LDA, LDZ, LWORK, M, N, LRWORK, LIWORK
       DOUBLE PRECISION   ABSTOL, VL, VU
       INTEGER            IWORK( * ), ISUPPZ( * )
       DOUBLE PRECISION   RWORK( * ), W( * )
       COMPLEX*16         A( LDA, * ), WORK( * ), Z( LDZ, * )
     end SUBROUTINE ZHEEVR
  end interface

  !> DC algorithm
  !> on exit the array (A) is overwritten by the unitary matrix Z which contains the eigenvectors of A.
  interface
     SUBROUTINE ZHEEVD( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, RWORK, LRWORK, IWORK, LIWORK, INFO )
       implicit none
       CHARACTER          JOBZ, UPLO
       INTEGER            INFO, LDA, LDZ, LWORK, N, LRWORK, LIWORK
       INTEGER            IWORK( * )
       DOUBLE PRECISION   RWORK( * ), W( * )
       COMPLEX*16         A( LDA, * ), WORK( * )
     end SUBROUTINE ZHEEVD
  end interface  
  
  !> Computes selected eigenvalues and, optionally, eigenvectors of a Hermitian matrix.
  interface
     SUBROUTINE HEEVX(A, W, UPLO, Z, VL, VU, IL, IU, M, IFAIL, ABSTOL, INFO )
       implicit none
       CHARACTER          UPLO
       INTEGER            INFO, LWORK, M, IL, IU, VL, VU
       INTEGER            IFAIL( * )
       DOUBLE PRECISION   ABSTOL
       DOUBLE PRECISION   W( * )
       COMPLEX*16         A( * ), Z( * )
     end SUBROUTINE HEEVX
  end interface

  interface
     SUBROUTINE DGETRF( M, N, A, LDA, IPIV, INFO )
       implicit none
       INTEGER            INFO, LDA, M, N
       INTEGER            IPIV( * )
       DOUBLE PRECISION   A( LDA, * )
     end SUBROUTINE DGETRF
  end interface

  interface
     SUBROUTINE DGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO )
       implicit none
       INTEGER            INFO, LDA, LWORK, N
       INTEGER            IPIV( * )
       DOUBLE PRECISION   A( LDA, * ), WORK( * )
     end SUBROUTINE DGETRI
  end interface

  interface
     SUBROUTINE ZGETRF( M, N, A, LDA, IPIV, INFO )
       implicit none
       INTEGER            INFO, LDA, M, N
       INTEGER            IPIV( * )
       COMPLEX*16         A( LDA, * )
     end SUBROUTINE ZGETRF
  end interface

  interface
     SUBROUTINE ZGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO )
       implicit none
       INTEGER            INFO, LDA, LWORK, N
       INTEGER            IPIV( * )
       COMPLEX*16         A( LDA, * ), WORK( * )
     end SUBROUTINE ZGETRI
  end interface

  interface
     SUBROUTINE ZHEEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, RWORK, INFO )
       implicit none
       CHARACTER          JOBZ, UPLO
       INTEGER            INFO, LDA, LWORK, N
       DOUBLE PRECISION   RWORK( * ), W( * )
       COMPLEX*16         A( LDA, * ), WORK( * )
     end SUBROUTINE ZHEEV
  end interface

  interface
     SUBROUTINE DSYEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO )
       implicit none
       CHARACTER          JOBZ, UPLO
       INTEGER            INFO, LDA, LWORK, N
       DOUBLE PRECISION   A( LDA, * ), W( * ), WORK( * )
     end SUBROUTINE DSYEV
  end interface

  interface
     SUBROUTINE DGEQRF( M, N, A, LDA, TAU, WORK, LWORK, INFO )
       implicit none
       INTEGER            INFO, LDA, LWORK, M, N
       DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
     end SUBROUTINE DGEQRF
  end interface

  interface
     SUBROUTINE ZGEQRF( M, N, A, LDA, TAU, WORK, LWORK, INFO )
       implicit none
       INTEGER            INFO, LDA, LWORK, M, N
       COMPLEX*16         A( LDA, * ), TAU( * ), WORK( * )
     end SUBROUTINE ZGEQRF
  end interface

  interface
     SUBROUTINE DORGQR( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
       implicit none
       INTEGER            INFO, K, LDA, LWORK, M, N
       DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
     end SUBROUTINE DORGQR
  end interface

  interface
     SUBROUTINE ZUNGQR( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
       implicit none
       INTEGER            INFO, K, LDA, LWORK, M, N
       COMPLEX*16         A( LDA, * ), TAU( * ), WORK( * )
     end SUBROUTINE ZUNGQR
  end interface

  interface
     SUBROUTINE DPOTRF( UPLO, N, A, LDA, INFO )
       implicit none
       CHARACTER          UPLO
       INTEGER            INFO, LDA, N
       DOUBLE PRECISION   A( LDA, * )
     end SUBROUTINE DPOTRF
  end interface

  interface
     DOUBLE PRECISION FUNCTION DLAMCH( CMACH )
       implicit none
       CHARACTER          CMACH
     END FUNCTION DLAMCH
  end interface

end module lapack_m

#define IMAG(x) AIMAG(x)
#define CONJG(x) CONJG(x)
#define CMPLX(x,y) CMPLX(x,y,kind=DPC)
#define ABS2(x) (DBLE(x)**2 + AIMAG(x)**2)
!The following macro puts any point/array in the [-0.5, 0.5) range:
#define MIN_RANGE(x) (x - FLOOR(x + 0.5d0))
!The following macro puts any point/array in the [0, 1) range:
#define UNIT_RANGE(x) (x - FLOOR(x))
!Integer division of a/b rounded up*/
#define DIVUP(a,b) ((a+b-1)/(b))
!Rounds a up to the smallest multiple of b*/
#define ROUNDUP(a,b) (DIVUP(a,b)*(b))

! disable Fortran OMP pragmas if not -DOMP*/
! note: C standard does not permit $ in identifiers, however this seems acceptable
!   as an extension, for all versions of cpp I tried. --DAS
#ifndef OMP
#define $OMP disabled
#define $omp disabled
#endif

#define MPI_REAL_DP MPI_DOUBLE_PRECISION
#define MPI_COMPLEX_DPC MPI_DOUBLE_COMPLEX

#define RFLAVOR "Real"
#define CFLAVOR "Complex"

#ifdef CPLX
#define MYFLAVOR CFLAVOR
#define SCALAR COMPLEX(DPC)
#define MPI_SCALAR MPI_COMPLEX_DPC
#define MYCONJG(x) CONJG(x)
#define MYREAL(x) REAL(x,kind=DP)
#define MYABS2(x) ABS2(x)
#define X(x) z ## x
#define pX(x) pz ## x
#define ZERO (0.0d0,0.0d0)
#define ONE (1.0d0,0.0d0)
#define SCALARIFY(x) (x)
#define SCALARIFY2(x,y) DCMPLX(x,y)
!#define SCALARIFY2(x,y) CMPLX(x,y)
#define COMPLEXIFY(x) (x)
#define SCALARSIZE 2
#define ONLYIFCPLX(x) x, 
#else
#define MYFLAVOR RFLAVOR
#define SCALAR REAL(DP)
#define MPI_SCALAR MPI_REAL_DP
#define MYCONJG(x) (x)
#define MYREAL(x) (x)
#define MYABS2(x) ((x)**2)
#define X(x) d ## x
#define pX(x) pd ## x
#define ZERO 0.0d0
#define ONE 1.0d0
!#define SCALARIFY(x) real(x,kind=DP)
#define SCALARIFY(x) DBLE(x,kind=DP)
!#define SCALARIFY2(x,y) real(x,kind=DP)
#define SCALARIFY2(x,y) DBLE(x,kind=DP)
!#define COMPLEXIFY(x) cmplx(x,kind=DPC)
#define COMPLEXIFY(x) DCMPLX(x,kind=DPC)
#define SCALARSIZE 1
#define ONLYIFCPLX(x)
#endif

! truncate spaces in string 
!#!define TRUNC(s) trim(adjustl(s))

! Sun compiler has a length limit of 132 characters and won`t support these macros 
#if defined(DEBUG) && (! defined(SUN) )
! Use this instead of the intrinsic 'allocate' 
#define SAFE_ALLOCATE(x,y) \
allocate(x y,stat=alc); \
call alloc_check(alc,SIZEOF(x),#x,__FILE__,__LINE__,.true.)

! Do not use this directly 
#define MY_DEALLOCATE(x) \
sz=SIZEOF(x);deallocate(x,stat=alc); \
call alloc_check(alc,sz,#x,__FILE__,__LINE__,.false.)

#else
! No checking for faster performance, if not in debug mode 
#define SAFE_ALLOCATE(x,y) allocate(x y)
#define MY_DEALLOCATE(x) deallocate(x)
#endif

! Use this instead of the intrinsic 'deallocate' for pointers 
#define SAFE_DEALLOCATE_P(x) \
if(associated(x))then;\
MY_DEALLOCATE(x);\
nullify(x);\
endif

! Use this instead of the intrinsic 'deallocate' for arrays 
#define SAFE_DEALLOCATE(x) \
if(allocated(x))then;\
MY_DEALLOCATE(x);\
endif

!the TOSTRING macro converts a macro into a string
#define STRINGIFY(x) #x
#define TOSTRING(x)  STRINGIFY(x)

#ifdef DEBUG
#define PUSH_SUB(routine) call push_sub(__FILE__+"."+TOSTRING(routine))
#define POP_SUB(routine) call pop_sub(__FILE__+"."+TOSTRING(routine))
#else
#define PUSH_SUB(routine) 
#define POP_SUB(routine) 
#endif

! ! deprecated identifiers 
! #define DFLOAT @@_use_dble_instead_@@
! #define DREAL @@_use_dble_instead_@@
! #define REAL @@_use_dble_instead_@@
! !#define DIMAG @@_use_IMAG_instead_@@
! #define DCONJG @@_use_CONJG_instead_@@
! !#define DCMPLX @@_use_CMPLX_instead_@@
! #define dfloat @@_use_dble_instead_@@
! #define dreal @@_use_dble_instead_@@
! !#define dimag @@_use_IMAG_instead_@@
! #define dconjg @@_use_CONJG_instead_@@
! #define dcmplx @@_use_CMPLX_instead_@@

#include "compiler.h"

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

!============================================================================
!
! MODULE: fftw_m, originally by DAS 1/14/2011
!
!>   Routines used with FFTW, as well as interfaces for library calls.
!
! DESCRIPTION:
!>   No FFTW calls should exist outside this routine: the wrapper routines
!!   here should be used everywhere.
!!
!!   Interfaces for FFTW2 functions are formulated from fftw-2.1.5/fftw/fftwf77.c
!!   http://www.fftw.org/fftw2_doc/fftw_5.html. Every FFTW2 function used should
!!   have an interface to ensure the argument types are correct.
!!   Contains include file fftw_f77.i which has parameters used in FFTW2 calls.
!!
!!   FFTW3 provides its own file of constants and interfaces, fftw3.f03.
!!   See http://fftw.org/doc/Overview-of-Fortran-interface.html
!
!============================================================================

#include "f_defs.h"

module fftw_m

#ifdef USEFFTW3
  !> use of this is recommended by FFTW3 documentation. It causes no harm for
  !! FFTW2, but some compilers (e.g. Open64) do not have it available, so we
  !! will keep it hidden for FFTW2 so as not to have to solve that problem yet.
  use, intrinsic :: iso_c_binding
#endif
  use global_m
  implicit none

  private

#ifdef USEFFTW3
!> It is better to use this one which has interfaces too, rather than fftw3.f which has only constants
  include 'fftw3.f'
#else
  include 'fftw_f77.i'
#endif

#ifdef USEFFTW3
    integer*8, private :: fft_plan = 0
    integer, private :: ifirst = 0
    integer, private :: num_fft_threads = 0
    integer :: iret
#else
!> fftw_plan type in C is recommended to be integer*8 by FFTW documentation
    integer*8, private :: plus_plan = 0
    integer*8, private :: minus_plan = 0
#endif
    integer, private :: Nfftold(3) = 0

#ifndef USEFFTW3

  interface
    subroutine fftwnd_f77_create_plan(p, rank, n, idir, flags)
      implicit none
      integer*8 :: p
      integer :: rank, n(*), idir, flags
    end subroutine fftwnd_f77_create_plan
  end interface

  interface
    subroutine fftwnd_f77_one(p, in, out)
      implicit none
      integer*8 :: p
      complex*16 :: in(*)
      integer :: out
!< The argument is really complex*16 out(*), but we only use in-place transforms,
!! in which case this argument is ignored. For simplicity we just pass it 0.    
    end subroutine fftwnd_f77_one
  end interface
  
  interface
    subroutine fftwnd_f77_destroy_plan(p)
      implicit none
      integer*8 :: p
    end subroutine fftwnd_f77_destroy_plan
  end interface

#endif

  public ::            &
    check_FFT_size,    &
    setup_FFT_sizes,   &
    gvec_to_fft_index, &
    put_into_fftbox,   &
    get_from_fftbox,   &
    do_FFT,            &
    conjg_fftbox,      &
    multiply_fftboxes, &
    destroy_fftw_plans

  interface put_into_fftbox
    module procedure dput_into_fftbox, zput_into_fftbox
  end interface put_into_fftbox

  interface get_from_fftbox
    module procedure dget_from_fftbox, zget_from_fftbox
  end interface get_from_fftbox

contains

!> Originally by gsm      Last Modified: 4/10/2010 (gsm)
!!     Best FFT grid dimension is given by 2^a*3^b*5^c*7^d*11^e*13^f
!!     where a,b,c,d are arbitrary and e,f are 0 or 1
!! Ref: http://www.fftw.org/fftw2_doc/fftw_3.html
!!     On entry
!!             Nfft = FFT grid dimension to test
!!             Nfac = number of factors to test
!!     On exit
!!             check_FFT_size = .true. if good FFT grid dimension
  logical function check_FFT_size(Nfft, Nfac) 
    integer, intent(in) :: Nfft, Nfac
  
    integer :: remainder, product, ifac, ipow, maxpow
    integer, parameter :: maxfac = 6
    integer :: pow(maxfac)
    integer, parameter :: fac(maxfac) = (/ 2, 3, 5, 7, 11, 13 /)
  
    PUSH_SUB(check_FFT_size)
  
    if(Nfft .lt. 1 .or. Nfac .lt. 1 .or. Nfac .gt. maxfac) then
      call die('check_FFT_size input')
    endif
  
    remainder = Nfft
    do ifac = 1, maxfac
      pow(ifac) = 0
    enddo
  
    do ifac = 1, Nfac
      maxpow = int(log(dble(remainder)) / log(dble(fac(ifac)))) + 1
      do ipow = 1, maxpow
        if (mod(remainder, fac(ifac)) .eq. 0) then
          remainder = remainder / fac(ifac)
          pow(ifac) = pow(ifac) + 1
        endif
      enddo
    enddo
  
    product = remainder
    do ifac = 1, Nfac
      do ipow = 1, pow(ifac)
        product = product * fac(ifac)
      enddo
    enddo
    if (product .ne. Nfft) then
      call die('Internal error in check_FFT_size; factorization failed')
    endif
  
    check_FFT_size = remainder .eq. 1 .and. pow(5) .le. 1 .and. pow(6) .le. 1
  
    POP_SUB(check_FFT_size)
  
    return
  end function check_FFT_size

!> The former "fft_routines.f90"
!! Sohrab Ismail-Beigi   Feb 28 2001
!!
!! There are a set of Fast Fourier-related routines that are used
!! to compute the matrix elements of the type <nk|e^(i*G.r)|mk`>.
!! For many G-vectors, FFTs will be the fastest way to compute them.
!!
!! The FFTW (http://www.fftw.org) suite of routines do the actual work.
!! Most of what is below is interfacing code and routines that simplify
!! small and useful tasks.
 !
!!
!! Given gvec%FFTgrid(1:3) values (in FFTgrid), finds appropriate FFT box
!! sizes to use in Nfft(1:3).  scale = 1/(Nfftx*Nffty*Nfftz).
!!
  subroutine setup_FFT_sizes(FFTgrid,Nfft,scale)
    integer, intent(in) :: FFTgrid(3)
    integer, intent(out) :: Nfft(3)
    real(DP), intent(out) :: scale

    integer, parameter :: Nfac = 3
    integer :: i

    PUSH_SUB(setup_FFT_sizes)

    do i=1,3
      Nfft(i) = FFTgrid(i)
      do while (.not. check_FFT_size(Nfft(i), Nfac))
        Nfft(i) = Nfft(i) + 1
      enddo
    enddo
    scale = 1.0d0/product(Nfft(1:3))
    
    POP_SUB(setup_FFT_sizes)

    return
  end subroutine setup_FFT_sizes

!> Takes the G-vector g(1:3) and FFT box size Nfft(1:3) and finds the
!! point idx(1:3) in the box corresponding to that G-vector.
!!
  subroutine gvec_to_fft_index(g,idx,Nfft)
    integer, intent(in) :: g(3), Nfft(3)
    integer, intent(out) :: idx(3)

! no push/pop since called too frequently.

    idx(1:3) = g(1:3) + 1
    
    if (g(1) < 0) idx(1) = Nfft(1) + idx(1)
    if (g(2) < 0) idx(2) = Nfft(2) + idx(2)
    if (g(3) < 0) idx(3) = Nfft(3) + idx(3)

    return
  end subroutine gvec_to_fft_index

!> Do an FFT on the fftbox in place:  destroys contents of fftbox
!! and replaces them by the Fourier transform.
!!
!! The FFT done is:
!!
!!   fftbox(p) <- sum_j { fftbox(j)*e^{sign*i*j.p} }
!!
!! where j and p are integer 3-vectors ranging over Nfft(1:3).
!!
  subroutine do_FFT(fftbox, Nfft, sign)
    complex(DPC), intent(inout) :: fftbox(:,:,:)
    integer, intent(in) :: Nfft(3)
    integer, intent(in) :: sign

    character(len=100) :: str

    PUSH_SUB(do_FFT)

#ifdef USEFFTW3

!JRD To be removed
    !complex(DPC), allocatable :: fftbox2(:,:,:)

#ifdef OMP
    if (ifirst .eq. 0) then
      call dfftw_init_threads(iret)
      ifirst = 1
      write(str,'(a,i0)') 'Setup threaded FFTs. Return: ', iret
      call logit(str)

!#!$OMP PARALLEL
!#      num_fft_threads = OMP_GET_NUM_THREADS()
!#!$OMP END PARALLEL

      num_fft_threads = peinf%nthreads
      call dfftw_plan_with_nthreads(num_fft_threads)
      write(str,'(a,i0)') 'Doing threaded FFTs. Num threads: ', num_fft_threads
      call logit(str)
    endif 
#endif

    if (peinf%verb_max) then
      write(str,'(a,2(i0," x "),i0,a)') 'Creating ', Nfft(1:3), ' FFTW plans.'
      call logit(str)
    endif
      
    if(peinf%inode.eq.0) call timacc(93,1)
    if (sign == 1) then
      call dfftw_plan_dft(fft_plan,3,Nfft,fftbox,fftbox,FFTW_BACKWARD,FFTW_ESTIMATE)
      !call dfftw_plan_dft(fft_plan,3,Nfft,fftbox,fftbox2,FFTW_BACKWARD,FFTW_MEASURE)
      !call dfftw_plan_dft_3d(fft_plan,Nfft(1),Nfft(2),Nfft(3),fftbox,fftbox2,FFTW_BACKWARD,FFTW_MEASURE)
    else if (sign == -1) then
      call dfftw_plan_dft(fft_plan,3,Nfft,fftbox,fftbox,FFTW_FORWARD,FFTW_ESTIMATE)
      !call dfftw_plan_dft(fft_plan,3,Nfft,fftbox,fftbox2,FFTW_FORWARD,FFTW_MEASURE)
      !call dfftw_plan_dft_3d(fft_plan,Nfft(1),Nfft(2),Nfft(3),fftbox,fftbox2,FFTW_FORWARD,FFTW_MEASURE)
    else
      call die('sign is not 1 or -1 in do_FFT')
    endif
    if(peinf%inode.eq.0) call timacc(93,2)

    if(peinf%inode.eq.0) call timacc(94,1)
    call dfftw_execute_dft(fft_plan,fftbox,fftbox)
    if(peinf%inode.eq.0) call timacc(94,2)

    ! otherwise there is a memory leak
    call dfftw_destroy_plan(fft_plan)
    !call dfftw_cleanup_threads()

    Nfftold(:) = -1

!JRD To be removed
    !fftbox(:,:,:)=fftbox2(:,:,:)
    !SAFE_DEALLOCATE(fftbox2)

#else

!JRD Determine dimensions

!    ndim = 0
!    do idim=1,3
!      if (Nfft(idim) .gt. 1) ndim = ndim + 1
!    enddo
!    if (ndim .lt. 1) call die('fft of dimension less than 1')
!    if (ndim .gt. 3) call die('fft of dimension more than 3')

    if(any(Nfftold(1:3) .ne. Nfft(1:3))) then
      write(str,'(a,2(i0," x "),i0,a)') 'Creating ', Nfft(1:3), ' FFTW plans.'
      call logit(str)
      
      ! otherwise there is a memory leak
      if(any(Nfftold(1:3) /= -1)) call destroy_fftw_plans()
      
      call fftwnd_f77_create_plan(plus_plan,3,Nfft,FFTW_BACKWARD, &
        FFTW_MEASURE+FFTW_IN_PLACE+FFTW_USE_WISDOM)
      call fftwnd_f77_create_plan(minus_plan,3,Nfft,FFTW_FORWARD, &
        FFTW_MEASURE+FFTW_IN_PLACE+FFTW_USE_WISDOM)
      Nfftold(1:3) = Nfft(1:3)
      call logit('Done creating plans')
    endif

    if (sign == 1) then
      call fftwnd_f77_one(plus_plan,fftbox,0)
    else if (sign == -1) then
      call fftwnd_f77_one(minus_plan,fftbox,0)
    else
      call die('sign is not 1 or -1 in do_FFT')
    endif

#endif

    POP_SUB(do_FFT)

    return
  end subroutine do_FFT

  subroutine destroy_fftw_plans()
    character(len=100) :: str

    ! FFTW plan was never created
    if(all(Nfftold(1:3) == 0)) return

    PUSH_SUB(destroy_fftw_plans)

    if(all(Nfftold(1:3) == -1)) then 
#ifdef USEFFTW3
#ifdef OMP
      ifirst = 0
      call dfftw_cleanup_threads()
#endif
#endif
      POP_SUB(destroy_fftw_plans)
      return
      ! call die("Cannot destroy FFTW plan for a second time.")
    endif

    write(str,'(a,2(i0," x "),i0,a)') 'Destroying ', Nfftold(1:3), ' FFTW plans.'
    call logit(str)
    Nfftold(1:3) = -1 ! make clear there is no plan anymore so we do not try to destroy twice

#ifdef USEFFTW3
    call dfftw_destroy_plan(fft_plan)
#ifdef OMP
    call dfftw_cleanup_threads()
#endif
#else
    call fftwnd_f77_destroy_plan(plus_plan)
    call fftwnd_f77_destroy_plan(minus_plan)
#endif
    ! should forget wisdom here, but I cannot figure out how... --DAS

    POP_SUB(destroy_fftw_plans)
    return
  end subroutine destroy_fftw_plans

!> Complex conjugate contents of FFT box
!
  subroutine conjg_fftbox(fftbox,Nfft)
    integer, intent(in) :: Nfft(3)
    complex(DPC), intent(inout) :: fftbox(:,:,:) !< (Nfft(1), Nfft(2), Nfft(3))
    ! for some reason, absoft segfaults if dims specified for fftbox as above right

    integer :: ix,iy,iz

    PUSH_SUB(conjg_fftbox)

    if(peinf%inode.eq.0) call timacc(96,1)

!$OMP PARALLEL DO PRIVATE (ix,iy,iz)
    do iz = 1, Nfft(3)
    do iy = 1, Nfft(2)
    do ix = 1, Nfft(1)
      fftbox(ix,iy,iz) = CONJG(fftbox(ix,iy,iz))
    enddo
    enddo
    enddo
!$OMP END PARALLEL DO

    if(peinf%inode.eq.0) call timacc(96,2)

    POP_SUB(conjg_fftbox)

    return
  end subroutine conjg_fftbox

!> Multiply contents of two fft boxes, result into fftbox2
!
  subroutine multiply_fftboxes(fftbox1, fftbox2, Nfft)
    integer, intent(in) :: Nfft(3)
    complex(DPC), intent(in) :: fftbox1(:,:,:) !< (Nfft(1), Nfft(2), Nfft(3))
    complex(DPC), intent(inout) :: fftbox2(:,:,:) !< (Nfft(1), Nfft(2), Nfft(3))

    integer :: ix,iy,iz

    PUSH_SUB(multiply_fftboxes)

    if(peinf%inode.eq.0) call timacc(95,1)

    !forall(ix=1:Nfft(1), iy=1:Nfft(2), iz=1:Nfft(3)) &
    !  fftbox2(ix,iy,iz) = fftbox1(ix,iy,iz) * fftbox2(ix,iy,iz)

!$OMP PARALLEL DO PRIVATE (ix,iy,iz)
    do iz = 1, Nfft(3)
    do iy = 1, Nfft(2)
    do ix = 1, Nfft(1)
      fftbox2(ix,iy,iz) = fftbox1(ix,iy,iz) * fftbox2(ix,iy,iz)
    enddo
    enddo
    enddo
!$OMP END PARALLEL DO

    if(peinf%inode.eq.0) call timacc(95,2)

    POP_SUB(multiply_fftboxes)

    return
  end subroutine multiply_fftboxes

#include "undef.h"
!overrules flavor.mk
#undef CPLX
#include "f_defs.h"
#include "fftw_inc.f90"

#include "undef.h"

#define CPLX
#include "f_defs.h"
#include "fftw_inc.f90"

end module fftw_m

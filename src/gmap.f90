#include "f_defs.h"

!!===================================================================
!!
!! Module:
!!
!!    gmap          Last Modified by Meng Wu (2020)
!!
!!     Find the index array and phases needed to map G-vectors for a
!!     wavefunction at one k-point to the corresponding wavefunction
!!     at a symmetry-related k-point.
!!
!!===================================================================

module gmap_m

  use global_m
  use misc_m

  implicit none

  private

  public :: gmap

  interface gmap
     module procedure dgmap, zgmap
  end interface gmap
  
contains

  subroutine dgmap(gvec, syms, ngk, itran, kgq, isortc, isorti, ind, phase)
    type (gspace), intent(in) :: gvec         ! uses gvec%components(:,isrtc(1:ngk)), and index_vec
    type (symmetry), intent(in) :: syms       ! uses syms%mtrx(:,:,itran) & syms%tnp(:,itran)
    integer, intent(in) :: ngk                ! number of g-vector entries in a wavefunction
    integer, intent(in) :: itran              ! index of transformation
    integer, intent(in) :: kgq(3)             ! an umklapp vector (i.e. integer 3-vector)
    integer, intent(in) :: isortc(:)          ! index array for R(q) (1:ngk)
    integer, intent(in) :: isorti(:)          ! inverse index array for q (1:gvec%ng)
    integer, intent(out) :: ind(:)            ! indices for the vectors inv(symm(itran))*(g+kgq) (ngk)
    real(DP), intent(out) :: phase(:)         ! exp(-i*(g+kgq).dot.syms%tnp(itran)) (ngk)

    PUSH_SUB(dgmap)

    call gmap_base( gvec, syms, ngk, itran, kgq, isortc, isorti, ind, dphase = phase)

    POP_SUB(dgmap)
    return
  end subroutine dgmap
  
  subroutine zgmap(gvec, syms, ngk, itran, kgq, isortc, isorti, ind, phase)
    type (gspace), intent(in) :: gvec         ! uses gvec%components(:,isrtc(1:ngk)), and index_vec
    type (symmetry), intent(in) :: syms       ! uses syms%mtrx(:,:,itran) & syms%tnp(:,itran)
    integer, intent(in) :: ngk                ! number of g-vector entries in a wavefunction
    integer, intent(in) :: itran              ! index of transformation
    integer, intent(in) :: kgq(3)             ! an umklapp vector (i.e. integer 3-vector)
    integer, intent(in) :: isortc(:)          ! index array for R(q) (1:ngk)
    integer, intent(in) :: isorti(:)          ! inverse index array for q (1:gvec%ng)
    integer, intent(out) :: ind(:)            ! indices for the vectors inv(symm(itran))*(g+kgq) (ngk)
    complex(DPC), intent(out) :: phase(:)     ! exp(-i*(g+kgq).dot.syms%tnp(itran)) (ngk)

    PUSH_SUB(zgmap)

    call gmap_base(gvec, syms, ngk, itran, kgq, isortc, isorti, ind, zphase = phase)

    POP_SUB(zgmap)
    return
  end subroutine zgmap
  
  !==================================================================

  !! gmap_base() allow ig_old exceed ngk, which is the largest index of ig_new
  subroutine gmap_base(gvec, syms, ngk, itran, kgq, isortc, isorti, ind, dphase, zphase)
    type (gspace), intent(in) :: gvec         ! uses gvec%components(:,isrtc(1:ngk)), and index_vec
    type (symmetry), intent(in) :: syms       ! uses syms%mtrx(:,:,itran) & syms%tnp(:,itran)
    integer, intent(in) :: ngk                ! number of g-vector entries in a wavefunction
    integer, intent(in) :: itran              ! index of transformation
    integer, intent(in) :: kgq(3)             ! an umklapp vector (i.e. integer 3-vector)
    integer, intent(in) :: isortc(:)          ! index array for R(q) (1:ngk)
    integer, intent(in) :: isorti(:)          ! inverse index array for q (1:gvec%ng)
    integer, intent(out) :: ind(:)            ! indices for the vectors inv(symm(itran))*(g+kgq) (ngk)
    real(DP),     optional, intent(out) :: dphase(:) ! exp(-i*(g+kgq).dot.syms%tnp(itran)) (ngk)
    complex(DPC), optional, intent(out) :: zphase(:) ! exp(-i*(g+kgq).dot.syms%tnp(itran)) (ngk)
    integer :: ig, kd(3), kadd, kgrad, kgrad1
    integer :: kg(3), kgr(3)
    integer, dimension(3,3) :: inv_mtrx_reci_int
    real(DP), dimension(3,3) :: inv_mtrx_reci, mtrx_reci, diff_inv_mtrx_reci
    logical :: flag_inv
    real(DP) :: fi
    PUSH_SUB(gmap_base)

    if(present(dphase) .and. present(zphase)) then
       call die("gmap: cannot pass both dphase and zphase")
    else if(.not. present(dphase) .and. .not. present(zphase)) then
       call die("gmap: must pass either dphase or zphase")
    endif

    if (ngk > gvec%ng) call die("gmap: ngk (wfn cutoff) is greater than gvec%ng (rho cutoff)")
    if (ubound(isorti, 1) < gvec%ng) call die("gmap: isorti size < gvec%ng")
    if (any(isorti(1:gvec%ng) > gvec%ng)) call die("gmap: isorti cannot be greater than gvec%ng.")
    if (ubound(isortc, 1) < ngk)      call die("gmap: isortc size < ngk")
    if (any(isortc(1:ngk) < 1))       call die("gmap: isortc cannot be less than 1.")
    if (any(isortc(1:ngk) > gvec%ng)) call die("gmap: isortc cannot be greater than ng.")
    if (ubound(gvec%index_vec, 1) /= gvec%nFFTgridpts) call die("gmap: gvec%index_vec has wrong size")
    if (any(gvec%index_vec(1:gvec%nFFTgridpts) < 0)) call die("gmap: index_vec cannot be less than 0")
    if (any(gvec%index_vec(1:gvec%nFFTgridpts) > gvec%ng)) call die("gmap: index_vec cannot be greater than ng")
    if(present(dphase)) then
       if(ubound(dphase, 1) < ngk) call die("gmap: dphase size < ngk")
    else
       if(ubound(zphase, 1) < ngk) call die("gmap: zphase size < ngk")
    endif
    if(ubound(ind, 1) < ngk) call die("gmap: ind size < ngk")

    mtrx_reci(1:3,1:3) = dble(syms%mtrx_reci(1:3,1:3,itran))
    call M33INV(mtrx_reci,inv_mtrx_reci,flag_inv)
    if (.not. flag_inv) then
       call die("gmap_2: mtrx_reci not invertible.", only_root_writes=.true.)
    endif

    inv_mtrx_reci_int(:,:) = NINT(inv_mtrx_reci(:,:))
    diff_inv_mtrx_reci(:,:) = DBLE(inv_mtrx_reci_int(:,:)) - inv_mtrx_reci(:,:)
    IF ( NORM2(diff_inv_mtrx_reci) > TOL_Zero ) THEN
       call die("gmap_2: inv_mtrx_reci not integer ", only_root_writes=.true.)
    ENDIF   
    
    !! Loop over g-vectors in new wave function
    do ig = 1, ngk
       kg(1:3) = gvec%components(1:3, isortc(ig)) + kgq(1:3)
       kgr(1:3) = MATMUL(inv_mtrx_reci_int(1:3, 1:3), kg(1:3))      
       kd(1:3) = kgr(1:3) + gvec%FFTgrid(1:3) / 2 + 1
       if (any(kd(1:3) .lt. 1 .or. kd(1:3) .gt. gvec%FFTgrid(1:3))) then
          call die('gmap: kd out of bounds')
       endif

       kadd = ((kd(1) - 1) * gvec%FFTgrid(2) + kd(2) - 1) * gvec%FFTgrid(3) + kd(3)
       !! gvec%index_vec(ig_FFT) = ig_gvec \in [1, gvec%ng]
       kgrad1 = gvec%index_vec(kadd)
       if (kgrad1 .lt. 1 .or. kgrad1 .gt. gvec%ng) then
          write(0,*) 'itran = ', itran, 'ig = ', ig, ', kadd = ', kadd, ', kgrad1 = ', kgrad1
          call die('gmap: G-vectors falling outside of the charge-density G-sphere')
       endif
       !! isorti(ig_gvec) = ig_old
       kgrad = isorti(kgrad1)
       !! kgrad can be 0 if kgrad1 is outside the cutoff of wfn_old!
       ind(ig) = kgrad
       fi = 2.0D0 * PI_D * DOT_PRODUCT(DBLE(gvec%components(1:3, isortc(ig))),syms%tnp(1:3,itran))
       if(present(zphase)) then
          zphase(ig) = DCMPLX(cos(fi), -sin(fi))
       else
          dphase(ig) = cos(fi)
          if(abs(abs(dphase(ig)) - 1.0D0) .gt. TOL_Small) then
             write(0,'(a,i8,a,f12.8,a)') 'phase(', ig, ') = ', dphase(ig), ' != +/- 1'
             call die("Illegal non-unity phase in gmap, error in fractional translation.")
          endif

          if(abs(sin(fi)) .gt. TOL_Small) then
             write(0,'(a,i8,a,f12.8)') 'Im phase(', ig, ') = ', -sin(fi)
             call die("Illegal complex phase in gmap, error in fractional translation.")
          endif         
       endif
    enddo

    POP_SUB(gmap_base)
    return
  end subroutine gmap_base
end module gmap_m

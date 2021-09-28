#include "f_defs.h"

!===========================================================================
!
! (1) mtxel()
!
!     Compute matrix elements (gme) for valence state iv with all
!     conduction bands and for all G-vectors.
!
!                <c,k,ispin|exp{i(q+G).r}|v,(k-q),ispin> = M_{vc}(k,q,G)
!
!     On exit,
!       pol%gme(band,g-vector,spin) = plane wave matrix elements
!       pol%isrtx   orders the |G(i)|^2   i=1,pol%nmtx
!       vwfn%isort  orders |qk+g|^2    (in vwfn type)
!
!       energies are apparently assumed in Rydbergs.
!
! after get_from_fftbox(), we  get < c,k | exp(i(q+G).r) | v,(k-q) >

!> (default) pol%k_plus_q = F ==> we calculate mdat = < c, k   | e^{i(q+G).r} | v, k-q >
!>           pol%k_plus_q = T ==> we calculate mdat = < v, k+q | e^{i(q+G).r} | c, k   >
!===========================================================================

module mtxel_m
  use global_m
  use fftw_m
  use misc_m
  implicit none
  private
  public :: mtxel
contains
  subroutine mtxel(gvec, pol, ispin, wfn_bra, wfn_ket, mdat_complex, mdat_real)
    type (gspace), intent(in) :: gvec
    type (polarizability), intent(in) :: pol
    integer, intent(in) :: ispin
    type (wavefunction), intent(in) :: wfn_ket
    type (wavefunction), intent(in) :: wfn_bra
    complex(DPC), optional, dimension(:,:,:,:), intent(out) :: mdat_complex
    real(DPC), optional, dimension(:,:,:,:,:), intent(out) :: mdat_real

    integer :: ic, iv, jsp, ig !, iv_abs, ic_abs
    integer, dimension(3) :: Nfft
    complex(DPC), dimension(:,:,:), allocatable :: fftbox1,fftbox2
    SCALAR, dimension(:), allocatable :: tmparray
    real(DP) :: scale
    PUSH_SUB(mtxel)

    if (present(mdat_complex) .and. present(mdat_real)) then
       call die("Only one mdat is allowed.", only_root_writes=.true.)
    elseif ( (.not. present(mdat_complex)) .and. (.not. present(mdat_real))) then
       call die("Must have one mdat", only_root_writes=.true.)
    endif

    if(peinf%inode .eq. 0) call timacc(26,1)
    !> Use FFTs to calculate matrix elements
    !> Compute size of FFT box we need
    ! write(*,*) "pol%FFTgrid = ", pol%FFTgrid    
    call setup_FFT_sizes(pol%FFTgrid, Nfft, scale)
    !> Allocate FFT boxes
    SAFE_ALLOCATE(fftbox2, (Nfft(1),Nfft(2),Nfft(3)))
    SAFE_ALLOCATE(fftbox1, (Nfft(1),Nfft(2),Nfft(3)))
    SAFE_ALLOCATE(tmparray, (pol%nmtx))

    !> < n^c k | Exp[ I (q+G1).x1] |n^v (k-q)> = < n^c k up | Exp[I (q+G1).x1] |n^v (k-q) up > + < n^c k down | Exp[I (q+G1).x1] |n^v (k-q) down >
    do jsp = ispin, ispin * wfn_bra%nspinor
       !> put | v,k-q,vspinor > into fftbox1
       !> Common/fftw.p.f:
       !> subroutine zput_into_fftbox(ndata, data, glist, gindex, fftbox, Nfft)
       do iv = 1, wfn_ket%nband
          !> in BSE/input_mq.f90: iv = iv_rela = kpq%ifmax(irkq,is) - iv_absolute +1
          call put_into_fftbox(wfn_ket%ng, wfn_ket%cg(:,iv,jsp), gvec%components, wfn_ket%isort, fftbox1, Nfft)
          !> Common/fftw.p.f:
          !> subroutine do_FFT(fftbox, Nfft, sign)
          !> sign = 1 ==> FFTW_BACKWARD
          !>      = 2 ==> FFTW_FORWARD
          call do_FFT(fftbox1, Nfft, 1)

          if (pol%k_plus_q) then
             !> < v, k+q |
             call conjg_fftbox(fftbox1, Nfft)
          endif

          ! Now we loop over the conduction states and get the matrix elements:
          ! 1. Get conduction wave function and put it into box 2,
          ! 2. do FFT, get u_{ck}(r)
          ! 3. multiply by box1 contents, get F(r) = [u_{vk+q)(r)]^* u_{ck}(r)
          ! 4. do FFT again, and extract the resulting matrix elements and put the into pol
          ! We conjugate the final result since we really want <ck|e^{-i(q+G).r}|vk+q>
          ! but we have calculated <vk+q|e^{i(q+G).r}|ck>.
          !< ic_loc is relative local index of conduction bands
          do ic = 1, wfn_bra%nband
             !> (default) if pol%skip_nvb == 0 and pol%skip_ncb == 0, do not skip any transition
             !> skip the transitions between all the valence bands to the lower part of conduction bands
             if ((pol%skip_nvb .eq. 0) .and. (pol%skip_ncb .ne. 0)) then
                if (ic .le. pol%skip_ncb) then
                   cycle
                endif
                !> skip the transitions between the upper part of valence bands to all the conduction bands
             elseif ((pol%skip_nvb .ne. 0) .and. (pol%skip_ncb .eq. 0)) then
                if (iv .le. pol%skip_nvb) then
                   cycle
                endif
                !> skip the transitions between the upper part of valence bands to the lower part of conduction bands
             elseif ((pol%skip_nvb .ne. 0) .and. (pol%skip_ncb .ne. 0)) then
                if ((iv .le. pol%skip_nvb) .and. (ic .le. pol%skip_ncb)) then
                   cycle
                endif
             endif

             call put_into_fftbox(wfn_bra%ng, wfn_bra%cg(:,ic,jsp), gvec%components, wfn_bra%isort, fftbox2, Nfft)
             call do_FFT(fftbox2, Nfft, 1)
             
             if (.not. pol%k_plus_q) then
                !< We need the complex conjugate of u_{ck)(r) for the cross correlation
                !< c,k,is |
                call conjg_fftbox(fftbox2, Nfft)
             endif
             
             !> Multiply contents of two fft boxes, result into fftbox2
             !> fftbox1 is not updated
             call multiply_fftboxes(fftbox1, fftbox2, Nfft)
             call do_FFT(fftbox2, Nfft, 1)
             !> we already divide the output by Nfft, scale = 1/Nfft
             !> subroutine dget_from_fftbox(ndata, data, glist, gindex, fftbox, Nfft, scale)
             !> data(j) = fftbox(bidx(1),bidx(2),bidx(3))*scale
             !> subroutine gvec_to_fft_index(g,idx,Nfft)
             !> pol%isrtx
             !> tmparray contains the matrix elements
             call get_from_fftbox(pol%nmtx, tmparray, gvec%components, pol%isrtx, fftbox2, Nfft, scale)

             !> This accumulation summation is for nspinor = 2 and jsp = 2 case!
             if (present(mdat_complex)) then
                do ig = 1, pol%nmtx
                   mdat_complex(ig, ispin, iv, ic) = mdat_complex(ig, ispin, iv, ic) + tmparray(ig)
                enddo
             else
                do ig = 1, pol%nmtx
                   mdat_real(1, ig, ispin, iv, ic) = mdat_real(1, ig, ispin, iv, ic) + DBLE(tmparray(ig))
                   mdat_real(2, ig, ispin, iv, ic) = mdat_real(2, ig, ispin, iv, ic) + DIMAG(tmparray(ig))
                enddo
             endif
          enddo ! ic
       enddo ! iv
    enddo ! jsp loop

    SAFE_DEALLOCATE(tmparray)
    SAFE_DEALLOCATE(fftbox1)
    SAFE_DEALLOCATE(fftbox2)

    if(peinf%inode.eq.0) call timacc(26,2)

    POP_SUB(mtxel)
    return
  end subroutine mtxel
end module mtxel_m

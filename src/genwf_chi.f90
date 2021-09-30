#include "f_defs.h"

!!===================================================================
!!
!! Module:
!!
!!    genwf_chi()          by Meng Wu (2020)
!!
!!    Generate valence wavefunctions at k \pm q (wfnv) and conduction
!!    wavefunctions at k (wfnc) with corresponding wavefunctions (intwfnc,
!!    intwfnvq) in the reduced Brillouin zone using crystal point group
!!    symmetry operations.
!!
!!===================================================================

subroutine genwf_chi(pol, ik, iq, crys, gvec, kg, kgq, syms, symsq, wfnc, wfnv, intwfnc, intwfnvq)
  use global_m
  use blas_m
  use gmap_m
  use input_utils_m
  use sort_m
  use misc_m
  use so32su2_m
  implicit none

  type (polarizability), intent(in) :: pol
  integer, intent(in) :: ik, iq
  type (crystal), intent(in) :: crys
  type (gspace), intent(in) :: gvec
  type (grid), intent(in) :: kg, kgq
  type (symmetry), intent(in) :: syms, symsq
  !! We must set intent(inout) here, otherwise wfnv%nband will be 0 when iq >= pol%nq0
  type (wavefunction), intent(inout) :: wfnc, wfnv
  type (int_wavefunction), intent(in) :: intwfnc, intwfnvq
  
  integer :: ig, ib, is, ispinor, irk_loc, irkq_loc
  integer, allocatable :: ind(:), isort_new2gvec(:), isort_old2gvec(:), isorti_gvec2old(:)
  real(DP), allocatable :: ekin(:)
  !! cg_(ig, ib, is)  
  SCALAR, allocatable :: cg_(:,:,:), cg_2d(:,:), cg_2d_(:,:)
  SCALAR, allocatable :: ph(:)
  complex(DPC) :: umtrx(2,2), umtrx_transpose(2,2)
  PUSH_SUB(genwf_chi)
  
  irk_loc  =  peinf%irk_g2l(kg%indr(ik),              peinf%inode+1)
  irkq_loc = peinf%irkq_g2l(kgq%indr(pol%indexq(ik)), peinf%inode+1)

  SAFE_ALLOCATE(isort_new2gvec, (gvec%ng))
  SAFE_ALLOCATE(isort_old2gvec, (gvec%ng))
  SAFE_ALLOCATE(isorti_gvec2old, (gvec%ng))

  if (iq <= pol%nq0) then
     !! Rotate valence states

     !! Since we use a rotation here R[rk]=fk, so |rk|=|fk|
     !! In this way, psi_{rk} and psi_{fk} should have the same number of Gvectors
     wfnv%ng      = intwfnvq%ng(irkq_loc)
     wfnv%nband   = pol%nvb
     wfnv%nspin   = intwfnvq%nspin
     wfnv%nspinor = intwfnvq%nspinor

     !! cg_(ig, ib, is)
     SAFE_ALLOCATE(cg_, (wfnv%ng, wfnv%nband, wfnv%nspin*wfnv%nspinor))
     SAFE_ALLOCATE(ind, (wfnv%ng))
     SAFE_ALLOCATE(ph,  (wfnv%ng))

     isort_old2gvec(:) = intwfnvq%isort(:, irkq_loc)
     !! isorti_gvec2old(ig_gvec) = ig_old
     isorti_gvec2old = 0
     !! Some elements in isorti_gvec2old are 0!
     !! Range of isorti_gvec2old : [0, 1, ..., wfnv%ng]
     !! ind(:) has the same range of isorti_gvec2old
     do ig = 1, wfnv%ng
        isorti_gvec2old(isort_old2gvec(ig)) = ig
     enddo

     SAFE_ALLOCATE(ekin, (gvec%ng))
     call kinetic_energies(gvec, crys%bdot, ekin, qvec = kgq%f(1:3, pol%indexq(ik)))
     call sortrx(gvec%ng, ekin, isort_new2gvec, gvec = gvec%components)
     SAFE_DEALLOCATE(ekin)

     SAFE_ALLOCATE(wfnv%isort, (gvec%ng))
     wfnv%isort(:) = isort_new2gvec(:)

     !! Find ind and ph relating wavefunctions in fk to rk-kpoint
     ind = 0
     ph  = ZERO
     
     call gmap(gvec, symsq, wfnv%ng, kgq%itran(pol%indexq(ik)), &
          kgq%kg0(:,pol%indexq(ik)), isort_new2gvec, isorti_gvec2old, ind, ph)

     !! intwfnvq%cgk(ig, ib, is, irk)
     cg_(1:wfnv%ng, 1:wfnv%nband, 1:wfnv%nspin*wfnv%nspinor) &
          = intwfnvq%cgk(1:wfnv%ng, 1:wfnv%nband, 1:wfnv%nspin*wfnv%nspinor, irkq_loc)

     SAFE_ALLOCATE(wfnv%cg, (wfnv%ng, wfnv%nband, wfnv%nspin*wfnv%nspinor))
     !$OMP PARALLEL DO collapse(4)
     do ib = 1, wfnv%nband
        do is = 1, wfnv%nspin
           do ispinor = 1, wfnv%nspinor
              do ig = 1, wfnv%ng
                 if (ind(ig) > 0) then
                    wfnv%cg(ig, ib, is*ispinor) = ph(ig) * cg_(ind(ig), ib, is*ispinor)
                 else
                    wfnv%cg(ig, ib, is*ispinor) = ZERO
                 endif
              enddo ! ig
           enddo ! ispinor
        enddo ! is
     enddo ! ib
     !$OMP END PARALLEL DO
     
     SAFE_DEALLOCATE(cg_)
     SAFE_DEALLOCATE(ind)
     SAFE_DEALLOCATE(ph)

#ifdef CPLX
     if (wfnv%nspinor .eq. 2) then
        SAFE_ALLOCATE(cg_2d,  (wfnv%ng*wfnv%nband, wfnv%nspin*wfnv%nspinor))
        SAFE_ALLOCATE(cg_2d_, (wfnv%ng*wfnv%nband, wfnv%nspin*wfnv%nspinor))

        call so32su2(symsq%mtrx_cart(1:3, 1:3, kgq%itran(pol%indexq(ik))), umtrx)
        umtrx_transpose = TRANSPOSE(umtrx)
        !$OMP PARALLEL WORKSHARE
        cg_2d   = RESHAPE(wfnv%cg, (/wfnv%ng*wfnv%nband, wfnv%nspin*wfnv%nspinor/))
        cg_2d_  = MATMUL(cg_2d, umtrx_transpose)
        wfnv%cg = RESHAPE(cg_2d_, (/wfnv%ng, wfnv%nband, wfnv%nspin*wfnv%nspinor/))
        !$OMP END PARALLEL WORKSHARE

        SAFE_DEALLOCATE(cg_2d)
        SAFE_DEALLOCATE(cg_2d_)
     endif
#endif

     !! Check norm after rotation
     do ib = 1, wfnv%nband
        call checknorm("wfnvq", ib, pol%indexq(ik), wfnv%nspin, wfnv%cg(1:wfnv%ng, ib, :))
     enddo
  endif !! iq <= pol%nq0

  !! Rotate conduction wavefunctions
  wfnc%ng      = intwfnc%ng(irk_loc)
  wfnc%nband   = pol%ncb
  wfnc%nspin   = intwfnc%nspin
  wfnc%nspinor = intwfnc%nspinor

  SAFE_ALLOCATE(cg_, (wfnc%ng, wfnc%nband, wfnc%nspin*wfnc%nspinor))
  SAFE_ALLOCATE(ind, (wfnc%ng))
  SAFE_ALLOCATE(ph,  (wfnc%ng))

  isort_old2gvec(:) = intwfnc%isort(:, irk_loc)
  isorti_gvec2old = 0
  do ig = 1, wfnc%ng
     isorti_gvec2old(isort_old2gvec(ig)) = ig
  enddo

  SAFE_ALLOCATE(ekin, (gvec%ng))
  call kinetic_energies(gvec, crys%bdot, ekin, qvec = kg%f(:, ik))
  call sortrx(gvec%ng, ekin, isort_new2gvec, gvec = gvec%components)
  SAFE_DEALLOCATE(ekin)

  SAFE_ALLOCATE(wfnc%isort, (gvec%ng))
  wfnc%isort(:) = isort_new2gvec(:)

  ind = 0
  ph  = ZERO
  call gmap(gvec, syms, wfnc%ng, kg%itran(ik), kg%kg0(:,ik), isort_new2gvec, &
       isorti_gvec2old, ind, ph)

  cg_(1:wfnc%ng, 1:wfnc%nband, 1:wfnc%nspin*wfnc%nspinor) &
       = intwfnc%cgk(1:wfnc%ng, 1:wfnc%nband, 1:wfnc%nspin*wfnc%nspinor, irk_loc)

  SAFE_ALLOCATE(wfnc%cg, (wfnc%ng, wfnc%nband, wfnc%nspin*wfnc%nspinor))
  !$OMP PARALLEL DO collapse(4)
  do ib = 1, wfnc%nband
     do is = 1, wfnc%nspin
        do ispinor = 1, wfnc%nspinor
           do ig = 1, wfnc%ng
              if (ind(ig) > 0) then
                 wfnc%cg(ig, ib, is*ispinor) = ph(ig) * cg_(ind(ig), ib, is*ispinor)
              else
                 wfnc%cg(ig, ib, is*ispinor) = ZERO
              endif
           enddo
        enddo ! ispinor
     enddo ! is
  enddo ! ib
  !$OMP END PARALLEL DO

  SAFE_DEALLOCATE(cg_)
  SAFE_DEALLOCATE(ind)
  SAFE_DEALLOCATE(ph)

#ifdef CPLX
  if (wfnc%nspinor .eq. 2) then
     SAFE_ALLOCATE(cg_2d,  (wfnc%ng*wfnc%nband, wfnc%nspin*wfnc%nspinor))
     SAFE_ALLOCATE(cg_2d_, (wfnc%ng*wfnc%nband, wfnc%nspin*wfnc%nspinor))
     call so32su2(syms%mtrx_cart(1:3, 1:3, kg%itran(ik)), umtrx)
     umtrx_transpose = TRANSPOSE(umtrx)

     !$OMP PARALLEL WORKSHARE
     cg_2d   = RESHAPE(wfnc%cg, (/wfnc%ng*wfnc%nband, wfnc%nspin*wfnc%nspinor/))
     cg_2d_  = MATMUL(cg_2d, umtrx_transpose)
     wfnc%cg = RESHAPE(cg_2d_, (/wfnc%ng, wfnc%nband, wfnc%nspin*wfnc%nspinor/))
     !$OMP END PARALLEL WORKSHARE

     SAFE_DEALLOCATE(cg_2d)
     SAFE_DEALLOCATE(cg_2d_)
  endif
#endif

  do ib = 1, wfnc%nband
     call checknorm('wfnc', ib, ik, wfnc%nspin, wfnc%cg(1:wfnc%ng, ib, :))
  enddo

  SAFE_DEALLOCATE(isort_old2gvec)
  SAFE_DEALLOCATE(isort_new2gvec)
  SAFE_DEALLOCATE(isorti_gvec2old)

  POP_SUB(genwf_chi)
  return
end subroutine genwf_chi

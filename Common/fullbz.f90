#include "f_defs.h"

!--------------------------------------------------------------------
!
! module fullbz_m
!
! Routines: fullbz
!
! If wigner_seitz = .true. (for BSE, PlotXct, NonLinearOptics)
!     Uses a Wigner-Seitz construction to define the Brillouin zone.
! If wigner_seitz = .false. (for Epsilon, Sigma)
!     Uses the usual "box" BZ.
!
!     input: crys,syms type
!            ntran (number of symmetry operations)
!            gr%nr
!            gr%rk
!
!     output: gr type (except gr%nr, gr%rk)
!
!--------------------------------------------------------------------

!> after fullbz() we will have
!> MATMUL(DBLE(syms%mtrx_reci(1:3,1:3,gr%itran(ifk))),gr%r(1:3,gr%indr(ifk))) + gr%kg0(1:3,ifk) = gr%f(1:3,ifk)
module fullbz_m
  use global_m
  use misc_m
  implicit none
  private
  public :: fullbz, dealloc_grid, dealloc_grid_simple, bz_r2f_copy
contains

  subroutine dealloc_grid(gr)
    type(grid), intent(inout) :: gr

    PUSH_SUB(dealloc_grid)

    SAFE_DEALLOCATE_P(gr%r)
    SAFE_DEALLOCATE_P(gr%f)
    SAFE_DEALLOCATE_P(gr%indr)
    SAFE_DEALLOCATE_P(gr%itran)
    SAFE_DEALLOCATE_P(gr%kg0)
    SAFE_DEALLOCATE_P(gr%dege_r)
    SAFE_DEALLOCATE_P(gr%findex)

    POP_SUB(dealloc_grid)
    return
  end subroutine dealloc_grid

  subroutine dealloc_grid_simple(gr)
    type(grid), intent(inout) :: gr

    PUSH_SUB(dealloc_grid_simple)

    SAFE_DEALLOCATE_P(gr%r)
    SAFE_DEALLOCATE_P(gr%f)
    SAFE_DEALLOCATE_P(gr%indr)
    SAFE_DEALLOCATE_P(gr%itran)
    SAFE_DEALLOCATE_P(gr%kg0)

    POP_SUB(dealloc_grid_simple)
    return
  end subroutine dealloc_grid_simple

  !> We use WS BZ in fullbz
  subroutine fullbz(crys,syms,gr,use_identity_only,nfix)
    type (crystal), intent(in) :: crys
    type (symmetry), intent(in) :: syms
    type (grid), intent(inout) :: gr
    logical, optional, intent(in) :: use_identity_only
    integer, optional, intent(in) :: nfix !< don`t unfold the first nfix points
    real(DP) :: fk_temp(3), rk_temp(3)
    real(DP), allocatable :: fk_list(:,:)
    integer, allocatable :: indr_(:), itran_(:), g0_(:,:)
    integer :: gumk(3), ntran_useful, ntran_, it,irk,ifk
    logical :: found
    PUSH_SUB(fullbz)

    if (present(use_identity_only)) then
       if (use_identity_only) then
          ntran_useful = 1
       else
          ntran_useful = syms%ntran
       endif
    else
       ntran_useful = syms%ntran
    endif

    SAFE_ALLOCATE(fk_list, (3, gr%nr * ntran_useful))
    SAFE_ALLOCATE(indr_, (gr%nr * ntran_useful))
    SAFE_ALLOCATE(itran_, (gr%nr * ntran_useful))
    SAFE_ALLOCATE(g0_, (3,gr%nr * ntran_useful))
    fk_list(:,:) = 0.0D0
    indr_(:) = 0
    itran_(:) = 0
    g0_(:,:) = 0
    gr%nf=0

    !> Make sure that all the RBZ points gr%r are in the WS BZ
    do irk = 1, gr%nr
       rk_temp(:) = gr%r(:, irk)
       call get_gumk3(crys%bdot, rk_temp, gumk)
       rk_temp(:) = rk_temp(:) - DBLE(gumk(:))
       if (NORM2(rk_temp(:) - gr%r(:,irk)) > TOL_ZERO) then
          if (peinf%inode .eq. 0) then
             write(*,'(1X,A)') "kpoints in wave function outside WS BZ."
          endif
       endif
    enddo

    do irk = 1, gr%nr
       ntran_ = ntran_useful
       if (present(nfix)) then
          if (irk <= nfix) ntran_ = 1
       endif

       sym_loop : do it = 1, ntran_
          fk_temp(:) = MATMUL(DBLE(syms%mtrx_reci(:,:,it)), gr%r(:,irk))
          call get_gumk3(crys%bdot, fk_temp, gumk)

          !> Move to WS BZ
          fk_temp(:) = fk_temp(:) - DBLE(gumk(:))

          found = .false.
          do ifk = 1, gr%nf
             if (NORM2(fk_temp(1:3)-fk_list(1:3, ifk)) .lt. TOL_SMALL) then
                found = .true.
                cycle sym_loop
             endif
          enddo
          
          if (.not. found) then
             ! write(*,'(A,3E20.12,A,3E20.12,A,ES30.23)') "rk = ", gr%r(:,irk), " fk = ", fk_temp(:), " len(fk) = ", DOT_PRODUCT(fk_temp-1.0D-8, MATMUL(crys%bdot, fk_temp-1.0D-8))
             ! write(*,'(A,3E20.12,A,3E20.12,A,ES30.23)') "rk = ", gr%r(:,irk), " fk = ", fk_temp(:), " len(fk) = ", DOT_PRODUCT(fk_temp, MATMUL(crys%bdot, fk_temp))             
             gr%nf = gr%nf + 1
             if (gr%nf > gr%nr * ntran_useful) then
                call die('fullbz internal error: gr%nf > gr%nr * syms%ntran')
             endif

             fk_list(1:3,gr%nf) = fk_temp(1:3)
             indr_(gr%nf) = irk
             itran_(gr%nf) = it
             g0_(:,gr%nf) = - gumk(:)
          else
             write(*,'(A)') "found"
          endif
       enddo sym_loop !end loop over symmetries
    enddo !end loop over the q-points from epsmat and eps0mat

    SAFE_ALLOCATE(gr%f, (3,gr%nf))
    SAFE_ALLOCATE(gr%indr, (gr%nf))
    SAFE_ALLOCATE(gr%itran, (gr%nf))
    SAFE_ALLOCATE(gr%kg0, (3,gr%nf))
    gr%f(:,:) = 0.0D0
    gr%indr(:) = 0
    gr%itran(:) = 0
    gr%kg0(:,:) = 0

    gr%f(1:3,1:gr%nf) = fk_list(1:3,1:gr%nf)
    gr%indr(1:gr%nf) = indr_(1:gr%nf)
    gr%itran(1:gr%nf) = itran_(1:gr%nf)
    gr%kg0(:,1:gr%nf) = g0_(:,1:gr%nf)

    SAFE_DEALLOCATE(fk_list)
    SAFE_DEALLOCATE(indr_)
    SAFE_DEALLOCATE(itran_)
    SAFE_DEALLOCATE(g0_)
    POP_SUB(fullbz)

    return
  end subroutine fullbz

  !> We use WS BZ in fullbz
  subroutine bz_r2f_copy(gr)
    type (grid), intent(inout) :: gr
    integer :: ik
    
    PUSH_SUB(bz_r2f_copy)

    gr%nf = gr%nr
    SAFE_ALLOCATE(gr%kg0,(3,gr%nf))
    SAFE_ALLOCATE(gr%f,(3,gr%nf))
    SAFE_ALLOCATE(gr%itran,(gr%nf))
    SAFE_ALLOCATE(gr%indr,(gr%nf))
    gr%kg0(:,:)=0
    gr%itran(:)=1
    gr%f(:,:)=gr%r(:,:)
    do ik = 1, gr%nf
       gr%indr(ik)=ik
    enddo

    POP_SUB(bz_r2f_copy)
    return
  end subroutine bz_r2f_copy

end module fullbz_m

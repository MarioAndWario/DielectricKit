#include "f_defs.h"

module irrbz_m
  use global_m
  use misc_m
  use fullbz_m
  implicit none
  private
  public :: irrbz_2, irrbz_1
contains

  !> We will use symmetry in the little group
  !> and FBZ q points in gr%f(:,1:gr%nf)
  !> to construct RBZ q points in gr_little_group%r(:, 1:gr_little_group%nr)

  subroutine irrbz_2(crys, syms, syms_little_group, gr, gr_little_group, nfix, gr_little_group_to_gr)
    type (crystal), intent(in) :: crys
    type (symmetry), intent(in) :: syms, syms_little_group
    type(grid), intent(in) :: gr
    type(grid), intent(out) :: gr_little_group
    integer, optional, intent(in) :: nfix
    type(grid_relation), optional, intent(out) :: gr_little_group_to_gr
    real(DP), allocatable :: r_(:,:)
    integer, allocatable :: dege_r_(:)
    integer :: irk,ifk,isym,nfix_, ntran_, irk_little_group
    real(DP) :: k_temp(3)
    integer :: gumk(3)
    logical :: found
    PUSH_SUB(irrbz_2)

    !> Initialize gr_little_group
    gr_little_group%nf = 0    
    gr_little_group%nr = 0

    !> Temporary variables
    SAFE_ALLOCATE(r_, (3, gr%nf))
    SAFE_ALLOCATE(dege_r_, (gr%nf))
    r_(:,:) = 0.0D0
    dege_r_(:) = 0
    nfix_ = 1
    if (present(nfix)) nfix_ = nfix

    !> initialize number of points in irr. BZ
    gr_little_group%nr = 0
    ik_loop: do ifk = 1, gr%nf
       if (ifk .le. nfix_) then
          ntran_ = 1
       else
          ntran_ = syms_little_group%ntran
       endif
       do isym = 1, ntran_
          !> qk(:) = matmul(dble(syms%mtrx(:,:,syms%indsub(it))),fk(:,ik))
          k_temp(:) = MATMUL(DBLE(syms_little_group%mtrx_reci(:,:,isym)),gr%f(:,ifk))
          call get_gumk3(crys%bdot, k_temp, gumk)
          k_temp(:) = k_temp(:) - DBLE(gumk(:))

          !> compare to other k-points in the irr. BZ with respect to qvec
          do irk = 1, gr_little_group%nr
             if (NORM2(k_temp(1:3) - r_(1:3,irk)) .lt. TOL_SMALL) then
                dege_r_(irk) = dege_r_(irk) + 1
                cycle ik_loop
             endif
          enddo
       enddo ! loop over ntranq tranformation

       gr_little_group%nr = gr_little_group%nr + 1
       r_(:,gr_little_group%nr) = gr%f(:,ifk)
       dege_r_(gr_little_group%nr) = 1
    enddo ik_loop !end loop over full BZ

    SAFE_ALLOCATE(gr_little_group%r, (3, gr_little_group%nr))
    SAFE_ALLOCATE(gr_little_group%dege_r, (gr_little_group%nr))

    gr_little_group%r(:,1:gr_little_group%nr) = r_(:,1:gr_little_group%nr)
    gr_little_group%dege_r(1:gr_little_group%nr) = dege_r_(1:gr_little_group%nr)

    SAFE_DEALLOCATE(r_)
    SAFE_DEALLOCATE(dege_r_)

    !> Find the relation between gr_little_group to gr
    !> MATMUL(DBLE(syms%mtrx_reci(:,:,gr_little_group_to_gr%symsindex(irk_little_group))),gr%r(:,irk)) + gr_little_group_to_gr%gumk(:,irk_little_group) = gr_little_group%r(:,irk_little_group)
    !> irk = gr_little_group_to_gr%rindex(irk_little_group)
    if (present(gr_little_group_to_gr)) then
       SAFE_ALLOCATE(gr_little_group_to_gr%rindex, (gr_little_group%nr))
       SAFE_ALLOCATE(gr_little_group_to_gr%symsindex, (gr_little_group%nr))
       SAFE_ALLOCATE(gr_little_group_to_gr%gumk, (3,gr_little_group%nr))

       little_group_loop: do irk_little_group = 1, gr_little_group%nr
          found = .false.
          do irk = 1, gr%nr
             do isym = 1, syms%ntran
                k_temp(:) = MATMUL(DBLE(syms%mtrx_reci(:,:,isym)),gr%r(:,irk))
                call get_gumk3(crys%bdot, k_temp, gumk)
                k_temp(:) = k_temp(:) - DBLE(gumk(:))
                if (NORM2(k_temp(:) - gr_little_group%r(:,irk_little_group)) < TOL_SMALL ) then
                   found = .true.
                   gr_little_group_to_gr%rindex(irk_little_group) = irk
                   gr_little_group_to_gr%symsindex(irk_little_group) = isym
                   gr_little_group_to_gr%gumk(:,irk_little_group) = - gumk(:)
                   cycle little_group_loop
                endif
             enddo
          enddo
          if (.not. found) then
             call die("irrbz_2: Cannot find corresponding coarse point", only_root_writes=.true.)
          endif
       enddo little_group_loop
    endif

    POP_SUB(irrbz_2)

    return
  end subroutine irrbz_2

  !> For epsilon
  subroutine irrbz_1(crys, syms_little_group, gr, gr_little_group)
    type (crystal), intent(in) :: crys
    type (symmetry), intent(in) :: syms_little_group
    type(grid), intent(in) :: gr
    type(grid), intent(out) :: gr_little_group
    real(DP), allocatable :: r_(:,:)
    integer, allocatable :: dege_r_(:), findex_(:)
    integer :: irk, ifk, isym, ntran_, gumk(3)
    real(DP) :: k_temp(3)
    PUSH_SUB(irrbz_1)

    !> Initialize gr_little_group
    gr_little_group%nf = gr%nf
    gr_little_group%nr = 0

    !> Temporary variables
    SAFE_ALLOCATE(r_, (3, gr%nf))
    SAFE_ALLOCATE(dege_r_, (gr%nf))
    SAFE_ALLOCATE(findex_, (gr%nf))

    r_(:,:) = 0.0D0
    dege_r_(:) = 0
    findex_(:) = 0

    !> initialize number of points in irr. BZ
    gr_little_group%nr = 0
    ntran_ = syms_little_group%ntran

    ik_loop: do ifk = 1, gr%nf
       do isym = 1, ntran_
          !> qk(:) = matmul(dble(syms%mtrx(:,:,syms%indsub(it))),fk(:,ik))
          k_temp(:) = MATMUL(DBLE(syms_little_group%mtrx_reci(:,:,isym)), gr%f(:,ifk))
          call get_gumk3(crys%bdot, k_temp, gumk)
          k_temp(:) = k_temp(:) - DBLE(gumk(:))
          !> compare to other k-points in the irr. BZ with respect to qvec
          do irk = 1, gr_little_group%nr
             if (NORM2(k_temp(1:3) - r_(1:3,irk)) .lt. TOL_SMALL) then
                dege_r_(irk) = dege_r_(irk) + 1
                cycle ik_loop
             endif
          enddo
       enddo ! loop over ntranq tranformation
       gr_little_group%nr = gr_little_group%nr + 1
       r_(:,gr_little_group%nr) = gr%f(:,ifk)
       dege_r_(gr_little_group%nr) = 1
       !> findex_(irk_little_group) = ifk
       findex_(gr_little_group%nr) = ifk
    enddo ik_loop !end loop over full BZ

    SAFE_ALLOCATE(gr_little_group%r, (3, gr_little_group%nr))
    SAFE_ALLOCATE(gr_little_group%dege_r, (gr_little_group%nr))
    SAFE_ALLOCATE(gr_little_group%findex, (gr_little_group%nr))

    gr_little_group%r(:,1:gr_little_group%nr) = r_(:,1:gr_little_group%nr)
    gr_little_group%dege_r(1:gr_little_group%nr) = dege_r_(1:gr_little_group%nr)
    gr_little_group%findex(1:gr_little_group%nr) = findex_(1:gr_little_group%nr)

    SAFE_DEALLOCATE(r_)
    SAFE_DEALLOCATE(dege_r_)
    SAFE_DEALLOCATE(findex_)

    POP_SUB(irrbz_1)
    return
  end subroutine irrbz_1
end module irrbz_m

#include "f_defs.h"

!==================================================================
!
! Routines:
!
! 1. subgrp()   Originally By ?                 Last Modified 03/01/2020 (MW)
!
!    Determines a subgroup of the symmetry group that preserves a q-vector.
!
!==================================================================

!> MATMUL(syms_little_group%mtrx_reci(1:3,1:3,isym),kshift(1:3)) - syms_little_group%kgzero(1:3,i) = kshift(1:3), i=1,syms_little_group%ntran
!> Store little group in syms_little_group, not syms
subroutine subgrp_2(crys, syms, kshift, syms_little_group)
  use global_m
  use misc_m
  implicit none
  type (crystal), intent(in) :: crys
  type (symmetry), intent(in) :: syms
  real(DP), intent(in) :: kshift(3)
  type (symmetry), intent(out) :: syms_little_group
  integer ::  isym, gumk(3)
  real(DP) :: k_temp(3)
  PUSH_SUB(subgrp_2)

  syms_little_group%ntran = 0
  syms_little_group%mtrx(:,:,:) = 0
  syms_little_group%mtrx_reci(:,:,:) = 0
  syms_little_group%mtrx_cart(:,:,:) = 0.0D0
  syms_little_group%tnp(:,:) = 0.0D0
  syms_little_group%cell_symmetry = syms%cell_symmetry
  syms_little_group%kgzero(:,:) = 0

  do isym = 1, syms%ntran
     k_temp(1:3) = MATMUL(DBLE(syms%mtrx_reci(1:3, 1:3, isym)), kshift(1:3)) - kshift(1:3)
     call get_gumk3(crys%bdot, k_temp, gumk)
     k_temp(:) = k_temp(:) - DBLE(gumk(:))
     if (NORM2(k_temp(1:3)) .lt. TOL_SMALL) then
        syms_little_group%ntran = syms_little_group%ntran + 1
        syms_little_group%mtrx(:,:,syms_little_group%ntran) = syms%mtrx(:,:,isym)
        syms_little_group%mtrx_reci(:,:,syms_little_group%ntran) = syms%mtrx_reci(:,:,isym)
        syms_little_group%mtrx_cart(:,:,syms_little_group%ntran) = syms%mtrx_cart(:,:,isym)
        syms_little_group%tnp(:,syms_little_group%ntran) = syms%tnp(:,isym)
        syms_little_group%kgzero(:,syms_little_group%ntran) = gumk(:)
     endif
  enddo

  POP_SUB(subgrp_2)

  return
end subroutine subgrp_2

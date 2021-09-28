#include "f_defs.h"

!=====================================================================
!
! Routines:
!
! (1) distrib_chi()       Originally By MLT    Last Modified by Meng Wu (2020)
!
!      nv = # of valence bands
!      nc = # of conduction bands
!      nk = # of k-points
!      np = # of MPI processes
!
!=====================================================================

!! If WFNmq is not used, pol%indexq is identity mapping
subroutine distrib_chi(pol, kg, kgq)
  use global_m
  use misc_m
  use scalapack_m
  implicit none

  type (polarizability), intent(in) :: pol
  type (grid), intent(in) :: kg, kgq

  integer :: ik, ik_loc, nmax, ipes
  real(DP) :: mem, rmem, rmem2, rmemtemp1, rmemtemp2
  integer, allocatable :: irk_(:), irkq_(:)
  integer, allocatable :: irk_unique_global(:,:), irkq_unique_global(:,:)
  integer :: irk, irkq, irk_loc, irkq_loc
  PUSH_SUB(distrib_chi)

  !! 1D-Block distribution of fk on peinf%npes procs
  !! Number of fk on 1st proc (also the maximum among all procs)
  peinf%nkpe = iceil(kg%nf, peinf%npes)
  !! Number of fk owned by each proc
  SAFE_ALLOCATE(peinf%ikt, (peinf%npes))
  SAFE_ALLOCATE(peinf%ik,  (peinf%npes, peinf%nkpe))
  peinf%ikt = 0
  peinf%ik = 0
  do ipes = 1, peinf%npes
     peinf%ikt(ipes) = NUMROC(kg%nf, peinf%nkpe, ipes-1, 0, peinf%npes)
     do ik_loc = 1, peinf%ikt(ipes)
        peinf%ik(ipes, ik_loc) = INDXL2G(ik_loc, peinf%nkpe, ipes-1, 0, peinf%npes)
     enddo
  enddo
  nmax = NUMROC(kg%nf, peinf%nkpe, 0, 0, peinf%npes)

  !! Build peinf%nrk, peinf%irk_l2g, peinf%irk_g2l mappings
  SAFE_ALLOCATE(peinf%nrk, (peinf%npes))  !! For intwfnc
  peinf%nrk = 0
  !! irk contains the indices of rk in kg%r(:,1:kg%nr) ==> used by conduction bands
  SAFE_ALLOCATE(irk_, (nmax))
  SAFE_ALLOCATE(irk_unique_global, (nmax,peinf%npes))
  irk_unique_global = 0

  SAFE_ALLOCATE(peinf%nrkq, (peinf%npes)) !! For intwfnvq
  peinf%nrkq = 0
  !! irkq contains the indices of rkq in kgq%r(:,1:kgq%nr) ==> used by valence bands
  SAFE_ALLOCATE(irkq_, (nmax))
  SAFE_ALLOCATE(irkq_unique_global, (nmax,peinf%npes))
  irkq_unique_global = 0

  !! loop over all procs
  do ipes = 1, peinf%npes
     if (peinf%ikt(ipes) .eq. 0) cycle
     !! reset irk and irkq
     irk_ = 0
     irkq_ = 0
     do ik_loc = 1, peinf%ikt(ipes)
        ik = INDXL2G(ik_loc, peinf%nkpe, ipes-1, 0, peinf%npes)
        irk_(ik_loc)  =  kg%indr(ik)
        irkq_(ik_loc) = kgq%indr(pol%indexq(ik))
     enddo

     !! Find (number of unique elements in irk or irkq), this is the number of rk that is to be stored in intwfnc and intwfnvq, respectively!
     call unique(irk_(1:peinf%ikt(ipes)), peinf%nrk(ipes), irk_unique_global(1:peinf%ikt(ipes), ipes))
     call unique(irkq_(1:peinf%ikt(ipes)), peinf%nrkq(ipes), irkq_unique_global(1:peinf%ikt(ipes), ipes))
     call Bubble_Sort(irk_unique_global(1:peinf%nrk(ipes), ipes))
     call Bubble_Sort(irkq_unique_global(1:peinf%nrkq(ipes), ipes))

     ! if (peinf%inode .eq. 0) then
     !    write(*,'(A,I5,A,I5,A,I5,A,I5,A,I5)') "ipes = ", ipes, " nrk = ", peinf%nrk(ipes), " nrkq = ", peinf%nrkq(ipes)
     !    write(*,'(A,20I5)') "irk_ = ", irk_(:)
     !    write(*,'(A,20I5)') "irkq_ = ", irkq_(:)
     !    write(*,'(A,20I5)') "irk_unique_global = ", irk_unique_global(:, ipes)
     !    write(*,'(A,20I5)') "irkq_unique_global = ", irkq_unique_global(:, ipes)
     !    write(*,'(A)') "============================="
     ! endif
  enddo

  !! Initialize irk_l2g, irkq_l2g, irk_g2l, irkq_g2l
  !! irk_l2g(irk_loc) = irk
  !! irkq_l2g(irkq_loc) = irkq
  SAFE_ALLOCATE( peinf%irk_l2g, (MAX( peinf%nrk(peinf%inode+1),1)))
  SAFE_ALLOCATE(peinf%irkq_l2g, (MAX(peinf%nrkq(peinf%inode+1),1)))
  peinf%irk_l2g  = 0
  peinf%irkq_l2g = 0

  !! irk_g2l(irk, ipes) = irk_loc
  !! irkq_g2l(irkq, ipes) = irkq_loc
  SAFE_ALLOCATE( peinf%irk_g2l, ( kg%nr, peinf%npes))
  SAFE_ALLOCATE(peinf%irkq_g2l, (kgq%nr, peinf%npes))
  peinf%irk_g2l  = 0
  peinf%irkq_g2l = 0

  !! Get irk_g2l
  !! if peinf%irk_g2l(irk,ipes) = 0, it means the ipes-th proc does not store irk-th kg%r
  do irk = 1, kg%nr
     do ipes = 1, peinf%npes
        loop_irk_loc: do irk_loc = 1, peinf%nrk(ipes)
           if (irk_unique_global(irk_loc, ipes) .eq. irk) then
              peinf%irk_g2l(irk, ipes) = irk_loc
              !! We can do this exit because there are no duplicate elements in irk_unique_global(:,ipes)
              exit loop_irk_loc
           endif
        enddo loop_irk_loc
     enddo
  enddo

  !! Get irkq_g2l
  do irkq = 1, kgq%nr
     do ipes = 1, peinf%npes
        loop_irkq_loc: do irkq_loc = 1, peinf%nrkq(ipes)
           if (irkq_unique_global(irkq_loc, ipes) .eq. irkq) then
              peinf%irkq_g2l(irkq, ipes) = irkq_loc
              !! We can do this exit because there are no duplicate elements in irkq_unique_global(:,ipes)
              exit loop_irkq_loc
           endif
        enddo loop_irkq_loc
     enddo
  enddo

  !! Get irk_l2g
  do irk_loc = 1, peinf%nrk(peinf%inode+1)
     peinf%irk_l2g(irk_loc) = irk_unique_global(irk_loc,peinf%inode+1)
  enddo
  !! Get irkq_l2g
  do irkq_loc = 1, peinf%nrkq(peinf%inode+1)
     peinf%irkq_l2g(irkq_loc) = irkq_unique_global(irkq_loc,peinf%inode+1)
  enddo

  SAFE_DEALLOCATE(irk_unique_global)
  SAFE_DEALLOCATE(irkq_unique_global)

  POP_SUB(distrib_chi)
  return

contains

  function number_unique_elements(list_int)
    integer, intent(in) :: list_int(:)
    logical :: mask(size(list_int))
    integer :: number_unique_elements, i
    mask = .true.

    do i = size(list_int), 2, -1
       mask(i) = .NOT. (ANY(list_int(:i-1) == list_int(i)))
    enddo
    number_unique_elements = COUNT(mask)

    return
  end function number_unique_elements

  subroutine get_list_unique(list_int, nunique, list_unique)
    integer, intent(in) :: list_int(:)
    integer, intent(in) :: nunique
    integer, intent(out) :: list_unique(nunique)
    logical :: mask(size(list_int))
    integer :: i
    mask = .true.
    do i = size(list_int), 2, -1
       mask(i) = .NOT. (ANY(list_int(:i-1) == list_int(i)))
    enddo
    list_unique = PACK(list_int, mask)
    call Shell_Sort(list_unique)
    return
  end subroutine get_list_unique

  SUBROUTINE Shell_Sort(a)
    IMPLICIT NONE
    INTEGER :: i, j, increment
    INTEGER :: temp
    INTEGER, INTENT(inout) :: a(:)

    increment = SIZE(a) / 2
    DO WHILE (increment > 0)
       DO i = increment+1, SIZE(a)
          j = i
          temp = a(i)
          DO WHILE (j >= increment+1 .AND. a(j-increment) > temp)
             a(j) = a(j-increment)
             j = j - increment
          ENDDO
          a(j) = temp
       ENDDO
       IF (increment == 2) THEN
          increment = 1
       ELSE
          increment = increment * 5 / 11
       ENDIF
    ENDDO
  END SUBROUTINE Shell_Sort

  subroutine unique(list_int, k, res)
    integer, intent(in) :: list_int(:)
    integer, intent(out) :: k
    integer, intent(out) :: res(size(list_int))
    integer :: i, j
    k = 1
    res(1) = list_int(1)
    outer: do i=2,size(list_int)
       do j=1,k
          if (res(j) == list_int(i)) then
             ! Found a match so start looking again
             cycle outer
          end if
       end do
       ! No match found so add it to the output
       k = k + 1
       res(k) = list_int(i)
    end do outer
  end subroutine unique

  SUBROUTINE Bubble_Sort(a)
    integer, INTENT(in out), DIMENSION(:) :: a
    integer :: temp
    INTEGER :: i, j
    LOGICAL :: swapped

    DO j = SIZE(a)-1, 1, -1
       swapped = .FALSE.
       DO i = 1, j
          IF (a(i) > a(i+1)) THEN
             temp = a(i)
             a(i) = a(i+1)
             a(i+1) = temp
             swapped = .TRUE.
          END IF
       END DO
       IF (.NOT. swapped) EXIT
    END DO
  END SUBROUTINE Bubble_Sort

end subroutine distrib_chi

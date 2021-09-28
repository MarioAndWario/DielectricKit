!=============================================================================
!
! Routines:
!
! (1) distrib()         Originally By MLT       Last Modified 6/2008 FJR
!
!     Distributes qpoints among the processors
!
!=============================================================================

#include "f_defs.h"

subroutine distrib(peps)
  use global_m
  use realspace_common_m
  use scalapack_m
  implicit none
  type (realspace_t), intent(in) :: peps
  integer :: ifq, ipe, modnfqnpes, nfq_loc_count, ifq_loc, inode, ipes, ik_loc
  integer, allocatable :: startindex(:)
  PUSH_SUB(distrib)
  
  !> Blocksize
  peinf%nkpe = iceil(peps%nfq, peinf%npes)
  !> (global) peinf%ik(iproc, ifq_loc) = ifq_global
  SAFE_ALLOCATE(peinf%ik, (peinf%npes, peinf%nkpe))
  !> (global) peinf%ikt(iproc) = # of ifq_global dealt by each proc
  SAFE_ALLOCATE(peinf%ikt, (peinf%npes))
  peinf%ik  = 0
  peinf%ikt = 0
  
  do ipes = 1, peinf%npes
     peinf%ikt(ipes) = NUMROC(peps%nfq, peinf%nkpe, ipes-1, 0, peinf%npes)
     do ik_loc = 1, peinf%ikt(ipes)
        peinf%ik(ipes, ik_loc) = INDXL2G(ik_loc, peinf%nkpe, ipes-1, 0, peinf%npes)
     enddo
  enddo

  if (peinf%inode .eq. 0) then
     write(6,'(1X,A)') "Distribute fq points:"
     write(6,'(1X,A,I0)') "- Maximal number of full BZ qpoint on a proc: ", peinf%nkpe
     write(6,'(1X,A,I0)') "- Number of idle procs : ", COUNT(peinf%ikt(:) == 0)
     write(6,'(A)')     
  endif
  
  ! if ( peinf%inode .eq. 0 ) then
  !    do inode = 0, peinf%npes - 1
  !       write(6,'(A,I5,A,I5,A)') "Process #", inode, " nfq_loc = ", peinf%ikt(inode+1), " ifq_global : "
  !       ! if (peinf%ikt(inode+1) .gt. 0 ) then
  !       !    write(6,'(10I5)') peinf%ik(inode+1,1:peinf%ikt(inode+1))
  !       ! endif
  !    enddo
  ! endif

  POP_SUB(distrib)
  return
end subroutine distrib

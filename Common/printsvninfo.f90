!================================================================================
!
! Program:
!
! printsvninfo      Originally By DAS
!
! Returns information on svn repository name, version, and revision number.
! For use by scripts to add version info to output.
!
!================================================================================

#include "f_defs.h"

program printsvninfo

  use global_m
  use svninfo_m
  implicit none

  character*256 :: string

  call getsvninfo(string)

  write(6,'(a)') trim(string)

end program printsvninfo

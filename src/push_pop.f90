!==========================================================================
!
! Module push_pop_m     Originally By DAS
!
!   Create a stack trace of routines entered and exited, for debugging.
!   This only takes effect if the code is compiled with -DDEBUG.
!   Enable by setting 'debug_level' below, and recompile:
!     0: no debugging trace
!     1: only node 0 writes trace
!     2: all nodes write trace. Very slow.
!   Inspired by Octopus messages.F90 (originally revision 6920)
!
!==========================================================================

#include "f_defs.h"

module push_pop_m

  use message_m
  use nrtype_m
  use peinfo_m

  implicit none

#ifdef DEBUG

  private

  public ::              &
    open_debug_trace,    &
    push_sub,            &
    pop_sub,             &
    set_debug_level,     &
    operator(+)
  ! be sure to make visible '+' from message_m as used in PUSH/POP statements
  
  integer, parameter, private :: offset = 1000
  integer, parameter, private :: namelength = 100
  integer, private            :: debug_level = 0
  !< 0: no debugging trace
  !! 1: only node 0 writes trace
  !! 2: all nodes write trace. Very slow.
  logical, private            :: first_write = .true.

 !> The stack.
  character(len=namelength), public :: sub_stack(50)
  real(DP), public                  :: time_stack(50)
  integer, public                   :: stacksize = 0

contains

  ! ---------------------------------------------------------
  subroutine set_debug_level(level)
    integer, intent(in) :: level

    if(level < 0 .or. level > 2) call die("Invalid value for debug_level")
    debug_level = level
    
    return
  end subroutine set_debug_level

  ! ---------------------------------------------------------

  subroutine open_debug_trace(iunit)
    integer, intent(out) :: iunit

    character(len=6) :: filenum, pos
    character(len=7) :: status

    iunit = peinf%inode + offset
    write(filenum, '(i6.6)') iunit - offset

    if(first_write) then
      pos = 'rewind'
      first_write = .false.
      status = 'replace'
    else
      pos = 'append'
      status = 'old'
    endif
    call open_file(iunit, 'debug_trace.node.'//filenum, position=pos, status=status)

  end subroutine open_debug_trace

  ! ---------------------------------------------------------

  subroutine push_sub(sub_name)
    character(len=*), intent(in) :: sub_name

    integer :: iunit
    character(len=namelength) :: sub_name2

    if(debug_level == 0) return

    stacksize = stacksize + 1
    if(stacksize > 49) then
      sub_stack(50) = 'push_sub'
      call die('Too many recursion levels (max=50)')
    end if

    call trim_sub_name(sub_name, sub_name2)
    sub_stack(stacksize) = trim(sub_name2)

    if(debug_level > 1) then
      call open_debug_trace(iunit)
      call push_sub_write(iunit)
      ! close file to ensure flushing
      call close_file(iunit)
    else if(debug_level == 1 .and. peinf%inode == 0) then
      ! write to stderr if we are node 0
      call push_sub_write(0)
    end if

  contains

    subroutine push_sub_write(iunit_out)
      integer,  intent(in) :: iunit_out

      integer :: ii

      write(iunit_out,'(a)', advance='no') "* I | "
      do ii = stacksize - 1, 1, -1
        write(iunit_out,'(a)', advance='no') "..|"
      end do
      write(iunit_out,'(a)') trim(sub_name2)

    end subroutine push_sub_write

  end subroutine push_sub

  ! ---------------------------------------------------------

  subroutine pop_sub(sub_name)
    character(len=*), intent(in) :: sub_name

    integer :: iunit
    character(len=namelength) :: sub_name2
    character(len=namelength*3) :: tmp_str

    if(debug_level == 0) return

    call trim_sub_name(sub_name, sub_name2)

    if(stacksize <= 0) then
      stacksize = 1
      sub_stack(1) = 'pop_sub'
      write(tmp_str,'(a)') 'pop ' // trim(sub_name2) // ': Too few recursion levels.'
      call die(tmp_str)
    end if

    if(sub_name2 .ne. sub_stack(stacksize)) then
      write (tmp_str,'(a,3x,a,3x,a)') 'Wrong sub name on pop_sub ', &
        trim(sub_name2), trim(sub_stack(stacksize))
      call die(tmp_str)
    end if

    if(debug_level > 1) then
      call open_debug_trace(iunit)
      call pop_sub_write(iunit)
      ! close file to ensure flushing
      call close_file(iunit)
    else if (debug_level == 1 .and. peinf%inode == 0) then
      ! write to stderr if we are node 0
      call pop_sub_write(0)
    end if
    
    stacksize = stacksize - 1

  contains

    subroutine pop_sub_write(iunit_out)
      integer, intent(in) :: iunit_out

      integer :: ii

      write(iunit_out,'(a)', advance='no') "* O | "
      do ii = stacksize - 1, 1, -1
        write(iunit_out,'(a)', advance='no') "..|"
      end do
      write(iunit_out,'(a)') trim(sub_stack(stacksize))

    end subroutine pop_sub_write

  end subroutine pop_sub

  ! ---------------------------------------------------------

  ! remove leading directory names and slashes from routine names
  subroutine trim_sub_name(sub_name, output)
    character(len=*), intent(in)  :: sub_name
    character(len=*), intent(out) :: output

    integer :: ic, slash

    slash = 0 ! if no slash, take whole name
    do ic = 1, len(sub_name)
      if(sub_name(ic:ic) == '/') slash = ic
    enddo

    if(slash >= len(sub_name) - 1) then
      output = ""
    else
      output = sub_name(slash + 1:len(sub_name))
    endif

  end subroutine trim_sub_name

#endif

end module push_pop_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

!===================================================================
!
! Modules:
!
! 1. io_utils_m      Originally By FHJ
!
!    Provides set of routines to standardize stdout operations.
!
!    progress_info: tracks the progress of a task
!
!===================================================================

#include "f_defs.h"

module io_utils_m

  use global_m

  implicit none

  private

  !> Progress information object. Use this object whenever you want to have
  !! a nice standardized progress report on a given task, including an
  !! estimate of the remaining time to complete the task.
  !!
  !! \sa progress_init(), progress_step() and progress_free routines().
  type progress_info
     real(DP) :: t_pred  !< A prediction of when task will be done.
     !> Time when progress_init was called (before 1st step got executed).
     real(DP) :: t_wall0
     !> Time elapsed, wrt t_wall0, when 1st step finished (before 2nd step got executed).
     real(DP) :: dt1
     real(DP) :: report_every !< How often to report progress to stdout?
     integer :: num_steps !< Max. number of steps (not number of reports!)
     integer :: cur_step  !< Keep track of the current step
     !> What is the task we are tracking? This should be a noun(-phrase). Eg:
     !! "reading wavefunctions", "calculation of matrix elements", etc.
     character(len=256) :: str_task
     !> What defines each step/iteration? This should be a noun(-phrase). Eg:
     !! "k-point", "transition", etc.
     character(len=64) :: str_step
     character(len=11) :: str_tot !< this is the string " / $num_steps"
     logical :: should_print !< Should I write out my progress? Default is .true. for inode==0
     integer :: iunit=6 !< Defaults to 6 (stdout)
  end type progress_info

  public ::              &
       centralize,          &
       progress_info,       &
       progress_init,       &
       progress_step,       &
       progress_free,       &
       print_dealing_with

contains

  !> Trim and centralize a string. Max width hard-coded to 256.
  subroutine centralize(text_in, text_out, width)
    character(len=*), intent(in) :: text_in
    character(len=*), intent(out) :: text_out
    integer, intent(in) :: width

    character(len=256) :: text
    integer :: text_len, margin, wid

    PUSH_SUB(centralize)

    wid = width
    if (wid>256) wid=256
    if (len(TRUNC(text_in))>wid) then
       ! FHJ: avoid buffer overflow.
       text_out = text_in(1:wid)
       POP_SUB(centralize)
       return
    endif

    text = TRUNC(text_in)
    text_len = len(TRUNC(text))
    margin = (wid-text_len)/2
    text_out = repeat(' ', len(text_out))
    text_out(margin+1:margin+text_len) = text(1:text_len)

    POP_SUB(centralize)
    return

  end subroutine centralize

  !----------------------------------------------------------------------------
  ! FHJ: Progress report stuff
  !----------------------------------------------------------------------------

  !> Initialize a ::progress_info structure.
  !!
  !! This subroutine will also print the following sentence:
  !! "Started ${str_task} with ${num_steps} ${str_step}s at ${TIME}."
  subroutine progress_init(prog_info, str_task, str_step, num_steps, num_reports, should_print, iunit)
    type(progress_info), intent(out) :: prog_info
    character(len=*), intent(in) :: str_task !< See: progress_info%str_task
    character(len=*), intent(in) :: str_step !< See: progress_info%str_step
    integer, intent(in) :: num_steps !< How many times will you call progress_step?
    !> What is the total number of times we want to write the progress to
    !! stdout? Defaults to 10.
    integer, intent(in), optional :: num_reports
    logical, intent(in), optional :: should_print
    integer, intent(in), optional :: iunit

    character(len=15) :: mydate, mytime
    character(len=256) :: str_tmp
    real(DP) :: t_cpu
    integer :: num_rep

    PUSH_SUB(progress_init)

    prog_info%dt1 = 0d0
    prog_info%t_pred = 0d0
    prog_info%str_task = str_task
    prog_info%str_step = str_step
    prog_info%num_steps = num_steps
    prog_info%should_print = peinf%inode==0
    prog_info%iunit = 6
    num_rep = 10
    if (present(num_reports)) num_rep = num_reports
    if (present(should_print)) prog_info%should_print = should_print
    if (present(iunit)) prog_info%iunit = iunit

    if (num_steps<num_rep) then
       prog_info%report_every = 1d0
    else
       prog_info%report_every = dble(num_steps)/num_rep
    endif
    prog_info%cur_step = 0
    call timget(t_cpu, prog_info%t_wall0)

    if (prog_info%should_print) then
       call date_and_time(mydate,mytime)
       write(prog_info%iunit,*)
       write(str_tmp,'(i8)') num_steps
       ! FHJ: eg: "Started reading WFNs with 65536 k-point(s) at 12:00:00."
       write(prog_info%iunit,'(1x,4a,1x,3a)') &
            'Started ',TRUNC(str_task), ' with ', TRUNC(str_tmp), TRUNC(str_step), &
            '(s) at ', mytime(1:2)//':'//mytime(3:4)//':'//mytime(5:6)//'.'
       call centralize(str_step, str_tmp, 19)
       write(str_tmp,'(i8)') num_steps
       write(prog_info%str_tot,'(a,a)') ' / ', TRUNC(str_tmp)
    endif

    POP_SUB(progress_init)

  end subroutine progress_init

  !> Update prediction info for current task.
  !!
  !! Call this function at the beginning of your loop. The variable
  !! prog_info%t_pred will be updated with the prediction of when the task will
  !! be finished. Depending on the current step, some info about the progress
  !! of the current task will be printed. Although this function is pretty
  !! lightweight, you might consider not placing the call in the innermost loop
  !! if only a very simple operation is performed each step.
  subroutine progress_step(prog_info, step_)
    type(progress_info), intent(inout) :: prog_info
    integer, intent(in), optional :: step_ !< Current step/iteration

    !> Contains an estimate for the time remaining to finish the task
    character(len=11) :: str_remaining_tmp
    character(len=64) :: str_remaining
    real(DP) :: dt, t_cpu, t_remain
    character(len=15) :: mydate, mytime
    integer :: step, idx1, idx2
    character(len=15) :: step_str

    PUSH_SUB(progress_step)

    if (present(step_)) then
       prog_info%cur_step = step_
    else
       prog_info%cur_step = prog_info%cur_step + 1
    endif
    step = prog_info%cur_step

    if (step<=3 .or. mod(dble(step), prog_info%report_every)<1d0) then
       call timget(t_cpu, dt)
       str_remaining = '.'
       dt = dt - prog_info%t_wall0
       if (step==2) then
          prog_info%dt1 = dt
       elseif (step>2) then
          ! FHJ: to get a better estimate, we compare the time wrt step #2, since
          ! step #1 is usually atypical (setup_FFTs, allocate buffers, etc.)
          t_remain = (dt - prog_info%dt1) / (step-2d0) * (prog_info%num_steps-step+1d0)
          prog_info%t_pred = prog_info%t_wall0 + t_remain
          write(str_remaining_tmp, '(i6)') nint(t_remain)
          write(str_remaining, '(3a)') ', remaining: ', TRUNC(str_remaining_tmp),' s.'
       endif

       if (prog_info%should_print .and. ((step==1.or.step==3).or.&
            mod(dble(step), prog_info%report_every)<1d0)) then
          call date_and_time(mydate,mytime)
          ! FHJ: This ensures that the two numbers in (xx/yy) have the same width
          write(step_str, '(i15)') step
          idx2 = LEN(TRUNC(prog_info%str_tot)) - 2
          idx1 = 16 - idx2
          idx2 = idx1 + idx2 - 1
          write(prog_info%iunit,'(1x,3a,i3,a)', advance='no') & ! [ 12:00:00 |  20% ]
               '[ ', mytime(1:2)//':'//mytime(3:4)//':'//mytime(5:6), ' | ', &
               nint((step - 1d0)/dble(prog_info%num_steps)*1d2), '% ]'
          write(prog_info%iunit,'(1x,a,1x,a,1x,a)', advance='no') & ! transition  16 / 128
               TRUNC(prog_info%str_step), step_str(idx1:idx2), TRUNC(prog_info%str_tot)
          write(prog_info%iunit,'(a)') TRUNC(str_remaining)    !, remaining:    10 s.
       endif
    endif

    POP_SUB(progress_step)

  end subroutine progress_step

  !> Finalize a progress information task.
  !!
  !! This will essentially print the following sentence:
  !! "Finished ${str_task} at ${TIME}."
  subroutine progress_free(prog_info)
    type(progress_info), intent(in) :: prog_info

    real(DP) :: dt, t_cpu
    character(len=15) :: mydate, mytime

    PUSH_SUB(progress_free)

    if (prog_info%should_print) then
       call date_and_time(mydate,mytime)
       write(prog_info%iunit,'(1x,4a)') 'Finished ', TRUNC(prog_info%str_task), &
            ' at ', mytime(1:2)//':'//mytime(3:4)//':'//mytime(5:6)//'.'
       call timget(t_cpu, dt)
       dt = dt - prog_info%t_wall0
       write(mytime,'(i8)') nint(dt)
       write(prog_info%iunit,'(1x,a,a,a)') 'Elapsed time: ', TRUNC(mytime), ' s.'
       write(prog_info%iunit,*)
    endif

    POP_SUB(progress_free)

  end subroutine progress_free

  !> Print a banner with the new q/k point we are calculating
  subroutine print_dealing_with(ik, ik_max, kk, label, iunit)
    integer, intent(in) :: ik !< Index of the current k-point
    integer, intent(in) :: ik_max !< Total numeber of k-points
    real(DP), intent(in) :: kk(:) !< (3) K-point as a vector
    character(len=*), intent(in) :: label !< Either "k" or "q"
    integer, optional, intent(in) :: iunit

    !> Contains an estimate for the time remaining to finish the task
    real(DP) :: t_cpu, t_wall
    character(len=15) :: mydate, mytime
    character(len=15) :: str_cur, str_max
    character(len=30) :: str_sulfix
    integer :: ii, idx1, idx2, iunit_

    PUSH_SUB(print_dealing_with)

    call timget(t_cpu, t_wall)
    call date_and_time(mydate,mytime)
    iunit_ = 6
    if (present(iunit)) iunit_ = iunit
    write(iunit_,'(A)') repeat('=', 80)
    write(iunit_,'(1X,A)', advance='no') mytime(1:2)//':'//mytime(3:4)//':'//mytime(5:6)
    write(iunit_,'(3X,A,3F20.10)', advance='no') 'Dealing with '//label(1:1)//' =', (kk(ii),ii=1,3)
    ! wrote 58 characters
    write(str_max,'(I15)') ik_max
    write(str_cur,'(I15)') ik
    idx2 = LEN(TRUNC(str_max))
    idx1 = 16 - idx2
    idx2 = idx1 + idx2 - 1
    write(str_sulfix,'(A)') str_cur(idx1:idx2)//' / '//TRUNC(str_max)
    idx1 = 80 - 58 - LEN(TRUNC(str_sulfix))
    write(iunit_,'(A)') repeat(' ', idx1)//TRUNC(str_sulfix)
    write(iunit_,'(A,/)') repeat('=', 80)

    POP_SUB(print_dealing_with)

  end subroutine print_dealing_with

end module io_utils_m

#include "f_defs.h"

!========================================================================
!
! Routines:
!
! (1) date_time()       Originally by ?         Last Modified: 5/12/2008 (JRD)
!
!     Gets current date and time.
!
! (2) timget()          Originally by gsm       Last Modified: 4/29/2010 (gsm)
!
!     Gets current cpu and wall time.
!
! (3) timacc(n,option,tsec,nslices)   Originally by ?
!                                               Last Modified: 6/17/2009 (PWD)
!
!     Timing subroutine.  Calls machine-dependent subroutine timget
!     which returns elapsed cpu and wall clock times in seconds
!     Also return the number of times the counter has been called
!
!     Depending on value of "option" routine will:
!       (0) zero all accumulators
!       (1) start with new incremental time slice for accumulator n
!           also increase by one the counter for this accumulator
!       (2) stop time slice; add time to accumlator n
!       (3) report accumulated time for accumulator n
!           and number of time that the routine has been called
!       (4) report time slice for accumulator n (not full time accumulated)
!
!       If, on first entry, subroutine is not being initialized, it
!       will automatically initialize as well as rezero accumulator n.
!       However, initialization SHOULD be done explicitly by the user
!       so that it can be done near the top of his/her main routine.
!
!       Input:
!         n=index of accumulator (distinguish what is being timed); not used if option=0
!         option=see comment above
!       Output:
!         on option=3:
!         tottim(2,n)=accumulated time for accumulator n; otherwise
!         tottim is a dummy variable.
!         nslices is optional variable that give number of slices collected
!
! (4) logit()    Originally By (SIB)     Last Modified 6/12/2008 (JRD)
!
!    Write out a debugging message with an inputed string and write time.
!
! (5) logitint() Originally By (SIB)     Last Modified 6/12/2008 (JRD)
!
!    Same as logit but with an integer constant.
!
!========================================================================

module timing_m

  use intrinsics_m
  use message_m
  use nrtype_m
  use peinfo_m
  use push_pop_m

  implicit none

  private

  public ::    &
       date_time, &
       timget,    &
       timacc,    &
       logit,     &
       logitint

  !> MTIM determines the maximum number of "timing slots" available
  integer, parameter, private :: MTIM=100
  real(DP), private, save :: acctim(2,MTIM),tzero(2,MTIM)
  integer, private, save :: ncount(MTIM)

contains

  subroutine date_time(bdate,btime)
    character, intent(out) :: bdate*11,btime*14

    integer :: lmonth
    integer :: idate (8)
    character :: day*2,year*4
    character :: adate*8,atime*10,azone*5
    character :: hour*2,min*2,sec*2
    character*3 :: month(12)

    DATA month/'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep', &
         'Oct','Nov','Dec'/

    PUSH_SUB(date_time)

    call date_and_time(adate,atime,azone,idate)

    read(adate,101) year,lmonth,day
101 format(a4,i2,a2)
    write(bdate,102) day,month(lmonth),year
102 format(a2,'-',a3,'-',a4)
    read(atime,201) hour,min,sec
201 format(a2,a2,a2,a4)
    write(btime,202) hour,min,sec,azone
202 format(a2,':',a2,':',a2,1x,a5)

    POP_SUB(date_time)

    return
  end subroutine date_time

  !================================================================================

  subroutine timget(cpu,wall)
    real(DP), intent(out) :: cpu,wall

    integer :: values(8)

    ! no push_sub, called too frequently

    TIMGET(cpu)

    call date_and_time(VALUES=values)
    wall=((values(3)*24.0d0+values(5))*60.0d0 &
         +values(6))*60.0d0+values(7)+values(8)*1.0d-3

    return
  end subroutine timget

  !================================================================================

  subroutine timacc(n,option,tottim,nslices)
    integer, intent(in) :: n !< not used for option = 0
    integer, intent(in) :: option !< 0, 1, 2, 3, 4
    real(DP), intent(out), optional :: tottim(2) !< should be present if option=3 or 4
    integer, intent(out), optional :: nslices !< only used if option=3, still optional in that case

    real(DP) :: cpu,wall
    character*100 :: tmpstr

    ! no push_sub, called too frequently

    ! Check that n lies in sensible bounds

    if (n .lt. 0 .or. n .gt. MTIM) then
       write(tmpstr,'(a,i6,a,i8)') 'timacc: dim MTIM = ', MTIM,' but input n =', n
       call die(tmpstr)
    end if

    if (option==0) then

       ! Zero out all accumulators of time and init timers

       acctim(:,:)=0.0d0
       tzero(:,:)=0.0d0
       ncount(:)=0

    else if (option==1) then

       ! Initialize timepw for n

       call timget(cpu,wall)
       tzero(1,n)=cpu
       tzero(2,n)=wall

    else if (option==2) then

       ! Accumulate time for n

       call timget(cpu,wall)
       acctim(1,n)=acctim(1,n)+cpu -tzero(1,n)
       acctim(2,n)=acctim(2,n)+wall-tzero(2,n)
       ncount(n)=ncount(n)+1

    else if (option==3) then

       ! Return accumulated time for n

       if(.not. present(tottim)) call die("timacc requires tottim for option 3.")

       tottim(1)=acctim(1,n)
       tottim(2)=acctim(2,n)
       if(present(nslices)) then
          nslices=ncount(n)
       end if

    else if (option==4) then

       ! Return elapsed time for n (do not accumulate)

       if(.not. present(tottim)) call die("timacc requires tottim for option 4.")

       call timget(cpu,wall)
       tottim(1)=cpu-tzero(1,n)
       tottim(2)=wall-tzero(2,n)

    else

       write(tmpstr,'(a,i10,a)') 'timacc: input option = ', option, 'not valid.'
       call die(tmpstr)

    end if

    return
  end subroutine timacc

  !=====================================================================

  subroutine logit(str, should_print, iunit)
    character (len=*), intent(in) :: str
    logical, intent(in), optional :: should_print
    integer, intent(in), optional :: iunit

    character*15 :: mydate,mytime,tmpstr
    logical :: should_print_
    integer :: iunit_

    if (.not.peinf%verb_log) return

    iunit_ = 6
    if (present(iunit)) iunit_ = iunit
    should_print_ = peinf%inode==0
    if (present(should_print)) should_print_ = should_print

    if (should_print_) then
       call date_and_time(mydate,mytime)
       tmpstr = mytime(1:2)//':'//mytime(3:4)//':'//mytime(5:6)//'.'//mytime(8:10)
       mytime = tmpstr
       write(iunit_,*) '*** LOG: ', TRUNC(str),'  time = ', TRUNC(mytime)
    endif

  end subroutine logit

  !=====================================================================

  subroutine logitint(str,i)
    character (len=*), intent(in) :: str
    integer, intent(in) :: i

    character*100 :: tmpstr
    
    if (.not.peinf%verb_log) return

    write(tmpstr,'(a,i5)') str(1:len_trim(str)),i
    call logit(tmpstr)

  end subroutine logitint

  !=====================================================================

end module timing_m

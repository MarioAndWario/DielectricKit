!>===================================================================
!!
!! Module misc_m
!!
!! Routines:
!!
!! 3. get_volume()       Originally By (SIB)     Last Modified 6/12/2008 (JRD)
!!
!!    This assumes that b is a symmetric matrix.  It sets
!!    vol = (2pi)^3 / square_root(|det(b)|)
!!    This makes sense if b is the matrix of dot products of the recip
!!    lattice vectors, so vol is the real space volume.
!!
!! 4. findvector()       Originally By (SIB)     Last Modified 6/12/2008 (JRD)
!!
!!    Looks for the vector in the list of vectors
!!    gvec%components(1:3,1:gvec%ng).  If found, iout is its index.  Otherwise
!!    iout is zero.
!!
!! 11. procmem()          Originally By (gsm)     Last Modified 4/14/2009 (gsm)
!!
!!     Determines the amount of free memory per processor
!!     from the proc file system
!!
!! 12. sizeof_scalar()    Originally By (DAS)     Last Modified 1/25/2011 (DAS)
!!
!!     Return the size of the SCALAR type, for memory estimation.
!!
!!
!!===================================================================

#include "f_defs.h"

module misc_m
  use global_m
  use blas_m
  use scalapack_m
  implicit none
  private
  public ::   checknorm, get_volume, findvector, procmem, sizeof_scalar, &
       get_gumk3, get_gumk4, M33INV, output_header, generate_full_WS_BZ
contains

  subroutine checknorm(filename, iband, ik, nspin, wfn, ispin, tol)
    character (len=*), intent(in) :: filename
    integer, intent(in) :: iband,ik, nspin
    SCALAR, intent(in) :: wfn(:,:) !< (ng,nspin*nspinor)
    integer, optional, intent(in) :: ispin
    real(DP), optional, intent(in) :: tol
    real(DP) :: tol_
    SCALAR :: norm
    integer :: wfnsize
    SCALAR, allocatable :: wfn_1d(:)

    PUSH_SUB(checknorm)

    if (present(tol)) then
       tol_ = tol
    else
       tol_ = TOL_SMALL
    endif
    !> wfnsize = ng * nspinor
    wfnsize = SIZE(wfn)/nspin
    SAFE_ALLOCATE(wfn_1d, (wfnsize))

    if (nspin .eq. 1) then
       !> Do a dot product, make sure that the result has 1.0D0 as real part and 0.0D0 as imaginary part
       !> DOT(CONJG(wfn),wfn)

       !> Flatten the 2D wfn into 1D wfn_1d
       wfn_1d = PACK(wfn, .true.)
       norm = blas_dot(wfnsize,wfn_1d,1,wfn_1d,1)

       if (ABS(norm - 1.0D0) > tol_) then
          write(0,555) TRUNC(filename),ABS(norm-1.0D0),iband,1,ik
          call die("Incorrect normalization.")
       endif
    elseif (nspin .eq. 2) then
       if (present(ispin)) then
          if (ispin .eq. 1) then
             wfn_1d = PACK(wfn(:,1), .true.)
             norm = blas_dot(wfnsize,wfn_1d,1,wfn_1d,1)

             if (ABS(norm - 1.0D0) > tol_) then
                write(0,555) TRUNC(filename),ABS(norm-1.0D0),iband,1,ik
                call die("Incorrect normalization.")
             endif
          elseif (ispin .eq. 2) then
             wfn_1d = PACK(wfn(:,2), .true.)
             norm = blas_dot(wfnsize,wfn_1d,1,wfn_1d,1)

             if (ABS(norm - 1.0D0) > tol_) then
                write(0,555) TRUNC(filename),ABS(norm-1.0D0),iband,2,ik
                call die("Incorrect normalization.")
             endif
          else
             call die("wrong input of ispin.")
          endif
       else
          wfn_1d = PACK(wfn(:,1), .true.)
          norm = blas_dot(wfnsize,wfn_1d,1,wfn_1d,1)

          if (ABS(norm - 1.0D0) > tol_) then
             write(0,555) TRUNC(filename),ABS(norm-1.0D0),iband,1,ik
             call die("Incorrect normalization.")
          endif

          wfn_1d = PACK(wfn(:,2), .true.)
          norm = blas_dot(wfnsize,wfn_1d,1,wfn_1d,1)

          if (ABS(norm - 1.0D0) > tol_) then
             write(0,555) TRUNC(filename),ABS(norm-1.0D0),iband,2,ik
             call die("Incorrect normalization.")
          endif
       endif
    else
       call die("checknorm_2: Check nspin",only_root_writes=.true.)
    endif

    SAFE_DEALLOCATE(wfn_1d)

555 format(1x,'Wavefunction is not normalized in file',1x,a,/,&
         3x,'abs(norm - 1) =',f10.7,/,&
         3x,'iband =',i6,1x,'ispin =',i2,1x,'ik =',i6,/)

    POP_SUB(checknorm)

    return
  end subroutine checknorm

  subroutine get_volume(vol,b)
    real(DP), intent(out) :: vol
    real(DP), intent(in)  :: b(3,3)

    PUSH_SUB(get_volume)

    vol = b(1,1)*(b(2,2)*b(3,3) - b(2,3)**2) &
         + 2*b(1,2)*b(2,3)*b(3,1) &
         - b(2,2)*b(1,3)**2 - b(3,3)*b(1,2)**2
    vol = sqrt(abs(vol))
    vol = ((2.0d0*PI_D)**3)/vol

    POP_SUB(get_volume)

    return
  end subroutine get_volume

  !=====================================================================

  subroutine findvector(iout,kk,gvec)
    integer, intent(out) :: iout
    integer, intent(in)  :: kk(3)
    type (gspace), intent(in) :: gvec

    ! no push/pop since called too frequently

    iout=((kk(1)+gvec%FFTgrid(1)/2)*gvec%FFTgrid(2)+kk(2)+gvec%FFTgrid(2)/2)* &
         gvec%FFTgrid(3)+kk(3)+gvec%FFTgrid(3)/2+1
    if (iout .ge. 1 .and. iout .le. gvec%nFFTgridpts) then
       iout=gvec%index_vec(iout)
       if (iout .ge. 1 .and. iout .le. gvec%ng) then
          if (any(kk(1:3) /= gvec%components(1:3, iout))) iout = 0
       else
          iout = 0
       endif
    else
       iout = 0
    endif

    return
  end subroutine findvector

  !***********************************************************************************************************************************
  !  M33INV  -  Compute the inverse of a 3x3 matrix.
  !
  !  A       = input 3x3 matrix to be inverted
  !  AINV    = output 3x3 inverse of matrix A
  !  OK_FLAG = (output) .TRUE. if the input matrix could be inverted, and .FALSE. if the input matrix is singular.
  !***********************************************************************************************************************************

  SUBROUTINE M33INV (A, AINV, OK_FLAG)

    IMPLICIT NONE

    DOUBLE PRECISION, DIMENSION(3,3), INTENT(IN)  :: A
    DOUBLE PRECISION, DIMENSION(3,3), INTENT(OUT) :: AINV
    LOGICAL, INTENT(OUT) :: OK_FLAG

    DOUBLE PRECISION, PARAMETER :: EPS = 1.0D-10
    DOUBLE PRECISION :: DET
    DOUBLE PRECISION, DIMENSION(3,3) :: COFACTOR


    DET =   A(1,1)*A(2,2)*A(3,3)  &
         - A(1,1)*A(2,3)*A(3,2)  &
         - A(1,2)*A(2,1)*A(3,3)  &
         + A(1,2)*A(2,3)*A(3,1)  &
         + A(1,3)*A(2,1)*A(3,2)  &
         - A(1,3)*A(2,2)*A(3,1)

    IF (ABS(DET) .LE. EPS) THEN
       AINV = 0.0D0
       OK_FLAG = .FALSE.
       RETURN
    END IF

    COFACTOR(1,1) = +(A(2,2)*A(3,3)-A(2,3)*A(3,2))
    COFACTOR(1,2) = -(A(2,1)*A(3,3)-A(2,3)*A(3,1))
    COFACTOR(1,3) = +(A(2,1)*A(3,2)-A(2,2)*A(3,1))
    COFACTOR(2,1) = -(A(1,2)*A(3,3)-A(1,3)*A(3,2))
    COFACTOR(2,2) = +(A(1,1)*A(3,3)-A(1,3)*A(3,1))
    COFACTOR(2,3) = -(A(1,1)*A(3,2)-A(1,2)*A(3,1))
    COFACTOR(3,1) = +(A(1,2)*A(2,3)-A(1,3)*A(2,2))
    COFACTOR(3,2) = -(A(1,1)*A(2,3)-A(1,3)*A(2,1))
    COFACTOR(3,3) = +(A(1,1)*A(2,2)-A(1,2)*A(2,1))

    AINV = TRANSPOSE(COFACTOR) / DET

    OK_FLAG = .TRUE.

    RETURN

  END SUBROUTINE M33INV

  ! subroutine procmem(mem, ranks_per_node, nfreq_group)
  subroutine procmem(mem, ranks_per_node)
    real(DP), intent(out) :: mem !< Memory per MPI rank (lower bound)
    integer, intent(out) :: ranks_per_node !< # of MPI ranks per node (upper bound)
    ! integer, intent(in), optional :: nfreq_group

    integer, parameter :: host_len=64
    integer :: ierr,eof,info,iunit,m,n,i,pagesize
    real(DP) :: x,y,mac_m,mac_n
    !    integer :: ntot
    !    real(DP) :: xtot,ytot ! we do not use the total memory actually
    character(len=256) :: s, filename, my_hostname
    character(len=host_len), allocatable :: hostnames(:)

    PUSH_SUB(procmem)

    !-----------------------------------------------------
    !> determine the amount of free memory per node in kB

    m=0
    iunit=14
    x=0

    call open_file(unit=iunit,file='/proc/meminfo',form='formatted',iostat=ierr,status='old')
    if (ierr.eq.0) then
       eof=0
       do while(eof.eq.0)
          read(iunit,'(a)',iostat=eof)s
          if (s(1:7).eq."MemFree") then
             read(s(9:),*)n
             m=m+n
          endif
          !        if (s(1:8).eq."MemTotal") then
          !          read(s(10:),*)ntot
          !        endif
          if (s(1:6).eq."Cached") then
             read(s(8:),*)n
             m=m+n
          endif
       enddo
       x=dble(m)/dble(peinf%npes)
       call close_file(iunit)
    endif

    if(m == 0) then
       !> this is for Mac OS
       !! total memory is accessible instead from sysctl -n hw.usermem
       write(filename,'(a,i9.9)') 'vm_stat_', peinf%inode
       SYSTEMCALL("vm_stat > " + TRUNC(filename) + " 2> /dev/null")
       !> Fortran 2008 would use execute_command_line instead
       !! even if the command failed, still open file in order to delete it
       call open_file(unit=iunit,file=TRUNC(filename),form='formatted',iostat=ierr,status='old')
       if (ierr.eq.0) then
          eof=0
          do while(eof.eq.0)
             read(iunit,'(a)',iostat=eof)s
             if (s(1:45).eq."Mach Virtual Memory Statistics: (page size of") then
                read(s(46:),*)pagesize ! in bytes
             endif
             if (s(1:11).eq."Pages free:") then
                read(s(12:),*) mac_n
                mac_m = mac_m + mac_n
             endif
             if (s(1:18).eq."Pages speculative:") then
                read(s(19:),*) mac_n
                mac_m = mac_m + mac_n
             endif
          enddo
          call close_file(iunit, delete = .true.)
          x = mac_m * dble(pagesize) / dble(peinf%npes * 1024) ! to kB
       endif
    endif

    !> === Example output from vm_stat ===
    !! Mach Virtual Memory Statistics: (page size of 4096 bytes)
    !! Pages free:                           2886.
    !! Pages active:                       139635.
    !! Pages inactive:                      66906.
    !! Pages speculative:                    2376.
    !! Pages wired down:                    50096.
    !! "Translation faults":            123564742.
    !! Pages copy-on-write:              10525831.
    !! Pages zero filled:                53274329.
    !! Pages reactivated:                  739514.
    !! Pageins:                           2282166.
    !! Pageouts:                           306727.
    !! Object cache: 25 hits of 522230 lookups (0% hit rate)

    if(m == 0 .and. mac_m == 0) then ! BSD
       !> http://mario79t.wordpress.com/2008/08/29/memory-usage-on-freebsd/
       !! -bash-2.05b$ sysctl vm.stats.vm.v_free_count
       !! vm.stats.vm.v_free_count: 29835
       !! -bash-2.05b$ sysctl vm.stats.vm.v_page_count
       !! vm.stats.vm.v_page_count: 124419
       !! -bash-2.05b$ sysctl hw.pagesize
       !! hw.pagesize: 4096
       write(filename,'(a,i9.9)') 'sysctl_', peinf%inode
       SYSTEMCALL("sysctl -a > " + TRUNC(filename) + " 2> /dev/null")
       !> Fortran 2008 would use execute_command_line instead
       !! even if the command failed, still open file in order to delete it
       call open_file(unit=iunit,file=TRUNC(filename),form='formatted',iostat=ierr,status='old')
       if (ierr.eq.0) then
          eof=0
          do while(eof.eq.0)
             read(iunit,'(a)',iostat=eof)s
             if (s(1:12).eq."hw.pagesize:") then
                read(s(13:),*)pagesize ! in bytes
             endif
             if (s(1:25).eq."vm.stats.vm.v_free_count:") then
                read(s(26:),*) mac_n
                mac_m = mac_m + mac_n
             endif
             if (s(1:26).eq."vm.stats.vm.v_cache_count:") then
                read(s(27:),*) mac_n
                mac_m = mac_m + mac_n
             endif
          enddo
          call close_file(iunit, delete = .true.)
          x = mac_m * dble(pagesize) / dble(peinf%npes * 1024) ! to kB
       endif
    endif

    !    xtot=dble(ntot)/dble(peinf%npes)

#ifdef MPI
    call MPI_Allreduce(x,y,1,MPI_REAL_DP,MPI_SUM,MPI_COMM_WORLD,mpierr)
    !    call MPI_Allreduce(xtot,ytot,1,MPI_REAL_DP,MPI_SUM,MPI_COMM_WORLD,mpierr)
#else
    y=x
    !    ytot=xtot
#endif

    !----------------------------------------------
    ! Determine the number of processors per node
    ! FHJ: we do this in parallel. Each MPI rank compares with all other hostnames,
    ! and figures out the number of ranks running on the same node. At the end, we
    ! pick the upper bound for ranks_per_node.

    ! if(present(nfreq_group)) then
    !    SAFE_ALLOCATE(hostnames, (peinf%npes_orig))
    ! else
    SAFE_ALLOCATE(hostnames, (peinf%npes))
    ! endif
    HOSTNAMECALL(my_hostname, info)
    hostnames(peinf%inode+1) = my_hostname(1:host_len)

#ifdef MPI
    call MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, hostnames(1), &
         host_len, MPI_BYTE, MPI_COMM_WORLD, mpierr)
    ranks_per_node = 0
    do i = 1, peinf%npes
       ! FHJ: no need to trim, all strings have the same length
       if (hostnames(peinf%inode+1)==hostnames(i)) ranks_per_node = ranks_per_node + 1
    enddo
    call MPI_Allreduce(MPI_IN_PLACE, ranks_per_node, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, mpierr)
#else
    ranks_per_node = 1
#endif

    !-----------------------------------
    !> report the available memory in B

    if (ranks_per_node>1) y = y / dble(ranks_per_node)
    mem=y*1024.0d0

    !-----------------------------------
    !> warn if zero memory

    if (mem .lt. TOL_Small .and. peinf%inode .eq. 0) then
       write(0,666)
666    format(1x,'WARNING: estimation of memory available is zero, probably failed.',/)
    endif

    POP_SUB(procmem)

    return
  end subroutine procmem

  !================================================================================

  !> for memory estimation, tell what size of SCALAR type is
  integer function sizeof_scalar()

#ifdef NOSIZEOF
    !> if no sizeof, take a guess

#ifdef CPLX
    sizeof_scalar = 16
#else
    sizeof_scalar = 8
#endif

#else

    SCALAR :: dummy
    sizeof_scalar = sizeof(dummy)

#endif

  end function sizeof_scalar

  !> (kvector(:) - dble(gumk_out(:))) falls into first WS BZ, having the shortest distance to Gamma
  subroutine get_gumk3(bdot, kvector, gumk_out)
    real(DP), intent(in) :: bdot(3,3)
    real(DP), intent(in) :: kvector(3)
    integer, intent(out) :: gumk_out(3)
    real(DP) :: len2_kvector_temp, len2_kvector_temp_min, kvector_temp(3), kvector_(3)
    integer :: b1, b2, b3, ib1, ib2, ib3, b_out(3), b1list(7), b2list(7), b3list(7), nb1, nb2, nb3, nb(3), id
    logical :: dim_flag(3)

    !> If kvector(id) = 0 and bdot(id,id) is block-diagonal, e.g. bdot(1,3)=DOT[b1,b3] = 0, bdot(2,3)=DOT[b2,b3]=0
    !> do not apply integer shift!
    if ( ABS(kvector(1)) < TOL_ZERO .and. ABS(bdot(1,2)) < TOL_ZERO .and. ABS(bdot(1,3)) < TOL_ZERO ) then
       nb1 = 1
       b1list = 0
    else
       nb1 = 7
       b1list(:) = (/0, -1, 1, -2, 2, -3, 3/)
    endif

    if ( ABS(kvector(2)) < TOL_ZERO .and. ABS(bdot(1,2)) < TOL_ZERO .and. ABS(bdot(2,3)) < TOL_ZERO ) then
       nb2 = 1
       b2list = 0
    else
       nb2 = 7
       b2list(:) = (/0, -1, 1, -2, 2, -3, 3/)
    endif

    if ( ABS(kvector(3)) < TOL_ZERO .and. ABS(bdot(1,3)) < TOL_ZERO .and. ABS(bdot(2,3)) < TOL_ZERO ) then
       nb3 = 1
       b3list = 0
    else
       nb3 = 7
       b3list(:) = (/0, -1, 1, -2, 2, -3, 3/)
    endif

    len2_kvector_temp_min = INF
    kvector_ = 0.0D0
    b_out = 0

    nb = (/ nb1, nb2, nb3 /)
    dim_flag = (nb > 1)

    do ib1 = 1, nb1
       b1 = b1list(ib1)
       do ib2 = 1, nb2
          b2 = b2list(ib2)
          do ib3 = 1, nb3
             b3 = b3list(ib3)

             kvector_temp(:) = kvector(:) - (/DBLE(b1), DBLE(b2), DBLE(b3)/)
             len2_kvector_temp = DOT_PRODUCT(kvector_temp, MATMUL(bdot, kvector_temp))

             !> If T, it means the kvector is at the edge of BZ
             if (ABS(len2_kvector_temp - len2_kvector_temp_min) < TOL_SMALL) then
                !> Determine whether we should replace b_out
                !> We try to make the first meaningful fractional coordinates >= 0 and as large as possible
                !> Since BZ has inversion symmetry, we can always find an BZ edge point with the first meaningful fractional coordinate >= 0
                loop_id: do id = 1, 3
                   if (.not. dim_flag(id)) cycle
                   !> kvector_temp(id)
                   !> Pick the larger fractional coordinate
                   if (kvector_temp(id) > kvector_(id)) then
                      len2_kvector_temp_min = len2_kvector_temp
                      kvector_(:) = kvector_temp(:)
                      b_out(:) = (/ b1, b2, b3 /)
                      !> If this fractional coordinate is the same, cycle to the next one
                   elseif ( ABS(kvector_temp(id) - kvector_(id)) < TOL_ZERO) then
                      cycle
                   endif
                   exit loop_id
                enddo loop_id
             elseif ( len2_kvector_temp < len2_kvector_temp_min ) then
                !> Determine whether this is at edge of BZ
                len2_kvector_temp_min = len2_kvector_temp
                kvector_(:) = kvector_temp(:)
                b_out(:) = (/ b1, b2, b3 /)
             endif
          enddo
       enddo
    enddo
    gumk_out(:) = b_out(:)

  end subroutine get_gumk3

  !> (kvector(:) - dble(gumk_out(:))) falls into first WS BZ, having the shortest distance to Gamma
  subroutine get_gumk4(bdot, kvector, gumk_out)
    real(DP), intent(in) :: bdot(3,3)
    real(DP), intent(in) :: kvector(3)
    integer, intent(out) :: gumk_out(3)
    real(DP) :: len2_kvector_temp, len2_kvector_temp_min, kvector_temp(3), kvector_(3)
    integer :: b1, b2, b3, ib1, ib2, ib3, b_out(3), b1list(9), b2list(9), b3list(9), nb1, nb2, nb3, nb(3), id
    logical :: dim_flag(3)

    !> If kvector(id) = 0 and bdot(id,id) is block-diagonal, e.g. bdot(1,3)=DOT[b1,b3] = 0, bdot(2,3)=DOT[b2,b3]=0
    !> do not apply integer shift!
    if ( ABS(kvector(1)) < TOL_ZERO .and. ABS(bdot(1,2)) < TOL_ZERO .and. ABS(bdot(1,3)) < TOL_ZERO ) then
       nb1 = 1
       b1list = 0
    else
       nb1 = 9
       b1list(:) = (/0, -1, 1, -2, 2, -3, 3, -4, 4/)
    endif

    if ( ABS(kvector(2)) < TOL_ZERO .and. ABS(bdot(1,2)) < TOL_ZERO .and. ABS(bdot(2,3)) < TOL_ZERO ) then
       nb2 = 1
       b2list = 0
    else
       nb2 = 9
       b2list(:) = (/0, -1, 1, -2, 2, -3, 3, -4, 4/)
    endif

    if ( ABS(kvector(3)) < TOL_ZERO .and. ABS(bdot(1,3)) < TOL_ZERO .and. ABS(bdot(2,3)) < TOL_ZERO ) then
       nb3 = 1
       b3list = 0
    else
       nb3 = 9
       b3list(:) = (/0, -1, 1, -2, 2, -3, 3, -4, 4/)
    endif

    len2_kvector_temp_min = INF
    kvector_ = 0.0D0
    b_out = 0

    nb = (/ nb1, nb2, nb3 /)
    dim_flag = (nb > 1)

    do ib1 = 1, nb1
       b1 = b1list(ib1)
       do ib2 = 1, nb2
          b2 = b2list(ib2)
          do ib3 = 1, nb3
             b3 = b3list(ib3)

             kvector_temp(:) = kvector(:) - (/DBLE(b1), DBLE(b2), DBLE(b3)/)
             len2_kvector_temp = DOT_PRODUCT(kvector_temp, MATMUL(bdot, kvector_temp))

             !> If T, it means the kvector is at the edge of BZ
             if (ABS(len2_kvector_temp - len2_kvector_temp_min) < TOL_SMALL) then
                !> Determine whether we should replace b_out
                !> We try to make the first meaningful fractional coordinates >= 0 and as large as possible
                !> Since BZ has inversion symmetry, we can always find an BZ edge point with the first meaningful fractional coordinate >= 0
                loop_id: do id = 1, 3
                   if (.not. dim_flag(id)) cycle
                   !> kvector_temp(id)
                   !> Pick the larger fractional coordinate
                   if (kvector_temp(id) > kvector_(id)) then
                      len2_kvector_temp_min = len2_kvector_temp
                      kvector_(:) = kvector_temp(:)
                      b_out(:) = (/ b1, b2, b3 /)
                      !> If this fractional coordinate is the same, cycle to the next one
                   elseif ( ABS(kvector_temp(id) - kvector_(id)) < TOL_ZERO) then
                      cycle
                   endif
                   exit loop_id
                enddo loop_id
             elseif ( len2_kvector_temp < len2_kvector_temp_min ) then
                !> Determine whether this is at edge of BZ
                len2_kvector_temp_min = len2_kvector_temp
                kvector_(:) = kvector_temp(:)
                b_out(:) = (/ b1, b2, b3 /)
             endif
          enddo
       enddo
    enddo
    ! write(*,'(A,3F15.8,A,3F15.8)') "k_in = ", kvector, "k_out = ", kvector_
    gumk_out(:) = b_out(:)

  end subroutine get_gumk4

  subroutine output_header(file_header, kp, gvec, syms, crys)
    character (len=*), intent(in) :: file_header
    type (kpoints), intent(in) :: kp
    type (crystal), intent(in) :: crys
    type (symmetry), intent(in) :: syms
    type (gspace), intent(in) :: gvec
    integer :: ascunit=2345, ib, ik, is, isym, ii, iatom
    integer :: temp_i3(3)
    real(DP) :: temp_r3(3)
    open( unit = ascunit, file = TRUNC(file_header), form = 'formatted', status = 'replace' )

    write(ascunit,'(A)') '====== kpoints:kp information ======'
    write(ascunit,'(A,I5)') 'kp%nspin = ', kp%nspin
    write(ascunit,'(A,I5)') 'kp%nspinor = ', kp%nspinor
    write(ascunit,'(A,I5)') 'kp%nrk = ', kp%nrk
    write(ascunit,'(A,I5)') 'kp%mnband = ', kp%mnband
    write(ascunit,'(A,I5)') 'kp%ngkmax = ', kp%ngkmax ! npwx_g : the maximum number of Gvector owned by one kpoints among all kpoints
    write(ascunit,'(A,G20.10)') 'kp%ecutwfc = ', kp%ecutwfc ! reals

    write(ascunit,'(A, 3I5)') 'kp%kgrid(1:3) = ', kp%kgrid(:) ! integers

    write(ascunit,'(A, 3G20.10)') 'kp%shift(1:3) = ', kp%shift(:) ! reals

    write(ascunit,'(A)') "------------------------------------"
    write(ascunit,'(A)') '   ik   kp%ngk(ik) :' ! integers
    do ik=1, kp%nrk
       write(ascunit,'(I5, I5)') ik, kp%ngk(ik)
    enddo

    write(ascunit,'(A)') "------------------------------------"
    write(ascunit,'(A)') '   ik   kp%w(ik) : ' ! real
    do ik=1, kp%nrk
       write(ascunit,'(I5, G20.10)') ik, kp%w(ik)
    enddo

    write(ascunit,'(A)') "------------------------------------"
    write(ascunit,'(A)') '   ik   kp%rk(1:3, ik) : ' ! real
    do ik=1, kp%nrk
       write(ascunit,'(I5, 3G20.10)') ik, kp%rk(:,ik)
    enddo

    write(ascunit,'(A)') "------------------------------------"
    write(ascunit,'(A)') '   ik   is   kp%ifmin(ik,is) kp%ifmax(ik,is):' ! integers
    do ik=1, kp%nrk
       do is=1, kp%nspin
          write(ascunit,'(I5, I5, 2I5)') ik, is, kp%ifmin(ik,is), kp%ifmax(ik,is)
       enddo
    enddo

    write(ascunit,'(A)') "------------------------------------"
    write(ascunit,'(A)') '   ib   ik   is   kp%el(ib,ik,is) : ' ! real
    do ib=1, kp%mnband
       do ik=1, kp%nrk
          do is=1, kp%nspin
             write(ascunit,'(I5, I5, I5, G20.10)') ib, ik, is, kp%el(ib,ik,is)
          enddo
       enddo
    enddo

    write(ascunit,'(A)') "------------------------------------"
    write(ascunit,'(A)') '   ib   ik   is   kp%occ(ib,ik,is) : ' ! real
    do ib=1, kp%mnband
       do ik=1, kp%nrk
          do is=1, kp%nspin
             write(ascunit,'(I5, I5, I5, G20.10)') ib, ik, is, kp%occ(ib,ik,is)
          enddo
       enddo
    enddo

    write(ascunit,'(A)') "====================================="

    write(ascunit,'(A)') "====== gspace:gvec information ======"

    write(ascunit,'(A, I5)') "gvec%ng = ", gvec%ng ! ngm_g : global number of Gvectors (summed over procs)
    write(ascunit,'(A, G20.10)') "gvec%ecutrho = ", gvec%ecutrho
    write(ascunit,'(A, 3I5)') "gvec%FFTgrid(1:3) = ", gvec%FFTgrid

    write(ascunit,'(A)') "====================================="


    write(ascunit,'(A)') "====== symmetry:syms information ======"
    write(ascunit,'(A, I5)') "syms%ntran = ", syms%ntran
    write(ascunit,'(A, I5)') "syms%cell_symmetry = ", syms%cell_symmetry

    write(ascunit,'(A)') "---------------------------------------"
    write(ascunit,'(A, I5)') "syms%mtrx(1:3,1:3,1:48) "
    do isym=1,syms%ntran
       write(ascunit,'("#", I5, ":")') isym
       do ii=1,3
          temp_i3(:) = syms%mtrx(ii, 1:3, isym)
          WRITE (ascunit,'(I8,1X,I8,1X,I8)') temp_i3(:)
       enddo
    enddo

    ! write(ascunit,'(A)') "---------------------------------------"
    ! write(ascunit,'(A, I5)') "syms%mtrx_reci(1:3,1:3,1:48) "
    ! do isym=1,syms%ntran
    !    write(ascunit,'("#", I5, ":")') isym
    !    do ii=1,3
    !       temp_i3(:) = syms%mtrx_reci(ii, 1:3, isym)
    !       WRITE (ascunit,'(I8,1X,I8,1X,I8)') temp_i3(:)
    !    enddo
    ! enddo

    ! write(ascunit,'(A)') "---------------------------------------"
    ! write(ascunit,'(A, I5)') "syms%mtrx_cart(1:3,1:3,1:48) "
    ! do isym=1,syms%ntran
    !    write(ascunit,'("#", I5, ":")') isym
    !    do ii=1,3
    !       temp_r3(:) = syms%mtrx_cart(ii, 1:3, isym)
    !       WRITE (ascunit,'(F12.5,1X,F12.5,1X,F12.5)') temp_r3(:)
    !    enddo
    ! enddo

    write(ascunit,'(A)') "---------------------------------------"
    write(ascunit,'(A, I5)') "   isym   syms%tnp(1:3,isym) "
    do isym=1,syms%ntran
       temp_r3(:) = syms%tnp(1:3, isym)
       WRITE (ascunit,'(I5, 3G20.10)') isym, temp_r3(:)
    enddo

    write(ascunit,'(A)') "======================================="

    write(ascunit,'(A)') "====== crystal:crys information ======"

    write(ascunit,'(A, I5)') "crys%nat = ", crys%nat
    write(ascunit,'(A, G20.10)') "crys%celvol = ", crys%celvol
    write(ascunit,'(A, G20.10)') "crys%alat = ", crys%alat

    write(ascunit,'(A)') 'Lattice vectors crys%avec(1:3,1:3) : '
    do ii=1, 3
       temp_r3(:) = crys%avec(ii, 1:3)
       write(ascunit,'(F15.8,1X,F15.8,1X,F15.8,1X)') temp_r3(:)
    enddo

    write(ascunit,'(A)') '------------------------------------'

    write(ascunit,'(A)') 'crys%adot : '

    do ii=1, 3
       temp_r3(:) = crys%adot(ii, 1:3)
       write(ascunit,'(F15.8,1X,F15.8,1X,F15.8,1X)') temp_r3(:)
    enddo

    write(ascunit,'(A)') '------------------------------------'

    write(ascunit,'(A, G20.10)') 'crys%recvol = ', crys%recvol
    write(ascunit,'(A, G20.10)')   'crys%blat = ', crys%blat

    write(ascunit,'(A)') 'Reciprocal lattice vectors crys%bvec(1:3,1:3): '
    do ii=1, 3
       temp_r3(:) = crys%bvec(ii, 1:3)
       write(ascunit,'(F15.8,1X,F15.8,1X,F15.8,1X)') temp_r3(:)
    enddo
    write(ascunit,'(A)') '------------------------------------'

    write(ascunit,'(A)') 'crys%bdot : '
    do ii = 1, 3
       temp_r3(:) = crys%bdot(ii, 1:3)
       write(ascunit,'(F15.8,1X,F15.8,1X,F15.8,1X)') temp_r3(:)
    enddo
    write(ascunit,'(A)') '------------------------------------'
    write(ascunit,'(A)') "   iatom   crys%atyp(iatom)   crys%apos(1:3, iatom)" ! Cartesian positions of atoms in units of alat
    do iatom=1, crys%nat
       write(ascunit,'(I5, 3X, I5, 3G20.10)') iatom, crys%atyp(iatom), crys%apos(1:3,iatom)
    enddo

    write(ascunit,'(A)') "====================================="

    close(ascunit)
  end subroutine output_header

  subroutine generate_full_WS_BZ(kgrid, qshift, bdot, nfk, fk)
    integer, intent(in) :: kgrid(3)
    real(DP), intent(in) :: qshift(3) !> in reciprocal space crystal coordinates [b1,b2,b3]
    real(DP), intent(in) :: bdot(3,3)
    integer, intent(out) :: nfk
    real(DP), allocatable, intent(out) :: fk(:,:)
    integer :: ifk, gumk(3), i1, i2, i3

    nfk = PRODUCT(kgrid(:))
    if (nfk .le. 0) then
       call die("kgrid mistake.", only_root_writes=.true.)
    endif
    SAFE_ALLOCATE(fk, (3,nfk))

    !> Generate kmesh + qshift points
    ifk = 0
    do i1 = 0, kgrid(1)-1
       do i2 = 0, kgrid(2)-1
          do i3 = 0, kgrid(3)-1
             ifk = ifk + 1
             !> kpoint in crystal coordinates [0,1]
             fk(1,ifk) = dble(i1) / dble(kgrid(1)) + qshift(1)
             fk(2,ifk) = dble(i2) / dble(kgrid(2)) + qshift(2)
             fk(3,ifk) = dble(i3) / dble(kgrid(3)) + qshift(3)
             !> Move kpoint to Wigner-Seitz BZ
             call get_gumk3(bdot, fk(:,ifk), gumk)
             !> move kpoint to 1st WS BZ
             fk(:,ifk) = fk(:,ifk) - dble(gumk(:))
          enddo
       enddo
    enddo

    if (ifk .ne. nfk) then
       call die("nfk mistach.", only_root_writes=.true.)
    endif
  end subroutine generate_full_WS_BZ

end module misc_m

!>===================================================================
!!
!! Module misc_m
!!
!! Routines:
!!
!! 1. checknorm()   By Meng Wu (2020)
!!
!!    Check the normalization condition of Bloch waves for one kpoint, one band
!!
!! 2. get_volume()       Originally By (SIB)     Last Modified 6/12/2008 (JRD)
!!
!!    This assumes that b is a symmetric matrix.  It sets
!!    vol = (2pi)^3 / square_root(|det(b)|)
!!    This makes sense if b is the matrix of dot products of the recip
!!    lattice vectors, so vol is the real space volume.
!!
!! 3. findvector()       Originally By (SIB)     Last Modified 6/12/2008 (JRD)
!!
!!    Given Gvectors as a 3D vector, return its index in the gvec%components(:,:) array
!!    If not found, return 0
!!
!! 4. M33INV()
!!
!!    Inverse of a 3x3 matrix
!!
!! 5. procmem()          Originally By (gsm)     Last Modified 4/14/2009 (gsm)
!!
!!    Determines the amount of free memory per processor
!!    from the proc file system
!!
!! 6. sizeof_scalar()    Originally By (DAS)     Last Modified 1/25/2011 (DAS)
!!
!!    Return the size of the SCALAR type, for memory estimation.
!!
!! 7. get_gumk3()       By Meng Wu (2020)
!!
!!    Given a kpoint (could be outside 1st BZ), return its Umklapp vector
!!
!! 8. get_gumk3()       By Meng Wu (2020)
!!
!!    Given a kpoint (could be outside 1st BZ), return its Umklapp vector
!!
!! 9. generate_full_WS_BZ()   By Meng Wu (2020)
!!
!!    Given kgrid, generate all the kpoints in the Wigner-Seitz BZ
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
       get_gumk3, M33INV, generate_full_WS_BZ, prepare_syms
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

    if(m == 0 .and. mac_m == 0) then ! BSD
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
#else
    y=x
    !    ytot=xtot
#endif

    !----------------------------------------------
    ! Determine the number of processors per node
    ! FHJ: we do this in parallel. Each MPI rank compares with all other hostnames,
    ! and figures out the number of ranks running on the same node. At the end, we
    ! pick the upper bound for ranks_per_node.

    SAFE_ALLOCATE(hostnames, (peinf%npes))
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

  subroutine prepare_syms(crys, syms)
    type(crystal), intent(in) :: crys
    type(symmetry), intent(inout) :: syms
    real(DP), dimension(3,3) :: atTat
    real(DP), dimension(3,3) :: inv_atTat
    real(DP), dimension(3,3) :: inv_at
    integer :: itran
    logical :: flag_inv

    ! crys%avec(:,:) = at(:,:)
    ! at(:,:) = [a1, a2, a3]
    atTat = MATMUL(TRANSPOSE(crys%avec),crys%avec)
    ! -------------------------------
    ! [a1]
    ! ------adot = [a1 a2 a3] \cdot [a2]
    ! [a3]
    ! -------------------------------
    call M33INV(atTat, inv_atTat, flag_inv)
    if (.not. flag_inv) then
       call die("prepare_syms: atTat not invertable", only_root_writes=.true.)
    endif
    call M33INV(crys%avec,inv_at, flag_inv)
    if (.not. flag_inv) then
       call die("prepare_syms: crys%avec not invertable", only_root_writes=.true.)
    endif
    syms%mtrx_reci(:,:,:) = 0.0D0
    syms%mtrx_cart(:,:,:) = 0.0D0
    do itran = 1, syms%ntran
       syms%mtrx_reci(1:3,1:3,itran) = nint(MATMUL(MATMUL(atTat, dble(syms%mtrx(1:3,1:3,itran))), inv_atTat))
       syms%mtrx_cart(1:3,1:3,itran) = MATMUL(MATMUL(crys%avec, dble(syms%mtrx(1:3,1:3,itran))), inv_at)
    enddo

  end subroutine prepare_syms

end module misc_m

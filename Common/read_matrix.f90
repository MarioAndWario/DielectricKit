#include "f_defs.h"

!=================================================================================
!
! Module read_matrix
!
! (1) read_matrix_d()           Originally by JRD       Last Modified 5/1/2008 (JRD)
!
! This program reads a distributed matrix like chimat or epsmat to file.
!
! (2) read_matrix_f()           Originally by JRD       Last Modified 9/10/2010 (gsm)
!
! Modification of read_matrix_d for full-frequency.
!
!=================================================================================

module read_matrix_m
  use global_m
  use scalapack_m
  use epsread_hdf5_m
  use hdf5
  implicit none
  private
  public :: read_matrix_d_hdf5, read_matrix_d_hdf5_2, read_matrix_f_hdf5 !, read_matrix_f_hdf5_2
  ! read_matrix_d, &
  ! read_matrix_d_hdf5, read_matrix_d_hdf5_2, &
  ! read_matrix_f, &
  ! read_matrix_f_hdf5

contains

  ! subroutine read_matrix_d(scal,matrix,nmtx,iunit)
  !   type (scalapack), intent(in) :: scal
  !   SCALAR, intent(out) :: matrix(:,:) !< (scal%npr,scal%npc)
  !   integer, intent(in) :: nmtx
  !   integer, intent(in) :: iunit
  !   PUSH_SUB(read_matrix_d)

  !   call read_matrix_d_(scal,matrix,nmtx,iunit=iunit)

  !   POP_SUB(read_matrix_d)
  ! end subroutine read_matrix_d

  subroutine read_matrix_d_hdf5(scal,matrix,nmtx,fname,iq,is)
    type (scalapack), intent(in) :: scal
    SCALAR, intent(out) :: matrix(:,:) !< (scal%npr,scal%npc)
    integer, intent(in) :: nmtx
    character(len=*), intent(in) :: fname
    integer, intent(in) :: iq
    integer, intent(in) :: is
    PUSH_SUB(read_matrix_d_hdf5)

    call read_matrix_d_(scal,matrix,nmtx,fname=fname,iq=iq,is=is)

    POP_SUB(read_matrix_d_hdf5)
  end subroutine read_matrix_d_hdf5

  subroutine read_matrix_d_(scal,matrix,nmtx,iunit,fname,iq,is)
    type (scalapack), intent(in) :: scal
    SCALAR, intent(out) :: matrix(:,:) !< (scal%npr,scal%npc)
    integer, intent(in) :: nmtx
    integer, intent(in), optional :: iunit
    character(len=*), intent(in), optional :: fname
    integer, intent(in), optional :: iq
    integer, intent(in), optional :: is
    integer :: ii, jj
#ifdef USESCALAPACK
    SCALAR, allocatable :: tempcol(:)
    integer :: irow, icol, irowm, icolm
    integer :: icurr
#endif
    logical :: use_hdf5
    PUSH_SUB(read_matrix_d_)

    if (.not.present(iunit).and..not.(present(fname).and.present(iq))) then
       call die("Not enough arguments to read_matrix_d_", only_root_writes=.true.)
    endif
    if (present(iunit).and.(present(fname).or.present(iq))) then
       call die("Too many arguments to read_matrix_d_", only_root_writes=.true.)
    endif
    if ((present(fname).or.present(iq)).and..not.(present(fname).and.present(iq))) then
       call die("Inconsistent arguments to read_matrix_d_", only_root_writes=.true.)
    endif
    use_hdf5 = present(fname).and.present(iq)
    if (.not. use_hdf5) then
       call die("must use hdf5.", only_root_writes=.true.)
    endif
    if (peinf%verb_debug .and. peinf%inode==0) then
       write(6,*) 'Reading matrix: ', nmtx, fname
       write(6,*)
    endif

#ifdef USESCALAPACK
    SAFE_ALLOCATE(tempcol, (nmtx))
    icurr=0
    do jj = 1, nmtx
       if (peinf%inode .eq. 0) then
          call read_eps_matrix_col_hdf5(tempcol,jj,nmtx,iq,is,fname)
       endif
       call MPI_BCAST(tempcol,nmtx,MPI_SCALAR,0, MPI_COMM_WORLD,mpierr)
       icol = MOD(INT(((jj-1)/scal%nbl)+TOL_SMALL),scal%npcol)
       if (icol .eq. scal%mypcol) then
          do ii = 1, nmtx
             irow = MOD(INT(((ii-1)/scal%nbl)+TOL_SMALL),scal%nprow)
             if (irow .eq. scal%myprow) then
                icurr = icurr+1
                icolm = INT((icurr-1)/scal%npr+TOL_SMALL)+1
                irowm = MOD((icurr-1),scal%npr)+1
                matrix(irowm,icolm) = tempcol(ii)
             endif
          enddo
       endif
       call MPI_barrier(MPI_COMM_WORLD,mpierr)
    enddo
    SAFE_DEALLOCATE(tempcol)
#else
    if (peinf%inode .eq. 0) then
       do jj = 1, nmtx
          call read_eps_matrix_col_hdf5(matrix(:,jj),jj,nmtx,iq,is,fname)
       enddo
    endif
#endif

    POP_SUB(read_matrix_d_)
    return
  end subroutine read_matrix_d_

  !> used to calculate 0.5 * (chi1 + chi2)
  subroutine read_matrix_d_hdf5_2(scal,matrix,nmtx,fname,iq,is)
    type (scalapack), intent(in) :: scal
    SCALAR, intent(out) :: matrix(:,:) !< (scal%npr,scal%npc)
    integer, intent(in) :: nmtx
    character(len=*), intent(in) :: fname
    integer, intent(in) :: iq
    integer, intent(in) :: is
    PUSH_SUB(read_matrix_d_hdf5_2)

    call read_matrix_d_2_(scal,matrix,nmtx,fname=fname,iq=iq,is=is)

    POP_SUB(read_matrix_d_hdf5_2)
  end subroutine read_matrix_d_hdf5_2

  subroutine read_matrix_d_2_(scal,matrix,nmtx,iunit,fname,iq,is)
    type (scalapack), intent(in) :: scal
    SCALAR, intent(out) :: matrix(:,:) !< (scal%npr,scal%npc)
    integer, intent(in) :: nmtx
    integer, intent(in), optional :: iunit
    character(len=*), intent(in), optional :: fname
    integer, intent(in), optional :: iq
    integer, intent(in), optional :: is
    integer :: ii, jj
#ifdef USESCALAPACK
    SCALAR, allocatable :: tempcol(:)
    integer :: irow, icol, irowm, icolm
    integer :: icurr
#endif
    logical :: use_hdf5
    PUSH_SUB(read_matrix_d_2_)

    if (.not.present(iunit).and..not.(present(fname).and.present(iq))) then
       call die("Not enough arguments to read_matrix_d_", only_root_writes=.true.)
    endif
    if (present(iunit).and.(present(fname).or.present(iq))) then
       call die("Too many arguments to read_matrix_d_", only_root_writes=.true.)
    endif
    if ((present(fname).or.present(iq)).and..not.(present(fname).and.present(iq))) then
       call die("Inconsistent arguments to read_matrix_d_", only_root_writes=.true.)
    endif
    use_hdf5 = present(fname).and.present(iq)
    if (.not. use_hdf5) then
       call die("must use hdf5.", only_root_writes=.true.)
    endif
    if (peinf%verb_debug .and. peinf%inode .eq. 0) then
       write(6,*) 'Reading matrix: ', nmtx, fname
       write(6,*)
    endif

#ifdef USESCALAPACK
    SAFE_ALLOCATE(tempcol, (nmtx))
    icurr=0
    do jj = 1, nmtx
       if (peinf%inode .eq. 0) then
          call read_eps_matrix_col_hdf5(tempcol,jj,nmtx,iq,is,fname)
       endif
       call MPI_BCAST(tempcol,nmtx,MPI_SCALAR,0,MPI_COMM_WORLD,mpierr)
       icol = MOD(INT(((jj-1)/scal%nbl)+TOL_SMALL),scal%npcol)
       if (icol .eq. scal%mypcol) then
          do ii = 1, nmtx
             irow=MOD(INT(((ii-1)/scal%nbl)+TOL_SMALL),scal%nprow)
             if (irow .eq. scal%myprow) then
                icurr=icurr+1
                icolm=INT((icurr-1)/scal%npr+TOL_SMALL)+1
                irowm=MOD((icurr-1),scal%npr)+1
                matrix(irowm,icolm)= 0.5 * (matrix(irowm,icolm) + tempcol(ii) )
                ! matrix(irowm,icolm) = tempcol(ii)
             endif
          enddo
       endif
       call MPI_barrier(MPI_COMM_WORLD,mpierr)
    enddo
#else
    if(peinf%inode .eq. 0) then
       do jj = 1, nmtx
          call read_eps_matrix_col_hdf5_2(matrix(:,jj),jj,nmtx,iq,is,fname)
       enddo
    endif
#endif
    SAFE_DEALLOCATE(tempcol)

    POP_SUB(read_matrix_d_2_)
    return
  end subroutine read_matrix_d_2_

  !=================================================================================

  ! !> [DEBUG HERE!!!]
  ! !> FHJ: Front end for read_matrix_f_ for Fortran binary files. See that routine for more info.
  ! subroutine read_matrix_f(scal, nfreq, nfreq_in_group, retarded, nmtx, nfreq_group, iunit, advanced)
  !   type(scalapack), intent(in) :: scal
  !   integer, intent(in) :: nfreq
  !   integer, intent(in) :: nfreq_in_group
  !   complex(DPC), intent(out) :: retarded(:,:,:) !< (scal%npr,scal%npc,pol%nfreq_in_group)
  !   integer, intent(in) :: nmtx
  !   integer, intent(in) :: nfreq_group
  !   integer, intent(in) :: iunit
  !   complex(DPC), optional, intent(out) :: advanced(:,:,:) !< (scal%npr,scal%npc,pol%nfreq_in_group)
  !   PUSH_SUB(read_matrix_f)

  !   call read_matrix_f_(scal, nfreq, nfreq_in_group, retarded, nmtx, nfreq_group, iunit=iunit, advanced=advanced)

  !   POP_SUB(read_matrix_f)
  ! end subroutine read_matrix_f

  !> FHJ: Front end for read_matrix_f_ for HDF5 files. See that routine for more info.
  subroutine read_matrix_f_hdf5(scal, nfreq, nfreq_in_group, retarded, nmtx, nfreq_group, fname, iq, is, advanced)
    type(scalapack), intent(in) :: scal
    integer, intent(in) :: nfreq
    integer, intent(in) :: nfreq_in_group
    complex(DPC), intent(out) :: retarded(:,:,:) !< (scal%npr,scal%npc,pol%nfreq_in_group)
    integer, intent(in) :: nmtx
    integer, intent(in) :: nfreq_group
    character(len=*), intent(in) :: fname
    integer, intent(in) :: iq
    integer, intent(in) :: is
    complex(DPC), optional, intent(out) :: advanced(:,:,:) !< (scal%npr,scal%npc,pol%nfreq_in_group)
    PUSH_SUB(read_matrix_f_hdf5)

    call read_matrix_f_(scal, nfreq, nfreq_in_group, retarded, nmtx, nfreq_group, fname=fname, iq=iq, is=is, advanced=advanced)

    POP_SUB(read_matrix_f_hdf5)
  end subroutine read_matrix_f_hdf5

  !> FHJ: This routines the full-frequency chiR/epsR matrix from a file, and
  !! optionally chiA/epsA (note: you shouldn`t really need chiA, ever...)
  !! If using HDF5, we only read the retarded part. If legacy
  !! Fortran binary, we read the retarded and skip the advanced. The final
  !! matrix will be distributed in a ScaLAPACK layout given by scal. Note that
  !! this routine is pretty innefficient, but this is not a core component
  !! of BGW as it`s only used if you read_chi or use the eps*omega utility.
  subroutine read_matrix_f_(scal, nfreq, nfreq_in_group, retarded, nmtx, nfreq_group, iunit, fname, iq, is, advanced)
    type(scalapack), intent(in) :: scal
    integer, intent(in) :: nfreq
    integer, intent(in) :: nfreq_in_group
    complex(DPC), intent(out) :: retarded(:,:,:) !> !< (scal%npr,scal%npc,pol%nfreq_in_group)
    integer, intent(in) :: nmtx
    integer, intent(in) :: nfreq_group
    integer, intent(in), optional :: iunit
    character(len=*), intent(in), optional :: fname
    integer, intent(in), optional :: iq
    integer, intent(in), optional :: is
    complex(DPC), intent(out), optional :: advanced(:,:,:) !> !< (scal%npr,scal%npc,pol%nfreq_in_group)
    integer :: ii, jj, ifreq,ifreq_para,freq_grp_ind
#ifdef USESCALAPACK
    complex(DPC), allocatable :: tempcolR(:,:)
    complex(DPC), allocatable :: tempcolA(:,:)
    integer :: irow, icol, irowm, icolm
    integer :: icurr
#endif
    logical :: use_hdf5, want_advanced
    PUSH_SUB(read_matrix_f_)

    want_advanced = .false.
#ifdef CPLX
    want_advanced = present(advanced)
#endif
    if (.not.present(iunit).and..not.(present(fname).and.present(iq))) then
       call die("Not enough arguments to read_matrix_f_", only_root_writes=.true.)
    endif
    if (present(iunit).and.(present(fname).or.present(iq))) then
       call die("Too many arguments to read_matrix_f_", only_root_writes=.true.)
    endif
    if ((present(fname).or.present(iq)).and..not.(present(fname).and.present(iq))) then
       call die("Inconsistent arguments to read_matrix_f_", only_root_writes=.true.)
    endif
    use_hdf5 = present(fname).and.present(iq)
    if (.not. use_hdf5) then
       call die("must use hdf5.", only_root_writes=.true.)
    endif

    if (peinf%verb_debug .and. peinf%inode==0) then
       write(6,*) 'Reading matrix: ', nmtx, fname
       write(6,*)
    endif
#ifdef USESCALAPACK

    SAFE_ALLOCATE(tempcolR, (nmtx,nfreq))
    if (want_advanced) then
       call die("read_matrix: Advanced chimat not supported", only_root_writes=.true.)
       SAFE_ALLOCATE(tempcolA, (nmtx,nfreq))
    endif

    icurr = 0
    !> Read the jj-th column of epsinv: epsinv(:,jj) at all frequencies
    do jj = 1, nmtx
       if (peinf%inode .eq. 0) then
          if (.not. want_advanced) then
             call read_eps_matrix_col_f_hdf5(tempcolR, nfreq, jj, nmtx, iq, is, fname)
          else
             call read_eps_matrix_col_f_hdf5(tempcolR, nfreq, jj, nmtx, iq, is, fname, advanced=tempcolA)
          endif
       endif
       if (peinf%npes > 1) then
          call MPI_BCAST(tempcolR, nfreq*nmtx, MPI_COMPLEX_DPC, 0, MPI_COMM_WORLD, mpierr)
          if (want_advanced) then
             call MPI_BCAST(tempcolA, nfreq*nmtx, MPI_COMPLEX_DPC, 0, MPI_COMM_WORLD, mpierr)
          endif
       endif
       icol = MOD(INT(((jj-1)/scal%nbl)+TOL_SMALL),scal%npcol)
       if (icol .eq. scal%mypcol) then
          do ii = 1, nmtx
             irow = MOD(INT(((ii-1)/scal%nbl)+TOL_SMALL),scal%nprow)
             if (irow .eq. scal%myprow) then
                icurr = icurr+1
                icolm = INT((icurr-1)/scal%npr+TOL_SMALL)+1
                irowm = MOD((icurr-1),scal%npr)+1
                do ifreq = 1, nfreq
                   freq_grp_ind = MOD(ifreq-1,nfreq_group)
                   ifreq_para = (ifreq+nfreq_group-1)/nfreq_group
                   if (freq_grp_ind .eq. peinf%igroup_f) then
                      retarded(irowm,icolm,ifreq_para) = tempcolR(ii,ifreq)
                      if (want_advanced) then
                         !> [BUGS HERE]
                         !> comment this line will lead to no segmentation fault!
                         !> probably because advanced is optional (in fact not used in epsilon.x), but here we just assigne a value to it!
                         ! advanced(irowm,icolm,ifreq_para) = tempcolA(ii,ifreq)
                      endif
                   endif
                enddo
             endif
          enddo
       endif
       call MPI_barrier(MPI_COMM_WORLD,mpierr)
    enddo

    SAFE_DEALLOCATE(tempcolR)
    if (want_advanced) then
       SAFE_DEALLOCATE(tempcolA)
    endif
#else
    !> not USESCALAPACK
    if (peinf%inode .eq. 0) then
       do jj = 1, nmtx
          if (.not. want_advanced) then
             call read_eps_matrix_col_f_hdf5(retarded(:,:,jj), nfreq, jj, nmtx, iq, is, fname)
          else
             call read_eps_matrix_col_f_hdf5(retarded(:,:,jj), nfreq, jj, nmtx, iq, is, fname, advanced=advanced(:,:,jj))
          endif
       enddo
    endif
#endif

    POP_SUB(read_matrix_f_)
    return
  end subroutine read_matrix_f_

  !   !> [WORKING]
  !   !> nfreq_in_group = nfreq
  !   !> nfreq_group = 1
  !   subroutine read_matrix_f_hdf5_2(scal, nfreq, retarded, nmtx, fname, iq, is)
  !     type(scalapack), intent(in) :: scal
  !     integer, intent(in) :: nfreq
  !     complex(DPC), intent(out) :: retarded(:,:,:) !> !< (scal%npr,scal%npc,pol%nfreq)
  !     integer, intent(in) :: nmtx
  !     character(len=*), intent(in) :: fname
  !     integer, intent(in) :: iq
  !     integer, intent(in) :: is

  !     integer :: ii, jj, ifreq,ifreq_para
  ! #ifdef USESCALAPACK
  !     complex(DPC), allocatable :: tempcolR(:,:)
  !     complex(DPC), allocatable :: tempcolA(:,:)
  !     integer :: irow, icol, irowm, icolm
  !     integer :: icurr
  ! #endif
  !     logical :: use_hdf5, want_advanced
  !     PUSH_SUB(read_matrix_f_)

  !     want_advanced = .false.
  ! #ifdef CPLX
  !     want_advanced = present(advanced)
  ! #endif

  !     if ((present(fname).or.present(iq)).and..not.(present(fname).and.present(iq))) then
  !        call die("Inconsistent arguments to read_matrix_f_", only_root_writes=.true.)
  !     endif
  !     use_hdf5 = present(fname).and.present(iq)
  !     if (.not. use_hdf5) then
  !        call die("must use hdf5.", only_root_writes=.true.)
  !     endif

  !     if (peinf%verb_debug .and. peinf%inode==0) then
  !        write(6,*) 'Reading matrix: ', nmtx, fname
  !        write(6,*)
  !     endif
  ! #ifdef USESCALAPACK

  !     SAFE_ALLOCATE(tempcolR, (nmtx, nfreq))
  !     if (want_advanced) then
  !        call die("read_matrix: Advanced chimat not supported", only_root_writes=.true.)
  !     endif

  !     icurr = 0
  !     !> Read the ig2-th column of epsinv: epsinv(:,ig2) at all frequencies
  !     do ig2 = 1, nmtx
  !        if (peinf%inode .eq. 0) then
  !           !> [WORKING]
  !           !> Read ig2-th column epsinv(ig1=1:nmtx,ig2,1:nfreq)
  !           call read_eps_matrix_col_f_hdf5(tempcolR, nfreq, ig2, nmtx, iq, is, fname)
  !        endif
  !        if (peinf%npes > 1) then
  !           call MPI_BCAST(tempcolR, nfreq*nmtx, MPI_COMPLEX_DPC, 0, MPI_COMM_WORLD, mpierr)
  !        endif

  !        target_mypcol = INDXG2P(ig2, scal%nbl, scal%mypcol, 0, scal%npcol)
  !        if (target_mypcol .ne. scal%mypcol) cycle

  !        do ig1 = 1, nmtx
  !           target_myprow = INDXG2P(ig1, scal%nbl, scal%myprow, 0, scal%nprow)
  !           if (target_myprow .ne. scal%myprow) cycle

  !           ig1_loc = INDXG2L(ig1, scal%nbl, scal%myprow, 0, scal%nprow)
  !           ig2_loc = INDXG2L(ig2, scal%nbl, scal%mypcol, 0, scal%npcol)
  !           do ifreq = 1, nfreq
  !              retarded(ig1_loc, ig2_loc, ifreq) = tempcolR(ig1, ifreq)
  !           enddo

  !        enddo

  !        call MPI_barrier(MPI_COMM_WORLD,mpierr)
  !     enddo

  !     SAFE_DEALLOCATE(tempcolR)
  ! #else
  !     !> not USESCALAPACK
  !     if (peinf%inode .eq. 0) then
  !        do jj = 1, nmtx
  !           !> [WORKING]
  !           !> [BUGS HERE]
  !           call read_eps_matrix_col_f_hdf5(retarded(:,:,jj), nfreq, jj, nmtx, iq, is, fname)
  !        enddo
  !     endif
  ! #endif

  !     POP_SUB(read_matrix_f_hdf5_2)
  !     return
  !   end subroutine read_matrix_f_hdf5_2

!   subroutine read_matrix_d_par_hdf5(scal, matrix, nmtx, iq, is, name)
!     type(scalapack), intent(in) :: scal
!     SCALAR, intent(out) :: matrix(:,:) !< (scal%npr,scal%npc)
!     integer, intent(in) :: nmtx
!     integer, intent(in) :: iq
!     integer, intent(in) :: is
!     character(len=*), intent(in) :: name
!     integer :: ii, jj, error, rank
!     real(DP), allocatable :: data(:,:,:,:,:,:)
!     integer(HID_T) :: file_id       ! File identifier
!     integer(HID_T) :: dset_id       ! Dataset identifier
!     integer(HID_T) :: filespace     ! Dataspace identifier in file
!     integer(HID_T) :: memspace      ! Dataspace identifier in mem
!     integer(HID_T) :: plist_id      ! Property list identifier for parallel IO
!     integer(HSIZE_T) :: count(6), countm(6), offset(6), offsetm(6)
!     integer :: info
!     PUSH_SUB(read_matrix_d_par_hdf5)

!     if (peinf%inode.eq.0) call timacc(47,1)

!     ! JRD: We need a barrier here or else parallel file opening gets mixed up with
!     ! peinf%inode 0 opening the file to write the diagonal (which is called first).
!     call MPI_barrier(MPI_COMM_WORLD, mpierr)
!     info = MPI_INFO_NULL

!     rank = 6
!     !> hyperslab in memory
!     countm(1) = SCALARSIZE
!     countm(2) = nmtx
!     countm(3) = 
!     countm(4) = 1
!     countm(5) = 1
!     countm(6) = 1
!     offsetm(:) = 0

!     !> Determine hyperslab in file
!     count(1) = 1
!     count(2) = scal%npr / scal%nbl
!     count(3) = scal%npc / scal%nbl
!     count(4) = 1
!     count(5) = 1
!     count(6) = 1

!     offset(1) = 0
!     offset(2) = scal%myprow * scal%nbr
!     offset(3) = scal%mypcol * scal%nbc
!     offset(4) = 0
!     offset(5) = is - 1
!     offset(6) = iq - 1

!     stride(1) = 1
!     stride(2) = scal%nprow * scal%nbr
!     stride(3) = scal%npcol * scal%nbc
!     stride(4) = 1
!     stride(5) = 1
!     stride(6) = 1

!     call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
!     call h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, info, error)

!     call h5fopen_f(trim(name), H5F_ACC_RDWR_F, file_id, error, access_prp = plist_id)
!     call h5pclose_f(plist_id,error)

!     SAFE_ALLOCATE(data, (countm(1),countm(2),countm(3),countm(4),countm(5),countm(6)))
    
! !     do jj = 1, scal%npc
! !        !do jj = 1, scal%npc - mod(scal%npc,scal%nbl)
! !        do ii = 1, scal%npr
! !           !do ii = 1, scal%npr - mod(scal%npr,scal%nbl)
! !           data(1,ii,jj,1,1,1) = dble(matrix(ii,jj))
! ! #ifdef CPLX
! !           data(2,ii,jj,1,1,1) = IMAG(matrix(ii,jj))
! ! #endif
! !        enddo
! !     enddo

!     call h5screate_simple_f(rank, countm, memspace, error)
!     call h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, offsetm, countm, error)

!     call h5dopen_f(file_id, 'mats/matrix', dset_id, error)
!     call h5dget_space_f(dset_id, filespace, error)

!     call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, count, error, stride, block_)

!     ! Add in remainders, in case scal%nbl doesnt perfectly divide nmtx
!     ! Bottom Rows
!     rowremainder = MOD(scal%npr, scal%nbr)
!     if (rowremainder .ne. 0) then
!        offsetr = offset
!        countr = count
!        blockr = block_
!        strider = stride
!        offsetr(2) = nmtx-rowremainder
!        countr(2) = rowremainder
!        blockr(2) = 1
!        strider(2) = 1
!        call h5sselect_hyperslab_f(filespace, H5S_SELECT_OR_F, offsetr, countr, error, strider, blockr)
!        !write(6,*) peinf%inode, "I have the bottom row", rowremainder, scal%npc
!     endif

!     ! Right Columns
!     colremainder = MOD(scal%npc, scal%nbc)
!     if (colremainder .ne. 0) then
!        offsetr=offset
!        countr=count
!        blockr=block_
!        strider=stride
!        offsetr(3)=nmtx-colremainder
!        countr(3)=colremainder
!        blockr(3)=1
!        strider(3)=1
!        call h5sselect_hyperslab_f(filespace, H5S_SELECT_OR_F, offsetr, countr, error, strider, blockr)
!        !write(6,*) peinf%inode, "I have the right column", colremainder, scal%npr
!        ! Bottom Corner of Matrix
!        if (rowremainder .ne. 0) then
!           offsetr=offset
!           countr=count
!           blockr=block_
!           strider=stride
!           offsetr(2)=nmtx-rowremainder
!           countr(2)=rowremainder
!           blockr(2)=1
!           strider(2)=1
!           offsetr(3)=nmtx-colremainder
!           countr(3)=colremainder
!           blockr(3)=1
!           strider(3)=1
!           call h5sselect_hyperslab_f(filespace, H5S_SELECT_OR_F, offsetr, countr, error, strider, blockr)
!           !write(6,*) peinf%inode, "I have bottom both"
!        endif
!     endif

!     call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
!     call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
!     if(peinf%inode.eq.0) call timacc(47,2)
!     if(peinf%inode.eq.0) call timacc(48,1)
!     call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, data, countm, error, memspace, filespace, xfer_prp = plist_id)
!     if(peinf%inode.eq.0) call timacc(48,2)
!     call h5pclose_f(plist_id, error)

!     SAFE_DEALLOCATE(data)
!     call h5dclose_f(dset_id, error)
!     call h5sclose_f(memspace, error)
!     call h5sclose_f(filespace, error)
!     call h5fclose_f(file_id, error)

!     POP_SUB(read_matrix_d_par_hdf5)
!     return
!   end subroutine write_matrix_d_par_hdf

end module read_matrix_m

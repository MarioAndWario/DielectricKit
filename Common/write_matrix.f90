!=================================================================================
!
! Module write_matrix_m
!
! (1) write_matrix_d()          Originally by JRD       Last Modified 5/1/2008 (JRD)
!
! This program writes a distributed matrix like chimat or epsmat to file.
!
! (2) write_matrix_f()          Originally by JRD       Last Modified 2/5/2009 (CHP)
!
! Modification of write_matrix_d for full-frequency.
!
!=================================================================================

#include "f_defs.h"

module write_matrix_m
  use global_m
#ifdef HDF5
  use hdf5
#endif
  use hdf5_io_m
  use scalapack_m
  use io_utils_m
  implicit none
  private
  ! public :: write_matrix_ser_hdf, write_matrix_f_ser_hdf, write_matrix_diagonal_hdf, write_gvec_indices_hdf, write_pol_fftgrid_hdf ! write_matrix_f_hdf, write_matrix_d_hdf, write_vcoul_hdf
  public :: write_matrix_ser_hdf, write_matrix_f_ser_hdf
  public :: write_gvec_indices_hdf, write_pol_fftgrid_hdf
#ifdef USESCALAPACK
  public :: write_matrix_d_par_hdf_2, write_matrix_f_par_hdf_2 ! write_matrix_d_par_hdf, write_matrix_f_par_hdf
#endif
  ! write_matrix_f
  ! write_matrix_d, &
  ! write_matrix_d_sub, &
contains

  !> [WORKING]
  !> Collect nc_write columns of the chimat matrix and output to chimat.h5 using ROOT
  subroutine write_matrix_ser_hdf(scal, matrix, nmtx_max, nc_write, ifreq, iq, is, filename)
    type(scalapack), intent(in) :: scal
    SCALAR, intent(in) :: matrix(:,:) !< (scal%npr,scal%npc)
    integer, intent(in) :: nmtx_max, nc_write, ifreq, iq, is
    character(len=*), intent(in) :: filename
    integer :: error, info
    INTEGER(HID_T) :: file_id       ! File identifier
    integer(HID_T) :: dset_id
    INTEGER(HID_T) :: filespace
    integer(HID_T) :: memspace
    integer :: rank_mat = 6
    integer(HSIZE_T) :: count_mat(6), offset_mat(6)
    integer :: desc_root(9), desc_2d(9)
    integer :: cntxt_root, info_blacs
    integer :: nprow, npcol, myprow, mypcol
    integer :: npr, npc, lld_root, lld_2d
    integer, parameter :: ctxt_ = 2
    SCALAR, allocatable :: temp_c(:,:)
    real(DP), allocatable :: temp_r(:,:,:)
    integer :: nc_write_, icol_start !, alloc_stat
    type(progress_info) :: prog_info
    PUSH_SUB(write_matrix_ser_hdf)    

    !> scal%nbr = ICEIL(nmtx_max, scal%nprow)
    !> scal%nbc = ICEIL(nmtx_max, scal%npcol)
    !> Use p?gemr2d to copy a submatrix from one general rectangular matrix to another
    !> Each time, collect nc_write (default, scal%nbr) columns with p?gemr2d to ROOT and then output it

    !> Initialize BLACS grid for 2d block-cyclic distributed matrix matrix(nmtx_max, nmtx_max)
    lld_2d = MAX(scal%npr,1)
    if ( scal%myprow .ne. -1) then
       call DESCINIT(desc_2d, nmtx_max, nmtx_max, scal%nbr, scal%nbc, 0, 0, scal%icntxt, lld_2d, info_blacs)
    else
       desc_2d(ctxt_) = -1
       info_blacs = 0
    endif
    if (info_blacs .ne. 0) then
       call die('DESCINIT failed for pol%chi', only_root_writes=.true.)
    endif

    !> Initialize BLACS grid for temp_c with dimension (nmtx_max, nc_write) stored on ROOT
    !> Blocksize=nmtx_max
    call blacs_get(0, 0, cntxt_root)
    !> This context only has ROOT
    call blacs_gridinit(cntxt_root, 'R', 1, 1)
    call blacs_gridinfo(cntxt_root, nprow, npcol, myprow, mypcol)
    if ( (myprow .ge. 1) .or. (mypcol .ge. 1) ) then
       myprow = -1;
       mypcol = -1;
    endif
    if (myprow .ne. -1) then
       npr = NUMROC(nmtx_max, nmtx_max, myprow, 0, nprow)
       npc = NUMROC(nc_write, nc_write, mypcol, 0, npcol)
       if (npr .ne. nmtx_max) then
          call die("Distribution failed 1.", only_root_writes=.true.)
       endif
       if (npc .ne. nc_write) then
          call die("Distribution failed 2.", only_root_writes=.true.)
       endif
    else
       npr = 0
       npc = 0
    endif
    lld_root = MAX(1, npr)

    if (myprow .ne. -1) then
       if (lld_root .ne. nmtx_max) then
          call die("lld_evecs should be nmtx_max.", only_root_writes=.true.)
       endif
       call DESCINIT(desc_root, nmtx_max, nc_write, nmtx_max, nc_write, 0, 0, cntxt_root, lld_root, info)
    else
       desc_root(ctxt_) = -1
       info = 0
    endif
    if (info .ne. 0) then
       call die('DESCINIT failed for tempcol', only_root_writes=.true.)
    endif

    if (myprow .ne. -1) then
       SAFE_ALLOCATE(temp_c, (npr, npc))
    else
       SAFE_ALLOCATE(temp_c, (1, 1))
    endif
    temp_c = ZERO

    !> ROOT outputs [nmtx_max, nc_write_] elements to file
    if (myprow .ne. -1) then
       ! write(*,*) "AAA peinf%inode = ", peinf%inode
       !> Open file for serial HDF5 write
       call h5fopen_f(TRUNC(filename), H5F_ACC_RDWR_F, file_id, error)
       if (error .ne. 0) then
          call die("HDF5 error", only_root_writes=.true.)
       endif
       call h5dopen_f(file_id, 'mats/matrix', dset_id, error)
       if (error .ne. 0) then
          call die("HDF5 error", only_root_writes=.true.)
       endif
       call h5dget_space_f(dset_id, filespace, error)
       if (error .ne. 0) then
          call die("HDF5 error", only_root_writes=.true.)
       endif
    endif

    call progress_init(prog_info, 'Output to '//TRUNC(filename), 'column-blocks', (nmtx_max+nc_write-1)/nc_write)    

    !> All procs in the context of scal%icntxt should communicate with ROOT
    if (scal%myprow .ne. -1) then
       !> Loop over column blocks
       do icol_start = 1, nmtx_max, nc_write
          call progress_step(prog_info, (icol_start+nc_write-1)/nc_write)
          nc_write_ = MIN(nc_write, nmtx_max - icol_start + 1)
          ! write(*,'(A,I5,A,I5)') "icol_start = ", icol_start, " nc_write_ = ", nc_write_
          call pX(gemr2d)(nmtx_max, nc_write_, matrix, 1, icol_start, desc_2d, temp_c, 1, 1, desc_root, scal%icntxt)

          if (myprow .ne. -1) then
             ! write(*,*) "AAA peinf%inode = ", peinf%inode
             SAFE_ALLOCATE(temp_r, (2, npr, npc))
             temp_r(1,:,:) = DBLE(temp_c(:,:))
             if (SCALARSIZE .eq. 2) then
                temp_r(2,:,:) = DIMAG(temp_c(:,:))
             endif

             !> Write mats/matrix
             if (npc > 0) then
                !> output one frequency at a time
                count_mat(:) = (/ SCALARSIZE, nmtx_max, nc_write_, 1, 1, 1 /)
                offset_mat(:) = (/ 0, 0, icol_start-1, ifreq-1, is - 1, iq - 1 /)
             else
                count_mat = 0
                offset_mat = 0
             endif

             call h5screate_simple_f(rank_mat, count_mat, memspace, error)
             if (error .ne. 0) then
                call die("HDF5 error", only_root_writes=.true.)
             endif

             call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset_mat, count_mat, error)
             if (error .ne. 0) then
                call die("HDF5 error", only_root_writes=.true.)
             endif

             call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, temp_r, count_mat, error, mem_space_id = memspace, file_space_id = filespace)
             if (error .ne. 0) then
                call die("HDF5 error", only_root_writes=.true.)
             endif

             call h5sclose_f(memspace, error)
             SAFE_DEALLOCATE(temp_r)
          endif
       enddo
    endif

    if (myprow .ne. -1) then
       call h5sclose_f(filespace, error)
       call h5dclose_f(dset_id, error)
       call h5fclose_f(file_id, error)
    endif

    SAFE_DEALLOCATE(temp_c)

    if (myprow .ne. -1) then
       CALL BLACS_GRIDEXIT(cntxt_root)
    endif

    POP_SUB(write_matrix_ser_hdf)
    return
  end subroutine write_matrix_ser_hdf

  subroutine write_matrix_f_ser_hdf(scal, matrix, nmtx_max, nfreq, nc_write, iq, is, filename)
    type(scalapack), intent(in) :: scal
    SCALAR, intent(in) :: matrix(:,:,:) !< (scal%npr,scal%npc, nfreq)
    integer, intent(in) :: nmtx_max, nfreq, nc_write, iq, is
    character(len=*), intent(in) :: filename
    integer :: error, info
    INTEGER(HID_T) :: file_id       ! File identifier
    integer(HID_T) :: dset_id
    INTEGER(HID_T) :: filespace
    integer(HID_T) :: memspace
    integer :: rank_mat = 6
    integer(HSIZE_T) :: count_mat(6), offset_mat(6)
    integer :: desc_root(9), desc_2d(9)
    integer :: cntxt_root, info_blacs
    integer :: nprow, npcol, myprow, mypcol
    integer :: npr, npc, lld_root, lld_2d
    integer, parameter :: ctxt_ = 2
    SCALAR, allocatable :: temp_c(:,:)
    real(DP), allocatable :: temp_r(:,:,:)
    integer :: nc_write_, icol_start, ifreq
    type(progress_info) :: prog_info    
    PUSH_SUB(write_matrix_f_ser_hdf)

    !> scal%nbr = ICEIL(nmtx_max, scal%nprow)
    !> scal%nbc = ICEIL(nmtx_max, scal%npcol)
    !> Use p?gemr2d to copy a submatrix from one general rectangular matrix to another
    !> Each time, collect nc_write (default, scal%nbr) columns with p?gemr2d to ROOT and then output it

    !> Initialize BLACS grid for 2d block-cyclic distributed matrix matrix(1:nmtx_max, 1:nmtx_max, ifreq)
    lld_2d = MAX(scal%npr,1)
    if ( scal%myprow .ne. -1) then
       call DESCINIT(desc_2d, nmtx_max, nmtx_max, scal%nbr, scal%nbc, 0, 0, scal%icntxt, lld_2d, info_blacs)
    else
       desc_2d(ctxt_) = -1
       info_blacs = 0
    endif
    if (info_blacs .ne. 0) then
       call die('DESCINIT failed for pol%chi', only_root_writes=.true.)
    endif

    !> Initialize BLACS grid for temp_c with dimension (nmtx_max, nc_write) stored on ROOT
    !> Blocksize=nmtx_max
    call blacs_get(0, 0, cntxt_root)
    !> This context only has ROOT
    call blacs_gridinit(cntxt_root, 'R', 1, 1)
    call blacs_gridinfo(cntxt_root, nprow, npcol, myprow, mypcol)
    if ( (myprow .ge. 1) .or. (mypcol .ge. 1) ) then
       myprow = -1;
       mypcol = -1;
    endif
    if (myprow .ne. -1) then
       npr = NUMROC(nmtx_max, nmtx_max, myprow, 0, nprow)
       npc = NUMROC(nc_write, nc_write, mypcol, 0, npcol)
       if (npr .ne. nmtx_max) then
          call die("Distribution failed 1.", only_root_writes=.true.)
       endif
       if (npc .ne. nc_write) then
          call die("Distribution failed 2.", only_root_writes=.true.)
       endif
    else
       npr = 0
       npc = 0
    endif
    lld_root = MAX(1, npr)

    if (myprow .ne. -1) then
       if (lld_root .ne. nmtx_max) then
          call die("lld_evecs should be nmtx_max.", only_root_writes=.true.)
       endif
       call DESCINIT(desc_root, nmtx_max, nc_write, nmtx_max, nc_write, 0, 0, cntxt_root, lld_root, info)
    else
       desc_root(ctxt_) = -1
       info = 0
    endif
    if (info .ne. 0) then
       call die('DESCINIT failed for tempcol', only_root_writes=.true.)
    endif

    if (myprow .ne. -1) then
       SAFE_ALLOCATE(temp_c, (npr, npc))
    else
       SAFE_ALLOCATE(temp_c, (1, 1))
    endif
    temp_c = ZERO

    !> ROOT outputs [nmtx_max, nc_write_] elements to file
    if (myprow .ne. -1) then
       ! write(*,*) "AAA peinf%inode = ", peinf%inode
       !> Open file for serial HDF5 write
       call h5fopen_f(TRUNC(filename), H5F_ACC_RDWR_F, file_id, error)
       if (error .ne. 0) then
          call die("HDF5 error", only_root_writes=.true.)
       endif
       call h5dopen_f(file_id, 'mats/matrix', dset_id, error)
       if (error .ne. 0) then
          call die("HDF5 error", only_root_writes=.true.)
       endif
       call h5dget_space_f(dset_id, filespace, error)
       if (error .ne. 0) then
          call die("HDF5 error", only_root_writes=.true.)
       endif
    endif

    call progress_init(prog_info, 'Output to '//TRUNC(filename), 'column-blocks', nfreq*((nmtx_max+nc_write-1)/nc_write))    
    
    !> All procs in the context of scal%icntxt should communicate with ROOT
    if (scal%myprow .ne. -1) then
       do ifreq = 1, nfreq
          !> Loop over column blocks
          do icol_start = 1, nmtx_max, nc_write
             call progress_step(prog_info, (icol_start+nc_write-1)/nc_write + (ifreq-1) * ((nmtx_max+nc_write-1)/nc_write) )
             nc_write_ = MIN(nc_write, nmtx_max - icol_start + 1)
             ! write(*,'(A,I5,A,I5)') "icol_start = ", icol_start, " nc_write_ = ", nc_write_
             call pX(gemr2d)(nmtx_max, nc_write_, matrix(:,:,ifreq), 1, icol_start, desc_2d, temp_c, 1, 1, desc_root, scal%icntxt)

             if (myprow .ne. -1) then
                ! write(*,*) "AAA peinf%inode = ", peinf%inode
                SAFE_ALLOCATE(temp_r, (2, npr, npc))
                temp_r(1,:,:) = DBLE(temp_c(:,:))
                if (SCALARSIZE .eq. 2) then
                   temp_r(2,:,:) = DIMAG(temp_c(:,:))
                endif

                !> Write mats/matrix
                if (npc > 0) then
                   !> output one frequency at a time
                   count_mat(:) = (/ SCALARSIZE, nmtx_max, nc_write_, 1, 1, 1 /)
                   offset_mat(:) = (/ 0, 0, icol_start-1, ifreq-1, is - 1, iq - 1 /)
                else
                   count_mat = 0
                   offset_mat = 0
                endif

                call h5screate_simple_f(rank_mat, count_mat, memspace, error)
                if (error .ne. 0) then
                   call die("HDF5 error", only_root_writes=.true.)
                endif

                call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset_mat, count_mat, error)
                if (error .ne. 0) then
                   call die("HDF5 error", only_root_writes=.true.)
                endif

                call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, temp_r, count_mat, error, mem_space_id = memspace, file_space_id = filespace)
                if (error .ne. 0) then
                   call die("HDF5 error", only_root_writes=.true.)
                endif

                call h5sclose_f(memspace, error)
                SAFE_DEALLOCATE(temp_r)
             endif
          enddo ! icol_start
       enddo ! ifreq
    endif

    if (myprow .ne. -1) then
       call h5sclose_f(filespace, error)
       call h5dclose_f(dset_id, error)
       call h5fclose_f(file_id, error)
    endif

    SAFE_DEALLOCATE(temp_c)

    if (myprow .ne. -1) then
       CALL BLACS_GRIDEXIT(cntxt_root)
    endif

    POP_SUB(write_matrix_f_ser_hdf)
    return
  end subroutine write_matrix_f_ser_hdf

  !===================================================================================

  !   subroutine write_matrix_d(scal,matrix,nmtx,iunit)
  !     type(scalapack), intent(in) :: scal
  !     SCALAR, intent(in) :: matrix(:,:) !< (scal%npr,scal%npc)
  !     integer, intent(in) :: nmtx
  !     integer, intent(in) :: iunit

  !     integer :: ii, jj
  ! #ifdef USESCALAPACK
  !     SCALAR, allocatable :: tempcol(:),tempcol2(:)
  !     integer :: irow, icol, irowm, icolm
  !     integer :: icurr
  ! #endif
  !     type(progress_info) :: prog_info !< a user-friendly progress report

  !     PUSH_SUB(write_matrix_d)

  !     if (peinf%verb_debug .and. peinf%inode==0) then
  !        write(6,*) 'Writing matrix: ', nmtx, iunit
  !        write(6,*)
  !     endif

  ! #ifdef USESCALAPACK
  !     SAFE_ALLOCATE(tempcol, (nmtx))
  !     SAFE_ALLOCATE(tempcol2, (nmtx))

  !     icurr=0

  !     call progress_init(prog_info, 'writing matrix', 'column', nmtx)
  !     do jj = 1, nmtx
  !        call progress_step(prog_info, jj)
  !        !        if (peinf%inode .eq. 0) then
  !        !          write(6,*) ' In loop: ', ii
  !        !        endif
  !        icol=MOD(INT(((jj-1)/scal%nbl)+TOL_SMALL),scal%npcol)
  !        tempcol=0d0
  !        if (icol .eq. scal%mypcol) then
  !           do ii = 1, nmtx
  !              irow=MOD(INT(((ii-1)/scal%nbl)+TOL_SMALL),scal%nprow)
  !              if (irow .eq. scal%myprow) then
  !                 icurr=icurr+1
  !                 icolm=INT((icurr-1)/scal%npr+TOL_SMALL)+1
  !                 irowm=MOD((icurr-1),scal%npr)+1
  !                 tempcol(ii)=matrix(irowm,icolm)

  !                 !                if (icolm .gt. scal%npc .or. irowm.gt.scal%npr) then
  !                 !                  write(6,*) 'Error: ', scal%npc,scal%npr,icolm,irowm
  !                 !                endif

  !              endif
  !           enddo
  !        endif
  !        if (peinf%inode .eq. 0) then
  !           tempcol2=0d0
  !        endif
  !        call MPI_REDUCE(tempcol,tempcol2,nmtx,MPI_SCALAR,MPI_SUM,0, &
  !             MPI_COMM_WORLD,mpierr)
  !        if (peinf%inode .eq. 0) then
  !           write(iunit) (tempcol2(ii),ii=1,nmtx)
  !        endif

  !        call MPI_barrier(MPI_COMM_WORLD,mpierr)

  !     enddo
  !     call progress_free(prog_info)

  !     SAFE_DEALLOCATE(tempcol)
  !     SAFE_DEALLOCATE(tempcol2)

  ! #else
  !     if (peinf%inode .eq. 0) then
  !        call progress_init(prog_info, 'writing matrix', 'column', nmtx)
  !        do jj = 1, nmtx
  !           call progress_step(prog_info, jj)
  !           write(iunit) (matrix(ii, jj), ii = 1, nmtx)
  !        enddo
  !        call progress_free(prog_info)
  !     endif
  ! #endif

  !     POP_SUB(write_matrix_d)
  !     return
  !   end subroutine write_matrix_d

  !   subroutine write_matrix_d_sub(scal,matrix,nmtx,iunit,neig_sub)
  !     type(scalapack), intent(in) :: scal
  !     complex(DPC), intent(in) :: matrix(:,:) !< (scal%npr,scal%npc)
  !     integer, intent(in) :: nmtx
  !     integer, intent(in) :: iunit
  !     integer, intent(in), optional :: neig_sub
  !     integer :: ii, jj, nmtx_col
  ! #ifdef USESCALAPACK
  !     complex(DPC), allocatable :: tempcol(:),tempcol2(:)
  !     integer :: irow, icol, irowm, icolm
  !     integer :: icurr
  ! #endif
  !     type(progress_info) :: prog_info !< a user-friendly progress report
  !     PUSH_SUB(write_matrix_d_sub)

  !     if (peinf%verb_debug .and. peinf%inode==0) then
  !        write(6,*) 'Writing matrix: ', nmtx, iunit
  !        write(6,*)
  !     endif

  !     ! neig_sub allows to write only neig_sub columns of the actual matrix
  !     nmtx_col = nmtx
  !     IF(PRESENT(neig_sub)) nmtx_col = neig_sub

  ! #ifdef USESCALAPACK
  !     SAFE_ALLOCATE(tempcol, (nmtx))
  !     SAFE_ALLOCATE(tempcol2, (nmtx))

  !     icurr=0

  !     call progress_init(prog_info, 'writing matrix', 'column', nmtx_col)
  !     do jj = 1, nmtx_col
  !        call progress_step(prog_info, jj)
  !        !        if (peinf%inode .eq. 0) then
  !        !          write(6,*) ' In loop: ', ii
  !        !        endif
  !        icol=MOD(INT(((jj-1)/scal%nbl)+TOL_SMALL),scal%npcol)
  !        tempcol=0d0
  !        if (icol .eq. scal%mypcol) then
  !           do ii = 1, nmtx
  !              irow=MOD(INT(((ii-1)/scal%nbl)+TOL_SMALL),scal%nprow)
  !              if (irow .eq. scal%myprow) then
  !                 icurr=icurr+1
  !                 icolm=INT((icurr-1)/scal%npr+TOL_SMALL)+1
  !                 irowm=MOD((icurr-1),scal%npr)+1
  !                 tempcol(ii)=matrix(irowm,icolm)

  !                 !                if (icolm .gt. scal%npc .or. irowm.gt.scal%npr) then
  !                 !                  write(6,*) 'Error: ', scal%npc,scal%npr,icolm,irowm
  !                 !                endif

  !              endif
  !           enddo
  !        endif
  !        if (peinf%inode .eq. 0) then
  !           tempcol2=0d0
  !        endif
  !        call MPI_REDUCE(tempcol,tempcol2,nmtx,MPI_COMPLEX_DPC,MPI_SUM,0, &
  !             MPI_COMM_WORLD,mpierr)
  !        if (peinf%inode .eq. 0) then
  !           write(iunit) (tempcol2(ii),ii=1,nmtx)
  !        endif

  !        call MPI_barrier(MPI_COMM_WORLD,mpierr)

  !     enddo
  !     call progress_free(prog_info)

  !     if (peinf%inode .eq. 0) then
  !        ! write empty rows for the missin column so sigma will work also with previously
  !        ! generated epsmat (this can go in the future)
  !        do jj = nmtx_col + 1, nmtx
  !           write(iunit)
  !        end do
  !     end if

  !     SAFE_DEALLOCATE(tempcol)
  !     SAFE_DEALLOCATE(tempcol2)

  ! #else
  !     if (peinf%inode .eq. 0) then
  !        call progress_init(prog_info, 'writing matrix', 'column', nmtx_col)
  !        do jj = 1, nmtx_col
  !           call progress_step(prog_info, jj)
  !           write(iunit) (matrix(ii, jj), ii = 1, nmtx)
  !        enddo
  !        call progress_free(prog_info)
  !        do jj = nmtx_col + 1, nmtx
  !           write(iunit)
  !        end do
  !     endif
  ! #endif

  !     POP_SUB(write_matrix_d_sub)
  !     return
  !   end subroutine write_matrix_d_sub

  !   !=================================================================================
  !> The indices of tempcolR and tempcolA are wrong!

  !   subroutine write_matrix_f(scal,nfreq,retarded,nmtx,iunit,nfreq_group,advanced)
  !     type(scalapack), intent(in) :: scal
  !     integer, intent(in) :: nfreq
  !     complex(DPC), intent(in) :: retarded(:,:,:) !< (scal%npr,scal%npc,nfreq_in_group)
  !     integer, intent(in) :: nmtx
  !     integer, intent(in) :: iunit
  !     integer, intent(in) :: nfreq_group
  !     complex(DPC), optional, intent(in) :: advanced(:,:,:) !< (scal%npr,scal%npc,nfreq_in_group)

  !     integer :: ii, jj, ifreq
  ! #ifdef USESCALAPACK
  !     complex(DPC), allocatable :: tempcolR(:,:),tempcolR2(:,:)
  !     complex(DPC), allocatable :: tempcolA(:,:),tempcolA2(:,:)
  !     integer :: irow, icol, irowm, icolm,freq_grp_ind,ifreq_para
  !     integer :: icurr
  ! #endif
  !     type(progress_info) :: prog_info !< a user-friendly progress report
  !     logical :: has_advanced

  !     PUSH_SUB(write_matrix_f)

  !     if (peinf%verb_debug .and. peinf%inode==0) then
  !        write(6,*) 'Writing matrix: ', nfreq, nmtx, iunit
  !        write(6,*)
  !     endif

  !     has_advanced = present(advanced)
  ! #ifdef USESCALAPACK
  !     SAFE_ALLOCATE(tempcolR, (nfreq,nmtx))
  !     SAFE_ALLOCATE(tempcolR2, (nfreq,nmtx))
  ! #ifdef CPLX
  !     SAFE_ALLOCATE(tempcolA, (nfreq,nmtx))
  !     SAFE_ALLOCATE(tempcolA2, (nfreq,nmtx))
  ! #endif

  !     icurr=0

  !     call progress_init(prog_info, 'writing matrix', 'column', nmtx)
  !     do jj = 1, nmtx
  !        call progress_step(prog_info, jj)
  !        !        if (peinf%inode .eq. 0) then
  !        !          write(6,*) ' In loop: ', ii
  !        !        endif
  !        icol=MOD(INT(((jj-1)/scal%nbl)+TOL_SMALL),scal%npcol)
  !        tempcolR=0d0
  ! #ifdef CPLX
  !        tempcolA=0d0
  ! #endif
  !        if (icol .eq. scal%mypcol) then
  !           do ii = 1, nmtx
  !              irow=MOD(INT(((ii-1)/scal%nbl)+TOL_SMALL),scal%nprow)
  !              if (irow .eq. scal%myprow) then
  !                 icurr=icurr+1
  !                 icolm=INT((icurr-1)/scal%npr+TOL_SMALL)+1
  !                 irowm=MOD((icurr-1),scal%npr)+1
  !                 do ifreq=1,nfreq
  !                    freq_grp_ind=mod(ifreq-1,nfreq_group)
  !                    ifreq_para=(ifreq+nfreq_group-1)/nfreq_group
  !                    if(freq_grp_ind .eq. peinf%rank_mtxel) then
  !                       tempcolR(ifreq,ii)=retarded(irowm,icolm,ifreq_para)
  ! #ifdef CPLX
  !                       if (has_advanced) then
  !                          tempcolA(ifreq,ii)=advanced(irowm,icolm,ifreq_para)
  !                       endif
  ! #endif
  !                    endif
  !                 enddo
  !              endif
  !           enddo
  !        endif
  !        if (peinf%inode .eq. 0) then
  !           tempcolR2=0d0
  ! #ifdef CPLX
  !           if (has_advanced) then
  !              tempcolA2=0d0
  !           endif
  ! #endif
  !        endif
  !        call MPI_REDUCE(tempcolR(1,1),tempcolR2(1,1),nfreq*nmtx,MPI_COMPLEX_DPC,MPI_SUM,0,MPI_COMM_WORLD,mpierr)
  ! #ifdef CPLX
  !        if (has_advanced) then
  !           call MPI_REDUCE(tempcolA(1,1),tempcolA2(1,1),nfreq*nmtx,MPI_COMPLEX_DPC,MPI_SUM,0,MPI_COMM_WORLD,mpierr)
  !        endif
  ! #endif
  !        if (peinf%inode .eq. 0) then
  !           do ii = 1, nmtx
  !              write(iunit) (tempcolR2(ifreq,ii),ifreq=1,nfreq)
  !           enddo
  ! #ifdef CPLX
  !           if (has_advanced) then
  !              do ii = 1, nmtx
  !                 write(iunit) (tempcolA2(ifreq,ii),ifreq=1,nfreq)
  !              enddo
  !           else
  !              do ii = 1, nmtx
  !                 write(iunit)
  !              enddo
  !           endif
  ! #endif
  !        endif

  !        call MPI_barrier(MPI_COMM_WORLD,mpierr)

  !     enddo
  !     call progress_free(prog_info)

  !     SAFE_DEALLOCATE(tempcolR)
  !     SAFE_DEALLOCATE(tempcolR2)
  ! #ifdef CPLX
  !     SAFE_DEALLOCATE(tempcolA)
  !     SAFE_DEALLOCATE(tempcolA2)
  ! #endif

  ! #else

  !     if(peinf%inode .eq. 0) then
  !        call progress_init(prog_info, 'writing matrix', 'column', nmtx)
  !        do jj = 1, nmtx
  !           call progress_step(prog_info, jj)
  !           do ii = 1, nmtx
  !              write(iunit) (retarded(ii, jj, ifreq), ifreq= 1, nfreq)
  !           enddo
  ! #ifdef CPLX
  !           if (has_advanced) then
  !              do ii = 1, nmtx
  !                 write(iunit) (advanced(ii, jj, ifreq),ifreq = 1, nfreq)
  !              enddo
  !           else
  !              do ii = 1, nmtx
  !                 write(iunit)
  !              enddo
  !           endif
  ! #endif
  !        enddo
  !        call progress_free(prog_info)
  !     endif

  ! #endif

  !     POP_SUB(write_matrix_f)

  !     return
  !   end subroutine write_matrix_f

#ifdef HDF5

  !========================================================================
  ! JRD: The HDF5 Equivalents of the above routines.
  !========================================================================

  !==========================================================================================

  ! subroutine write_matrix_diagonal_hdf(epsdiag,nmtx,iq,isize,name)
  !   real(DP), intent(in) :: epsdiag(:,:,:) !< (isize,nmtx,1)
  !   integer, intent(in) :: nmtx
  !   integer, intent(in) :: iq
  !   integer, intent(in) :: isize
  !   character(len=*), intent(in) :: name

  !   integer(HID_T) :: file_id       ! File identifier
  !   integer(HID_T) :: dset_id       ! Dataset identifier
  !   integer(HID_T) :: filespace     ! Dataspace identifier in file
  !   integer(HID_T) :: memspace      ! Dataspace identifier in mem

  !   integer(HSIZE_T) :: dims(3), offset(3), offsetm(3)

  !   integer :: error, rank

  !   PUSH_SUB(write_matrix_diagonal_hdf)

  !   call open_file(99, trim(name), status='old')
  !   call close_file(99)

  !   call h5fopen_f(trim(name), H5F_ACC_RDWR_F, file_id, error)

  ! ! Write Array

  !   rank = 3
  !   dims(1) = isize
  !   dims(2) = nmtx
  !   dims(3) = 1
  !   offset(1) = 0
  !   offset(2) = 0
  !   offset(3) = iq - 1
  !   offsetm(:) = 0

  !   call h5dopen_f(file_id, 'mats/matrix-diagonal', dset_id, error)
  !   call h5screate_simple_f(rank, dims, memspace, error)
  !   call h5dget_space_f(dset_id,filespace,error)
  !   call h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, offsetm, dims, error)
  !   call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, dims, error)
  !   call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, epsdiag, dims, error, memspace, filespace)
  !   call h5dclose_f(dset_id, error)
  !   call h5sclose_f(memspace, error)
  !   call h5sclose_f(filespace, error)

  !   call h5fclose_f(file_id, error)

  !   POP_SUB(write_matrix_diagonal_hdf)

  ! end subroutine write_matrix_diagonal_hdf

  ! subroutine write_matrix_diagonal_hdf(epsdiag,nmtx,nfreq,iq,isize,name)
  !   real(DP), intent(in) :: epsdiag(:,:,:) !< (isize,nmtx,nfreq) for one rq (iq-th rq)
  !   integer, intent(in) :: nmtx, nfreq
  !   integer, intent(in) :: iq
  !   integer, intent(in) :: isize
  !   character(len=*), intent(in) :: name

  !   ! integer(HSIZE_T) :: countf(4), offsetf(4)
  !   integer :: countf(4), offsetf(4)
  !   integer(HID_T) :: file_id       ! File identifier

  !   integer :: error

  !   PUSH_SUB(write_matrix_diagonal_hdf)

  !   call h5fopen_f(trim(name), H5F_ACC_RDWR_F, file_id, error)
  !   if (error .ne. 0) then
  !      call die("HDF5 error", only_root_writes=.true.)
  !   endif

  !   countf(:) = (/ isize, nmtx, nfreq, 1/)
  !   offsetf(:) = (/0, 0, 0, iq-1/)

  !   call hdf5_write_double_hyperslab(file_id, 'mats/matrix-diagonal', countf, offsetf, epsdiag, error)

  !   if (error .ne. 0) then
  !      call die("HDF5 error", only_root_writes=.true.)
  !   endif

  !   call h5fclose_f(file_id, error)
  !   if (error .ne. 0) then
  !      call die("HDF5 error", only_root_writes=.true.)
  !   endif

  !   POP_SUB(write_matrix_diagonal_hdf)

  ! end subroutine write_matrix_diagonal_hdf


  !===================================================================================

  !   subroutine write_matrix_ser_hdf(matrix,nmtx,iq,is,name)
  !     SCALAR, intent(in) :: matrix(:,:) !< (nmtx,nmtx)
  !     integer, intent(in) :: nmtx
  !     integer, intent(in) :: iq
  !     integer, intent(in) :: is
  !     character(len=*), intent(in) :: name

  !     integer :: error, rank, ii, jj

  !     integer(HID_T) :: file_id       ! File identifier
  !     integer(HID_T) :: dset_id       ! Dataset identifier
  !     integer(HID_T) :: filespace     ! Dataspace identifier in file
  !     integer(HID_T) :: memspace      ! Dataspace identifier in mem

  !     integer(HSIZE_T) :: count(6), offset(6), offsetm(6)

  !     real(DP), allocatable :: data(:,:,:,:,:,:)

  !     PUSH_SUB(write_matrix_ser_hdf)

  !     call open_file(99, trim(name), status='old')
  !     call close_file(99)

  !     call h5fopen_f(trim(name), H5F_ACC_RDWR_F, file_id, error)

  !     rank=6
  !     count(1) = SCALARSIZE
  !     count(2) = nmtx
  !     count(3) = nmtx
  !     count(4) = 1
  !     count(5) = 1
  !     count(6) = 1

  !     offset(:) = 0
  !     offset(5) = is - 1
  !     offset(6) = iq - 1

  !     offsetm(:) = 0

  !     SAFE_ALLOCATE(data,(count(1),count(2),count(3),count(4),count(5),count(6)))

  !     do jj = 1, nmtx
  !        do ii = 1, nmtx
  !           data(1,ii,jj,1,1,1) = dble(matrix(ii,jj))
  ! #ifdef CPLX
  !           data(2,ii,jj,1,1,1) = IMAG(matrix(ii,jj))
  ! #endif
  !        enddo
  !     enddo

  !     call h5dopen_f(file_id, 'mats/matrix', dset_id, error)
  !     call h5screate_simple_f(rank, count, memspace, error)
  !     call h5dget_space_f(dset_id,filespace,error)
  !     call h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, offsetm, count, error)
  !     call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, count, error)
  !     call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, data, count, error, memspace, filespace)

  !     SAFE_DEALLOCATE(data)

  !     call h5dclose_f(dset_id, error)
  !     call h5sclose_f(memspace, error)
  !     call h5sclose_f(filespace, error)

  !     call h5fclose_f(file_id, error)

  !     POP_SUB(write_matrix_ser_hdf)

  !   end subroutine write_matrix_ser_hdf

  !   !========================================================================

  !   subroutine write_matrix_f_ser_hdf(nfreq, retarded, nmtx, iq, is, name)
  !     integer, intent(in) :: nfreq
  !     complex(DPC), intent(in) :: retarded(:,:,:) !< (nmtx,nmtx,nfreq)
  !     integer, intent(in) :: nmtx
  !     integer, intent(in) :: iq
  !     integer, intent(in) :: is
  !     character(len=*), intent(in) :: name
  !     integer :: jj, error, rank
  !     real(DP), allocatable :: data(:,:,:,:,:,:)
  !     type(progress_info) :: prog_info !< a user-friendly progress report
  !     integer(HID_T) :: file_id       ! File identifier
  !     integer(HID_T) :: dset_id       ! Dataset identifier
  !     integer(HID_T) :: filespace     ! Dataspace identifier in file
  !     integer(HID_T) :: memspace      ! Dataspace identifier in mem
  !     integer(HSIZE_T) :: count(6), offset(6), offsetm(6)
  !     PUSH_SUB(write_matrix_f_ser_hdf)

  !     ! DVF: this routine was built off of write_matrix_f_hdf to do the serial
  !     ! writing of an hdf format matrix. This is needed for epsmat_old2hdf5.f90

  !     rank=6
  !     count(1) = 2
  !     count(2) = nmtx
  !     count(3) = 1
  !     count(4) = nfreq
  !     count(5) = 1
  !     count(6) = 1

  !     SAFE_ALLOCATE(data, (count(1),count(2),count(3),count(4),count(5),count(6)))

  !     call open_file(99, trim(name), status='old')
  !     call close_file(99)

  !     call h5fopen_f(trim(name), H5F_ACC_RDWR_F, file_id, error)
  !     call h5dopen_f(file_id, 'mats/matrix', dset_id, error)
  !     call h5screate_simple_f(rank, count, memspace, error)
  !     call h5dget_space_f(dset_id,filespace,error)

  !     call progress_init(prog_info, 'writing matrix', 'column', nmtx)
  !     ! JRD XXX the fact that jj is not outer loop presents a bit of challenge
  !     ! but this serial routine is just for legacy support
  !     do jj = 1, nmtx
  !        call progress_step(prog_info, jj)
  !        data(1,:,1,:,1,1)=dble(retarded(:,jj,:))
  !        data(2,:,1,:,1,1)=IMAG(retarded(:,jj,:))

  !        offsetm(:) = 0
  !        call h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, offsetm, count, error)

  !        offset(1)=0
  !        offset(2)=0
  !        offset(3)=jj-1
  !        offset(4)=0
  !        offset(5)=is-1
  !        offset(6)=iq-1

  !        call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, count, error)
  !        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, data, count, error, memspace, filespace)
  !     enddo
  !     call progress_free(prog_info)

  !     SAFE_DEALLOCATE(data)
  !     call h5dclose_f(dset_id, error)
  !     call h5sclose_f(memspace, error)
  !     call h5sclose_f(filespace, error)
  !     call h5fclose_f(file_id, error)

  !     POP_SUB(write_matrix_f_ser_hdf)

  !     return

  !   end subroutine write_matrix_f_ser_hdf

  !========================================================================

  !   subroutine write_matrix_d_hdf(scal,matrix,nmtx,iq,is,name)
  !     type(scalapack), intent(in) :: scal
  !     SCALAR, intent(in) :: matrix(:,:) !< (scal%npr,scal%npc)
  !     integer, intent(in) :: nmtx
  !     integer, intent(in) :: iq
  !     integer, intent(in) :: is
  !     character(len=*), intent(in) :: name
  !     integer :: ii, jj, error, size, rank
  ! #ifdef USESCALAPACK
  !     real(DP), allocatable :: datatmp(:,:,:,:,:,:)
  !     integer :: irow, icol, irowm, icolm
  !     integer :: icurr
  ! #endif
  !     real(DP), allocatable :: data(:,:,:,:,:,:)
  !     type(progress_info) :: prog_info !< a user-friendly progress report
  !     integer(HID_T) :: file_id       ! File identifier
  !     integer(HID_T) :: dset_id       ! Dataset identifier
  !     integer(HID_T) :: filespace     ! Dataspace identifier in file
  !     integer(HID_T) :: memspace      ! Dataspace identifier in mem
  ! #ifdef USESCALAPACK
  !     !  integer(HID_T) :: plist_id      ! Property list identifier for parallel IO
  !     !                                 ! Not used yet...
  ! #endif

  !     integer(HSIZE_T) :: count(6), offset(6), offsetm(6)

  !     PUSH_SUB(write_matrix_d_hdf)

  !     if (peinf%verb_debug .and. peinf%inode==0) then
  !        write(6,*) 'Writing matrix: ', nmtx
  !        write(6,*)
  !     endif

  !     ! XXX: For now, we will still have only proc 0 write...
  !     ! We should changes this to parallel writes. But doing
  !     ! this effectively from the scalapack, block cyclic layout
  !     ! seems a bit tricky. So, ignoring for now...

  !     rank=6
  !     count(1) = SCALARSIZE
  !     count(2) = nmtx
  !     count(3) = 1
  !     count(4) = 1
  !     count(5) = 1
  !     count(6) = 1

  !     if (peinf%inode .eq. 0) then
  !        SAFE_ALLOCATE(data,(count(1),count(2),count(3),count(4),count(5),count(6)))

  !        call open_file(99, trim(name), status='old')
  !        call close_file(99)

  !        call h5fopen_f(trim(name), H5F_ACC_RDWR_F, file_id, error)
  !        call h5dopen_f(file_id, 'mats/matrix', dset_id, error)
  !        call h5screate_simple_f(rank, count, memspace, error)
  !        call h5dget_space_f(dset_id,filespace,error)
  !     endif

  ! #ifdef USESCALAPACK
  !     SAFE_ALLOCATE(datatmp, (count(1),count(2),count(3),count(4),count(5),count(6)))
  !     icurr=0
  ! #endif

  !     call progress_init(prog_info, 'writing matrix', 'column', nmtx)
  !     do jj = 1, nmtx

  !        call progress_step(prog_info, jj)
  ! #ifdef USESCALAPACK

  !        if(peinf%inode.eq.0) call timacc(47,1)

  !        icol=MOD(INT(((jj-1)/scal%nbl)+TOL_SMALL),scal%npcol)
  !        datatmp=0d0
  !        if (icol .eq. scal%mypcol) then
  !           do ii = 1, nmtx
  !              irow=MOD(INT(((ii-1)/scal%nbl)+TOL_SMALL),scal%nprow)
  !              if (irow .eq. scal%myprow) then
  !                 icurr=icurr+1
  !                 icolm=INT((icurr-1)/scal%npr+TOL_SMALL)+1
  !                 irowm=MOD((icurr-1),scal%npr)+1
  !                 datatmp(1,ii,1,1,1,1)=dble(matrix(irowm,icolm))
  ! #ifdef CPLX
  !                 datatmp(2,ii,1,1,1,1)=IMAG(matrix(irowm,icolm))
  ! #endif
  !              endif
  !           enddo
  !        endif
  !        if (peinf%inode .eq. 0) then
  !           data=0d0
  !        endif

  !        ! XXX This is a big waste of communication. Should be fixed when do
  !        ! parallel IO.

  !        size = nmtx * SCALARSIZE

  !        call MPI_REDUCE(datatmp,data,size,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
  !             MPI_COMM_WORLD,mpierr)

  !        if(peinf%inode.eq.0) call timacc(47,2)

  ! #else

  !        if (peinf%inode .eq. 0) then
  !           do ii = 1, nmtx
  !              data(1,ii,1,1,1,1) = dble(matrix(ii,jj))
  ! #ifdef CPLX
  !              data(2,ii,1,1,1,1) = IMAG(matrix(ii,jj))
  ! #endif
  !           enddo
  !        endif

  ! #endif

  !        if(peinf%inode.eq.0) call timacc(48,1)

  !        if (peinf%inode .eq. 0) then

  !           offsetm(:) = 0
  !           call h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, offsetm, count, error)

  !           offset(1)=0
  !           offset(2)=0
  !           offset(3)=jj-1
  !           offset(4)=0
  !           offset(5)=is-1
  !           offset(6)=iq-1

  !           call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, count, error)

  !           call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, data, count, error, memspace, filespace)

  !        endif

  !        if(peinf%inode.eq.0) call timacc(48,2)

  ! #ifdef USESCALAPACK
  !        call MPI_barrier(MPI_COMM_WORLD,mpierr)
  ! #endif

  !     enddo
  !     call progress_free(prog_info)

  ! #ifdef USESCALAPACK
  !     SAFE_DEALLOCATE(datatmp)
  ! #endif

  !     if (peinf%inode .eq. 0) then
  !        SAFE_DEALLOCATE(data)
  !        call h5dclose_f(dset_id, error)
  !        call h5sclose_f(memspace, error)
  !        call h5sclose_f(filespace, error)
  !        call h5fclose_f(file_id, error)
  !     endif

  !     POP_SUB(write_matrix_d_hdf)

  !     return
  !   end subroutine write_matrix_d_hdf

  !========================================================================

#ifdef USESCALAPACK

  !   subroutine write_matrix_d_par_hdf(scal, matrix, nmtx, iq, is, name)
  !     type(scalapack), intent(in) :: scal
  !     SCALAR, intent(in) :: matrix(:,:) !< (scal%npr,scal%npc)
  !     integer, intent(in) :: nmtx, iq, is
  !     character(len=*), intent(in) :: name
  !     integer :: ii, jj, error, rank
  !     real(DP), allocatable :: data(:,:,:,:,:,:)
  !     integer(HID_T) :: file_id       ! File identifier
  !     integer(HID_T) :: dset_id       ! Dataset identifier
  !     integer(HID_T) :: filespace     ! Dataspace identifier in file
  !     integer(HID_T) :: memspace      ! Dataspace identifier in mem
  !     integer(HID_T) :: plist_id      ! Property list identifier for parallel IO
  !     integer(HSIZE_T) :: count(6), countm(6), offset(6), offsetm(6), stride(6), block_(6)
  !     integer(HSIZE_T) :: countr(6), offsetr(6), strider(6), blockr(6)
  !     integer :: comm, info, rowremainder, colremainder
  !     PUSH_SUB(write_matrix_d_par_hdf)

  !     if (peinf%inode .eq. 0) call timacc(47,1)

  !     ! JRD: We need a barrier here or else parallel file opening gets mixed up with
  !     ! peinf%inode 0 opening the file to write the diagonal (which is called first).
  !     call MPI_barrier(MPI_COMM_WORLD,mpierr)

  !     comm = MPI_COMM_WORLD
  !     info = MPI_INFO_NULL

  !     if (peinf%verb_debug .and. peinf%inode==0) then
  !        write(6,*) 'Writing matrix: ', nmtx
  !        write(6,*)
  !     endif

  !     ! JRD Should be ok with npr and npc = 0
  !     !if (scal%npr .eq. 0 .or. scal%npc .eq. 0) then
  !     !  write(6,*) peinf%inode,"Zero npr or npc!!", scal%npr, scal%npc
  !     !endif

  !     rank=6
  !     countm(1) = SCALARSIZE
  !     countm(2) = scal%npr
  !     countm(3) = scal%npc
  !     countm(4) = 1
  !     countm(5) = 1
  !     countm(6) = 1

  !     offsetm(:) = 0

  !     count(1) = 1
  !     count(2) = scal%npr/scal%nbl
  !     count(3) = scal%npc/scal%nbl
  !     count(4) = 1
  !     count(5) = 1
  !     count(6) = 1

  !     block_(1) = SCALARSIZE
  !     block_(2) = scal%nbl
  !     block_(3) = scal%nbl
  !     block_(4) = 1
  !     block_(5) = 1
  !     block_(6) = 1

  !     offset(1) = 0
  !     offset(2) = scal%myprow*scal%nbl
  !     offset(3) = scal%mypcol*scal%nbl
  !     offset(4) = 0
  !     offset(5) = is-1
  !     offset(6) = iq-1

  !     stride(1) = 1
  !     stride(2) = scal%nprow*scal%nbl
  !     stride(3) = scal%npcol*scal%nbl
  !     stride(4) = 1
  !     stride(5) = 1
  !     stride(6) = 1

  !     call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
  !     call h5pset_fapl_mpio_f(plist_id, comm, info, error)

  !     call open_file(99, trim(name), status='old')
  !     call close_file(99)

  !     call h5fopen_f(trim(name), H5F_ACC_RDWR_F, file_id, error, access_prp = plist_id)
  !     call h5pclose_f(plist_id,error)

  !     SAFE_ALLOCATE(data, (countm(1),countm(2),countm(3),countm(4),countm(5),countm(6)))
  !     do jj = 1, scal%npc
  !        do ii = 1, scal%npr
  !           data(1,ii,jj,1,1,1) = dble(matrix(ii,jj))
  ! #ifdef CPLX
  !           data(2,ii,jj,1,1,1) = IMAG(matrix(ii,jj))
  ! #endif
  !        enddo
  !     enddo

  !     call h5screate_simple_f(rank, countm, memspace, error)
  !     call h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, offsetm, countm, error)

  !     call h5dopen_f(file_id, 'mats/matrix', dset_id, error)
  !     call h5dget_space_f(dset_id,filespace,error)
  !     call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, count, error, stride, block_)

  !     ! Add in remainders, in case scal%nbl doesnt perfectly divide nmtx

  !     ! Bottom Rows
  !     rowremainder = mod(scal%npr, scal%nbl)
  !     if (rowremainder .ne. 0) then
  !        offsetr=offset
  !        countr=count
  !        blockr=block_
  !        strider=stride
  !        offsetr(2)=nmtx-rowremainder
  !        countr(2)=rowremainder
  !        blockr(2)=1
  !        strider(2)=1
  !        call h5sselect_hyperslab_f(filespace, H5S_SELECT_OR_F, offsetr, countr, error, strider, blockr)
  !        !write(6,*) peinf%inode, "I have the bottom row", rowremainder, scal%npc
  !     endif

  !     ! Right Columns
  !     colremainder = mod(scal%npc,scal%nbl)
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
  !     call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, data, countm, error, memspace, filespace, &
  !          xfer_prp = plist_id)
  !     if(peinf%inode.eq.0) call timacc(48,2)
  !     call h5pclose_f(plist_id, error)

  !     SAFE_DEALLOCATE(data)
  !     call h5dclose_f(dset_id, error)
  !     call h5sclose_f(memspace, error)
  !     call h5sclose_f(filespace, error)
  !     call h5fclose_f(file_id, error)

  !     POP_SUB(write_matrix_d_par_hdf)
  !     return
  !   end subroutine write_matrix_d_par_hdf

  !> Use nbr&nbc instead of nbl
  subroutine write_matrix_d_par_hdf_2(scal, matrix, nmtx_max, iq, is, name)
    type(scalapack), intent(in) :: scal
    SCALAR, intent(in) :: matrix(:,:) !< (scal%npr,scal%npc)
    integer, intent(in) :: nmtx_max, iq, is
    character(len=*), intent(in) :: name
    integer :: ii, jj, error, rank
    real(DP), allocatable :: data(:,:,:,:,:,:)
    integer(HID_T) :: file_id       ! File identifier
    integer(HID_T) :: dset_id       ! Dataset identifier
    integer(HID_T) :: filespace     ! Dataspace identifier in file
    integer(HID_T) :: memspace      ! Dataspace identifier in mem
    integer(HID_T) :: plist_id      ! Property list identifier for parallel IO
    integer(HSIZE_T) :: count(6), countm(6), offset(6), offsetm(6), stride(6), block_(6)
    integer(HSIZE_T) :: countr(6), offsetr(6), strider(6), blockr(6)
    integer :: comm, info, rowremainder, colremainder
    PUSH_SUB(write_matrix_d_par_hdf)

    if (peinf%inode .eq. 0) call timacc(47,1)
    ! JRD: We need a barrier here or else parallel file opening gets mixed up with
    ! peinf%inode 0 opening the file to write the diagonal (which is called first).
    call MPI_barrier(MPI_COMM_WORLD,mpierr)
    comm = MPI_COMM_WORLD
    info = MPI_INFO_NULL

    if (peinf%verb_debug .and. peinf%inode==0) then
       write(6,*) 'Writing matrix: ', nmtx_max
       write(6,*)
    endif

    ! JRD Should be ok with npr and npc = 0
    !if (scal%npr .eq. 0 .or. scal%npc .eq. 0) then
    !  write(6,*) peinf%inode,"Zero npr or npc!!", scal%npr, scal%npc
    !endif

    rank = 6
    countm(1) = SCALARSIZE
    countm(2) = scal%npr
    countm(3) = scal%npc
    countm(4) = 1
    countm(5) = 1
    countm(6) = 1

    offsetm(:) = 0

    ! https://ecpannualmeeting.com/assets/overview/sessions/HDF5-ECP-AM-2020-Tutorial.pdf
    !> If MOD[scal%npr, scal%nbr] !=0, then there are some residue elements
    !> If MOD[scal%npc, scal%nbc] !=0, then there are some residue elements
    count(1) = 1
    count(2) = scal%npr/scal%nbr
    count(3) = scal%npc/scal%nbc
    count(4) = 1
    count(5) = 1
    count(6) = 1

    block_(1) = SCALARSIZE
    block_(2) = scal%nbr
    block_(3) = scal%nbc
    block_(4) = 1
    block_(5) = 1
    block_(6) = 1

    offset(1) = 0
    offset(2) = scal%myprow*scal%nbr
    offset(3) = scal%mypcol*scal%nbc
    offset(4) = 0
    offset(5) = is-1
    offset(6) = iq-1

    stride(1) = 1
    stride(2) = scal%nprow*scal%nbr
    stride(3) = scal%npcol*scal%nbc
    stride(4) = 1
    stride(5) = 1
    stride(6) = 1

    call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif
    call h5pset_fapl_mpio_f(plist_id, comm, info, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    ! call open_file(99, trim(name), status='old')
    ! call close_file(99)

    call h5fopen_f(trim(name), H5F_ACC_RDWR_F, file_id, error, access_prp = plist_id)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif
    call h5pclose_f(plist_id,error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    SAFE_ALLOCATE(data, (countm(1),countm(2),countm(3),countm(4),countm(5),countm(6)))
    do jj = 1, scal%npc
       do ii = 1, scal%npr
          data(1,ii,jj,1,1,1) = dble(matrix(ii,jj))
#ifdef CPLX
          data(2,ii,jj,1,1,1) = IMAG(matrix(ii,jj))
#endif
       enddo
    enddo

    call h5screate_simple_f(rank, countm, memspace, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif
    call h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, offsetm, countm, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call h5dopen_f(file_id, 'mats/matrix', dset_id, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif
    call h5dget_space_f(dset_id,filespace,error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif
    call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, count, error, stride, block_)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    ! Add in remainders, in case scal%nbr doesnt perfectly divide nmtx_max

    ! Bottom Rows
    rowremainder = mod(scal%npr, scal%nbr)
    if (rowremainder .ne. 0) then
       offsetr = offset
       countr = count
       blockr = block_
       strider = stride
       offsetr(2) = nmtx_max - rowremainder
       countr(2) = rowremainder
       blockr(2) = 1
       strider(2) = 1
       call h5sselect_hyperslab_f(filespace, H5S_SELECT_OR_F, offsetr, countr, error, strider, blockr)
       if (error .ne. 0) then
          call die("HDF5 error", only_root_writes=.true.)
       endif
       !write(6,*) peinf%inode, "I have the bottom row", rowremainder, scal%npc
    endif

    ! Right Columns
    colremainder = mod(scal%npc, scal%nbc)
    if (colremainder .ne. 0) then
       offsetr = offset
       countr = count
       blockr = block_
       strider = stride
       offsetr(3) = nmtx_max - colremainder
       countr(3) = colremainder
       blockr(3) = 1
       strider(3) = 1
       call h5sselect_hyperslab_f(filespace, H5S_SELECT_OR_F, offsetr, countr, error, strider, blockr)
       if (error .ne. 0) then
          call die("HDF5 error", only_root_writes=.true.)
       endif
       !write(6,*) peinf%inode, "I have the right column", colremainder, scal%npr
       ! Bottom Corner of Matrix
       if (rowremainder .ne. 0) then
          offsetr = offset
          countr = count
          blockr = block_
          strider = stride
          offsetr(2) = nmtx_max - rowremainder
          countr(2) = rowremainder
          blockr(2) = 1
          strider(2) = 1
          offsetr(3) = nmtx_max - colremainder
          countr(3) = colremainder
          blockr(3) = 1
          strider(3) = 1
          call h5sselect_hyperslab_f(filespace, H5S_SELECT_OR_F, offsetr, countr, error, strider, blockr)
          if (error .ne. 0) then
             call die("HDF5 error", only_root_writes=.true.)
          endif
          !write(6,*) peinf%inode, "I have bottom both"
       endif
    endif

    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif
    if (peinf%inode .eq. 0) call timacc(47,2)
    if (peinf%inode .eq. 0) call timacc(48,1)
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, data, countm, error, memspace, filespace, xfer_prp = plist_id)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif
    if (peinf%inode .eq. 0) call timacc(48,2)
    call h5pclose_f(plist_id, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    SAFE_DEALLOCATE(data)
    call h5dclose_f(dset_id, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif
    call h5sclose_f(memspace, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif
    call h5sclose_f(filespace, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    ! !> Flushes the internal HDF5 buffers then asks the operating system (the OS) to flush the system buffers for the open files
    ! call h5fflush_f(file_id, H5F_SCOPE_GLOBAL_F, error)
    ! if (error .ne. 0) then
    !    call die("HDF5 error", only_root_writes=.true.)
    ! endif

    call h5fclose_f(file_id, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    POP_SUB(write_matrix_d_par_hdf_2)
    return
  end subroutine write_matrix_d_par_hdf_2

  !========================================================================================

  ! subroutine write_matrix_f_par_hdf(scal, nfreq, retarded, nmtx, iq, is, name, nfreq_group)
  !   type(scalapack), intent(in) :: scal
  !   integer, intent(in) :: nfreq_in_group
  !   complex(DPC), intent(in) :: retarded(:,:,:) !< (scal%npr,scal%npc,nfreq_in_group)
  !   integer, intent(in) :: nmtx
  !   integer, intent(in) :: iq
  !   integer, intent(in) :: is
  !   character(len=*), intent(in) :: name
  !   integer, intent(in) :: nfreq_group
  !   integer :: error, rank
  !   real(DP), allocatable :: data(:,:,:,:,:,:)
  !   integer(HID_T) :: file_id       ! File identifier
  !   integer(HID_T) :: dset_id       ! Dataset identifier
  !   integer(HID_T) :: filespace     ! Dataspace identifier in file
  !   integer(HID_T) :: memspace      ! Dataspace identifier in mem
  !   integer(HID_T) :: plist_id      ! Property list identifier for parallel IO
  !   integer(HSIZE_T) :: count(6), countm(6), offset(6), offsetm(6), stride(6), block_(6)
  !   integer(HSIZE_T) :: countr(6), offsetr(6), strider(6), blockr(6)
  !   integer :: comm, info, rowremainder, colremainder
  !   PUSH_SUB(write_matrix_f_par_hdf)

  !   if (peinf%inode .eq. 0) call timacc(47,1)

  !   ! JRD: We need a barrier here or else parallel file opening gets mixed up with
  !   ! peinf%inode 0 opening the file to write the diagonal (which is called first).
  !   call MPI_barrier(MPI_COMM_WORLD, mpierr)

  !   comm = MPI_COMM_WORLD
  !   info = MPI_INFO_NULL

  !   if (peinf%verb_debug .and. peinf%inode .eq. 0) then
  !      write(6,*) 'Writing matrix: ', nmtx
  !      write(6,*)
  !   endif

  !   ! JRD Should be ok with npr and npc = 0
  !   rank = 6
  !   countm(1) = 2
  !   countm(2) = MAX(1, scal%npr)
  !   countm(3) = MAX(1, scal%npc)
  !   countm(4) = MAX(1, nfreq_in_group)
  !   countm(5) = 1
  !   countm(6) = 1
  !   offsetm(:) = 0

  !   count(1) = 1
  !   count(2) = scal%npr/scal%nbl
  !   count(3) = scal%npc/scal%nbl
  !   if (nfreq_group .gt. 1) then
  !      count(4) = nfreq_in_group
  !   else
  !      count(4) = 1
  !   endif
  !   count(5) = 1
  !   count(6) = 1

  !   block_(1) = 2
  !   block_(2) = scal%nbl
  !   block_(3) = scal%nbl
  !   if (nfreq_group .gt. 1) then
  !      block_(4) = 1
  !   else
  !      block_(4) = nfreq_in_group
  !   endif
  !   block_(5) = 1
  !   block_(6) = 1

  !   offset(1) = 0
  !   offset(2) = scal%myprow*scal%nbl
  !   offset(3) = scal%mypcol*scal%nbl
  !   offset(4) = peinf%igroup_f
  !   offset(5) = is-1
  !   offset(6) = iq-1

  !   stride(1) = 1
  !   stride(2) = scal%nprow*scal%nbl
  !   stride(3) = scal%npcol*scal%nbl
  !   stride(4) = nfreq_group
  !   stride(5) = 1
  !   stride(6) = 1

  !   call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
  !   call h5pset_fapl_mpio_f(plist_id, comm, info, error)

  !   call open_file(99, trim(name), status='old')
  !   call close_file(99)

  !   call h5fopen_f(trim(name), H5F_ACC_RDWR_F, file_id, error, access_prp = plist_id)
  !   call h5pclose_f(plist_id, error)

  !   SAFE_ALLOCATE(data, (countm(1),countm(2),countm(3),countm(4),countm(5),countm(6)))
  !   !XXX create data can we avoid duplication?
  !   !XXX THREAD? Yes we should to get better bandwidth
  !   if (scal%npr/=0 .and. scal%npc/=0 .and. nfreq_in_group/=0) then
  !      data(1,:,:,:,1,1) = dble(retarded(:,:,:))
  !      data(2,:,:,:,1,1) = IMAG(retarded(:,:,:))
  !   endif
  !   call h5screate_simple_f(rank, countm, memspace, error)

  !   if (scal%npr*scal%npc .ne. 0) then
  !      call h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, offsetm, countm, error)
  !   else
  !      call H5sselect_none_f(memspace,error)
  !   endif

  !   call h5dopen_f(file_id, 'mats/matrix', dset_id, error)
  !   call h5dget_space_f(dset_id, filespace, error)

  !   if (scal%npr*scal%npc .ne. 0) then
  !      call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, count, error, stride, block_)
  !   else
  !      call H5sselect_none_f(filespace, error)
  !   endif

  !   ! Add in remainders, in case scal%nbl doesnt perfectly divide nmtx

  !   ! The condition that nfreq_in_group .ne. 0 is here because this
  !   ! condition will only be met by the processors that are excluded
  !   ! from the calculation when using parallel frequencies with a
  !   ! number of processors not divisible by the number of frequency
  !   ! groups (for the paranoid: this condition that some processors don`t own any frequencies
  !   ! could also be met if the number of frequency groups
  !   ! requested by the user is greater than the total number of frequencies calculated,
  !   ! but we don`t allow that, i.e. the code dies if the user makes such a request).
  !   ! They will have a row/colremainder = 0 because they own no part of the dielectric
  !   ! matrix, but we don`t want them to be involved with the hyperslab operations because
  !   ! they have no data and are spectators.

  !   ! Bottom Rows
  !   rowremainder = mod(scal%npr, scal%nbl)
  !   if (rowremainder .ne. 0 .and. nfreq_in_group .ne. 0) then
  !      offsetr=offset
  !      countr=count
  !      blockr=block_
  !      strider=stride
  !      offsetr(2)=nmtx-rowremainder
  !      countr(2)=rowremainder
  !      blockr(2)=1
  !      strider(2)=1
  !      call h5sselect_hyperslab_f(filespace, H5S_SELECT_OR_F, offsetr, countr, error, strider, blockr)
  !      !write(6,*) peinf%inode, "I have the bottom row", rowremainder, scal%npc
  !   endif

  !   ! Right Columns
  !   colremainder = mod(scal%npc, scal%nbl)
  !   if (colremainder .ne. 0 .and. nfreq_in_group .ne. 0) then
  !      offsetr=offset
  !      countr=count
  !      blockr=block_
  !      strider=stride
  !      offsetr(3)=nmtx-colremainder
  !      countr(3)=colremainder
  !      blockr(3)=1
  !      strider(3)=1
  !      call h5sselect_hyperslab_f(filespace, H5S_SELECT_OR_F, offsetr, countr, error, strider, blockr)
  !      !write(6,*) peinf%inode, "I have the right column", colremainder, scal%npr
  !      ! Bottom Corner of Matrix
  !      if (rowremainder .ne. 0) then
  !         offsetr=offset
  !         countr=count
  !         blockr=block_
  !         strider=stride
  !         offsetr(2)=nmtx-rowremainder
  !         countr(2)=rowremainder
  !         blockr(2)=1
  !         strider(2)=1
  !         offsetr(3)=nmtx-colremainder
  !         countr(3)=colremainder
  !         blockr(3)=1
  !         strider(3)=1
  !         call h5sselect_hyperslab_f(filespace, H5S_SELECT_OR_F, offsetr, countr, error, strider, blockr)
  !         !write(6,*) peinf%inode, "I have bottom both"
  !      endif
  !   endif

  !   call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
  !   !call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_INDEPENDENT_F, error)
  !   call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
  !   if(peinf%inode.eq.0) call timacc(47,2)
  !   if(peinf%inode.eq.0) call timacc(48,1)
  !   call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, data, countm, error, memspace, filespace, xfer_prp = plist_id)
  !   if(peinf%inode.eq.0) call timacc(48,2)
  !   call h5pclose_f(plist_id, error)

  !   SAFE_DEALLOCATE(data)
  !   call h5dclose_f(dset_id, error)
  !   call h5sclose_f(memspace, error)
  !   call h5sclose_f(filespace, error)
  !   call h5fclose_f(file_id, error)

  !   POP_SUB(write_matrix_f_par_hdf)
  !   return
  ! end subroutine write_matrix_f_par_hdf

  !> [WORKING]
  subroutine write_matrix_f_par_hdf_2(scal, nfreq, retarded, nmtx_max, iq, is, name)
    type(scalapack), intent(in) :: scal
    integer, intent(in) :: nfreq
    complex(DPC), intent(in) :: retarded(:,:,:) !< (scal%npr,scal%npc,nfreq)
    integer, intent(in) :: nmtx_max
    integer, intent(in) :: iq
    integer, intent(in) :: is
    character(len=*), intent(in) :: name
    integer :: error, rank
    real(DP), allocatable :: data(:,:,:,:,:,:)
    integer(HID_T) :: file_id       ! File identifier
    integer(HID_T) :: dset_id       ! Dataset identifier
    integer(HID_T) :: filespace     ! Dataspace identifier in file
    integer(HID_T) :: memspace      ! Dataspace identifier in mem
    integer(HID_T) :: plist_id      ! Property list identifier for parallel IO
    integer(HSIZE_T) :: count(6), countm(6), offset(6), offsetm(6), stride(6), block_(6)
    integer(HSIZE_T) :: countr(6), offsetr(6), strider(6), blockr(6)
    integer :: comm, info, rowremainder, colremainder
    PUSH_SUB(write_matrix_f_par_hdf_2)

    if (peinf%inode .eq. 0) call timacc(47,1)

    ! JRD: We need a barrier here or else parallel file opening gets mixed up with
    ! peinf%inode 0 opening the file to write the diagonal (which is called first).
    call MPI_barrier(MPI_COMM_WORLD, mpierr)

    comm = MPI_COMM_WORLD
    info = MPI_INFO_NULL

    if (peinf%verb_debug .and. peinf%inode .eq. 0) then
       write(6,*) 'Writing matrix: ', nmtx_max
       write(6,*)
    endif

    ! JRD Should be ok with npr and npc = 0
    rank = 6
    countm(1) = 2
    countm(2) = scal%npr
    countm(3) = scal%npc
    countm(4) = nfreq
    countm(5) = 1
    countm(6) = 1
    offsetm(:) = 0

    count(1) = 1
    count(2) = scal%npr/scal%nbr
    count(3) = scal%npc/scal%nbc
    count(4) = 1
    count(5) = 1
    count(6) = 1

    block_(1) = 2
    block_(2) = scal%nbr
    block_(3) = scal%nbc
    block_(4) = nfreq
    block_(5) = 1
    block_(6) = 1

    offset(1) = 0
    offset(2) = scal%myprow*scal%nbr
    offset(3) = scal%mypcol*scal%nbc
    offset(4) = 0
    offset(5) = is-1
    offset(6) = iq-1

    stride(1) = 1
    stride(2) = scal%nprow*scal%nbr
    stride(3) = scal%npcol*scal%nbc
    stride(4) = 1
    stride(5) = 1
    stride(6) = 1

    call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif
    call h5pset_fapl_mpio_f(plist_id, comm, info, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    ! call open_file(99, trim(name), status='old')
    ! call close_file(99)

    call h5fopen_f(trim(name), H5F_ACC_RDWR_F, file_id, error, access_prp = plist_id)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif
    call h5pclose_f(plist_id, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    SAFE_ALLOCATE(data, (countm(1),countm(2),countm(3),countm(4),countm(5),countm(6)))

    !XXX create data can we avoid duplication?
    !XXX THREAD? Yes we should to get better bandwidth
    if (scal%npr*scal%npc .ne. 0) then
       data(1,:,:,:,1,1) = dble(retarded(:,:,:))
       data(2,:,:,:,1,1) = IMAG(retarded(:,:,:))
    endif
    call h5screate_simple_f(rank, countm, memspace, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif
    if (scal%npr*scal%npc .ne. 0) then
       call h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, offsetm, countm, error)
       if (error .ne. 0) then
          call die("HDF5 error", only_root_writes=.true.)
       endif
    else
       call H5sselect_none_f(memspace, error)
       if (error .ne. 0) then
          call die("HDF5 error", only_root_writes=.true.)
       endif
    endif

    call h5dopen_f(file_id, 'mats/matrix', dset_id, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif
    call h5dget_space_f(dset_id, filespace, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    if (scal%npr*scal%npc .ne. 0) then
       call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, count, error, stride, block_)
       if (error .ne. 0) then
          call die("HDF5 error", only_root_writes=.true.)
       endif
    else
       call H5sselect_none_f(filespace, error)
       if (error .ne. 0) then
          call die("HDF5 error", only_root_writes=.true.)
       endif
    endif

    ! Add in remainders, in case scal%nbr&nbc doesnt perfectly divide nmtx_max
    ! The condition that nfreq_in_group .ne. 0 is here because this
    ! condition will only be met by the processors that are excluded
    ! from the calculation when using parallel frequencies with a
    ! number of processors not divisible by the number of frequency
    ! groups (for the paranoid: this condition that some processors don`t own any frequencies
    ! could also be met if the number of frequency groups
    ! requested by the user is greater than the total number of frequencies calculated,
    ! but we don`t allow that, i.e. the code dies if the user makes such a request).
    ! They will have a row/colremainder = 0 because they own no part of the dielectric
    ! matrix, but we don`t want them to be involved with the hyperslab operations because
    ! they have no data and are spectators.

    ! Bottom Rows
    rowremainder = mod(scal%npr, scal%nbr)
    if (rowremainder .ne. 0) then
       offsetr = offset
       countr = count
       blockr = block_
       strider = stride
       offsetr(2) = nmtx_max - rowremainder
       countr(2) = rowremainder
       blockr(2) = 1
       strider(2) = 1
       call h5sselect_hyperslab_f(filespace, H5S_SELECT_OR_F, offsetr, countr, error, strider, blockr)
       if (error .ne. 0) then
          call die("HDF5 error", only_root_writes=.true.)
       endif
       !write(6,*) peinf%inode, "I have the bottom row", rowremainder, scal%npc
    endif

    ! Right Columns
    colremainder = mod(scal%npc, scal%nbc)
    if (colremainder .ne. 0) then
       offsetr = offset
       countr = count
       blockr = block_
       strider = stride
       offsetr(3) = nmtx_max - colremainder
       countr(3) = colremainder
       blockr(3) = 1
       strider(3) = 1
       call h5sselect_hyperslab_f(filespace, H5S_SELECT_OR_F, offsetr, countr, error, strider, blockr)
       if (error .ne. 0) then
          call die("HDF5 error", only_root_writes=.true.)
       endif
       !write(6,*) peinf%inode, "I have the right column", colremainder, scal%npr
       ! Bottom Corner of Matrix
       if (rowremainder .ne. 0) then
          offsetr = offset
          countr = count
          blockr = block_
          strider = stride
          offsetr(2) = nmtx_max - rowremainder
          countr(2) = rowremainder
          blockr(2) = 1
          strider(2) = 1
          offsetr(3) = nmtx_max - colremainder
          countr(3) = colremainder
          blockr(3) = 1
          strider(3) = 1
          call h5sselect_hyperslab_f(filespace, H5S_SELECT_OR_F, offsetr, countr, error, strider, blockr)
          if (error .ne. 0) then
             call die("HDF5 error", only_root_writes=.true.)
          endif
          !write(6,*) peinf%inode, "I have bottom both"
       endif
    endif

    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif
    !call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_INDEPENDENT_F, error)
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif
    if(peinf%inode.eq.0) call timacc(47,2)
    if(peinf%inode.eq.0) call timacc(48,1)
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, data, countm, error, memspace, filespace, xfer_prp = plist_id)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif
    if(peinf%inode.eq.0) call timacc(48,2)
    call h5pclose_f(plist_id, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    SAFE_DEALLOCATE(data)
    call h5dclose_f(dset_id, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif
    call h5sclose_f(memspace, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif
    call h5sclose_f(filespace, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif
    call h5fclose_f(file_id, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    POP_SUB(write_matrix_f_par_hdf_2)
    return
  end subroutine write_matrix_f_par_hdf_2

#endif

  !   !==========================================================================================

  !   subroutine write_matrix_f_hdf(scal, nfreq, retarded, nmtx, iq, is, name)
  !     type(scalapack), intent(in) :: scal
  !     integer, intent(in) :: nfreq
  !     complex(DPC), intent(in) :: retarded(:,:,:) !< (scal%npr,scal%npc,nfreq)
  !     integer, intent(in) :: nmtx
  !     integer, intent(in) :: iq
  !     integer, intent(in) :: is
  !     character(len=*), intent(in) :: name
  !     integer :: ii, jj, error, size, rank
  ! #ifdef USESCALAPACK
  !     real(DP), allocatable :: datatmp(:,:,:,:,:,:)
  !     integer :: irow, icol, irowm, icolm
  !     integer :: icurr
  ! #endif
  !     real(DP), allocatable :: data(:,:,:,:,:,:)
  !     type(progress_info) :: prog_info !< a user-friendly progress report
  !     integer(HID_T) :: file_id       ! File identifier
  !     integer(HID_T) :: dset_id       ! Dataset identifier
  !     integer(HID_T) :: filespace     ! Dataspace identifier in file
  !     integer(HID_T) :: memspace      ! Dataspace identifier in mem
  !     integer(HSIZE_T) :: count(6), offset(6), offsetm(6)
  !     PUSH_SUB(write_matrix_f_hdf)

  !     if (peinf%verb_debug .and. peinf%inode==0) then
  !        write(6,*) 'Writing matrix: ', nmtx, nfreq
  !        write(6,*)
  !     endif

  !     ! XXX: For now, we will still have only proc 0 write...
  !     ! We should changes this to parallel writes. But doing
  !     ! this effectively from the scalapack, block cyclic layout
  !     ! seems a bit tricky. So, ignoring for now...

  !     rank=6
  !     count(1) = 2
  !     count(2) = nmtx
  !     count(3) = 1
  !     count(4) = nfreq
  !     count(5) = 1
  !     count(6) = 1

  !     if (peinf%inode .eq. 0) then
  !        SAFE_ALLOCATE(data, (count(1),count(2),count(3),count(4),count(5),count(6)))

  !        call open_file(99, trim(name), status='old')
  !        call close_file(99)

  !        call h5fopen_f(trim(name), H5F_ACC_RDWR_F, file_id, error)
  !        call h5dopen_f(file_id, 'mats/matrix', dset_id, error)
  !        call h5screate_simple_f(rank, count, memspace, error)
  !        call h5dget_space_f(dset_id,filespace,error)
  !     endif

  ! #ifdef USESCALAPACK
  !     SAFE_ALLOCATE(datatmp, (count(1),count(2),count(3),count(4),count(5),count(6)))
  !     icurr=0
  ! #endif

  !     call progress_init(prog_info, 'writing matrix', 'column', nmtx)
  !     do jj = 1, nmtx
  !        call progress_step(prog_info, jj)

  ! #ifdef USESCALAPACK

  !        !    if(peinf%inode.eq.0) call timacc(51,1)

  !        icol=MOD(INT(((jj-1)/scal%nbl)+TOL_SMALL),scal%npcol)
  !        datatmp=0d0
  !        ! JRD XXX The below is going to be incredibly slow. Hoping all over memory.
  !        if (icol .eq. scal%mypcol) then
  !           do ii = 1, nmtx
  !              irow=MOD(INT(((ii-1)/scal%nbl)+TOL_SMALL),scal%nprow)
  !              if (irow .eq. scal%myprow) then
  !                 icurr=icurr+1
  !                 icolm=INT((icurr-1)/scal%npr+TOL_SMALL)+1
  !                 irowm=MOD((icurr-1),scal%npr)+1
  !                 datatmp(1,ii,1,:,1,1)=dble(retarded(irowm,icolm,:))
  !                 datatmp(2,ii,1,:,1,1)=IMAG(retarded(irowm,icolm,:))
  !              endif
  !           enddo
  !        endif
  !        if (peinf%inode .eq. 0) then
  !           data=0d0
  !        endif
  !        ! XXX This is a big waste of communication. Should be fixed when do
  !        ! parallel IO.

  !        size = nmtx*nfreq*2

  !        !    if(peinf%inode.eq.0) call timacc(51,2)
  !        !    if(peinf%inode.eq.0) call timacc(52,1)

  !        call MPI_REDUCE(datatmp,data,size,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
  !             MPI_COMM_WORLD,mpierr)

  !        !    if(peinf%inode.eq.0) call timacc(52,2)

  ! #else
  !        ! JRD XXX The below is horrendously slow. Jumping all over memory
  !        if (peinf%inode .eq. 0) then
  !           do ii = 1, nmtx
  !              data(1,ii,1,:,1,1)=dble(retarded(ii,jj,:))
  !              data(2,ii,1,:,1,1)=IMAG(retarded(ii,jj,:))
  !           enddo
  !        endif

  ! #endif

  !        !    if(peinf%inode.eq.0) call timacc(53,1)

  !        if (peinf%inode .eq. 0) then

  !           offsetm(:) = 0
  !           call h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, offsetm, count, error)

  !           offset(1)=0
  !           offset(2)=0
  !           offset(3)=jj-1
  !           offset(4)=0
  !           offset(5)=is-1
  !           offset(6)=iq-1

  !           call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, count, error)

  !           call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, data, count, error, memspace, filespace)

  !        endif

  !        !    if(peinf%inode.eq.0) call timacc(53,2)
  !        !    if(peinf%inode.eq.0) call timacc(54,1)

  ! #ifdef USESCALAPACK
  !        call MPI_barrier(MPI_COMM_WORLD,mpierr)
  ! #endif

  !        !    if(peinf%inode.eq.0) call timacc(54,2)


  !     enddo
  !     call progress_free(prog_info)

  ! #ifdef USESCALAPACK
  !     SAFE_DEALLOCATE(datatmp)
  ! #endif

  !     if (peinf%inode .eq. 0) then
  !        SAFE_DEALLOCATE(data)
  !     endif

  !     if (peinf%inode .eq. 0) then
  !        call h5dclose_f(dset_id, error)
  !        call h5sclose_f(memspace, error)
  !        call h5sclose_f(filespace, error)
  !        call h5fclose_f(file_id, error)
  !     endif

  !     POP_SUB(write_matrix_f_hdf)

  !     return
  !   end subroutine write_matrix_f_hdf

  ! =======================================================================
  ! ======
  ! From epsilon/epsinv.f90:
  ! call write_gvec_indices_hdf(gvec%ng,pol%isrtx,isorti,ekin,my_iq,filename)
  ! ------
  ! ng = gvec%ng = ngm_g
  ! gind_eps2rho(1:ng) = pol%isrtx ==> sort gvec%components(1:3,1:gvec%ng) = g_global(1:3,1:ng=ngm_g)
  ! gind_rho2eps(1:ng) = isorti ==> inverse sorting index of pol%isrtx
  ! ------
  ! from epsilon_main.f90:
  ! ======
  ! Common/sort.p.f:
  ! subroutine sortrx_gvec(NVAL, AA, ord, gvec)
  ! the rank is written into pol%isrtx in increasing order
  ! ekin(pol%isrtx(1)) < ekin(pol%isrtx(2)) in this way, we can get the value of the i-th smallest in ekin
  ! call sortrx(gvec%ng, ekin, pol%isrtx, gvec = gvec%components)
  ! ------
  ! from epsilon_main.f90:
  ! do i=1,gvec%ng
  !    isorti(pol%isrtx(i)) = i
  ! end do
  ! ------
  ! ekin(1:ng) = ekin
  ! iq = my_iq

  subroutine write_gvec_indices_hdf(ng, gind_eps2rho, gind_rho2eps, ekin, iq, name)
    integer, intent(in) :: ng
    integer, intent(in) :: gind_eps2rho(:) !< (ng)
    integer, intent(in) :: gind_rho2eps(:) !< (ng)
    real(DP), intent(in) :: ekin(:) !< (ng)
    integer, intent(in) :: iq
    character(len=*), intent(in) :: name
    integer(HID_T) :: file_id
    integer :: error, countf(2), offsetf(2)
    PUSH_SUB(write_gvec_indices_hdf)

    ! call open_file(99, trim(name), status='old')
    ! call close_file(99)

    call h5fopen_f(trim(name), H5F_ACC_RDWR_F, file_id, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif
    countf(:) = (/ng, 1/)
    offsetf(:) = (/0, iq-1/)
    call hdf5_write_int_hyperslab(file_id, 'eps_header/gspace/gind_eps2rho', countf, offsetf, gind_eps2rho, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif
    call hdf5_write_int_hyperslab(file_id, 'eps_header/gspace/gind_rho2eps', countf, offsetf, gind_rho2eps, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif
    call hdf5_write_double_hyperslab(file_id, 'eps_header/gspace/ekin', countf, offsetf, ekin, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    ! !> Flushes the internal HDF5 buffers then asks the operating system (the OS) to flush the system buffers for the open files
    ! call h5fflush_f(file_id, H5F_SCOPE_GLOBAL_F, error)
    ! if (error .ne. 0) then
    !    call die("HDF5 error", only_root_writes=.true.)
    ! endif

    call h5fclose_f(file_id, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    POP_SUB(write_gvec_indices_hdf)
  end subroutine write_gvec_indices_hdf

  subroutine write_pol_fftgrid_hdf(pol, name)
    type (polarizability), intent(in) :: pol
    character(len=*), intent(in) :: name
    integer(HID_T) :: file_id
    integer :: error
    PUSH_SUB(write_pol_fftgrid_hdf)

    call h5fopen_f(trim(name), H5F_ACC_RDWR_F, file_id, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif
    call hdf5_write_int_array(file_id, 'eps_header/params/polFFTgrid', (/3/), pol%FFTgrid(:), error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif
    call h5fclose_f(file_id, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    POP_SUB(write_pol_fftgrid_hdf)
  end subroutine write_pol_fftgrid_hdf

#endif

end module write_matrix_m

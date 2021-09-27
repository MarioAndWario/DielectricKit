!=================================================================================
!
! Module:
!
! write_matrix_m          Originally by JRD       Last Modified 5/1/2008 (JRD)
!
! This program writes a distributed matrix like chimat or epsmat to file.
!
!=================================================================================

#include "f_defs.h"

module write_matrix_m
  use global_m
  use hdf5
  use hdf5_io_m
  use scalapack_m
  use io_utils_m
  implicit none
  private
  public :: write_matrix_ser_hdf, write_matrix_f_ser_hdf, &
       write_gvec_indices_hdf, write_matrix_d_par_hdf, write_matrix_f_par_hdf !, write_pol_fftgrid_hdf
contains

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

  !> Use nbr&nbc instead of nbl
  subroutine write_matrix_d_par_hdf(scal, matrix, nmtx_max, iq, is, name)
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

    call h5fclose_f(file_id, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    POP_SUB(write_matrix_d_par_hdf)
    return
  end subroutine write_matrix_d_par_hdf

  !========================================================================================

  subroutine write_matrix_f_par_hdf(scal, nfreq, retarded, nmtx_max, iq, is, name)
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
    PUSH_SUB(write_matrix_f_par_hdf)

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

    call h5fopen_f(trim(name), H5F_ACC_RDWR_F, file_id, error, access_prp = plist_id)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif
    call h5pclose_f(plist_id, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    SAFE_ALLOCATE(data, (countm(1),countm(2),countm(3),countm(4),countm(5),countm(6)))

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

    POP_SUB(write_matrix_f_par_hdf)
    return
  end subroutine write_matrix_f_par_hdf

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

  ! subroutine write_pol_fftgrid_hdf(pol, name)
  !   type (polarizability), intent(in) :: pol
  !   character(len=*), intent(in) :: name
  !   integer(HID_T) :: file_id
  !   integer :: error
  !   PUSH_SUB(write_pol_fftgrid_hdf)

  !   call h5fopen_f(trim(name), H5F_ACC_RDWR_F, file_id, error)
  !   if (error .ne. 0) then
  !      call die("HDF5 error", only_root_writes=.true.)
  !   endif
  !   call hdf5_write_int_array(file_id, 'eps_header/params/polFFTgrid', (/3/), pol%FFTgrid(:), error)
  !   if (error .ne. 0) then
  !      call die("HDF5 error", only_root_writes=.true.)
  !   endif
  !   call h5fclose_f(file_id, error)
  !   if (error .ne. 0) then
  !      call die("HDF5 error", only_root_writes=.true.)
  !   endif

  !   POP_SUB(write_pol_fftgrid_hdf)
  ! end subroutine write_pol_fftgrid_hdf

end module write_matrix_m

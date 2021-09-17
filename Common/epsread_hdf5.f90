#include "f_defs.h"

!>=========================================================================
!!
!!  Module:
!!
!!  epsread_hdf5_m     Originally by JRD     Last Modified 12/2014 (FHJ)
!!
!!    Routines to read header info and matrices from epsmat files in
!!    HDF5 format.
!!
!!=========================================================================

module epsread_hdf5_m
#ifdef HDF5
  use global_m
  use hdf5
  use hdf5_io_m
  implicit none
  private
  public :: &
       read_eps_params_hdf5, read_eps_params_hdf5_2, read_eps_matrix_flavor_hdf5, &
       read_eps_grid_sizes_hdf5, &
       read_eps_matrix_diagonal_hdf5, read_eps_matrix_diagonal_allq_hdf5, &
       read_eps_qgrid_hdf5, &
       read_eps_freqgrid_hdf5, &
       read_eps_old_gvecs_hdf5, &
       read_eps_gvecsofq_hdf5, &
       read_eps_matrix_col_f_hdf5, &
       read_eps_matrix_col_hdf5, read_eps_matrix_col_hdf5_2, &
       read_eps_matrix_ser_hdf5, read_eps_matrix_ser_hdf5_2, read_eps_matrix_ser_allq_hdf5, &
       read_eps_matrix_par_hdf5, read_eps_matrix_par_hdf5_2, read_eps_matrix_par_hdf5_3, read_eps_matrix_par_hdf5_3_, read_eps_matrix_par_hdf5_4, read_eps_matrix_par_allq_hdf5, &
       read_eps_matrix_par_f_hdf5, read_eps_matrix_par_f_hdf5_2, read_eps_matrix_par_f_hdf5_3, read_eps_matrix_par_f_hdf5_3_, read_eps_matrix_par_f_hdf5_4, read_eps_matrix_par_distribute_rq_hdf5
  ! read_vcoul_hdf5, &
contains

  subroutine read_eps_params_hdf5(pol, name, nband)
    type(polarizability), intent(inout) :: pol
    character(len=*), intent(in) :: name
    integer, intent(out), optional :: nband
    integer(HID_T) :: file_id       ! File identifier
    integer :: error
    logical :: exists
    PUSH_SUB(read_eps_params_hdf5)

    call h5fopen_f(trim(name), H5F_ACC_RDONLY_F, file_id, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call hdf5_read_int(file_id, 'eps_header/params/matrix_type', pol%matrix_type, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call hdf5_read_logical(file_id, 'eps_header/params/has_advanced', pol%has_advanced, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call hdf5_read_int(file_id, 'eps_header/params/nmatrix', pol%nmatrix, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call hdf5_read_int(file_id, 'eps_header/params/matrix_flavor', pol%matrix_flavor, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call hdf5_read_int(file_id, 'eps_header/params/icutv', pol%icutv, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call hdf5_read_double(file_id, 'eps_header/params/ecuts', pol%ecuts, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    if (present(nband)) then
       call hdf5_read_int(file_id, 'eps_header/params/nband', nband, error)
       if (error .ne. 0) then
          call die("HDF5 error", only_root_writes=.true.)
       endif
    endif
    call hdf5_read_double(file_id, 'eps_header/params/efermi', pol%efermi, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    pol%efermi = pol%efermi*ryd
    call hdf5_read_int(file_id, 'eps_header/params/intraband_flag', pol%intraband_flag, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call hdf5_read_double(file_id, 'eps_header/params/intraband_overlap_min', pol%intraband_overlap_min, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call hdf5_read_logical(file_id, 'eps_header/params/subsample', pol%subsample, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call h5lexists_f(file_id, 'eps_header/params/timeordered', exists, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif
    if (exists .and. error == 0) then
       call hdf5_read_logical(file_id, 'eps_header/params/timeordered', pol%timeordered, error)
       if (error .ne. 0) then
          call die("HDF5 error", only_root_writes=.true.)
       endif
    else
       pol%timeordered = .true.
    endif

    call h5fclose_f(file_id, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    POP_SUB(read_eps_params_hdf5)
  end subroutine read_eps_params_hdf5

  !> Read eps_header/params into pol
  subroutine read_eps_params_hdf5_2(name, pol)
    character(len=*), intent(in) :: name
    type(polarizability), intent(inout) :: pol
    integer(HID_T) :: file_id       ! File identifier
    integer :: error
    logical :: exists
    PUSH_SUB(read_eps_params_hdf5_2)

    call h5fopen_f(TRIM(name), H5F_ACC_RDONLY_F, file_id, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call hdf5_read_int_array(file_id, 'eps_header/params/FFTgrid', (/3/), pol%FFTgrid, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call hdf5_read_double(file_id, 'eps_header/params/ecuts', pol%ecuts, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call hdf5_read_double(file_id, 'eps_header/params/efermi', pol%efermi, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif
    pol%efermi = pol%efermi*ryd

    call hdf5_read_logical(file_id, 'eps_header/params/has_advanced', pol%has_advanced, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call hdf5_read_int(file_id, 'eps_header/params/icutv', pol%icutv, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call hdf5_read_int(file_id, 'eps_header/params/intraband_flag', pol%intraband_flag, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call hdf5_read_double(file_id, 'eps_header/params/intraband_overlap_min', pol%intraband_overlap_min, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call hdf5_read_int(file_id, 'eps_header/params/matrix_flavor', pol%matrix_flavor, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call hdf5_read_int(file_id, 'eps_header/params/matrix_type', pol%matrix_type, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call hdf5_read_int(file_id, 'eps_header/params/nband', pol%nband, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call hdf5_read_int(file_id, 'eps_header/params/ncb', pol%ncb, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call hdf5_read_int(file_id, 'eps_header/params/nmatrix', pol%nmatrix, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call hdf5_read_int(file_id, 'eps_header/params/nvb', pol%nvb, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call hdf5_read_int(file_id, 'eps_header/params/skip_ncb', pol%skip_ncb, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call hdf5_read_int(file_id, 'eps_header/params/skip_nvb', pol%skip_nvb, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call hdf5_read_logical(file_id, 'eps_header/params/subsample', pol%subsample, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call hdf5_read_logical(file_id, 'eps_header/params/subspace', pol%subspace, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call h5lexists_f(file_id, 'eps_header/params/timeordered', exists, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif
    if (exists .and. error == 0) then
       call hdf5_read_logical(file_id, 'eps_header/params/timeordered', pol%timeordered, error)
       if (error .ne. 0) then
          call die("HDF5 error", only_root_writes=.true.)
       endif
    else
       pol%timeordered = .true.
    endif

    call h5fclose_f(file_id, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    POP_SUB(read_eps_params_hdf5_2)
  end subroutine read_eps_params_hdf5_2

  subroutine read_eps_matrix_flavor_hdf5(matrix_flavor, name)
    integer, intent(out) :: matrix_flavor
    character(len=*), intent(in) :: name
    integer(HID_T) :: file_id       ! File identifier
    integer :: error
    PUSH_SUB(read_eps_matrix_flavor_hdf5)

    call h5fopen_f(trim(name), H5F_ACC_RDONLY_F, file_id, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call hdf5_read_int(file_id, 'eps_header/params/matrix_flavor', matrix_flavor, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call h5fclose_f(file_id, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    POP_SUB(read_eps_matrix_flavor_hdf5)
  end subroutine read_eps_matrix_flavor_hdf5

  !===================================================================================
  ! From Sigma/epscopy.f90: epscopy_init(sig,neps)
  ! ------
  ! call read_eps_grid_sizes_hdf5(ngarbage1, sig%nq0, dgarbage1, sig%nFreq, sig%nfreq_imag, nmtxmax, ngarbage3, freq_dep_flag, 'eps0mat.h5')
  ! ------
  ! ng = ngarbage1
  ! nq = sig%nq0
  ! ecuts = dgarbage1
  ! nfreq = sig%nFreq
  ! nfreq_imag = sig%nfreq_imag
  ! nmtxmax = nmtxmax
  ! qgrid = ngarbage3
  ! freq_dep_flag = freq_dep
  ! name = 'eps0mat.h5'
  ! ------
  ! From sigma/epscopy.f90: epscopy(...)
  ! call read_eps_grid_sizes_hdf5(ng_old, nq_tmp, ecuts, nfreq, nfreq_imag, nmtxmax, qgrid, freq_dep_flag, TRUNC(fname))
  ! ------
  ! [INPUT]
  ! name = TRUNC(fname)
  ! ------
  ! [OUTPUT]
  ! ng = ng_old
  ! nq_tmp = nq
  !> [serial]
  subroutine read_eps_grid_sizes_hdf5(ng, nq, ecuts, nfreq, nfreq_imag, nmtxmax, qgrid, freq_dep, name)
    integer, intent(out) :: ng
    integer, intent(out) :: nq
    real(DP), intent(out) :: ecuts
    integer, intent(out) :: nfreq
    integer, intent(out) :: nfreq_imag
    integer, intent(out) :: nmtxmax
    integer, intent(out) :: qgrid(3)
    integer, intent(out) :: freq_dep
    character(len=*), intent(in) :: name
    integer(HID_T) :: file_id       ! File identifier
    integer :: error
    logical :: exists
    PUSH_SUB(read_eps_grid_sizes_hdf5)

    ! call open_file(99, trim(name), status='old')
    ! call close_file(99)
    call h5fopen_f(trim(name), H5F_ACC_RDONLY_F, file_id, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    ! FHJ: We support version 2 onwards
    call hdf5_require_version(file_id, 'eps_header/versionnumber', 2, trim(name))
    call hdf5_require_version(file_id, 'mf_header/versionnumber', VER_WFN_HDF5, trim(name))
    call hdf5_require_flavor(file_id, 'eps_header/flavor', SCALARSIZE, trim(name))
    call hdf5_require_version(file_id, 'mf_header/flavor', SCALARSIZE, trim(name))

    call hdf5_read_int(file_id, 'eps_header/qpoints/nq', nq, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call hdf5_read_int(file_id, 'mf_header/gspace/ng', ng, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call hdf5_read_int(file_id, 'eps_header/freqs/nfreq', nfreq, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call h5lexists_f(file_id, 'eps_header/freqs/nfreq_imag', exists, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    if (exists.and.error==0) then
       call hdf5_read_int(file_id, 'eps_header/freqs/nfreq_imag', nfreq_imag, error)
       if (error .ne. 0) then
          call die("HDF5 error", only_root_writes=.true.)
       endif
    else
       nfreq_imag = 0
    endif

    call hdf5_read_int(file_id, 'eps_header/freqs/freq_dep', freq_dep, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call hdf5_read_int(file_id, 'eps_header/gspace/nmtx_max', nmtxmax, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call hdf5_read_int_array(file_id, 'eps_header/qpoints/qgrid', (/3/), qgrid, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call hdf5_read_double(file_id, 'eps_header/params/ecuts', ecuts, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call h5fclose_f(file_id,error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    POP_SUB(read_eps_grid_sizes_hdf5)
  end subroutine read_eps_grid_sizes_hdf5

  !====================================================================================

  ! subroutine read_eps_matrix_diagonal_hdf5(nmtx, iq, epsdiag, name)
  !   integer, intent(in) :: nmtx
  !   integer, intent(in) :: iq
  !   real(DP), intent(out) :: epsdiag(:,:) !< (matrix_flavor,nmtx)
  !   character(len=*), intent(in) :: name

  !   integer(HID_T) :: file_id       ! File identifier
  !   integer :: error, matrix_flavor

  !   PUSH_SUB(read_eps_matrix_diagonal_hdf5)

  !   call h5fopen_f(trim(name), H5F_ACC_RDONLY_F, file_id, error)
  !   if (error .ne. 0) then
  !      call die("HDF5 error", only_root_writes=.true.)
  !   endif

  !   call hdf5_read_int(file_id, 'eps_header/params/matrix_flavor', matrix_flavor, error)
  !   if (error .ne. 0) then
  !      call die("HDF5 error", only_root_writes=.true.)
  !   endif

  !   if (matrix_flavor/=size(epsdiag,1)) then
  !      write(0,*) 'ERROR: Got size(epsdiag,1)=',size(epsdiag,1),', but matrix_flavor=',matrix_flavor
  !      call die('Internal error in read_eps_matrix_diagonal_hdf5', &
  !           only_root_writes=.true.)
  !   endif
  !   call hdf5_read_double_hyperslab(file_id, 'mats/matrix-diagonal', &
  !        (/matrix_flavor,nmtx,1/), (/0,0,iq-1/), epsdiag, error)
  !   if (error .ne. 0) then
  !      call die("HDF5 error", only_root_writes=.true.)
  !   endif

  !   call h5fclose_f(file_id,error)
  !   if (error .ne. 0) then
  !      call die("HDF5 error", only_root_writes=.true.)
  !   endif

  !   POP_SUB(read_eps_matrix_diagonal_hdf5)

  ! end subroutine read_eps_matrix_diagonal_hdf5

  subroutine read_eps_matrix_diagonal_hdf5(nmtx, ifreq, iq, epsdiag, name)
    integer, intent(in) :: nmtx
    integer, intent(in) :: ifreq, iq
    real(DP), intent(out) :: epsdiag(:,:) !< (matrix_flavor,nmtx)
    character(len=*), intent(in) :: name
    integer(HID_T) :: file_id       ! File identifier
    integer :: error, matrix_flavor
    integer(HID_T) :: dset_id
    integer(HID_T) :: dataspace
    integer :: rank
    PUSH_SUB(read_eps_matrix_diagonal_hdf5)

    call h5fopen_f(trim(name), H5F_ACC_RDONLY_F, file_id, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call hdf5_read_int(file_id, 'eps_header/params/matrix_flavor', matrix_flavor, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    if (matrix_flavor/=size(epsdiag,1)) then
       write(0,*) 'ERROR: Got size(epsdiag,1)=',size(epsdiag,1),', but matrix_flavor=',matrix_flavor
       call die('Internal error in read_eps_matrix_diagonal_hdf5', only_root_writes=.true.)
    endif

    !> Check the dimension of mats/matrix-diagonal
    call h5dopen_f(file_id, 'mats/matrix-diagonal', dset_id, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call h5dget_space_f(dset_id, dataspace, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call H5Sget_simple_extent_ndims_f(dataspace,rank,error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif
    if (rank .ne. 4) then
       call die("Check version of BGW, matrix-diagonal should have rank=4 to deal with FF.")
    endif

    call h5dclose_f(dset_id, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call hdf5_read_double_hyperslab(file_id, 'mats/matrix-diagonal', (/matrix_flavor,nmtx,1,1/), (/0,0,ifreq-1,iq-1/), epsdiag, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call h5fclose_f(file_id, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    POP_SUB(read_eps_matrix_diagonal_hdf5)
  end subroutine read_eps_matrix_diagonal_hdf5

  !> Serial code
  !> For epsmat for all q, one ifreq
  subroutine read_eps_matrix_diagonal_allq_hdf5(epsdiag, nmtxmax, nq, ifreq, name)
    SCALAR, intent(out) :: epsdiag(:,:) !< (nmtxmax,nq)
    integer, intent(in) :: nmtxmax, nq
    integer, intent(in) :: ifreq
    character(len=*), intent(in) :: name
    integer(HID_T) :: file_id       ! File identifier
    integer(HID_T) :: dset_id
    integer(HID_T) :: dataspace
    integer :: error
    integer :: rank
    !> data(imatrix_flavor,ig,iq)
    real(DP), allocatable :: data(:,:,:)
    PUSH_SUB(read_eps_matrix_diagonal_allq_hdf5)

    call h5fopen_f(trim(name), H5F_ACC_RDONLY_F, file_id, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    !> Check the dimension of mats/matrix-diagonal
    call h5dopen_f(file_id, 'mats/matrix-diagonal', dset_id, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call h5dget_space_f(dset_id, dataspace, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call H5Sget_simple_extent_ndims_f(dataspace, rank, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif
    if (rank .ne. 4) then
       call die("Check version of BGW, matrix-diagonal should have rank=4 to deal with FF.")
    endif

    call h5dclose_f(dset_id, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    SAFE_ALLOCATE(data, (SCALARSIZE, nmtxmax, nq))
    call hdf5_read_double_hyperslab(file_id, 'mats/matrix-diagonal', (/SCALARSIZE, nmtxmax, 1, nq/), (/0, 0, ifreq-1, 0/), data, error)

    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif
    epsdiag(1:nmtxmax,1:nq) = SCALARIFY2(data(1,1:nmtxmax,1:nq), data(2,1:nmtxmax,1:nq))
    SAFE_DEALLOCATE(data)

    call h5fclose_f(file_id,error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    POP_SUB(read_eps_matrix_diagonal_allq_hdf5)
  end subroutine read_eps_matrix_diagonal_allq_hdf5

  !====================================================================================
  ! sigma/epscopy.f90:
  ! call read_eps_qgrid_hdf5(nq_tmp, sig%qpt(:,qoffset+1:), nmtx_of_q, TRUNC(fname))
  ! ------
  ! nq = nq_tmp
  ! qpts = sig%qpt(:,qoffset+1:), qoffset is 0 for q0, nq0 for q1
  ! nmtx = nmtx_of_q
  subroutine read_eps_qgrid_hdf5(nq, qpts, nmtx, name)
    integer, intent(in) :: nq
    real(DP), intent(inout) :: qpts(:,:) !< (3,nq)
    integer, intent(out) :: nmtx(:) !< (nq)
    character(len=*), intent(in) :: name
    integer(HID_T) :: file_id       ! File identifier
    integer :: error
    PUSH_SUB(read_eps_qgrid_hdf5)

    call h5fopen_f(trim(name), H5F_ACC_RDONLY_F, file_id, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call hdf5_read_double_array(file_id, 'eps_header/qpoints/qpts', (/3,nq/), qpts, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call hdf5_read_int_array(file_id, 'eps_header/gspace/nmtx', (/nq/), nmtx, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call h5fclose_f(file_id, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    POP_SUB(read_eps_qgrid_hdf5)
  end subroutine read_eps_qgrid_hdf5

  !===================================================================================

  subroutine read_eps_freqgrid_hdf5(nfreq, dFreqGrid, dFreqBrd, name)
    integer, intent(in) :: nfreq
    real(DP), intent(out) :: dFreqGrid(:) !< (nfreq)
    complex(DPC), intent(out) :: dFreqBrd(:) !< (nfreq)
    character(len=*), intent(in) :: name
    real(DP) :: freqs_tmp(2,nfreq)
    integer :: iw
    integer(HID_T) :: file_id       ! File identifier
    integer :: error
    PUSH_SUB(read_eps_freqgrid_hdf5)

    call h5fopen_f(trim(name), H5F_ACC_RDONLY_F, file_id, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call hdf5_read_double_array(file_id, 'eps_header/freqs/freqs', (/2,nfreq/), freqs_tmp, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    do iw=1,nfreq
       dFreqGrid(iw) = freqs_tmp(1,iw)
       dFreqBrd(iw) = DCMPLX(0,freqs_tmp(2,iw))
    enddo
    call h5fclose_f(file_id,error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    POP_SUB(read_eps_freqgrid_hdf5)
  end subroutine read_eps_freqgrid_hdf5

  !=================================================================================
  !> sigma/epscopy.f90:
  !> call read_eps_old_gvecs_hdf5(ng_old, gvecs_old, TRUNC(fname))
  !> ng = ng_old
  !> gvecs = gvecs_old
  subroutine read_eps_old_gvecs_hdf5(ng, gvecs, name)
    integer, intent(in) :: ng
    integer, intent(out) :: gvecs(:,:) !< (3,ng)
    character(len=*), intent(in) :: name
    integer(HID_T) :: file_id       ! File identifier
    integer :: error
    PUSH_SUB(read_eps_old_gvecs_hdf5)

    call h5fopen_f(trim(name), H5F_ACC_RDONLY_F, file_id, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call hdf5_read_int_array(file_id, 'mf_header/gspace/components', (/3,ng/), gvecs, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call h5fclose_f(file_id,error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    POP_SUB(read_eps_old_gvecs_hdf5)
  end subroutine read_eps_old_gvecs_hdf5

  !====================================================================================
  ! =====
  ! from sigma/epscopy.f90:
  ! call read_eps_gvecsofq_hdf5(ng_old,isrtold,isrtinvdummy,ekold,iq,TRUNC(fname))
  ! ------
  ! ng = ng_old
  ! gind_eps2rho = isrtold
  ! gind_rho2eps = isrtinvdummy
  ! ekin = ekold
  ! iq = iq
  subroutine read_eps_gvecsofq_hdf5(ng, gind_eps2rho, gind_rho2eps, ekin, iq, name)
    integer, intent(in) :: ng !< Number of G-vectors
    integer, intent(out) :: gind_eps2rho(:) !< (ng)
    integer, intent(out) :: gind_rho2eps(:) !< (ng)
    real(DP), intent(out) :: ekin(:) !< (ng)
    integer, intent(in) :: iq
    character(len=*), intent(in) :: name
    integer(HID_T) :: file_id
    integer :: error
    integer :: countf(2), offsetf(2)
    PUSH_SUB(read_eps_gvecsofq_hdf5)

    call h5fopen_f(trim(name), H5F_ACC_RDONLY_F, file_id, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    countf(:) = (/ng, 1/)
    offsetf(:) = (/0, iq-1/)

    call hdf5_read_int_hyperslab(file_id, 'eps_header/gspace/gind_eps2rho', countf, offsetf, gind_eps2rho, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call hdf5_read_int_hyperslab(file_id, 'eps_header/gspace/gind_rho2eps', countf, offsetf, gind_rho2eps, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call hdf5_read_double_hyperslab(file_id, 'eps_header/gspace/ekin', countf, offsetf, ekin, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call h5fclose_f(file_id, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    POP_SUB(read_eps_gvecsofq_hdf5)
  end subroutine read_eps_gvecsofq_hdf5

  !===========================================================================================

  subroutine read_eps_matrix_col_f_hdf5(retarded, nFreq, igp, nmtx, iq, is, name, advanced)
    integer, intent(in) :: iq
    integer, intent(in) :: is
    integer, intent(in) :: nFreq
    integer, intent(in) :: nmtx
    integer, intent(in) :: igp
    complex(DPC), intent(out) :: retarded(:,:) !< (ig1=1,nmtx,nFreq)
    character(len=*), intent(in) :: name
    complex(DPC), optional, intent(out) :: advanced(:,:) !< (nmtx,nFreq)
    integer(HID_T) :: file_id       ! File identifier
    integer(HID_T) :: data_id       ! Property list identifier
    integer(HID_T) :: dataspace        ! Property list identifier
    integer(HID_T) :: memspace        ! Property list identifier
    integer :: error
    integer(HSIZE_T) :: count(6), offset(6)
    real(DP), allocatable :: data(:,:,:,:,:,:)
    integer :: nmatrix_per_spin, nspin, buf_sz, version
    logical :: has_advanced
    PUSH_SUB(read_eps_matrix_col_f_hdf5)

    call h5fopen_f(trim(name), H5F_ACC_RDONLY_F, file_id, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    ! FHJ: the default is never to read the advanced matrix, unless the file
    ! version is <3 (on which case we didn`t store the Coulomb interaction)
    call hdf5_read_int(file_id, 'eps_header/versionnumber', version, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call hdf5_read_int(file_id, 'eps_header/params/nmatrix', nmatrix_per_spin, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call hdf5_read_int(file_id, 'mf_header/kpoints/nspin', nspin, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    nmatrix_per_spin = nmatrix_per_spin / nspin
    has_advanced = .false.
    buf_sz = 1
    if (version < 3) then
       call hdf5_read_logical(file_id, 'eps_header/params/has_advanced', has_advanced, error)
       if (error .ne. 0) then
          call die("HDF5 error", only_root_writes=.true.)
       endif
       if (has_advanced) then
          call die("Advanced chimat not supported", only_root_writes=.true.)
       endif
       if (present(advanced) .and. .not.has_advanced) then
          call die('Inconsistent epsmat file: version<3, but no advanced matrix', only_root_writes=.true.)
       endif
    endif

    if (has_advanced) then
       call die("Advanced chimat not supported", only_root_writes=.true.)
       buf_sz = 2
    endif

    call h5dopen_f(file_id, 'mats/matrix', data_id, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call h5dget_space_f(data_id,dataspace,error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    !> [WORKING]
    !> Should use nmtx_max?

    !> Read epsinv(:,igp,ifreq) for one q
    count(:) = (/2, nmtx, 1, nFreq, buf_sz, 1/)
    call h5screate_simple_f(6, count, memspace, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    offset(:) = (/0, 0, igp-1, 0, nmatrix_per_spin*(is-1), iq-1/)
    call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, offset, count, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    SAFE_ALLOCATE(data,(count(1),count(2),count(3),count(4),count(5),count(6)))
    call h5dread_f(data_id, H5T_NATIVE_DOUBLE, data, count, error, memspace, dataspace)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    retarded(1:nmtx,1:nFreq) = DCMPLX(data(1,1:nmtx,1,1:nFreq,1,1),data(2,1:nmtx,1,1:nFreq,1,1))
    if (has_advanced .and. present(advanced)) then
       call die("Advanced chimat not supported", only_root_writes=.true.)
       advanced(1:nmtx,1:nFreq) = DCMPLX(data(1,1:nmtx,1,1:nFreq,2,1),data(2,1:nmtx,1,1:nFreq,2,1))
    endif

    SAFE_DEALLOCATE(data)
    call h5sclose_f(memspace, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    if (.not.has_advanced .and. present(advanced)) then
       call die("epsread: Advanced chimat not supported", only_root_writes=.true.)
       ! call get_advanced_from_retarded()
    endif

    call h5sclose_f(dataspace, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call h5dclose_f(data_id,error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call h5fclose_f(file_id,error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    POP_SUB(read_eps_matrix_col_f_hdf5)

  end subroutine read_eps_matrix_col_f_hdf5

  !===========================================================================================

  subroutine read_eps_matrix_col_hdf5(eps,j,nmtx,iq,is,name)
    integer, intent(in) :: iq
    integer, intent(in) :: is
    integer, intent(in) :: nmtx
    integer, intent(in) :: j
    SCALAR, intent(out) :: eps(:) !< (nmtx)
    character(len=*), intent(in) :: name
    integer(HID_T) :: file_id       ! File identifier
    integer(HID_T) :: data_id       ! Property list identifier
    integer(HID_T) :: dataspace        ! Property list identifier
    integer(HID_T) :: memspace        ! Property list identifier
    integer :: error, rank
    integer(HSIZE_T) :: count(6), offset(6)
    real(DP), allocatable :: data(:,:,:,:,:,:)
    PUSH_SUB(read_eps_matrix_col_hdf5)

    call h5fopen_f(trim(adjustl(name)), H5F_ACC_RDONLY_F, file_id, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    rank = 6
    count(1) = SCALARSIZE
    count(2) = nmtx
    count(3) = 1
    count(4) = 1
    count(5) = 1 !mat
    count(6) = 1 !iq

    offset(:) = 0
    offset(3) = j - 1
    offset(5) = is - 1
    offset(6) = iq - 1

    call h5screate_simple_f(rank, count, memspace, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call h5dopen_f(file_id, '/mats/matrix', data_id, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call h5dget_space_f(data_id,dataspace,error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, offset, count, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    SAFE_ALLOCATE(data,(count(1),count(2),count(3),count(4),count(5),count(6)))
    call h5dread_f(data_id, H5T_NATIVE_DOUBLE, data, count, error, memspace, dataspace)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif
    eps(1:nmtx) = SCALARIFY2(data(1,1:nmtx,1,1,1,1),data(2,1:nmtx,1,1,1,1))
    SAFE_DEALLOCATE(data)
    call h5sclose_f(memspace, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call h5sclose_f(dataspace, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call h5dclose_f(data_id,error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call h5fclose_f(file_id,error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    POP_SUB(read_eps_matrix_col_hdf5)

  end subroutine read_eps_matrix_col_hdf5

  !> used calcualte 0.5 * (chi1 + chi2)
  subroutine read_eps_matrix_col_hdf5_2(eps,j,nmtx,iq,is,name)
    integer, intent(in) :: iq
    integer, intent(in) :: is
    integer, intent(in) :: nmtx
    integer, intent(in) :: j
    SCALAR, intent(out) :: eps(:) !< (nmtx)
    character(len=*), intent(in) :: name
    integer(HID_T) :: file_id       ! File identifier
    integer(HID_T) :: data_id       ! Property list identifier
    integer(HID_T) :: dataspace        ! Property list identifier
    integer(HID_T) :: memspace        ! Property list identifier
    integer :: error, rank
    integer(HSIZE_T) :: count(6), offset(6)
    real(DP), allocatable :: data(:,:,:,:,:,:)
    PUSH_SUB(read_eps_matrix_col_hdf5_2)

    call h5fopen_f(trim(name), H5F_ACC_RDONLY_F, file_id, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call h5dopen_f(file_id, 'mats/matrix', data_id, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call h5dget_space_f(data_id,dataspace,error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    rank = 6
    count(1) = SCALARSIZE
    count(2) = nmtx
    count(3) = 1
    count(4) = 1
    count(5) = 1 !mat
    count(6) = 1 !iq

    call h5screate_simple_f(rank, count, memspace, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    offset(:) = 0
    offset(3) = j - 1
    offset(5) = is - 1
    offset(6) = iq - 1

    call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, offset, count, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    SAFE_ALLOCATE(data,(count(1),count(2),count(3),count(4),count(5),count(6)))
    call h5dread_f(data_id, H5T_NATIVE_DOUBLE, data, count, error, memspace, dataspace)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    eps(1:nmtx) = 0.5 * ( eps(1:nmtx) + SCALARIFY2(data(1,1:nmtx,1,1,1,1),data(2,1:nmtx,1,1,1,1)) )

    SAFE_DEALLOCATE(data)

    call h5sclose_f(memspace, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call h5sclose_f(dataspace, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call h5dclose_f(data_id,error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call h5fclose_f(file_id,error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    POP_SUB(read_eps_matrix_col_hdf5_2)
  end subroutine read_eps_matrix_col_hdf5_2

  !===========================================================================================

  !> Serial code
  subroutine read_eps_matrix_ser_hdf5(eps,nmtx,iq,is,name,ifreq_zero)
    integer, intent(in) :: iq
    integer, intent(in) :: is
    integer, intent(in) :: nmtx
    SCALAR, intent(out) :: eps(:,:) !< (nmtx,nmtx)
    character(len=*), intent(in) :: name
    integer, optional, intent(in) :: ifreq_zero
    integer(HID_T) :: file_id       ! File identifier
    integer(HID_T) :: data_id       ! Property list identifier
    integer(HID_T) :: dataspace        ! Property list identifier
    integer(HID_T) :: memspace        ! Property list identifier
    integer :: error, rank
    integer(HSIZE_T) :: count(6), offset(6)
    real(DP), allocatable :: data(:,:,:,:,:,:)
    PUSH_SUB(read_eps_matrix_ser_hdf5)

    call h5fopen_f(trim(name), H5F_ACC_RDONLY_F, file_id, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call h5dopen_f(file_id, 'mats/matrix', data_id, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call h5dget_space_f(data_id,dataspace,error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    rank = 6
    count(:) = (/ SCALARSIZE, nmtx, nmtx, 1, 1, 1 /)
    call h5screate_simple_f(rank, count, memspace, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    offset(:) = 0
    if (present(ifreq_zero)) then
       offset(4) = ifreq_zero-1
    endif
    offset(5) = is - 1
    offset(6) = iq - 1

    call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, offset, count, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    SAFE_ALLOCATE(data,(count(1),count(2),count(3),count(4),count(5),count(6)))
    call h5dread_f(data_id, H5T_NATIVE_DOUBLE, data, count, error, memspace, dataspace)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    eps(1:nmtx,1:nmtx) = SCALARIFY2(data(1,1:nmtx,1:nmtx,1,1,1),data(2,1:nmtx,1:nmtx,1,1,1))
    SAFE_DEALLOCATE(data)

    call h5sclose_f(memspace, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call h5sclose_f(dataspace, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call h5dclose_f(data_id,error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call h5fclose_f(file_id,error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    POP_SUB(read_eps_matrix_ser_hdf5)
  end subroutine read_eps_matrix_ser_hdf5

  subroutine read_eps_matrix_ser_hdf5_2(eps, nmtxmax, iq, is, ifreq, name)
    SCALAR, intent(out) :: eps(:,:) !< (nmtxmax,nmtxmax)
    integer, intent(in) :: nmtxmax, iq, is, ifreq
    character(len=*), intent(in) :: name
    integer(HID_T) :: file_id       ! File identifier
    integer(HID_T) :: data_id       ! Property list identifier
    integer(HID_T) :: dataspace        ! Property list identifier
    integer(HID_T) :: memspace        ! Property list identifier
    integer :: error, rank
    integer(HSIZE_T) :: count(6), offset(6)
    real(DP), allocatable :: data(:,:,:,:,:,:)
    PUSH_SUB(read_eps_matrix_ser_hdf5_2)

    call h5fopen_f(trim(name), H5F_ACC_RDONLY_F, file_id, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call h5dopen_f(file_id, 'mats/matrix', data_id, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call h5dget_space_f(data_id,dataspace,error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    rank = 6
    count(:) = (/ SCALARSIZE, nmtxmax, nmtxmax, 1, 1, 1 /)
    call h5screate_simple_f(rank, count, memspace, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    offset(:) = (/ 0, 0, 0, ifreq-1, is-1, iq-1 /)
    call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, offset, count, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    SAFE_ALLOCATE(data,(count(1),count(2),count(3),count(4),count(5),count(6)))
    call h5dread_f(data_id, H5T_NATIVE_DOUBLE, data, count, error, memspace, dataspace)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    eps(:,:) = SCALARIFY2(data(1,:,:,1,1,1),data(2,:,:,1,1,1))
    SAFE_DEALLOCATE(data)

    call h5sclose_f(memspace, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call h5sclose_f(dataspace, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call h5dclose_f(data_id,error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call h5fclose_f(file_id,error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    POP_SUB(read_eps_matrix_ser_hdf5_2)
  end subroutine read_eps_matrix_ser_hdf5_2

  !> Serial code
  !> Read epsmat for all q, one freq
  subroutine read_eps_matrix_ser_allq_hdf5(eps, nmtxmax, nq, is, ifreq, name)
    SCALAR, intent(out) :: eps(:,:,:) !< (nmtxmax,nmtxmax,nq)
    integer, intent(in) :: nmtxmax, nq, is
    integer, intent(in) :: ifreq
    character(len=*), intent(in) :: name
    integer(HID_T) :: file_id       ! File identifier
    integer(HID_T) :: data_id       ! Property list identifier
    integer(HID_T) :: dataspace        ! Property list identifier
    integer(HID_T) :: memspace        ! Property list identifier
    integer :: error, rank
    integer(HSIZE_T) :: count(6), offset(6)
    !> data(imatrix_flavor, ig1, ig2, ifreq, is, iq)
    real(DP), allocatable :: data(:,:,:,:,:,:)
    PUSH_SUB(read_eps_matrix_ser_allq_hdf5)

    call h5fopen_f(trim(name), H5F_ACC_RDONLY_F, file_id, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call h5dopen_f(file_id, 'mats/matrix', data_id, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call h5dget_space_f(data_id,dataspace,error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    rank = 6
    count(:) = (/ SCALARSIZE, nmtxmax, nmtxmax, 1, 1, nq /)

    call h5screate_simple_f(rank, count, memspace, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    offset(:) = (/ 0, 0, 0, ifreq-1, is-1, 0 /)
    call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, offset, count, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    SAFE_ALLOCATE(data,(count(1),count(2),count(3),count(4),count(5),count(6)))
    call h5dread_f(data_id, H5T_NATIVE_DOUBLE, data, count, error, memspace, dataspace)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif
    eps(:,:,:) = SCALARIFY2(data(1,:,:,1,1,:),data(2,:,:,1,1,:))
    SAFE_DEALLOCATE(data)

    call h5sclose_f(memspace, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call h5sclose_f(dataspace, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call h5dclose_f(data_id,error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call h5fclose_f(file_id,error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    POP_SUB(read_eps_matrix_ser_allq_hdf5)
  end subroutine read_eps_matrix_ser_allq_hdf5

  !===========================================================================================
  ! from Sigma/epscopy.f90:
  ! call read_eps_matrix_par_hdf5(epsmpi%eps(:,:,iq+qoffset), epsmpi%nb, peinf%pool_rank, peinf%npes_pool, nmtx, iq, 1, fname)
  ! ------
  ! eps(1:npes,1:ngpown) = epsmpi%eps(:,:,iq+qoffset)
  ! nb = epsmpi%nb = 1
  ! rank = peinf%pool_rank = The proc index in a pool for current proc, from 0 to peinf%npes_pool-1
  ! npes = peinf%npes_pool = peinf%npes_pool = LowerBound[peinf%npes/peinf%npools]
  ! nmtx = nmtx = nmtx_of_q(iq) : number of Gvectors for epsinv(:,:,iq)
  ! iq = iq
  ! is = 1
  ! ------
  ! Recall
  subroutine read_eps_matrix_par_hdf5(eps, nb, rank, npes, nmtx, iq, is, name, ifreq_zero)
    SCALAR, intent(inout) :: eps(:,:) !< (neps, ngpown)
    integer, intent(in) :: nb !< block size
    integer, intent(in) :: rank !< processor rank for column distribution
    integer, intent(in) :: npes !< number of processors over which we distribute
    integer, intent(in) :: nmtx
    integer, intent(in) :: iq
    integer, intent(in) :: is
    character(len=*), intent(in) :: name
    integer, optional, intent(in) :: ifreq_zero !> use to specify which frequency is the zero frequency

    integer(HID_T) :: file_id       ! File identifier
    integer(HID_T) :: data_id       ! Property list identifier
    integer(HID_T) :: plist_id       ! Property list identifier
    integer(HID_T) :: dataspace        ! Property list identifier
    integer(HID_T) :: memspace        ! Property list identifier
    integer :: error
    integer :: comm, info
    integer(HSIZE_T) :: count(6), offset(6), countm(6)
    real(DP), allocatable :: data(:,:,:,:,:,:)
    integer :: ngpown_max, igp, igp_loc
    PUSH_SUB(read_eps_matrix_par_hdf5)

    !> Common/scalapack.f90:
    !> Each pool has the same epsmpi
    !> epsmpi(:,:,iq) is distributed within a pool (npes procs)
    !> Here nb = 1, nmtx = nmtx_of_q(iq), npes = peinf%npes_pool
    !> (ngpown_max) is a proc-local quantities, which is the number of Gvectors of epsinv for current proc.
    !> nmtx is Maximal number of Gvectors for the current q vector
    !> so here ngpown_max is q-dependent
    ngpown_max = NUMROC(nmtx, nb, 0, 0, npes)

    !> We only distribute ig2 of epsinv(ig1,ig2,iq) over procs within a pool
    SAFE_ALLOCATE(data,(SCALARSIZE,nmtx,1,1,1,1))

#ifdef MPI
    comm = MPI_COMM_WORLD
    info = MPI_INFO_NULL

    call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call h5pset_fapl_mpio_f(plist_id, comm, info, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call h5fopen_f(trim(name), H5F_ACC_RDONLY_F, file_id, error, access_prp = plist_id)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call h5pclose_f(plist_id,error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

#else
    call h5fopen_f(trim(name), H5F_ACC_RDONLY_F, file_id, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

#endif

    !> Note that ngpown_max is q-dependent
    do igp_loc = 1, ngpown_max
       ! ======
       ! igp_loc is the proc-local index of a Gvector of epsinv
       ! igp is its global index
       ! ------
       ! rank = peinf%pool_rank
       ! npes = peinf%npes_pool
       ! ------
       ! INDXLOC starts from 1
       ! IPROC starts from 0
       ! ------
       ! INDXL2G( INDXLOC, NB, IPROC, ISRCPROC, NPROCS )
       ! INDXL2G = NPROCS*NB*((INDXLOC-1)/NB) + MOD(INDXLOC-1,NB) + MOD(NPROCS+IPROC-ISRCPROC, NPROCS)*NB + 1
       ! ------
       ! suppose in a pool, npes_pool = 4, there are 9 points
       ! proc # 0 : 1 5 9
       ! proc # 1 : 2 6
       ! proc # 2 : 3 7
       ! proc # 3 : 4 8
       ! -------
       igp = INDXL2G(igp_loc, nb, rank, 0, npes)

       call h5dopen_f(file_id, 'mats/matrix', data_id, error)
       if (error .ne. 0) then
          call die("HDF5 error", only_root_writes=.true.)
       endif

       call h5dget_space_f(data_id,dataspace,error)
       if (error .ne. 0) then
          call die("HDF5 error", only_root_writes=.true.)
       endif

       ! JRD: The commented code in this routine represents efforts to use a single HDF5 read call
       ! with an appropriate block and stride for each proc. In general, it currently appears to
       ! perform worse (though this may be related to size of matrix. So, it is commented until
       ! further investigation.

       !> Recall Common/write_matrix.f90: write_matrix_d_par_hdf(scal,matrix,nmtx,iq,is,name)
       countm(1) = SCALARSIZE ! icomplex
       countm(2) = nmtx       ! ig1, q-dependent size
       countm(3) = 1          ! ig2
       countm(4) = 1          ! ifreq
       countm(5) = 1          ! is
       countm(6) = 1          ! iq

       call h5screate_simple_f(6, countm, memspace, error)
       if (error .ne. 0) then
          call die("HDF5 error", only_root_writes=.true.)
       endif

       count(1) = SCALARSIZE
       count(2) = nmtx
       count(3) = 1
       count(4) = 1
       count(5) = 1
       count(6) = 1

       !> Construct data and offset
       if (igp <= nmtx) then

          offset(1)=0     ! icomplex
          offset(2)=0     ! ig1
          offset(3)=igp-1 ! ig2
          if (present(ifreq_zero)) then
             offset(4)=ifreq_zero-1
          else
             offset(4)=0     ! ifreq
          endif
          offset(5)=is-1  ! is
          offset(6)=iq-1  ! iq

          !> Select hyperslab
          call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, offset, count, error)
          if (error .ne. 0) then
             call die("HDF5 error", only_root_writes=.true.)
          endif
       else
          call H5sselect_none_f(memspace,error);
          if (error .ne. 0) then
             call die("HDF5 error", only_root_writes=.true.)
          endif

          call H5sselect_none_f(dataspace,error);
          if (error .ne. 0) then
             call die("HDF5 error", only_root_writes=.true.)
          endif
       endif

#ifdef MPI
       call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
       if (error .ne. 0) then
          call die("HDF5 error", only_root_writes=.true.)
       endif

       call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
       if (error .ne. 0) then
          call die("HDF5 error", only_root_writes=.true.)
       endif

       !> Collectively read the file
       call h5dread_f(data_id, H5T_NATIVE_DOUBLE, data, countm, error, memspace, dataspace, xfer_prp = plist_id)
       if (error .ne. 0) then
          call die("HDF5 error", only_root_writes=.true.)
       endif
# else
       call h5dread_f(data_id, H5T_NATIVE_DOUBLE, data, countm, error, memspace, dataspace)
       if (error .ne. 0) then
          call die("HDF5 error", only_root_writes=.true.)
       endif
#endif
       if (igp <= nmtx) then
          !    if (ngpown .ne. 0 .and. my_igp .le. nmtx) then
          !XXX PROBABLY NEED THREADED LOOP HERE
          eps(1:nmtx,igp_loc) = SCALARIFY2(data(1,1:nmtx,1,1,1,1),data(2,1:nmtx,1,1,1,1))
       endif

       call h5sclose_f(memspace, error)
       if (error .ne. 0) then
          call die("HDF5 error", only_root_writes=.true.)
       endif

       call h5sclose_f(dataspace, error)
       if (error .ne. 0) then
          call die("HDF5 error", only_root_writes=.true.)
       endif

       call h5dclose_f(data_id, error)
       if (error .ne. 0) then
          call die("HDF5 error", only_root_writes=.true.)
       endif
    enddo ! enddo igp_loc

    call h5fclose_f(file_id,error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    SAFE_DEALLOCATE(data)

    POP_SUB(read_eps_matrix_par_hdf5)
  end subroutine read_eps_matrix_par_hdf5

  subroutine read_eps_matrix_par_hdf5_2(eps, nb, pool_rank, npes_pool, my_pool, nmtx, nmtxmax, iq, is, ifreq, name)
    SCALAR, intent(inout) :: eps(:,:) !> (nmtxmax, MAX(ngpown,1))
    integer, intent(in) :: nb !> block size, not necessarily 1
    integer, intent(in) :: pool_rank, npes_pool, my_pool !> not necessarily peinf%pool_rank, peinf%npes_pool
    integer, intent(in) :: nmtx, nmtxmax
    integer, intent(in) :: iq, is
    integer, intent(in) :: ifreq !> use to specify which frequency is the zero frequency
    character(len=*), intent(in) :: name

    integer(HID_T) :: file_id       ! File identifier
    integer(HID_T) :: data_id       ! Property list identifier
    integer(HID_T) :: plist_id       ! Property list identifier
    integer(HID_T) :: dataspace        ! Property list identifier
    integer(HID_T) :: memspace        ! Property list identifier
    integer :: error
    integer(HSIZE_T) :: count(6), offset(6) !, countm(6)
    real(DP), allocatable :: data(:,:,:,:,:,:)
    integer :: ngpown, ngpown_max, igp, igp_loc
    PUSH_SUB(read_eps_matrix_par_hdf5_2)

    !> We only distribute ig2 of epsinv(ig1,ig2,iq) over procs within a pool
    SAFE_ALLOCATE(data, (SCALARSIZE, nmtxmax, 1, 1, 1, 1))

#ifdef MPI
    ngpown     = NUMROC(nmtxmax, nb, pool_rank, 0, npes_pool)
    ngpown_max = NUMROC(nmtxmax, nb,         0, 0, npes_pool)

    call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call h5fopen_f(TRUNC(name), H5F_ACC_RDONLY_F, file_id, error, access_prp = plist_id)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call h5pclose_f(plist_id,error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif
#else
    ngpown = nmtxmax
    ngpown_max = nmtxmax

    call h5fopen_f(TRUNC(name), H5F_ACC_RDONLY_F, file_id, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif
#endif

    do igp_loc = 1, ngpown_max
       if (my_pool>=0) then
          igp = INDXL2G(igp_loc, nb, pool_rank, 0, npes_pool)
          !> <- idle proc ->
       else
          igp = nmtx + 1
       endif

       call h5dopen_f(file_id, 'mats/matrix', data_id, error)
       if (error .ne. 0) then
          call die("HDF5 error", only_root_writes=.true.)
       endif

       call h5dget_space_f(data_id,dataspace,error)
       if (error .ne. 0) then
          call die("HDF5 error", only_root_writes=.true.)
       endif

       count(:) = (/ SCALARSIZE, nmtxmax, 1, 1, 1, 1 /)

       call h5screate_simple_f(6, count, memspace, error)
       if (error .ne. 0) then
          call die("HDF5 error", only_root_writes=.true.)
       endif

       !> Construct data and offset
       ! if (igp <= nmtx) then
       if ((igp_loc <= ngpown) .and. (igp <= nmtx)) then
          offset(:) = (/ 0, 0, igp-1, ifreq-1, is-1, iq-1 /)
          !> Select hyperslab
          call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, offset, count, error)
          if (error .ne. 0) then
             call die("HDF5 error", only_root_writes=.true.)
          endif
       else
          call H5sselect_none_f(memspace,error);
          if (error .ne. 0) then
             call die("HDF5 error", only_root_writes=.true.)
          endif
          call H5sselect_none_f(dataspace,error);
          if (error .ne. 0) then
             call die("HDF5 error", only_root_writes=.true.)
          endif
       endif

#ifdef MPI
       call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
       if (error .ne. 0) then
          call die("HDF5 error", only_root_writes=.true.)
       endif

       call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
       if (error .ne. 0) then
          call die("HDF5 error", only_root_writes=.true.)
       endif

       !> Collectively read the file
       call h5dread_f(data_id, H5T_NATIVE_DOUBLE, data, count, error, memspace, dataspace, xfer_prp = plist_id)
       if (error .ne. 0) then
          call die("HDF5 error", only_root_writes=.true.)
       endif
# else
       call h5dread_f(data_id, H5T_NATIVE_DOUBLE, data, count, error, memspace, dataspace)
       if (error .ne. 0) then
          call die("HDF5 error", only_root_writes=.true.)
       endif
#endif
       if ((igp_loc <= ngpown) .and. (igp <= nmtx)) then
          eps(1:nmtx, igp_loc) = SCALARIFY2(data(1,1:nmtx,1,1,1,1),data(2,1:nmtx,1,1,1,1))
       endif

       call h5sclose_f(memspace, error)
       if (error .ne. 0) then
          call die("HDF5 error", only_root_writes=.true.)
       endif

       call h5sclose_f(dataspace, error)
       if (error .ne. 0) then
          call die("HDF5 error", only_root_writes=.true.)
       endif

       call h5dclose_f(data_id, error)
       if (error .ne. 0) then
          call die("HDF5 error", only_root_writes=.true.)
       endif
    enddo !> igp_loc

    call h5fclose_f(file_id,error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif
    SAFE_DEALLOCATE(data)

    POP_SUB(read_eps_matrix_par_hdf5_2)
  end subroutine read_eps_matrix_par_hdf5_2

  subroutine read_eps_matrix_par_hdf5_3(eps, nb, ngpown, ngpown_max, pool_rank, npes_pool, my_pool, nmtx, nmtxmax, iq, is, ifreq, name)
    SCALAR, intent(inout) :: eps(:,:) !> (nmtxmax, MAX(ngpown,1))
    integer, intent(in) :: nb, ngpown, ngpown_max !> block size, not necessarily 1
    integer, intent(in) :: pool_rank, npes_pool, my_pool !> not necessarily peinf%pool_rank, peinf%npes_pool
    integer, intent(in) :: nmtx, nmtxmax, iq, is
    integer, intent(in) :: ifreq !> use to specify which frequency is the zero frequency
    character(len=*), intent(in) :: name
    integer(HID_T) :: file_id       ! File identifier
    integer(HID_T) :: data_id       ! Property list identifier
    integer(HID_T) :: plist_id       ! Property list identifier
    integer(HID_T) :: dataspace        ! Property list identifier
    integer(HID_T) :: memspace        ! Property list identifier
    integer :: error
    integer(HSIZE_T) :: count(6), offset(6) !, countm(6)
    real(DP), allocatable :: data(:,:,:,:,:,:)
    integer :: igp, igp_loc
    PUSH_SUB(read_eps_matrix_par_hdf5_3)

    !> We only distribute ig2 of epsinv(ig1,ig2,iq) over procs within a pool
    SAFE_ALLOCATE(data, (SCALARSIZE, nmtxmax, 1, 1, 1, 1))
    data = 0.0D0

#ifdef MPI
    call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call h5fopen_f(TRUNC(name), H5F_ACC_RDONLY_F, file_id, error, access_prp = plist_id)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call h5pclose_f(plist_id,error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif
#else
    call h5fopen_f(TRUNC(name), H5F_ACC_RDONLY_F, file_id, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif
#endif

    call h5dopen_f(file_id, 'mats/matrix', data_id, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call h5dget_space_f(data_id,dataspace,error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    do igp_loc = 1, ngpown_max
       if (my_pool>=0) then
          igp = INDXL2G(igp_loc, nb, pool_rank, 0, npes_pool)
          !> <- idle proc ->
       else
          igp = nmtx + 1
       endif

       !> Read one column at a time!
       count(:) = (/ SCALARSIZE, nmtxmax, 1, 1, 1, 1 /)

       call h5screate_simple_f(6, count, memspace, error)
       if (error .ne. 0) then
          call die("HDF5 error", only_root_writes=.true.)
       endif

       !> Construct data and offset
       ! if (igp <= nmtx) then
       if ((igp_loc <= ngpown) .and. (igp <= nmtx)) then
          offset(:) = (/ 0, 0, igp-1, ifreq-1, is-1, iq-1 /)
          !> Select hyperslab
          call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, offset, count, error)
          if (error .ne. 0) then
             call die("HDF5 error", only_root_writes=.true.)
          endif
       else
          call H5sselect_none_f(memspace,error);
          if (error .ne. 0) then
             call die("HDF5 error", only_root_writes=.true.)
          endif
          call H5sselect_none_f(dataspace,error);
          if (error .ne. 0) then
             call die("HDF5 error", only_root_writes=.true.)
          endif
       endif

#ifdef MPI
       call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
       if (error .ne. 0) then
          call die("HDF5 error", only_root_writes=.true.)
       endif

       call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
       if (error .ne. 0) then
          call die("HDF5 error", only_root_writes=.true.)
       endif

       !> Collectively read the file
       call h5dread_f(data_id, H5T_NATIVE_DOUBLE, data, count, error, memspace, dataspace, xfer_prp = plist_id)
       if (error .ne. 0) then
          call die("HDF5 error", only_root_writes=.true.)
       endif
#else
       call h5dread_f(data_id, H5T_NATIVE_DOUBLE, data, count, error, memspace, dataspace)
       if (error .ne. 0) then
          call die("HDF5 error", only_root_writes=.true.)
       endif
#endif
       if ((igp_loc <= ngpown) .and. (igp <= nmtx)) then
          eps(1:nmtx, igp_loc) = SCALARIFY2(data(1,1:nmtx,1,1,1,1), data(2,1:nmtx,1,1,1,1))
       endif

       call h5sclose_f(memspace, error)
       if (error .ne. 0) then
          call die("HDF5 error", only_root_writes=.true.)
       endif
    enddo !> igp_loc

    call h5sclose_f(dataspace, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call h5dclose_f(data_id, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call h5fclose_f(file_id,error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif
    SAFE_DEALLOCATE(data)

    POP_SUB(read_eps_matrix_par_hdf5_3)
  end subroutine read_eps_matrix_par_hdf5_3

  subroutine read_eps_matrix_par_hdf5_3_(eps, nb, ngpown_, ngpown_max, pool_rank, npes_pool, my_pool, nmtx, nmtxmax, iq, is, ifreq_zero, name)
    SCALAR, intent(inout) :: eps(:,:) !> (nmtxmax, MAX(ngpown,1))
    integer, intent(in) :: nb, ngpown_, ngpown_max !> block size, not necessarily 1
    integer, intent(in) :: pool_rank, npes_pool, my_pool !> not necessarily peinf%pool_rank, peinf%npes_pool
    integer, intent(in) :: nmtx, nmtxmax, iq, is
    integer, intent(in) :: ifreq_zero !> use to specify which frequency is the zero frequency
    character(len=*), intent(in) :: name
    integer(HID_T) :: file_id       ! File identifier
    integer(HID_T) :: data_id       ! Property list identifier
    integer(HID_T) :: plist_id      ! Property list identifier
    integer(HID_T) :: dataspace     ! Property list identifier
    integer(HID_T) :: memspace      ! Property list identifier
    integer :: error
    integer(HSIZE_T) :: count(6), offset(6) !, countm(6)
    real(DP), allocatable :: data(:,:,:,:,:,:)
    integer :: igp_offset, igp_loc, igp_start, ngpown
    PUSH_SUB(read_eps_matrix_par_hdf5_3_)

    !> We only distribute ig2 of epsinv(ig1,ig2,iq) over procs within a pool
    igp_start = INDXL2G(1, nb, pool_rank, 0, npes_pool)
    if (nmtxmax - igp_start < 0) then
       ngpown = 0
    else
       ngpown = MIN(ngpown_, nmtxmax - igp_start + 1)
    endif

    if (ngpown > 0) then
       SAFE_ALLOCATE(data, (SCALARSIZE, nmtxmax, ngpown, 1, 1, 1))
       igp_offset = INDXL2G(1, nb, pool_rank, 0, npes_pool) - 1
    else
       SAFE_ALLOCATE(data, (1, 1, 1, 1, 1, 1))
       igp_offset = -1
    endif
    data = 0.0D0

#ifdef MPI
    call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call h5fopen_f(TRUNC(name), H5F_ACC_RDONLY_F, file_id, error, access_prp = plist_id)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call h5pclose_f(plist_id,error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif
#else
    call h5fopen_f(TRUNC(name), H5F_ACC_RDONLY_F, file_id, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif
#endif

    call h5dopen_f(file_id, 'mats/matrix', data_id, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call h5dget_space_f(data_id,dataspace,error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    count(:) = (/ SCALARSIZE, nmtxmax, MAX(ngpown, 1), 1, 1, 1 /)

    call h5screate_simple_f(6, count, memspace, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    !> Construct data and offset
    if (ngpown > 0) then
       offset(:) = (/ 0, 0, igp_offset, ifreq_zero-1, is-1, iq-1 /)
       !> Select hyperslab
       call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, offset, count, error)
       if (error .ne. 0) then
          call die("HDF5 error", only_root_writes=.true.)
       endif
    else
       call H5sselect_none_f(memspace,error);
       if (error .ne. 0) then
          call die("HDF5 error", only_root_writes=.true.)
       endif
       call H5sselect_none_f(dataspace,error);
       if (error .ne. 0) then
          call die("HDF5 error", only_root_writes=.true.)
       endif
    endif

#ifdef MPI
    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    !> Collectively read the file
    call h5dread_f(data_id, H5T_NATIVE_DOUBLE, data, count, error, memspace, dataspace, xfer_prp = plist_id)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif
#else
    call h5dread_f(data_id, H5T_NATIVE_DOUBLE, data, count, error, memspace, dataspace)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif
#endif

    if (ngpown > 0) then
       eps(1:nmtxmax, 1:ngpown) = SCALARIFY2(data(1,1:nmtxmax,1:ngpown,1,1,1), data(2,1:nmtxmax,1:ngpown,1,1,1))
    endif

    call h5sclose_f(memspace, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call h5sclose_f(dataspace, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call h5dclose_f(data_id, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call h5fclose_f(file_id,error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif
    SAFE_DEALLOCATE(data)

    POP_SUB(read_eps_matrix_par_hdf5_3_)
  end subroutine read_eps_matrix_par_hdf5_3_

  !> Read and distribute epsmat.h5 to all procs (1:peinf%npes)
  subroutine read_eps_matrix_par_hdf5_4(eps, nb, ngpown_, nmtxmax, iq, is, ifreq_zero, name)
    SCALAR, intent(inout) :: eps(:,:) !> (nmtxmax, MAX(ngpown,1))
    integer, intent(in) :: nb, ngpown_ !> block size, not necessarily 1
    integer, intent(in) :: nmtxmax, iq, is, ifreq_zero
    character(len=*), intent(in) :: name

    integer(HID_T) :: file_id       ! File identifier
    integer(HID_T) :: data_id       ! Property list identifier
    integer(HID_T) :: plist_id       ! Property list identifier
    integer(HID_T) :: dataspace        ! Property list identifier
    integer(HID_T) :: memspace        ! Property list identifier
    integer :: error
    integer(HSIZE_T) :: count(6), offset(6) !, countm(6)
    real(DP), allocatable :: data(:,:,:,:,:,:)
    integer :: igp_offset, igp_start, ngpown
    PUSH_SUB(read_eps_matrix_par_hdf5_4)

    !> epsmpi%nb = iceil(neps, peinf%npes)
    !> epsmpi%ngpown = NUMROC(neps, epsmpi%nb, peinf%inode, 0, peinf%npes)
    !> However, nmtxmax can be < neps = MAX(nmtx_eps0mat, nmtx_epsmat)
    !> When nmtxmax < neps, the task cannot read ngpown elements (because the file just does not have the elements)
    igp_start = INDXL2G(1, nb, peinf%inode, 0, peinf%npes)
    if (nmtxmax - igp_start < 0) then
       ngpown = 0
    else
       ngpown = MIN(ngpown_, nmtxmax - igp_start + 1)
    endif

    !> We only distribute ig2 of epsinv(ig1,ig2,iq) over procs within a pool
    if (ngpown > 0) then
       SAFE_ALLOCATE(data, (SCALARSIZE, nmtxmax, ngpown, 1, 1, 1))
       igp_offset = INDXL2G(1, nb, peinf%inode, 0, peinf%npes) - 1
    else
       SAFE_ALLOCATE(data, (1, 1, 1, 1, 1, 1))
       igp_offset = -1
    endif
    data = 0.0D0

#ifdef MPI
    call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call h5fopen_f(TRUNC(name), H5F_ACC_RDONLY_F, file_id, error, access_prp = plist_id)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call h5pclose_f(plist_id,error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif
#else
    call h5fopen_f(TRUNC(name), H5F_ACC_RDONLY_F, file_id, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif
#endif

    call h5dopen_f(file_id, 'mats/matrix', data_id, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call h5dget_space_f(data_id, dataspace, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    count(:) = (/ SCALARSIZE, nmtxmax, MAX(ngpown, 1), 1, 1, 1 /)

    call h5screate_simple_f(6, count, memspace, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    !> Construct data and offset
    if (ngpown > 0) then
       offset(:) = (/ 0, 0, igp_offset, ifreq_zero-1, is-1, iq-1 /)
       !> Select hyperslab
       call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, offset, count, error)
       if (error .ne. 0) then
          call die("HDF5 error", only_root_writes=.true.)
       endif
    else
       call H5sselect_none_f(memspace,error);
       if (error .ne. 0) then
          call die("HDF5 error", only_root_writes=.true.)
       endif
       call H5sselect_none_f(dataspace,error);
       if (error .ne. 0) then
          call die("HDF5 error", only_root_writes=.true.)
       endif
    endif

#ifdef MPI
    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    !> Collectively read the file
    call h5dread_f(data_id, H5T_NATIVE_DOUBLE, data, count, error, memspace, dataspace, xfer_prp = plist_id)
    if (error .ne. 0) then
       write(*,'(I5,A,6I8,A,6I8,A,I8,A,I8)') peinf%inode, " count = ", count, "offset = ", offset, " nmtxmax = ", nmtxmax, " ngpown = ", ngpown
       call die("HDF5 error I", only_root_writes=.true.)
    endif
#else
    call h5dread_f(data_id, H5T_NATIVE_DOUBLE, data, count, error, memspace, dataspace)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif
#endif
    if (ngpown > 0) then
       eps(1:nmtxmax, 1:ngpown) = SCALARIFY2(data(1,1:nmtxmax,1:ngpown,1,1,1), data(2,1:nmtxmax,1:ngpown,1,1,1))
    endif
    SAFE_DEALLOCATE(data)

    call h5sclose_f(memspace, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call h5sclose_f(dataspace, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call h5dclose_f(data_id, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call h5fclose_f(file_id,error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    POP_SUB(read_eps_matrix_par_hdf5_4)
  end subroutine read_eps_matrix_par_hdf5_4

  !> Parallel code, all procs participate
  !> Read epsmat for all q, one freq
  subroutine read_eps_matrix_par_allq_hdf5(eps, nmtxmax, nq, is, igp_start, ngp_loc, ifreq, name)
    SCALAR, intent(inout) :: eps(:,:,:) !< (nmtxmax,MAX(ngp_loc,1),nq)
    integer, intent(in) :: nmtxmax, nq, is, igp_start, ngp_loc
    integer, intent(in) :: ifreq !> use to specify which frequency is the zero frequency
    character(len=*), intent(in) :: name

    integer(HID_T) :: file_id       ! File identifier
    integer(HID_T) :: data_id       ! Property list identifier
    integer(HID_T) :: plist_id       ! Property list identifier
    integer(HID_T) :: dataspace        ! Property list identifier
    integer(HID_T) :: memspace        ! Property list identifier
    integer :: error
    integer(HSIZE_T) :: count(6), offset(6) !, countm(6)
    !> data(imatrix_flavor, ig1, ig2, ifreq, is, iq)
    real(DP), allocatable :: data(:,:,:,:,:,:)
    PUSH_SUB(read_eps_matrix_par_allq_hdf5)

    !> We only distribute ig2 of epsinv(ig1,ig2,iq) over all procs
    SAFE_ALLOCATE(data, (SCALARSIZE, nmtxmax, MAX(ngp_loc,1), 1, 1, nq))

#ifdef MPI
    call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call h5fopen_f(trim(name), H5F_ACC_RDONLY_F, file_id, error, access_prp = plist_id)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call h5pclose_f(plist_id,error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif
#else
    call h5fopen_f(trim(name), H5F_ACC_RDONLY_F, file_id, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif
#endif

    call h5dopen_f(file_id, 'mats/matrix', data_id, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call h5dget_space_f(data_id,dataspace,error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    count(:) = (/ SCALARSIZE, nmtxmax, MAX(ngp_loc,1), 1, 1, nq /)

    call h5screate_simple_f(6, count, memspace, error)
    if (error .ne. 0) then
       call die("1 HDF5 error", only_root_writes=.true.)
    endif

    !> Construct data and offset
    if (ngp_loc > 0) then
       offset(:) = (/ 0, 0, igp_start - 1, ifreq - 1, is - 1, 0 /)
       !> Select hyperslab
       call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, offset, count, error)
       if (error .ne. 0) then
          call die("2 HDF5 error", only_root_writes=.true.)
       endif
    else
       call H5sselect_none_f(memspace,error);
       if (error .ne. 0) then
          call die("3 HDF5 error", only_root_writes=.true.)
       endif
       call H5sselect_none_f(dataspace,error);
       if (error .ne. 0) then
          call die("4 HDF5 error", only_root_writes=.true.)
       endif
    endif

#ifdef MPI
    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
    if (error .ne. 0) then
       call die("5 HDF5 error", only_root_writes=.true.)
    endif

    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
    if (error .ne. 0) then
       call die("6 HDF5 error", only_root_writes=.true.)
    endif

    !> Collectively read the file
    call h5dread_f(data_id, H5T_NATIVE_DOUBLE, data, count, error, memspace, dataspace, xfer_prp = plist_id)
    if (error .ne. 0) then
       call die("7 HDF5 error", only_root_writes=.true.)
    endif
#else
    call h5dread_f(data_id, H5T_NATIVE_DOUBLE, data, count, error, memspace, dataspace)
    if (error .ne. 0) then
       call die("8 HDF5 error", only_root_writes=.true.)
    endif
#endif

    if (ngp_loc > 0) then
       eps(:, :, :) = SCALARIFY2(data(1, :, :, 1, 1, :),data(2, :, :, 1, 1, :))
    endif

    call h5sclose_f(memspace, error)
    if (error .ne. 0) then
       call die("9 HDF5 error", only_root_writes=.true.)
    endif

    call h5sclose_f(dataspace, error)
    if (error .ne. 0) then
       call die("10 HDF5 error", only_root_writes=.true.)
    endif

    call h5dclose_f(data_id, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call h5fclose_f(file_id,error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    SAFE_DEALLOCATE(data)

    POP_SUB(read_eps_matrix_par_allq_hdf5)
  end subroutine read_eps_matrix_par_allq_hdf5

  !> Partition over rq
  subroutine read_eps_matrix_par_distribute_rq_hdf5(eps, nmtxmax, is, irq_start, nrq_loc, ifreq, name)
    SCALAR, intent(inout) :: eps(:,:,:) !< (nmtxmax,nmtxmax,nrq_loc)
    integer, intent(in) :: nmtxmax, is, irq_start, nrq_loc
    integer, intent(in) :: ifreq !> use to specify which frequency is the zero frequency
    character(len=*), intent(in) :: name

    integer(HID_T) :: file_id       ! File identifier
    integer(HID_T) :: data_id       ! Property list identifier
    integer(HID_T) :: plist_id       ! Property list identifier
    integer(HID_T) :: dataspace        ! Property list identifier
    integer(HID_T) :: memspace        ! Property list identifier
    integer :: error
    integer(HSIZE_T) :: count(6), offset(6) !, countm(6)
    !> data(imatrix_flavor, ig1, ig2, ifreq, is, iq)
    real(DP), allocatable :: data(:,:,:,:,:,:)
    PUSH_SUB(read_eps_matrix_par_distribute_rq_hdf5)

    !> We only distribute rq of epsinv(ig1,ig2,iq) over all procs
    SAFE_ALLOCATE(data, (SCALARSIZE, nmtxmax, nmtxmax, 1, 1, MAX(nrq_loc,1)))

#ifdef MPI
    call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call h5fopen_f(trim(name), H5F_ACC_RDONLY_F, file_id, error, access_prp = plist_id)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call h5pclose_f(plist_id,error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif
#else
    call h5fopen_f(trim(name), H5F_ACC_RDONLY_F, file_id, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif
#endif

    call h5dopen_f(file_id, 'mats/matrix', data_id, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call h5dget_space_f(data_id,dataspace,error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    count(:) = (/ SCALARSIZE, nmtxmax, nmtxmax, 1, 1, MAX(nrq_loc,1) /)

    call h5screate_simple_f(6, count, memspace, error)
    if (error .ne. 0) then
       call die("1 HDF5 error", only_root_writes=.true.)
    endif

    !> Construct data and offset
    if (nrq_loc > 0) then
       offset(:) = (/ 0, 0, 0, ifreq - 1, is - 1, irq_start - 1 /)
       !> Select hyperslab
       call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, offset, count, error)
       if (error .ne. 0) then
          call die("2 HDF5 error", only_root_writes=.true.)
       endif
    else
       call H5sselect_none_f(memspace,error);
       if (error .ne. 0) then
          call die("3 HDF5 error", only_root_writes=.true.)
       endif
       call H5sselect_none_f(dataspace,error);
       if (error .ne. 0) then
          call die("4 HDF5 error", only_root_writes=.true.)
       endif
    endif

#ifdef MPI
    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
    if (error .ne. 0) then
       call die("5 HDF5 error", only_root_writes=.true.)
    endif

    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
    if (error .ne. 0) then
       call die("6 HDF5 error", only_root_writes=.true.)
    endif

    !> Collectively read the file
    call h5dread_f(data_id, H5T_NATIVE_DOUBLE, data, count, error, memspace, dataspace, xfer_prp = plist_id)
    if (error .ne. 0) then
       call die("7 HDF5 error", only_root_writes=.true.)
    endif
#else
    call h5dread_f(data_id, H5T_NATIVE_DOUBLE, data, count, error, memspace, dataspace)
    if (error .ne. 0) then
       call die("8 HDF5 error", only_root_writes=.true.)
    endif
#endif

    if (nrq_loc > 0) then
       eps(:, :, :) = SCALARIFY2(data(1, :, :, 1, 1, :),data(2, :, :, 1, 1, :))
    endif

    call h5sclose_f(memspace, error)
    if (error .ne. 0) then
       call die("9 HDF5 error", only_root_writes=.true.)
    endif

    call h5sclose_f(dataspace, error)
    if (error .ne. 0) then
       call die("10 HDF5 error", only_root_writes=.true.)
    endif

    call h5dclose_f(data_id, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call h5fclose_f(file_id,error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    SAFE_DEALLOCATE(data)
    POP_SUB(read_eps_matrix_par_distribute_rq_hdf5)
  end subroutine read_eps_matrix_par_distribute_rq_hdf5

  !==========================================================================================

  subroutine read_eps_matrix_par_f_hdf5(retarded, nb, pool_comm, my_pool, npools, nmtx, nFreq, iq, is, name, advanced)
    complex(DPC), intent(inout) :: retarded(:,:,:) !< (neps,ngpown,nFreq)
    integer, intent(in) :: nb !< block size
    integer, intent(in) :: pool_comm !< MPI comm for each pool
    integer, intent(in) :: my_pool !< my pool, starting from 0 (0=no pools)
    integer, intent(in) :: npools !< number of pools (1=no pools).
    integer, intent(in) :: nmtx
    integer, intent(in) :: nFreq
    integer, intent(in) :: iq
    integer, intent(in) :: is
    character(len=*), intent(in) :: name
    complex(DPC), optional, intent(inout) :: advanced(:,:,:) !< (nFreq,neps,ngpown)
    integer(HID_T) :: file_id       ! File identifier
    integer(HID_T) :: data_id       ! Property list identifier
    integer(HID_T) :: plist_id       ! Property list identifier
    integer(HID_T) :: dataspace        ! Property list identifier
    integer(HID_T) :: memspace        ! Property list identifier
    integer :: error
    integer(HSIZE_T) :: count(6), offset(6)
    real(DP), allocatable :: data(:,:,:,:,:,:)
    integer :: pool_rank !< processor rank for column distribution
    integer :: npes_pool !< number of processors over which we distribute
    integer :: ngpown_max, igp, igp_loc, nmatrix_per_spin, nspin, buf_sz, version
    logical :: want_advanced, read_advanced
    PUSH_SUB(read_eps_matrix_par_f_hdf5)

#ifdef MPI
    call MPI_Comm_rank(pool_comm, pool_rank, mpierr)
    call MPI_Comm_size(pool_comm, npes_pool, mpierr)
    ! FHJ: We need the following BCast for processors left over out of the pool
    call MPI_Bcast(npes_pool, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
    call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call h5fopen_f(trim(name), H5F_ACC_RDONLY_F, file_id, error, access_prp = plist_id)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call h5pclose_f(plist_id,error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif
#else
    pool_rank = 0
    npes_pool = 1
    call h5fopen_f(trim(name), H5F_ACC_RDONLY_F, file_id, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif
#endif
    want_advanced = present(advanced)
    if (want_advanced) then
       call die("epsread: Advanced eps not supported.", only_root_writes=.true.)
    endif
    ! FHJ: the default is never to read the advanced matrix, unless the file
    ! version is <3 (on which case we didn`t store the Coulomb interaction)
    call hdf5_read_int(file_id, 'eps_header/versionnumber', version, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call hdf5_read_int(file_id, 'eps_header/params/nmatrix', nmatrix_per_spin, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call hdf5_read_int(file_id, 'mf_header/kpoints/nspin', nspin, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    nmatrix_per_spin = nmatrix_per_spin / nspin
    read_advanced = .false.
    buf_sz = 1
    if (version<3) then
       call hdf5_read_logical(file_id, 'eps_header/params/has_advanced', read_advanced, error)
       if (error .ne. 0) then
          call die("HDF5 error", only_root_writes=.true.)
       endif
       if (SCALARSIZE==2 .and. .not.read_advanced) then
          call die('Inconsistent epsmat file: version<3, but no advanced matrix', only_root_writes=.true.)
       endif
    endif
    if (read_advanced) buf_sz = 2

    ngpown_max = NUMROC(nmtx, nb, 0, 0, npes_pool)
    SAFE_ALLOCATE(data,(2,nmtx,1,nFreq,buf_sz,1))

    call logit('Reading HDF5 file')
    do igp_loc = 1, ngpown_max
       if (my_pool>=0) then
          igp = INDXL2G(igp_loc, nb, pool_rank, 0, npes_pool)
          !> <- idle proc ->
          !> [DEBUG]
       else
          igp = nmtx + 1
       endif
       call h5dopen_f(file_id, 'mats/matrix', data_id, error)
       if (error .ne. 0) then
          call die("HDF5 error", only_root_writes=.true.)
       endif

       call h5dget_space_f(data_id,dataspace,error)
       if (error .ne. 0) then
          call die("HDF5 error", only_root_writes=.true.)
       endif

       count(:) = (/2, nmtx, 1, nFreq, buf_sz, 1/)
       call h5screate_simple_f(6, count, memspace, error)
       if (error .ne. 0) then
          call die("HDF5 error", only_root_writes=.true.)
       endif

       ! Construct data and offset
       !> [DEBUG]
       !> Use my_pool instead of igp here!
       if (igp <= nmtx) then
          offset = (/0, 0, igp-1, 0, nmatrix_per_spin*(is-1), iq-1/)
          call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, offset, count, error)
          if (error .ne. 0) then
             call die("HDF5 error", only_root_writes=.true.)
          endif
       else
          call H5sselect_none_f(memspace,error);
          if (error .ne. 0) then
             call die("HDF5 error", only_root_writes=.true.)
          endif
          call H5sselect_none_f(dataspace,error);
          if (error .ne. 0) then
             call die("HDF5 error", only_root_writes=.true.)
          endif
       endif

       !> Create property list for collective dataset read
       !> Read is serial for now
#ifdef MPI
       call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
       if (error .ne. 0) then
          call die("HDF5 error", only_root_writes=.true.)
       endif
       call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
       if (error .ne. 0) then
          call die("HDF5 error", only_root_writes=.true.)
       endif
       call h5dread_f(data_id, H5T_NATIVE_DOUBLE, data, count, error, memspace, dataspace, xfer_prp = plist_id)
       if (error .ne. 0) then
          call die("HDF5 error", only_root_writes=.true.)
       endif
#else
       call h5dread_f(data_id, H5T_NATIVE_DOUBLE, data, count, error, memspace, dataspace)
       if (error .ne. 0) then
          call die("HDF5 error", only_root_writes=.true.)
       endif
#endif
       if (igp <= nmtx) then
          retarded(1:nmtx,igp_loc,1:nFreq) = DCMPLX(data(1,1:nmtx,1,1:nFreq,1,1),data(2,1:nmtx,1,1:nFreq,1,1))
          if (want_advanced .and. read_advanced) then
             advanced(1:nmtx,igp_loc,1:nFreq) = DCMPLX(data(1,1:nmtx,1,1:nFreq,2,1),data(2,1:nmtx,1,1:nFreq,2,1))
          endif
       endif

       call h5sclose_f(memspace, error)
       if (error .ne. 0) then
          call die("HDF5 error", only_root_writes=.true.)
       endif

       call h5sclose_f(dataspace, error)
       if (error .ne. 0) then
          call die("HDF5 error", only_root_writes=.true.)
       endif

       call h5dclose_f(data_id, error)
       if (error .ne. 0) then
          call die("HDF5 error", only_root_writes=.true.)
       endif
    enddo

    SAFE_DEALLOCATE(data)
    if (want_advanced .and. .not.read_advanced) then
       call die("epsread: Advanced eps not supported.", only_root_writes=.true.)
    endif
    call h5fclose_f(file_id,error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    POP_SUB(read_eps_matrix_par_f_hdf5)
  end subroutine read_eps_matrix_par_f_hdf5

  subroutine read_eps_matrix_par_f_hdf5_2(retarded, nb, pool_rank, npes_pool, my_pool, nmtx, nmtxmax, iq, is, nFreq, name, advanced)
    complex(DPC), intent(inout) :: retarded(:,:,:) !< (neps,ngpown,nFreq)
    integer, intent(in) :: nb !< block size
    ! integer, intent(in) :: pool_comm !< MPI comm for each pool
    !> my_pool = -1 ==> idle proc
    integer, intent(in) :: pool_rank, npes_pool, my_pool
    integer, intent(in) :: nmtx, nmtxmax
    integer, intent(in) :: iq, is, nFreq
    character(len=*), intent(in) :: name
    complex(DPC), optional, intent(inout) :: advanced(:,:,:) !< (nFreq,neps,ngpown)
    integer(HID_T) :: file_id       ! File identifier
    integer(HID_T) :: data_id       ! Property list identifier
    integer(HID_T) :: plist_id       ! Property list identifier
    integer(HID_T) :: dataspace        ! Property list identifier
    integer(HID_T) :: memspace        ! Property list identifier
    integer :: error
    integer(HSIZE_T) :: count(6), offset(6)
    real(DP), allocatable :: data(:,:,:,:,:,:)
    integer :: ngpown, ngpown_max, igp, igp_loc, nmatrix_per_spin, nspin, buf_sz, version
    logical :: want_advanced, read_advanced
    PUSH_SUB(read_eps_matrix_par_f_hdf5_2)

#ifdef MPI
    ngpown     = NUMROC(nmtxmax, nb, pool_rank, 0, npes_pool)
    ngpown_max = NUMROC(nmtxmax, nb,         0, 0, npes_pool)

    call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call h5fopen_f(trim(name), H5F_ACC_RDONLY_F, file_id, error, access_prp = plist_id)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call h5pclose_f(plist_id,error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif
#else
    ngpown = nmtxmax
    ngpown_max = nmtxmax
    call h5fopen_f(trim(name), H5F_ACC_RDONLY_F, file_id, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif
#endif

    want_advanced = present(advanced)
    if (want_advanced) then
       call die("epsread: Advanced eps not supported.", only_root_writes=.true.)
    endif
    ! FHJ: the default is never to read the advanced matrix, unless the file
    ! version is <3 (on which case we didn`t store the Coulomb interaction)
    call hdf5_read_int(file_id, 'eps_header/versionnumber', version, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call hdf5_read_int(file_id, 'eps_header/params/nmatrix', nmatrix_per_spin, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call hdf5_read_int(file_id, 'mf_header/kpoints/nspin', nspin, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    nmatrix_per_spin = nmatrix_per_spin / nspin
    read_advanced = .false.
    buf_sz = 1
    if (version<3) then
       call hdf5_read_logical(file_id, 'eps_header/params/has_advanced', read_advanced, error)
       if (error .ne. 0) then
          call die("HDF5 error", only_root_writes=.true.)
       endif
       if (SCALARSIZE==2 .and. .not.read_advanced) then
          call die('Inconsistent epsmat file: version<3, but no advanced matrix', only_root_writes=.true.)
       endif
    endif
    if (read_advanced) buf_sz = 2

    SAFE_ALLOCATE(data, (2, nmtxmax, 1, nFreq, buf_sz, 1))

    !> [Q] Is npes_pool a global variable?
    !> [A] Yes. peinf%npes_pool = peinf%npes/peinf%npools
    !> [Q] Is it really necessary to include my_pool in the argument list?
    !> [A] Yes. Idle procs do not belong to any pool, and therefore should not read epsmat.
    call logit('Reading HDF5 file')

    do igp_loc = 1, ngpown_max
       if (my_pool>=0) then
          igp = INDXL2G(igp_loc, nb, pool_rank, 0, npes_pool)
          !> <- idle proc ->
       else
          igp = nmtx + 1
       endif

       call h5dopen_f(file_id, 'mats/matrix', data_id, error)
       if (error .ne. 0) then
          call die("HDF5 error", only_root_writes=.true.)
       endif

       call h5dget_space_f(data_id, dataspace, error)
       if (error .ne. 0) then
          call die("HDF5 error", only_root_writes=.true.)
       endif

       count(:) = (/2, nmtxmax, 1, nFreq, buf_sz, 1/)
       call h5screate_simple_f(6, count, memspace, error)
       if (error .ne. 0) then
          call die("HDF5 error", only_root_writes=.true.)
       endif

       if ( (igp_loc <= ngpown) .and. (igp <= nmtx)) then
          offset = (/0, 0, igp-1, 0, nmatrix_per_spin*(is-1), iq-1/)
          call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, offset, count, error)
          if (error .ne. 0) then
             call die("HDF5 error", only_root_writes=.true.)
          endif
       else
          call H5sselect_none_f(memspace,error);
          if (error .ne. 0) then
             call die("HDF5 error", only_root_writes=.true.)
          endif
          call H5sselect_none_f(dataspace,error);
          if (error .ne. 0) then
             call die("HDF5 error", only_root_writes=.true.)
          endif
       endif

       !> Create property list for collective dataset read
       !> Read is serial for now
#ifdef MPI
       call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
       if (error .ne. 0) then
          call die("HDF5 error", only_root_writes=.true.)
       endif
       call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
       if (error .ne. 0) then
          call die("HDF5 error", only_root_writes=.true.)
       endif
       call h5dread_f(data_id, H5T_NATIVE_DOUBLE, data, count, error, memspace, dataspace, xfer_prp = plist_id)
       if (error .ne. 0) then
          call die("HDF5 error", only_root_writes=.true.)
       endif
#else
       call h5dread_f(data_id, H5T_NATIVE_DOUBLE, data, count, error, memspace, dataspace)
       if (error .ne. 0) then
          call die("HDF5 error", only_root_writes=.true.)
       endif
#endif
       if ( (igp_loc <= ngpown) .and. (igp <= nmtx)) then
          retarded(1:nmtx,igp_loc,1:nFreq) = SCALARIFY2(data(1,1:nmtx,1,1:nFreq,1,1),data(2,1:nmtx,1,1:nFreq,1,1))
          if (want_advanced .and. read_advanced) then
             call die("advanced not supported.", only_root_writes=.true.)
             advanced(1:nmtx,igp_loc,1:nFreq) = SCALARIFY2(data(1,1:nmtx,1,1:nFreq,2,1),data(2,1:nmtx,1,1:nFreq,2,1))
          endif
       endif

       call h5sclose_f(memspace, error)
       if (error .ne. 0) then
          call die("HDF5 error", only_root_writes=.true.)
       endif

       call h5sclose_f(dataspace, error)
       if (error .ne. 0) then
          call die("HDF5 error", only_root_writes=.true.)
       endif

       call h5dclose_f(data_id, error)
       if (error .ne. 0) then
          call die("HDF5 error", only_root_writes=.true.)
       endif
    enddo

    SAFE_DEALLOCATE(data)
    if (want_advanced .and. .not.read_advanced) then
       call die("epsread: Advanced eps not supported.", only_root_writes=.true.)
       ! call get_advanced_from_retarded()
    endif
    call h5fclose_f(file_id,error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    POP_SUB(read_eps_matrix_par_f_hdf5_2)
  end subroutine read_eps_matrix_par_f_hdf5_2

  subroutine read_eps_matrix_par_f_hdf5_3(retarded, nb, ngpown, ngpown_max, pool_rank, npes_pool, my_pool, nmtx, nmtxmax, iq, is, nFreq, name)
    complex(DPC), intent(inout) :: retarded(:,:,:) !< (neps,ngpown,nFreq)
    integer, intent(in) :: nb, ngpown, ngpown_max !< block size
    ! integer, intent(in) :: pool_comm !< MPI comm for each pool
    !> my_pool = -1 ==> idle proc
    integer, intent(in) :: pool_rank, npes_pool, my_pool
    integer, intent(in) :: nmtx, nmtxmax, iq, is, nFreq
    character(len=*), intent(in) :: name
    integer(HID_T) :: file_id       ! File identifier
    integer(HID_T) :: data_id       ! Property list identifier
    integer(HID_T) :: plist_id       ! Property list identifier
    integer(HID_T) :: dataspace        ! Property list identifier
    integer(HID_T) :: memspace        ! Property list identifier
    integer :: error
    integer(HSIZE_T) :: count(6), offset(6)
    real(DP), allocatable :: data(:,:,:,:,:,:)
    integer :: igp, igp_loc, nmatrix_per_spin, nspin, buf_sz, version
    PUSH_SUB(read_eps_matrix_par_f_hdf5_3)

#ifdef MPI
    call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call h5fopen_f(trim(name), H5F_ACC_RDONLY_F, file_id, error, access_prp = plist_id)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call h5pclose_f(plist_id,error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif
#else
    call h5fopen_f(trim(name), H5F_ACC_RDONLY_F, file_id, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif
#endif

    call hdf5_read_int(file_id, 'eps_header/versionnumber', version, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call hdf5_read_int(file_id, 'eps_header/params/nmatrix', nmatrix_per_spin, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call hdf5_read_int(file_id, 'mf_header/kpoints/nspin', nspin, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    nmatrix_per_spin = nmatrix_per_spin / nspin
    buf_sz = 1
    SAFE_ALLOCATE(data, (2, nmtxmax, 1, nFreq, buf_sz, 1))
    data = 0.0D0

    !> [Q] Is npes_pool a global variable?
    !> [A] Yes. peinf%npes_pool = peinf%npes/peinf%npools
    !> [Q] Is it really necessary to include my_pool in the argument list?
    !> [A] Yes. Idle procs do not belong to any pool, and therefore should not read epsmat.
    call logit('Reading HDF5 file')

    call h5dopen_f(file_id, 'mats/matrix', data_id, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call h5dget_space_f(data_id, dataspace, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    do igp_loc = 1, ngpown_max
       if (my_pool>=0) then
          igp = INDXL2G(igp_loc, nb, pool_rank, 0, npes_pool)
          !> <- idle proc ->
       else
          igp = nmtx + 1
       endif

       count(:) = (/2, nmtxmax, 1, nFreq, buf_sz, 1/)
       call h5screate_simple_f(6, count, memspace, error)
       if (error .ne. 0) then
          call die("HDF5 error", only_root_writes=.true.)
       endif

       if ( (igp_loc <= ngpown) .and. (igp <= nmtx)) then
          offset = (/0, 0, igp-1, 0, nmatrix_per_spin*(is-1), iq-1/)
          call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, offset, count, error)
          if (error .ne. 0) then
             call die("HDF5 error", only_root_writes=.true.)
          endif
       else
          call H5sselect_none_f(memspace,error);
          if (error .ne. 0) then
             call die("HDF5 error", only_root_writes=.true.)
          endif
          call H5sselect_none_f(dataspace,error);
          if (error .ne. 0) then
             call die("HDF5 error", only_root_writes=.true.)
          endif
       endif

       !> Create property list for collective dataset read
       !> Read is serial for now
#ifdef MPI
       call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
       if (error .ne. 0) then
          call die("HDF5 error", only_root_writes=.true.)
       endif
       call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
       if (error .ne. 0) then
          call die("HDF5 error", only_root_writes=.true.)
       endif
       call h5dread_f(data_id, H5T_NATIVE_DOUBLE, data, count, error, memspace, dataspace, xfer_prp = plist_id)
       if (error .ne. 0) then
          call die("HDF5 error", only_root_writes=.true.)
       endif
#else
       call h5dread_f(data_id, H5T_NATIVE_DOUBLE, data, count, error, memspace, dataspace)
       if (error .ne. 0) then
          call die("HDF5 error", only_root_writes=.true.)
       endif
#endif
       if ( (igp_loc <= ngpown) .and. (igp <= nmtx)) then
          retarded(1:nmtx, igp_loc, 1:nFreq) = SCALARIFY2(data(1, 1:nmtx, 1, 1:nFreq, 1, 1), data(2, 1:nmtx, 1, 1:nFreq, 1, 1))
       endif

       call h5sclose_f(memspace, error)
       if (error .ne. 0) then
          call die("HDF5 error", only_root_writes=.true.)
       endif
    enddo

    SAFE_DEALLOCATE(data)

    call h5sclose_f(dataspace, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call h5dclose_f(data_id, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call h5fclose_f(file_id,error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    POP_SUB(read_eps_matrix_par_f_hdf5_3)
  end subroutine read_eps_matrix_par_f_hdf5_3

  subroutine read_eps_matrix_par_f_hdf5_3_(epsDyn, nb, ngpown_, ngpown_max, pool_rank, npes_pool, my_pool, nmtx, nmtxmax, iq, is, nFreq, name)
    complex(DPC), intent(inout) :: epsDyn(:,:,:) !< (neps,ngpown,nFreq)
    integer, intent(in) :: nb, ngpown_, ngpown_max !< block size
    integer, intent(in) :: pool_rank, npes_pool, my_pool
    integer, intent(in) :: nmtx, nmtxmax, iq, is, nFreq
    character(len=*), intent(in) :: name
    integer(HID_T) :: file_id       ! File identifier
    integer(HID_T) :: data_id       ! Property list identifier
    integer(HID_T) :: plist_id       ! Property list identifier
    integer(HID_T) :: dataspace        ! Property list identifier
    integer(HID_T) :: memspace        ! Property list identifier
    integer :: error
    integer(HSIZE_T) :: count(6), offset(6)
    real(DP), allocatable :: data(:,:,:,:,:,:)
    integer :: igp_offset, igp_start, ngpown
    PUSH_SUB(read_eps_matrix_par_f_hdf5_3_)

    igp_start = INDXL2G(1, nb, pool_rank, 0, npes_pool)
    if (nmtxmax - igp_start < 0) then
       ngpown = 0
    else
       ngpown = MIN(ngpown_, nmtxmax - igp_start + 1)
    endif

#ifdef MPI
    call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call h5fopen_f(trim(name), H5F_ACC_RDONLY_F, file_id, error, access_prp = plist_id)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call h5pclose_f(plist_id,error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif
#else
    call h5fopen_f(trim(name), H5F_ACC_RDONLY_F, file_id, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif
#endif

    if (ngpown > 0) then
       SAFE_ALLOCATE(data, (2, nmtxmax, ngpown, nfreq, 1, 1))
       igp_offset = INDXL2G(1, nb, pool_rank, 0, npes_pool) - 1
    else
       SAFE_ALLOCATE(data, (1, 1, 1, 1, 1, 1))
       igp_offset = -1
    endif
    data = 0.0D0

    call logit('Reading HDF5 file')

    call h5dopen_f(file_id, 'mats/matrix', data_id, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call h5dget_space_f(data_id, dataspace, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    count(:) = (/2, nmtxmax, MAX(ngpown, 1), nfreq, 1, 1/)
    call h5screate_simple_f(6, count, memspace, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    if (ngpown > 0) then
       offset = (/0, 0, igp_offset, 0, is-1, iq-1/)
       call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, offset, count, error)
       if (error .ne. 0) then
          call die("HDF5 error", only_root_writes=.true.)
       endif
    else
       call H5sselect_none_f(memspace,error);
       if (error .ne. 0) then
          call die("HDF5 error", only_root_writes=.true.)
       endif
       call H5sselect_none_f(dataspace,error);
       if (error .ne. 0) then
          call die("HDF5 error", only_root_writes=.true.)
       endif
    endif

    !> Create property list for collective dataset read
    !> Read is serial for now
#ifdef MPI
    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif
    call h5dread_f(data_id, H5T_NATIVE_DOUBLE, data, count, error, memspace, dataspace, xfer_prp = plist_id)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif
#else
    call h5dread_f(data_id, H5T_NATIVE_DOUBLE, data, count, error, memspace, dataspace)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif
#endif

    if (ngpown > 0) then
       epsDyn(1:nmtxmax,1:ngpown,1:nfreq) = SCALARIFY2(data(1,1:nmtxmax,1:ngpown,1:nfreq,1,1), data(2,1:nmtxmax,1:ngpown,1:nfreq,1,1))
    endif
    SAFE_DEALLOCATE(data)

    call h5sclose_f(memspace, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call h5sclose_f(dataspace, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call h5dclose_f(data_id, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call h5fclose_f(file_id,error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    POP_SUB(read_eps_matrix_par_f_hdf5_3_)
  end subroutine read_eps_matrix_par_f_hdf5_3_

  !> Read and distribute epsmat.h5 using all procs (1:peinf%npes)
  subroutine read_eps_matrix_par_f_hdf5_4(epsDyn, nb, ngpown_, nmtxmax, iq, is, nFreq, name)
    complex(DPC), intent(inout) :: epsDyn(:,:,:) !< (neps,ngpown,nFreq)
    integer, intent(in) :: nb, ngpown_, nmtxmax, iq, is, nFreq
    character(len=*), intent(in) :: name
    integer(HID_T) :: file_id       ! File identifier
    integer(HID_T) :: data_id       ! Property list identifier
    integer(HID_T) :: plist_id       ! Property list identifier
    integer(HID_T) :: dataspace        ! Property list identifier
    integer(HID_T) :: memspace        ! Property list identifier
    integer :: error, igp_offset
    integer(HSIZE_T) :: count(6), offset(6)
    real(DP), allocatable :: data(:,:,:,:,:,:)
    integer :: igp_start, ngpown
    PUSH_SUB(read_eps_matrix_par_f_hdf5_4)

    igp_start = INDXL2G(1, nb, peinf%inode, 0, peinf%npes)
    if (nmtxmax - igp_start < 0) then
       ngpown = 0
    else
       ngpown = MIN(ngpown_, nmtxmax - igp_start + 1)
    endif

#ifdef MPI
    call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call h5fopen_f(trim(name), H5F_ACC_RDONLY_F, file_id, error, access_prp = plist_id)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call h5pclose_f(plist_id,error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif
#else
    call h5fopen_f(trim(name), H5F_ACC_RDONLY_F, file_id, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif
#endif

    if (ngpown > 0) then
       SAFE_ALLOCATE(data, (2, nmtxmax, ngpown, nFreq, 1, 1))
       igp_offset = INDXL2G(1, nb, peinf%inode, 0, peinf%npes) - 1
    else
       SAFE_ALLOCATE(data, (1, 1, 1, 1, 1, 1))
       igp_offset = -1
    endif
    data = 0.0D0

    call logit('Reading HDF5 file')

    call h5dopen_f(file_id, 'mats/matrix', data_id, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call h5dget_space_f(data_id, dataspace, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    count(:) = (/2, nmtxmax, MAX(ngpown, 1), nFreq, 1, 1/)
    call h5screate_simple_f(6, count, memspace, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    if ( ngpown > 0 ) then
       offset = (/0, 0, igp_offset, 0, is-1, iq-1/)
       call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, offset, count, error)
       if (error .ne. 0) then
          call die("HDF5 error", only_root_writes=.true.)
       endif
    else
       call H5sselect_none_f(memspace, error);
       if (error .ne. 0) then
          call die("HDF5 error", only_root_writes=.true.)
       endif
       call H5sselect_none_f(dataspace, error);
       if (error .ne. 0) then
          call die("HDF5 error", only_root_writes=.true.)
       endif
    endif

    !> Create property list for collective dataset read
#ifdef MPI
    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif
    call h5dread_f(data_id, H5T_NATIVE_DOUBLE, data, count, error, memspace, dataspace, xfer_prp = plist_id)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif
#else
    call h5dread_f(data_id, H5T_NATIVE_DOUBLE, data, count, error, memspace, dataspace)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif
#endif

    if ( ngpown > 0) then
       epsDyn(1:nmtxmax,1:ngpown,1:nfreq) = SCALARIFY2(data(1,1:nmtxmax,1:ngpown,1:nfreq,1,1), data(2,1:nmtxmax,1:ngpown,1:nfreq,1,1))
    endif
    SAFE_DEALLOCATE(data)

    call h5sclose_f(memspace, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call h5sclose_f(dataspace, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call h5dclose_f(data_id, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call h5fclose_f(file_id,error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    POP_SUB(read_eps_matrix_par_f_hdf5_4)
  end subroutine read_eps_matrix_par_f_hdf5_4

#endif
end module epsread_hdf5_m

#include "f_defs.h"

!>=========================================================================
!!
!!  Module:
!!
!!  epsread_hdf5_m     Originally by JRD     Last Modified by Meng Wu
!!
!!    Routines to read header info and matrices from epsmat files in
!!    HDF5 format.
!!
!!=========================================================================

module epsread_hdf5_m
  use global_m
  use hdf5
  use hdf5_io_m
  implicit none
  private
  public :: read_eps_grid_sizes_hdf5, read_eps_freqgrid_hdf5, &
       read_eps_matrix_flavor_hdf5, read_eps_params_hdf5, &
       read_eps_qgrid_hdf5, read_eps_matrix_ser_allq_hdf5, &
       read_eps_matrix_par_distribute_rq_hdf5

contains

  !> Read eps_header/params into pol
  subroutine read_eps_params_hdf5(name, pol)
    character(len=*), intent(in) :: name
    type(polarizability), intent(inout) :: pol
    integer(HID_T) :: file_id       ! File identifier
    integer :: error
    logical :: exists
    PUSH_SUB(read_eps_params_hdf5)

    call h5fopen_f(TRIM(name), H5F_ACC_RDONLY_F, file_id, error)

    call hdf5_read_double(file_id, 'eps_header/params/ecuts', pol%ecuts, error)
    call hdf5_read_double(file_id, 'eps_header/params/efermi', pol%efermi, error)
    pol%efermi = pol%efermi*ryd
    call hdf5_read_logical(file_id, 'eps_header/params/has_advanced', pol%has_advanced, error)
    call hdf5_read_int(file_id, 'eps_header/params/icutv', pol%icutv, error)
    call hdf5_read_int(file_id, 'eps_header/params/intraband_flag', pol%intraband_flag, error)
    call hdf5_read_double(file_id, 'eps_header/params/intraband_overlap_min', pol%intraband_overlap_min, error)
    call hdf5_read_int(file_id, 'eps_header/params/matrix_flavor', pol%matrix_flavor, error)
    call hdf5_read_int(file_id, 'eps_header/params/matrix_type', pol%matrix_type, error)
    call hdf5_read_int(file_id, 'eps_header/params/nband', pol%nband, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call hdf5_read_int(file_id, 'eps_header/params/nmatrix', pol%nmatrix, error)
    call hdf5_read_logical(file_id, 'eps_header/params/subsample', pol%subsample, error)
    call hdf5_read_logical(file_id, 'eps_header/params/subspace', pol%subspace, error)
    call h5lexists_f(file_id, 'eps_header/params/timeordered', exists, error)
    ! call hdf5_read_int(file_id, 'eps_header/params/ncb', pol%ncb, error)    
    ! call hdf5_read_int(file_id, 'eps_header/params/nvb', pol%nvb, error)
    ! call hdf5_read_int(file_id, 'eps_header/params/skip_ncb', pol%skip_ncb, error)
    ! call hdf5_read_int(file_id, 'eps_header/params/skip_nvb', pol%skip_nvb, error)
    ! call hdf5_read_int_array(file_id, 'eps_header/params/FFTgrid', (/3/), pol%FFTgrid, error)
    
    if (exists .and. error == 0) then
       call hdf5_read_logical(file_id, 'eps_header/params/timeordered', pol%timeordered, error)
    else
       pol%timeordered = .true.
    endif

    call h5fclose_f(file_id, error)

    POP_SUB(read_eps_params_hdf5)
  end subroutine read_eps_params_hdf5

  subroutine read_eps_matrix_flavor_hdf5(matrix_flavor, name)
    integer, intent(out) :: matrix_flavor
    character(len=*), intent(in) :: name
    integer(HID_T) :: file_id       ! File identifier
    integer :: error
    PUSH_SUB(read_eps_matrix_flavor_hdf5)

    call h5fopen_f(trim(name), H5F_ACC_RDONLY_F, file_id, error)
    call hdf5_read_int(file_id, 'eps_header/params/matrix_flavor', matrix_flavor, error)
    call h5fclose_f(file_id, error)

    POP_SUB(read_eps_matrix_flavor_hdf5)
  end subroutine read_eps_matrix_flavor_hdf5

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

    call h5fopen_f(trim(name), H5F_ACC_RDONLY_F, file_id, error)

    call hdf5_require_version(file_id, 'eps_header/versionnumber', 2, trim(name))
    call hdf5_require_version(file_id, 'mf_header/versionnumber', VER_WFN_HDF5, trim(name))
    call hdf5_require_flavor(file_id, 'eps_header/flavor', SCALARSIZE, trim(name))
    call hdf5_require_version(file_id, 'mf_header/flavor', SCALARSIZE, trim(name))
    call hdf5_read_int(file_id, 'eps_header/qpoints/nq', nq, error)
    call hdf5_read_int(file_id, 'mf_header/gspace/ng', ng, error)
    call hdf5_read_int(file_id, 'eps_header/freqs/nfreq', nfreq, error)
    call h5lexists_f(file_id, 'eps_header/freqs/nfreq_imag', exists, error)

    if (exists.and.error==0) then
       call hdf5_read_int(file_id, 'eps_header/freqs/nfreq_imag', nfreq_imag, error)
    else
       nfreq_imag = 0
    endif

    call hdf5_read_int(file_id, 'eps_header/freqs/freq_dep', freq_dep, error)
    call hdf5_read_int(file_id, 'eps_header/gspace/nmtx_max', nmtxmax, error)
    call hdf5_read_int_array(file_id, 'eps_header/qpoints/qgrid', (/3/), qgrid, error)
    call hdf5_read_double(file_id, 'eps_header/params/ecuts', ecuts, error)
    call h5fclose_f(file_id,error)

    POP_SUB(read_eps_grid_sizes_hdf5)
  end subroutine read_eps_grid_sizes_hdf5

  subroutine read_eps_qgrid_hdf5(nq, qpts, nmtx, name)
    integer, intent(in) :: nq
    real(DP), intent(inout) :: qpts(:,:) !< (3,nq)
    integer, intent(out) :: nmtx(:) !< (nq)
    character(len=*), intent(in) :: name
    integer(HID_T) :: file_id       ! File identifier
    integer :: error
    PUSH_SUB(read_eps_qgrid_hdf5)

    call h5fopen_f(trim(name), H5F_ACC_RDONLY_F, file_id, error)
    call hdf5_read_double_array(file_id, 'eps_header/qpoints/qpts', (/3,nq/), qpts, error)
    call hdf5_read_int_array(file_id, 'eps_header/gspace/nmtx', (/nq/), nmtx, error)
    call h5fclose_f(file_id, error)

    POP_SUB(read_eps_qgrid_hdf5)
  end subroutine read_eps_qgrid_hdf5

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
    call hdf5_read_double_array(file_id, 'eps_header/freqs/freqs', (/2,nfreq/), freqs_tmp, error)
    do iw=1,nfreq
       dFreqGrid(iw) = freqs_tmp(1,iw)
       dFreqBrd(iw) = DCMPLX(0,freqs_tmp(2,iw))
    enddo
    call h5fclose_f(file_id,error)

    POP_SUB(read_eps_freqgrid_hdf5)
  end subroutine read_eps_freqgrid_hdf5

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
    call h5dopen_f(file_id, 'mats/matrix', data_id, error)
    call h5dget_space_f(data_id,dataspace,error)
    rank = 6
    count(:) = (/ SCALARSIZE, nmtxmax, nmtxmax, 1, 1, nq /)

    call h5screate_simple_f(rank, count, memspace, error)
    offset(:) = (/ 0, 0, 0, ifreq-1, is-1, 0 /)
    call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, offset, count, error)

    SAFE_ALLOCATE(data,(count(1),count(2),count(3),count(4),count(5),count(6)))
    call h5dread_f(data_id, H5T_NATIVE_DOUBLE, data, count, error, memspace, dataspace)
    eps(:,:,:) = SCALARIFY2(data(1,:,:,1,1,:),data(2,:,:,1,1,:))
    SAFE_DEALLOCATE(data)

    call h5sclose_f(memspace, error)
    call h5sclose_f(dataspace, error)
    call h5dclose_f(data_id,error)
    call h5fclose_f(file_id,error)

    POP_SUB(read_eps_matrix_ser_allq_hdf5)
  end subroutine read_eps_matrix_ser_allq_hdf5

  
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
    call h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, error)
    call h5fopen_f(trim(name), H5F_ACC_RDONLY_F, file_id, error, access_prp = plist_id)
    call h5pclose_f(plist_id,error)
#else
    call h5fopen_f(trim(name), H5F_ACC_RDONLY_F, file_id, error)
#endif

    call h5dopen_f(file_id, 'mats/matrix', data_id, error)
    call h5dget_space_f(data_id,dataspace,error)

    count(:) = (/ SCALARSIZE, nmtxmax, nmtxmax, 1, 1, MAX(nrq_loc,1) /)

    call h5screate_simple_f(6, count, memspace, error)

    !> Construct data and offset
    if (nrq_loc > 0) then
       offset(:) = (/ 0, 0, 0, ifreq - 1, is - 1, irq_start - 1 /)
       !> Select hyperslab
       call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, offset, count, error)
    else
       call H5sselect_none_f(memspace,error);
       call H5sselect_none_f(dataspace,error);
    endif

#ifdef MPI
    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
    !> Collectively read the file
    call h5dread_f(data_id, H5T_NATIVE_DOUBLE, data, count, error, memspace, dataspace, xfer_prp = plist_id)
#else
    call h5dread_f(data_id, H5T_NATIVE_DOUBLE, data, count, error, memspace, dataspace)
#endif

    if (nrq_loc > 0) then
       eps(:, :, :) = SCALARIFY2(data(1, :, :, 1, 1, :),data(2, :, :, 1, 1, :))
    endif

    call h5sclose_f(memspace, error)
    call h5sclose_f(dataspace, error)
    call h5dclose_f(data_id, error)
    call h5fclose_f(file_id,error)

    SAFE_DEALLOCATE(data)
    POP_SUB(read_eps_matrix_par_distribute_rq_hdf5)
  end subroutine read_eps_matrix_par_distribute_rq_hdf5

end module epsread_hdf5_m

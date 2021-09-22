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
#ifdef HDF5
  use global_m
  use hdf5
  use hdf5_io_m
  implicit none
  private
  public :: read_eps_grid_sizes_hdf5, &
       read_eps_freqgrid_hdf5, &
       read_eps_matrix_flavor_hdf5, &
       read_eps_params_hdf5, &
       read_eps_qgrid_hdf5

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

    call hdf5_read_int_array(file_id, 'eps_header/params/FFTgrid', (/3/), pol%FFTgrid, error)
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
    call hdf5_read_int(file_id, 'eps_header/params/ncb', pol%ncb, error)
    if (error .ne. 0) then
       call die("HDF5 error", only_root_writes=.true.)
    endif

    call hdf5_read_int(file_id, 'eps_header/params/nmatrix', pol%nmatrix, error)
    call hdf5_read_int(file_id, 'eps_header/params/nvb', pol%nvb, error)
    call hdf5_read_int(file_id, 'eps_header/params/skip_ncb', pol%skip_ncb, error)
    call hdf5_read_int(file_id, 'eps_header/params/skip_nvb', pol%skip_nvb, error)
    call hdf5_read_logical(file_id, 'eps_header/params/subsample', pol%subsample, error)
    call hdf5_read_logical(file_id, 'eps_header/params/subspace', pol%subspace, error)
    call h5lexists_f(file_id, 'eps_header/params/timeordered', exists, error)

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

#endif
end module epsread_hdf5_m

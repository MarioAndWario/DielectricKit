!>=========================================================================
!!
!!  Module:
!!
!!  epswrite_hdf5_m     Originally by JRD     Last Modified by Meng Wu
!!
!!    Routines to write header info for epsmat files in HDF5 format.
!!
!!=========================================================================

#include "f_defs.h"

module epswrite_hdf5_m
  use h5lt
  use hdf5
  use global_m
  use hdf5_io_m
  use wfn_io_hdf5_m
  implicit none
  private
  public :: set_qpt_done, is_qpt_done, eps_hdf5_setup
contains
  
  subroutine set_qpt_done(fname, iq)
    character(len=*), intent(in) :: fname
    integer, intent(in) :: iq

    integer(HID_T) :: file_id
    integer :: nq, error
    logical, allocatable :: qpt_done(:)
    PUSH_SUB(set_qpt_done)

    call open_file(99, trim(fname), status='old')
    call close_file(99)

    call h5fopen_f(trim(fname), H5F_ACC_RDWR_F, file_id, error)
    call hdf5_read_int(file_id, 'eps_header/qpoints/nq', nq, error)
    SAFE_ALLOCATE(qpt_done, (nq))
    call hdf5_read_logical_array(file_id, 'eps_header/qpoints/qpt_done', (/nq/), qpt_done, error)
    qpt_done(iq) = .true.
    call hdf5_write_logical_array(file_id, 'eps_header/qpoints/qpt_done', (/nq/), qpt_done, error)

    call h5fclose_f(file_id, error)

    SAFE_DEALLOCATE(qpt_done)
    POP_SUB(set_qpt_done)
  end subroutine set_qpt_done

  logical function is_qpt_done(fname, iq)
    character(len=*), intent(in) :: fname
    integer, intent(in) :: iq
    integer(HID_T) :: file_id
    integer :: nq, error
    logical, allocatable :: qpt_done(:)
    PUSH_SUB(is_qpt_done)

    call open_file(99, trim(fname), status='old')
    call close_file(99)

    call h5fopen_f(trim(fname), H5F_ACC_RDONLY_F, file_id, error)
    call hdf5_read_int(file_id, 'eps_header/qpoints/nq', nq, error)
    SAFE_ALLOCATE(qpt_done, (nq))
    call hdf5_read_logical_array(file_id, 'eps_header/qpoints/qpt_done', (/nq/), qpt_done, error)
    is_qpt_done = qpt_done(iq)
    call h5fclose_f(file_id, error)
    
    SAFE_DEALLOCATE(qpt_done)
    POP_SUB(is_qpt_done)
  end function is_qpt_done

  subroutine eps_hdf5_setup(kp, gvec, syms, crys, pol, name, restart)
    type(kpoints), intent(in) :: kp
    type(gspace), intent(in) :: gvec
    type(symmetry), intent(in) :: syms
    type(crystal), intent(in) :: crys
    type(polarizability), intent(in) :: pol
    character(len=*), intent(in) :: name
    logical, intent(inout), optional :: restart
    integer :: qgrid(3)
    integer :: nq
    real(DP), allocatable :: qpts(:,:) !< (3,nq)
    integer, allocatable :: nmtx(:) !< (nq)
    integer :: nmtx_max
    integer(HID_T) :: file_id
    integer :: error, ii
    logical, allocatable :: qpt_done(:)
    real(DP) :: freqs_tmp(2,pol%nfreq)
    logical :: restart_, file_exists, file_ok
    character(len=3) :: sheader='WFN'
    PUSH_SUB(eps_hdf5_setup)

    qgrid(:) = pol%qgrid(:)
    nq = pol%nq
    nmtx_max = pol%nmtx

    SAFE_ALLOCATE(qpts, (3, nq))
    SAFE_ALLOCATE(nmtx, (nq))
    SAFE_ALLOCATE(qpt_done, (nq))

    qpts(:,1:nq) = pol%qpt(:,1:nq)
    nmtx(1:nq) = pol%nmtx_of_q(1:nq)

    restart_=.false.
    if (present(restart)) restart_ = restart

    ! FHJ: Set up file: write MF header and create groups
    write(6,'(1x,2a)') "Initializing ", trim(name)
    call setup_hdf5_mf_file(trim(name))
    call write_hdf5_header_type(trim(name), sheader, SCALARSIZE, kp, gvec, syms, crys)
    call write_hdf5_gvectors(trim(name), gvec%ng, gvec%components)
    call h5fopen_f(trim(name), H5F_ACC_RDWR_F, file_id, error)

    call hdf5_create_group(file_id, 'eps_header', error)
    call hdf5_create_group(file_id, 'eps_header/params', error)
    call hdf5_create_group(file_id, 'eps_header/qpoints', error)
    call hdf5_create_group(file_id, 'eps_header/freqs', error)
    call hdf5_create_group(file_id, 'eps_header/gspace', error)
    call hdf5_create_group(file_id, 'mats', error)
    if( pol%subspace .and. (.not. pol%use_hdf5) ) then
       call hdf5_create_group(file_id, 'eps_header/subspace', error)
    endif
    call hdf5_write_int(file_id, 'eps_header/versionnumber', VER_EPS_HDF5, error)
    call hdf5_write_int(file_id, 'eps_header/flavor', SCALARSIZE, error)
    call hdf5_write_int(file_id, 'eps_header/params/matrix_type', pol%matrix_type, error)

    call hdf5_write_logical(file_id, 'eps_header/params/has_advanced', pol%has_advanced, error)
    call hdf5_write_logical(file_id, 'eps_header/params/timeordered', pol%timeordered, error)
    call hdf5_write_int(file_id, 'eps_header/params/nmatrix', pol%nmatrix, error)
    call hdf5_write_int(file_id, 'eps_header/params/matrix_flavor', pol%matrix_flavor, error)
    call hdf5_write_int(file_id, 'eps_header/params/icutv', pol%icutv, error)
    call hdf5_write_double(file_id, 'eps_header/params/ecuts', pol%ecuts, error)
    call hdf5_write_int(file_id, 'eps_header/params/nband', pol%nband, error)
    call hdf5_write_int(file_id, 'eps_header/params/skip_nvb', pol%skip_nvb, error)
    call hdf5_write_int(file_id, 'eps_header/params/skip_ncb', pol%skip_ncb, error)
    call hdf5_write_int(file_id, 'eps_header/params/nvb', pol%nvb, error)
    call hdf5_write_int(file_id, 'eps_header/params/ncb', pol%ncb, error)
    call hdf5_write_logical(file_id, 'eps_header/params/correcthead', .false., error)    
    ! call hdf5_write_int_array(file_id, 'eps_header/params/FFTgrid', (/3/), pol%FFTgrid, error)
   
    call hdf5_write_double(file_id, 'eps_header/params/efermi', pol%efermi/ryd, error)
    call hdf5_write_int(file_id, 'eps_header/params/intraband_flag', pol%intraband_flag, error)
    call hdf5_write_double(file_id, 'eps_header/params/intraband_overlap_min', pol%intraband_overlap_min, error)
    call hdf5_write_logical(file_id, 'eps_header/params/subsample', pol%subsample, error)
    call hdf5_write_logical(file_id, 'eps_header/params/subspace', pol%subspace, error)

    qpt_done(:) = .false.
    call hdf5_write_int(file_id, 'eps_header/qpoints/nq', nq, error)
    call hdf5_write_double_array(file_id, 'eps_header/qpoints/qpts', (/3,nq/), qpts, error)
    call hdf5_write_int_array(file_id, 'eps_header/qpoints/qgrid', (/3/), qgrid, error)
    call hdf5_write_logical_array(file_id, 'eps_header/qpoints/qpt_done', (/nq/), qpt_done, error)
    call hdf5_write_int(file_id, 'eps_header/freqs/freq_dep', pol%freq_dep, error)
    call hdf5_write_int(file_id, 'eps_header/freqs/nfreq', pol%nfreq, error)
    call hdf5_write_int(file_id, 'eps_header/freqs/nfreq_imag', pol%nfreq_imag, error)

    do ii = 1, pol%nfreq
       freqs_tmp(1,ii) = pol%dFreqGrid(ii) + dble(pol%dFreqBrd(ii))
       freqs_tmp(2,ii) = IMAG(pol%dFreqBrd(ii))
    enddo
    call hdf5_write_double_array(file_id, 'eps_header/freqs/freqs', (/2, pol%nfreq/), freqs_tmp, error)

    !! G-vectors-related datasets
    call hdf5_write_int_array(file_id, 'eps_header/gspace/nmtx', (/nq/), nmtx, error)
    call hdf5_write_int(file_id, 'eps_header/gspace/nmtx_max',  nmtx_max, error)
    call hdf5_create_dset(file_id, 'eps_header/gspace/ekin', H5T_NATIVE_DOUBLE, (/gvec%ng,nq/), error)
    call hdf5_create_dset(file_id, 'eps_header/gspace/gind_eps2rho', H5T_NATIVE_INTEGER, (/gvec%ng,nq/), error)
    call hdf5_create_dset(file_id, 'eps_header/gspace/gind_rho2eps', H5T_NATIVE_INTEGER, (/gvec%ng,nq/), error)
    call hdf5_create_dset(file_id, 'mats/matrix', H5T_NATIVE_DOUBLE, (/pol%matrix_flavor, nmtx_max, nmtx_max, pol%nfreq, pol%nmatrix, nq/), error)
    call hdf5_create_dset(file_id, 'mats/matrix-diagonal', H5T_NATIVE_DOUBLE, (/pol%matrix_flavor, nmtx_max, pol%nfreq, nq/), error)
    call h5fclose_f(file_id, error)

    SAFE_DEALLOCATE(qpts)
    SAFE_DEALLOCATE(nmtx)
    SAFE_DEALLOCATE(qpt_done)
    POP_SUB(eps_hdf5_setup_2)

  end subroutine eps_hdf5_setup
  
end module epswrite_hdf5_m

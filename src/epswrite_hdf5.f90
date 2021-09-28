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
  public :: set_qpt_done, is_qpt_done, eps_hdf5_setup, eps_hdf5_setup_part
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

  subroutine eps_hdf5_setup_part(kp, gvec, syms, crys, pol, name, iq_start, iq_end, restart)
    type(kpoints), intent(in) :: kp
    type(gspace), intent(in) :: gvec
    type(symmetry), intent(in) :: syms
    type(crystal), intent(in) :: crys
    type(polarizability), intent(in) :: pol
    character(len=*), intent(in) :: name
    integer, intent(in) :: iq_start, iq_end    
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
    PUSH_SUB(eps_hdf5_setup_part)

    qgrid(:) = pol%qgrid(:)
    nq = iq_end - iq_start + 1
    if (nq <= 0) then
       call die("Invalid iq_start and iq_end.", only_root_writes=.true.)
    endif
    nmtx_max = MAXVAL(pol%nmtx_of_q(iq_start:iq_end))

    SAFE_ALLOCATE(qpts, (3, nq))
    SAFE_ALLOCATE(nmtx, (nq))
    SAFE_ALLOCATE(qpt_done, (nq))

    qpts(:,1:nq) = pol%qpt(:,iq_start:iq_end)
    nmtx(1:nq) = pol%nmtx_of_q(iq_start:iq_end)

    restart_=.false.
    if (present(restart)) restart_ = restart

    ! FHJ: try to restart the calculation, if possible and desired.
    ! We ignore the restart flags if the file doesn`t exist. However, the code
    ! dies if the file exists and is incompatible or looks corrupted.
    if (restart_) then
       call try_restart()
       if (file_ok) return
    endif

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
    ! FHJ: General datasets
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
    call hdf5_write_int_array(file_id, 'eps_header/params/FFTgrid', (/3/), pol%FFTgrid, error)       
    call hdf5_write_double(file_id, 'eps_header/params/efermi', pol%efermi/ryd, error)
    call hdf5_write_int(file_id, 'eps_header/params/intraband_flag', pol%intraband_flag, error)
    call hdf5_write_double(file_id, 'eps_header/params/intraband_overlap_min', pol%intraband_overlap_min, error)
    call hdf5_write_logical(file_id, 'eps_header/params/subsample', pol%subsample, error)
    call hdf5_write_logical(file_id, 'eps_header/params/subspace', pol%subspace, error)
    ! FHJ: Q-points-related datasets
    qpt_done(:) = .false.
    call hdf5_write_int(file_id, 'eps_header/qpoints/nq', nq, error)
    call hdf5_write_double_array(file_id, 'eps_header/qpoints/qpts', (/3,nq/), qpts, error)
    call hdf5_write_int_array(file_id, 'eps_header/qpoints/qgrid', (/3/), qgrid, error)
    call hdf5_write_logical_array(file_id, 'eps_header/qpoints/qpt_done', (/nq/), qpt_done, error)
    ! FHJ: Frequency-related datasets
    call hdf5_write_int(file_id, 'eps_header/freqs/freq_dep', pol%freq_dep, error)
    call hdf5_write_int(file_id, 'eps_header/freqs/nfreq', pol%nfreq, error)
    call hdf5_write_int(file_id, 'eps_header/freqs/nfreq_imag', pol%nfreq_imag, error)
    do ii=1,pol%nfreq !FHJ: TODO - have an unique complex freq grid in the future!
       freqs_tmp(1,ii) = pol%dFreqGrid(ii) + dble(pol%dFreqBrd(ii))
       freqs_tmp(2,ii) = IMAG(pol%dFreqBrd(ii))
    enddo
    call hdf5_write_double_array(file_id, 'eps_header/freqs/freqs', (/2, pol%nfreq/), freqs_tmp, error)

    ! FHJ: G-vectors-related datasets
    call hdf5_write_int_array(file_id, 'eps_header/gspace/nmtx', (/nq/), nmtx, error)
    call hdf5_write_int(file_id, 'eps_header/gspace/nmtx_max',  nmtx_max, error)
    call hdf5_create_dset(file_id, 'eps_header/gspace/ekin', H5T_NATIVE_DOUBLE, (/gvec%ng,nq/), error)
    call hdf5_create_dset(file_id, 'eps_header/gspace/gind_eps2rho', H5T_NATIVE_INTEGER, (/gvec%ng,nq/), error)
    call hdf5_create_dset(file_id, 'eps_header/gspace/gind_rho2eps', H5T_NATIVE_INTEGER, (/gvec%ng,nq/), error)
    call hdf5_create_dset(file_id, 'mats/matrix', H5T_NATIVE_DOUBLE, (/ pol%matrix_flavor, nmtx_max, nmtx_max, pol%nfreq, pol%nmatrix, nq /), error)
    call hdf5_create_dset(file_id, 'mats/matrix-diagonal', H5T_NATIVE_DOUBLE, (/ pol%matrix_flavor, nmtx_max, pol%nfreq, nq /), error)
    call h5fclose_f(file_id, error)
    SAFE_DEALLOCATE(qpts)
    SAFE_DEALLOCATE(nmtx)
    SAFE_DEALLOCATE(qpt_done)
    POP_SUB(eps_hdf5_setup_part)
  contains

    subroutine try_restart()
      integer :: nspin_old, nq_old
      real(DP) :: qpts_old(3,nq)
      integer :: freq_dep_old, nfreq_old
      integer :: ng_old, gvecs_old(3,gvec%ng), nmtx_old(nq), matrix_flavor_old, nmatrix_old

      write(6,'(1x,2a)') "Trying to restart file ", trim(name)
      inquire(file=trim(name), exist=file_exists)
      if (file_exists) then
         call h5fopen_f(trim(name), H5F_ACC_RDONLY_F, file_id, error)
         if (error==0) then
            ! FHJ: Consistency check.
            call h5lexists_f(file_id, 'eps_header/qpoints/qpt_done', file_ok, error)

            if (file_ok) then
               call hdf5_read_int(file_id, 'mf_header/kpoints/nspin', nspin_old, error)

               if (error==0) call hdf5_read_int(file_id, 'eps_header/qpoints/nq', nq_old, error)
               if (error==0) call hdf5_read_double_array(file_id, 'eps_header/qpoints/qpts', (/3,nq/), qpts_old, error)
               if (error==0) call hdf5_read_int(file_id, 'eps_header/freqs/freq_dep', freq_dep_old, error)
               if (error==0) call hdf5_read_int(file_id, 'eps_header/freqs/nfreq', nfreq_old, error)
               if (error==0) call hdf5_read_int(file_id, 'mf_header/gspace/ng', ng_old, error)
               if (error==0) call hdf5_read_int_array(file_id, 'mf_header/gspace/components', (/3,gvec%ng/), gvecs_old, error)
               if (error==0) call hdf5_read_int_array(file_id, 'eps_header/gspace/nmtx', (/nq/), nmtx_old, error)
               if (error==0) call hdf5_read_int(file_id, 'eps_header/params/matrix_flavor', matrix_flavor_old, error)
               if (error==0) call hdf5_read_int(file_id, 'eps_header/params/nmatrix', nmatrix_old, error)
               file_ok = (error==0) .and. (nspin_old==kp%nspin .and. nq_old==nq .and. &
                    all(dabs(qpts_old-qpts)<TOL_ZERO) .and. freq_dep_old==pol%freq_dep .and. &
                    nfreq_old==pol%nfreq .and. ng_old==gvec%ng .and. all(gvecs_old==gvec%components) .and. &
                    all(nmtx_old==nmtx) .and. matrix_flavor_old==pol%matrix_flavor .and. &
                    nmatrix_old==pol%nmatrix)
            endif
            call h5fclose_f(file_id, error)

            if (file_ok) then
               ! FHJ: File *seems* alright, we don`t have to initialize it
               write(6,'(1x,2a)') "Everything looks ok: restarting file ", trim(name)
               return
            endif
            write(0,*)
            write(0,'(3a)') "ERROR: file ", trim(name), " is incompatible with the current calculation."
            write(0,*) 'Values from file vs. calculation:'
            write(0,*) 'nspin:', nspin_old, kp%nspin
            write(0,*) 'nq:', nq_old, nq
            write(0,*) 'qpts (same?):', all(dabs(qpts_old-qpts)<TOL_ZERO)
            write(0,*) 'freq_dep:', freq_dep_old, pol%freq_dep
            write(0,*) 'nfreq:', nfreq_old, pol%nfreq

            write(0,*) 'ng:', ng_old, gvec%ng
            write(0,*) 'gvecs (same?):', all(gvecs_old==gvec%components)
            write(0,*) 'nmtx (same?):', all(nmtx_old==nmtx)
            write(0,*) 'matrix_flavor:', matrix_flavor_old, pol%matrix_flavor
            write(0,*) 'nmatrix:', nmatrix_old, pol%nmatrix
            write(0,*)
            write(0,*) 'NOTE: you should only trust the first pair of values that disagree.'
            write(0,*)
            call die("file "//trim(name)//" is incompatible with the current calculation.", only_root_writes=.true.)
         else
            write(0,'(3a,i0,a)') "ERROR: ", trim(name), " is not a valid HDF5 file (error code: ", error," )."
            write(0,'(a)') 'Make sure the file is not corrupted'
            call die("file "//trim(name)//" looks corrupted", only_root_writes=.true.)
         endif
      endif
      file_ok = .false.
      write(0,'(3a)') "WARNING: file ", trim(name), " doesn`t exist. We`ll start the calculation from scratch."
      restart_ = .false.
      if (present(restart)) restart = .false.
    end subroutine try_restart
  end subroutine eps_hdf5_setup_part
  
end module epswrite_hdf5_m

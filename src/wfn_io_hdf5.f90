#include "f_defs.h"

!!===================================================================
!!
!! Module:
!!
!! (1) wfn_io_hdf5_m Originally by JIM Last Modified 4/25/2012 (JIM)
!!
!! Routines to read and write wavefunctions in 1 format.
!! The code is generated through repeated inclusion of a file with
!! different preprocessor definitions each time. You are not expected to
!! understand this. Consult the resulting .p.f file for clarity.
!!
!!=================================================================
module wfn_io_hdf5_m
  use global_m
  use hdf5_io_m
  use hdf5
  implicit none
  private
  public :: &
       read_hdf5_header_type , write_hdf5_header_type , &
       read_hdf5_gvectors , write_hdf5_gvectors , &
       setup_hdf5_mf_file , &
       setup_hdf5_wfn_file , &
       read_hdf5_mf_header , &
       write_hdf5_mf_header , &
       read_hdf5_gvec_kpt, read_info, read_kpoints, read_gspace, read_symmetry, read_crystal
contains

  !=========================================================================
  ! begin read/write header
  subroutine read_hdf5_header_type(sFileName, sheader, iflavor, kp, gvec, syms, crys)
    character(len=*), intent(in) :: sFileName
    character(len=3), intent(inout) :: sheader
    integer, intent(inout) :: iflavor
    type(kpoints), intent(out) :: kp
    type(gspace), intent(out) :: gvec
    type(symmetry), intent(out) :: syms
    type(crystal), intent(out) :: crys
    ! set values based on epsilon calculation
    logical :: is_get=.false.
    logical :: wfnflag=.true.

    if (peinf%inode == 0) then
       call read_info(TRUNC(sFileName),iflavor)
       call read_kpoints(TRUNC(sFileName),kp)
       call read_gspace(TRUNC(sFileName),gvec)
       call read_symmetry(TRUNC(sFileName),syms)
       call read_crystal(TRUNC(sFileName),crys)
    endif
    gvec%nFFTgridpts = product(gvec%FFTgrid(1:3))
    if (peinf%npes > 1) then
       call MPI_BCAST(is_get, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, mpierr)
       call MPI_BCAST(kp%nspin, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
       call MPI_BCAST(kp%nspinor, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
       call MPI_BCAST(gvec%ng, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
       call MPI_BCAST(syms%ntran, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
       call MPI_BCAST(syms%cell_symmetry, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
       call MPI_BCAST(crys%nat, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
       call MPI_BCAST(gvec%ecutrho, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
       call MPI_BCAST(gvec%FFTgrid, 3, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
       call MPI_BCAST(gvec%nFFTgridpts, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
       call MPI_BCAST(crys%celvol, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
       call MPI_BCAST(crys%alat, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
       call MPI_BCAST(crys%avec, 9, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
       call MPI_BCAST(crys%adot, 9, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
       call MPI_BCAST(crys%recvol, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
       call MPI_BCAST(crys%blat, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
       call MPI_BCAST(crys%bvec, 9, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
       call MPI_BCAST(crys%bdot, 9, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
       call MPI_BCAST(syms%mtrx(1,1,1), 3*3*48, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
       ! call MPI_BCAST(syms%mtrx_reci(1,1,1), 3*3*48, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
       ! call MPI_BCAST(syms%mtrx_cart(1,1,1), 3*3*48, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
       call MPI_BCAST(syms%tnp, 3*48, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
       if (wfnflag) then
          call MPI_BCAST(kp%nrk, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
          call MPI_BCAST(kp%mnband, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
          call MPI_BCAST(kp%ngkmax, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
          call MPI_BCAST(kp%ecutwfc, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
          call MPI_BCAST(kp%kgrid, 3, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
          call MPI_BCAST(kp%shift, 3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
       endif
    endif
    if(peinf%inode > 0 ) then
       allocate(crys%atyp (crys%nat))
       allocate(crys%apos (3, crys%nat))
       if (wfnflag) then
          allocate(kp%ngk (kp%nrk))
          allocate(kp%w (kp%nrk))
          allocate(kp%rk (3, kp%nrk))
          allocate(kp%ifmin (kp%nrk, kp%nspin))
          allocate(kp%ifmax (kp%nrk, kp%nspin))
          allocate(kp%el (kp%mnband, kp%nrk, kp%nspin))
          allocate(kp%occ (kp%mnband, kp%nrk, kp%nspin))
       endif
    endif
    if (peinf%npes > 1) then
       call MPI_BCAST(crys%atyp, crys%nat, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
       call MPI_BCAST(crys%apos, 3*crys%nat, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
       if (wfnflag) then
          call MPI_BCAST(kp%ngk, kp%nrk, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
          call MPI_BCAST(kp%w, kp%nrk, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
          call MPI_BCAST(kp%rk(1,1), 3*kp%nrk, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
          call MPI_BCAST(kp%ifmin(1,1), size(kp%ifmin), MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
          call MPI_BCAST(kp%ifmax(1,1), size(kp%ifmax), MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
          call MPI_BCAST(kp%el(1,1,1), size(kp%el), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
          call MPI_BCAST(kp%occ(1,1,1), size(kp%occ), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
       endif
    endif

  end subroutine read_hdf5_header_type
  subroutine read_info(sFileName, iflavor)
    character(len=*), intent(in) :: sFileName
    integer, intent(out) :: iflavor
    integer(HID_T) :: hidFile
    integer :: iError
    logical :: exists

    write(*,*) "sFileName = ", sFileName
    call h5fopen_f(sFileName, H5F_ACC_RDONLY_F, hidFile, iError)
    !FHJ: Keep backwards compatibility with older WFN.h5 files
    call h5lexists_f(hidFile, '/mf_header/flavor', exists, iError)
    if (exists .and. iError==0) then
       call hdf5_read_int(hidFile, '/mf_header/flavor', iflavor, iError)
    else
       call hdf5_read_int(hidFile, '/info', iflavor, iError)
    endif
    call h5fclose_f(hidFile, iError)

  end subroutine read_info
  subroutine read_gspace(sFileName,gvec)
    character(len=*), intent(in) :: sFileName
    type(gspace), intent(out) :: gvec
    integer(HID_T) :: hidFile
    integer :: iError

    call h5fopen_f(sFileName, H5F_ACC_RDONLY_F, hidFile, iError)
    call hdf5_read_int(hidFile, '/mf_header/gspace/ng', gvec%ng, iError)
    call hdf5_read_double(hidFile, '/mf_header/gspace/ecutrho', gvec%ecutrho, iError)
    call hdf5_read_int_array(hidFile, '/mf_header/gspace/FFTgrid', (/3/), gvec%FFTgrid, iError)
    call h5fclose_f(hidFile, iError)

  end subroutine read_gspace
  subroutine read_hdf5_gvectors(sFileName, ng, gvec)
    character(len=*), intent(in) :: sFileName
    integer, intent(in) :: ng !< used size of array
    integer, intent(out) :: gvec(:, :) !< (3, ng_bound)
    integer(HID_T) :: hidFile
    integer :: iError
    logical :: bcast_, dont_read_

    dont_read_=.false.
    bcast_=.not. dont_read_
    if(peinf%inode == 0) then
       call h5fopen_f(sFileName, H5F_ACC_RDONLY_F, hidFile, iError)
       call hdf5_read_int_array(hidFile, '/mf_header/gspace/components', (/3,ng/), gvec, iError)
       call h5fclose_f(hidFile, iError)
    endif
    if(peinf%npes > 1) then
       if(bcast_) then
          call MPI_BCAST(gvec(1,1), 3 * ng, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
       endif
    endif

  end subroutine read_hdf5_gvectors
  subroutine read_symmetry(sFileName,syms)
    character(len=*), intent(in) :: sFileName
    type(symmetry), intent(out) :: syms
    integer(HID_T) :: hidFile
    integer :: iError

    call h5fopen_f(sFileName, H5F_ACC_RDONLY_F, hidFile, iError)
    call hdf5_read_int(hidFile, '/mf_header/symmetry/ntran', syms%ntran, iError)
    call hdf5_read_int(hidFile, '/mf_header/symmetry/cell_symmetry', syms%cell_symmetry, iError)
    call hdf5_read_int_array(hidFile, '/mf_header/symmetry/mtrx', (/3, 3, 48/), syms%mtrx, iError)
    ! call hdf5_read_int_array(hidFile, '/mf_header/symmetry/mtrx_reci', (/3, 3, 48/), syms%mtrx_reci, iError)
    ! call hdf5_read_double_array(hidFile, '/mf_header/symmetry/mtrx_cart', (/3, 3, 48/), syms%mtrx_cart, iError)
    call hdf5_read_double_array(hidFile, '/mf_header/symmetry/tnp', (/3, 48/), syms%tnp, iError)
    call h5fclose_f(hidFile, iError)

  end subroutine read_symmetry
  subroutine read_crystal(sFileName,crys)
    character(len=*), intent(in) :: sFileName
    type(crystal), intent(out) :: crys
    integer(HID_T) :: hidFile
    integer :: iError

    call h5fopen_f(sFileName, H5F_ACC_RDONLY_F, hidFile, iError)
    call hdf5_read_double(hidFile, '/mf_header/crystal/celvol', crys%celvol, iError)
    call hdf5_read_double(hidFile, '/mf_header/crystal/recvol', crys%recvol, iError)
    call hdf5_read_double(hidFile, '/mf_header/crystal/alat', crys%alat, iError)
    call hdf5_read_double(hidFile, '/mf_header/crystal/blat', crys%blat, iError)
    call hdf5_read_int(hidFile, '/mf_header/crystal/nat', crys%nat, iError)
    call hdf5_read_double_array(hidFile, '/mf_header/crystal/avec', (/3, 3/), crys%avec, iError)
    call hdf5_read_double_array(hidFile, '/mf_header/crystal/bvec', (/3, 3/), crys%bvec, iError)
    call hdf5_read_double_array(hidFile, '/mf_header/crystal/adot', (/3, 3/), crys%adot, iError)
    call hdf5_read_double_array(hidFile, '/mf_header/crystal/bdot', (/3, 3/), crys%bdot, iError)
    allocate(crys%atyp (crys%nat))
    allocate(crys%apos (3,crys%nat))
    call hdf5_read_int_array(hidFile, '/mf_header/crystal/atyp', (/crys%nat/), crys%atyp, iError)
    call hdf5_read_double_array(hidFile, '/mf_header/crystal/apos', (/3, crys%nat/), crys%apos, iError)
    call h5fclose_f(hidFile, iError)

  end subroutine read_crystal
  subroutine read_kpoints(sFileName,kp)
    character(len=*), intent(in) :: sFileName
    type(kpoints), intent(out) :: kp
    integer(HID_T) :: hidFile
    integer :: iError

    call h5fopen_f(sFileName, H5F_ACC_RDONLY_F, hidFile, iError)
    call hdf5_read_int(hidFile, '/mf_header/kpoints/nspin', kp%nspin, iError)
    call hdf5_read_int(hidFile, '/mf_header/kpoints/nspinor', kp%nspinor, iError)
    call hdf5_read_int(hidFile, '/mf_header/kpoints/nrk', kp%nrk, iError)
    call hdf5_read_int(hidFile, '/mf_header/kpoints/mnband', kp%mnband, iError)
    call hdf5_read_int(hidFile, '/mf_header/kpoints/ngkmax', kp%ngkmax, iError)
    call hdf5_read_double(hidFile, '/mf_header/kpoints/ecutwfc', kp%ecutwfc, iError)
    call hdf5_read_int_array(hidFile, '/mf_header/kpoints/kgrid', (/3/), kp%kgrid, iError)
    call hdf5_read_double_array(hidFile, '/mf_header/kpoints/shift', (/3/), kp%shift, iError)
    allocate(kp%ngk (kp%nrk))
    allocate(kp%ifmin (kp%nrk, kp%nspin))
    allocate(kp%ifmax (kp%nrk, kp%nspin))
    allocate(kp%w (kp%nrk))
    allocate(kp%rk (3,kp%nrk))
    allocate(kp%el (kp%mnband,kp%nrk,kp%nspin))
    allocate(kp%occ (kp%mnband,kp%nrk,kp%nspin))
    call hdf5_read_int_array(hidFile, '/mf_header/kpoints/ngk', &
         (/kp%nrk/), kp%ngk, iError)
    call hdf5_read_int_array(hidFile, '/mf_header/kpoints/ifmin', &
         (/kp%nrk, kp%nspin/), kp%ifmin, iError)
    call hdf5_read_int_array(hidFile, '/mf_header/kpoints/ifmax', &
         (/kp%nrk, kp%nspin/), kp%ifmax, iError)
    call hdf5_read_double_array(hidFile, '/mf_header/kpoints/w', &
         (/kp%nrk/), kp%w, iError)
    call hdf5_read_double_array(hidFile, '/mf_header/kpoints/rk', &
         (/3, kp%nrk/), kp%rk, iError)
    call hdf5_read_double_array(hidFile, '/mf_header/kpoints/el', &
         (/kp%mnband, kp%nrk, kp%nspin/), kp%el, iError)
    call hdf5_read_double_array(hidFile, '/mf_header/kpoints/occ', &
         (/kp%mnband, kp%nrk, kp%nspin/), kp%occ, iError)
    call h5fclose_f(hidFile, iError)

  end subroutine read_kpoints
  ! end read/write header

  ! begin read/write wfn gvectors
  subroutine read_hdf5_wfn_gvectors(fname, gvec, ngktot)
    character(len=*) :: fname
    integer, intent(inout) :: gvec(:,:)
    integer, intent(in) :: ngktot
    integer(HID_T) :: file_id
    integer :: error

    if(peinf%inode == 0) then
       call h5fopen_f(fname, H5F_ACC_RDONLY_F, file_id, error)
       call hdf5_read_int_array(file_id, 'wfns/gvecs', (/3, ngktot/), gvec, error)
       call h5fclose_f(file_id, error)
    endif
    if(peinf%npes > 1) then
       call MPI_Bcast(gvec(1,1), 3*ngktot, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
    endif

  end subroutine read_hdf5_wfn_gvectors
  ! end read/write wfn gvectors

  ! begin read/write wfn data
  subroutine read_hdf5_band_complex(fname, wfn, ngk, nstot, ioffsetk, ioffsetb)
    character(len=*) :: fname
    complex(DPC), intent(out) :: wfn(:,:) !< (ngk,nstot=kp%ns*kp%nspinor)
    integer, intent(in) :: ngk
    integer, intent(in) :: nstot
    integer, intent(in) :: ioffsetk
    integer, intent(in) :: ioffsetb
    real(DP) :: dwfn(2,ngk,nstot,1)
    integer(HID_T) :: file_id
    integer(HID_T) :: dataset_id
    integer(HID_T) :: dataspace_id
    integer(HID_T) :: memspace_id
    integer(HSIZE_T) :: a3(4), offset(4), count(4)
    integer :: error

    a3(1) = 2
    a3(2) = ngk
    a3(3) = nstot
    a3(4) = 1
    offset(1) = 0
    offset(2) = ioffsetk
    offset(3) = 0
    offset(4) = ioffsetb
    count(1) = 2
    count(2) = ngk
    count(3) = nstot
    count(4) = 1
    if(peinf%inode == 0) then
       call h5fopen_f(fname, H5F_ACC_RDONLY_F, file_id, error)
       call h5dopen_f(file_id, 'wfns/coeffs', dataset_id, error)
       CALL h5screate_simple_f(4, count, memspace_id, error)
       call h5dget_space_f(dataset_id, dataspace_id, error)
       call h5sselect_hyperslab_f(dataspace_id, H5S_SELECT_SET_F, offset, count, error)
       call h5dread_f(dataset_id, H5T_NATIVE_DOUBLE, dwfn, a3, error, file_space_id = dataspace_id, mem_space_id = memspace_id)
       call h5dclose_f(dataset_id, error)
       call h5sclose_f(dataspace_id, error)
       call h5sclose_f(memspace_id, error)
       call h5fclose_f(file_id, error)
    endif
    wfn(:,:) = DCMPLX(dwfn(1,:,:,1), dwfn(2,:,:,1))

  end subroutine read_hdf5_band_complex
  ! end read/write wfn data

  ! begin read/write wfn data
  subroutine read_hdf5_band_real(fname, wfn, ngk, nstot, ioffsetk, ioffsetb)
    character(len=*) :: fname
    real(DPC), intent(out) :: wfn(:,:) !< (ngk,nstot=kp%ns*kp%nspinor)
    integer, intent(in) :: ngk
    integer, intent(in) :: nstot
    integer, intent(in) :: ioffsetk
    integer, intent(in) :: ioffsetb
    real(DP) :: dwfn(1,ngk,nstot,1)
    integer(HID_T) :: file_id
    integer(HID_T) :: dataset_id
    integer(HID_T) :: dataspace_id
    integer(HID_T) :: memspace_id
    integer(HSIZE_T) :: a3(4), offset(4), count(4)
    integer :: error

    a3(1) = 1
    a3(2) = ngk
    a3(3) = nstot
    a3(4) = 1
    offset(1) = 0
    offset(2) = ioffsetk
    offset(3) = 0
    offset(4) = ioffsetb
    count(1) = 1
    count(2) = ngk
    count(3) = nstot
    count(4) = 1
    if(peinf%inode == 0) then
       call h5fopen_f(fname, H5F_ACC_RDONLY_F, file_id, error)
       call h5dopen_f(file_id, 'wfns/coeffs', dataset_id, error)
       CALL h5screate_simple_f(4, count, memspace_id, error)
       call h5dget_space_f(dataset_id, dataspace_id, error)
       call h5sselect_hyperslab_f(dataspace_id, H5S_SELECT_SET_F, offset, count, error)
       call h5dread_f(dataset_id, H5T_NATIVE_DOUBLE, dwfn, a3, error, file_space_id = dataspace_id, mem_space_id = memspace_id)
       call h5dclose_f(dataset_id, error)
       call h5sclose_f(dataspace_id, error)
       call h5sclose_f(memspace_id, error)
       call h5fclose_f(file_id, error)
    endif
    wfn(:,:) = dwfn(1,:,:,1)

  end subroutine read_hdf5_band_real
  ! end read/write wfn data

  ! begin read/write header
  subroutine write_hdf5_header_type(sFileName, sheader, iflavor, kp, gvec, syms, crys)
    character(len=*), intent(in) :: sFileName
    character(len=3), intent(inout) :: sheader
    integer, intent(in) :: iflavor
    type(kpoints), intent(in) :: kp
    type(gspace), intent(in) :: gvec
    type(symmetry), intent(in) :: syms
    type(crystal), intent(in) :: crys
    ! set values based on epsilon calculation
    logical :: is_get=.false.
    logical :: wfnflag=.true.

    if (peinf%inode == 0) then
       call write_info(TRUNC(sFileName),iflavor)
       call write_kpoints(TRUNC(sFileName),kp)
       call write_gspace(TRUNC(sFileName),gvec)
       call write_symmetry(TRUNC(sFileName),syms)
       call write_crystal(TRUNC(sFileName),crys)
    endif

  end subroutine write_hdf5_header_type
  subroutine write_info(sFileName, iflavor)
    character(len=*), intent(in) :: sFileName
    integer, intent(in) :: iflavor
    integer(HID_T) :: hidFile
    integer :: iError

    write(*,*) "sFileName = ", sFileName
    call h5fopen_f(sFileName, H5F_ACC_RDWR_F, hidFile, iError)
    call hdf5_write_int(hidFile, '/mf_header/versionnumber', VER_WFN_HDF5, iError)
    call hdf5_write_int(hidFile, '/mf_header/flavor', iflavor, iError)
    call h5fclose_f(hidFile, iError)

  end subroutine write_info
  subroutine write_gspace(sFileName,gvec)
    character(len=*), intent(in) :: sFileName
    type(gspace), intent(in) :: gvec
    integer(HID_T) :: hidFile
    integer :: iError

    call h5fopen_f(sFileName, H5F_ACC_RDWR_F, hidFile, iError)
    call hdf5_write_int(hidFile, '/mf_header/gspace/ng', gvec%ng, iError)
    call hdf5_write_double(hidFile, '/mf_header/gspace/ecutrho', gvec%ecutrho, iError)
    call hdf5_write_int_array(hidFile, '/mf_header/gspace/FFTgrid', (/3/), gvec%FFTgrid, iError)
    call h5fclose_f(hidFile, iError)

  end subroutine write_gspace
  subroutine write_hdf5_gvectors(sFileName, ng, gvec)
    character(len=*), intent(in) :: sFileName
    integer, intent(in) :: ng !< used size of array
    integer, intent(in) :: gvec(:, :) !< (3, ng_bound)
    integer(HID_T) :: hidFile
    integer :: iError
    logical :: bcast_, dont_read_

    dont_read_=.false.
    bcast_=.not. dont_read_
    if(peinf%inode == 0) then
       call h5fopen_f(sFileName, H5F_ACC_RDWR_F, hidFile, iError)
       call hdf5_write_int_array(hidFile, '/mf_header/gspace/components', (/3,ng/), gvec, iError)
       call h5fclose_f(hidFile, iError)
    endif

  end subroutine write_hdf5_gvectors
  subroutine write_symmetry(sFileName,syms)
    character(len=*), intent(in) :: sFileName
    type(symmetry), intent(in) :: syms
    integer(HID_T) :: hidFile
    integer :: iError

    call h5fopen_f(sFileName, H5F_ACC_RDWR_F, hidFile, iError)
    call hdf5_write_int(hidFile, '/mf_header/symmetry/ntran', syms%ntran, iError)
    call hdf5_write_int(hidFile, '/mf_header/symmetry/cell_symmetry', syms%cell_symmetry, iError)
    call hdf5_write_int_array(hidFile, '/mf_header/symmetry/mtrx', (/3, 3, 48/), syms%mtrx, iError)
    ! call hdf5_write_int_array(hidFile, '/mf_header/symmetry/mtrx_reci', (/3, 3, 48/), syms%mtrx_reci, iError)
    ! call hdf5_write_double_array(hidFile, '/mf_header/symmetry/mtrx_cart', (/3, 3, 48/), syms%mtrx_cart, iError)
    call hdf5_write_double_array(hidFile, '/mf_header/symmetry/tnp', (/3, 48/), syms%tnp, iError)
    call h5fclose_f(hidFile, iError)

  end subroutine write_symmetry
  subroutine write_crystal(sFileName,crys)
    character(len=*), intent(in) :: sFileName
    type(crystal), intent(in) :: crys
    integer(HID_T) :: hidFile
    integer :: iError

    call h5fopen_f(sFileName, H5F_ACC_RDWR_F, hidFile, iError)
    call hdf5_write_double(hidFile, '/mf_header/crystal/celvol', crys%celvol, iError)
    call hdf5_write_double(hidFile, '/mf_header/crystal/recvol', crys%recvol, iError)
    call hdf5_write_double(hidFile, '/mf_header/crystal/alat', crys%alat, iError)
    call hdf5_write_double(hidFile, '/mf_header/crystal/blat', crys%blat, iError)
    call hdf5_write_int(hidFile, '/mf_header/crystal/nat', crys%nat, iError)
    call hdf5_write_double_array(hidFile, '/mf_header/crystal/avec', (/3, 3/), crys%avec, iError)
    call hdf5_write_double_array(hidFile, '/mf_header/crystal/bvec', (/3, 3/), crys%bvec, iError)
    call hdf5_write_double_array(hidFile, '/mf_header/crystal/adot', (/3, 3/), crys%adot, iError)
    call hdf5_write_double_array(hidFile, '/mf_header/crystal/bdot', (/3, 3/), crys%bdot, iError)
    call hdf5_write_int_array(hidFile, '/mf_header/crystal/atyp', (/crys%nat/), crys%atyp, iError)
    call hdf5_write_double_array(hidFile, '/mf_header/crystal/apos', (/3, crys%nat/), crys%apos, iError)
    call h5fclose_f(hidFile, iError)

  end subroutine write_crystal
  subroutine write_kpoints(sFileName,kp)
    character(len=*), intent(in) :: sFileName
    type(kpoints), intent(in) :: kp
    integer(HID_T) :: hidFile
    integer :: iError

    call h5fopen_f(sFileName, H5F_ACC_RDWR_F, hidFile, iError)
    call hdf5_write_int(hidFile, '/mf_header/kpoints/nspin', kp%nspin, iError)
    call hdf5_write_int(hidFile, '/mf_header/kpoints/nspinor', kp%nspinor, iError)
    call hdf5_write_int(hidFile, '/mf_header/kpoints/nrk', kp%nrk, iError)
    call hdf5_write_int(hidFile, '/mf_header/kpoints/mnband', kp%mnband, iError)
    call hdf5_write_int(hidFile, '/mf_header/kpoints/ngkmax', kp%ngkmax, iError)
    call hdf5_write_double(hidFile, '/mf_header/kpoints/ecutwfc', kp%ecutwfc, iError)
    call hdf5_write_int_array(hidFile, '/mf_header/kpoints/kgrid', (/3/), kp%kgrid, iError)
    call hdf5_write_double_array(hidFile, '/mf_header/kpoints/shift', (/3/), kp%shift, iError)
    call hdf5_write_int_array(hidFile, '/mf_header/kpoints/ngk', &
         (/kp%nrk/), kp%ngk, iError)
    call hdf5_write_int_array(hidFile, '/mf_header/kpoints/ifmin', &
         (/kp%nrk, kp%nspin/), kp%ifmin, iError)
    call hdf5_write_int_array(hidFile, '/mf_header/kpoints/ifmax', &
         (/kp%nrk, kp%nspin/), kp%ifmax, iError)
    call hdf5_write_double_array(hidFile, '/mf_header/kpoints/w', &
         (/kp%nrk/), kp%w, iError)
    call hdf5_write_double_array(hidFile, '/mf_header/kpoints/rk', &
         (/3, kp%nrk/), kp%rk, iError)
    call hdf5_write_double_array(hidFile, '/mf_header/kpoints/el', &
         (/kp%mnband, kp%nrk, kp%nspin/), kp%el, iError)
    call hdf5_write_double_array(hidFile, '/mf_header/kpoints/occ', &
         (/kp%mnband, kp%nrk, kp%nspin/), kp%occ, iError)
    call h5fclose_f(hidFile, iError)

  end subroutine write_kpoints
  ! end read/write header

  ! begin read/write wfn gvectors
  subroutine write_hdf5_wfn_gvectors(fname, gvec, ngktot)
    character(len=*) :: fname
    integer, intent(inout) :: gvec(:,:)
    integer, intent(in) :: ngktot
    integer(HID_T) :: file_id
    integer :: error

    if(peinf%inode == 0) then
       call h5fopen_f(fname, H5F_ACC_RDWR_F, file_id, error)
       call hdf5_write_int_array(file_id, 'wfns/gvecs', (/3, ngktot/), gvec, error)
       call h5fclose_f(file_id, error)
    endif

  end subroutine write_hdf5_wfn_gvectors
  ! end read/write wfn gvectors

  ! begin read/write wfn data
  subroutine write_hdf5_band_complex(fname, wfn, ngk, nstot, ioffsetk, ioffsetb)
    character(len=*) :: fname
    complex(DPC), intent(in) :: wfn(:,:) !< (ngk,nstot=kp%ns*kp%nspinor)
    integer, intent(in) :: ngk
    integer, intent(in) :: nstot
    integer, intent(in) :: ioffsetk
    integer, intent(in) :: ioffsetb
    real(DP) :: dwfn(2,ngk,nstot,1)
    integer(HID_T) :: file_id
    integer(HID_T) :: dataset_id
    integer(HID_T) :: dataspace_id
    integer(HID_T) :: memspace_id
    integer(HSIZE_T) :: a3(4), offset(4), count(4)
    integer :: error

    a3(1) = 2
    a3(2) = ngk
    a3(3) = nstot
    a3(4) = 1
    offset(1) = 0
    offset(2) = ioffsetk
    offset(3) = 0
    offset(4) = ioffsetb
    count(1) = 2
    count(2) = ngk
    count(3) = nstot
    count(4) = 1
    dwfn(1,:,:,1) = real(wfn(:,:))
    dwfn(2,:,:,1) = AIMAG(wfn(:,:))
    if(peinf%inode == 0) then
       call h5fopen_f(fname, H5F_ACC_RDWR_F, file_id, error)
       call h5dopen_f(file_id, 'wfns/coeffs', dataset_id, error)
       CALL h5screate_simple_f(4, count, memspace_id, error)
       call h5dget_space_f(dataset_id, dataspace_id, error)
       call h5sselect_hyperslab_f(dataspace_id, H5S_SELECT_SET_F, offset, count, error)
       call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, dwfn, a3, error, file_space_id = dataspace_id, mem_space_id = memspace_id)
       call h5dclose_f(dataset_id, error)
       call h5sclose_f(dataspace_id, error)
       call h5sclose_f(memspace_id, error)
       call h5fclose_f(file_id, error)
    endif

  end subroutine write_hdf5_band_complex
  ! end read/write wfn data

  ! begin read/write wfn data
  subroutine write_hdf5_band_real(fname, wfn, ngk, nstot, ioffsetk, ioffsetb)
    character(len=*) :: fname
    real(DPC), intent(in) :: wfn(:,:) !< (ngk,nstot=kp%ns*kp%nspinor)
    integer, intent(in) :: ngk
    integer, intent(in) :: nstot
    integer, intent(in) :: ioffsetk
    integer, intent(in) :: ioffsetb
    real(DP) :: dwfn(1,ngk,nstot,1)
    integer(HID_T) :: file_id
    integer(HID_T) :: dataset_id
    integer(HID_T) :: dataspace_id
    integer(HID_T) :: memspace_id
    integer(HSIZE_T) :: a3(4), offset(4), count(4)
    integer :: error

    a3(1) = 1
    a3(2) = ngk
    a3(3) = nstot
    a3(4) = 1
    offset(1) = 0
    offset(2) = ioffsetk
    offset(3) = 0
    offset(4) = ioffsetb
    count(1) = 1
    count(2) = ngk
    count(3) = nstot
    count(4) = 1
    dwfn(1,:,:,1) = wfn(:,:)
    if(peinf%inode == 0) then
       call h5fopen_f(fname, H5F_ACC_RDWR_F, file_id, error)
       call h5dopen_f(file_id, 'wfns/coeffs', dataset_id, error)
       CALL h5screate_simple_f(4, count, memspace_id, error)
       call h5dget_space_f(dataset_id, dataspace_id, error)
       call h5sselect_hyperslab_f(dataspace_id, H5S_SELECT_SET_F, offset, count, error)
       call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, dwfn, a3, error, file_space_id = dataspace_id, mem_space_id = memspace_id)
       call h5dclose_f(dataset_id, error)
       call h5sclose_f(dataspace_id, error)
       call h5sclose_f(memspace_id, error)
       call h5fclose_f(file_id, error)
    endif

  end subroutine write_hdf5_band_real
  ! end read/write wfn data

  ! end read/write wfn data
  subroutine read_hdf5_bands_block(file_id, kp, nbownmax, nbownactual, does_it_ownb, ib_first, wfnsout, ioffset)
    integer(HID_T), intent(in) :: file_id
    type(kpoints), intent(in) :: kp
    integer, intent(in) :: nbownmax
    integer, intent(in) :: nbownactual !< how many bands I own
    logical, intent(in) :: does_it_ownb(:,:)
    integer, intent(in) :: ib_first !< first band of the set of bands I own
    COMPLEX(DPC), intent(out) :: wfnsout(:,:,:) !< (SUM(kp%ngk), kp%nspin*kp%nspinor, nbownactual)
    integer, optional, intent(in) :: ioffset
    real(DP), allocatable :: wfndata(:,:,:,:)
    integer(HID_T) :: plist_id
    integer(HID_T) :: dset_id
    integer(HID_T) :: dataspace
    integer(HID_T) :: memspace
    integer(HSIZE_T) :: count(4), offset(4)
    integer :: error
    integer :: ipe, reader
    integer :: nread
    integer :: ngktot
    integer :: ioffset_
    integer, allocatable :: ranks(:)
    integer :: icount, is
    integer :: mpiworldgrp, mpigrp, bandcomm, ib, ib_, max_bands_read, bands_read, bands_read_max
    ! integer, parameter :: max_bytes_read = 536870912 ! don`t read/send more than 1/2 GB at a time
    integer(kind=8), parameter :: max_bytes_read = 1073741824 ! don`t read/send more than 1 GB at a time
    ! integer(kind=8), parameter :: max_bytes_read = 2147483648 ! don`t read/send more than 2 GB at a time
    ! integer, parameter :: max_bytes_read = 4294967296 ! don`t read/send more than 4 GB at a time
    ! integer, parameter :: max_bytes_read = 8589934592 ! don`t read/send more than 8 GB at a time
    ! 0=native 1; 1=manual group comm; 2=manual send/recvs
    integer, parameter :: comm_style=0
    logical :: do_read

    call logit('Reading HDF5 WFNs by blocks')
    allocate(ranks (peinf%npes))
    ioffset_ = 0
    if(present(ioffset)) then
       ioffset_ = ioffset
    endif
    ngktot = SUM(kp%ngk)
    nread = peinf%npools
    call h5dopen_f(file_id, 'wfns/coeffs', dset_id, error)
    call h5dget_space_f(dset_id, dataspace, error)
    ! find lowest rank PE that owns the first band of this group of bands
    ! for example, for valence bands, all the procs in a pool owns the same range of vbands
    ! in this case reader should be the lowest global rank of all these procs.
    ! And only the reader-th proc will perform the read operation and broadcast the wavefunctions to the remaining procs in this pool.
    reader = -1
    if (nbownactual>0) then
       do ipe = 1, peinf%npes
          if(does_it_ownb(ib_first,ipe)) then
             reader = ipe - 1
             exit ! when we find one, we exit, so this is the lowest rank !!
          endif
       enddo
       if (reader==-1) call die("Cannot find reader in read_hdf5_bands_block", only_root_writes=.true.)
       if (comm_style==1) then
          ! We use MPI_Bcast with 1 groups
          icount = 0
          do ipe = 1, peinf%npes
             if(does_it_ownb(ib_first,ipe)) then
                icount = icount + 1
                ranks(icount) = ipe - 1
             endif
          enddo
          call MPI_Comm_Group(MPI_COMM_WORLD, mpiworldgrp, mpierr)
          call MPI_Group_Incl(mpiworldgrp, icount, ranks, mpigrp, mpierr)
          call MPI_Comm_Create(MPI_COMM_WORLD, mpigrp, bandcomm, mpierr)
       endif
    else
       if (comm_style==1) then
          ! FHJ: Note that MPI_Comm_Create must be called by everyone in MPI_COMM_WORLD!
          call MPI_Comm_Create(MPI_COMM_WORLD, MPI_GROUP_EMPTY, bandcomm, mpierr)
       endif
    endif
    ! FHJ: read at most max_bands_read bands to avoid 1/1 buffer overflow.
    ! Here, 2*kp%nspin*kp%nspinor*dble(ngktot)*8d0 is the size of a
    ! single band, including G-vectors from all k-points.

    ! max_bands_read < nbownmax = peinf%nvownmax
    max_bands_read = min(nbownmax, &
         int(max_bytes_read/(2*kp%nspin*kp%nspinor*dble(ngktot)*8d0)))
    if (max_bands_read==0) then
       max_bands_read = 1
       if (peinf%inode==0) then
          write(0,*)
          write(0,'(a)') 'WARNING: could not honor limit of '
          write(0,'(f0.3,a)') max_bytes_read/1024d0/1024d0,' MB per chunk when'
          write(0,'(a)') 'reading HDF5 WFN file. Using chunks of '
          write(0,'(f0.3,a)') (kp%nspin*kp%nspinor*2*dble(ngktot)*8d0)/1024d0/1024d0,' MB.'
          write(0,*)
       endif
    endif
    !> wfndata(2,ig/ik,is,ib)
    allocate(wfndata (2, ngktot, kp%nspin*kp%nspinor, max_bands_read))
    ib = 1
    do while (ib<=nbownmax) ! nbownmax = peinf%nvownmax
       ! bands_read is local
       ! bands_read_max is global, in other word, MAX_ACROSS_PROCS[bands_read] = bands_read_max
       ! There must be at least one procs which actually reads bands_read_max bands in one iteration, which means bands_read = bands_read_max for these procs.
       ! if for current proc, we have bands_read<bands_read_max, it means this is the last reading operation, after this iteration, this LOOP will end.
       bands_read = max(min(nbownactual, ib-1 + max_bands_read) - ib + 1, 0)
       bands_read_max = max(min(nbownmax, ib-1 + max_bands_read) - ib + 1, 0)
       count(1) = 2
       count(2) = ngktot
       count(3) = kp%nspin*kp%nspinor
       count(4) = bands_read
       call h5screate_simple_f(4, count, memspace, error)
       ! if comm_style .ne. 0, then only the reader-th procs will have do_read=.true.
       ! <->
       ! if comm_style == 0, all the procs will participate in the reading
       ! ==> do_read always be .true., and all the procs in the same pool will read in the
       ! same wavefunctions
       do_read = bands_read>0.and.(peinf%inode==reader.or.comm_style==0)
       if (do_read) then
          offset(1) = 0
          offset(2) = 0
          offset(3) = 0
          offset(4) = (ib_first-1)+ioffset_+(ib-1)
          if (peinf%verb_debug .and. peinf%inode==reader) then
             write(6,'(4(a,i0,1x))') 'ib=',ib,'bands_read=',bands_read,'offset=',offset(3),'ib_first=',ib_first
          endif
       else
          offset(:) = 0
          call H5sselect_none_f(memspace,error)
       endif
       call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, offset, count, error)
       if (.not.do_read) then
          call H5sselect_none_f(dataspace,error)
       endif
       if (peinf%verb_debug .and. peinf%inode==reader) then
          write(6,'(a,i0,a)') 'ib=',ib,' before read!'
       endif
       call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
       call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
       call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, wfndata(:,:,:,:), count, error, memspace, dataspace, xfer_prp=plist_id)
       call h5pclose_f(plist_id, error)
       if (peinf%verb_debug .and. peinf%inode==reader) then
          write(6,'(a,i0,a)') 'ib=',ib,' after read!'
       endif
       call h5sclose_f(memspace, error)
       ! for those reader-th procs
       if (do_read) then
          do is = 1, kp%nspin*kp%nspinor
             wfnsout(:,is,ib:ib+bands_read-1) = DCMPLX(wfndata(1,:,is,1:bands_read),wfndata(2,:,is,1:bands_read))
          enddo
       endif
       if (bands_read>0) then
          ! FHJ: No manual distribution is necessary for comm_style==0
          if (comm_style>0) call logit('Sending bands')
          if (comm_style==2) then
             if (peinf%inode==reader) then
                do ipe = 1, peinf%npes
                   if(does_it_ownb(ib_first,ipe) .and. ipe-1 .ne. peinf%inode) then
                      call MPI_Send(wfnsout(1,1,ib), ngktot*bands_read*kp%nspin*kp%nspinor, MPI_DOUBLE_COMPLEX, ipe-1, 0, MPI_COMM_WORLD, mpierr)
                   endif
                enddo
             else
                call MPI_Recv(wfnsout(1,1,ib), ngktot*bands_read*kp%nspin*kp%nspinor, MPI_DOUBLE_COMPLEX, reader, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierr)
             endif
          elseif (comm_style==1) then
             call MPI_Bcast(wfnsout(1,1,ib), ngktot*bands_read*kp%nspin*kp%nspinor, MPI_DOUBLE_COMPLEX, 0, bandcomm, mpierr)
          endif
       endif
       ib = ib + bands_read_max
    enddo
    if (comm_style==1.and.nbownactual>0) then
       call MPI_Comm_free(bandcomm, mpierr)
       call MPI_Group_free(mpigrp, mpierr)
    endif
    if (allocated(ranks)) then
       deallocate(ranks)
    endif
    if(allocated(wfndata)) then
       deallocate(wfndata)
    endif
    call h5sclose_f(dataspace, error)
    call h5dclose_f(dset_id, error)

  end subroutine read_hdf5_bands_block
  ! read/write other

  !===============================================================================
  !> Create the appropriate structures that hold info about the mean-field
  subroutine setup_hdf5_mf_file(fname, create_file)
    character(len=*), intent(in) :: fname
    logical, intent(in), optional :: create_file
    integer(HID_T) :: file_id
    integer(HID_T) :: group_id
    integer :: error
    logical :: create_file_

    create_file_ = .true.
    if (present(create_file)) create_file_ = create_file
    if (create_file_) then
       call h5fcreate_f(fname, H5F_ACC_TRUNC_F, file_id, error)
    else
       call h5fopen_f(fname, H5F_ACC_RDWR_F, file_id, error)
    endif
    call h5gcreate_f(file_id, '/mf_header', group_id, error)
    call h5gclose_f(group_id, error)
    call h5gcreate_f(file_id, '/mf_header/gspace', group_id, error)
    call h5gclose_f(group_id, error)
    call h5gcreate_f(file_id, '/mf_header/symmetry', group_id, error)
    call h5gclose_f(group_id, error)
    call h5gcreate_f(file_id, '/mf_header/crystal', group_id, error)
    call h5gclose_f(group_id, error)
    call h5gcreate_f(file_id, '/mf_header/kpoints', group_id, error)
    call h5gclose_f(group_id, error)
    call h5fclose_f(file_id, error)

  end subroutine setup_hdf5_mf_file
  !> Create the appropriate structures that hold info about the WFN
  !! coefficients in an 1 file. No data is actually written.
  subroutine setup_hdf5_wfn_file(fname, iflavor, kp)
    character(len=*), intent(in) :: fname
    integer, intent(in) :: iflavor
    type(kpoints), intent(in) :: kp
    integer(HID_T) :: file_id
    integer(HID_T) :: group_id
    integer(HID_T) :: dataspace_id
    integer(HID_T) :: dataset_id
    integer :: error
    integer, parameter :: rank_wfn = 5
    integer(HSIZE_T) :: dataspace_wfn(5)
    ! Data chunking
    integer(HID_T) :: dcpl
    integer, parameter :: rank_chunk = 5
    integer(HSIZE_T) :: chunk(5)

    if (iflavor == 1) then
       call die("Only support iflavor = 2 case")
    endif
    call setup_hdf5_mf_file(fname)
    call h5fopen_f(fname, H5F_ACC_RDWR_F, file_id, error)
    call h5gcreate_f(file_id, '/wfns', group_id, error)
    call h5gclose_f(group_id, error)
    if (iflavor/=1 .and. iflavor/=2) then
       write(0,*) 'ERROR: got iflavor=', iflavor
       call die('Internal error: invalid flavor in setup_hdf5_wfn_file.', only_root_writes=.true.)
    endif
    dataspace_wfn(1) = iflavor ! complex DP ic
    dataspace_wfn(2) = kp%ngkmax
    dataspace_wfn(3) = kp%mnband ! ib
    dataspace_wfn(4) = kp%nspinor*kp%nspin ! is
    dataspace_wfn(5) = kp%nrk ! ik
    call h5screate_simple_f(rank_wfn, dataspace_wfn, dataspace_id, error)
    call h5pcreate_f(H5P_DATASET_CREATE_F, dcpl, error)
    chunk(1) = 2
    chunk(2) = kp%ngkmax
    chunk(3) = 1
    chunk(4) = kp%nspinor*kp%nspin
    chunk(5) = 1
    call h5pset_chunk_f(dcpl,rank_chunk, chunk, error)
    call h5dcreate_f(file_id, '/wfns/coeffs', H5T_NATIVE_DOUBLE, dataspace_id, dataset_id, error, dcpl)
    call h5sclose_f(dataspace_id, error)
    call h5dclose_f(dataset_id, error)
    call h5pclose_f(dcpl , error)
    call h5fclose_f(file_id, error)

  end subroutine setup_hdf5_wfn_file
  !> A high-level wrapper for write_*_header* functions
  subroutine write_hdf5_mf_header(fname, mf)
    character(len=*), intent(in) :: fname
    type(mf_header_t), intent(in) :: mf
    character(len=3) :: sheader

    sheader = mf%sheader
    call write_hdf5_header_type(fname, sheader, mf%iflavor, mf%kp, mf%gvec, mf%syms, mf%crys)

  end subroutine write_hdf5_mf_header
  !> A high-level wrapper for read_*_header* functions
  !! Note: the optional fields are ignored for now.
  subroutine read_hdf5_mf_header(fname, mf, iflavor, sheader, warn, dont_warn_kgrid)
    character(len=*), intent(in) :: fname
    type(mf_header_t), intent(out) :: mf
    integer, intent(in), optional :: iflavor
    character(len=3), intent(in), optional :: sheader
    logical, intent(in), optional :: warn
    logical, intent(in), optional :: dont_warn_kgrid

    if (present(sheader)) then
       mf%sheader = sheader
    else
       mf%sheader = 'GET'
    endif
    if (present(iflavor)) then
       mf%iflavor = iflavor
    else
       mf%iflavor = -1
    endif
    call read_hdf5_header_type(fname, mf%sheader, mf%iflavor, mf%kp, mf%gvec, mf%syms, mf%crys)
    !mf%kp, mf%gvec, mf%syms, mf%crys, version=mf%version, sdate=mf%sdate, stime=mf%stime)
    !FHJ: FIXME - Implement in WFN file
    mf%sheader = 'WFN'
    mf%stime = 'N/A'
    mf%sdate = 'N/A'

  end subroutine read_hdf5_mf_header
  subroutine read_hdf5_gvec_kpt(filename, ngk, gvec_kpt_offset, components)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: ngk
    integer, intent(in) :: gvec_kpt_offset
    integer, intent(out) :: components(:,:) !< (id=1:nd,ig=1:ngk)
    integer(hid_t) :: file_id_wfn
    integer(hid_t) :: dataset_id_gvec
    integer(hid_t) :: dataspace_gvec
    integer(hid_t) :: memspace_gvec
    integer, parameter :: rank_gvec=2
    integer(hsize_t) :: count_gvec(2), offset_gvec(2)
    integer :: error
    if ( peinf%inode .eq. 0 ) then
       call h5fopen_f(trim(adjustl(filename)), H5F_ACC_RDONLY_F, file_id_wfn, error)
       call h5dopen_f(file_id_wfn, '/wfns/gvecs', dataset_id_gvec, error)
       call h5dget_space_f(dataset_id_gvec, dataspace_gvec, error)
       count_gvec(1) = 3
       count_gvec(2) = ngk
       call h5screate_simple_f(rank_gvec, count_gvec, memspace_gvec, error)
       offset_gvec(1) = 0
       offset_gvec(2) = gvec_kpt_offset
       components(:,:) = 0
       call h5sselect_hyperslab_f(dataspace_gvec, h5s_select_set_f, offset_gvec, count_gvec, error)
       call h5dread_f(dataset_id_gvec, h5t_native_integer, components(1:3,1:ngk), count_gvec, error, memspace_gvec, dataspace_gvec)
       call h5sclose_f(memspace_gvec, error)
       call h5sclose_f(dataspace_gvec, error)
       call h5dclose_f(dataset_id_gvec, error)
       call h5fclose_f(file_id_wfn, error)
    endif
    if ( peinf%npes > 1 ) then
       call mpi_bcast(components(1,1), 3*ngk, mpi_integer, 0, mpi_comm_world, mpierr)
    endif
  end subroutine read_hdf5_gvec_kpt
end module wfn_io_hdf5_m

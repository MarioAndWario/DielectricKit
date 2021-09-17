!=========================================================================
!
! Included from file wfn_io_hdf5.f90.
! You are not expected to understand this. --JIM
!
!=========================================================================

#ifdef READ
#define HDF5_READ_WRITE(x) hdf5_read ## x
#define READ_WRITE(x) read ## x
#define INTENT out
#define FLAVOR_INTENT inout
#define H5D_READ_WRITE call h5dread_f
#define H5F_FILE_ACCESS H5F_ACC_RDONLY_F
#else
#define HDF5_READ_WRITE(x) hdf5_write ## x
#define READ_WRITE(x) write ## x
#define INTENT in
#define FLAVOR_INTENT in
#define H5D_READ_WRITE call h5dwrite_f
#define H5F_FILE_ACCESS H5F_ACC_RDWR_F
#endif

#ifdef HDF5
#define NAME(x) READ_WRITE(_hdf5 ## x)
#endif

#ifdef TEMP_WFN_DATA
#ifdef TEMP_COMPLEX
#define TEMP_SCALAR complex(DPC)
#define MPI_TEMP_SCALAR MPI_COMPLEX_DPC
#define LONGNAME(x) NAME(x ## _complex)
#else
#define TEMP_SCALAR real(DPC)
#define MPI_TEMP_SCALAR MPI_REAL_DP
#define LONGNAME(x) NAME(x ## _real)
#endif
#endif

! begin read/write header
#ifdef TEMP_HEADER
subroutine NAME(_header_type)(sFileName, sheader, iflavor, kp, gvec, syms, crys)
  character(len=*), intent(in) :: sFileName
  character(len=3), intent(inout) :: sheader
  integer, intent(FLAVOR_INTENT) :: iflavor
  type(kpoints), intent(INTENT) :: kp
  type(gspace), intent(INTENT) :: gvec
  type(symmetry), intent(INTENT) :: syms
  type(crystal), intent(INTENT) :: crys

  ! set values based on epsilon calculation
  logical :: is_get=.false.
  logical :: wfnflag=.true.

  PUSH_SUB(NAME(_header_type))

  if (peinf%inode == 0) then
     call READ_WRITE(_info)(TRUNC(sFileName),iflavor)
     call READ_WRITE(_kpoints)(TRUNC(sFileName),kp)
     call READ_WRITE(_gspace)(TRUNC(sFileName),gvec)
     call READ_WRITE(_symmetry)(TRUNC(sFileName),syms)
     call READ_WRITE(_crystal)(TRUNC(sFileName),crys)
  endif

#ifdef READ
  ! ======
  gvec%nFFTgridpts = product(gvec%FFTgrid(1:3))
  ! ------
#ifdef MPI

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
     ! =====
     call MPI_BCAST(gvec%nFFTgridpts, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
     ! -----
     call MPI_BCAST(crys%celvol, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
     call MPI_BCAST(crys%alat, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
     call MPI_BCAST(crys%avec, 9, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
     call MPI_BCAST(crys%adot, 9, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
     call MPI_BCAST(crys%recvol, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
     call MPI_BCAST(crys%blat, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
     call MPI_BCAST(crys%bvec, 9, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
     call MPI_BCAST(crys%bdot, 9, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
     call MPI_BCAST(syms%mtrx(1,1,1), 3*3*48, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
     ! ======
     call MPI_BCAST(syms%mtrx_reci(1,1,1), 3*3*48, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
     call MPI_BCAST(syms%mtrx_cart(1,1,1), 3*3*48, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
     ! ======
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
#endif
  if(peinf%inode > 0 ) then
     SAFE_ALLOCATE(crys%atyp, (crys%nat))
     SAFE_ALLOCATE(crys%apos, (3, crys%nat))
     if (wfnflag) then
        SAFE_ALLOCATE(kp%ngk, (kp%nrk))
        SAFE_ALLOCATE(kp%w, (kp%nrk))
        SAFE_ALLOCATE(kp%rk, (3, kp%nrk))
        SAFE_ALLOCATE(kp%ifmin, (kp%nrk, kp%nspin))
        SAFE_ALLOCATE(kp%ifmax, (kp%nrk, kp%nspin))
        SAFE_ALLOCATE(kp%el, (kp%mnband, kp%nrk, kp%nspin))
        SAFE_ALLOCATE(kp%occ, (kp%mnband, kp%nrk, kp%nspin))
     endif
  endif
#endif

#if defined READ && defined MPI
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
#endif
  POP_SUB(NAME(_header_type))
end subroutine NAME(_header_type)

subroutine READ_WRITE(_info)(sFileName, iflavor)
  character(len=*), intent(in) :: sFileName
  integer, intent(INTENT) :: iflavor

  integer(HID_T) :: hidFile
  integer :: iError
#ifdef READ
  logical :: exists
#endif

  PUSH_SUB(READ_WRITE(_info))

  call h5fopen_f(sFileName, H5F_FILE_ACCESS, hidFile, iError)
#ifndef READ
  call HDF5_READ_WRITE(_int)(hidFile, '/mf_header/versionnumber', VER_WFN_HDF5, iError)
#endif
#ifdef READ
  !FHJ: Keep backwards compatibility with older WFN.h5 files
  call h5lexists_f(hidFile, '/mf_header/flavor', exists, iError)
  if (exists .and. iError==0) then
     call HDF5_READ_WRITE(_int)(hidFile, '/mf_header/flavor', iflavor, iError)
  else
     call HDF5_READ_WRITE(_int)(hidFile, '/info', iflavor, iError)
  endif
#else
  call HDF5_READ_WRITE(_int)(hidFile, '/mf_header/flavor', iflavor, iError)
#endif
  call h5fclose_f(hidFile, iError)

  POP_SUB(READ_WRITE(_info))
end subroutine READ_WRITE(_info)

subroutine READ_WRITE(_gspace)(sFileName,gvec)
  character(len=*), intent(in) :: sFileName
  type(gspace), intent(INTENT) :: gvec

  integer(HID_T) :: hidFile
  integer :: iError

  PUSH_SUB(READ_WRITE(_gspace))

  call h5fopen_f(sFileName, H5F_FILE_ACCESS, hidFile, iError)
  call HDF5_READ_WRITE(_int)(hidFile, '/mf_header/gspace/ng', gvec%ng, iError)
  call HDF5_READ_WRITE(_double)(hidFile, '/mf_header/gspace/ecutrho', gvec%ecutrho, iError)
  call HDF5_READ_WRITE(_int_array)(hidFile, '/mf_header/gspace/FFTgrid', (/3/), gvec%FFTgrid, iError)
  call h5fclose_f(hidFile, iError)

  POP_SUB(READ_WRITE(_gspace))
end subroutine READ_WRITE(_gspace)

subroutine NAME(_gvectors)(sFileName, ng, gvec)
  character(len=*), intent(in) :: sFileName
  integer, intent(in) :: ng !< used size of array
  integer, intent(INTENT) :: gvec(:, :) !< (3, ng_bound)

  integer(HID_T) :: hidFile
  integer :: iError
  logical :: bcast_, dont_read_

  PUSH_SUB(NAME(_gvectors))

  dont_read_=.false.
  bcast_=.not. dont_read_

  if(peinf%inode == 0) then
     call h5fopen_f(sFileName, H5F_FILE_ACCESS, hidFile, iError)
     call HDF5_READ_WRITE(_int_array)(hidFile, '/mf_header/gspace/components', (/3,ng/), gvec, iError)
     call h5fclose_f(hidFile, iError)
  endif

#ifdef READ
#ifdef MPI
  if(peinf%npes > 1) then
     if(bcast_) then
        call MPI_BCAST(gvec(1,1), 3 * ng, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
     endif
  endif
#endif
#endif

  POP_SUB(NAME(_gvectors))
end subroutine NAME(_gvectors)

subroutine READ_WRITE(_symmetry)(sFileName,syms)
  character(len=*), intent(in) :: sFileName
  type(symmetry), intent(INTENT) :: syms

  integer(HID_T) :: hidFile
  integer :: iError

  PUSH_SUB(READ_WRITE(_symmetry))

  call h5fopen_f(sFileName, H5F_FILE_ACCESS, hidFile, iError)
  call HDF5_READ_WRITE(_int)(hidFile, '/mf_header/symmetry/ntran', syms%ntran, iError)
  call HDF5_READ_WRITE(_int)(hidFile, '/mf_header/symmetry/cell_symmetry', syms%cell_symmetry, iError)
  call HDF5_READ_WRITE(_int_array)(hidFile, '/mf_header/symmetry/mtrx', (/3, 3, 48/), syms%mtrx, iError)
  ! ======
  call HDF5_READ_WRITE(_int_array)(hidFile, '/mf_header/symmetry/mtrx_reci', (/3, 3, 48/), syms%mtrx_reci, iError)
  call HDF5_READ_WRITE(_double_array)(hidFile, '/mf_header/symmetry/mtrx_cart', (/3, 3, 48/), syms%mtrx_cart, iError)
  ! ======
  call HDF5_READ_WRITE(_double_array)(hidFile, '/mf_header/symmetry/tnp', (/3, 48/), syms%tnp, iError)
  call h5fclose_f(hidFile, iError)

  POP_SUB(READ_WRITE(_symmetry))
end subroutine READ_WRITE(_symmetry)

subroutine READ_WRITE(_crystal)(sFileName,crys)
  character(len=*), intent(in) :: sFileName
  type(crystal), intent(INTENT) :: crys

  integer(HID_T) :: hidFile
  integer :: iError

  PUSH_SUB(READ_WRITE(_crystal))

  call h5fopen_f(sFileName, H5F_FILE_ACCESS, hidFile, iError)
  call HDF5_READ_WRITE(_double)(hidFile, '/mf_header/crystal/celvol', crys%celvol, iError)
  call HDF5_READ_WRITE(_double)(hidFile, '/mf_header/crystal/recvol', crys%recvol, iError)
  call HDF5_READ_WRITE(_double)(hidFile, '/mf_header/crystal/alat', crys%alat, iError)
  call HDF5_READ_WRITE(_double)(hidFile, '/mf_header/crystal/blat', crys%blat, iError)
  call HDF5_READ_WRITE(_int)(hidFile, '/mf_header/crystal/nat', crys%nat, iError)
  call HDF5_READ_WRITE(_double_array)(hidFile, '/mf_header/crystal/avec', (/3, 3/), crys%avec, iError)
  call HDF5_READ_WRITE(_double_array)(hidFile, '/mf_header/crystal/bvec', (/3, 3/), crys%bvec, iError)
  call HDF5_READ_WRITE(_double_array)(hidFile, '/mf_header/crystal/adot', (/3, 3/), crys%adot, iError)
  call HDF5_READ_WRITE(_double_array)(hidFile, '/mf_header/crystal/bdot', (/3, 3/), crys%bdot, iError)

#ifdef READ
  SAFE_ALLOCATE(crys%atyp, (crys%nat))
  SAFE_ALLOCATE(crys%apos, (3,crys%nat))
#endif

  call HDF5_READ_WRITE(_int_array)(hidFile, '/mf_header/crystal/atyp', (/crys%nat/), crys%atyp, iError)
  call HDF5_READ_WRITE(_double_array)(hidFile, '/mf_header/crystal/apos', (/3, crys%nat/), crys%apos, iError)
  call h5fclose_f(hidFile, iError)

  POP_SUB(READ_WRITE(_crystal))
end subroutine READ_WRITE(_crystal)

subroutine READ_WRITE(_kpoints)(sFileName,kp)
  character(len=*), intent(in) :: sFileName
  type(kpoints), intent(INTENT) :: kp

  integer(HID_T) :: hidFile
  integer :: iError

  PUSH_SUB(READ_WRITE(_kpoints))

  call h5fopen_f(sFileName, H5F_FILE_ACCESS, hidFile, iError)
  call HDF5_READ_WRITE(_int)(hidFile, '/mf_header/kpoints/nspin', kp%nspin, iError)
  call HDF5_READ_WRITE(_int)(hidFile, '/mf_header/kpoints/nspinor', kp%nspinor, iError)
  call HDF5_READ_WRITE(_int)(hidFile, '/mf_header/kpoints/nrk', kp%nrk, iError)
  call HDF5_READ_WRITE(_int)(hidFile, '/mf_header/kpoints/mnband', kp%mnband, iError)
  call HDF5_READ_WRITE(_int)(hidFile, '/mf_header/kpoints/ngkmax', kp%ngkmax, iError)
  call HDF5_READ_WRITE(_double)(hidFile, '/mf_header/kpoints/ecutwfc', kp%ecutwfc, iError)
  call HDF5_READ_WRITE(_int_array)(hidFile, '/mf_header/kpoints/kgrid', (/3/), kp%kgrid, iError)
  call HDF5_READ_WRITE(_double_array)(hidFile, '/mf_header/kpoints/shift', (/3/), kp%shift, iError)

#ifdef READ
  SAFE_ALLOCATE(kp%ngk, (kp%nrk))
  SAFE_ALLOCATE(kp%ifmin, (kp%nrk, kp%nspin))
  SAFE_ALLOCATE(kp%ifmax, (kp%nrk, kp%nspin))
  SAFE_ALLOCATE(kp%w, (kp%nrk))
  SAFE_ALLOCATE(kp%rk, (3,kp%nrk))
  SAFE_ALLOCATE(kp%el, (kp%mnband,kp%nrk,kp%nspin))
  SAFE_ALLOCATE(kp%occ, (kp%mnband,kp%nrk,kp%nspin))
#endif

  call HDF5_READ_WRITE(_int_array)(hidFile, '/mf_header/kpoints/ngk', &
       (/kp%nrk/), kp%ngk, iError)
  call HDF5_READ_WRITE(_int_array)(hidFile, '/mf_header/kpoints/ifmin', &
       (/kp%nrk, kp%nspin/), kp%ifmin, iError)
  call HDF5_READ_WRITE(_int_array)(hidFile, '/mf_header/kpoints/ifmax', &
       (/kp%nrk, kp%nspin/), kp%ifmax, iError)
  call HDF5_READ_WRITE(_double_array)(hidFile, '/mf_header/kpoints/w', &
       (/kp%nrk/), kp%w, iError)
  call HDF5_READ_WRITE(_double_array)(hidFile, '/mf_header/kpoints/rk', &
       (/3, kp%nrk/), kp%rk, iError)
  call HDF5_READ_WRITE(_double_array)(hidFile, '/mf_header/kpoints/el', &
       (/kp%mnband, kp%nrk, kp%nspin/), kp%el, iError)
  call HDF5_READ_WRITE(_double_array)(hidFile, '/mf_header/kpoints/occ', &
       (/kp%mnband, kp%nrk, kp%nspin/), kp%occ, iError)
  call h5fclose_f(hidFile, iError)

  POP_SUB(READ_WRITE(_kpoints))
end subroutine READ_WRITE(_kpoints)
#endif
! end read/write header

!#BEGIN_INTERNAL_ONLY
! begin read/write wfn gvectors
#ifdef TEMP_WFN_GVEC
subroutine NAME(_wfn_gvectors)(fname, gvec, ngktot)
  character(len=*) :: fname
  integer, intent(inout) :: gvec(:,:)
  integer, intent(in) :: ngktot

  integer(HID_T) :: file_id
  integer :: error

  PUSH_SUB(NAME(_wfn_gvectors))

  if(peinf%inode == 0) then
     call h5fopen_f(fname, H5F_FILE_ACCESS, file_id, error)
     call HDF5_READ_WRITE(_int_array)(file_id, 'wfns/gvecs', (/3, ngktot/), gvec, error)
     call h5fclose_f(file_id, error)
  endif

#if defined READ && defined MPI
  if(peinf%npes > 1) then
     call MPI_Bcast(gvec(1,1), 3*ngktot, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
  endif
#endif
  POP_SUB(NAME(_wfn_gvectors))
end subroutine NAME(_wfn_gvectors)
#endif
! end read/write wfn gvectors

! begin read/write wfn data
#ifdef TEMP_WFN_DATA
subroutine LONGNAME(_band)(fname, wfn, ngk, nstot, ioffsetk, ioffsetb)
  character(len=*) :: fname
  TEMP_SCALAR, intent(INTENT) :: wfn(:,:) !< (ngk,nstot=kp%ns*kp%nspinor)
  integer, intent(in) :: ngk
  integer, intent(in) :: nstot
  integer, intent(in) :: ioffsetk
  integer, intent(in) :: ioffsetb

#ifdef TEMP_COMPLEX
  real(DP) :: dwfn(2,ngk,nstot,1)
#else
  real(DP) :: dwfn(1,ngk,nstot,1)
#endif
  integer(HID_T) :: file_id
  integer(HID_T) :: dataset_id
  integer(HID_T) :: dataspace_id
  integer(HID_T) :: memspace_id
  integer(HSIZE_T) :: a3(4), offset(4), count(4)
  integer :: error

  PUSH_SUB(LONGNAME(_band))

#ifdef TEMP_COMPLEX
  a3(1) = 2
#else
  a3(1) = 1
#endif
  a3(2) = ngk
  a3(3) = nstot
  a3(4) = 1
  offset(1) = 0
  offset(2) = ioffsetk
  offset(3) = 0
  offset(4) = ioffsetb

#ifdef TEMP_COMPLEX
  count(1) = 2
#else
  count(1) = 1
#endif
  count(2) = ngk
  count(3) = nstot
  count(4) = 1

#ifndef READ
#ifdef TEMP_COMPLEX
  dwfn(1,:,:,1) = real(wfn(:,:))
  dwfn(2,:,:,1) = IMAG(wfn(:,:))
#else
  dwfn(1,:,:,1) = wfn(:,:)
#endif
#endif

  if(peinf%inode == 0) then
     call h5fopen_f(fname, H5F_FILE_ACCESS, file_id, error)
     call h5dopen_f(file_id, 'wfns/coeffs', dataset_id, error)
     CALL h5screate_simple_f(4, count, memspace_id, error)
     call h5dget_space_f(dataset_id, dataspace_id, error)

     call h5sselect_hyperslab_f(dataspace_id, H5S_SELECT_SET_F, offset, count, error)

     H5D_READ_WRITE(dataset_id, H5T_NATIVE_DOUBLE, dwfn, a3, error, file_space_id = dataspace_id, mem_space_id = memspace_id)
     call h5dclose_f(dataset_id, error)
     call h5sclose_f(dataspace_id, error)
     call h5sclose_f(memspace_id, error)
     call h5fclose_f(file_id, error)
  endif

#ifdef READ
#ifdef TEMP_COMPLEX
  wfn(:,:) = DCMPLX(dwfn(1,:,:,1), dwfn(2,:,:,1))
#else
  wfn(:,:) = dwfn(1,:,:,1)
#endif
#endif

  POP_SUB(LONGNAME(_band))

end subroutine LONGNAME(_band)
#endif
! end read/write wfn data

#ifdef TEMP_OTHER
! ======
! from Epsilon/input.f90:
!   call read_hdf5_bands_block(file_id, kp, peinf%nvownmax, peinf%nvownactual, &
!                   peinf%does_it_ownv, ib_first, wfns, ioffset=vwfn%ncore_excl)
! ------
! nbownmax = peinf%nvownmax
! nbownactual = peinf%nvownactual
! does_it_ownb = peinf%does_it_ownv
! ib_first = peinf%invindexv(1) !peinf%invindexc(icb_relative_local) = icb_relative_global
! wfns = wfnsout
! ioffset = vwfn%ncore_excl
subroutine read_hdf5_bands_block(file_id, kp, nbownmax, nbownactual, does_it_ownb, ib_first, wfnsout, ioffset)
  integer(HID_T), intent(in) :: file_id
  type(kpoints), intent(in) :: kp
  integer, intent(in) :: nbownmax
  integer, intent(in) :: nbownactual !< how many bands I own
  logical, intent(in) :: does_it_ownb(:,:)
  integer, intent(in) :: ib_first !< first band of the set of bands I own
  SCALAR, intent(out) :: wfnsout(:,:,:) !< (SUM(kp%ngk), kp%nspin*kp%nspinor, nbownactual)
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
  ! 0=native HDF5; 1=manual group comm; 2=manual send/recvs
  integer, parameter :: comm_style=0
  logical :: do_read

  PUSH_SUB(read_hdf5_bands_block)

  call logit('Reading HDF5 WFNs by blocks')

  SAFE_ALLOCATE(ranks,(peinf%npes))
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
  ! ======
  ! The following IF is useless when comm_style = 0
  ! <->
  if (nbownactual>0) then
     do ipe = 1, peinf%npes
        !
        if(does_it_ownb(ib_first,ipe)) then
           reader = ipe - 1
           exit ! when we find one, we exit, so this is the lowest rank !!

        endif
     enddo

     if (reader==-1) call die("Cannot find reader in read_hdf5_bands_block", only_root_writes=.true.)

#ifdef MPI
     ! >-<
     if (comm_style==1) then
        ! We use MPI_Bcast with MPI groups
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
#endif
     ! >-<
  else
#ifdef MPI
     ! >-<
     if (comm_style==1) then
        ! FHJ: Note that MPI_Comm_Create must be called by everyone in MPI_COMM_WORLD!
        call MPI_Comm_Create(MPI_COMM_WORLD, MPI_GROUP_EMPTY, bandcomm, mpierr)
     endif
#endif
  endif

  ! FHJ: read at most max_bands_read bands to avoid MPI/HDF5 buffer overflow.
  ! Here, SCALARSIZE*kp%nspin*kp%nspinor*dble(ngktot)*8d0 is the size of a
  ! single band, including G-vectors from all k-points.
  ! ======
  ! Common/f_defs.h
  ! #ifdef CPLX
  ! #define SCALARSIZE 2
  ! #else
  ! #define SCALARSIZE 1
  ! #endif
  ! ------
  ! max_bands_read < nbownmax = peinf%nvownmax
  max_bands_read = min(nbownmax, &
       int(max_bytes_read/(SCALARSIZE*kp%nspin*kp%nspinor*dble(ngktot)*8d0)))

  ! >-<
  if (max_bands_read==0) then
     max_bands_read = 1
     if (peinf%inode==0) then
        write(0,*)
        write(0,'(a)') 'WARNING: could not honor limit of '
        write(0,'(f0.3,a)') max_bytes_read/1024d0/1024d0,' MB per chunk when'
        write(0,'(a)') 'reading HDF5 WFN file. Using chunks of '
        write(0,'(f0.3,a)') (kp%nspin*kp%nspinor*SCALARSIZE*dble(ngktot)*8d0)/1024d0/1024d0,' MB.'
        write(0,*)
     endif
  endif

  !> wfndata(2,ig/ik,is,ib)
  SAFE_ALLOCATE(wfndata, (SCALARSIZE, ngktot, kp%nspin*kp%nspinor, max_bands_read))

  ib = 1
  do while (ib<=nbownmax) ! nbownmax = peinf%nvownmax
     ! bands_read is local
     ! bands_read_max is global, in other word, MAX_ACROSS_PROCS[bands_read] = bands_read_max
     ! There must be at least one procs which actually reads bands_read_max bands in one iteration, which means bands_read = bands_read_max for these procs.
     ! if for current proc, we have bands_read<bands_read_max, it means this is the last reading operation, after this iteration, this LOOP will end.
     bands_read = max(min(nbownactual, ib-1 + max_bands_read) - ib + 1, 0)
     bands_read_max = max(min(nbownmax, ib-1 + max_bands_read) - ib + 1, 0)
     count(1) = SCALARSIZE
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

#ifdef MPI
     call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
     call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
     call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, wfndata(:,:,:,:), count, error, memspace, dataspace, xfer_prp=plist_id)
     call h5pclose_f(plist_id, error)
#else
     call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, wfndata(:,:,:,:), count, error, memspace, dataspace)
#endif

     if (peinf%verb_debug .and. peinf%inode==reader) then
        write(6,'(a,i0,a)') 'ib=',ib,' after read!'
     endif

     call h5sclose_f(memspace, error)

     ! for those reader-th procs
     if (do_read) then
#ifdef DEBUG
        if (peinf%verb_max .and. peinf%inode==reader) then
           do ib_ = 1, bands_read
              write(6,'(" band = ",i6,"; avg(norm) = ", f12.6)') &
                   offset(3) + ib_, sum(wfndata(:,:,:,ib_)**2)/dble(kp%nrk)
           enddo
           FLUSH(6)
        endif
#endif
        do is = 1, kp%nspin*kp%nspinor
           wfnsout(:,is,ib:ib+bands_read-1) = SCALARIFY2(wfndata(1,:,is,1:bands_read),wfndata(2,:,is,1:bands_read))
        enddo
     endif

#ifdef MPI
     if (bands_read>0) then
        ! FHJ: No manual distribution is necessary for comm_style==0
        ! >-<
        if (comm_style>0) call logit('Sending bands')
        ! >-<
        if (comm_style==2) then
           if (peinf%inode==reader) then
              do ipe = 1, peinf%npes
                 if(does_it_ownb(ib_first,ipe) .and. ipe-1 .ne. peinf%inode) then
                    call MPI_Send(wfnsout(1,1,ib), ngktot*bands_read*kp%nspin*kp%nspinor, MPI_SCALAR, ipe-1, 0, MPI_COMM_WORLD, mpierr)
                 endif
              enddo
           else
              call MPI_Recv(wfnsout(1,1,ib), ngktot*bands_read*kp%nspin*kp%nspinor, MPI_SCALAR, reader, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierr)
           endif
           ! >-<
        elseif (comm_style==1) then
           call MPI_Bcast(wfnsout(1,1,ib), ngktot*bands_read*kp%nspin*kp%nspinor, MPI_SCALAR, 0, bandcomm, mpierr)
        endif
     endif
#endif
     ib = ib + bands_read_max
  enddo

#ifdef MPI
  if (comm_style==1.and.nbownactual>0) then
     call MPI_Comm_free(bandcomm, mpierr)
     call MPI_Group_free(mpigrp, mpierr)
  endif
#endif
  SAFE_DEALLOCATE(ranks)
  SAFE_DEALLOCATE(wfndata)

  call h5sclose_f(dataspace, error)
  call h5dclose_f(dset_id, error)
  POP_SUB(read_hdf5_bands_block)
end subroutine read_hdf5_bands_block
#endif
! read/write other
!#END_INTERNAL_ONLY


#undef READ_WRITE
#undef INTENT
#undef FLAVOR_INTENT
#undef NAME
#undef TEMP_SCALAR
#undef MPI_TEMP_SCALAR
#undef LONGNAME

#undef HDF5_READ_WRITE
#undef H5D_READ_WRITE
#undef H5F_FILE_ACCESS

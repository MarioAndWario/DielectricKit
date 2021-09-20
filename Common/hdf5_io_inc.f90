!
!==============================================================================
! Included from file hdf5_io_hdf5.f90. Originally by JIM.
! Last modified 12/2014 (FHJ)
!==============================================================================
#ifdef READ
#define READ_WRITE(x) hdf5_read ## x
#define INTENT inout
#define H5F_FILE_ACCESS H5F_ACC_RDONLY_F
#else
#define READ_WRITE(x) hdf5_write ## x
#define INTENT in
#define H5F_FILE_ACCESS H5F_ACC_RDWR_F
#endif

#ifdef TYPE_INT
#define THE_TYPE integer
#define SUBROUTINE_NAME_SCALAR READ_WRITE(_int)
#define SUBROUTINE_NAME_ARRAY READ_WRITE(_int_array)
#define SUBROUTINE_NAME_HYPERSLAB READ_WRITE(_int_hyperslab)
#define TYPE_SUFFIX(x) x ## _int_f
#define H5TYPE H5T_NATIVE_INTEGER
#endif
#ifdef TYPE_DOUBLE
#define THE_TYPE real(DP)
#define SUBROUTINE_NAME_SCALAR READ_WRITE(_double)
#define SUBROUTINE_NAME_ARRAY READ_WRITE(_double_array)
#define SUBROUTINE_NAME_HYPERSLAB READ_WRITE(_double_hyperslab)
#define TYPE_SUFFIX(x) x ## _double_f
#define H5TYPE H5T_NATIVE_DOUBLE
#endif
#ifdef TYPE_LOGICAL
#define THE_TYPE logical
#define SUBROUTINE_NAME_SCALAR READ_WRITE(_logical)
#define SUBROUTINE_NAME_ARRAY READ_WRITE(_logical_array)
#define SUBROUTINE_NAME_HYPERSLAB READ_WRITE(_logical_hyperslab)
#define TYPE_SUFFIX(x) x ## _int_f
#define H5TYPE H5T_NATIVE_INTEGER
#endif

#ifdef READ
#define STR_READ_WRITE Reads
#define H5LT_READ_WRITE TYPE_SUFFIX(h5ltread_dataset)
#else
#define STR_READ_WRITE Writes
#define H5LT_READ_WRITE TYPE_SUFFIX(h5ltmake_dataset)
#endif
!/* The following routine reads or writes a scalar value, of type THE_TYPE, from the HDF5 file */
!> STR_READ_WRITE a rank-0 THE_TYPE scalar.
subroutine SUBROUTINE_NAME_SCALAR(loc_id, dset_name, buf, errcode)
  integer(HID_T), intent(in) :: loc_id !< HDF5 file id
  character(LEN=*), intent(in) :: dset_name !< HDF5 dataset name
  THE_TYPE, intent(INTENT) :: buf !< data buffer
  integer, intent(out) :: errcode !< HDF5 error code

  integer(HSIZE_T), dimension(1) :: dims_h5type(1)
#ifdef READ
  THE_TYPE, dimension(1) :: buf_1(1)
#else
  THE_TYPE :: buf_(1) !< data buffer
  logical :: exists_
  integer(HID_T) :: dset_id
#endif
#ifdef TYPE_LOGICAL
  integer :: ibuf, ibuf_1(1), ibuf_(1)
#endif

  PUSH_SUB(SUBROUTINE_NAME_SCALAR)

  dims_h5type(1) = 0

#ifdef READ
#ifdef TYPE_LOGICAL
  call H5LT_READ_WRITE(loc_id, dset_name, ibuf_1, dims_h5type, errcode)
  if (errcode .ne. 0) then
     call die("HDF5 error", only_root_writes=.true.)
  endif

  select case (ibuf_1(1))
  case (0)
     buf = .false.
  case (1)
     buf = .true.
  case default
     call die("Invalid value for dset '"//dset_name//"'." , only_root_writes=.true.)
  endselect
#else
  call H5LT_READ_WRITE(loc_id, dset_name, buf_1, dims_h5type, errcode)
  if (errcode .ne. 0) then
     call die("HDF5 error", only_root_writes=.true.)
  endif

  buf = buf_1(1)
#endif
#else
  !FHJ: We can`t use H5LTmake_dataset_* if the dataset already exists!
  call h5lexists_f(loc_id, dset_name, exists_, errcode)
  if (errcode .ne. 0) then
     call die("HDF5 error", only_root_writes=.true.)
  endif

  if (exists_) then
     call h5dopen_f(loc_id, dset_name, dset_id, errcode)
     if (errcode .ne. 0) then
        call die("HDF5 error", only_root_writes=.true.)
     endif

#ifdef TYPE_LOGICAL
     ibuf = 0
     if (buf) ibuf = 1
     call h5dwrite_f(dset_id, H5TYPE, (/ibuf/), dims_h5type, errcode)
     if (errcode .ne. 0) then
        call die("HDF5 error", only_root_writes=.true.)
     endif

#else
     call h5dwrite_f(dset_id, H5TYPE, (/buf/), dims_h5type, errcode)
     if (errcode .ne. 0) then
        call die("HDF5 error", only_root_writes=.true.)
     endif

#endif
     call h5dclose_f(dset_id, errcode)
     if (errcode .ne. 0) then
        call die("HDF5 error", only_root_writes=.true.)
     endif

  else
#ifdef TYPE_LOGICAL
     ibuf = 0
     if (buf) ibuf = 1
     ibuf_(1) = ibuf
     call H5LT_READ_WRITE(loc_id, dset_name, 0, dims_h5type, ibuf_, errcode)

     if (errcode .ne. 0) then
        call die("HDF5 error", only_root_writes=.true.)
     endif

#else
     buf_(1) = buf
     call H5LT_READ_WRITE(loc_id, dset_name, 0, dims_h5type, buf_, errcode)
     if (errcode .ne. 0) then
        call die("HDF5 error", only_root_writes=.true.)
     endif
#endif
  endif
#endif

  POP_SUB(SUBROUTINE_NAME_SCALAR)

end subroutine SUBROUTINE_NAME_SCALAR
!/* The following routine reads or writes an array, of type THE_TYPE, from the HDF5 file */
!> STR_READ_WRITE a(n) THE_TYPE array.
subroutine SUBROUTINE_NAME_ARRAY(loc_id, dset_name, dims, buf, errcode)
  integer(HID_T), intent(in) :: loc_id !< HDF5 file id
  character(LEN=*), intent(in) :: dset_name !< HDF5 dataset name
  integer, intent(in), dimension(:) :: dims !< size of the buffer buf
  THE_TYPE, intent(INTENT), dimension(*) :: buf !< data buffer
  integer, intent(out) :: errcode !< error code

  integer :: rank
  integer(HSIZE_T), dimension(size(dims)) :: dims_h5type
#ifndef READ
  logical :: exists_
  integer(HID_T) :: dset_id
#endif
#ifdef TYPE_LOGICAL
  integer :: ibuf(product(dims))
#ifdef READ
  integer :: ii
#endif
#endif

  PUSH_SUB(SUBROUTINE_NAME_ARRAY)

  rank = size(dims)
  dims_h5type = dims

#ifdef READ
#ifdef TYPE_LOGICAL
  call H5LT_READ_WRITE(loc_id, dset_name, ibuf, dims_h5type, errcode)
  if (errcode .ne. 0) then
     call die("HDF5 error", only_root_writes=.true.)
  endif

  if (any(ibuf<0).or.any(ibuf>1)) then
     call die("Invalid value for dset '"//dset_name//"'." , only_root_writes=.true.)
  endif
  buf(1:product(dims)) = .false.
  do ii=1,product(dims)
     if (ibuf(ii)/=0) buf(ii) = .true.
  enddo
#else
  call H5LT_READ_WRITE(loc_id, dset_name, buf, dims_h5type, errcode)
  if (errcode .ne. 0) then
     call die("HDF5 error", only_root_writes=.true.)
  endif

#endif
#else
  !FHJ: We can`t use H5LTmake_dataset_* if the dataset already exists!
  call h5lexists_f(loc_id, dset_name, exists_, errcode)
  if (errcode .ne. 0) then
     call die("HDF5 error", only_root_writes=.true.)
  endif

  if (exists_) then
     call h5dopen_f(loc_id, dset_name, dset_id, errcode)
     if (errcode .ne. 0) then
        call die("HDF5 error", only_root_writes=.true.)
     endif

#ifdef TYPE_LOGICAL
     ibuf(:) = 0
     where(buf(1:product(dims)))
        ibuf = 1
     endwhere
     call h5dwrite_f(dset_id, H5TYPE, ibuf, dims_h5type, errcode)
     if (errcode .ne. 0) then
        call die("HDF5 error", only_root_writes=.true.)
     endif

#else
     call h5dwrite_f(dset_id, H5TYPE, buf, dims_h5type, errcode)
     if (errcode .ne. 0) then
        call die("HDF5 error", only_root_writes=.true.)
     endif

#endif
     call h5dclose_f(dset_id, errcode)
     if (errcode .ne. 0) then
        call die("HDF5 error", only_root_writes=.true.)
     endif

  else
#ifdef TYPE_LOGICAL
     ibuf(:) = 0
     where(buf(1:dims(1)))
        ibuf = 1
     endwhere
     call H5LT_READ_WRITE(loc_id, dset_name, rank, dims_h5type, ibuf, errcode)
     if (errcode .ne. 0) then
        call die("HDF5 error", only_root_writes=.true.)
     endif

#else
     call H5LT_READ_WRITE(loc_id, dset_name, rank, dims_h5type, buf, errcode)
     if (errcode .ne. 0) then
        call die("HDF5 error", only_root_writes=.true.)
     endif

#endif
  endif
#endif

  POP_SUB(SUBROUTINE_NAME_ARRAY)

end subroutine SUBROUTINE_NAME_ARRAY
!/* The following routine reads or writes a hyperslab from an array, of type THE_TYPE, from the HDF5 file */
!> STR_READ_WRITE a portion of a(n) THE_TYPE array.
subroutine SUBROUTINE_NAME_HYPERSLAB(loc_id, dset_name, countf, offsetf, buf, errcode)
  integer(HID_T), intent(in) :: loc_id !< HDF5 file id
  character(LEN=*), intent(in) :: dset_name !< HDF5 dataset name
  !> Number of elements to read from the dataset for each dimention
  integer, intent(in) :: countf(:)
  !> Offset when reading dataset from file.
  integer, intent(in) :: offsetf(:)
  !> Data buffer. We treat it as a flat contiguous 1D array.
  THE_TYPE, intent(INTENT), dimension(*) :: buf
  integer, intent(out) :: errcode !< error code

  integer(HSIZE_T) :: hcountf(size(countf)) !< Count for file dataspace
  integer(HSIZE_T) :: hcountm(1) !< Count for memory dataspace
  integer(HSIZE_T) :: hoffsetf(size(offsetf)) !< Offset for file dataspace
  integer(HID_T) :: dset_id
  integer(HID_T) :: dataspace
  integer(HID_T) :: memspace
  logical :: exists_
#ifdef TYPE_LOGICAL
  integer :: ibuf(product(countf))
  integer :: ii
#endif

  PUSH_SUB(SUBROUTINE_NAME_HYPERSLAB)

  call h5lexists_f(loc_id, dset_name, exists_, errcode)
  if (errcode .ne. 0) then
     call die("HDF5 1 error", only_root_writes=.true.)
  endif

  if (.not. exists_) call die('Cannot read/write hyperslap from "'//dset_name//&
       "': dataset doesn`t exists!", only_root_writes=.true.)
  call h5dopen_f(loc_id, dset_name, dset_id, errcode)
  if (errcode .ne. 0) then
     call die("HDF5 2 error", only_root_writes=.true.)
  endif

  ! FHJ: Get 2D file dataspace and set selection mask
  call h5dget_space_f(dset_id, dataspace, errcode)
  if (errcode .ne. 0) then
     call die("HDF5 3 error", only_root_writes=.true.)
  endif

  hcountf(:) = countf(:)
  hoffsetf(:) = offsetf(:)
  call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, hoffsetf, hcountf, errcode)
  if (errcode .ne. 0) then
     call die("HDF5 4 error", only_root_writes=.true.)
  endif

  ! FHJ: Create flat memory dataspace
  hcountm(1) = product(countf)
  call h5screate_simple_f(1, hcountm, memspace, errcode)
  if (errcode .ne. 0) then
     call die("HDF5 5 error", only_root_writes=.true.)
  endif

  ! FHJ: Read dataspace
#ifdef TYPE_LOGICAL
#ifdef READ
  call h5dread_f(dset_id, H5TYPE, ibuf, hcountm, errcode, memspace, dataspace)
  if (errcode .ne. 0) then
     call die("HDF5 6 error", only_root_writes=.true.)
  endif

  buf(1:product(countf)) = .false.
  do ii=1,product(countf)
     if (ibuf(ii)/=0) buf(ii) = .true.
  enddo
#else
  ibuf(:) = 0
  do ii=1,product(countf)
     if (buf(ii)) ibuf(ii) = 1
  enddo
  call h5dwrite_f(dset_id, H5TYPE, ibuf, hcountm, errcode, memspace, dataspace)
  if (errcode .ne. 0) then
     call die("HDF5 7 error", only_root_writes=.true.)
  endif

#endif
#else
#ifdef READ
  call h5dread_f(dset_id, H5TYPE, buf, hcountm, errcode, memspace, dataspace)
  if (errcode .ne. 0) then
     call die("HDF5 8 error", only_root_writes=.true.)
  endif

#else
  call h5dwrite_f(dset_id, H5TYPE, buf, hcountm, errcode, memspace, dataspace)
  if (errcode .ne. 0) then
     call die("HDF5 9 error", only_root_writes=.true.)
  endif

#endif
#endif
  call h5sclose_f(memspace, errcode)
  if (errcode .ne. 0) then
     call die("HDF5 10 error", only_root_writes=.true.)
  endif

  call h5sclose_f(dataspace, errcode)
  if (errcode .ne. 0) then
     call die("HDF5 11 error", only_root_writes=.true.)
  endif

  call h5dclose_f(dset_id, errcode)
  if (errcode .ne. 0) then
     call die("HDF5 12 error", only_root_writes=.true.)
  endif

  POP_SUB(SUBROUTINE_NAME_HYPERSLAB)

end subroutine SUBROUTINE_NAME_HYPERSLAB

#undef READ_WRITE
#undef INTENT
#undef THE_TYPE
#undef TYPE_SUFFIX
#undef H5TYPE
#undef SUBROUTINE_NAME_SCALAR
#undef SUBROUTINE_NAME_ARRAY
#undef SUBROUTINE_NAME_HYPERSLAB
#undef STR_READ_WRITE
#undef H5LT_READ_WRITE
#undef H5F_FILE_ACCESS

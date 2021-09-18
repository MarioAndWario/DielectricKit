#include "f_defs.h"

!>===================================================================
!!
!! Module:
!!
!! (1) wfn_io_hdf5_m     Originally by JIM     Last Modified 4/25/2012 (JIM)
!!
!!     Routines to read and write wavefunctions in HDF5 format.
!!     The code is generated through repeated inclusion of a file with
!!     different preprocessor definitions each time. You are not expected to
!!     understand this. Consult the resulting .p.f file for clarity.
!!
!!=================================================================

module wfn_io_hdf5_m
  use global_m
#ifdef HDF5
  use hdf5_io_m
  use hdf5

  implicit none

  private

  public :: &
       read_hdf5_header_type   , write_hdf5_header_type  , &
       read_hdf5_gvectors      , write_hdf5_gvectors     , &
       setup_hdf5_mf_file      , &
       setup_hdf5_wfn_file     , &
       read_hdf5_mf_header     , &
       write_hdf5_mf_header    , &
       read_hdf5_gvec_kpt, read_info, read_kpoints, read_gspace, read_symmetry, read_crystal

contains

  !===============================================================================
#define READ
  !
#define TEMP_HEADER
#include "wfn_io_hdf5_inc.f90"
#undef TEMP_HEADER
  !
#define TEMP_WFN_GVEC
#include "wfn_io_hdf5_inc.f90"
#undef TEMP_WFN_GVEC
  !
#define TEMP_WFN_DATA
#define TEMP_COMPLEX
#include "wfn_io_hdf5_inc.f90"
#undef TEMP_COMPLEX
#include "wfn_io_hdf5_inc.f90"
#undef TEMP_WFN_DATA
  !===============================================================================
#undef READ
  !
#define TEMP_HEADER
#include "wfn_io_hdf5_inc.f90"
#undef TEMP_HEADER
  !
#define TEMP_WFN_GVEC
#include "wfn_io_hdf5_inc.f90"
#undef TEMP_WFN_GVEC
  !
#define TEMP_WFN_DATA
#define TEMP_COMPLEX
#include "wfn_io_hdf5_inc.f90"
#undef TEMP_COMPLEX
#include "wfn_io_hdf5_inc.f90"
#undef TEMP_WFN_DATA
  !
#define TEMP_OTHER
#include "wfn_io_hdf5_inc.f90"
#undef TEMP_OTHER
  !===============================================================================


  !> Create the appropriate structures that hold info about the mean-field
  !! calculation in an HDF5 file. No data is actually written.
  subroutine setup_hdf5_mf_file(fname, create_file)
    character(len=*), intent(in) :: fname
    logical, intent(in), optional :: create_file

    integer(HID_T) :: file_id
    integer(HID_T) :: group_id
    integer :: error
    logical :: create_file_

    PUSH_SUB(setup_hdf5_mf_file)

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

    POP_SUB(setup_hdf5_mf_file)

  end subroutine setup_hdf5_mf_file


  !> Create the appropriate structures that hold info about the WFN
  !! coefficients in an HDF5 file. No data is actually written.
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
    PUSH_SUB(setup_hdf5_wfn_file)

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

    POP_SUB(setup_hdf5_wfn_file)
  end subroutine setup_hdf5_wfn_file

  !> A high-level wrapper for write_*_header* functions
  subroutine write_hdf5_mf_header(fname, mf)
    character(len=*), intent(in) :: fname
    type(mf_header_t), intent(in) :: mf
    character(len=3) :: sheader
    PUSH_SUB(write_hdf5_mf_header)

    sheader = mf%sheader
    call write_hdf5_header_type(fname, sheader, mf%iflavor, mf%kp, mf%gvec, mf%syms, mf%crys)

    POP_SUB(write_hdf5_mf_header)
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
    PUSH_SUB(read_hdf5_mf_header)

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

    POP_SUB(read_hdf5_mf_header)
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

#endif

end module wfn_io_hdf5_m

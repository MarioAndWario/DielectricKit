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
       read_hdf5_header_type   , &
       write_hdf5_header_type  , &
       read_hdf5_gvectors      , &
       write_hdf5_gvectors     , &
                                !#BEGIN_INTERNAL_ONLY
       read_hdf5_wfn_gvectors  , &
       write_hdf5_wfn_gvectors , &
       read_hdf5_band_real     , &
       write_hdf5_band_real    , &
       read_hdf5_band_complex  , &
       write_hdf5_band_complex , &
       read_hdf5_band          , &
       write_hdf5_band         , &
       read_hdf5_bands_block   , &
                                !#END_INTERNAL_ONLY
       setup_hdf5_mf_file      , &
       setup_hdf5_wfn_file     , &
       read_hdf5_mf_header     , &
       write_hdf5_mf_header    , &
       read_hdf5_bands_block_sigma_cplx , read_hdf5_vxc, read_hdf5_gvec_kpt, read_hdf5_wfn_bse, read_hdf5_wfn_complex_hdf2wfn, &
       read_info, read_kpoints, read_gspace, read_symmetry, read_crystal, write_hdf5_band_complex_wfn2hdf

  !#BEGIN_INTERNAL_ONLY
  interface read_hdf5_band
     module procedure read_hdf5_band_real, read_hdf5_band_complex
  end interface read_hdf5_band

  interface write_hdf5_band
     module procedure write_hdf5_band_real, write_hdf5_band_complex
  end interface write_hdf5_band
  !#END_INTERNAL_ONLY

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

  !> Sigma/input.f90:
  ! call read_hdf5_bands_block_sigma_cplx(file_id, kp, peinf%ntband_max, peinf%ntband_node, peinf%does_it_own_inner, ib_first, wfnkqmpi%cg, ioffset=sig%ncore_excl)
  ! ------
  ! nbownmax = peinf%ntband_max
  ! nbownactual = peinf%ntband_node
  ! does_it_ownb = peinf%does_it_own_inner
  ! wfnsout = wfnkqmpi%cg
  ! ioffset = sig%ncore_excl
  ! ------
  subroutine read_hdf5_bands_block_sigma_cplx(file_id, kp, nbownmax, nbownactual, does_it_ownb, ib_first, wfnsout, ioffset)
    integer(HID_T), intent(in) :: file_id
    type(kpoints), intent(in) :: kp
    integer, intent(in) :: nbownmax
    integer, intent(in) :: nbownactual !< how many bands I own
    logical, intent(in) :: does_it_ownb(:,:)
    integer, intent(in) :: ib_first !< first band of the set of bands I own
    complex(DPC), intent(out) :: wfnsout(:,:,:,:) !< (ig,ib,is,ik)
    integer, optional, intent(in) :: ioffset
    real(DP), allocatable :: wfndata(:,:,:,:,:) !< (ic,ig,ib,is,ik)
    integer(HID_T) :: plist_id
    integer(HID_T) :: dset_id
    integer(HID_T) :: dataspace
    integer(HID_T) :: memspace
    integer, parameter :: rank_hdf5_wfn=5
    integer(HSIZE_T) :: count(5), offset(5)
    integer :: error, ipe, reader, nread, ngktot, ioffset_, icount
    integer, allocatable :: ranks(:)
    integer :: mpiworldgrp, mpigrp, bandcomm, ib, ib_, max_bands_read, bands_read, bands_read_max
    integer(kind=8), parameter :: max_bytes_read = 1073741824 ! don`t read/send more than 1 GB at a time
    ! 0=native 1; 1=manual group comm; 2=manual send/recvs
    integer, parameter :: comm_style=0
    logical :: do_read
    character(LEN=30) :: file_wfn
    character(LEN=20) :: ik_string, ip_string
    character(LEN=20) :: ib_string
    integer :: unit_wfn, ib_proc, ib_global, ig, irk, bands_read_old
    call logit('Reading HDF5 WFNs by blocks')

    ioffset_ = 0
    if(present(ioffset)) then
       ioffset_ = ioffset
    endif
    
    ngktot = SUM(kp%ngk)
    nread = peinf%npools
    call h5dopen_f(file_id, '/wfns/coeffs', dset_id, error)
    call h5dget_space_f(dset_id, dataspace, error)

    ! ======
    ! find lowest rank PE that owns the first band of this group of bands
    ! ------
    ! for example, for valence bands, all the procs in a pool owns the same range of vbands
    ! in this case reader should be the lowest global rank of all these procs.
    ! And only the reader-th proc will perform the read operation and broadcast the wavefunctions to the remaining procs in this pool.
    reader = -1
    if (nbownactual > 0) then
       do ipe = 1, peinf%npes
          if (does_it_ownb(ib_first,ipe)) then
             reader = ipe - 1
             exit ! when we find one, we exit, so this is the lowest rank !!
          endif
       enddo
       if (reader==-1) call die("Cannot find reader in read_hdf5_bands_block", only_root_writes=.true.)
    endif

    ! FHJ: read at most max_bands_read bands to avoid 1/1 buffer overflow.
    ! Here, 2*kp%nspin*kp%nspinor*dble(ngktot)*8d0 is the size of a
    ! single band, including G-vectors from all k-points.
    ! max_bands_read < nbownmax = peinf%nvownmax
    ! [ IMPORTANT ]
    max_bands_read = min(nbownmax, int(max_bytes_read/(2*kp%nspin*kp%nspinor*dble(ngktot)*8d0)))
    ! max_bands_read = nbownmax

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

    ib = 1
    bands_read_old = 0
    do while (ib <= nbownmax) ! nbownmax = peinf%nvownmax
       ! bands_read is local
       ! bands_read_max is global, in other word, MAX_ACROSS_PROCS[bands_read] = bands_read_max
       ! There must be at least one procs which actually reads bands_read_max bands in one iteration, which means bands_read = bands_read_max for these procs.
       ! if for current proc, we have bands_read<bands_read_max, it means this is the last reading operation, after this iteration, this LOOP will end.
       bands_read = max(min(nbownactual, ib-1 + max_bands_read) - ib + 1, 0)
       bands_read_max = max(min(nbownmax, ib-1 + max_bands_read) - ib + 1, 0)

       count(1) = 2 ! ic
       count(2) = kp%ngkmax ! ig
       count(3) = bands_read ! ib
       count(4) = kp%nspin*kp%nspinor ! is
       count(5) = kp%nrk ! ik

       ! We will recycle wfndata is the dimension is the same with last one
       if (bands_read > 0) then
          if (bands_read .ne. bands_read_old) then

             if (allocated(wfndata)) then
                deallocate(wfndata)
             endif

             allocate(wfndata (2, kp%ngkmax, bands_read, kp%nspin*kp%nspinor, kp%nrk))
          endif
       else
          if (allocated(wfndata)) then
             deallocate(wfndata)
          endif

          allocate(wfndata (1, 1, 1, 1, 1))
       endif

       bands_read_old = bands_read
       call h5screate_simple_f(rank_hdf5_wfn, count, memspace, error)
       ! <->
       ! if comm_style == 0, all the procs will participate in the reading
       ! ==> do_read always be .true., and all the procs in the same pool will read in the
       ! same wavefunctions
       do_read = (bands_read>0) .and. (peinf%inode==reader.or.comm_style==0)
       if (do_read) then
          offset(1) = 0 ! ic
          offset(2) = 0 ! ig
          offset(3) = (ib_first-1)+ioffset_+(ib-1) ! ib
          offset(4) = 0 ! is
          offset(5) = 0 ! ik
       else
          offset(:) = 0
          call H5sselect_none_f(memspace,error)
       endif

       call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, offset, count, error)

       if (.not.do_read) then
          call H5sselect_none_f(dataspace,error)
       endif

       call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
       call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)       
       call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, wfndata, count, error, memspace, dataspace, xfer_prp=plist_id)
       call h5pclose_f(plist_id, error)
       call h5sclose_f(memspace, error)

       !> for those reader-th procs
       if (do_read) then
          wfnsout(:,ib:ib+bands_read-1,:,:) = DCMPLX(wfndata(1,:,1:bands_read,:,:),wfndata(2,:,1:bands_read,:,:))
       endif
       ib = ib + bands_read_max
    enddo

    if (allocated(wfndata)) then
       deallocate(wfndata)
    endif

    call h5sclose_f(dataspace, error)
    call h5dclose_f(dset_id, error)
  end subroutine read_hdf5_bands_block_sigma_cplx

  subroutine read_hdf5_vxc(sfilename, sig, kp, vxc)
    character(len=*), intent(in) :: sfilename
    type (vxcinfo), intent(inout) :: vxc
    type (kpoints), intent(in) :: kp
    type (siginfo), intent(in) :: sig
    logical :: found
    integer :: nspin, nspinor, nrk, irk, ik_outer, ik_vxc
    integer :: ierror
    integer(hid_t) :: hidfile
    real(dp), pointer :: mtxeld_temp(:,:,:)
    real(dp), pointer :: mtxelo_temp(:,:,:,:)

    if (peinf%inode .eq. 0) then

       call h5fopen_f(sfilename, h5f_acc_rdonly_f, hidfile, ierror)
       call hdf5_read_int(hidfile, '/mf_header/kpoints/nspin', vxc%nspin, ierror)
       call hdf5_read_int(hidfile, '/mf_header/kpoints/nspinor', vxc%nspinor, ierror)
       call hdf5_read_int(hidfile, '/mf_header/kpoints/nrk', vxc%nrk, ierror)

       if (vxc%nrk <= 0) then
          call die('please check the number of kvectors in vxc.h5')
       endif

       if (sig%nkn <= 0) then
          call die('please check the number of kvectors in sigma.inp')
       endif

       allocate(vxc%rk(3,vxc%nrk))
       allocate(vxc%kmap_outer2vxc(sig%nkn))
       allocate(vxc%kmap_vxc2outer(vxc%nrk))

       vxc%kmap_outer2vxc = 0
       vxc%kmap_vxc2outer = 0

       !> read the kvectors
       call hdf5_read_double_array(hidfile, '/mf_header/kpoints/rk', (/3, vxc%nrk/), vxc%rk, ierror)

       do irk = 1, vxc%nrk
          write(*,*) "irk = ", irk, " (", vxc%rk(:,irk) ,")"
       enddo

       !> check the kvectors from 'vxc.h5' with k_outer from 'sig%k'
       !> we have to make sure that eack k_outer can be found in vxc%rk
       write(*,'(a)') "============"
       write(*,'(a)') "checking compatibility between kvectors from outerbands and those from vxc matrix elements ..."

       do ik_outer = 1, sig%nkn
          found = .false.
          do ik_vxc = 1, vxc%nrk
             if(all(abs(sig%kpt(1:3,ik_outer)-vxc%rk(1:3,ik_vxc)) .lt. tol_small)) then
                found = .true.
                !> indice mapping
                vxc%kmap_vxc2outer(ik_vxc) = ik_outer
                vxc%kmap_outer2vxc(ik_outer) = ik_vxc
                exit
             endif
          enddo

          if (found .eqv. .false.) then
             call die('k-point in sigma.inp k_points block not available in vxc.h5 file')
          endif

       enddo

       write(*,'(a)') "done. they are compatible."
       write(*,'(a)') "============"

       !> read vxc matrix elements
       call hdf5_read_int(hidfile, '/vxc/ndiag', vxc%ndiag, ierror)
       call hdf5_read_int(hidfile, '/vxc/diag_nmin', vxc%diag_nmin, ierror)
       call hdf5_read_int(hidfile, '/vxc/diag_nmax', vxc%diag_nmax, ierror)
       call hdf5_read_int(hidfile, '/vxc/noffdiag', vxc%noffdiag, ierror)
       call hdf5_read_int(hidfile, '/vxc/offdiag_nmin', vxc%offdiag_nmin, ierror)
       call hdf5_read_int(hidfile, '/vxc/offdiag_nmax', vxc%offdiag_nmax, ierror)

       if ((vxc%ndiag .eq. 0) .and. (vxc%noffdiag .eq. 0)) then
          call die('ndiag = 0 and noffdiag = 0 in vxc.h5 file')
       endif

       !> check that the ib_outer is contained in [diag_nmin,diag_nmax] range
       ! <- offdiag ->
       if ((vxc%offdiag_nmin .ne. 0) .and. (vxc%offdiag_nmin .ne. 0) ) then
          if ((sig%ib_outer_min_abso_global < vxc%offdiag_nmin) .or. (sig%ib_outer_max_abso_global > vxc%offdiag_nmax)) then
             call die('offdiag : outerbands not available in vxc.h5 file')
          endif
          ! <- diag ->
       else
          if ( (sig%ib_outer_min_abso_global < vxc%diag_nmin) .or. (sig%ib_outer_max_abso_global > vxc%diag_nmax)) then
             call die('diag : outerbands not available in vxc.h5 file')
          endif
       endif

       write(*,'(a)') "============"

       if (vxc%ndiag > 0) then
          write(*,'(a)') "reading diagonal matrix elements from vxc.h5"
          allocate(mtxeld_temp(2,vxc%ndiag,vxc%nrk))
          allocate(vxc%mtxeld(vxc%ndiag,vxc%nrk))
          call hdf5_read_double_array(hidfile, '/vxc/mtxeld', (/2,vxc%ndiag,vxc%nrk/), mtxeld_temp, ierror)
          vxc%mtxeld(:,:) = DCMPLX(mtxeld_temp(1,:,:),mtxeld_temp(2,:,:))
          deallocate(mtxeld_temp)
       endif

       if (vxc%noffdiag > 0) then
          write(*,'(a)') "reading off-diagonal matrix elements from vxc.h5"
          allocate(mtxelo_temp(2,vxc%noffdiag,vxc%noffdiag,vxc%nrk))
          allocate(vxc%mtxelo(vxc%noffdiag,vxc%noffdiag,vxc%nrk))
          call hdf5_read_double_array(hidfile, '/vxc/mtxelo', (/2,vxc%noffdiag,vxc%noffdiag,vxc%nrk/), mtxelo_temp, ierror)
          vxc%mtxelo(:,:,:) = DCMPLX(mtxelo_temp(1,:,:,:),mtxelo_temp(2,:,:,:))
          deallocate(mtxelo_temp)
       endif

       write(*,'(a)') "============"

       call h5fclose_f(hidfile, ierror)

    endif ! endif ionode == 0

    !> broadcast size and range of bands
    if (peinf%npes > 1) then
       call mpi_bcast(vxc%nspin, 1, mpi_integer, 0, mpi_comm_world, mpierr)
       call mpi_bcast(vxc%nspinor, 1, mpi_integer, 0, mpi_comm_world, mpierr)
       call mpi_bcast(vxc%nrk, 1, mpi_integer, 0, mpi_comm_world, mpierr)
       call mpi_bcast(vxc%ndiag, 1, mpi_integer, 0, mpi_comm_world, mpierr)
       call mpi_bcast(vxc%diag_nmin, 1, mpi_integer, 0, mpi_comm_world, mpierr)
       call mpi_bcast(vxc%diag_nmax, 1, mpi_integer, 0, mpi_comm_world, mpierr)
       call mpi_bcast(vxc%noffdiag, 1, mpi_integer, 0, mpi_comm_world, mpierr)
       call mpi_bcast(vxc%offdiag_nmin, 1, mpi_integer, 0, mpi_comm_world, mpierr)
       call mpi_bcast(vxc%offdiag_nmax, 1, mpi_integer, 0, mpi_comm_world, mpierr)
    endif

    !> allocate space on non-root procs
    if (peinf%inode > 0) then
       allocate(vxc%rk(3,vxc%nrk))
       allocate(vxc%kmap_outer2vxc(sig%nkn))
       allocate(vxc%kmap_vxc2outer(vxc%nrk))
       if (vxc%ndiag > 0) then
          allocate(vxc%mtxeld(vxc%ndiag,vxc%nrk))
       endif
       if (vxc%noffdiag > 0) then
          allocate(vxc%mtxelo(vxc%noffdiag,vxc%noffdiag,vxc%nrk))
       endif
    endif

    !> broadcast vxc matrix elements
    if (peinf%npes > 1) then
       call mpi_bcast(vxc%rk,3*vxc%nrk, mpi_double, 0, mpi_comm_world, mpierr)
       call mpi_bcast(vxc%kmap_outer2vxc, sig%nkn, mpi_integer, 0, mpi_comm_world, mpierr)
       call mpi_bcast(vxc%kmap_vxc2outer, vxc%nrk, mpi_integer, 0, mpi_comm_world, mpierr)
       if (vxc%ndiag > 0) then
          call mpi_bcast(vxc%mtxeld,vxc%ndiag*vxc%nrk, mpi_double_complex, 0, mpi_comm_world, mpierr)
       endif
       if (vxc%noffdiag > 0) then
          call mpi_bcast(vxc%mtxelo,vxc%noffdiag*vxc%noffdiag*vxc%nrk, mpi_double_complex, 0, mpi_comm_world, mpierr)
       endif
    endif
  end subroutine read_hdf5_vxc

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

  subroutine read_hdf5_wfn_bse(filename, kp, ib, irk, cg_real, cgarray)
    character(len=*), intent(in) :: filename
    type(kpoints), intent(in) :: kp
    integer, intent(in) :: ib, irk
    real(DP), intent(out) :: cg_real(:,:,:)
    SCALAR, intent(out) :: cgarray(:,:)
    integer(hid_t) :: file_id_wfn
    integer(hid_t) :: dataset_id_wfn
    integer(hid_t) :: dataspace_wfn
    integer(hid_t) :: memspace_wfn
    integer, parameter :: rank_wfn=5
    integer(hsize_t) :: count_wfn(5), offset_wfn(5)
    integer :: error

    !> On ROOT, use hdf5 to read the (ib=ib_global,irk=irk) wavefunction
    if ( peinf%inode .eq. 0 ) then
       call h5fopen_f('WFN_fi.h5', H5F_ACC_RDONLY_F, file_id_wfn,error)
       call h5dopen_f (file_id_wfn, '/wfns/coeffs', dataset_id_wfn, error)
       call h5dget_space_f(dataset_id_wfn, dataspace_wfn, error)
       count_wfn(1) = 2 ! ic
       count_wfn(2) = kp%ngkmax ! ig
       count_wfn(3) = 1 ! ib
       count_wfn(4) = kp%nspin*kp%nspinor ! is
       count_wfn(5) = 1 ! ik

       call h5screate_simple_f(rank_wfn, count_wfn, memspace_wfn,error)
       offset_wfn(1) = 0 ! ic
       offset_wfn(2) = 0 ! ig
       offset_wfn(3) = ib-1 ! ib
       offset_wfn(4) = 0 ! is
       offset_wfn(5) = irk-1 ! ik

       call h5sselect_hyperslab_f(dataspace_wfn, H5S_SELECT_SET_F, offset_wfn, count_wfn, error)

       call h5dread_f(dataset_id_wfn, H5T_NATIVE_DOUBLE, cg_real(1:2,1:kp%ngkmax,1:kp%nspin*kp%nspinor), count_wfn, error, memspace_wfn, dataspace_wfn)
       ! REAL ==> COMPLEX
       cgarray(1:kp%ngk(irk),1:kp%nspin*kp%nspinor) = DCMPLX(cg_real(1,1:kp%ngk(irk),1:kp%nspin*kp%nspinor),cg_real(2,1:kp%ngk(irk),1:kp%nspin*kp%nspinor))

       call h5sclose_f(memspace_wfn, error)
       call h5sclose_f(dataspace_wfn, error)
       call h5dclose_f(dataset_id_wfn, error)
       call h5fclose_f(file_id_wfn, error)
    endif ! ROOT

  end subroutine read_hdf5_wfn_bse

  !> Used in Meanfield/Utilities/hdf2wfn.x
  subroutine read_hdf5_wfn_complex_hdf2wfn(filename, kp, ib, irk, cgarray)
    character(len=*), intent(in) :: filename
    type(kpoints), intent(in) :: kp
    integer, intent(in) :: ib, irk
    SCALAR, intent(out) :: cgarray(:,:)
    real(DP), allocatable :: cg_real(:,:,:)    
    integer(hid_t) :: file_id_wfn
    integer(hid_t) :: dataset_id_wfn
    integer(hid_t) :: dataspace_wfn
    integer(hid_t) :: memspace_wfn
    integer, parameter :: rank_wfn=5
    integer(hsize_t) :: count_wfn(5), offset_wfn(5)
    integer :: error

    !> On ROOT, use hdf5 to read the (ib=ib_global,irk=irk) wavefunction
    if ( peinf%inode .eq. 0 ) then
       call h5fopen_f(filename, H5F_ACC_RDONLY_F, file_id_wfn,error)
       call h5dopen_f (file_id_wfn, '/wfns/coeffs', dataset_id_wfn, error)
       call h5dget_space_f(dataset_id_wfn, dataspace_wfn, error)
       count_wfn(1) = 2 ! ic
       count_wfn(2) = kp%ngkmax ! ig
       count_wfn(3) = 1 ! ib
       count_wfn(4) = kp%nspin*kp%nspinor ! is
       count_wfn(5) = 1 ! ik

       call h5screate_simple_f(rank_wfn, count_wfn, memspace_wfn,error)
       offset_wfn(1) = 0 ! ic
       offset_wfn(2) = 0 ! ig
       offset_wfn(3) = ib-1 ! ib
       offset_wfn(4) = 0 ! is
       offset_wfn(5) = irk-1 ! ik

       call h5sselect_hyperslab_f(dataspace_wfn, H5S_SELECT_SET_F, offset_wfn, count_wfn, error)

       allocate(cg_real(2,kp%ngkmax,kp%nspin*kp%nspinor))
       call h5dread_f(dataset_id_wfn, H5T_NATIVE_DOUBLE, cg_real(1:2,1:kp%ngkmax,1:kp%nspin*kp%nspinor), count_wfn, error, memspace_wfn, dataspace_wfn)
       ! REAL ==> COMPLEX
       cgarray(1:kp%ngk(irk),1:kp%nspin*kp%nspinor) = DCMPLX(cg_real(1,1:kp%ngk(irk),1:kp%nspin*kp%nspinor),cg_real(2,1:kp%ngk(irk),1:kp%nspin*kp%nspinor))

       deallocate(cg_real)
       call h5sclose_f(memspace_wfn, error)
       call h5sclose_f(dataspace_wfn, error)
       call h5dclose_f(dataset_id_wfn, error)
       call h5fclose_f(file_id_wfn, error)
    endif ! ROOT

  end subroutine read_hdf5_wfn_complex_hdf2wfn
  
  subroutine write_hdf5_band_complex_wfn2hdf(fname, wfn, ngk, nstot, ik, ib)
    character(len=*) :: fname
    complex(DPC), intent(in) :: wfn(:,:) !< (ngk,nstot=kp%ns*kp%nspinor)
    integer, intent(in) :: ngk
    integer, intent(in) :: nstot
    integer, intent(in) :: ik
    integer, intent(in) :: ib

    real(DP) :: dwfn(2,ngk,nstot)

    integer(HID_T) :: file_id
    integer(HID_T) :: dataset_id
    integer(HID_T) :: dataspace_id
    integer(HID_T) :: memspace_id
    integer(HSIZE_T) :: offset(5), count(5)
    integer :: error
    integer, parameter :: rank_wfn = 5

    count(1) = 2 ! ic complex = 2
    count(2) = ngk ! ig
    count(3) = 1 ! ib
    count(4) = nstot ! is
    count(5) = 1 ! ik

    offset(1) = 0
    offset(2) = 0
    offset(3) = ib-1
    offset(4) = 0
    offset(5) = ik-1

    dwfn(1,:,:) = real(wfn(:,:))
    dwfn(2,:,:) = aimag(wfn(:,:))

    if(peinf%inode == 0) then
       call h5fopen_f(fname, H5F_ACC_RDWR_F, file_id, error)
       call h5dopen_f(file_id, 'wfns/coeffs', dataset_id, error)
       CALL h5screate_simple_f(rank_wfn, count, memspace_id, error)
       call h5dget_space_f(dataset_id, dataspace_id, error)

       call h5sselect_hyperslab_f(dataspace_id, H5S_SELECT_SET_F, offset, count, error)

       call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, dwfn(1:2, 1:ngk, 1:nstot), count, error, file_space_id = dataspace_id, mem_space_id = memspace_id)

       call h5dclose_f(dataset_id, error)
       call h5sclose_f(dataspace_id, error)
       call h5sclose_f(memspace_id, error)
       call h5fclose_f(file_id, error)
    endif

  end subroutine write_hdf5_band_complex_wfn2hdf

#endif

end module wfn_io_hdf5_m

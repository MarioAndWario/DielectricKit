#include "f_defs.h"

!====================================================================
!
! Routines:
!
! (1) ploteps
!
! input files: epsmat.h5 or chimat.h5, ploteps.inp
!
! output files:
!
! 1. Now each proc reads in the whole eps[0]mat, but if epsinv is too
!    large, we might need to read eps[0]mat within the ifq loop
! 3. Check multiple r2 points
!====================================================================

program ploteps
  use global_m
  use fftw_m
  use fullbz_m
  use inread_m
  use misc_m
  use sort_m
  use write_eps_m
  use ploteps_common_m
  use wfn_io_hdf5_m
  use input_utils_m
  use gmap_m
  use epsread_hdf5_m
  use io_utils_m
  use hdf5
  use h5lt
  use hdf5_io_m
  use scalapack_m
  implicit none

  type (crystal) :: crys
  type (symmetry) :: syms
  type (gspace) :: gvec
  type (kpoints) :: kp
  type (grid) :: qg
  type (ploteps_t) :: peps
  integer :: ii, igrid1_super, igrid2_super, igrid3_super, igrid1_primitive, igrid2_primitive, igrid3_primitive, igrid1_super_downsample, igrid2_super_downsample, igrid3_super_downsample, ir2_loop
  integer :: ncount, ntim
  integer :: ngrid_super(3), ngrid_super_downsample(3), nfft(3)
  real(DP) :: r1(3), scale, tsec(2)
  character :: filename*20
  character*16, allocatable :: routnam(:)
  complex(DPC), allocatable :: ucfft(:,:,:) !< (nfft1, nfft2, nfft3)
  complex(DPC), allocatable :: scfft(:,:,:) !< (ngrid_super_downsample1, ngrid_super_downsample2, ngrid_super_downsample3)
  real(DP) :: r2(3), qr2_, qr1_
  complex(DPC) :: phmqr2

  !> ffeps(ig=1,pol%nmtx, 1: peinf%ikt(peinf%inode+1))
  complex(DPC), allocatable :: ffeps_fq(:)
  character(len=30) :: filename_eps_hdf5 !, filename_eps0_hdf5
  integer :: ifreq, ifq_loc, irq
  logical :: skip_checkbz
  !> integer(HID_T) :: file_id
  integer :: error
  logical :: file_exists
  integer :: ifq_global, irq_global, itran, kg0(3), ig, nmtx, irq_global_last_step
  real(DP) :: rq(3), fq(3)
  real(DP), allocatable :: ekin(:) , freqs_(:,:), q0(:,:), dFreqGrid(:)
  integer, allocatable :: isrtx(:), ind(:), nmtx_(:), nmtx_file(:)
  complex(DPC), allocatable :: ph(:), phr2(:), epsinv(:,:,:), epsinv_rq(:,:), epsinv_fq(:,:), dFreqBrd(:) !, eps0diag(:,:)
  integer :: ig1, ig2, ifreq_target, ifreq_target_, matrix_flavor, matrix_flavor_
  complex(DPC) :: prefactor, temp ! eps0head
  integer :: ng, ng_, nrq0, nfreq_, nfreq_imag_, nmtx_max_file, nmtx_max_file_, qgrid_(3), freq_dep_, total_gvec_grid, total_peps_grid, box_min(3), box_max(3)
  real(DP) :: ecuts_, ecuts, fi, temp_ph, error_max
  integer :: R_vector_frac(3)
  real(DP) :: xi1_frac(3)
  type(progress_info) :: prog_info
  integer :: request, irq_start, irq_end, irq_need, irq_loc, irq_target, nfq_loc_max, irq_eps_start, ipes, ik_loc

  call peinfo_init()
  call h5open_f(error)
  call timacc(0,0)
  call timacc(1,1)

  ! call write_program_header('PlotEps', .false.)
  !> Read ploteps.inp
  call inread(peps)

  !> epsinv in {q G G'} space
  filename_eps_hdf5=TRUNC(peps%filename)

  if (peinf%inode .eq. 0) then
     write(6,'(1X,A)')
     write(6,'(1X,A)') "Reading"//filename_eps_hdf5//" file"
     INQUIRE(FILE=TRUNC(filename_eps_hdf5), EXIST=file_exists)
     if (.not. file_exists) then
        call die("ploteps: "//TRUNC(peps%filename)//" file not exist", only_root_writes=.TRUE.)
     endif

     call read_eps_grid_sizes_hdf5(ng, peps%nrq, ecuts, peps%nfreq, peps%nfreq_imag, nmtx_max_file, peps%qgrid, peps%freq_dep, TRUNC(filename_eps_hdf5))
     peps%nfq = PRODUCT(peps%qgrid(:))
     peps%nsuper(:) = peps%qgrid(:)
     if (peps%nfq .le. 0) then
        call die("Number of qpoints <= 0, please check qgrid", only_root_writes = .true.)
     endif
     !> Check that all the r2 are within the peps%nsuper(3) range
     if ( (any(peps%r2(1,:) > peps%nsuper(1))) .or. (any(peps%r2(2,:) > peps%nsuper(2))) .or. (any(peps%r2(3,:) > peps%nsuper(3))) ) then
        call die("r2 beyond range", only_root_writes=.TRUE.)
     elseif ( (any(peps%r2(1,:) < 0)) .or. (any(peps%r2(2,:) < 0)) .or. (any(peps%r2(3,:) < 0)) ) then
        call die("r2 beyond range", only_root_writes=.TRUE.)
     endif
     if ( ANY( (peps%nsuper(:) - peps%qgrid(:)) > 0 ) ) then
        call die("nsuper cannot exceed qgrid", only_root_writes = .TRUE.)
     endif
     if ( ANY( (peps%nsuper(:) .le. 0 ) ) ) then
        call die("Please check supercell_size", only_root_writes = .TRUE.)
     endif
     if ( ANY( (peps%qgrid(:) .le. 0 ) ) ) then
        call die("Please check qgrid", only_root_writes = .TRUE.)
     endif

     write(6,'(1X,A,3(i4,1x))') '- Supercell size:', peps%nsuper
     write(6,'(1X,A,3(i4,1x),A,I4,A)') '- Qgrid :', peps%qgrid(:), " with ", peps%nfq, " qpoints in full BZ"
     write(*,'(1X,A,F10.3,A)') "- G-cutoff of epsinv read from "//TRUNC(filename_eps_hdf5)//" : ", ecuts, " Ry"

     if ((peps%freq_dep .eq. 0) .and. (peps%nfreq .gt. 1)) then
        call die("ploteps: GPP epsmat with > 1 frequencies.", only_root_writes=.true.)
     endif
     !> if FF but with 1 frequency
     if ((peps%freq_dep .ne. 0) .and. (peps%nfreq .le. 1)) then
        call die("ploteps: FF epsmat with <= 1 frequencies.", only_root_writes=.true.)
     endif
     if ((peps%freq_dep .ne. 0) .and. (SCALARSIZE .eq. 1)) then
        call die("ploteps: FF not compatible with complex-flavor.", only_root_writes=.true.)
     endif

     if (peps%ecut > ecuts) then
        write(*,'(A,F10.3,A,F10.3)') "peps%ecut = ", peps%ecut, " ecuts = ", ecuts
        call die("ploteps: peps%ecut > ecuts.", only_root_writes=.true.)
     endif

     if (peps%freq_dep .ne. 0) then
        SAFE_ALLOCATE(dFreqGrid, (peps%nfreq))
        SAFE_ALLOCATE(dFreqBrd,  (peps%nfreq))
        call read_eps_freqgrid_hdf5(peps%nfreq, dFreqGrid, dFreqBrd, TRUNC(filename_eps_hdf5))
        ifreq_target = 0
        do ifreq = 1, peps%nfreq
           if (ABS(dFreqGrid(ifreq) + dFreqBrd(ifreq) - peps%freq_target) < TOL_SMALL) then
              ifreq_target = ifreq
              exit
           endif
        enddo
        if (ifreq_target .eq. 0) then
           call die("epscopy: ifreq_zero not found.")
        endif
        SAFE_DEALLOCATE(dFreqGrid)
        SAFE_DEALLOCATE(dFreqBrd)
        !> if GPP epsmat, ifreq_target = 1
     else
        if (ABS(peps%freq_target) > TOL_ZERO) then
           call die("GPP only has zero frequency.", only_root_writes=.true.)
        endif
        ifreq_target = 1
     endif

     write(6,'(1X,A,F12.5,",",F12.5,A)') "- Target frequency : (", peps%freq_target, ")."
     write(6,'(1X,A,I5)') "- ifreq_target in "//TRUNC(filename_eps_hdf5)//" = ", ifreq_target

     !> GPP + real-flavor
     matrix_flavor = SCALARSIZE
     !> Check matrix_flavor in epsmat
     call read_eps_matrix_flavor_hdf5(matrix_flavor_, TRUNC(filename_eps_hdf5))
     if (matrix_flavor .ne. matrix_flavor_) then
        call die("ploteps: matrix_flavor mismatch betwen code and epsmat.")
     endif

     call read_gspace(TRUNC(filename_eps_hdf5), gvec)
     call read_symmetry(TRUNC(filename_eps_hdf5), syms)
     call read_crystal(TRUNC(filename_eps_hdf5), crys)

  endif

  if (peinf%npes > 1) then
     call MPI_BCAST(ifreq_target, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
     call MPI_BCAST(peps%nfq, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
     call MPI_BCAST(peps%qgrid, 3, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
     call MPI_BCAST(peps%nsuper, 3, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)

     call MPI_BCAST(peps%nrq, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
     call MPI_BCAST(nmtx_max_file, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
     call MPI_BCAST(peps%freq_dep, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
     call MPI_BCAST(peps%nfreq, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
     call MPI_BCAST(peps%nfreq_imag, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)

     call MPI_BCAST(gvec%ng, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
     call MPI_BCAST(gvec%ecutrho, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
     call MPI_BCAST(gvec%FFTgrid, 3, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)

     call MPI_BCAST(crys%nat, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
     call MPI_BCAST(crys%celvol, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
     call MPI_BCAST(crys%alat, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
     call MPI_BCAST(crys%avec, 9, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
     call MPI_BCAST(crys%adot, 9, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
     call MPI_BCAST(crys%recvol, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
     call MPI_BCAST(crys%blat, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
     call MPI_BCAST(crys%bvec, 9, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
     call MPI_BCAST(crys%bdot, 9, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)

     call MPI_BCAST(syms%ntran, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
     call MPI_BCAST(syms%cell_symmetry, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
     call MPI_BCAST(syms%mtrx(1,1,1), 3*3*48, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
     call MPI_BCAST(syms%mtrx_reci(1,1,1), 3*3*48, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
     call MPI_BCAST(syms%mtrx_cart(1,1,1), 3*3*48, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
     call MPI_BCAST(syms%tnp, 3*48, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
  endif

  SAFE_ALLOCATE(gvec%components, (3, gvec%ng))
  call read_hdf5_gvectors(TRUNC(filename_eps_hdf5), gvec%ng, gvec%components)
  !> Common/input_utils.f90:
  !> Compute index_vec indices relating G-vectors in fractional coordinates to positions in the FFT grid
  !> The address of gvec%components(:,i) is given by
  !> address=((gx_i+gxmax)*(2*gymax+1)+gy_i+gymax)*(2*gzmax+1)+gz_i+gzmax+1
  !> gvec%index_vec(iadd) = ig
  call gvec_index(gvec)

  SAFE_ALLOCATE(peps%rq, (3,peps%nrq))
  SAFE_ALLOCATE(peps%nmtx, (peps%nrq))
  SAFE_ALLOCATE(peps%isrtxi, (gvec%ng, peps%nrq))
  peps%isrtxi = 0
  SAFE_ALLOCATE(isrtx, (gvec%ng))
  SAFE_ALLOCATE(ekin,  (gvec%ng))

  if (peinf%inode .eq. 0) then
     ! write(6,'(A)') "Reading data from eps[0]mat.h5 files"
     ! write(6,'(A)')

     SAFE_ALLOCATE(nmtx_file, (peps%nrq))

     !> Read header of epsmat.h5 and check input parameters
     call read_eps_qgrid_hdf5(peps%nrq, peps%rq, nmtx_file, TRUNC(filename_eps_hdf5))
     if (NORM2(peps%rq(:,1)) > TOL_SMALL) then
        call die("ploteps: first rq must be Gamma.", only_root_writes=.true.)
     endif

     !> Note that nmtx_max_file is for epsinv_rq read from epsmat.h5
     !> while peps%nmtx is the number of Gvectors with peps%ecuts
     if (nmtx_max_file .ne. MAXVAL(nmtx_file(:))) then
        call die("nmtx_max_file .ne. MAXVAL(nmtx_file(:))",only_root_writes=.true.)
     endif

     ! write(6,'(1X,A,I5,A)') "We will consider ", peps%nrq," RBZ qpoints"
     ! do irq = 1, peps%nrq
     !    write(6,'(1X , "(", 3F12.5 , ")", A, I5 )') peps%rq(:, irq), " nmtx_file = ", nmtx_file(irq)
     ! enddo

     !> Use peps%ecut to determine a new peps%nmtx, instead of using that from epsmat.h5
     !> Try to reduce the size of nfft using get_eps_fftgrid(...)
     box_max(:) = 0
     box_min(:) = 0
     do irq = 1, peps%nrq
        call kinetic_energies(gvec, crys%bdot, ekin, qvec = peps%rq(:, irq))
        call sortrx(gvec%ng, ekin, isrtx, gvec = gvec%components)
        peps%nmtx(irq) = gcutoff(gvec%ng, ekin, isrtx, peps%ecut)
        if (peps%nmtx(irq) > nmtx_file(irq)) then
           call die("ploteps: nmtx(irq) > nmtx_file(irq).", only_root_writes=.true.)
        endif
        do ig = 1, gvec%ng
           peps%isrtxi(isrtx(ig), irq) = ig
        enddo
        do ig = 1, peps%nmtx(irq)
           box_min(1:3) = MIN(box_min(1:3), gvec%components(1:3, isrtx(ig)))
           box_max(1:3) = MAX(box_max(1:3), gvec%components(1:3, isrtx(ig)))
        enddo
     enddo
     peps%nmtx_max = MAXVAL(peps%nmtx(:))
     SAFE_DEALLOCATE(nmtx_file)

     if (peps%high_resolution) then
        peps%FFTgrid(:) = gvec%FFTgrid(:)
     else
        peps%FFTgrid(1:3) = peps%FFTfactor * ( box_max(1:3) - box_min(1:3) + 1 )
        !> Use gvec%FFTgrid to correct pol%WFN_FFTgrid, to make sure that the ratio betwen gvec%FFTgrid(:) is preseved
        total_gvec_grid = SUM(gvec%FFTgrid(:))
        total_peps_grid = SUM(peps%FFTgrid(:))
        peps%FFTgrid(1) = NINT(DBLE(gvec%FFTgrid(1)) / DBLE(Total_gvec_grid) * DBLE(total_peps_grid))
        peps%FFTgrid(2) = NINT(DBLE(gvec%FFTgrid(2)) / DBLE(Total_gvec_grid) * DBLE(total_peps_grid))
        peps%FFTgrid(3) = NINT(DBLE(gvec%FFTgrid(3)) / DBLE(Total_gvec_grid) * DBLE(total_peps_grid))
     endif
  endif

  if (peinf%npes > 1) then
     call MPI_BCAST(      peps%rq,       3*peps%nrq, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
     call MPI_BCAST(    peps%nmtx,         peps%nrq,          MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
     call MPI_BCAST(peps%nmtx_max,                1,          MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
     call MPI_BCAST(  peps%isrtxi, gvec%ng*peps%nrq,          MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
     call MPI_BCAST( peps%FFTgrid,                3,          MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
  endif

  qg%nr = peps%nrq
  SAFE_ALLOCATE(qg%r, (3,qg%nr))
  qg%r(:,:) = peps%rq(:,:)
  !> Use fullbz to generate from peps%rq all the qpoints in FBZ
  if (peps%unfold) then
     call fullbz_2(crys, syms, qg)
  else
     call fullbz_2(crys, syms, qg, use_identity_only=.true.)
  endif

  !> Block-distribute all fq points over peinf%npes procs
  !> Blocksize
  peinf%nkpe = iceil(peps%nfq, peinf%npes)
  !> (global) peinf%ik(iproc, ifq_loc) = ifq_global
  SAFE_ALLOCATE(peinf%ik, (peinf%npes, peinf%nkpe))
  !> (global) peinf%ikt(iproc) = # of ifq_global dealt by each proc
  SAFE_ALLOCATE(peinf%ikt, (peinf%npes))
  peinf%ik  = 0
  peinf%ikt = 0
  nfq_loc_max = NUMROC(peps%nfq, peinf%nkpe, 0, 0, peinf%npes)
  if (nfq_loc_max .ne. peinf%nkpe) then
     call die("nfq_loc_max .ne. nkpe.", only_root_writes=.true.)
  endif

  do ipes = 1, peinf%npes
     peinf%ikt(ipes) = NUMROC(peps%nfq, peinf%nkpe, ipes-1, 0, peinf%npes)
     do ik_loc = 1, peinf%ikt(ipes)
        peinf%ik(ipes, ik_loc) = INDXL2G(ik_loc, peinf%nkpe, ipes-1, 0, peinf%npes)
     enddo
  enddo

  if (peinf%inode .eq. 0) then
     write(6,'(1X,A)') "Distribute fq points:"
     write(6,'(1X,A,I0)') "- Maximal number of full BZ qpoint on a proc: ", nfq_loc_max
     write(6,'(1X,A,I0)') "- Number of idle procs : ", COUNT(peinf%ikt(:) == 0)
     write(6,'(A)')
  endif

  call MPI_barrier(MPI_COMM_WORLD,mpierr)

  !> low_comm
  !> ROOT reads in epsinv matrix from epsmat.h5 and broadcast to all procs
  if (peps%low_comm) then
     peps%rq_blocksize = peps%nrq
     peps%nrq_loc = peps%nrq
     peps%nrq_loc_max = peps%nrq
     SAFE_ALLOCATE(epsinv, (nmtx_max_file, nmtx_max_file, peps%nrq))
     if (peinf%inode .eq. 0) then
        call read_eps_matrix_ser_allq_hdf5(epsinv, nmtx_max_file, peps%nrq, 1, ifreq_target, TRUNC(filename_eps_hdf5))
        !> Correct the head using eps0head
        epsinv(1,1,1) = DCMPLX(peps%epsinvhead)
     endif
     if (peinf%npes > 1) then
        call MPI_Bcast(epsinv, nmtx_max_file*nmtx_max_file*peps%nrq, MPI_COMPLEX_DPC, 0, MPI_COMM_WORLD, mpierr)
     endif
  else
     SAFE_ALLOCATE(peinf%irq, (peinf%npes, nfq_loc_max))
     DO ipes = 1, peinf%npes
        DO ifq_loc = 1, peinf%ikt(ipes)
           peinf%irq(ipes, ifq_loc) = qg%indr(peinf%ik(ipes, ifq_loc))
        ENDDO
     ENDDO
     !> mid_comm
     !> Distribute nrq
     peps%rq_blocksize = MAX(CEILING(DBLE(peps%nrq) / DBLE(peinf%npes)), 1)
     peps%nrq_loc      = NUMROC(peps%nrq, peps%rq_blocksize, peinf%inode, 0, peinf%npes)
     peps%nrq_loc_max  = NUMROC(peps%nrq, peps%rq_blocksize,           0, 0, peinf%npes)
     SAFE_ALLOCATE(epsinv, (nmtx_max_file, nmtx_max_file, peps%nrq_loc_max))
     if (peinf%inode .eq. 0) then
        write(6,'(1X,A,I5)') "Maximal number of rq points stored in a proc : ", peps%nrq_loc_max
     endif
     !> Parallel HDF5 read
     if (peps%nrq_loc > 0) then
        irq_eps_start = INDXL2G(1, peps%rq_blocksize, peinf%inode, 0, peinf%npes)
     else
        irq_eps_start = 0
     endif
     !> we can use MPI_Gather to assemble an epsinv with all G2, for one iq.
     call read_eps_matrix_par_distribute_rq_hdf5(epsinv, nmtx_max_file, 1, irq_eps_start, peps%nrq_loc, ifreq_target, TRUNC(filename_eps_hdf5))
     if (peinf%inode .eq. 0) then
        ! epsinv(1,1,1) = eps0head
        epsinv(1,1,1) = DCMPLX(peps%epsinvhead)
     endif
  endif

  SAFE_ALLOCATE(epsinv_rq, (nmtx_max_file, nmtx_max_file))
  !> prefactor = 1/(N \Omega), crys%celvol in units of Bohr^3, prefactor in units of Ang^(-3)

  prefactor = 1.0D0 / (DBLE(qg%nf) * crys%celvol * BOHR**3)
  call setup_FFT_sizes(peps%FFTgrid, nfft, scale)
  !> The number of real-space points within a primitive cell is PRODUCT(nfft)
  SAFE_ALLOCATE(ucfft, (nfft(1), nfft(2), nfft(3)))
  !> ngrid_super defines the real-space mesh in the BvO supercell
  ngrid_super(:) = nfft(:) * peps%nsuper(:)
  ngrid_super_downsample(:) = (ngrid_super(:) + peps%downsample(:) - 1) / peps%downsample(:) !> CEIL[ngrid_super/downsample]
  SAFE_ALLOCATE(scfft, (ngrid_super_downsample(1), ngrid_super_downsample(2), ngrid_super_downsample(3)))
  !> The number of real-space points within the large cell is PRODUCT(ngrid_super_downsample)
  if (peinf%inode .eq. 0) then
     write(6,'(1X,a)') 'Grids information:'
     write(6,'(1X,a,3(1x,i0))') '- FFT box size:', nfft
     write(6,'(1X,a,3(1x,i0))') '- Supercell grid:', ngrid_super
     write(6,'(1X,a,3(1x,i0))') '- Effective supercell grid:', ngrid_super_downsample
  endif

  !> Loop over r2 positions within a unit cell, r2 positions read from ploteps.inp
  DO ir2_loop = 1, peps%nr2
     !> r2 in fractional coordinates
     r2(:) = peps%r2(:, ir2_loop)

     if (peinf%inode .eq. 0) then
        write(6,'(1X,A,3ES30.23,A)') "r2 = (", r2(:), ")"
     endif

     scfft = ZERO

     DO ifq_loc = 1, nfq_loc_max
        ! call progress_step(prog_info, ifq_loc)
        if (ifq_loc <= peinf%ikt(peinf%inode+1)) then
           ifq_global = peinf%ik(peinf%inode+1, ifq_loc)
           !> Find corresponding iq_RBZ
           irq_global = qg%indr(ifq_global)
           rq(:)  = qg%r(:, irq_global)
           fq(:)  = qg%f(:, ifq_global)
           itran  = qg%itran(ifq_global)
           kg0(:) = qg%kg0(:, ifq_global)
        else
           ifq_global = 0
           irq_global = 0
           rq = 0.0D0
           fq = 0.0D0
           itran = 0
           kg0 = 0
        endif

        if (peps%low_comm) then
           !> low_comm
           if (ifq_loc <= peinf%ikt(peinf%inode+1)) then
              epsinv_rq(:, :) = epsinv(:, :, irq_global)
           else
              epsinv_rq = ZERO
           endif
        else
           !> mid_comm
           if (peps%nrq_loc > 0) then
              irq_start = INDXL2G(1,              peps%rq_blocksize, peinf%inode, 0, peinf%npes)
              irq_end   = INDXL2G(peps%nrq_loc, peps%rq_blocksize, peinf%inode, 0, peinf%npes)
           else
              irq_start = 0
              irq_end = 0
           endif

           !> NON-BLOCKING RECEIVE
           !> current proc is not idle
           if (ifq_loc <= peinf%ikt(peinf%inode+1)) then
              !> if the rq needed by current proc is elsewhere, MPI_IRECV
              !> if the rq needed by current proc is just here, copy epsinv to epsinv_loc_rq
              irq_need = peinf%irq(peinf%inode+1, ifq_loc)
              !> Current proc has the irq_need
              if ((irq_need <= irq_end) .and. (irq_need >= irq_start)) then
                 irq_loc = INDXG2L(irq_need, peps%rq_blocksize, peinf%inode, 0, peinf%npes)
                 epsinv_rq(:, :) = epsinv(:, :, irq_loc)
              else
                 call MPI_IRECV( epsinv_rq(1,1), nmtx_max_file*nmtx_max_file, MPI_SCALAR, MPI_ANY_SOURCE, irq_need, MPI_COMM_WORLD, request, mpierr)
              endif
           endif

           !> loop over all procs, see if current proc needs to send out epscol to a target proc
           !> All procs take care of the ipes-th procs
           do ipes = 1, peinf%npes
              !> ipes-th proc is not idle
              if (peinf%ik(ipes, ifq_loc) .ne. 0) then
                 irq_target = peinf%irq(ipes, ifq_loc)
                 !> See if ipes-th proc needs epscol from current procs
                 if ((irq_target <= irq_end) .and. (irq_target >= irq_start)) then
                    irq_loc = INDXG2L(irq_target, peps%rq_blocksize, 0, 0, peinf%npes)
                    if (ipes .ne. (peinf%inode+1)) then
                       call MPI_SEND(epsinv(1,1,irq_loc), nmtx_max_file*nmtx_max_file, MPI_SCALAR, ipes-1, irq_target, MPI_COMM_WORLD, mpierr)
                    endif
                 endif
              endif
           enddo
           !> Once current procs receive the epsinv_loc_rq it needs, it is free to move on
           if (ifq_loc <= peinf%ikt(peinf%inode+1)) then
              !> Current proc does not have the irq_need
              if ( .not. ((irq_need <= irq_end) .and. (irq_need >= irq_start)) ) then
                 call MPI_WAIT(request, MPI_STATUS_IGNORE, mpierr)
              endif
           endif
        endif

        if (ifq_loc <= peinf%ikt(peinf%inode+1)) then
           call kinetic_energies(gvec, crys%bdot, ekin, qvec = fq)
           call sortrx(gvec%ng, ekin, isrtx, gvec = gvec%components)
           nmtx = gcutoff(gvec%ng, ekin, isrtx, peps%ecut)

           !> Calculate phr2(:) = e^{-i ( q + G ) \dot r2} by rewriting part of gmap(...)
           !> Here G is order by |fq+G|^2
           SAFE_ALLOCATE(phr2, (nmtx))
           do ig = 1, nmtx
              fi = 2.0D0 * PI_D * DOT_PRODUCT((DBLE(gvec%components(1:3, isrtx(ig))) + fq(1:3)), r2(1:3))
              phr2(ig) = DCMPLX(cos(fi), -sin(fi))
           enddo

           SAFE_ALLOCATE(ind, (nmtx))
           SAFE_ALLOCATE(ph, (nmtx))
           !> MATMUL(syms%mtrx_reci(1:3,1:3,itran),rq(1:3)) + kg0(1:3) = fq(1:3)
           !> isorti: map from old file (rq) to gvec%components
           !> isort: map from gvec%components to new file (fq)
           call gmap_2(gvec, syms, nmtx, itran, kg0, isrtx, peps%isrtxi(:, irq_global), ind, ph)

           !> Unfold epsinv_rq to epsinv_fq, with extra CONJG(phr2(ig1)) = e^{i G \cdot r2} phase
           SAFE_ALLOCATE(epsinv_fq, (nmtx, nmtx))
           !$OMP PARALLEL DO collapse(2)
           do ig1 = 1, nmtx
              do ig2 = 1, nmtx
                 epsinv_fq(ig1, ig2) = ph(ig1) * CONJG(ph(ig2)) * epsinv_rq(ind(ig1),ind(ig2))
              enddo
           enddo
           !$OMP END PARALLEL DO
           SAFE_DEALLOCATE(ind)
           SAFE_DEALLOCATE(ph)
           
           !> Take zgemv of epsinv_fq(ig1,ig2) and phr2(ig2) to calculate ffeps_fq(ig1)
           !> SUBROUTINE ZGEMV(TRANS,M,N,ALPHA,A,LDA,x,INCX,BETA,Y,INCY)
           SAFE_ALLOCATE(ffeps_fq, (nmtx))
           call zgemv('n', nmtx, nmtx, ONE, epsinv_fq(:,:), nmtx, phr2(:), 1, ZERO, ffeps_fq(:), 1)

           SAFE_DEALLOCATE(phr2)
           SAFE_DEALLOCATE(epsinv_fq)

           !> Take FFT of f(iq_FBZ) to calculate ucfft (FFT grid = gvec%FFTgrid)
           ucfft = ZERO
           !> gvec_index
           call put_into_fftbox(nmtx, ffeps_fq(:), gvec%components, isrtx, ucfft, nfft)

           SAFE_DEALLOCATE(ffeps_fq)
           call do_FFT(ucfft, nfft, 1)

           !> r_frac = [      0,       0,       0] ==> ucfft(1,1,1)
           !> r_frac = [  delta,   delta,   delta] ==> ucfft(2,2,2)
           !> r_frac = [1-delta, 1-delta, 1-delta] ==> ucfft(Nfft(1),Nfft(2),Nfft(3))

           !> Add e^{iq.r1} phase to ucfft and add it to scfft, of which the size is ngrid_super_downsample(1:3)
           !> [Important] OMP does not work with loops with steps (e.g., do i = 1, N, m)
           !> !$OMP PARALLEL DO collapse(3) private(igrid1_super_downsample,igrid2_super_downsample,igrid3_super_downsample,igrid1_primitive,igrid2_primitive,igrid3_primitive,r1,qr1_) reduction(+:scfft)
           do igrid3_super = 1, ngrid_super(3), peps%downsample(3)
              do igrid2_super = 1, ngrid_super(2), peps%downsample(2)
                 do igrid1_super = 1, ngrid_super(1), peps%downsample(1)
                    !> igrid_super_downsample from [1, 1, 1] --> [ngrid_super_downsample(1), ngrid_super_downsample(2), ngrid_super_downsample(3)]
                    !> igrid_super_downsample = CEIL(igrid_super, peps%downsample)
                    igrid3_super_downsample = ((igrid3_super-1) / peps%downsample(3)) + 1
                    igrid2_super_downsample = ((igrid2_super-1) / peps%downsample(2)) + 1
                    igrid1_super_downsample = ((igrid1_super-1) / peps%downsample(1)) + 1

                    !> igrid_primitive from [1, 1, 1] --> [nfft(1), nfft(2), nfft(3)]
                    igrid3_primitive = MOD(igrid3_super-1, nfft(3)) + 1 !> the corresponding FFT grid index in the central primitive cell
                    igrid2_primitive = MOD(igrid2_super-1, nfft(2)) + 1
                    igrid1_primitive = MOD(igrid1_super-1, nfft(1)) + 1

                    R_vector_frac(:) = (/ (igrid1_super-1)/nfft(1), (igrid2_super-1)/nfft(2), (igrid3_super-1)/nfft(3) /)
                    xi1_frac(:) = (/ DBLE(igrid1_primitive-1)/DBLE(nfft(1)), DBLE(igrid2_primitive-1)/DBLE(nfft(2)), DBLE(igrid3_primitive-1)/DBLE(nfft(3))  /)

                    !> r1 ==> fractional coordinates (can be > 1) in the BvO cell
                    !> r1(:) within a primitive cell from [0, 0, 0] --> [1-delta, 1-delta, 1-delta]
                    r1(3) = DBLE(igrid3_super-1)/DBLE(nfft(3))
                    r1(2) = DBLE(igrid2_super-1)/DBLE(nfft(2))
                    r1(1) = DBLE(igrid1_super-1)/DBLE(nfft(1))

                    qr1_ = 2.0D0 * PI_D * DOT_PRODUCT(fq, (xi1_frac(:) + DBLE(R_vector_frac(:))))
                    scfft(igrid1_super_downsample, igrid2_super_downsample, igrid3_super_downsample) = scfft(igrid1_super_downsample, igrid2_super_downsample, igrid3_super_downsample) + prefactor * DCMPLX(cos(qr1_),sin(qr1_)) * ucfft(igrid1_primitive, igrid2_primitive, igrid3_primitive)

                 enddo ! igrid1_super
              enddo ! igrid2_super
           enddo ! igrid3_super
           !> !$OMP END PARALLEL DO
        endif !> if not idle
     ENDDO ! ifq_loc

     if (peinf%npes > 1) then
        !> Sum over all the procs (all fq)
        if (peinf%inode .eq. 0) then
           call MPI_REDUCE(MPI_IN_PLACE, scfft(1,1,1), ngrid_super_downsample(1)*ngrid_super_downsample(2)*ngrid_super_downsample(3), MPI_COMPLEX_DPC, MPI_SUM, 0, MPI_COMM_WORLD, mpierr)
        else
           call MPI_REDUCE(scfft(1,1,1), scfft(1,1,1), ngrid_super_downsample(1)*ngrid_super_downsample(2)*ngrid_super_downsample(3), MPI_COMPLEX_DPC, MPI_SUM, 0, MPI_COMM_WORLD, mpierr)
        endif
     endif

     !> ROOT writes out epsinv for ir2_loop
     if (peinf%inode.eq.0) then
        call write_eps(peps, crys, nfft, scfft, ir2_loop)
     endif
  ENDDO ! ir2_loop
  call progress_free(prog_info)

  call destroy_fftw_plans()

  if (peps%mid_comm) then
     SAFE_DEALLOCATE_P(peinf%irq)
  endif
  SAFE_DEALLOCATE(epsinv_rq)
  SAFE_DEALLOCATE(epsinv)
  SAFE_DEALLOCATE_P(gvec%components)
  SAFE_DEALLOCATE_P(gvec%index_vec)
  SAFE_DEALLOCATE_P(peps%rq)
  SAFE_DEALLOCATE_P(peps%nmtx)
  SAFE_DEALLOCATE(ekin)
  SAFE_DEALLOCATE(isrtx)
  SAFE_DEALLOCATE(ucfft)
  SAFE_DEALLOCATE(scfft)

  call dealloc_grid_simple(qg)

  !> Time accounting
  ntim = 1
  SAFE_ALLOCATE(routnam, (ntim))
  routnam(1) = 'TOTAL:'
  call timacc(1,2)
  if (peinf%inode .eq. 0) then
     write(6,*)
     write(6,9000) 'CPU (s)','WALL (s)','#'
     write(6,*)
     do ii=2,ntim
        call timacc(ii,3,tsec,ncount)
        write(6,9001) routnam(ii),tsec(1),tsec(2),ncount
     enddo
     call timacc(1,3,tsec,ncount)
     write(6,9004) routnam(1),tsec(1),tsec(2)
     write(6,*)
9000 format(22x,a13,  3x,a13,  3x,a8)
9001 format(1x,a16,'      ',f13.3,3x,f13.3,3x,i8)
9004 format(1x,a16,'      ',f13.3,3x,f13.3)
  endif
  call h5close_f(error)
#ifdef MPI
  call MPI_FINALIZE(mpierr)
#endif

end program ploteps

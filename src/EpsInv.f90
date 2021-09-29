#include "f_defs.h"

!=====================================================================
!
! Program
!
!   EpsInv      by Meng Wu (2020)
!
! DESCRIPTION:
!
!   EpsInv.x calculates inverse dielectric response function $\epsilon^{-1}_{G_1 G_2}(q, \omega)$
!   with input of irreducible polarizability $\chi^{\star}_{G_1 G_2}(q, \omega)$
!   by inverting $\delta_{G_1 G_2} - v_{G_1}(q) chi^{\star}_{G_1 G_2}(q, \omega)$
!
!   We calculate each q-point sequentially. For each q-point, We parallel over
!   G1-vectors and G2-vectors using SCALAPACK 2D block-cyclic layout.
!
! USAGE:
!
!   1. For q1 on a uniform grid: EpsInv.x chimat.h5
!   2. For q0 on a shifted grid: EpsInv.x chi0mat.h5
!
!=====================================================================

program EpsInv
  use global_m
  use symmetries_m
  use blas_m
  use misc_m
  use input_utils_m
  use wfn_io_hdf5_m
  use sort_m
  use hdf5
  use h5lt
  use hdf5_io_m
  use epswrite_hdf5_m
  use epsread_hdf5_m
  use write_matrix_m
  use vcoul_generator_m
  use inversion_m
  use scalapack_m
  use lapack_m
  implicit none

  type (crystal) :: crys
  type (polarizability) :: pol
  type (symmetry) :: syms
  type (gspace) :: gvec
  type (kpoints) :: kp
  type (scalapack) :: scal

  real(DP) :: tsec(2),tmin(2),tmax(2)
  character*16, allocatable :: routnam(:)
  integer, allocatable :: routsrt(:)
  integer :: error, info, ii, ncount
  real(DP), allocatable :: ekin(:)

  !! HDF5 identifiers
  integer(HID_T) :: file_id
  integer(HID_T) :: dset_id
  integer(HID_T) :: plist_id
  integer(HID_T) :: group_id
  integer(HID_T) :: filespace
  integer(HID_T) :: memspace
  integer :: rank_mat = 6, rank_matdiag = 4
  integer(HSIZE_T) :: count_mat(6), offset_mat(6), count_matdiag(4), offset_matdiag(4)

  integer :: iq, ig
  logical :: file_exists
  character(len=3) :: sheader
  character(len=50) :: filename_chi_hdf5, filename_eps_hdf5
  integer :: iflavor, ng_, nfq, matrix_flavor_
  integer :: nprow_, iprow, npr_max, npc_max, nroutnam
  logical :: is_q1, qpt_done
  !! dmat_1d_block(icomplex, ig1, ig2, ifreq, is, iq)
  !! dmat_1d_block_diag(icomplex, ig1=ig2, ifreq, iq)
  real(DP), allocatable :: dmat_1d_block(:,:,:,:,:), dmat_1d_block_diag(:,:,:)
  complex(DPC), allocatable :: zmat_1d_block(:,:)
  real(DP) :: fact
  integer :: desc_1d(9), desc_2d(9)
  integer :: cntxt_1d, info_blacs
  integer :: block_size_col, nprow, npcol, myprow, mypcol
  integer :: npr, npc, lld_1d, lld_2d
  integer :: ig2_offset, is, ifreq, ig1_loc, ig2_loc, ig1, ig2
  integer, parameter :: ctxt_ = 2
  real(DP), allocatable :: vcoul(:), vcoul_temp(:)
  SCALAR, allocatable :: eps(:,:)
  complex(DPC), allocatable :: epsRDyn(:,:,:)
  complex(DPC), allocatable :: eps_temp(:,:)
  integer :: target_myprow, target_mypcol, ipes, ig1_gvec, ig2_gvec, irow, icol
  integer :: ntran_, cell_symmetry_
  real(DP) :: tnp_(3,48)
  integer :: mtrx_(3,3,48), nfft_(3)=(/0,0,0/), celltype
  integer :: spacegroup_, i
  character :: symbol_*21, temp_char*20
  real(DP), allocatable :: apos_frac(:,:)
  real(DP) :: q_(3), qlen2_, ekinx_min

  info = MPI_INFO_NULL
  call peinfo_init()
  call h5open_f(error)
  call timacc(0,0)
  call timacc(1,1)

  filename_chi_hdf5 = "chimat.h5"
  if (peinf%inode .eq. 0) then
     if ((command_argument_count() .eq. 0)) then
        write(6,'(1x,a)') "usage: EpsInv.x chi[0]mat.h5"
        call h5close_f(error)
#ifdef MPI
        call MPI_FINALIZE(mpierr)
#endif
        stop
     endif
     call get_command_argument(1, filename_chi_hdf5)

     if (TRIM(filename_chi_hdf5) .eq. "chimat.h5") then
        is_q1 = .true.
     elseif (TRIM(filename_chi_hdf5) .eq. "chi0mat.h5") then
        is_q1 = .false.
     else
        call die(TRUNC(filename_chi_hdf5)//" is not supported", only_root_writes=.true.)
     endif

     INQUIRE(FILE=TRUNC(filename_chi_hdf5), EXIST=file_exists)
     if (.not. file_exists) then
        call die("Cannot find "//TRUNC(filename_chi_hdf5), only_root_writes=.true.)
     endif
  endif

  if (peinf%npes > 1) then
     call MPI_BCAST( filename_chi_hdf5, 50, MPI_CHARACTER, 0, MPI_COMM_WORLD, mpierr)
     call MPI_BCAST( is_q1, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, mpierr)
  endif

  call timacc(2,1)

  sheader = 'CHI'
  iflavor = 0
  call read_hdf5_header_type(TRUNC(filename_chi_hdf5), sheader, iflavor, kp, gvec, syms, crys)
  SAFE_ALLOCATE(apos_frac, (3, crys%nat))
  apos_frac = MATMUL(TRANSPOSE(crys%bvec), crys%apos)

  !! Use SPGlib to get the crystal symmetry
  call get_symmetries(crys%nat, crys%atyp, apos_frac, crys%avec, nfft_, cell_symmetry_, ntran_, mtrx_, tnp_, spacegroup_, symbol_)

  if (peinf%inode .eq. 0) then
     write(6,'(A)') "Atoms:"
     do i = 1, crys%nat
        write(6,'(A,I5,A,I5,A,3F12.9,A)') "Atom #", i, " Z = ", crys%atyp(i), " f = (", apos_frac(:,i), " )"
     enddo
     write(6,*) 'Lattice vectors = [a1, a2, a3] :'
     do i = 1, 3
        write(6,'(F15.8,1X,F15.8,1X,F15.8,1X)') crys%avec(i,:)
     enddo
     write(6,*) 'Reciprocal vectors = [b1, b2, b3] :'
     do i = 1, 3
        write(6,'(F15.8,1X,F15.8,1X,F15.8,1X)') crys%bvec(i,:)
     enddo
     write(6,*) 'bdot = [b1, b2, b3]^T \cdot [b1, b2, b3] :'
     do i = 1, 3
        write(6,'(F15.8,1X,F15.8,1X,F15.8,1X)') crys%bdot(i,:)
     enddo

     write(6,'(A,I10)') "cell_symmetry_ = ", cell_symmetry_
     write(6,'(A,I10)') "ntran_ = ", ntran_
     write(6,'(A,I10)') "spacegroup_ = ", spacegroup_
     write(6,'(A,A)') "symbol_ = ", symbol_
     write(6,'(A)') '-----------------------'
  endif
  SAFE_DEALLOCATE(apos_frac)

  !! https://homepage.univie.ac.at/michael.leitner/lattice/pearson/ctype.html#cftype
  !! FCC
  if ((spacegroup_ .eq. 216) .or. (spacegroup_ .eq. 225) .or. (spacegroup_ .eq. 227)) then
     celltype = 1
     !! BCC
  elseif (spacegroup_ .eq. 229) then
     celltype = 2
     !! simple cubic
  elseif (spacegroup_ .eq. 221) then
     celltype = 3
     !! Hexagonal + triangular cell
  elseif ((spacegroup_ .ge. 143) .and. (spacegroup_ .le. 194)) then
     celltype = 4
     !! Other
  else
     celltype = 0
  endif

  if (peinf%inode .eq. 0) then
     if (celltype .eq. 1) then
        write(6,'(1X,A)') "FCC structure"
     elseif (celltype .eq. 2) then
        write(6,'(1X,A)') "BCC structure"
     elseif (celltype .eq. 3) then
        write(6,'(1X,A)') "Simple cubic structure"
     elseif (celltype .eq. 4) then
        write(6,'(1X,A)') "Triangular structure"
     else
        write(6,'(1X,A)') "Other structure"
     endif
  endif

  pol%timeordered = .true.

  if (peinf%inode .eq. 0) then
     call read_eps_grid_sizes_hdf5(ng_, pol%nq, pol%ecuts, pol%nfreq, pol%nfreq_imag, pol%nmtx, pol%qgrid, pol%freq_dep, TRUNC(filename_chi_hdf5))
     if (is_q1) then
        pol%nq0 = 0
        pol%nq1 = pol%nq
        write(6,'(1X,A)') 'For q1 vectors:'
     else
        pol%nq0 = pol%nq
        pol%nq1 = 0
        write(6,'(1X,A)') 'For q0 vectors:'
     endif
     nfq = PRODUCT(pol%qgrid(:))
     if (nfq .le. 0) then
        call die("Number of qpoints <= 0, please check qgrid", only_root_writes = .true.)
     endif
     if ( ANY( (pol%qgrid(:) .le. 0 ) ) ) then
        call die("Please check qgrid", only_root_writes = .TRUE.)
     endif

     write(6,'(1X,A)') "Reciprocal-space information:"
     write(6,'(1X,A,3(i4,1x),A,I4,A)') '- Q-grid :', pol%qgrid(:), " with ", nfq, " qpoints in full BZ."
     write(6,'(1X,A,F10.3,A)') "- G-cutoff of eps read from "//TRUNC(filename_chi_hdf5)//" : ", pol%ecuts, " Ryd."

     if ((pol%freq_dep .eq. 0) .and. (pol%nfreq .gt. 1)) then
        call die("EpsInv: GPP epsmat with > 1 frequencies.", only_root_writes=.true.)
     endif

     if ((pol%freq_dep .ne. 0) .and. (pol%nfreq .le. 1)) then
        call die("EpsInv: FF epsmat with <= 1 frequencies.", only_root_writes=.true.)
     endif

     if ((pol%freq_dep .ne. 0) .and. (SCALARSIZE .eq. 1)) then
        call die("EpsInv: FF not compatible with complex-flavor.", only_root_writes=.true.)
     endif

     SAFE_ALLOCATE(pol%dFreqGrid, (pol%nfreq))
     SAFE_ALLOCATE(pol%dFreqBrd,  (pol%nfreq))
     call read_eps_freqgrid_hdf5(pol%nfreq, pol%dFreqGrid, pol%dFreqBrd, TRUNC(filename_chi_hdf5))
     pol%nfreq_real = pol%nfreq - pol%nfreq_imag
     !! GPP + real-flavor
     pol%matrix_flavor = SCALARSIZE
     !! Check matrix_flavor in epsmat
     call read_eps_matrix_flavor_hdf5(matrix_flavor_, TRUNC(filename_chi_hdf5))
     if (pol%matrix_flavor .ne. matrix_flavor_) then
        call die("EpsInv: matrix_flavor mismatch betwen code and epsmat.")
     endif

     call read_eps_params_hdf5(TRUNC(filename_chi_hdf5), pol)
     if (pol%timeordered) then
        write(6,'(1X,A)') "Time-ordered <=> varepsilon is non-zero for real-axis frequencies"
        if (pol%nfreq_real .gt. 0) then
           write(6,'(1X,A,I10)') "Real-axis frequencies: nfreq_real = ", pol%nfreq_real
           do ifreq = 1, pol%nfreq_real
              !! pol%dFreqBrd(ifreq) is frequency-dependent varepsilon in this case
              write(6,'(1X,A,I5,A,F12.5,A,F12.5,A)') "#", ifreq, " omega = (", pol%dFreqGrid(ifreq), " ) eV varepsilon = ", DIMAG(pol%dFreqBrd(ifreq)), " eV"
           enddo
        endif

        if (pol%nfreq_imag .gt. 0) then
           write(6,'(1X,A,I10)') "Imaginary-axis frequencies: nfreq_imag = ", pol%nfreq_imag
           do ifreq = pol%nfreq_real+1, pol%nfreq
              write(6,'(1X,A,I5,A,"(",F12.5, " + i ", F12.5,")",A)') "#", ifreq, " omega = ", pol%dFreqGrid(ifreq), DIMAG(pol%dFreqBrd(ifreq)), " eV varepsilon = 0 eV"
           enddo
        endif
     else
        write(6,'(1X,A)') "Retarded <==> varepsilon = 0 eV"
        if (pol%nfreq_real .gt. 0) then
           write(6,'(1X,A,I10)') "Real-axis frequencies: nfreq_real = ", pol%nfreq_real
           do ifreq = 1, pol%nfreq_real
              write(6,'(1X,A,I5,A,"(",F12.5, " + i ", F12.5,")",A)') "#", ifreq, " omega = ", pol%dFreqGrid(ifreq), DIMAG(pol%dFreqBrd(ifreq)), " eV"
           enddo
        endif

        if (pol%nfreq_imag .gt. 0) then
           write(6,'(1X,A,I10)') "Imaginary-axis frequencies: nfreq_imag = ", pol%nfreq_imag
           do ifreq = pol%nfreq_real+1, pol%nfreq
              write(6,'(1X,A,I5,A,"(",F12.5, " + i ", F12.5,")",A)') "#", ifreq, " omega = ", pol%dFreqGrid(ifreq), DIMAG(pol%dFreqBrd(ifreq)), " eV"
           enddo
        endif
     endif
     write(6,'(1X,A)')

     !! Coulomb truncation schemes
     if (peinf%inode .eq. 0) then
        select case (pol%icutv)
        case (0)
           write(6,'(1X,A)') 'We are using no truncation'
        case (2)
           call die("Spherical truncation not supported.", only_root_writes=.true.)
        case (4)
           write(6,'(1X,A)') 'We are using a truncated Coulomb interaction: Cell Wire'
        case (5)
           write(6,'(1X,A)') 'We are using a truncated Coulomb interaction: Cell Box'
        case (6)
           write(6,'(1X,A)') 'We are using a truncated Coulomb interaction: Cell Slab'
        case (7)
           write(6,'(1X,A)') 'We are using a truncated Coulomb interaction: Supercell Box'
        end select
     endif

     if (pol%matrix_type .ne. 2) then
        call die("You should read a chimat.h5 file containing polarizability.", only_root_writes=.true.)
     endif
     !! pol%matrix_type = 0 means we will output a epsmat file
     pol%matrix_type = 0
     ! write(6,'(/1X,A)') 'More job parameters:'
     write(6,'(1X,A,I0)') '- Number of valence bands: ', pol%nvb
     ! write(6,'(1X,A,I0)') '- Number of valence bands to skip: ', pol%skip_nvb
     write(6,'(1X,A,I0)') '- Number of conduction bands: ', pol%ncb
     ! write(6,'(1X,A,I0)') '- Number of conduction bands to skip: ', pol%skip_ncb
     write(6,'(1X,A,I0)') '- Number of spins: ', kp%nspin
     ! write(6,'(1X,A,3I10)') "- FFTgrid = ", pol%FFTgrid
     write(6,'(1X,A,I0)') '- Max. number of G1/G2 vectors for one q: ', pol%nmtx
     write(6,'()')
  endif

  !! BCAST pol to all procs
  if (peinf%npes > 1) then
     call MPI_BCAST(pol%freq_dep, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
     call MPI_BCAST(pol%nfreq, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
     call MPI_BCAST(pol%nfreq_imag, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
     call MPI_BCAST(pol%nq, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
     call MPI_BCAST(pol%nq0, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
     call MPI_BCAST(pol%nq1, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
     call MPI_BCAST(pol%nband, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
     call MPI_BCAST(pol%nvb, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
     call MPI_BCAST(pol%ncb, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
     call MPI_BCAST(pol%nmtx, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
     ! call MPI_BCAST(pol%skip_nvb, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
     ! call MPI_BCAST(pol%skip_ncb, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
     call MPI_BCAST(pol%matrix_flavor, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
     call MPI_BCAST(pol%matrix_type, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
     call MPI_BCAST(pol%nmatrix, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
     call MPI_BCAST(pol%icutv, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
     call MPI_BCAST(pol%intraband_flag, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
     call MPI_BCAST(pol%intraband_overlap_min, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
     call MPI_BCAST(pol%has_advanced, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, mpierr)
     call MPI_BCAST(pol%subsample, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, mpierr)
     call MPI_BCAST(pol%subspace, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, mpierr)
     call MPI_BCAST(pol%timeordered, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, mpierr)
     ! call MPI_BCAST(pol%FFTgrid, 3, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
     call MPI_BCAST(pol%qgrid, 3, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
     call MPI_BCAST(pol%ecuts, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
     call MPI_BCAST(pol%efermi, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
     if (peinf%inode .ne. 0) then
        SAFE_ALLOCATE(pol%dFreqGrid, (pol%nfreq))
        SAFE_ALLOCATE(pol%dFreqBrd,  (pol%nfreq))
     endif
     call MPI_BCAST(pol%dFreqGrid, pol%nfreq, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
     call MPI_BCAST(pol%dFreqBrd, pol%nfreq, MPI_COMPLEX_DPC, 0, MPI_COMM_WORLD, mpierr)
  endif

  !! Read Gvectors
  SAFE_ALLOCATE(gvec%components, (3, gvec%ng))
  call read_hdf5_gvectors(TRUNC(filename_chi_hdf5), gvec%ng, gvec%components)
  call gvec_index(gvec)

  nfq = PRODUCT(pol%qgrid(:))
  SAFE_ALLOCATE(pol%qpt, (3, pol%nq))
  SAFE_ALLOCATE(pol%nmtx_of_q, (pol%nq))
  call read_eps_qgrid_hdf5(pol%nq, pol%qpt, pol%nmtx_of_q, TRUNC(filename_chi_hdf5))

  !! Setup eps[0]mat.h5 file
  if (is_q1) then
     filename_eps_hdf5 = "epsmat.h5"
  else
     filename_eps_hdf5 = "eps0mat.h5"
  endif

  call MPI_BARRIER(MPI_COMM_WORLD, mpierr)

  if (peinf%inode .eq. 0) then
     call eps_hdf5_setup(kp, gvec, syms, crys, pol, TRUNC(filename_eps_hdf5), restart = pol%restart)
  endif
  
  !! Scalapack setup
  !! 2D block-cyclic distribution of chimat
  !! Find a grid such that nprow * npcol = peinf%npes
  nprow_ = CEILING(SQRT(DBLE(peinf%npes) + 1.0D-6))
  do iprow = nprow_, 1, -1
     if(mod(peinf%npes, iprow) .eq. 0) exit
  enddo
  scal%nprow = iprow
  scal%npcol = peinf%npes / scal%nprow
  !! All procs will be in the proc grid!
  !! But there could be idle proc with no data to work on.

  !! Blocksize
  !! Here pol%nmtx = MAXVAL(pol%nmtx_of_q)
  scal%nbr = MIN(64, pol%nmtx / MIN(scal%nprow, scal%npcol))
  scal%nbc = MIN(64, pol%nmtx / MIN(scal%nprow, scal%npcol))

  scal%myprow = -1
  scal%mypcol = -1
  call blacs_get(-1, 0, scal%icntxt)
  !! proc grid is row-major
  !! The data grid is column-major.
  call blacs_gridinit(scal%icntxt, 'r', scal%nprow, scal%npcol)
  call blacs_gridinfo(scal%icntxt, scal%nprow, scal%npcol, scal%myprow, scal%mypcol)

  !! Number of G1 owned by this proc, could be 0
  scal%npr = NUMROC(pol%nmtx, scal%nbr, scal%myprow, 0, scal%nprow)
  !! Number of G2 owned by this proc, could be 0
  scal%npc = NUMROC(pol%nmtx, scal%nbc, scal%mypcol, 0, scal%npcol)
  !! Number of (G1,G2) pair owned by this proc
  peinf%myown = scal%npr * scal%npc
  !! The 1st proc must have the heaviest work load
  npr_max = NUMROC(pol%nmtx, scal%nbr, 0, 0, scal%nprow)
  npc_max = NUMROC(pol%nmtx, scal%nbc, 0, 0, scal%npcol)
  peinf%nckpe = npr_max * npc_max
  if (peinf%inode .eq. 0) then
     write(6,'(A)')
     write(6,'(1X,A,I8,A,I8)') "Process grid : ", scal%nprow ," x ", scal%npcol
     write(6,'(1X,A,I8,A,I8)') "Blocksize of (ik, ikp) : ", scal%nbr, " x ", scal%nbc
     write(6,'(1X,A,I8,A,I8,A,I8)') "Max. number of (ik, ikp) on a proc : ", npr_max, " x ", npc_max, " = ", peinf%nckpe
     write(6,'(A)')
  endif

  fact = 4.0D0 / (DBLE(nfq) * crys%celvol * DBLE(kp%nspin*kp%nspinor))

  !! Initialize BLACS grid for 2d block-cyclic distributed matrix pol%chi(ig1_loc,ig2_loc,is)
  lld_2d = MAX(scal%npr,1)
  if ( scal%myprow .ne. -1) then
     call DESCINIT(desc_2d, pol%nmtx, pol%nmtx, scal%nbr, scal%nbc, 0, 0, scal%icntxt, lld_2d, info_blacs)
  else
     desc_2d(ctxt_) = -1
     info_blacs = 0
  endif
  if (info_blacs .ne. 0) then
     call die('DESCINIT failed for pol%chi', only_root_writes=.true.)
  endif

  !! Initialize BLACS grid for 1d block distributed column matrix dmat_1d_block(ig1,ig2_loc,ifreq,is)
  block_size_col = ICEIL(pol%nmtx, peinf%npes)

  call blacs_get(0, 0, cntxt_1d)
  call blacs_gridinit(cntxt_1d, 'R', 1, peinf%npes)
  call blacs_gridinfo(cntxt_1d, nprow, npcol, myprow, mypcol)
  if ( (myprow .ge. 1) .or. (mypcol .ge. peinf%npes) ) then
     myprow = -1;
     mypcol = -1;
  endif
  if (myprow .ne. -1) then
     npr = NUMROC(pol%nmtx, pol%nmtx, myprow, 0, nprow)
     npc = NUMROC(pol%nmtx, block_size_col, mypcol, 0, npcol)
  else
     npr = 0
     npc = 0
  endif
  lld_1d = MAX(1, npr)
  if (npr .ne. pol%nmtx) then
     call die("Distribution failed 2.", only_root_writes=.true.)
  endif
  !! This is for evecs_r
  if ( myprow .ne. -1) then
     if (lld_1d .ne. pol%nmtx) then
        call die("lld_id should be pol%nmtx.", only_root_writes=.true.)
     endif
     call DESCINIT(desc_1d, pol%nmtx, pol%nmtx, pol%nmtx, block_size_col, 0, 0, cntxt_1d, lld_1d, info_blacs)
  else
     desc_1d(ctxt_) = -1
     info_blacs = 0
  endif
  if (info_blacs .ne. 0) then
     call die('DESCINIT failed for dmat_1d_block', only_root_writes=.true.)
  endif

  SAFE_ALLOCATE(zmat_1d_block, (pol%nmtx, MAX(npc, 1)))
  SAFE_ALLOCATE(pol%isrtx, (gvec%ng))
  SAFE_ALLOCATE(pol%isrtxi, (gvec%ng))
  SAFE_ALLOCATE(ekin,  (gvec%ng))
  SAFE_ALLOCATE(vcoul, (pol%nmtx))

  call timacc(2,2)

  !! Main loop over all q-points
  iq_loop: do iq = 1, pol%nq
     if (peinf%inode .eq. 0) then
        write(6,'(1X,A)') "======================================================================"
        write(6,'(1X,I5,A,I5,5X,A,3F20.10)') iq, " /", pol%nq, " q = ", pol%qpt(:, iq)
        write(6,'(1X,A)') "======================================================================"
        write(6,'(1X,A)')
        write(6,'(1X,A,F15.8,A)') "|q| = ", SQRT(DOT_PRODUCT(pol%qpt(:,iq), MATMUL(crys%bdot, pol%qpt(:,iq)))), " Bohr^{-1}"
     endif

     !! |q+G|^2
     if (.not. is_q1) then
        if (peinf%inode .eq. 0) then
           write(6,'(1X,A)') 'This is the special q->0 point.'
           write(6,'(1X,A)')
        endif
        filename_eps_hdf5 = 'eps0mat.h5'
        call kinetic_energies(gvec, crys%bdot, ekin)
     else
        if (peinf%inode .eq. 0) then
           write(6,'(1X,A)') 'This is a regular non-zero q-point.'
           write(6,'(1X,A)')
        endif
        filename_eps_hdf5 = 'epsmat.h5'
        call kinetic_energies(gvec, crys%bdot, ekin, qvec = pol%qpt(:, iq))
     endif

     if (pol%restart) then
        if (peinf%inode .eq. 0) qpt_done = is_qpt_done(TRUNC(filename_eps_hdf5), iq)
#ifdef MPI
        if (peinf%npes > 1) then
           call MPI_BCAST(qpt_done, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, mpierr)
        endif
#endif
        if (qpt_done) then
           if (peinf%inode == 0) then
              write(6,'(/,1X,A,/)') 'This q-point was already calculated: skipping.'
              write(6,'(1X,A)')
           endif
           cycle iq_loop
        endif
     endif

     call sortrx(gvec%ng, ekin, pol%isrtx, gvec = gvec%components)
     pol%isrtxi = 0
     do ig = 1, gvec%ng
        pol%isrtxi(pol%isrtx(ig)) = ig
     enddo

     if (peinf%inode .eq. 0) then
        call write_gvec_indices_hdf(gvec%ng, pol%isrtx, pol%isrtxi, ekin, iq, TRUNC(filename_eps_hdf5))
        write(6, '(1X,A,I0)') 'Rank of the polarizability matrix: ', pol%nmtx_of_q(iq)
        write(6,'(1X,A)')
     endif

     call timacc(3,1)

     if (peinf%inode .eq. 0) then
        write(6,'(1X,A)') 'Reading '//TRUNC(filename_chi_hdf5)
     endif

     !! PHDF5 reads chimat.h5 for one q, all frequencies into 1D block-cyclic distributed dmat_1d_block(icomplex,ig1,ig2_loc,iomega,is)
     !! use pzgemr2d to load dmat_1d_block to pol%chi(:, :, is)
     ig2_offset = INDXL2G( 1, block_size_col, peinf%inode, 0, peinf%npes)
     SAFE_ALLOCATE(dmat_1d_block, (SCALARSIZE, pol%nmtx, MAX(npc, 1), pol%nfreq, kp%nspin))
     dmat_1d_block = 0.0D0

     !! Open filename_chi_hdf5 for collective read
     call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
     call h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, info, error)
     call h5fopen_f(TRUNC(filename_chi_hdf5), H5F_ACC_RDONLY_F, file_id, error, access_prp = plist_id)
     call h5pclose_f(plist_id, error)
     CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
     CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)

     !! pol%nmtx = MAXVAL(pol%nmtx_of_q) here
     if (npc > 0) then
        count_mat(:) = (/ SCALARSIZE, pol%nmtx, npc, pol%nfreq, kp%nspin, 1 /)
        offset_mat(:) = (/ 0, 0, ig2_offset - 1, 0, 0, iq - 1 /)
     else
        count_mat = 0
        offset_mat = 0
     endif

     call h5screate_simple_f(rank_mat, count_mat, memspace, error)
     if (error .ne. 0) then
        call die("HDF5 error", only_root_writes=.true.)
     endif

     !! Select hyperslab in the file
     call h5dopen_f(file_id, 'mats/matrix', dset_id, error)
     if (error .ne. 0) then
        call die("HDF5 error", only_root_writes=.true.)
     endif

     call h5dget_space_f(dset_id, filespace, error)
     if (error .ne. 0) then
        call die("HDF5 error", only_root_writes=.true.)
     endif

     call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset_mat, count_mat, error)
     if (error .ne. 0) then
        call die("HDF5 error", only_root_writes=.true.)
     endif

     if ( npc <= 0) then
        call H5sselect_none_f(filespace, error)
        if (error .ne. 0) then
           call die("HDF5 error", only_root_writes=.true.)
        endif
        call H5sselect_none_f(memspace, error)
        if (error .ne. 0) then
           call die("HDF5 error", only_root_writes=.true.)
        endif
     endif

     call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, dmat_1d_block, count_mat, error, mem_space_id = memspace, file_space_id = filespace, xfer_prp=plist_id)
     if (error .ne. 0) then
        call die("HDF5 error", only_root_writes=.true.)
     endif

     call h5sclose_f(memspace, error)
     call h5sclose_f(filespace, error)
     call h5dclose_f(dset_id, error)
     call h5pclose_f(plist_id, error)
     call h5fclose_f(file_id, error)

     !! Assemble pol%chi or pol%chiRDyn
     if (pol%freq_dep .eq. 0) then
        SAFE_ALLOCATE(pol%chi, (MAX(scal%npr,1), MAX(scal%npc,1), kp%nspin))
        pol%chi = ZERO
     else
        SAFE_ALLOCATE(pol%chiRDyn, (MAX(scal%npr,1), MAX(scal%npc,1), pol%nfreq, kp%nspin))
        pol%chiRDyn = (0.0D0, 0.0D0)
     endif

     do is = 1, kp%nspin
        do ifreq = 1, pol%nfreq
           zmat_1d_block = ZERO
           !! Trim extra G1 and G2 vector outside pol%nmtx_of_q(iq)
           g2_loop: do ig2_loc = 1, npc
              ig2 = INDXL2G(ig2_loc, block_size_col, mypcol, 0, npcol)
              if (ig2 > pol%nmtx_of_q(iq)) then
                 cycle g2_loop
              endif
              do ig1 = 1, pol%nmtx_of_q(iq)
                 zmat_1d_block(ig1, ig2_loc) = DCMPLX(dmat_1d_block(1, ig1, ig2_loc, ifreq, is), dmat_1d_block(2, ig1, ig2_loc, ifreq, is))
              enddo
           enddo g2_loop

           if (pol%freq_dep .eq. 0) then
              call pX(gemr2d)(pol%nmtx, pol%nmtx, zmat_1d_block, 1, 1, desc_1d, pol%chi(:, :, is), 1, 1, desc_2d, cntxt_1d)
           else
              call pX(gemr2d)(pol%nmtx, pol%nmtx, zmat_1d_block, 1, 1, desc_1d, pol%chiRDyn(:, :, ifreq, is), 1, 1, desc_2d, cntxt_1d)
           endif
        enddo
     enddo
     SAFE_DEALLOCATE(dmat_1d_block)
     call timacc(3,2)

     !! Calculate vcoul
     call timacc(4,1)
     vcoul = 0.0D0

     if (peinf%inode .eq. 0) then
        write(6,'(1X,A)') "call vcoul_generator_zerovq0(...)"
     endif
     call vcoul_generator(pol%icutv, gvec, crys, pol%nmtx, pol%isrtx, pol%qpt(:, iq), vcoul)
     call timacc(4,2)

     if (pol%freq_dep .eq. 0) then
        SAFE_ALLOCATE(eps, (MAX(scal%npr,1), MAX(scal%npc,1)))
        eps = ZERO
     else
        SAFE_ALLOCATE(epsRDyn, (MAX(scal%npr,1), MAX(scal%npc,1), pol%nfreq))
        SAFE_ALLOCATE(eps_temp, (MAX(scal%npr,1), MAX(scal%npc,1)))
        epsRDyn = (0.0D0, 0.0D0)
        eps_temp = (0.0D0, 0.0D0)
     endif

     call timacc(5,1)

     !! Construct eps

     !! Generalized plasmon-pole (GPP)
     if (pol%freq_dep .eq. 0) then
        !$OMP PARALLEL DO collapse(2) private(target_myprow, target_mypcol, ig1_loc, ig2_loc)
        do ig2 = 1, pol%nmtx
           do ig1 = 1, pol%nmtx
              target_myprow = INDXG2P(ig1, scal%nbr, scal%myprow, 0, scal%nprow)
              target_mypcol = INDXG2P(ig2, scal%nbc, scal%mypcol, 0, scal%npcol)
              if ( (target_myprow .ne. scal%myprow) .or. (target_mypcol .ne. scal%mypcol) ) cycle
              ig1_loc = INDXG2L(ig1, scal%nbr, scal%myprow, 0, scal%nprow)
              ig2_loc = INDXG2L(ig2, scal%nbc, scal%mypcol, 0, scal%npcol)
              if (ig1 .eq. ig2) then
                 eps(ig1_loc, ig2_loc) = ONE - vcoul(ig1) * pol%chi(ig1_loc, ig2_loc, 1)
              else
                 eps(ig1_loc, ig2_loc) =     - vcoul(ig1) * pol%chi(ig1_loc, ig2_loc, 1)
              endif
           enddo
        enddo
        !$OMP END PARALLEL DO
        !! Full-frequency (FF)
     else
        SAFE_ALLOCATE(vcoul_temp, (scal%npr))
        eps_temp = (0.0D0, 0.0D0)
        vcoul_temp = 0.0D0
        do ig1 = 1, pol%nmtx
           ig2 = ig1
           !! Here the third argument is dummy, not used by the function
           target_myprow = INDXG2P(ig1, scal%nbr, scal%myprow, 0, scal%nprow)
           target_mypcol = INDXG2P(ig2, scal%nbc, scal%mypcol, 0, scal%npcol)
           ig1_loc = INDXG2L(ig1, scal%nbr, scal%myprow, 0, scal%nprow)
           ig2_loc = INDXG2L(ig2, scal%nbc, scal%mypcol, 0, scal%npcol)
           if (target_myprow .eq. scal%myprow) then
              vcoul_temp(ig1_loc) = vcoul(ig1)
           endif
           if ( (target_myprow .ne. scal%myprow) .or. (target_mypcol .ne. scal%mypcol) ) then
              cycle
           endif
           !! Initialize the diagonal elements
           eps_temp(ig1_loc, ig2_loc) = (1.0D0, 0.0D0)
        enddo
        !$OMP PARALLEL DO collapse(3)
        do ifreq = 1, pol%nfreq
           do ig2_loc = 1, scal%npc
              do ig1_loc = 1, scal%npr
                 epsRDyn(ig1_loc, ig2_loc, ifreq) = eps_temp(ig1_loc, ig2_loc) - vcoul_temp(ig1_loc) * pol%chiRDyn(ig1_loc, ig2_loc, ifreq, 1)
              enddo
           enddo
        enddo
        !$OMP END PARALLEL DO
        SAFE_DEALLOCATE(vcoul_temp)
     endif

     if (pol%freq_dep .eq. 0) then
        SAFE_DEALLOCATE(pol%chi)
     else
        SAFE_DEALLOCATE(pol%chiRDyn)
     endif

     if (peinf%inode .eq. 0) then
        !! GPP
        if (pol%freq_dep .eq. 0) then
           write(6,'(1X,A,I6,A,2F15.8,A,I5)') 'q-pt ', iq, ':    Eps(G=0,Gp=0) = ', eps(pol%isrtxi(1), pol%isrtxi(1)), " isrtxi(1) = ", pol%isrtxi(1)
           !! FF
        else
           do ifreq = 1, pol%nfreq
              write(6,'(1X,A,I6,A,I6,A,2F15.8,A,I5)') 'q-pt ', iq, " freq ", ifreq, ':    Eps(G=0,Gp=0) = ', epsRDyn(pol%isrtxi(1),pol%isrtxi(1),ifreq), " isrtxi(1) = ", pol%isrtxi(1)
           enddo
        endif
     endif
     call timacc(5,2)

     call timacc(6,1)
     if (peinf%inode .eq. 0) then
        write(6,'(1X,A)')
        write(6,'(1X,A)') 'Inverting epsilon matrix'
     endif
     !! Invert eps to get epsinv using SCALAPACK
     if (pol%freq_dep .eq. 0) then
        call X(invert_matrix)(scal, desc_2d, pol%nmtx, eps)
     else
        do ifreq = 1, pol%nfreq
           eps_temp(:, :) = epsRDyn(:, :, ifreq)
           call zinvert_matrix(scal, desc_2d, pol%nmtx, eps_temp)
           epsRDyn(:, :, ifreq) = eps_temp(:, :)
        enddo
        SAFE_DEALLOCATE(eps_temp)
     endif

     if (peinf%inode .eq. 0) then
        !! GPP: only static frequency, omega = 0
        if (pol%freq_dep .eq. 0) then
           write(6,'(1X,A,I6,A,2F15.8,A,I5)') 'q-pt ', iq, ': Epsinv(G=0,Gp=0) = ', eps(pol%isrtxi(1), pol%isrtxi(1)), " isrtxi(1) = ", pol%isrtxi(1)
           !! FF
        else
           do ifreq = 1, pol%nfreq
              write(6,'(1X,A,I6,A,I6,A,2F15.8,A,I5)') 'q-pt ', iq, " ifreq ", ifreq, ': Epsinv(G=0,Gp=0) = ', epsRDyn(pol%isrtxi(1), pol%isrtxi(1), ifreq), " isrtxi(1) = ", pol%isrtxi(1)
           enddo
        endif
     endif

     SAFE_ALLOCATE(dmat_1d_block, (SCALARSIZE, pol%nmtx, MAX(npc, 1), pol%nfreq, 1))
     dmat_1d_block = 0.0D0
     SAFE_ALLOCATE(dmat_1d_block_diag, (SCALARSIZE, MAX(npc, 1), pol%nfreq))
     dmat_1d_block_diag = 0.0D0

     !! Copy epsRDyn to matrix_1d_block using pX(gemr2d)
     do ifreq = 1, pol%nfreq
        if (pol%freq_dep .eq. 0) then
           call pX(gemr2d)(pol%nmtx, pol%nmtx, eps, 1, 1, desc_2d, zmat_1d_block, 1, 1, desc_1d, cntxt_1d)
        else
           call pX(gemr2d)(pol%nmtx, pol%nmtx, epsRDyn(:,:,ifreq), 1, 1, desc_2d, zmat_1d_block, 1, 1, desc_1d, cntxt_1d)
        endif
        do ig2_loc = 1, npc
           do ig1 = 1, pol%nmtx
              !! Here the last index (spin index) is not used
              dmat_1d_block(1, ig1, ig2_loc, ifreq, 1) =  DBLE(zmat_1d_block(ig1, ig2_loc))
              dmat_1d_block(2, ig1, ig2_loc, ifreq, 1) = DIMAG(zmat_1d_block(ig1, ig2_loc))
           enddo
           !! Assign diagonal terms to dmat_1d_block_diag
           ig1 = INDXL2G(ig2_loc, block_size_col, peinf%inode, 0, peinf%npes)
           dmat_1d_block_diag(1, ig2_loc, ifreq) =  DBLE(zmat_1d_block(ig1, ig2_loc))
           dmat_1d_block_diag(2, ig2_loc, ifreq) = DIMAG(zmat_1d_block(ig1, ig2_loc))
        enddo
     enddo

     if (pol%freq_dep .eq. 0) then
        SAFE_DEALLOCATE(eps)
     else
        SAFE_DEALLOCATE(epsRDyn)
     endif
     call timacc(6,2)

     call timacc(7,1)
     if (peinf%inode .eq. 0) then
        write(6,'(1X,A)')
        write(6,'(1X,A)') 'Outputing '//TRUNC(filename_eps_hdf5)
        write(6,'(1X,A)')
     endif

     !! Open filename_eps_hdf5 for collective write
     call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
     call h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, info, error)
     call h5fopen_f(TRUNC(filename_eps_hdf5), H5F_ACC_RDWR_F, file_id, error, access_prp = plist_id)
     call h5pclose_f(plist_id, error)
     CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
     CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)

     !! Write mats/matrix
     if (npc > 0) then
        count_mat(:) = (/ SCALARSIZE, pol%nmtx, npc, pol%nfreq, 1, 1 /)
        offset_mat(:) = (/ 0, 0, ig2_offset - 1, 0, 0, iq - 1 /)
     else
        count_mat = 0
        offset_mat = 0
     endif

     call h5screate_simple_f(rank_mat, count_mat, memspace, error)
     if (error .ne. 0) then
        call die("HDF5 error", only_root_writes=.true.)
     endif

     call h5dopen_f(file_id, 'mats/matrix', dset_id, error)
     if (error .ne. 0) then
        call die("HDF5 error", only_root_writes=.true.)
     endif

     call h5dget_space_f(dset_id, filespace, error)
     if (error .ne. 0) then
        call die("HDF5 error", only_root_writes=.true.)
     endif

     call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset_mat, count_mat, error)
     if (error .ne. 0) then
        call die("HDF5 error", only_root_writes=.true.)
     endif

     if ( npc <= 0) then
        call H5sselect_none_f(filespace, error)
        if (error .ne. 0) then
           call die("HDF5 error", only_root_writes=.true.)
        endif
        call H5sselect_none_f(memspace, error)
        if (error .ne. 0) then
           call die("HDF5 error", only_root_writes=.true.)
        endif
     endif

     call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, dmat_1d_block, count_mat, error, mem_space_id = memspace, file_space_id = filespace, xfer_prp=plist_id)
     if (error .ne. 0) then
        call die("HDF5 error", only_root_writes=.true.)
     endif

     call h5sclose_f(memspace, error)
     call h5sclose_f(filespace, error)
     call h5dclose_f(dset_id, error)

     !! Write mats/matrix_diagonal
     if (npc > 0) then
        count_matdiag(:) = (/ SCALARSIZE, npc, pol%nfreq, 1 /)
        offset_matdiag(:) = (/ 0, ig2_offset - 1, 0, iq - 1 /)
     else
        count_matdiag = 0
        offset_matdiag = 0
     endif

     call h5screate_simple_f(rank_matdiag, count_matdiag, memspace, error)
     if (error .ne. 0) then
        call die("HDF5 error", only_root_writes=.true.)
     endif

     !! Select hyperslab in the file
     call h5dopen_f(file_id, 'mats/matrix-diagonal', dset_id, error)
     if (error .ne. 0) then
        call die("HDF5 error", only_root_writes=.true.)
     endif

     call h5dget_space_f(dset_id, filespace, error)
     if (error .ne. 0) then
        call die("HDF5 error", only_root_writes=.true.)
     endif

     call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset_matdiag, count_matdiag, error)
     if (error .ne. 0) then
        call die("HDF5 error", only_root_writes=.true.)
     endif

     if ( npc <= 0) then
        call H5sselect_none_f(filespace, error)
        if (error .ne. 0) then
           call die("HDF5 error", only_root_writes=.true.)
        endif
        call H5sselect_none_f(memspace, error)
        if (error .ne. 0) then
           call die("HDF5 error", only_root_writes=.true.)
        endif
     endif

     call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, dmat_1d_block_diag, count_matdiag, error, mem_space_id = memspace, file_space_id = filespace, xfer_prp=plist_id)
     if (error .ne. 0) then
        call die("HDF5 error", only_root_writes=.true.)
     endif

     call h5sclose_f(memspace, error)
     call h5sclose_f(filespace, error)
     call h5dclose_f(dset_id, error)

     call h5pclose_f(plist_id, error)
     call h5fclose_f(file_id, error)

     if (peinf%inode .eq. 0) then
        call set_qpt_done(TRUNC(filename_eps_hdf5), iq)
     endif

     SAFE_DEALLOCATE(dmat_1d_block)
     SAFE_DEALLOCATE(dmat_1d_block_diag)

     call timacc(7,2)
  enddo iq_loop !! iq

  SAFE_DEALLOCATE(vcoul)
  SAFE_DEALLOCATE(zmat_1d_block)
  SAFE_DEALLOCATE_P(pol%isrtx)
  SAFE_DEALLOCATE_P(pol%isrtxi)
  SAFE_DEALLOCATE_P(pol%nmtx_of_q)
  SAFE_DEALLOCATE_P(pol%qpt)
  SAFE_DEALLOCATE_P(pol%dFreqGrid)
  SAFE_DEALLOCATE_P(pol%dFreqBrd)
  SAFE_DEALLOCATE_P(gvec%components)
  SAFE_DEALLOCATE(ekin)

#ifdef MPI
  call MPI_BARRIER(MPI_COMM_WORLD, mpierr)
#endif

  !! Time accounting
  nroutnam = 7
  SAFE_ALLOCATE(routnam, (nroutnam))

  routnam(1)='TOTAL:'
  routnam(2)='INIT:'
  routnam(3)='READ CHI:'
  routnam(4)='CALC VCOUL:'
  routnam(5)='GET EPS:'
  routnam(6)='INVERT EPS:'
  routnam(7)='OUTPUT:'

  SAFE_ALLOCATE(routsrt, (nroutnam))
  routsrt=(/ (ii, ii=2,nroutnam), 1 /)

  call timacc(1,2)

  if (peinf%inode .eq. 0) then
     write(6,*)
     write(6,9000) 'CPU (s)','WALL (s)','#'
     write(6,*)
  endif

  do ii = 1, nroutnam
     call timacc(routsrt(ii), 3, tsec, ncount)
#ifdef MPI
     call MPI_ALLREDUCE(tsec, tmin, 2, MPI_REAL_DP, MPI_MIN, MPI_COMM_WORLD, mpierr)
     call MPI_ALLREDUCE(tsec, tmax, 2, MPI_REAL_DP, MPI_MAX, MPI_COMM_WORLD, mpierr)
#else
     tmin = tsec
     tmax = tsec
#endif
     if (peinf%inode .eq. 0) then
        write(6,9001) routnam(routsrt(ii)),tmin(1),tmin(2),ncount
        write(6,9002) tsec(1),tsec(2)
        write(6,9003) tmax(1),tmax(2)
     endif
  enddo

9000 format(23x,a13,3x,a13,3x,a8)
9001 format(1X,A16,'(min.)',f13.3,3x,f13.3,3x,i8)
9002 format(   17x,'(PE 0)',f13.3,3x,f13.3)
9003 format(   17x,'(max.)',f13.3,3x,f13.3)

  call h5close_f(error)

#ifdef MPI
  call MPI_FINALIZE(mpierr)
#endif
contains

  subroutine dinvert_matrix(scal, desc_2d, nmtx, matrix)
    type (scalapack), intent(in) :: scal
    integer, intent(in) :: desc_2d(9)
    integer, intent(in) :: nmtx
    real(DP), intent(inout) :: matrix(scal%npr,scal%npc)
    integer :: info, lwork, liwork, ipiv(scal%npr+scal%nbr)
    integer, allocatable :: iwork(:)
    real(DP), allocatable :: work(:)

    !! LU factorization of a general m x n distributed matrix
    call pdgetrf(nmtx, nmtx, matrix, 1, 1, desc_2d, ipiv, info)
    if (info /= 0) then
       if (peinf%inode == 0) write(0,*) 'ERROR: got info = ', info, ' in p?getrf'
       call die('p?getrf failed')
    endif

    !! Computes the inverse of a LU-factored distributed matrix
    SAFE_ALLOCATE(work, (10))
    SAFE_ALLOCATE(iwork, (10))
    !! http://www.netlib.org/scalapack/explore-html/d2/d63/pzgetri_8f_source.html
    !! lwork = -1 : workspace query
    !! liwork = -1 : workspace query
    !! on exit, matrix contains the inverse of the original matrix
    call pdgetri(nmtx, matrix, 1, 1, desc_2d, ipiv, work, -1, iwork, -1, info)
    if (info /= 0) then
       if (peinf%inode == 0) write(0,*) 'ERROR: got info = ', info, ' in p?getri'
       call die('p?getri failed for query mode')
    endif
    lwork = MAX(NINT(work(1)),1)
    liwork = MAX(iwork(1),1)
    SAFE_DEALLOCATE(work)
    SAFE_DEALLOCATE(iwork)

    SAFE_ALLOCATE(work, (lwork))
    SAFE_ALLOCATE(iwork, (liwork))
    call pdgetri(nmtx, matrix, 1, 1, desc_2d, ipiv, work, lwork, iwork, liwork, info)
    if (info /= 0) then
       if (peinf%inode == 0) write(0,*) 'ERROR: got info = ', info, ' in p?getri'
       call die('p?getri failed')
    endif
    SAFE_DEALLOCATE(iwork)
    SAFE_DEALLOCATE(work)
  end subroutine dinvert_matrix

  subroutine zinvert_matrix(scal, desc_2d, nmtx, matrix)
    type (scalapack), intent(in) :: scal
    integer, intent(in) :: desc_2d(9)
    integer, intent(in) :: nmtx
    complex(DPC), intent(inout) :: matrix(scal%npr,scal%npc)
    integer :: info, lwork, liwork, ipiv(scal%npr+scal%nbr)
    integer, allocatable :: iwork(:)
    complex(DPC), allocatable :: work(:)

    !! LU factorization of a general m x n distributed matrix
    call pzgetrf(nmtx, nmtx, matrix, 1, 1, desc_2d, ipiv, info)
    if (info /= 0) then
       if (peinf%inode == 0) write(0,*) 'ERROR: got info = ', info, ' in p?getrf'
       call die('p?getrf failed')
    endif

    !! Computes the inverse of a LU-factored distributed matrix
    SAFE_ALLOCATE(work, (10))
    SAFE_ALLOCATE(iwork, (10))
    !! http://www.netlib.org/scalapack/explore-html/d2/d63/pzgetri_8f_source.html
    !! lwork = -1 : workspace query
    !! liwork = -1 : workspace query
    !! on exit, matrix contains the inverse of the original matrix
    call pzgetri(nmtx, matrix, 1, 1, desc_2d, ipiv, work, -1, iwork, -1, info)
    if (info /= 0) then
       if (peinf%inode == 0) write(0,*) 'ERROR: got info = ', info, ' in p?getri'
       call die('p?getri failed for query mode')
    endif
    lwork = MAX(NINT(DBLE(work(1))),1)
    liwork = MAX(iwork(1),1)
    SAFE_DEALLOCATE(work)
    SAFE_DEALLOCATE(iwork)

    SAFE_ALLOCATE(work, (lwork))
    SAFE_ALLOCATE(iwork, (liwork))
    call pzgetri(nmtx, matrix, 1, 1, desc_2d, ipiv, work, lwork, iwork, liwork, info)
    if (info /= 0) then
       if (peinf%inode == 0) write(0,*) 'ERROR: got info = ', info, ' in p?getri'
       call die('p?getri failed')
    endif
    SAFE_DEALLOCATE(iwork)
    SAFE_DEALLOCATE(work)
  end subroutine zinvert_matrix

end program EpsInv

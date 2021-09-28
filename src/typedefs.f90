#include "f_defs.h"

!==========================================================================
!
! Modules:
!
! (1) typedefs      Originally By GMR      Last Modified Meng Wu
!
!    Derived types that are used throughout the code.
!
!==========================================================================

module typedefs_m

  use nrtype_m
  implicit none
  public

  !---------------------------

  type crystal
     real(DP) :: celvol !< cell volume in real space (a.u.)
     real(DP) :: recvol !< cell volume in reciprocal space (a.u.)
     real(DP) :: alat !< lattice constant in real space (a.u.)
     real(DP) :: blat !< lattice constant in reciprocal space (a.u.)
     real(DP) :: avec(3,3) !< lattice vectors in real space (alat) = [a1,a2,a3]
     real(DP) :: bvec(3,3) !< lattice vectors in reciprocal space (blat) = [b1,b2,b3]
     real(DP) :: adot(3,3) !< metric tensor in real space (a.u.)
     real(DP) :: bdot(3,3) !< metric tensor in reciprocal space (a.u.)
     integer :: nat !< number of atoms
     integer, pointer :: atyp(:) !< atomic species, atyp(1:nat)
     real(DP), pointer :: apos(:,:) !< atomic positions, apos(1:3,1:nat) (alat)
  end type crystal

  !---------------------------

  type kpoints
     integer :: nspinor = 1 !< nspinor = 2 if doing two-component spinor calculation; 1 is default
     integer :: nspin   !< nspin = 1 or 2; nspin = 1 when npsinor = 2
     integer :: nrk     !< number of k-points
     integer :: mnband  !< max number of bands
     integer :: nvband  !< number of valence bands
     integer :: ncband  !< number of conduction bands
     integer  :: kgrid(3) !< Monkhorst-Pack number of k-points in each direction
     real(DP) :: shift(3) !< Monkhorst-Pack shift of grid
     real(DP) :: ecutwfc            !< wave-function cutoff, in Ry
     integer, pointer :: ngk(:)     !< number of g-vectors for each k-point
     integer :: ngkmax              !< max(ngk(:))
     integer, pointer :: ifmin(:,:) !< lowest occupied band (kpoint,spin)
     integer, pointer :: ifmax(:,:) !< highest occupied band (kpoint,spin)
     real(DP), pointer :: w(:)      !< weights (kpoint) (between 0 and 1)
     real(DP), pointer :: rk(:,:)   !< k-vector (3, kpoint) in crystal coords
     real(DP), pointer :: el(:,:,:) !< band energies (band, kpoint, spin)
     real(DP), pointer :: elda(:,:,:) !< band energies before eqp correction
     real(DP), pointer :: occ(:,:,:)  !< occupations (between 0 and 1)
     integer, pointer :: degeneracy(:,:,:) !< size of deg. subspace for (band, kpoint, spin)
  end type kpoints

  !---------------------------

  type symmetry
     integer :: ntran         !< number of operations in full group
     integer :: mtrx(3,3,48)  !< symmetry matrix
     integer :: mtrx_reci(3,3,48) ! < symmetry matrix in basis of reciprocal space lattice vectors, b1, b2, b3
     real(DP) :: mtrx_cart(3,3,48) ! < symmetry matrix in cartesian coordinates
     real(DP) :: tnp(3,48)    !< fractional translations
     integer :: kgzero(3,48)  !< Umklapp vectors for subgroup symmetry operations
     integer :: cell_symmetry !< 0 = cubic, 1 = hexagonal
  end type symmetry

  !---------------------------

  type grid
     integer :: nr  !< number in reduced zone
     integer :: nf  !< number in full zone
     real(DP) :: sz !< radius of a spherical subzone equivalent to
     integer, pointer :: itran(:) !< sym op to go from irrbz to fullbz
     integer, pointer :: indr(:)  !< irrbz k/q-point mapped to fullbz
     integer, pointer :: kg0(:,:) !< Umklapp vectors (for Wigner-Seitz cell)
     real(DP), pointer :: r(:,:)  !< k/q-points in reduced zone
     real(DP), pointer :: f(:,:)  !< k/q-points in full zone
     integer, pointer :: dege_r(:)
     integer, pointer :: findex(:) !< from little group irk_little_group to FBZ kpoint ifk
  end type grid

  type grid_relation
     integer, pointer :: rindex(:)
     integer, pointer :: symsindex(:)
     integer, pointer :: gumk(:,:)
  end type grid_relation

  !-----------------------------------

  type gspace
     integer :: ng       !< number of G-vectors
     integer :: nFFTgridpts !< number in FFT grid = product(FFTgrid(1:3))
     real(DP) :: ecutrho !< charge-density cutoff, in Ry
     integer, pointer :: components(:,:) !< the G-vectors, in units of 2 pi / a
     integer :: FFTgrid(3)  !< gsm: FFTgrid is the size of the FFT grid, not the maximum G-vector
     integer, pointer :: index_vec(:) ! mapping to FFT grid
     real(DP), pointer :: ekin(:) !< kinetic energy of G-vectors
  end type gspace

  !---------------------------

  type polarizability
     logical :: serial_output
     integer :: serial_output_nc
     logical :: low_mem
     integer :: eqp_start, eqp_end
     integer :: freq_dep        !> frequency dependence of the inverse dielectric matrix
     !> 0: static calculation 2: full frequency 3: two imaginary frequencies
     integer :: freq_dep_method !< full frequency calculation. 0: Adler-Wiser; 1: Shishkin and Kresse 2006
     integer :: nFreq           !< number of frequencies used in full frequency calculation
     integer :: nfreq_imag, nfreq_real      !< number of imaginary freqs for CD (also 1 for GN GPP)
     real(DP) :: dInitFreq      !< initial frequency (eV) for polarizability energy denominator
     real(DP) :: dDeltaFreq     !< frequency increment (eV) for polarizability energy denominator
     real(DP) :: dBrdning       !< Lorentzian broadening (eV) for polarizability energy denominator
     integer :: nBrdning
     real(DP) :: Brdning_stepsize
     real(DP), pointer :: dFreqGrid(:) !< Grid of Frequencies for Full Frequency
     real(DP) :: dFreqCutoff
     real(DP) :: delta_freq_imag !> first non-zero imaginary frequency

     integer :: nSFreq    !< number of frequencies used in spectral function
     real(DP) :: dInitSFreq  !< initial frequency (eV) for polarizability spectral function
     real(DP) :: dDeltaSFreq !< frequency increment (eV) for polarizability spectral function
     real(DP), pointer :: dSFreqGrid(:) !< Grid of Frequencies for spectral function
     real(DP) :: dSFreqStepIncrease
     real(DP) :: dSFreqCutoff1
     real(DP) :: dSFreqCutoff2

     logical :: has_advanced !< Do we store eps_A or just eps_R?
     integer :: matrix_type !< 0 to write epsilon^{-1}, 1 for epsilon, 2 for chi0.
     integer :: nmatrix !< has_advanced+1. Multiply by nspin if matrix_type==2
     integer :: matrix_flavor !< 2 (=CMPLX), unless we have freq_dep==0 and SCALARSIZE==1.

     logical :: eqp_corrections !< are we using eqp.dat and eqp_q.dat files
     complex(DPC), pointer :: dFreqBrd(:)  !< Corresponding Broadenings for Full Frequency
     integer :: fullConvLog !< logging pol matrix head & tail convergence
     integer :: nmtx
     integer, pointer :: nmtx_of_q(:)
     integer :: qgrid(3)
     integer :: nq0, nq1, nq !< Number of q->0 points, q/=0 points, and total number of q-points
     logical :: subsample !< whether we have more than one q0 point (used in subsampled calculation)
     logical :: non_uniform !< do non-uniform sampling using Voronoi decomposition of BZ
     integer :: gcomm
     logical :: min_fftgrid   !< use the smallest possible fftbox
     integer :: os_opt_ffts       !< optimizes calculation/reuse of FFTs (real-space WFNs)
     logical :: os_hdf5           !< use parallel IO?
     logical :: restart        !< are we restarting the calculation? Only ok with HDF5
     integer :: stop_after_qpt !< pretend the calculation was prematurely killed after this qpt (-1=don`t kill)
     integer :: intraband_flag !< 0=regular calculation, 1=only include intraband, 2=only interband
     real(DP) :: intraband_overlap_min !< a transition is intraband if |<uvk|uvk+q>| is larger than this
     logical :: patched_sampling !< Do we have only a patch in the BZ?
     integer :: WFN_FFTgrid(3)!< max. size FFTgrid that holds all WFNs
     integer :: FFTgrid(3)    !< FFT grid to use (RHO or economical one)
     logical :: skip_epsilon
     logical :: skip_chi
     logical :: use_hdf5      !< with -DHDF5, whether or not we actually use hdf5
     logical :: need_WFNq     !< will we need the WFNq file? (nq0>0.and.valueq0==1.and.iqexactlyzero==0)
     integer, pointer :: irow(:)
     integer, pointer :: isrtx(:)
     integer, pointer :: isrtxi(:)
     integer :: icutv               !< icutv encodes presence and type of truncation
     real(DP) :: truncval(3)   !< in Bohr (au)
     real(DP), pointer :: qpt(:,:)
     SCALAR, allocatable :: gme(:,:,:,:,:)
     SCALAR, allocatable :: gme2(:,:,:,:)
     SCALAR, allocatable :: chi(:,:,:)

     integer :: ncrit
     real(DP) :: efermi
     real(DP) :: efermi_input
     logical :: rfermi
     real(DP) :: ecuts    !< energy cutoff of screened coulomb interaction in Ry
     !> Reference regarding retarded/advanced functions: Catalin`s thesis, Eq. (1.44)
     complex(DPC), allocatable :: chiRDyn(:,:,:,:) !< Retarded polarizability
     complex(DPC), allocatable :: chiTDyn(:,:,:,:) !< Spectral function of polarizability
     real(DP), pointer :: edenDyn(:,:,:,:) !< Dynamic energy denominator
     real(DP), pointer :: edenDyn2(:,:,:) !< Dynamic energy denominator
     logical :: degeneracy_check_override
     real(DP) :: lin_denominator !< energy threshold below which to activate lin_denominator
     real(DP) :: de_min, de_max
     ! velocities for calculating linearized denominator in dynamic case
     real(DP) :: imaginary_frequency  !< purely imaginary frequency used in Godby-Needs GPP
     ! variables for subspace truncation method in epsilon
     logical  :: subspace
     real(DP) :: chi_eigenvalue_cutoff
     integer  :: neig_sub_input
     logical  :: use_elpa
     logical  :: need_full_chi
     logical  :: keep_full_eps_static
     logical  :: matrix_in_subspace_basis
     integer  :: nrow_local_sub, ncol_local_sub, neig_sub
     complex(DPC), allocatable :: chiRDyn_sym_omega0(:,:)
     complex(DPC), allocatable :: eigenvect_omega0(:,:)
     real(DP), allocatable :: eigenval_omega0(:)
     real(DP), allocatable :: vcoul_sub(:)
     ! variables for nonblocking scheme
     logical :: nonblocking_cyclic
     logical :: mtxel_2
     integer :: skip_nvb, nvb
     integer :: skip_ncb, ncb
     integer :: nband
     logical :: zero_vq0
     logical :: timeordered, resetrealfreq, k_plus_q, time_reversal, full_axis_frequency
     integer, allocatable :: indexq(:)
  end type polarizability

  !> FHJ: mean-field header
  type mf_header_t
     integer :: version
     character(len=3) :: sheader
     character(len=32) :: sdate
     character(len=32) :: stime
     integer :: iflavor
     type(crystal) :: crys
     type(kpoints) :: kp
     type(symmetry) :: syms
     type(gspace):: gvec
  end type mf_header_t

  type wavefunction
     integer :: ng
     integer :: nband
     integer :: nspin
     integer :: nspinor = 1 !< nspinor = 2 if doing two-component spinor calculation; 1 is default
     integer, pointer :: isort(:)
     ! cg(ig, ib, is)     
     SCALAR, pointer :: cg(:,:,:)
  end type wavefunction
  
  type int_wavefunction
     integer :: nspin
     integer :: nspinor = 1 !< nspinor = 2 if doing two-component spinor calculation; 1 is default
     integer, pointer ::  ng(:)     !< (nk)
     integer, pointer :: isort(:,:) !< (ngmax, nk)
     integer, pointer :: cbi(:)
     !> I think this can be decommissioned if we use kp%rk instead
     real(DP), pointer :: qk(:,:)
     ! cg(ig, ib, is)
     SCALAR, pointer :: cg(:,:,:)
     ! cgk(ig, ib, is, ik)     
     SCALAR, pointer :: cgk(:,:,:,:)
  end type int_wavefunction
  
end module typedefs_m

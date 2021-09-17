#include "f_defs.h"

!==========================================================================
!
! Modules:
!
! (1) typedefs      Originally By GMR      Last Modified 7/8/2008 (JRD)
!
!>    Derived types that are used throughout the code.
!
!==========================================================================

module typedefs_m

  use nrtype_m

  implicit none

  public ! only types in this module

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
     ! integer :: ntranq        !< number of operations in small group of q
     ! real(DP) :: rq(3)        !< The q-point this ntranq belongs to
     
     integer :: mtrx(3,3,48)  !< symmetry matrix
     ! <<<<<<
     ! mtrx_reci(1:3,1:3,itran) = (A \cdot A^{T}) \cdot mtrx(1:3,1:3,itran) \cdot (A \cdot A^{T})^{-1}
     ! ------------------------
     ! at = A^{T}
     ! bg = B
     ! ------------------------
     integer :: mtrx_reci(3,3,48) ! < symmetry matrix in basis of reciprocal space lattice vectors, b1, b2, b3
     real(DP) :: mtrx_cart(3,3,48) ! < symmetry matrix in cartesian coordinates
     ! >>>>>>
     real(DP) :: tnp(3,48)    !< fractional translations
     ! integer :: indsub(48)    !< symmetry operations in subgroup of q
     integer :: kgzero(3,48)  !< Umklapp vectors for subgroup symmetry operations
     integer :: cell_symmetry !< 0 = cubic, 1 = hexagonal
  end type symmetry

  !---------------------------

  type grid
     integer :: nr  !< number in reduced zone
     integer :: nf  !< number in full zone
     real(DP) :: sz !< radius of a spherical subzone equivalent to
     ! !! one point in the set f
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
  !> Parameters for scissors operators
  !> e_cor = e_in + es + edel * (e_in - e0)
  !> es and e0 are in eV. edel is a dimensionless slope.
  type sub_scissors_t
     real(DP) :: es
     real(DP) :: edel
     real(DP) :: e0
  end type sub_scissors_t

  type scissors_t
     type(sub_scissors_t) :: val
     type(sub_scissors_t) :: cond
  end type scissors_t

  !---------------------------

  type wavefunction
     integer :: ng
     integer :: nband
     integer :: nspin
     integer :: nspinor = 1 !< nspinor = 2 if doing two-component spinor calculation; 1 is default
     integer, pointer :: isort(:)
     SCALAR, pointer :: cg(:,:,:)
  end type wavefunction

  type wavefunction_kerdiag
     integer :: ng
     integer :: nband
     integer :: nspin
     integer :: nspinor = 1 !< nspinor = 2 if doing two-component spinor calculation; 1 is default
     integer, allocatable :: isort(:)
     SCALAR, allocatable :: cg(:,:,:)
  end type wavefunction_kerdiag
  
  !---------------------------
  !> For Epsilon: this is the wavefunction before unfolding the irr. BZ.
  type int_wavefunction
     integer :: nspin
     integer :: nspinor = 1 !< nspinor = 2 if doing two-component spinor calculation; 1 is default
     ! integer, pointer ::  ng(:)     !< (nk)
     ! integer, pointer :: isort(:,:) !< (ngmax, nk)
     ! ! integer, pointer :: cbi(:)
     ! ! [NEW] cg(ig,ib,is,ik)
     ! SCALAR, pointer :: cg(:,:,:,:)
     ! SCALAR, pointer :: cgk(:,:,:,:)

     integer, allocatable ::  ng(:)     !< (nk)
     integer, allocatable :: isort(:,:) !< (ngmax, nk)
     ! integer, pointer :: cbi(:)
     ! [NEW] cg(ig,ib,is,ik)
     SCALAR, allocatable :: cg(:,:,:,:)
     SCALAR, allocatable :: cgk(:,:,:,:)
     
  end type int_wavefunction

  ! For BSE
  type int_wavefunction_bse
     integer :: nspin
     integer :: nspinor = 1 !< nspinor = 2 if doing two-component spinor calculation; 1 is default
     integer, pointer ::  ng(:)     !< (nk)
     integer, pointer :: isort(:,:) !< (ngmax, nk)
     integer, pointer :: cbi(:)
     !> I think this can be decommissioned if we use kp%rk instead
     real(DP), pointer :: qk(:,:)
     ! [OLD] cg(ig,ik/ib,is)
     SCALAR, pointer :: cg(:,:,:)
     ! [NEW] cg(ig,ib,is,ik)
     ! SCALAR, pointer :: cg(:,:,:,:)
     SCALAR, pointer :: cgk(:,:,:,:)
  end type int_wavefunction_bse

  !------------------------------------

  !> FHJ: valence WFNs for 1 particular kpt and 1 band
  !! It stores all bands in the case of real-space WFNs
  type valence_wfns
     integer :: nband  !< This is indeed the number of valence bands
     ! integer :: ncore_excl
     !>This is the number of core states that are not included in the valence in polarizability calculation.
     integer :: ngv    !< Number of G-vectors
     integer :: idx_kp !< Idx of current kpt in kp/kpq structure
     ! integer, pointer :: isort(:)
     ! !> (nband+ncrit,spin). Note: the nband+ncrit index is actually useless!
     ! real(DP), pointer :: ev(:,:)
     ! SCALAR, pointer :: zv(:,:)   !< (ngv,spin)
     integer, allocatable :: isort(:)
     !> (nband+ncrit,spin). Note: the nband+ncrit index is actually useless!
     real(DP), allocatable :: ev(:,:)
     SCALAR, allocatable :: zv(:,:)   !< (ngv,spin)
     
     ! !> real-space wavefunction for all "local" val. bands (fft1,fft2,fft3,band)
     ! complex(DPC), pointer :: wfn_fft(:,:,:,:)
  end type valence_wfns

  !-------------------------------

  !> FHJ: conduction WFNs for 1 particular kpt and all bands (!) the processor owns
  type conduction_wfns
     integer :: nband  !< This is actually the number of valence+conduction bands!
     integer :: ngc    !< Number of G-vectors
     integer :: idx_kp !< Idx of current kpt in kp structure
     ! integer, pointer :: isort(:)
     ! real(DP), pointer :: ec(:,:) !< (nband,nspin)
     ! SCALAR, pointer :: zc(:,:)   !< (ngc*ncownactual,spin)
     integer, allocatable :: isort(:)
     real(DP), allocatable :: ec(:,:) !< (nband,nspin)
     SCALAR, allocatable :: zc(:,:)   !< (ngc*ncownactual,spin)
     
     ! !> real-space wavefunction for all "local" cond. bands (fft1,fft2,fft3,band)
     ! complex(DPC), pointer :: wfn_fft(:,:,:,:)
  end type conduction_wfns

  !----------------------------

  !> splines knots and coefficients
  type spline_tck
     integer :: n              !< number of knots
     integer :: k              !< degree of spline (1=linear, etc.)
     real(DP), pointer :: t(:) !< position of knots
     real(DP), pointer :: c(:) !< splines coefficient
  end type spline_tck

  !-------------------------------------

  type coulomb_modifier_t

     real(DP) :: short_range_frac_fock  !< Short range exchange fraction
     real(DP) :: long_range_frac_fock   !< Long range exchange fraction
     real(DP) :: screening_length       !< Screening length
     !< The above 3 parameters are used
     !< only for TDDFT and sigma calculations.

  end type coulomb_modifier_t

  !----------------------------

  type siginfo
     logical :: calc_vcoul0 !> If true, calculate vcoul0 instead of wcoul0
     real(DP) :: qlen_threshold !> control the range of q to do weighted linear regression
     real(DP) :: epsinvhead !> used in truncated calculations
     real(DP) :: finiteq(3) !> used in calc_wcoul0_pchip.f90
     !> Deprecated variables, only for Sigma
     logical :: low_io
     logical :: need_advanced
     logical :: gamma_only       !< nonzero if Gamma is the only q-point, 0 otherwise
     real(DP) :: avgcut         !< Cut in which we do cell averages on
     SCALAR :: wcoul0
     logical :: freplacebz
     logical :: fwritebz
     real(DP) :: q0vec(3)
     logical :: averagew
     type(coulomb_modifier_t) :: coulomb_mod
     logical :: coul_mod_flag   !< Flag which tells if the coulomb interaction has been
     ! modified (for hybrid functional like calculation in sigma)
     !> -----------------
     integer :: eqp_start, eqp_end
     logical :: average_vq     
     real(DP) :: varepsilon
     logical :: arem !> adaptive remainder
     logical :: real_gpp !> whether use complex \hbar \tilde{\omega} in GPP model, complex_gpp
     logical :: full_axis_frequency
     integer :: signmc, max_iteration !, signmc_mbz
     real(DP) :: mc_diff_threshold
     integer :: freq_dep    !< frequency dependence of the inverse dielectric matrix
     integer :: freq_dep_method    !< frequency dependence method of the inverse dielectric matrix
     integer :: nFreq       !< number of frequencies used in full frequency calculation
     real(DP) :: dDeltaFreq !< frequency increment (eV) for polarizability energy denominator
     real(DP) :: dBrdning   !< Lorentzian broadening (eV) for polarizability energy denominator
     real(DP), pointer :: dFreqGrid(:)   !< Grid of Frequencies for Full Frequency
     complex(DPC), pointer :: dFreqBrd(:) !< Corresponding Broadenings for Full Frequency
     integer :: nSFreq    !< number of frequencies used in spectral functions
     real(DP) :: dDeltaSFreq !< frequency increment (eV) for spectral functions
     real(DP), pointer :: dSFreqGrid(:)   !< Grid of Frequencies for spectral functions
     integer :: exact_ch    !< compute the exact static CH
     integer :: fullConvLog !< logging CH convergence
     ! real(DP) :: tol        !< energy tolerance for degeneracy
     logical :: use_hdf5     !< with -DHDF5, whether or not we actually use hdf5
     logical :: wfn_hdf5     !< with -DHDF5, use HDF5 for WFN files?
     logical :: use_xdat     !< use saved exchange matrix elements from file x.dat
     logical :: use_vxcdat   !< use saved exchange-correlation matrix elements from file vxc.dat
     logical :: use_vxc2dat  !< use saved exchange-correlation matrix elements from file vxc2.dat
     integer :: nkn          !< number of k-points on which to calculate Sigma (from sigma.inp)
     integer :: nfreq_imag
     integer :: cd_int_method !< Integration method for CD calculations
     logical :: cspline
     integer :: invalid_gpp_mode !< How to treat invalid GPP modes? See sigma.inp
     integer :: nq0, nq1, nq !< Number of q->0 points, q/=0 points, and total number of q-points
     logical :: subsample !< whether we perform a subsampling of the voronoi cell of q=0
     real(DP), pointer :: subweights(:) !< (nq0) weight for each subsampling q-point
     integer :: nvband       !< number of bands in bare exchange
     integer :: ntband       !< number of bands in dynamical sigma summation
     !> number of core states that are not included in both the bare exchange and dynamical sigma summations
     integer :: nspin
     integer :: spin_index(2)
     integer :: icutv               !< icutv encodes presence and type of truncation
     integer :: iuseremainder
     integer :: qgrid(3)
     integer :: iscreen !< what type of screening is present. 0 = semiconductor, 1 = graphene, 2 = metal
     integer :: fdf !< finite difference form for numerical derivative of Sigma
     integer, pointer :: indkn(:) !< mapping of k-points from sigma.inp to those in kp%rk from WFN files
     integer, pointer :: diag(:) !< energy bands for which Sigma diagonal matrix elements are calculated
     integer :: ndiag    !< number of bands contained in diag(:)
     integer :: ib_outer_min_abso_global, ib_outer_max_abso_global
     logical :: vxc_hdf5
     integer :: noffdiag !< offdiag
     integer :: loff
     integer :: toff
     integer :: bmin
     integer :: bmax
     integer, pointer :: off1(:) !< offdiag <bra|
     integer, pointer :: off2(:) !< offdiag |ket>
     integer, pointer :: off3(:) !< offdiag energy at which to evaluate
     integer, pointer :: offmap(:,:)  !< sig%off1(ii) = sig%diag(sig%offmap(ii,1))
     real(DP) :: dw       !< finite difference spacing for numerical derivative of Sigma in eV
     real(DP) :: ecutb    !< energy cutoff of bare coulomb interaction in Ry
     real(DP) :: ecuts    !< energy cutoff of screened coulomb interaction in Ry
     real(DP) :: xfrac    !< fraction of bare exchange
     real(DP) :: gamma    !< GPP broadening
     real(DP) :: sexcut   !< GPP SX cutoff
     integer :: freq_grid_shift !< How to shift the requency grid. See sigma.inp for more info.
     integer :: nfreqeval
     real(DP) :: freqevalmin
     real(DP) :: freqevalstep
     logical :: eqp_corrections !< are we using eqp.dat
     logical :: eqp_outer_corrections !< are we using eqp_outer.dat
     type(scissors_t) :: scis
     type(scissors_t) :: scis_outer
     logical  :: wfn_outer_present
     type(spline_tck) :: spl_tck       !< scissors b-spline coefficients
     type(spline_tck) :: spl_tck_outer !< scissors b-spline coeffs for WFN_outer
     real(DP) :: truncval(3)    !< in Bohr (au)
     real(DP), pointer :: kpt(:,:)
     real(DP), pointer :: qpt(:,:) !< (3,nq) q-points in eps0mat+epsmat files, or provided in sigma.inp
     integer :: ncrit           !< number of partially occupied bands
     real(DP) :: efermi         !< Fermi level
     real(DP) :: efermi_input   !< The value to set E_Fermi in the input file, in eV
     logical :: rfermi          !< Measure the new Fermi level relative to that of the neutral system
     SCALAR, pointer :: vxc(:,:) !< Vxc(G)
     SCALAR, pointer :: vxc2(:,:) !< Vxc(G) for hybrid functional type calculations
     logical :: degeneracy_check_override
     logical :: offdiagsym
     logical :: qgridsym
     logical :: die_outside_sphere
     logical :: sigma_correction !< if .true., compute only a correction to the QP self energy, ie,
     !! don`t subtract vxc and don`t compute the bare exchange.
     ! variables for subspace truncation
     logical :: subspace_q0, matrix_in_subspace_basis_q0, keep_full_eps_static_q0
     logical :: subspace, matrix_in_subspace_basis, keep_full_eps_static
     integer :: neig_sub_max
     real(DP) :: eqpmin_outer, eqpmax_inner
     integer :: cd_ny
     real(DP) :: ekinx_min !> |q|^2 of the smallest non-zero q1
     integer :: celltype !> FCC = 1, hexagional=2, other=0
  end type siginfo

  type vxcinfo
     integer :: nspin
     integer :: nspinor
     integer :: nrk
     integer :: ndiag, diag_nmin, diag_nmax
     integer :: noffdiag, offdiag_nmin, offdiag_nmax
     real(DP), pointer :: rk(:,:)   !< k-vector (3, kpoint) in crystal coords
     SCALAR, pointer :: mtxeld(:,:) ! (ib_diag_rela_global,ik)
     SCALAR, pointer :: mtxelo(:,:,:) ! (ib1_offdiag_rela_global,ib2_offdiag_rela_global,ik_inner)
     real(DP), pointer :: mtxeld_real(:,:,:) ! (1:2,ib_diag_rela_global,ik)
     real(DP), pointer :: mtxelo_real(:,:,:,:) ! (1:2,ib1_offdiag_rela_global,ib2_offdiag_rela_global,ik_inner)
     integer, pointer :: kmap_outer2vxc(:), kmap_vxc2outer(:)
  end type vxcinfo

  !---------------------------
  !> Dielectric matrix info using comm_mpi
  !! In BLACS terms, we distribute the columns of the epsinv matrix in a
  !! 1D-block-cyclic fashion. Currently, the block size nb is 1
  !! The following old quantities were removed, as they can be calculated on the fly:
  !! epsmpi%igp_owner(igp) = INDXG2P(igp, epsmpi%nb, peinf%pool_rank, 0, peinf%npes_pool)
  !! epsmpi%igp_index(igp) = INDXG2L(igp, epsmpi%nb, peinf%pool_rank, 0, peinf%npes_pool)
  !! The following array should be removed soon:
  !! epsmpi%inv_igp_index(igp_loc) = INDXL2G(igp_loc, epsmpi%nb, peinf%pool_rank, 0, peinf%npes_pool)
  type epsmpiinfo
     integer :: nb !< block size. Currently set to 1 (round robin)
     integer :: nrq !< number of q vectors in eps0mat.h5 + epsmat.h5
     integer :: ngpown !< number of columns I own. Same as numroc(neps, nb, pool_rank, 0, npes_pool)
     integer :: ngpown_max !< number of columns owned by pool_rank==0.
     integer, pointer :: isrtq(:,:)   !< These 3 arrays have a dimension of (1:gvec%ng,1:(nq+1)) where
     integer, pointer :: isrtqi(:,:)  !! (nq+1) is the total # of q`s including q0.

     integer, pointer :: inv_igp_index(:)
     integer, pointer :: nmtx(:)      !< dimension of eps(q) for each q
     real(DP), pointer :: rq(:,:)
     SCALAR, pointer :: eps(:,:,:)    !< eps(1:gvec%ng,1:ngpown,1:(nq+1))

     !> dimension of epsR and epsA (1:gvec%ng,1:ngpown,1:sig%nFreq,1:(nq+1))
     complex(DPC), pointer :: epsR(:,:,:,:)
     complex(DPC), pointer :: epsA(:,:,:,:)
  end type epsmpiinfo

  !---------------------------

  type wfnkqmpiinfo
     integer, pointer ::  nkptotal(:)
     integer, pointer :: isort(:,:)
     integer, pointer :: band_index(:,:)
     real(DP), pointer :: qk(:,:)
     real(DP), pointer :: el(:,:,:)   !> energies used in calculating Sigma (mtxel_cor), can use eqp.dat if eqp_correction=T
     !> [important]
     real(DP), pointer :: elda(:,:,:) !> to store DFT band energies
     SCALAR, pointer :: cg(:,:,:,:)
  end type wfnkqmpiinfo

  ! ======
  ! used to store outerband wavefunction in HDF5 case
  type OuterWfnInfo
     integer, pointer ::  ngk(:)
     integer, pointer :: isort(:,:)
     integer, pointer :: band_index(:,:)
     real(DP), pointer :: qk(:,:)
     real(DP), pointer :: el(:,:,:)
     SCALAR, pointer :: cg(:,:,:,:) ! ig, ib, is, ik
  end type OuterWfnInfo
  ! ------

  !---------------------------

  type wfnkmpiinfo
     integer, pointer ::  nkptotal(:)
     integer, pointer :: isort(:,:)
     real(DP), pointer :: qk(:,:)
     real(DP), pointer :: el(:,:,:)
     real(DP), pointer :: elda(:,:,:)
     ! ======
     ! cg(ibig,is,ik)
     SCALAR, pointer :: cg(:,:,:)
     ! ------
     ! cg(ig,ib,is,ik)
     ! SCALAR, pointer :: cg(:,:,:,:)
  end type wfnkmpiinfo

  !---------------------------

  type wpgen
     real(DP) :: wpsq(2) !< square of free el plasma freq for each spin
     real(DP) :: nelec(2) !< number of electrons for each spin per cell
     SCALAR, pointer :: rho(:,:) !< density, (ig, ispin)
  end type wpgen

  !---------------------------

  type polarizability
     ! logical :: no_screening
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
     ! real(DP) :: dFreqStepIncrease
     ! real(DP) :: dFreqCutoff1
     ! real(DP) :: dFreqCutoff2
     real(DP) :: dFreqCutoff     
     ! real(DP) :: omegappoverDeltaE !> determine the largest imaginary frequency
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

     type(scissors_t) :: scis
     logical :: eqp_corrections !< are we using eqp.dat and eqp_q.dat files
     complex(DPC), pointer :: dFreqBrd(:)  !< Corresponding Broadenings for Full Frequency
     integer :: fullConvLog !< logging pol matrix head & tail convergence
     ! integer :: iwritecoul !< flag to write vcoul
     integer :: nmtx
     integer, pointer :: nmtx_of_q(:)
     integer :: qgrid(3)
     ! integer, pointer :: qflags(:)
     integer :: nq0, nq1, nq !< Number of q->0 points, q/=0 points, and total number of q-points
     logical :: subsample !< whether we have more than one q0 point (used in subsampled calculation)
     logical :: non_uniform !< do non-uniform sampling using Voronoi decomposition of BZ
     integer :: gcomm
     logical :: min_fftgrid   !< use the smallest possible fftbox
     ! FHJ: These flags control some experimental optimizations
     integer :: os_opt_ffts       !< optimizes calculation/reuse of FFTs (real-space WFNs)
     ! integer :: nfreq_group     !< num. of frequencies to calculate in parallel
     ! integer :: nfreq_in_group  !< num. of epsilon frequencies held by any processor
     ! integer :: os_nsfreq_para !< num. of spectral frequencies held by any processor
     
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
     ! integer :: iqexactlyzero !< 1 if the q->0 point is *exactly* zero and will be read from WFN; 0 otherwise
     ! integer :: valueq0       !< 1=semiconductor (read from WFNq); 2=metal (read from WFN)
     integer, pointer :: irow(:)
     integer, pointer :: isrtx(:)
     integer, pointer :: isrtxi(:)
     integer :: icutv               !< icutv encodes presence and type of truncation
     real(DP) :: truncval(3)   !< in Bohr (au)
     real(DP), pointer :: qpt(:,:)
     !> FHJ: gme = <c,k|e^(-i(q+G).r)|v,k+q>, and the indices are:
     !! (nmtx, ncownactual, nvownactual, nspin, nrk, nfreq_group)
     ! SCALAR, pointer :: gme(:,:,:,:,:,:)
     ! SCALAR, allocatable :: gme(:,:,:,:,:,:)     
     SCALAR, allocatable :: gme(:,:,:,:,:)
     SCALAR, allocatable :: gme2(:,:,:,:)
     ! SCALAR, pointer :: chi(:,:,:)
     SCALAR, allocatable :: chi(:,:,:)

     integer :: ncrit
     real(DP) :: efermi
     real(DP) :: efermi_input
     logical :: rfermi
     real(DP) :: ecuts    !< energy cutoff of screened coulomb interaction in Ry
                          !> Reference regarding retarded/advanced functions: Catalin`s thesis, Eq. (1.44)
     ! complex(DPC), pointer :: chiRDyn(:,:,:,:) !< Retarded polarizability
     ! complex(DPC), pointer :: chiTDyn(:,:,:,:) !< Spectral function of polarizability
     complex(DPC), allocatable :: chiRDyn(:,:,:,:) !< Retarded polarizability
     complex(DPC), allocatable :: chiTDyn(:,:,:,:) !< Spectral function of polarizability
     
     ! real(DP), pointer :: edenDyn(:,:,:,:,:) !< Dynamic energy denominator     
     real(DP), pointer :: edenDyn(:,:,:,:) !< Dynamic energy denominator     
     real(DP), pointer :: edenDyn2(:,:,:) !< Dynamic energy denominator     
     logical :: degeneracy_check_override
     real(DP) :: lin_denominator !< energy threshold below which to activate lin_denominator
     type(cvpair_info), pointer :: lin_edenDyn(:,:,:,:,:) !< energies and
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

  !--------------------------------

  type cvpair_info
     real(DP) :: vc(2) !< conduction band velocity
     real(DP) :: vv(2) !< valence band velocity
     real(DP) :: ec !< conduction band energy
     real(DP) :: ev !< valence band energy
     integer :: idx_kp !< kpoint index
     logical :: vltc !< ev - ec < TOL_Degeneracy
  end type cvpair_info

  !--------------------------------

  type wfnkstates
     integer :: nkpt
     integer :: ndv
     integer, pointer :: isrtk(:)
     real(DP) :: k(3)
     real(DP), pointer :: ek(:,:)
     real(DP), pointer :: elda(:,:)
     SCALAR, pointer :: zk(:,:)
  end type wfnkstates

  !---------------------------

  type wfnkqstates
     integer :: nkpt
     integer, pointer :: isrtkq(:)
     real(DP), pointer :: ekq(:,:)
     SCALAR, pointer :: zkq(:,:)
  end type wfnkqstates

  !---------------------------------

  !> Used in haydock/diag only (see epsdiag.f90)
  type epsinfo
     integer :: nq !< number of q-vectors stored
     real(DP) :: emax !< maximum length of the stored q-vectors
     real(DP), pointer :: q(:,:) !< (3, nq) coordinates of q-vectors
     real(DP), pointer :: eps(:) !< (nq) head of dielectric matrix at each q-vector
     SCALAR :: epshead !< head of dielectric matrix at q->0
     real(DP) :: q0vec(3) !< coordinates of the q->0 vector
  end type epsinfo

  !------------------------------------

  !> Used in haydock/diag only
  !! Note that the bands in eqpv/eqpc are indexed with respect to the Fermi
  !! level, i.e., eqpv(1,:,:) is the VBM, eqpv(2,:,:) is VMB-1, etc.
  type eqpinfo
     type(scissors_t) :: scis
     type(spline_tck) :: spl_tck       !< scissors spline coefficients
     real(DP), pointer :: evqp(:,:,:)
     real(DP), pointer :: ecqp(:,:,:)
     real(DP), pointer :: evqp_co(:,:,:)
     real(DP), pointer :: ecqp_co(:,:,:)
     real(DP), pointer :: evqp_co_q(:,:,:)
     real(DP), pointer :: ecqp_co_q(:,:,:)
     real(DP), pointer :: evlda(:,:,:)
     real(DP), pointer :: eclda(:,:,:)
     real(DP), pointer :: evlda_co(:,:,:)
     real(DP), pointer :: eclda_co(:,:,:)
     real(DP), pointer :: evlda_co_q(:,:,:)
     real(DP), pointer :: eclda_co_q(:,:,:)
     real(DP), pointer :: evshift(:,:,:)
     real(DP), pointer :: ecshift(:,:,:)
     real(DP), pointer :: evshift_co(:,:,:)
     real(DP), pointer :: ecshift_co(:,:,:)
     real(DP), pointer :: evshift_co_q(:,:,:)
     real(DP), pointer :: ecshift_co_q(:,:,:)
  end type eqpinfo

  !------------------------------------

  !> moments for Haydock
  type mmtsinfo
     integer :: nmax
     integer :: nmaxp
     real(DP) :: norm
     real(DP) :: vol
     real(DP), pointer :: an(:)
     real(DP), pointer :: bn(:)
  end type mmtsinfo

  !------------------------------------

  type xctinfo     
     logical :: half_kernel, output_chi, use_3D_exchange
     integer, allocatable :: ik2ikm_index(:), ikq2ikqm_index(:), nrq_loc(:), irq_l2g(:,:), irq_g2l(:)
     real(DP) :: q2min
     integer :: celltype
     logical :: vcoul_mc              !< If vcoul_mc = T, calculate vcoul using MC average
     logical :: long_range_exchange   !< If long_range_exchange=T, include the long-range exchange kernel
     logical :: is_absorption         !< whether we are running the absorption code
     integer :: algo                  !< algorithm to use to solve BSE. See Common/nrtype.f90
     logical :: inteqp                !< whether we are interpolating
     logical :: is_periodic(3)        !< which dimensions are periodic
     integer :: idimensions           !< how many total periodic dimensions
     integer :: nkpt_co               !< number of kpts in the coarse grid
     integer :: nvb_co                !< number of valence bands in the coarse grid
     integer :: ncb_co                !< number of conduction bands in the coarse grid
     integer :: n1b_co                !< nvb_co for TDA calculations, nvb_co + ncb_co for non-TDA
     integer :: n2b_co                !< ncb_co for TDA calculations, nvb_co + ncb_co for non-TDA
     integer :: nspin
     integer :: nspinor = 1           !< nspinor = 2 if doing two-component spinor calculation; 1 is default
     integer :: qflag                 !< =0 for finite Q calculation with arbitrary Q (deprecated)
     !< =1 for Q=0 calculation (default)
     !< =2 use Q commensurate with WFN_co kgrid (under construction)
     integer :: iscreen !< what type of screening is present. 0 = semiconductor, 1 = graphene, 2 = metal
     logical :: renorm_transf         !< renormalize the dcc/dvv interpolation transformation?
     !> Calculate kernel blocks other than (vc)->(v'c') transitions? This will
     !! even include transitions such as (e,e)->(e',e'). In principle, we should
     !! always include these blocks if we are interpolating later on, but they
     !! are typically not important for semiconductors within TDA.
     logical :: extended_kernel
     !> If true, we extend the co/fi transf. to WFNs of different characters:
     !! |v> = \sum_n` d_vn`|n`>, where |n`> can be |v`> or |c`>
     !! If false, we restrict the character of the expansion coefficients:
     !! |v> = \sum_v` d_vv`|v`>
     logical :: unrestricted_transf
     !> Zero out dvc/dcv coefficients
     logical :: zero_unrestricted_contrib
     logical :: patched_sampling      !< simplest case of non-uniform sampling. See absorption.inp.
     logical :: patched_sampling_co   !< Use non-uniform sampling for coarse grid for Kernel. See absorption.inp.
     logical :: tda                   !< use Tamm-Dancoff approximation? (Absorption only)
     logical :: zero_coupling_block   !< If true, zero hbse_b before calling p*bseig
     integer :: iabsorp0              !< 1 means noeh_only, 0 otherwise
     logical :: rpa_only = .false.
     integer :: iwriteint = 1         !< = 0 for comm_disk, = 1 for comm_mpi
     logical :: eqp_corrections       !< do we use eqp.dat and eqp_q.dat
     logical :: eqp_co_corrections    !< do we use eqp_co.dat
     logical :: eqp_co_q_corrections    !< do we use eqp_co_q.dat

     !> For Coulomb interaction truncation
     ! integer :: iwritecoul
     integer :: icutv              !< icutv encodes presence and type of truncation
     real(DP) :: truncval(3)       !< in Bohr (au)
     !! double integral truncated_factor
     logical :: use_hdf5        !< with -DHDF5, whether or not we actually use hdf5
     logical :: bNoComm         !< If this is true, each processor will store the entire epsilon matrix
     logical :: bLittleComm     !< If this is true, each processor will store its own epsilon matrix          
     logical :: bLowComm        !< If this is true, ROOT will store the entire epsilon matrix
     logical :: bMidComm        !< If this is true, each processor will store part of rq for the epsilon matrix
     logical :: bHighComm       !< If this is true, each processor will store columns of the epsilon matrix     
     logical :: delaunay_interp !< use Delaunay interpolation?
     integer :: neps            !< Number of G vectors to capture the dielectric cutoff
     integer :: ilowmem
     logical :: skipinterp
     integer :: ivpar, icpar
     integer :: nn             !< PlotXct restrict_kpoints
     integer :: ng             !> size of mccp, mvvp
     integer :: ng_x           !> size of mvc, mvpcp
     integer :: nktotal        !< total number of unit cells
     !> Number of vertices in co k-grid that are used to expand each k-point in
     !! the fine grid for the **kernel** interpolation. This is 1 for the
     !! greedy interpolation (previous behaviour of the code), and ndims+1
     !! if we are performing Delaunay interpolation.
     integer :: npts_intp_kernel
     real(DP) :: eta           !< energy resolution
     real(DP) :: sigma         !< (used to calculate the optical spectrum)
     real(DP) :: gamma         !< (used to calculate the optical spectrum)
     real(DP) :: qshift
     real(DP) :: shift(3)      !< shift vector (this is the small shift,
     !< used to generate WFNq_fi, referenced only if xct%read_kpoints)
     real(DP) :: finiteq(3)    !< center-of-mass momentum of exciton
     logical :: energy_loss   !< calculate energy loss spectrum
     real(DP) :: lpol          !< norm of pol
     real(DP) :: pol(3)        !< light polarization for transition matrix elements
     integer :: npol           !< number of polarizations we have. Either 1 or 3.
     integer :: nmtxmax        !< max. number of columns in epsmat or eps0mat
     integer :: theory         !< theory level in kernel calculation
     !< 0 - GW-BSE, 1 - TDDFT
     real(DP) :: q0vec(3)      ! This is a hack for passing q0vec for
     ! TDDFT calculations (never used otherwise)
     ! when there is no epsilon
     type(coulomb_modifier_t) :: coulomb_mod
     logical :: coul_mod_flag   !< Flag which tells if the coulomb interaction has been
     !< modified (for TDDFT calculations)
     integer, pointer :: indexq(:), indexq_fi(:), umk(:,:) !< When exciton has finite center-of-mass momentum,
     !< maps between valence states at k+Q and conduction states at k
     integer, pointer :: isrtqi(:,:), nmtxa(:)
     integer, pointer :: ifmax(:,:), ifmaxq(:,:)
     real(DP) :: ecute         !< energy cutoff used in dielectric matrix
     real(DP) :: ecutg         !< energy cutoff used in wavefunctions and interaction
     !< kernel, see Rohlfing & Louie, PRB 62(8),p. 4938
     !< (must be slightly longer than xct%ecute because of umklapp vectors)
     real(DP) :: efermi        !< computed efermi
     real(DP) :: efermi_input  !< as set in input file
     logical :: rfermi         !< relative or absolute Fermi level
     SCALAR, pointer :: epsdiag(:,:) !< (nmtxmax, nq+1)
     type (wpgen) :: wpg !< charge density/fxc for TDDFT
     ! FHJ: TODO - move the following quantities to a separate derived type,
     ! ir reuse epsmpi
     !> Regular comm: (nmtxmax, ngpown_max, nq+1). The local processor stores row ig
     !! and a "local column" igp_l from epsilon in epscol(ig, igp_l, ik).
     !! Low comm: (nmtxmax, nmtxmax, nq+1). Each processor stores all eps(0)mat.
     !! local column = global epsilon column.
     SCALAR, pointer :: epscol(:,:,:)
     SCALAR :: epshead
     integer :: ngpown_max !< max. number of eps columns a PE can have
     integer :: ngpown !< number of columns of epsinv I own, for largest matrix
     integer :: nb !< block size for column distribution

     integer :: nrqown, nrqown_max, nb_rq     
     ! The arrays epsown and epsowni were removed, as they can be calculated on the fly:
     ! iowner = INDXG2P(igp, xct%nb, peinf%inode, 0, peinf%npes)
     ! xct%epsown(igp) = INDXG2P(igp, xct%nb, peinf%inode, 0, peinf%npes)
     ! xct%epsowni(igp_loc, iowner+1) = INDXG2L(igp, xct%nb, peinf%inode, 0, peinf%npes)
     logical :: output_hbse
     real (DP), dimension(3) :: finiteq_cart_unit
     ! logical :: long_range_exchange
     logical :: calcM, calcChi, calcAbs
     
     !> Used in haydock/diag only
     integer :: nkpt_fi         !< number of kpts in the fine grid
     integer :: nkptq_fi        !< number of kpts in the q-shifted fine grid
     integer :: nvb_fi          !< number of valence bands in the fine grid
     integer :: ncb_fi          !< number of conduction bands in the fine grid
     real(DP) :: avgcut         !< Cut in which we do cell averages on
     real(DP) :: wplasmon
     integer :: vmin,vmax
     integer :: rgrid(3) !< regular grid used to calculate qpt_averages (default is kgrid)
     logical :: freplacebz
     logical :: fwritebz
     logical :: degeneracy_check_override
     logical :: die_outside_sphere
     logical :: averagew
     logical :: subsample_line    !< during kernel interpolation, replace interpolated matrix
     !< elements with matrix elements from a precalculated bsemat file
     !< on a subsampled grid
     real(DP) :: subsample_cutoff !< use subsampled BSE matrix elements when |q| is less than cutoff
     real(DP) :: delta_frequency  !< Frequency step for absorption spectrum

     logical :: save_memory
     ! logical :: speed_eps
     logical :: zero_vq0 !< Use V_G(q) = 1/|q+G|^2
     integer :: diag_algo !> = 0: Use pdsyevx/pzheevx, parallel bisection followed by inverse iteration
                          !> = 1: Use pdsyevr/pzheevr, multiple relatively robust representations
                          !> = 2: Use pdsyevd/pzheevd, divide and conquer
     integer :: ng_in_d   !> Number of different G vectors used in constructing D array for interpolation
     integer, allocatable :: gindex_in_d_list(:) !> Index of those Gvectors in gvec%components(1:3,1:gvec%ng)
     logical :: eps_use_sym
     logical :: use_sym_wfn_co
     logical :: use_sym_wfnmq_co
     logical :: use_sym_wfn_fi     
     logical :: use_sym_wfnmq_fi
     integer :: kernel_type
     integer :: neigenvector
     logical :: save_memory_mdat
     real(DP) :: vcoul0, vcoul0_q0 !> store MC averaged vcoul(q=G=0) calculated from calc_vcoul0
     SCALAR :: wcoul0, wcoul0_q0
     ! real(DP) :: wcoul0 !> store MC averaged epsinv(q=G=0)*vcoul(q=G=0) calculated from calc_wcoul0
     integer :: qgrid(3)
     integer :: broadeningtype !> 0: lorentzian, 1: gaussian, 2: voigt
     integer :: niter !> number of iterations in Lanczos algorithm
     logical :: use_elpa
  end type xctinfo

  !------------------------------------

  type flags
     !>
     !> Used in haydock, diag, nonlinearoptics
     !>
     !> Defined flags:
     !>
     !>  eig = 0  --> do not write eigenvectors (default)
     !>      < 0 --> write all eigenvectors
     !>      > 0 --> write the first flag%eig eigenvectors
     !>
     !>  read_epsdiag = false --> read files 'eps0mat'/'epsmat' (default)
     !>                 true  --> read file 'epsdiag.dat'
     !>
     !>  krnl = 0 --> spin triplet kernel, direct kernel only (only allow for nspin = 1)
     !>         1 --> spin singlet kernel (default)
     !>         2 --> local-fields + RPA, exchange kernel only
     !>         3 --> spinor kernel
     !>
     !>  opr = 0 --> use velocity operator
     !>        1 --> use momentum operator
     !>        2 --> use JDOS operator (Haydock only)
     !>
     !>  lor = 0 --> use Lorentzian broadening
     !>        1 --> use Gaussian broadening
     !>        2 --> use Voigt broadening
     !>
     !>  spec = 0 --> go through the whole exciton calculation (default)
     !>         1 --> calculate only absorption spectrum (this option skips
     !>               all calculation and goes right to the end of the code)
     !>
     !>  vm = 0 --> calculate velocity/momentum matrix elements (default)
     !>       1 --> read velocity/momentum matrix elements from file vmtxel
     !>       2 --> use vectors from previous iteration (Haydock only!)
     !>
     !>  job = 0 --> ultrafast calculation
     !>  job = 1 --> two-photon calculation
     !>
     integer :: lor
     integer :: eig
     integer :: krnl
     integer :: opr
     integer :: job
     !> Use averaged Gauss quadrature in Lanczos algorithm? Default is true.
     logical :: lanczos_gauss_quad
     logical :: debug_lanczos    
     ! !> when we read from vmt.dat, we should set flag%calculate_vmt = .false.
     ! logical :: calculate_vmt
     logical :: bare_vmt     
  end type flags

  !---------------------------------

  type otherinfo
     integer :: ithreeD
     integer :: knx
     integer :: kny
     integer :: knz
     real(DP) :: keta
  end type otherinfo

  !---------------------------------

  type windowinfo
     real(DP) :: evalue
     real(DP) :: emax
     real(DP), pointer :: cstates(:)
     real(DP), pointer :: estates(:)
     integer, pointer :: istates(:)
     integer :: nstates
  end type windowinfo

  !-----------------------------
  !> coarse-grid wavefunctions for diag/haydock
  type tdistgwf

     integer :: block_sz !< block size for BLACS distribution = DIVUP(ngm, npes)
     integer :: ngm !< maximum number of G-vectors = kp*%ngkmax
     integer :: ngl !< local number of G-vectors = NUNROC(...). At most block_sz.
     integer :: tgl !< local to global translation = block_sz * peinf%inode

     !> local to global index translation : ig_g = ig_l + tgl
     !> ig_g = 1 ... ng(ik) is the global G-index
     !> ig_l = 1 ... ngl is the local G-index

     integer :: nk !< number of k-points
     integer :: ns !< number of spin components
     integer :: nspinor = 1 !< nspinor = 2 if doing two-component spinor calculation; 1 is default
     integer :: nv !< number of valence bands
     integer :: nc !< number of conduction bands

     integer, pointer :: ng(:) !< (nk)
     integer, pointer :: isort(:,:) !< (ngl,nk)
     SCALAR, pointer :: zv(:,:,:,:) !< (ngl,nv,ns*nspinor,nk)
     SCALAR, pointer :: zc(:,:,:,:) !< (ngl,nc,ns*nspinor,nk)

  end type tdistgwf

  !-----------------------------
  !> MJ: work arrays - getting rid of save statements

  type work_genwf
     integer :: ikold = 0
     integer :: nb
     integer :: ng
     integer :: ns
     integer :: nspinor = 1 !< nspinor = 2 if doing two-component spinor calculation; 1 is default
     ! ======
     integer :: nspin
     ! ------
     SCALAR, pointer :: cg(:,:,:)
     SCALAR, pointer :: ph(:)
     integer, pointer :: ind(:)
     integer, pointer :: isort(:)
  end type work_genwf

  !-----------------------------
  !> (gsm) work arrays - getting rid of save statements

  type twork_scell
     integer :: dNfft(3)
     integer :: Nplane
     integer :: Nrod
     complex(DPC), pointer :: fftbox_1D(:,:)
  end type twork_scell

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

  !> FHJ: header information for kernel files (bsedmat, bsexmat, bsemat.h5)
  type kernel_header_t

     ! Mean-field and other general information
     type(mf_header_t) :: mf !< mf header containing number of k-points, WFN cutoff, etc.

     integer :: version
     integer :: iflavor

     integer :: iscreen !< screening flag
     integer :: icutv   !< truncation flag
     real(DP) :: ecuts  !< epsilon cutoff
     real(DP) :: ecutg  !< WFN cutoff
     real(DP) :: efermi !< Fermi energy found by the code after any shift
     integer :: theory  !< 0 for GW-BSE, 1 for TD-HF, 2 for TD-DFT
     !> How many transitions blocks are there in the kernel matrix?
     !! 1 for restricted TDA kernel: vc -> v`c`
     !! 2 for restricted non-TDA kernel: {vc,cv} -> {v`c`,c`v`}  [not implemented]
     !! 4 for extended kernel: {n1,n2} -> {n1`,n2`}
     integer :: nblocks
     integer :: storage !< 0 if storing full matrix (only option supported now)
     integer :: nmat    !< number of matrices in the file (1 for bsexmat, 3 for bsedmat)
     logical :: energy_loss !< is this an energy-loss calculation?

     integer :: nvb     !< number of valence bands in the coarse grid
     integer :: ncb     !< number of conduction bands in the coarse grid
     integer :: n1b     !< nvb_co if kernel_sz==1; nvb_co + ncb_co if kernel_sz=4
     integer :: n2b     !< ncb_co if kernel_sz==1; nvb_co + ncb_co if kernel_sz=4
     integer :: ns      !< number of spins
     integer :: nspinor = 1 !< nspinor = 2 if doing two-component spinor calculation; 1 is default
     logical :: patched_sampling !< are we doing a calculation on a patch?

     ! Variables specific to kernel files
     integer :: nk      !< number of k-points
     real(DP), pointer :: kpts(:,:)
     integer :: kgrid(3)
     !> 0 for finite Q calculation with arbitrary Q (deprecated)
     !! 1 for Q=0 calculation (default)
     !! 2 use Q commensurate with WFN_co kgrid (under construction)
     integer :: qflag
     real(DP) :: center_mass_q(3)

  end type kernel_header_t

end module typedefs_m

!!================================================================================
!!
!! Modules:
!!
!! (1) nrtype_m      Originally By ?      Last Modified 8/22/2010 (gsm)
!!
!!     Global constants and parameters.
!!
!!================================================================================

#include "f_defs.h"

#define AUTO_VER(x) VER_ ## x ## _HDF5

module nrtype_m

  implicit none

  public ! only public parameters here

  ! Below are the version number for each BerkeleyGW file. These numbers should
  ! be changed whenever the structure of the file is altered and there`s either
  ! incompatibility with the previous version or a new feature. The version
  ! should be -1 if that file is not versioned yet, which corresponds to
  ! the formats used in the Berkeley 1.0.x family.
  integer, parameter :: VER_WFN_FORT = -1
  integer, parameter :: VER_WFN_HDF5 = 1
  integer, parameter :: VER_WFN = AUTO_VER(WFN)
  integer, parameter :: VER_EPS_FORT = -1
  integer, parameter :: VER_EPS_HDF5 = 3
  integer, parameter :: VER_EPS = AUTO_VER(EPS)
  integer, parameter :: VER_BSE_FORT = 1
  integer, parameter :: VER_BSE_HDF5 = 2
  integer, parameter :: VER_BSE = AUTO_VER(BSE)

  !! Maximum number of bands supported by the *inread* routines. This sets the
  !! size of arrays such as "occupations". These arrays should all be allocated
  !! dynamically in the future.
  integer, parameter :: MAX_BANDS = 100000 ! "occupations" array => 7MB
  !! Maximum number of {k,q}-points supported by the *inread* routines.
  !! The actual number of k-points/q-points in the WFN/bsemat/epsmat files
  !! can be larger.
  integer, parameter :: MAX_KPTS = 100000 ! "kpt_read" array => 0.8 MB

  !! parameters for real-space resolution in cell-truncation schemes
  integer, parameter :: n_in_box = 2
  integer, parameter :: n_in_wire = 4

  !! parameter for construction of Wigner-Seitz cell
  integer, parameter :: ncell = 3

  !! number of Monte-Carlo integration points
  integer, parameter :: nmc_coarse = 250000
  integer, parameter :: nmc_fine = 2500000
  integer, parameter :: nmc = nmc_fine

  !! type definitions following the convention of Numerical Recipes
  !! do not ever use single-precision!!
  !  integer, parameter :: SP = kind(1.0)
  integer, parameter :: DP = kind(1.0d0)
  !  integer, parameter :: SPC = kind((1.0,1.0))
  integer, parameter :: DPC = kind((1.0d0,1.0d0))

  !! a shift on the grid in order to avoid the singularity for truncation
  real(DP), parameter :: trunc_shift(3) = (/0.5d0, 0.5d0, 0.5d0/)

  !! physical constants
  !!
  !! These are the "2010 CODATA recommended values" taken from
  !! "The NIST Reference on Constants, Units, and Uncertainty"
  !! http://physics.nist.gov/cuu/
  !!
  !! The following variables are used throughout the package:
  !!     'BOHR', 'bohr' is Bohr radius, in Angstrom
  !!     'RYD', 'ryd2eV', 'rydberg' is Rydberg constant times hc, in eV
  !!     'HARTREE', 'hartree' is Hartree energy, in eV
  !!     'LIGHTSPEED' is inverse alpha (fine-structure constant)
  !!
  real(DP), parameter :: BOHR = 0.52917721092_dp
  real(DP), parameter :: RYD = 13.60569253_dp
  real(DP), parameter :: LIGHTSPEED = 137.035999074_dp

  !! Mathematical constants
  !! real(SP), parameter :: PI_S = 3.1415926535897932384626433832795_sp
  real(DP), parameter :: PI_D = 3.1415926535897932384626433832795_dp
  !! https://en.wikipedia.org/wiki/Euler%E2%80%93Mascheroni_constant
  real(DP), parameter :: gamma_D = 0.5772156649015328606065120900824024310421_dp
  real(DP), parameter :: TOL_SMALL = 1.0D-6
  real(DP), parameter :: TOL_ZERO = 1.0D-12
  real(DP), parameter :: TOL_Degeneracy = 1.0D-6
  real(DP), parameter :: INF = 1.0d12

end module nrtype_m

#undef AUTO_VER

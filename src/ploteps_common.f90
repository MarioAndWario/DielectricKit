#include "f_defs.h"

module ploteps_common_m

  use global_m
  implicit none

  private

  public :: ploteps_t
  
  type ploteps_t
     character*20 :: filename         !< either epsmat.h5 or chimat.h5
     real(DP) :: epsinvhead           !< need to be set for filename = epsmat.h5
     integer :: nrq                   !< number of qpoints in RBZ
     real(DP), pointer :: rq(:,:)     !< qpoints in RBZ
     integer :: nfq                   !< number of qpoints in FBZ
     real(DP), pointer :: fq(:,:)     !< qpoints in FBZ
     integer, pointer :: nmtx(:)      !< number of Gvectors for each RBZ qpoint
     integer :: nmtx_max              !< nmtx_max = MAXVAL(nmtx(:))
     integer :: nmtx_used             !< nmtx_used <= MINVAL(nmtx(:))
     real(DP) :: ecut                 !< nmtx_used is determined from ecut
     integer :: qgrid(3)              !< qgrid
     integer :: ispin                 !< which spin component to plot
     integer :: nr2                   !< number of r2 points
     real(DP), allocatable :: r2(:,:) !< coordinate of r2 points, with respect to lattice vectors
     logical :: unfold                !< use symmetries to unfold grid
     integer :: downsample(3)         !< downsample each dimention of the real-space grid by this factor. Defaults to 2.
     integer :: nsuper(3)             !< number of unit-cell repetitions
     integer, pointer :: isrtx(:,:)   !< gvec%components(:,isrtx(ig)) is the ig-th Gvector in epsmat
     integer, pointer :: isrtxi(:,:)  !< ig_gvec = isrtx(ig), isrtxi(ig_gvec) = ig
     logical :: realpart              !< output real part of scfft 
     logical :: imagpart              !< output imaginary part of scfft
     logical :: abs                   !< output abs(scfft)    
     logical :: abs2                  !< output abs(scfft)**2

     !> Frequency dependence
     integer :: freq_dep        !< frequency dependence of the inverse dielectric matrix
                                ! 0: static calculation 2: full frequency 3: two imaginary frequencies
     integer :: nfreq           !< number of frequencies used in full frequency calculation
     integer :: nfreq_imag      !< number of imaginary freqs for CD (also 1 for GN GPP)
     ! complex(DP), pointer :: freqs(:) !< complex frequencies
     complex(DP) :: freq_target
     integer :: FFTgrid(3), FFTfactor
     logical :: high_resolution
     logical :: low_comm, mid_comm
     integer :: rq_blocksize, nrq_loc, nrq_loc_max
  end type ploteps_t

end module ploteps_common_m

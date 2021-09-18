!==============================================================================
!
! Module:
!
! (1) vcoul_generator()     Originally by JRD, modifed by Meng Wu
!
!     Generates the (Truncated) Coulomb Interaction for all G at a particular
!     q.  Outputs what would be 8Pi/q^2 if not for truncation.
!
!==============================================================================

#include "f_defs.h"

!! Truncation flag. The current supported options are:
!! 0: No truncation (for 3D systems)
!! 6: Slab truncation (for 2D systems)
module vcoul_generator_m
  use global_m
  ! use minibzaverage_m
  implicit none
  private :: length2
  public :: vcoul_generator, destroy_qran, destroy_qran_mBZ
  !! qran is fractional coordinates
  !! qran_shell_cart in cartesian coordinates
  real(DP), allocatable, private :: qran(:,:), random_vector_cart_mBZ(:,:)
  real(DP), private :: q2_inner = 0d0, q2_outer = 0d0
  integer, private :: ifirst = 1, size_random_vector_mBZ
  real(dP), allocatable, private :: vcoul1(:)
  logical, allocatable, private :: vcoul1_done(:)
contains
  !! Calculates the (average) bare Coulomb interaction v(q+G) in Rydberg.
  !! This function uses different forms of v(q+G) (see below) depending on the
  !! dimensionality of the system and whether the potential is truncated.

  !! ROOT calculate vcoul and BCAST to all procs
  subroutine vcoul_generator(itruncflag, gvec, crys, ncoul, isrtq, qvec, vcoul)
    !! Truncation flag. The current supported options are:
    !! 0: No truncation (for 3D systems)
    !! 6: Slab truncation (for 2D systems)
    integer, intent(in) :: itruncflag
    !! G space containing the G-vectors used to calculate v(q+G)
    type (gspace), intent(in) :: gvec
    type (crystal), intent(in) :: crys
    !! Number of G vectors to calculate v(q+G)
    integer, intent(in) :: ncoul
    !! (ncoul) Indices of the G vectors from gvec which we are using to build v(q+G)
    integer, intent(in) :: isrtq(:)
    !! The q vector that we use to calculate v(q+G).
    !! The null vector q=(/0,0,0/) is perfectly valid here.
    real(DP), intent(in) :: qvec(3)
    !! (ncoul) The output, 8*pi*<1/(q+G)^2> if there`s no truncation.
    real(DP), intent(out) :: vcoul(:)
    integer :: ig
    real(DP) :: q_plus_G(3), q_plus_G_xy(3), qlen, qG2, kxy
    real(DP) :: kz, half_Lz, q_plus_G_z(3)
    logical :: verbose
    real(DP) :: U_(3,3)
    integer :: info_
    PUSH_SUB(vcoul_generator)

    verbose = peinf%verb_debug
    if(size(isrtq) < ncoul) call die("vcoul_generator: isrtq not allocated to size ncoul")
    if(size(vcoul) < ncoul) call die("vcoul_generator: vcoul not allocated to size ncoul")
    U_(1:3, 1:3) = crys%bdot(1:3, 1:3)
    call dpotrf('U', 3, U_, 3, info_)
    vcoul = 0.0D0
    !! <- Epsilon: for q0, qlen = |q0| > 0 ->
    !! for q0 vectors, qlen = 0
    ! qlen = sqrt(DOT_PRODUCT(qvec, MATMUL(crys%bdot, qvec)))
    qlen = SQRT(length2(qvec, U_))

    !! No Truncation
    if (itruncflag .eq. 0 ) then
       if (peinf%inode .eq. 0) then
          !! Since we order |rq+G|^2 using isrtq,
          !! rq+G(isrtq(1)) must has the smallest |rq+G|^2
          q_plus_G(:) = DBLE(gvec%components(:, isrtq(1))) + qvec(:)
          ! qG2 = DOT_PRODUCT(q_plus_G, MATMUL(crys%bdot, q_plus_G))
          qG2 = length2(q_plus_G, U_)
          if (qG2 .lt. Tol_zero) then
             vcoul(1) = 0.0D0
          else
             vcoul(1) = 8.0D0 * PI_D / qG2
          endif
          !! ncoul = max(ncouls,ncoulb)
          do ig = 2, ncoul
             q_plus_G(:) = DBLE(gvec%components(:, isrtq(ig))) + qvec(:)
             ! qG2 = DOT_PRODUCT(q_plus_G, MATMUL(crys%bdot, q_plus_G))
             qG2 = length2(q_plus_G, U_)
             vcoul(ig) = 8.0D0 * PI_D / qG2
          enddo
       endif
#ifdef MPI
       !! call MPI_ALLREDUCE(MPI_IN_PLACE,vcoul,ncoul,MPI_REAL_DP,MPI_SUM,MPI_COMM_WORLD,mpierr)
       if (peinf%npes > 1) then
          call MPI_Bcast(vcoul, ncoul, MPI_REAL_DP, 0, MPI_COMM_WORLD, mpierr)
       endif
#endif
       !! Slab truncation
    elseif (itruncflag .eq. 6) then
       if (ABS(qvec(3)) .gt. Tol_ZERO) then
          write(0,'(a)') 'You asked for cell slab truncation but have more q-points in z direction than qz=0!'
          call die('Bad Truncation', only_root_writes = .true.)
       endif

       if (peinf%inode .eq. 0) then
          q_plus_G(:) = DBLE(gvec%components(:, isrtq(1))) + qvec(:)
          ! qG2 = DOT_PRODUCT(q_plus_G, MATMUL(crys%bdot, q_plus_G))
          qG2 = length2(q_plus_G, U_)
          !! (G+q) in xy plane
          q_plus_G_xy(1:2) = q_plus_G(1:2)
          q_plus_G_xy(3) = 0.0D0
          !! kxy = |q+G|_{\parallel}
          kxy = SQRT(length2(q_plus_G_xy, U_))
          q_plus_G_z(1:2) = 0.0D0
          q_plus_G_z(3) = q_plus_G(3)
          !! kz = |q+G|_z = |G|_z, since we require that qvec(3) = 0
          kz = SQRT(length2(q_plus_G_z, U_))
          if (kz > TOL_ZERO) then
             call die("kz > 0", only_root_writes=.true.)
          endif
          half_Lz = SQRT(crys%adot(3,3)) / 2.0D0

          if (qG2 .lt. TOL_ZERO) then
             vcoul(1) = 0.0D0
             !<- q0 + G/=0 or q1 + any G ->
             !! no random integral method here
          else
             !! here we assume R=Lz/2 as the slab truncation length
             !! Also note that here we don't use random-integration method to calculate <- q0 + G/= 0 or q1 + any G ->
             vcoul(1) = 8.0D0 * PI_D / qG2 * (1.0D0 - EXP(-kxy*half_Lz)*COS(kz*half_Lz))
          endif

          do ig = 2, ncoul
             q_plus_G(:) = DBLE(gvec%components(:, isrtq(ig))) + qvec(:)
             ! qG2 = DOT_PRODUCT(q_plus_G, MATMUL(crys%bdot, q_plus_G))
             qG2 = length2(q_plus_G, U_)
             !! (G+q) in xy plane
             q_plus_G_xy(1:2) = q_plus_G(1:2)
             q_plus_G_xy(3) = 0.0D0
             !! kxy = |q+G|_{\parallel}
             kxy = SQRT(length2(q_plus_G_xy, U_))
             q_plus_G_z(1:2) = 0.0D0
             q_plus_G_z(3) = q_plus_G(3)
             !! kz = |q+G|_z = |G|_z, since we require that qvec(3) = 0
             kz = SQRT(length2(q_plus_G_z, U_))
             ! half_Lz = 2D0*PI_D/(sqrt(crys%bdot(3,3))*2D0)
             vcoul(ig) = 8.0D0 * PI_D / qG2 * (1.0D0 - EXP(-kxy*half_Lz)*COS(kz*half_Lz))
          enddo
       endif
#ifdef MPI
       if (peinf%npes > 1) then
          call MPI_Bcast(vcoul, ncoul, MPI_REAL_DP, 0, MPI_COMM_WORLD, mpierr)
       endif
#endif
    endif

    POP_SUB(vcoul_generator)
    return
  end subroutine vcoul_generator

  !-----------------------------------------------------------------
  subroutine destroy_qran()

    PUSH_SUB(destroy_qran)

    ifirst = 2
    SAFE_DEALLOCATE(qran)

    call logit('Deallocated random numbers.')

    POP_SUB(destroy_qran)
    return
  end subroutine destroy_qran

  subroutine destroy_qran_mBZ()
    PUSH_SUB(destroy_qran_mBZ)
    ifirst = 2
    ! if (peinf%inode .eq. 0) then
    SAFE_DEALLOCATE(random_vector_cart_mBZ)
    SAFE_DEALLOCATE(vcoul1)
    ! endif
    SAFE_DEALLOCATE(vcoul1_done)
    call logit('Deallocated random numbers in mBZ.')
    POP_SUB(destroy_qran_mBZ)
    return
  end subroutine destroy_qran_mBZ

  !! Based on uniform random numbers between [0,1], generate uniform random points (in cartesian coordinates) within a spherical shell,
  !! defined by the inner radius r_inner and outer radius r_outer
  !! If r_inner = 0, we have uniform distributed points within the entire sphere with radius r_outer
  subroutine spherical_uniform_distribution(nn, ran, r_inner, r_outer, num_random_points, random_points_cart, if_antithetic)
    real(DP), intent(in) :: r_inner, r_outer
    integer, intent(in) :: nn
    real(DP), intent(in) :: ran(3,nn)
    integer, intent(out) :: num_random_points
    real(DP), intent(out) :: random_points_cart(:,:)
    logical, intent(in), optional :: if_antithetic
    logical :: if_antithetic_ = .false.
    real(DP), allocatable :: theta(:), phi(:), r(:), r_antithetic(:)
    PUSH_SUB(spherical_uniform_distribution)

    if (present(if_antithetic)) then
       if_antithetic_ = if_antithetic
    endif

    SAFE_ALLOCATE(theta, (nn))
    SAFE_ALLOCATE(phi, (nn))
    SAFE_ALLOCATE(r,(nn))
    theta = 0.0D0
    phi = 0.0D0
    r = 0.0D0

    if (if_antithetic_) then
       num_random_points = 2 * nn
       SAFE_ALLOCATE(r_antithetic,(nn))
       r_antithetic = 0.0D0
    else
       num_random_points = nn
    endif

    if (size(random_points_cart,DIM=2) .ne. num_random_points) then
       call die("Check size of random_points_cart", only_root_writes=.true.)
    endif

    ! SAFE_ALLOCATE(random_points_cart, (3,num_random_points))
    random_points_cart = 0.0D0
    !! Spherical coordinates
    theta(:) = 2.0D0 * PI_D * ran(1,:)
    phi(:) = ACOS(2.0D0 * ran(2,:)-1.0D0)
    r(:) = (r_inner**3 + ran(3,:) * ( r_outer**3 - r_inner**3 ) )**(1.0D0/3.0D0)
    !! cartesian coordinates
    !! x
    random_points_cart(1,1:nn) = r(:) * SIN(phi(:)) * COS(theta(:))
    !! y
    random_points_cart(2,1:nn) = r(:) * SIN(phi(:)) * SIN(theta(:))
    !! z
    random_points_cart(3,1:nn) = r(:) * COS(phi(:))
    if (if_antithetic_) then
       r_antithetic(:) = (r_inner**3 + (1.0D0-ran(3,:)) * ( r_outer**3 - r_inner**3 ) )**(1.0D0/3.0D0)
       !! cartesian coordinates
       !! x
       random_points_cart(1,nn+1:2*nn) = r_antithetic(:) * SIN(phi(:)) * COS(theta(:))
       !! y
       random_points_cart(2,nn+1:2*nn) = r_antithetic(:) * SIN(phi(:)) * SIN(theta(:))
       !! z
       random_points_cart(3,nn+1:2*nn) = r_antithetic(:) * COS(phi(:))
    endif
    SAFE_DEALLOCATE(theta)
    SAFE_DEALLOCATE(phi)
    SAFE_DEALLOCATE(r)
    if (if_antithetic_) then
       SAFE_DEALLOCATE(r_antithetic)
    endif

    POP_SUB(spherical_uniform_distribution)
    return
  end subroutine spherical_uniform_distribution

  !! Based on uniform random numbers between [0,1], generate uniform random points (in cartesian coordinates) within a circular annulus,
  !! defined by the inner radius r_inner and outer radius r_outer
  !! If r_inner = 0, we have uniform distributed points within the entire circle with radius r_outer
  subroutine circular_uniform_distribution(nn, ran, r_inner, r_outer, num_random_points, random_points_cart, if_antithetic)
    real(DP), intent(in) :: r_inner, r_outer
    integer, intent(in) :: nn
    real(DP), intent(in) :: ran(2,nn)
    integer, intent(out) :: num_random_points
    real(DP), intent(out) :: random_points_cart(:,:)
    logical, intent(in), optional :: if_antithetic
    logical :: if_antithetic_ = .false.
    real(DP), allocatable :: theta(:), r(:), r_antithetic(:)
    PUSH_SUB(circular_uniform_distribution)

    if (present(if_antithetic)) then
       if_antithetic_ = if_antithetic
    endif

    SAFE_ALLOCATE(theta, (nn))
    SAFE_ALLOCATE(r,(nn))
    theta = 0.0D0
    r = 0.0D0

    if (if_antithetic_) then
       num_random_points = 2 * nn
       SAFE_ALLOCATE(r_antithetic,(nn))
       r_antithetic = 0.0D0
    else
       num_random_points = nn
    endif

    if (size(random_points_cart,DIM=2) .ne. num_random_points) then
       call die("Check size of random_points_cart", only_root_writes=.true.)
    endif

    random_points_cart = 0.0D0
    !! Spherical coordinates
    theta(:) = 2.0D0 * PI_D * ran(1,:)
    r(:) = SQRT(r_inner**2 + ran(2,:) * ( r_outer**2 - r_inner**2 ) )
    !! cartesian coordinates
    !! x
    random_points_cart(1,1:nn) = r(:) * COS(theta(:))
    !! y
    random_points_cart(2,1:nn) = r(:) * SIN(theta(:))

    if (if_antithetic_) then
       r_antithetic(:) = SQRT(r_inner**2 + (1.0D0-ran(2,:)) * ( r_outer**2 - r_inner**2 ) )
       !! cartesian coordinates
       !! x
       random_points_cart(1,nn+1:2*nn) = r_antithetic(:) * COS(theta(:))
       !! y
       random_points_cart(2,nn+1:2*nn) = r_antithetic(:) * SIN(theta(:))
    endif
    SAFE_DEALLOCATE(theta)
    SAFE_DEALLOCATE(r)
    if (if_antithetic_) then
       SAFE_DEALLOCATE(r_antithetic)
    endif

    POP_SUB(circular_uniform_distribution)
    return
  end subroutine circular_uniform_distribution

  ! subroutine generate_random_number(ran_seed, ran)
  !   integer, intent(in) :: ran_seed(2)
  !   real(DP), intent(out) :: ran(:,:)
  !   PUSH_SUB(generate_random_number)

  !   call random_seed(put=ran_seed)
  !   !! ran will be used in two places:
  !   !! 1. used to generate random points in minibz for vcoul with nonzero |q+G|
  !   !! 2. used to generate random points in spherical shell for vcoul0
  !   call random_number(ran)

  !   POP_SUB(spherical_uniform_distribution)
  !   return
  ! end subroutine generate_random_number

  subroutine generate_random_number(ran_seed, ran)
    integer, intent(in) :: ran_seed(:)
    real(DP), intent(out) :: ran(:,:)
    PUSH_SUB(generate_random_number)

    call random_seed(put=ran_seed)
    !! ran will be used in two places:
    !! 1. used to generate random points in minibz for vcoul with nonzero |q+G|
    !! 2. used to generate random points in spherical shell for vcoul0
    call random_number(ran)

    POP_SUB(spherical_uniform_distribution)
    return
  end subroutine generate_random_number

  subroutine pick_points_within_mBZ(crys, qgrid, size_random_vector_pool, random_vector_pool, size_random_vector_mBZ, random_vector_cart_mBZ)
    type(crystal), intent(in) :: crys
    integer, intent(in) :: qgrid(3), size_random_vector_pool
    real(DP), intent(in) :: random_vector_pool(:,:) !< (3, size_random_vector_pool)
    integer, intent(out) :: size_random_vector_mBZ
    real(DP), intent(out) :: random_vector_cart_mBZ(:,:) !< (3, size_random_vector_pool)
    integer :: ii, ib_combination
    real(DP) :: b_combination_cart(3,13), b_frac_(3)
    integer :: b_aux_(3,13)
    PUSH_SUB(pick_points_within_mBZ)

    b_aux_(:,1) = (/ -1,-1,-1 /)
    b_aux_(:,2) = (/ -1,-1, 0 /)
    b_aux_(:,3) = (/ -1,-1, 1 /)
    b_aux_(:,4) = (/ -1, 0,-1 /)
    b_aux_(:,5) = (/ -1, 0, 0 /)
    b_aux_(:,6) = (/ -1, 0, 1 /)
    b_aux_(:,7) = (/ -1, 1,-1 /)
    b_aux_(:,8) = (/ -1, 1, 0 /)
    b_aux_(:,9) = (/ -1, 1, 1 /)
    b_aux_(:,10) =(/  0,-1,-1 /)
    b_aux_(:,11) =(/  0,-1, 0 /)
    b_aux_(:,12) =(/  0,-1, 1 /)
    b_aux_(:,13) =(/  0, 0,-1 /)

    do ib_combination = 1, 13
       b_frac_(1) = DBLE(b_aux_(1,ib_combination))/DBLE(qgrid(1))
       b_frac_(2) = DBLE(b_aux_(2,ib_combination))/DBLE(qgrid(2))
       b_frac_(3) = DBLE(b_aux_(3,ib_combination))/DBLE(qgrid(3))
       b_combination_cart(:,ib_combination) = 2.0D0 * PI_D / crys%alat * MATMUL(crys%bvec, b_frac_)
    enddo

    !! size_random_vector_mBZ <= size_random_vector_pool
    random_vector_cart_mBZ = 0.0D0
    size_random_vector_mBZ = 0
    iqran_loop : do ii = 1, size_random_vector_pool
       do ib_combination = 1, 13
          if (ABS(DOT_PRODUCT(random_vector_pool(:,ii),b_combination_cart(:,ib_combination))) > 0.5D0 * NORM2(b_combination_cart(:,ib_combination))**2 ) then
             cycle iqran_loop
          endif
       enddo
       size_random_vector_mBZ = size_random_vector_mBZ + 1
       random_vector_cart_mBZ(:,size_random_vector_mBZ) = random_vector_pool(:,ii)
    enddo iqran_loop
    ! write(*,'(A,F10.6)') "size_random_vector_mBZ / size_random_vector_pool = ", DBLE(size_random_vector_mBZ) / DBLE(size_random_vector_pool)

    POP_SUB(pick_points_within_mBZ)
    return
  end subroutine pick_points_within_mBZ

  subroutine pick_q2_inner(crys, qgrid, q2_inner)
    type(crystal), intent(in) :: crys
    integer, intent(in) :: qgrid(3)
    real(DP), intent(out) :: q2_inner
    integer :: ib_combination, b_aux_(3,13)
    real(DP) :: b_frac_(3), q2_list(13) ! , b_combination_cart(3,13)
    PUSH_SUB(pick_q2_inner)

    b_aux_(:,1) = (/ -1,-1,-1 /)
    b_aux_(:,2) = (/ -1,-1, 0 /)
    b_aux_(:,3) = (/ -1,-1, 1 /)
    b_aux_(:,4) = (/ -1, 0,-1 /)
    b_aux_(:,5) = (/ -1, 0, 0 /)
    b_aux_(:,6) = (/ -1, 0, 1 /)
    b_aux_(:,7) = (/ -1, 1,-1 /)
    b_aux_(:,8) = (/ -1, 1, 0 /)
    b_aux_(:,9) = (/ -1, 1, 1 /)
    b_aux_(:,10) =(/  0,-1,-1 /)
    b_aux_(:,11) =(/  0,-1, 0 /)
    b_aux_(:,12) =(/  0,-1, 1 /)
    b_aux_(:,13) =(/  0, 0,-1 /)

    do ib_combination = 1, 13
       b_frac_(1) = DBLE(b_aux_(1,ib_combination))/DBLE(qgrid(1)) / 2.0D0
       b_frac_(2) = DBLE(b_aux_(2,ib_combination))/DBLE(qgrid(2)) / 2.0D0
       b_frac_(3) = DBLE(b_aux_(3,ib_combination))/DBLE(qgrid(3)) / 2.0D0
       q2_list(ib_combination) = DOT_PRODUCT(b_frac_, MATMUL(crys%bdot, b_frac_))
    enddo

    q2_inner = MINVAL(q2_list)

    POP_SUB(pick_q2_inner)
    return
  end subroutine pick_q2_inner

  subroutine pick_points_within_mBZ_2D(crys, qgrid, size_random_vector_pool, random_vector_pool, size_random_vector_mBZ, random_vector_cart_mBZ)
    type(crystal), intent(in) :: crys
    integer, intent(in) :: qgrid(3), size_random_vector_pool
    real(DP), intent(in) :: random_vector_pool(:,:) !< (2, size_random_vector_pool)
    integer, intent(out) :: size_random_vector_mBZ
    real(DP), intent(out) :: random_vector_cart_mBZ(:,:) !< (2, size_random_vector_pool)
    integer :: ii, ib_combination
    real(DP) :: b_combination_cart(3,4), b_frac_(3)
    integer :: b_aux_(3,4)
    PUSH_SUB(pick_points_within_mBZ_2D)

    !! Here we already use TRS to reduce number of b_aux_
    b_aux_(:,1) = (/ -1,-1, 0 /)
    b_aux_(:,2) = (/ -1, 0, 0 /)
    b_aux_(:,3) = (/ -1, 1, 0 /)
    b_aux_(:,4) =(/  0,-1, 0 /)

    do ib_combination = 1, 4
       b_frac_(1) = DBLE(b_aux_(1,ib_combination))/DBLE(qgrid(1))
       b_frac_(2) = DBLE(b_aux_(2,ib_combination))/DBLE(qgrid(2))
       b_frac_(3) = 0.0D0
       b_combination_cart(:,ib_combination) = 2.0D0 * PI_D / crys%alat * MATMUL(crys%bvec, b_frac_)
    enddo

    !! size_random_vector_mBZ <= size_random_vector_pool
    random_vector_cart_mBZ = 0.0D0
    size_random_vector_mBZ = 0
    iqran_loop : do ii = 1, size_random_vector_pool
       do ib_combination = 1, 4
          if (ABS(DOT_PRODUCT(random_vector_pool(:,ii),b_combination_cart(1:2,ib_combination))) > 0.5D0 * NORM2(b_combination_cart(1:2,ib_combination))**2 ) then
             cycle iqran_loop
          endif
       enddo
       size_random_vector_mBZ = size_random_vector_mBZ + 1
       random_vector_cart_mBZ(:,size_random_vector_mBZ) = random_vector_pool(:,ii)
    enddo iqran_loop
    ! write(*,'(A,F10.6)') "size_random_vector_mBZ / size_random_vector_pool = ", DBLE(size_random_vector_mBZ) / DBLE(size_random_vector_pool)

    POP_SUB(pick_points_within_mBZ_2D)
    return
  end subroutine pick_points_within_mBZ_2D

  subroutine pick_q2_inner_2D(crys, qgrid, q2_inner)
    type(crystal), intent(in) :: crys
    integer, intent(in) :: qgrid(3)
    real(DP), intent(out) :: q2_inner
    integer :: ib_combination
    real(DP) :: b_frac_(3), q2_list(4)
    integer :: b_aux_(3,4)
    PUSH_SUB(pick_q2_inner_2D)

    !! Here we already use TRS to reduce number of b_aux_
    b_aux_(:,1) = (/ -1,-1, 0 /)
    b_aux_(:,2) = (/ -1, 0, 0 /)
    b_aux_(:,3) = (/ -1, 1, 0 /)
    b_aux_(:,4) =(/  0,-1, 0 /)

    do ib_combination = 1, 4
       b_frac_(1) = DBLE(b_aux_(1,ib_combination))/DBLE(qgrid(1)) / 2.0D0
       b_frac_(2) = DBLE(b_aux_(2,ib_combination))/DBLE(qgrid(2)) / 2.0D0
       b_frac_(3) = 0.0D0
       q2_list(ib_combination) = DOT_PRODUCT(b_frac_, MATMUL(crys%bdot, b_frac_))
    enddo

    q2_inner = MINVAL(q2_list)

    POP_SUB(pick_q2_inner_2D)
    return
  end subroutine pick_q2_inner_2D

  pure function length2(x, U) result(r2)
    real(DP), intent(in)  :: x(3), U(3,3)
    real(DP) :: vmid(3)
    real(DP) :: r2
    vmid(1) = U(1,1)*x(1) + U(1,2)*x(2) + U(1,3)*x(3)
    vmid(2) =               U(2,2)*x(2) + U(2,3)*x(3)
    vmid(3) =                             U(3,3)*x(3)
    r2 = vmid(1)**2 + vmid(2)**2 + vmid(3)**2
  end function length2

end module vcoul_generator_m

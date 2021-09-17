#include "f_defs.h"

!=======================================================================
!
! Routines:
!
! (1) minibzaverage_3d_oneoverq2() Originally by JRD/MJ Last Modified: 8/27/2009 (MJ/JRD)
!
! (2) minibzaverage_3d_oneoverq() Originally by JRD/MJ  Last Modified: 8/27/2009 (MJ/JRD)
!
! (3) minibzaverage_2d_oneoverq2() Originally by JRD/MJ Last Modified: 9/15/2009 (MJ/JRD)
!
!  Output: average of <V_q> on the mini-BZ for a 1-D system.
!  output units: units equivalent to 8Pi/q^2
!
!=======================================================================

module minibzaverage_m
  use global_m
  ! use bessel_m
  use misc_m
  use slatec_m 
  implicit none
  private
  public :: minibzaverage_3d_oneoverq2_spherical_shell, minibzaverage_3d_oneoverq2_simple, minibzaverage_2d_oneoverq2_circular_shell, minibzaverage_2d_oneoverq2_simple, minibzaverage_3d_oneoverq2, minibzaverage_3d_oneoverq, minibzaverage_2d_oneoverq2,  minibzaverage_3d_oneoverq2_mod
  
contains

  subroutine minibzaverage_3d_oneoverq2(nn, bdot, integral, qran, qk, averagew, epshead, wcoul0, q0sph2, celvol, nfk)
    integer, intent(in) :: nn
    real(DP), intent(in) :: bdot(3,3), qran(:,:) !< (3, nn)
    real(DP), intent(out) :: integral
    real(DP), intent(in) ::  qk(3)
    logical, intent(in)  :: averagew
    SCALAR, intent(in) :: epshead
    real(DP), intent(in) ::  q0sph2
    real(DP), intent(in) :: celvol
    integer, intent(in) :: nfk
    real(DP) :: gkq(3), length,length_qk
    SCALAR, intent(inout)  :: wcoul0
    integer :: ii, nn2
    PUSH_SUB(minibzaverage_3d_oneoverq2)

    integral = 0D0
    !> In fact, we have already calculated ekinx = length_qk before calling this subroutine
    length_qk = DOT_PRODUCT(qk,MATMUL(bdot,qk))
    !> q0=0 + (G=0)
    if( length_qk < TOL_Zero ) then
       nn2 = nn ! nmc
       do ii = 1, nn2 ! Loop over Monte Carlo points
          gkq(:) = qk(:) + qran(:,ii) ! = qran(1:3,ii)
          length = DOT_PRODUCT(gkq,MATMUL(bdot,gkq))
          ! Skip the small value of q that will be integrated analytically later on
          if ( length > q0sph2 ) integral = integral + 1D0/length ! 1D0/length = 1/q^2
       enddo
       !> q1 + any G or q0 + G/=0
    else
       ! FHJ: for spherical integration regions, one can make the error per MC
       ! integration ~const. by choosing the number of points such that N ~ 1/ekinx.
       ! This is because, in 3D, error = sigma/N^{3/2}, and sigma ~ 1/ekinx^{3/2}
       ! If we fix the number of points such that N(ekinx=4*q0sph2) = nmc_coarse,
       ! Common/nrtype.f90:
       ! nmc_coarse = 250000
       ! nmc_fine = 2500000
       ! nmc = nmc_fine
       nn2 = nint(nmc_coarse * 4d0 * q0sph2 / length_qk)
       nn2 = max(1, min(nn2, nn))
       do ii = 1, nn2
          gkq(:) = qk(:) + qran(:,ii)
          length = DOT_PRODUCT(gkq,MATMUL(bdot,gkq))
          integral = integral + 1D0/length
       enddo
    endif
    integral = integral * 8D0 * PI_D / dble(nn2)
    !> q0 + G=0
    if( length_qk < TOL_Zero ) then
       !> Add analytical results within the sphere q0sph2
       !> nfk = number of qvectors in FBZ
       integral = integral + 32.0D0 * PI_D**2 * SQRT(q0sph2) / ( 8.0D0 * PI_D**3 / (celvol * dble(nfk)) )
    endif

    !> q0 + G = 0
    !> averagew = .true.
    !> epshead = epsinv(1,1) <== q0, G1=0, G2=0 component
    if (length_qk .lt. TOL_Zero .and. averagew) then
       wcoul0 = integral * epshead
    endif

    POP_SUB(minibzaverage_3d_oneoverq2)
    return
  end subroutine minibzaverage_3d_oneoverq2

  !> MC within the spherical shell + analytical integral
  !> See TMISC Sec. 15.14
  subroutine minibzaverage_3d_oneoverq2_spherical_shell(size_random_vector, crys, qgrid, integral, random_vector_cart, q2_inner, q2_outer)
    integer, intent(in) :: size_random_vector
    type (crystal), intent(in) :: crys
    integer, intent(in) :: qgrid(3)
    real(DP), intent(in) :: random_vector_cart(3,size_random_vector) !< (3, size_random_vector)
    real(DP), intent(out) :: integral
    real(DP), intent(in) ::  q2_inner, q2_outer
    integer :: nfk, ii, ib_combination, b_aux_(3,13)
    real(DP) :: prefactor1, prefactor2, b_combination_cart(3,13), b_frac_(3)
    logical :: outside_mBZ
    real(DP) :: integral_B, integral_C
    ! integer :: unit_qran_shell=234, counter_B, counter_C
    PUSH_SUB(minibzaverage_3d_oneoverq2_spherical_shell)

    integral = 0.0D0
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

    nfk = PRODUCT(qgrid)
    prefactor1 = 4.0D0 / 3.0D0 * PI_D * ( SQRT(q2_outer)**3 - SQRT(q2_inner)**3 ) / DBLE(size_random_vector)
    prefactor2 = DBLE(nfk) * crys%celvol / ( 16.0D0 * PI_D**3 )

    do ib_combination = 1, 13
       b_frac_(1) = DBLE(b_aux_(1,ib_combination))/DBLE(qgrid(1))
       b_frac_(2) = DBLE(b_aux_(2,ib_combination))/DBLE(qgrid(2))
       b_frac_(3) = DBLE(b_aux_(3,ib_combination))/DBLE(qgrid(3))
       b_combination_cart(:,ib_combination) = 2.0D0 * PI_D / crys%alat * MATMUL(crys%bvec, b_frac_)
    enddo

    ! if (peinf%inode .eq. 0) then
    !    do ib_combination = 1, 13
    !       write(*,'(A,I5,A,3F15.8)') "# ", ib_combination, " b = ", b_combination_cart(:,ib_combination)
    !    enddo
    ! endif

    !> MC part: (B-C)
    !> Loop over Monte Carlo points
    integral_B = 0.0D0
    integral_C = 0.0D0
    !> OMP will not speed up the code
    !> !$OMP PARALLEL DO default(shared) private(outside_mBZ) reduction(+:integral_B) reduction(-:integral_C)
    do ii = 1, size_random_vector
       outside_mBZ = .false.
       do ib_combination = 1, 13
          if (ABS(DOT_PRODUCT(random_vector_cart(:,ii),b_combination_cart(:,ib_combination))) > 0.5D0 * NORM2(b_combination_cart(:,ib_combination))**2) then
             outside_mBZ = .true.
             exit
          endif
       enddo
       if (.not. outside_mBZ) then
          ! integral = integral + 8.0D0 * PI_D/(NORM2(random_vector_cart(:,ii))**2)
          integral_B = integral_B + 8.0D0 * PI_D/(NORM2(random_vector_cart(:,ii))**2)
       else
          ! integral = integral - 8.0D0 * PI_D/(NORM2(random_vector_cart(:,ii))**2)
          integral_C = integral_C - 8.0D0 * PI_D/(NORM2(random_vector_cart(:,ii))**2)
       endif
    enddo
    !> !$OMP END PARALLEL DO
    integral = integral_B + integral_C
    !> Add analytical part: (A+D)
    integral = prefactor2 * ( prefactor1 * integral + 32.0D0 * PI_D**2 * ( SQRT(q2_inner) + SQRT(q2_outer) ) )

    POP_SUB(minibzaverage_3d_oneoverq2_spherical_shell)
    return
  end subroutine minibzaverage_3d_oneoverq2_spherical_shell

  ! !> [WORKING]
  ! subroutine minibzaverage_2d_oneoverq2_circular_shell(size_random_vector, crys, qgrid, integral, random_vector_cart, q2_inner, q2_outer)
  !   integer, intent(in) :: size_random_vector
  !   type (crystal), intent(in) :: crys
  !   integer, intent(in) :: qgrid(3)
  !   real(DP), intent(in) :: random_vector_cart(2,size_random_vector) !< (2, size_random_vector)
  !   real(DP), intent(out) :: integral
  !   real(DP), intent(in) ::  q2_inner, q2_outer
  !   integer :: nfk, ii, ib_combination, b_aux_(3,4)
  !   real(DP) :: qkxy(3), qkz(3), kxy, Lz, half_Lz
  !   real(DP) :: prefactor1, prefactor2, b_combination_cart(3,4), b_frac_(3)
  !   logical :: outside_mBZ
  !   real(DP) :: integral_B, integral_C
  !   complex(DPC) :: two_gamma
  !   ! integer :: unit_qran_shell=234, counter_B, counter_C
  !   PUSH_SUB(minibzaverage_2d_oneoverq2_circular_shell)

  !   integral = 0.0D0
  !   b_aux_(:,1) = (/ -1,-1, 0 /)
  !   b_aux_(:,2) = (/ -1, 0, 0 /)
  !   b_aux_(:,3) = (/ -1, 1, 0 /)
  !   b_aux_(:,4) =(/  0,-1, 0 /)
  !   nfk = PRODUCT(qgrid)
  !   Lz = SQRT(crys%adot(3,3))
  !   half_Lz = Lz/2.0D0

  !   !> [WORKING]
  !   prefactor1 = PI_D * (q2_outer - q2_inner) / DBLE(size_random_vector)
  !   !> prefactor2 = N_q Omega / (2 * (2 \pi)^2 Lz)
  !   prefactor2 = DBLE(nfk) * crys%celvol / ( 8.0D0 * PI_D**2 * Lz)

  !   ! ! !> Check area of mBZ [PASS]
  !   ! write(*,'(A,ES20.12,A,ES20.12,A,ES20.12)') "S_mBZ = ", 4.0D0*PI_D**2*Lz/DBLE(nfk)/crys%celvol, " |b1|^2 = ", crys%bdot(1,1), " S_mBZp = ", SQRT(3.0D0)/2.0D0*crys%bdot(1,1)/DBLE(nfk)

  !   do ib_combination = 1, 4
  !      b_frac_(1) = DBLE(b_aux_(1,ib_combination))/DBLE(qgrid(1))
  !      b_frac_(2) = DBLE(b_aux_(2,ib_combination))/DBLE(qgrid(2))
  !      b_frac_(3) = 0.0D0
  !      b_combination_cart(:, ib_combination) = 2.0D0 * PI_D / crys%alat * MATMUL(crys%bvec, b_frac_)
  !   enddo

  !   !> MC part: (B-C)
  !   !> Loop over Monte Carlo points
  !   integral_B = 0.0D0
  !   integral_C = 0.0D0
  !   !> OMP will not speed up the code
  !   !> !$OMP PARALLEL DO default(shared) private(kxy, outside_mBZ) reduction(+:integral_B) reduction(-:integral_C)
  !   do ii = 1, size_random_vector
  !      kxy = NORM2(random_vector_cart(1:2, ii))
  !      !> Here kz = 0
  !      outside_mBZ = .false.
  !      do ib_combination = 1, 4
  !         if (ABS(DOT_PRODUCT(random_vector_cart(:, ii),b_combination_cart(1:2, ib_combination))) > 0.5D0 * NORM2(b_combination_cart(1:2, ib_combination))**2) then
  !            outside_mBZ = .true.
  !            exit
  !         endif
  !      enddo
  !      if (.not. outside_mBZ) then
  !         integral_B = integral_B + 8.0D0 * PI_D/(NORM2(random_vector_cart(:, ii))**2) * (1.0D0-EXP(-kxy*half_Lz))
  !      else
  !         integral_C = integral_C - 8.0D0 * PI_D/(NORM2(random_vector_cart(:, ii))**2) * (1.0D0-EXP(-kxy*half_Lz))
  !      endif
  !   enddo
  !   !> !$OMP END PARALLEL DO
  !   integral = integral_B + integral_C
  !   two_gamma = cdig(ZERO, DCMPLX(Lz/2.0D0 * SQRT(q2_inner))) + cdig(ZERO, DCMPLX(Lz/2.0D0 * SQRT(q2_outer)))
  !   ! write(*,'(A,2ES20.12)') "two_gamma = ", two_gamma   
  !   !> [WORKING]
  !   !> Add analytical part: (A+D)
  !   integral = prefactor2 * ( prefactor1 * integral + 32.0D0 * PI_D**2 * (gamma_D + LOG(Lz/2.0D0)) + 16.0D0 * PI_D**2 * ( DBLE(two_gamma) + 0.5D0 * LOG(q2_inner * q2_outer) ) )

  !   POP_SUB(minibzaverage_2d_oneoverq2_circular_shell)
  !   return
  ! end subroutine minibzaverage_2d_oneoverq2_circular_shell

  !> [WORKING]
  subroutine minibzaverage_2d_oneoverq2_circular_shell(size_random_vector, crys, qgrid, integral, random_vector_cart, q2_inner, q2_outer)
    integer, intent(in) :: size_random_vector
    type (crystal), intent(in) :: crys
    integer, intent(in) :: qgrid(3)
    real(DP), intent(in) :: random_vector_cart(2,size_random_vector) !< (2, size_random_vector)
    real(DP), intent(out) :: integral
    real(DP), intent(in) ::  q2_inner, q2_outer
    integer :: nfk, ii, ib_combination, b_aux_(3,4)
    real(DP) :: kxy, Lz, half_Lz
    real(DP) :: prefactor1, prefactor2, b_combination_cart(3,4), b_frac_(3)
    logical :: outside_mBZ
    real(DP) :: integral_B, integral_C
    real(DP) :: two_gamma
    ! integer :: unit_qran_shell=234, counter_B, counter_C
    PUSH_SUB(minibzaverage_2d_oneoverq2_circular_shell)

    integral = 0.0D0
    b_aux_(:,1) = (/ -1,-1, 0 /)
    b_aux_(:,2) = (/ -1, 0, 0 /)
    b_aux_(:,3) = (/ -1, 1, 0 /)
    b_aux_(:,4) =(/  0,-1, 0 /)
    nfk = PRODUCT(qgrid)
    Lz = SQRT(crys%adot(3,3))
    half_Lz = Lz/2.0D0

    !> [WORKING]
    prefactor1 = PI_D * (q2_outer - q2_inner) / DBLE(size_random_vector)
    !> prefactor2 = N_q Omega / (2 * (2 \pi)^2 Lz)
    prefactor2 = DBLE(nfk) * crys%celvol / ( 8.0D0 * PI_D**2 * Lz)

    ! ! !> Check area of mBZ [PASS]
    ! write(*,'(A,ES20.12,A,ES20.12,A,ES20.12)') "S_mBZ = ", 4.0D0*PI_D**2*Lz/DBLE(nfk)/crys%celvol, " |b1|^2 = ", crys%bdot(1,1), " S_mBZp = ", SQRT(3.0D0)/2.0D0*crys%bdot(1,1)/DBLE(nfk)

    do ib_combination = 1, 4
       b_frac_(1) = DBLE(b_aux_(1,ib_combination))/DBLE(qgrid(1))
       b_frac_(2) = DBLE(b_aux_(2,ib_combination))/DBLE(qgrid(2))
       b_frac_(3) = 0.0D0
       b_combination_cart(:, ib_combination) = 2.0D0 * PI_D / crys%alat * MATMUL(crys%bvec, b_frac_)
    enddo

    !> MC part: (B-C)
    !> Loop over Monte Carlo points
    integral_B = 0.0D0
    integral_C = 0.0D0
    !> OMP will not speed up the code
    !> !$OMP PARALLEL DO default(shared) private(kxy, outside_mBZ) reduction(+:integral_B) reduction(-:integral_C)
    do ii = 1, size_random_vector
       kxy = NORM2(random_vector_cart(1:2, ii))
       !> Here kz = 0
       outside_mBZ = .false.
       do ib_combination = 1, 4
          if (ABS(DOT_PRODUCT(random_vector_cart(:, ii),b_combination_cart(1:2, ib_combination))) > 0.5D0 * NORM2(b_combination_cart(1:2, ib_combination))**2) then
             outside_mBZ = .true.
             exit
          endif
       enddo
       if (.not. outside_mBZ) then
          integral_B = integral_B + 8.0D0 * PI_D/(NORM2(random_vector_cart(:, ii))**2) * (1.0D0-EXP(-kxy*half_Lz))
       else
          integral_C = integral_C - 8.0D0 * PI_D/(NORM2(random_vector_cart(:, ii))**2) * (1.0D0-EXP(-kxy*half_Lz))
       endif
    enddo
    !> !$OMP END PARALLEL DO
    integral = integral_B + integral_C
    ! two_gamma = cdig(ZERO, DCMPLX(Lz/2.0D0 * SQRT(q2_inner))) + cdig(ZERO, DCMPLX(Lz/2.0D0 * SQRT(q2_outer)))
    !> [WORKING]
    !> DGAMIC is the complementary incomplete gamma function (https://people.math.sc.edu/Burkardt/f_src/slatec/slatec.html), also called the upper incomplete gamma function (https://en.wikipedia.org/wiki/Incomplete_gamma_function)
    two_gamma = DGAMIC(0.0D0, Lz/2.0D0 * SQRT(q2_inner)) + DGAMIC(0.0D0, Lz/2.0D0 * SQRT(q2_outer))
    
    ! write(*,'(A,2ES20.12)') "two_gamma = ", two_gamma   
    !> [WORKING]
    !> Add analytical part: (A+D)
    integral = prefactor2 * ( prefactor1 * integral + 32.0D0 * PI_D**2 * (gamma_D + LOG(Lz/2.0D0)) + 16.0D0 * PI_D**2 * ( two_gamma + 0.5D0 * LOG(q2_inner * q2_outer) ) )

    POP_SUB(minibzaverage_2d_oneoverq2_circular_shell)
    return
  end subroutine minibzaverage_2d_oneoverq2_circular_shell
  
  !> random_vector_cart(1:3,1:nrandom_vector_cart) are uniformly distributed points within mBZ
  subroutine minibzaverage_3d_oneoverq2_simple(crys, q_plus_G, size_random_vector, random_vector_cart, q2_min, integral)
    type (crystal), intent(in) :: crys
    real(DP), intent(in) :: q_plus_G(3) !> fractional coordinates
    integer, intent(in) :: size_random_vector
    real(DP), intent(in) :: random_vector_cart(:,:) !< (3, size_random_vector)
    real(DP), intent(out) :: integral
    real(DP), intent(in) :: q2_min
    real(DP) :: vector_cart(3), vector_random_cart(3), len2
    integer :: ii, size_random_vector_effective
    PUSH_SUB(minibzaverage_3d_oneoverq2_simple)

    len2 = DOT_PRODUCT(q_plus_G, MATMUL(crys%bdot, q_plus_G))
    ! size_random_vector_effective = MAX(nmc_low_threshold, NINT(size_random_vector * (q2_min**3) / (len2**3) ))
    !> error is a constant
    ! size_random_vector_effective = NINT(DBLE(size_random_vector) * (q2_min / len2)**3 + TOL_SMALL)
    !> error \propto 1/|q+G|^2
    
    if (len2 > TOL_ZERO) then
       if (len2 > q2_min) then
          size_random_vector_effective = NINT(DBLE(size_random_vector) * (q2_min / len2) + TOL_SMALL)
       else
          size_random_vector_effective = size_random_vector          
       endif
    else
       size_random_vector_effective = size_random_vector       
    endif

    ! write(*,'(A,3F12.6,A,F15.8)') "q_plus_G = ", q_plus_G(:), " size_random_vector_effective / size_random_vector = ", DBLE(size_random_vector_effective)/DBLE(size_random_vector)
    ! !> error \propto 1/|q+G|^3
    ! size_random_vector_effective = size_random_vector
    !> error \propto 1/|q+G|^2.5
    ! size_random_vector_effective = NINT(DBLE(size_random_vector) * SQRT(q2_min / len2) + TOL_SMALL)
    if (size_random_vector_effective .le. 1) then
       if (len2 < TOL_ZERO) then
          call die("len2 = 0", only_root_writes=.true.)
       endif       
       ! integral = 8.0D0 * PI_D/(NORM2(vector_cart)**2)
       integral = 8.0D0 * PI_D / len2
       return
    endif

    vector_cart = 2.0D0 * PI_D / crys%alat * MATMUL(crys%bvec, q_plus_G)
    ! write(*,'(A,ES20.12,A,ES20.12)') "len2 = ", len2, " NORM2(vector_cart)**2 = ", NORM2(vector_cart)**2
    integral = 0.0D0
    !> OMP here will not speed up the calculation
    !$OMP PARALLEL private(vector_random_cart) IF(size_random_vector_effective > 10000)
    !$OMP DO reduction(+:integral)
    do ii = 1, size_random_vector_effective
       vector_random_cart(:) = vector_cart(:) + random_vector_cart(:,ii)
       integral = integral + 1.0D0/(NORM2(vector_random_cart)**2)
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
    integral = integral * 8.0D0 * PI_D / DBLE(size_random_vector_effective)

    POP_SUB(minibzaverage_3d_oneoverq2_simple)
    return
  end subroutine minibzaverage_3d_oneoverq2_simple

  !========================================================================
  !> This is for Slab Truncation
  subroutine minibzaverage_2d_oneoverq2(nn, bdot, integral, qran, qk, kz, zc, epshead, q0len, averagew, wcoul0)
    integer, intent(in) :: nn
    real(DP), intent(in) :: bdot(3,3), qran(:,:) !< (3,nn)
    real(DP), intent(out) :: integral
    real(DP), intent(in) ::  qk(3)
    logical, intent(in)  :: averagew
    SCALAR, intent(in) :: epshead
    SCALAR, intent(inout)  :: wcoul0
    real(DP), intent(in) ::  zc, q0len
    real(DP), intent(out) :: kz
    integer ::  ii
    real(DP) :: gkq(3), length, kxy, gkqxy(3),lengthqk
    real(DP) :: gkqz(3),epsmodel,gamma,alpha,vc,vc_qtozero
    SCALAR  :: integralW
    PUSH_SUB(minibzaverage_2d_oneoverq2)

    ! Sahar:
    ! Define Gamma parameter for model epsilon (see Sohrab, PRB 2006)
    ! Extract the quadratic dependence of 1/epsinv(00)
    ! 1/epsinv(q;0,0) = 1 + q^2*vc(q)*gamma

    ! q0len = sqrt(DOT_PRODUCT(q0vec,MATMUL(bdot,q0vec))) in units of [Bohr]^{-1}
    ! for q0, q0len is finite
    ! zc=2D0*PI_D/(sqrt(bdot(3,3))*2D0) = a3/2 in units of [Bohr]
    !get Vc (q -> 0, G=0)
    vc_qtozero=((1.0d0 - exp(-q0len*zc))/q0len**2)
    ! Define Gamma
    gamma = (1.0d0/epshead-1.0d0)/((q0len**2)*vc_qtozero)
    ! Define alpha
    ! Set to zero for now
    alpha = 0.0d0
    ! length of q + G
    lengthqk = sqrt(DOT_PRODUCT(qk, MATMUL(bdot, qk)))
    integral = 0D0
    integralW = 0D0

    do ii = 1, nn ! nmc Monte Carlo points
       gkq(:) = qk(:)
       gkq(1:2) = gkq(1:2) + qran(1:2,ii)
       gkqxy(1:2) = gkq(1:2)
       gkqxy(3) = 0D0
       kxy=sqrt(DOT_PRODUCT(gkqxy, MATMUL(bdot, gkqxy))) ! in units of [Bohr]^{-1}
       length = DOT_PRODUCT(gkq, MATMUL(bdot, gkq))

       ! This is Temporary??
       gkqz(:)=gkq(:)
       gkqz(1)=0D0
       gkqz(2)=0D0
       kz=sqrt(DOT_PRODUCT(gkqz, MATMUL(bdot, gkqz)))

       ! First average v
       ! This is the most general expression of slab truncated Coulomb interaction, we haven't assume R=Lz/2 here
       integral = integral + (1.0d0+exp(-kxy*zc) * ((kz/kxy)*sin(kz*zc) - cos(kz*zc))) / length

       ! Do we also want to average W?
       ! This is a waste of time if we are not qk=0
       ! ======
       ! <- q0 + G=0 ->
       if (lengthqk.lt.TOL_zero.and.averagew) then
          ! Use model epsilon here
          ! Normalize integral by head of epsilon
          vc = ((1.0d0 - exp(-kxy*zc))/kxy**2)
          epsmodel=1.0d0 + vc * kxy**2 * gamma*exp(-alpha*kxy)
          integralW = integralW + (vc/epsmodel)
          !      write(6,*)  'USING MODEL EPSILON FOR AVERAGING OF W'
          !      write(6,*)  'gamma: ', gamma, 'alpha: ', alpha, 'qk', qk
          !      write(6,*)   'qk', qk
          ! No model epsilon here
       endif
    enddo
    ! Convert integral to Ry
    integral = integral * 8D0 * PI_D / dble(nn)
    if (lengthqk.lt.TOL_zero.and.averagew) then
       wcoul0 = integralW * 8D0 * PI_D / dble(nn)
    endif

    POP_SUB(minibzaverage_2d_oneoverq2)
    return
  end subroutine minibzaverage_2d_oneoverq2

  !> only for averagew = T and |q_plus_G|<TOL_zero
  subroutine minibzaverage_2d_oneoverq2_simple(crys, q_plus_G, size_random_vector, random_vector_cart, ekinx_min, integral)
    type (crystal), intent(in) :: crys
    real(DP), intent(in) ::  q_plus_G(3)
    integer, intent(in) :: size_random_vector
    real(DP), intent(in) :: random_vector_cart(:,:) !< (2, size_random_vector)
    real(DP), intent(in) :: ekinx_min
    real(DP), intent(out) :: integral
    real(DP) :: q_plus_G_xy(3), q_plus_G_z(3), kxy, kz, half_Lz, len2 !, kxy_, kz_
    real(DP) :: vector_cart(3), vector_random_cart(3)
    integer :: ii, size_random_vector_effective
    PUSH_SUB(minibzaverage_2d_oneoverq2_simple)

    !> (G+q) in xy plane
    q_plus_G_xy(1:2) = q_plus_G(1:2)
    q_plus_G_xy(3) = 0.0D0
    !> kxy = |q+G|_{\parallel}
    kxy = SQRT(DOT_PRODUCT(q_plus_G_xy, MATMUL(crys%bdot, q_plus_G_xy)))
    q_plus_G_z(1:2) = 0.0D0
    q_plus_G_z(3) = q_plus_G(3)
    !> kz = |q+G|_z = |G|_z, since we require that qvec(3) = 0
    kz = SQRT(DOT_PRODUCT(q_plus_G_z, MATMUL(crys%bdot, q_plus_G_z)))
    half_Lz = SQRT(crys%adot(3,3))/2.0D0

    !> |q+G|^2
    len2 = DOT_PRODUCT(q_plus_G, MATMUL(crys%bdot, q_plus_G))
    !> error \propto 1/|q+G|^2
    if (len2 > TOL_ZERO) then
       if (len2 > ekinx_min) then
          size_random_vector_effective = NINT(DBLE(size_random_vector) * (ekinx_min / len2) + TOL_SMALL)
       else
          size_random_vector_effective = size_random_vector          
       endif
    else
       size_random_vector_effective = size_random_vector       
    endif
    
    ! write(*,'(I5,3ES20.12)') size_random_vector, ekinx_min, len2, DBLE(size_random_vector) * (ekinx_min / len2) + TOL_SMALL
    
    if (size_random_vector_effective .le. 1) then
       if (len2 < TOL_ZERO) then
          call die("len2 = 0", only_root_writes=.true.)
       endif
       integral = 8.0D0 * PI_D / len2 * (1.0D0-EXP(-kxy*half_Lz)*COS(kz*half_Lz))
       return
    endif

    vector_cart = 2.0D0 * PI_D / crys%alat * MATMUL(crys%bvec, q_plus_G)
    ! !> [DEBUGGING]
    ! kxy_ = NORM2(vector_cart(1:2))
    ! kz_ = ABS(vector_cart(3))
    ! write(*,'(4(A,ES20.12))') "kxy_ = ", kxy_, " kz_ = ", kz_, " kxy = ", kxy, " kz = ", kz

    integral = 0.0D0
    ! write(*,'(A,I10,A,I10,A,ES20.12,A,ES20.12)') "size_effective = ", size_random_vector_effective, " size = ", size_random_vector, "ekinx_min = ", ekinx_min, "len2 = ", len2

    !$OMP PARALLEL private(vector_random_cart) IF(size_random_vector_effective > 10000)
    !$OMP DO reduction(+:integral)
    do ii = 1, size_random_vector_effective
       vector_random_cart(1:2) = vector_cart(1:2) + random_vector_cart(1:2, ii)
       vector_random_cart(3) = vector_cart(3)
       kxy = NORM2(vector_random_cart(1:2))
       kz = ABS(vector_random_cart(3))
       integral = integral + 1.0D0/(NORM2(vector_random_cart)**2) * (1.0D0-EXP(-kxy*half_Lz)*COS(kz*half_Lz))
    enddo
    !$OMP END DO
    !$OMP END PARALLEL    
    integral = integral * 8.0D0 * PI_D / DBLE(size_random_vector_effective)

    POP_SUB(minibzaverage_2d_oneoverq2_simple)
    return
  end subroutine minibzaverage_2d_oneoverq2_simple

  !========================================================================
  ! call from vcoul_generator.f90:
  ! call minibzaverage_3d_oneoverq(nn,bdot,dvalue,qran,qvec_mod)
  ! nn = nmc = Number of Monte Carlo points
  ! [Output] integral = dvalue
  ! qran(1:3,1:nmc) = random small qvector around Gamma
  ! qvec_mod = sig%qvec (for q0 point, qvec_mod = 0.0, for q1 point, it is just q1(1:3))
  subroutine minibzaverage_3d_oneoverq(nn, bdot, integral, qran, qk)
    integer, intent(in) :: nn
    real(DP), intent(in) :: bdot(3,3), qran(:,:) !< (3,nn)
    real(DP), intent(out) :: integral
    real(DP), intent(in) ::  qk(3)
    integer :: ii
    real(DP) :: gkq(3), length
    PUSH_SUB(minibzaverage_3d_oneoverq)
    
    integral = 0D0
    !> Loop over Monte Carlo points
    do ii = 1, nn ! nmc
       !> qk(1:3) = qvec_mod(1:3)
       gkq(:) = qk(:) + qran(:,ii) ! gkq(1:3) is around qk(1:3)=rq(1:3,iq_qRBZ_k_outer)
       length = DOT_PRODUCT(gkq,MATMUL(bdot,gkq))
       length = sqrt(length)
       !> integral is just the integral of 1/q within a qpoint cube (or grid), we assume epsinv is a constant within this cube and only integrate 1/q
       integral = integral + 1D0/length
    enddo
    !> [IMPORTANT]
    !> integral = 8 * Pi * 1/Vq * \int_{ around q_i} dq 1/q = 8 * Pi * <1/q> within q_i cube
    integral = integral * 8D0 * PI_D / dble(nn)

    POP_SUB(minibzaverage_3d_oneoverq)
    return
  end subroutine minibzaverage_3d_oneoverq

  !===========================================================================

  ! For modified coulomb interaction
  subroutine minibzaverage_3d_oneoverq2_mod(nn,bdot,integral,qran,qk,coulomb_mod)
    integer, intent(in) :: nn
    real(DP), intent(in) :: bdot(3,3), qran(:,:) !< (3, nn)
    real(DP), intent(out) :: integral
    real(DP), intent(in) ::  qk(3)
    type(coulomb_modifier_t), intent(in)  :: coulomb_mod
    real(DP) :: gkq(3), screeninv, temp_exp, length
    integer :: ii
    PUSH_SUB(minibzaverage_3d_oneoverq2_mod)

    integral = 0D0
    screeninv = 1.0D0/(4.0D0 * coulomb_mod%screening_length *coulomb_mod%screening_length)
    ! convert screening_length from A^{-1} to Bohr
    screeninv = screeninv/(BOHR*BOHR)
    do ii = 1, nn
       gkq(:) = qk(:) + qran(:,ii)
       length = DOT_PRODUCT(gkq,MATMUL(bdot,gkq))
       temp_exp = exp(-screeninv*length)
       integral = integral + 1D0/length*(temp_exp*coulomb_mod%long_range_frac_fock + &
            (1.0D0 - temp_exp)*coulomb_mod%short_range_frac_fock)
    enddo
    integral = integral * 8D0 * PI_D / dble(nn)

    POP_SUB(minibzaverage_3d_oneoverq2_mod)
    return
  end subroutine minibzaverage_3d_oneoverq2_mod
  
end module minibzaverage_m

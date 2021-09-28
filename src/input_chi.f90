#include "f_defs.h"
!> Read parameters from file WFN, Initialize k-points sampling, kg type, Initialize G-space, gvec
!> output: crys,gvec,syms,kg types

subroutine input_chi(pol, crys, gvec, kp, kg, syms, syms_wfn, intwfnv, intwfnc, kpq, kgq, symsq_wfn, intwfnvq)
  use global_m
  use fullbz_m
  use input_utils_m
  use misc_m
  use sort_m
  use wfn_rho_vxc_io_m
  ! use read_rho_vxc_m
  use scalapack_m
  use hdf5
  use eqpcor_m
  implicit none

  type (polarizability), intent(inout) :: pol
  type (crystal), intent(out) :: crys
  type (gspace), intent(out) :: gvec
  type (grid), intent(out) :: kg, kgq
  type (kpoints), intent(out) :: kp, kpq
  type (symmetry), intent(out) :: syms, syms_wfn, symsq_wfn
  !> intwfnc and intwfnv stores WFN, intwfnvq stores WFNmq or WFNq
  type (int_wavefunction), intent(out) :: intwfnc, intwfnv, intwfnvq

  type (symmetry) :: symsq
  type(gspace) :: gvec_kpt, gvecq
  type (crystal) :: crysq
  integer :: irk, irkq, irk_loc, irk_loc_, irkq_loc
  integer :: is, ig, ikq, ik
  real(DP) :: delta,qq_temp(3)
  SCALAR, allocatable :: cg(:,:), cg_c_(:,:,:), cg_v_(:,:,:)
  character(len=3) :: sheader
  integer :: iflavor
  integer :: ib, itran, umk(3), gumk(3), icb_relative, ivb_relative
  real(DP) :: dq(3), q_temp(3), q0_vector(3)
  integer, allocatable :: isort_(:)
  integer :: nvb_diff
  integer :: request_v, request_c, irk_loc_target, irkq_loc_target
  integer :: irk_start, irk_end, nrk_loc, nrk_loc_max, max1_nrk_loc, rk_blocksize
  integer :: ipes, ipes_, irk_start_, irk_end_, nrk_loc_, ib_min, ib_max
  character(LEN=50) :: file_eqp

  call logit('input_chi:  reading WFN')
  if (peinf%inode .eq. 0) call open_file(25, file='WFN', form='unformatted', status='old')
  sheader = 'WFN'
  iflavor = 0
  call read_binary_header_type(25, sheader, iflavor, kp, gvec, syms, crys)
  !> prepare syms%mtrx_reci and syms%mtrx_cart
  call prepare_syms(crys, syms)
  !> Common/input_utils.f90:
  !> Write a warning if any k-point is nonzero in a truncated direction
  ! call check_trunc_kpts(pol%icutv, kp)
  
  pol%nband = pol%nvb + pol%ncb

  ib_min = MINVAL(kp%ifmax(:,:)) - pol%nvb + 1
  ib_max = MAXVAL(kp%ifmax(:,:)) + pol%ncb
  if (ib_min <= 0) then
     call die("Invalid number of valence bands.", only_root_writes=.true.)
  endif
  if (ib_max > kp%mnband) then
     if (peinf%inode .eq. 0) then
        write(*,'(A,I5,A,I5,A,I5,A,I5)') "min_ifmax", MINVAL(kp%ifmax(:,:)), " max_ifmax", MAXVAL(kp%ifmax(:,:)), "nvb = ", pol%nvb, " ncb = ", pol%ncb        
        write(*,'(A,I5,A,I5)') "ib_max = ", ib_max, " kp%mnband = ", kp%mnband
        call die("Not enough bands in WFN.", only_root_writes=.true.)
     endif
  endif

  SAFE_ALLOCATE(kp%elda, (kp%mnband, kp%nrk, kp%nspin))
  kp%elda(:,:,:) = kp%el(:,:,:)

  if (pol%eqp_corrections) then
     file_eqp = 'eqp.dat'
     if (peinf%inode .eq. 0) then
        write(6,'(1X,A,I5,A,I5)') "Use "//TRUNC(file_eqp)//" for bands from ", ib_min, " to ", ib_max
     endif
     call eqpcor(crys, kp, ib_min, ib_max, file_eqp)
  endif
  call calc_efermi(kp, pol%efermi, ib_min, ib_max)

  !> eqp correction
  if (MINVAL(kp%ifmax(:,:)) .ne. MAXVAL(kp%ifmax(:,:))) then
     call die("Not a semiconductor", only_root_writes=.true.)
  endif

  call logit('input_kernel:  reading gvec info')
  SAFE_ALLOCATE(gvec%components, (3, gvec%ng))
  call read_binary_gvectors(25, gvec%ng, gvec%ng, gvec%components)

  if (ANY(kp%ifmax(:,:) .eq. 0)) then
     call die("BSE codes cannot handle a system where some k-points have no occupied bands.", only_root_writes = .true.)
  endif

  !> Number of valence bands
  nvb_diff = MAXVAL(kp%ifmax(:,:)) - MINVAL(kp%ifmax(:,:))
  if (nvb_diff > 0) then
     call die("Metal not supported.", only_root_writes=.true.)
  else
     kp%nvband = MAXVAL(kp%ifmax(:,:))
  endif
  !> Number of available conduction bands
  kp%ncband = kp%mnband - kp%nvband

  !> inread_kernel.f90: pol%nvb = number_val_bands
  if (pol%nvb .gt. kp%nvband) then
     call die("The requested number of valence bands is not available in WFN.", only_root_writes=.true.)
  endif
  !> inread_kernel.f90: pol%ncb = number_cond_bands
  if (pol%ncb .gt. kp%ncband) then
     call die("The requested number of conduction bands is not available in WFN.", only_root_writes=.true.)
  endif

  !> Common/input_utils.f90:
  !> magnitude index <=> FFT index
  call gvec_index(gvec)

  !> Generate full brillouin zone from irreducible wedge, rk -> fk
  if (peinf%inode .eq. 0) then
     write(*,'(1X,A)') "Use symmetries to expand the coarse grid sampling."
     write(*,'(1X,A)') "Use symmetries to expand the shifted coarse grid sampling."
  endif
  call timacc(7,1)

  kg%nr = kp%nrk
  SAFE_ALLOCATE(kg%r, (3, kg%nr))
  kg%r(:, 1:kg%nr) = kp%rk(:, 1:kp%nrk)
  call subgroup(syms_wfn, syms, crys, kg%r(:,1))
  call fullbz(crys, syms_wfn, kg)

  call timacc(7,2)

  if (peinf%verb_high .and. peinf%inode .eq. 0) then
     write(6,'(/1x,a6,14x,a7,12x,2(1x,a6),3x,a3)') 'i', 'k-point', 'indr', 'itran', 'kg0'
     write(6,'(1x,6("-"),1x,32("-"),2(1x,6("-")),1x,8("-"))')
     do ik = 1, kg%nf
        write(6,'(1x,i6,3(1x,f10.6),2(1x,i6),3(1x,i2))') ik, kg%f(:,ik), kg%indr(ik), kg%itran(ik), kg%kg0(:,ik)
     enddo
  endif

  !> If there is a finite center-of-mass momentum, Q, find mapping between k and k+Q
  SAFE_ALLOCATE(pol%indexq, (kg%nf))
  do ik = 1, kg%nf
     pol%indexq(ik) = ik
  enddo
  pol%qgrid = kp%kgrid

  if (pol%nq0 > 0) then
     !> [WORKING]
     !> We only allow one q0 in calc_chi.f90 calculation
     q0_vector(:) = pol%qpt(:,1)
     !> [WORKING]
     !> Broken TRS
     if (pol%k_plus_q) then
        if (peinf%inode .eq. 0) call open_file(26, file='WFNq', form='unformatted', status='old')
     else
        if (peinf%inode .eq. 0) call open_file(26, file='WFNmq', form='unformatted', status='old')
     endif

     sheader = 'WFN'
     iflavor = 0
     call read_binary_header_type(26, sheader, iflavor, kpq, gvecq, symsq, crysq)
     call prepare_syms(crysq, symsq)
     call check_header('WFN', kp, gvec, syms, crys, 'WFNmq', kpq, gvecq, symsq, crysq, is_wfn = .true.)
     if (ANY(kp%kgrid(1:3) .ne. kpq%kgrid(1:3))) then
        if(peinf%inode .eq. 0) then
           write(0,'(A,3I5)') 'WFN  kgrid = ', kp%kgrid(1:3)
           write(0,'(A,3I5)') 'WFNmq kgrid = ', kpq%kgrid(1:3)
        endif
        call die('kgrids for WFN and WFNmq must be the same', only_root_writes = .true.)
     endif

     SAFE_ALLOCATE(kpq%elda, (kpq%mnband, kpq%nrk, kpq%nspin))
     kpq%elda(:,:,:) = kpq%el(:,:,:)

     !> For WFNmq case, we need to call eqpcor_2 again for kpq
     !> kpq is for valence states
     if (pol%eqp_corrections) then
        file_eqp = 'eqpmq.dat'
        !> Here we only correct part of valence states
        ib_min = MINVAL(kpq%ifmax(:,:)) - pol%nvb + 1
        ib_max = MAXVAL(kpq%ifmax(:,:))
        call eqpcor(crys, kpq, ib_min, ib_max, file_eqp)
     endif

     SAFE_ALLOCATE(gvecq%components, (3, gvecq%ng))
     call read_binary_gvectors(26, gvecq%ng, gvecq%ng, gvecq%components)

     if (ANY(kpq%ifmax(:,:) .eq. 0)) then
        call die("BSE codes cannot handle a system where some k-points have no occupied bands.", only_root_writes = .true.)
     endif

     nvb_diff = MAXVAL(kpq%ifmax(:,:)) - MINVAL(kpq%ifmax(:,:))
     if (nvb_diff > 0) then
        call die("Metal not supported II.", only_root_writes=.true.)
     else
        kpq%nvband = MAXVAL(kpq%ifmax(:,:))
     endif
     kpq%ncband = kpq%mnband - kpq%nvband

     if (kpq%nvband .ne. kp%nvband) then
        call die("Number of valence bands mismatch between WFN and WFNmq.", only_root_writes=.true.)
     endif

     kgq%nr = kpq%nrk
     SAFE_ALLOCATE(kgq%r, (3,kgq%nr))
     kgq%r(1:3,1:kgq%nr) = kpq%rk(1:3,1:kpq%nrk)
     !> Reset symsq here, original symsq = syms = crystal symmetry
     if (peinf%inode .eq. 0) then
        write(*,*) "construct little group of q"
     endif

     !> Check that the q0 and kgq%r are commensurate
     !> That is, the small q used in WFNmq (for valence states) should be -q0
     !> [WORKING] broken TRS
     if (pol%k_plus_q) then
        if ( any(abs(kgq%r(:,1) - kg%r(:,1) - q0_vector(:)) > TOL_Small )  ) then
           if (peinf%inode .eq. 0) then
              write(6,'(A, 3F15.7)') "kgq%r(:,1) - kg%r(:,1) - q0(:) = ", kgq%r(:,1) - kg%r(:,1) - q0_vector(:)
           endif
           call die("q0 and WFN grids are not compatible", only_root_writes = .true.)
        endif
     else
        if ( any(abs(kgq%r(:,1) - kg%r(:,1) + q0_vector(:)) > TOL_Small )  ) then
           if (peinf%inode .eq. 0) then
              write(6,'(A, 3F15.7)') "kgq%r(:,1) - kg%r(:,1) + q0(:) = ", kgq%r(:,1) - kg%r(:,1) + q0_vector(:)
           endif
           call die("q0 and WFN grids are not compatible", only_root_writes = .true.)
        endif
     endif

     call subgroup(symsq, syms, crys, kgq%r(:,1)-kg%r(:,1))

     ! !> Determine the relative shift using the first RBZ kpoint in WFN and WFNmq
     ! q_temp(:) = kgq%r(:,1) - kg%r(:,1)
     ! ! if (NORM2(q0_vector - q_temp) > TOL_ZERO) then
     ! !    call die("q0_vector .ne. q_temp.", only_root_writes=.true.)
     ! ! endif
     ! symsq%ntran = 0
     ! symsq%mtrx(:,:,:) = 0
     ! symsq%mtrx_reci(:,:,:) = 0.D0
     ! symsq%mtrx_cart(:,:,:) = 0.D0
     ! symsq%kgzero(:,:) = 0.D0
     ! symsq%tnp(:,:) = 0.D0
     ! do itran = 1, syms%ntran
     !    dq(1:3) = MATMUL(DBLE(syms%mtrx_reci(1:3, 1:3, itran)), q_temp(1:3)) - q_temp(1:3)
     !    call get_gumk3(crys%bdot, dq, gumk)
     !    dq(:) = dq(:) - DBLE(gumk(:))
     !    if (all(abs(dq(1:3)) .lt. TOL_Small)) then
     !       !> Store index of element of subgroup
     !       symsq%ntran = symsq%ntran + 1
     !       symsq%mtrx(1:3,1:3,symsq%ntran) = syms%mtrx(1:3,1:3,itran)
     !       symsq%mtrx_reci(1:3,1:3,symsq%ntran) = syms%mtrx_reci(1:3,1:3,itran)
     !       symsq%mtrx_cart(1:3,1:3,symsq%ntran) = syms%mtrx_cart(1:3,1:3,itran)
     !       symsq%tnp(1:3,symsq%ntran) = syms%tnp(1:3,itran)
     !       !> matmul(symsq%mtrx_reci(1:3, 1:3, itranq), q_temp(1:3)) + symsq%kgzero(1:3,itranq) = q_temp(1:3), itranq = 1, symsq%ntran
     !       symsq%kgzero(1:3, symsq%ntran) = - gumk(1:3)
     !    endif
     ! enddo

     call subgroup(symsq_wfn, syms, crys, kgq%r(:,1))
     call fullbz(crysq, symsq_wfn, kgq)

     !> [WORKING]
     !> This only works for WFNmq, need to modify for WFNq
     !> Find mapping between kgq and kg: indexq
     do ik = 1, kg%nf
        !> For each kc, find the corresponding kv, such that, kv = kc - Q
        do ikq = 1, kgq%nf
           !> [WORKING] broken TRS
           if (pol%k_plus_q) then
              !> kv = kc + Q - umk, kv in FBZ, kc in FBZ, Q in FBZ, kc + Q could be outside of FBZ
              !> kv <-- ikq
              !> kv <-- ik
              qq_temp(:) = kg%f(:,ik) + q0_vector(:) - kgq%f(:,ikq)
           else
              !> kv = kc - Q - umk, kv in FBZ, kc in FBZ, Q in FBZ, kc - Q could be outside of FBZ
              !> kv <-- ikq
              !> kv <-- ik
              qq_temp(:) = kg%f(:,ik) - q0_vector(:) - kgq%f(:,ikq)
           endif
           call get_gumk3(crys%bdot, qq_temp, umk)
           qq_temp(:) = qq_temp(:) - dble(umk(:))
           delta = DOT_PRODUCT(qq_temp,MATMUL(crys%bdot,qq_temp))
           if (delta < TOL_ZERO) then
              exit
           endif
        enddo
        if (delta .gt. TOL_Zero) then
           if(peinf%inode.eq.0) then
              write(0,*) '  Could not find point equivalent to ', kg%f(:,ik), " ik = ", ik
           endif
           call die('q0 not commensurate with kgrid of WFNmq',only_root_writes = .true.)
        endif
        !> kg%f(:,ik) - q0(:) = kgq%f(:,indexq(ik)) + umk(:,ik)
        pol%indexq(ik) = ikq
        kgq%f(:,ikq) = kgq%f(:,ikq) + DBLE(umk(:))
        !> When generating |v (fk-q) >, we will allow fk-q outside FBZ
        kgq%kg0(:,ikq) = kgq%kg0(:,ikq) + umk(:)
     enddo ! ik
     SAFE_DEALLOCATE_P(gvecq%components)
  else
     !> Copy kp into kpq, copy kg into kgq, copy syms into symsq
     !> WFNmq = WFN
     !> Initialize symsq with syms
     symsq%ntran            = syms%ntran
     symsq%mtrx(:,:,:)      = syms%mtrx(:,:,:)
     symsq%mtrx_reci(:,:,:) = syms%mtrx_reci(:,:,:)
     symsq%mtrx_cart(:,:,:) = syms%mtrx_cart(:,:,:)
     symsq%kgzero(:,:)      = syms%kgzero(:,:)
     symsq%tnp(:,:)         = syms%tnp(:,:)

     symsq_wfn%ntran            = syms_wfn%ntran
     symsq_wfn%mtrx(:,:,:)      = syms_wfn%mtrx(:,:,:)
     symsq_wfn%mtrx_reci(:,:,:) = syms_wfn%mtrx_reci(:,:,:)
     symsq_wfn%mtrx_cart(:,:,:) = syms_wfn%mtrx_cart(:,:,:)
     symsq_wfn%kgzero(:,:)      = syms_wfn%kgzero(:,:)
     symsq_wfn%tnp(:,:)         = syms_wfn%tnp(:,:)

     !> Initialize kpq with kp
     kpq%nspinor = kp%nspinor
     kpq%nspin   = kp%nspin
     kpq%nrk     = kp%nrk
     kpq%mnband  = kp%mnband
     kpq%nvband  = kp%nvband
     kpq%ncband  = kp%ncband
     kpq%kgrid   = kp%kgrid
     kpq%shift   = kp%shift
     kpq%ecutwfc = kp%ecutwfc
     kpq%ngkmax  = kp%ngkmax

     SAFE_ALLOCATE(kpq%ngk, (kpq%nrk))
     SAFE_ALLOCATE(kpq%ifmin, (kpq%nrk,kpq%nspin))
     SAFE_ALLOCATE(kpq%ifmax, (kpq%nrk,kpq%nspin))
     SAFE_ALLOCATE(kpq%w, (kpq%nrk))
     SAFE_ALLOCATE(kpq%rk, (3,kpq%nrk))
     SAFE_ALLOCATE(kpq%el, (kpq%mnband,kpq%nrk,kpq%nspin))
     SAFE_ALLOCATE(kpq%occ, (kpq%mnband,kpq%nrk,kpq%nspin))

     kpq%ngk(:)     = kp%ngk(:)
     kpq%ifmin(:,:) = kp%ifmin(:,:)
     kpq%ifmax(:,:) = kp%ifmax(:,:)
     kpq%w(:)       = kp%w(:)
     kpq%rk(:,:)    = kp%rk(:,:)
     kpq%el(:,:,:)  = kp%el(:,:,:)
     kpq%occ(:,:,:) = kp%occ(:,:,:)

     !> initialize kgq with kpq
     kgq%nr = kg%nr
     kgq%nf = kg%nf
     SAFE_ALLOCATE(kgq%r, (3,kgq%nr))
     SAFE_ALLOCATE(kgq%f, (3,kgq%nf))
     SAFE_ALLOCATE(kgq%indr, (kgq%nf))
     SAFE_ALLOCATE(kgq%itran, (kgq%nf))
     SAFE_ALLOCATE(kgq%kg0, (3,kgq%nf))
     kgq%r(:,:) = kg%r(:,:)
     kgq%f(:,:) = kg%f(:,:)
     kgq%indr(:) = kg%indr(:)
     kgq%itran(:) = kg%itran(:)
     kgq%kg0(:,:) = kg%kg0(:,:)
  endif

  !> Distribute conduction bands of WFN and valence bands of WFN[m]q
  call distrib_chi(pol, kg, kgq)

  !> Distribute valence bands of WFN
  rk_blocksize = iceil(kp%nrk, peinf%npes)
  nrk_loc      = NUMROC(kp%nrk, rk_blocksize, peinf%inode, 0, peinf%npes)
  nrk_loc_max  = NUMROC(kp%nrk, rk_blocksize,           0, 0, peinf%npes)
  max1_nrk_loc = MAX(1, nrk_loc)
  if (nrk_loc >= 1) then
     irk_start = INDXL2G(      1, rk_blocksize, peinf%inode, 0, peinf%npes)
     irk_end   = INDXL2G(nrk_loc, rk_blocksize, peinf%inode, 0, peinf%npes)
  else
     irk_start = 0
     irk_end = 0
  endif

  ! write(*,'(A,I5,A,I5,A,I5,A,I5,A,I5,A,I5)') "rk_blocksize = ", rk_blocksize, " nrk_loc = ", nrk_loc, " nrk_loc_max = ", nrk_loc_max, " max1_nrk_loc = ", max1_nrk_loc, " irk_start = ", irk_start, " irk_end = ", irk_end

  if (pol%nq0 > 0) then
     SAFE_ALLOCATE(intwfnvq%cgk, (kpq%ngkmax, pol%nvb, kpq%nspin*kpq%nspinor, peinf%nrkq(peinf%inode+1)))
     SAFE_ALLOCATE(intwfnvq%isort, (gvec%ng, peinf%nrkq(peinf%inode+1)))
     SAFE_ALLOCATE(intwfnvq%ng, (peinf%nrkq(peinf%inode+1)))
     intwfnvq%cgk = ZERO
     intwfnvq%isort = 0
     intwfnvq%nspin = kpq%nspin
     intwfnvq%nspinor = kpq%nspinor
  endif

  if (pol%nq1 > 0) then
     !> intwfnv%cgk are distributed among all procs
     !> intwfnv%ng can be replaced by kp%ngk
     SAFE_ALLOCATE(intwfnv%ng, (max1_nrk_loc))
     SAFE_ALLOCATE(intwfnv%isort, (gvec%ng, max1_nrk_loc))
     SAFE_ALLOCATE(intwfnv%cgk,  (kp%ngkmax, pol%nvb, kp%nspin*kp%nspinor,  max1_nrk_loc))
     intwfnv%ng = 0
     intwfnv%isort = 0
     intwfnv%cgk = ZERO
     intwfnv%nspin = kp%nspin
     intwfnv%nspinor = kp%nspinor
  endif

  SAFE_ALLOCATE(intwfnc%cgk,  (kp%ngkmax, pol%ncb,   kp%nspin*kp%nspinor,  peinf%nrk(peinf%inode+1)))
  SAFE_ALLOCATE(intwfnc%isort, (gvec%ng, peinf%nrk(peinf%inode+1)))
  SAFE_ALLOCATE(intwfnc%ng, (peinf%nrk(peinf%inode+1)))
  intwfnc%cgk = ZERO
  intwfnc%isort = 0
  intwfnc%nspin = kp%nspin
  intwfnc%nspinor = kp%nspinor

  SAFE_ALLOCATE(isort_, (gvec%ng))
  SAFE_ALLOCATE(gvec_kpt%components, (3, kp%ngkmax))
  SAFE_ALLOCATE(cg, (kp%ngkmax, kp%nspin*kp%nspinor))
  SAFE_ALLOCATE(cg_c_, (kp%ngkmax, pol%ncb, kp%nspin*kp%nspinor))
  SAFE_ALLOCATE(cg_v_, (kp%ngkmax, pol%nvb, kp%nspin*kp%nspinor))

  do irk = 1, kp%nrk
     call read_binary_gvectors(25, kp%ngk(irk), kp%ngkmax, gvec_kpt%components)
     do ig = 1, kp%ngk(irk)
        call findvector(isort_(ig), gvec_kpt%components(:, ig), gvec)
        if (isort_(ig) .eq. 0) call die('could not find gvec', only_root_writes=.true.)
     enddo

     do ib = 1, kp%mnband
        !> ROOT reads wavefunction, bcase=.false. means no broadcast to all proc
        call read_binary_data(25, kp%ngk(irk), kp%ngkmax, kp%nspin*kp%nspinor, cg, bcast=.false.)
        if (peinf%inode .eq. 0) then
           call checknorm('WFN', ib, irk, kp%nspin, cg(1:kp%ngk(irk),:))
           do is = 1, kp%nspin
              !> Conduction bands
              if ( (ib .ge. (kp%ifmax(irk,is) + 1)) .and. (ib .le. (kp%ifmax(irk,is) + pol%ncb)) ) then
                 icb_relative = ib - kp%ifmax(irk,is)
                 !> Here only 1:kp%ngk(irk) are meaningful, but we copy the entire array for simplicity
                 cg_c_(:, icb_relative, is:is*kp%nspinor) = cg(:, is:is*kp%nspinor)
              endif
              !> Valence bands
              if (pol%nq1 > 0) then
                 if ( (ib .le. kpq%ifmax(irk,is)) .and. (ib .ge. (kpq%ifmax(irk,is) - pol%nvb + 1)) ) then
                    ivb_relative = kpq%ifmax(irk,is) - ib + 1
                    !> Here only 1:kp%ngk(irk) are meaningful, but we copy the entire array for simplicity
                    cg_v_(:, ivb_relative, is:is*kpq%nspinor) = cg(:, is:is*kpq%nspinor)
                 endif
              endif
           enddo
        endif
     enddo !> ib
     ! write(*,'(A,I5,A,2ES20.12)') "irk = ", irk, " cg_v_(1,1,1) = ", cg_v_(1,1,1)

     !> Old way
     ! call MPI_BCAST(cg_c_, kp%ngkmax*kp%nspin*kp%nspinor*pol%ncb, MPI_SCALAR, 0, MPI_COMM_WORLD, mpierr)
     ! irk_loc = peinf%irk_g2l(irk, peinf%inode + 1)
     ! if (irk_loc .ne. 0) then
     !    intwfnc%ng(irk_loc) = kp%ngk(irk)
     !    intwfnc%isort(1:kp%ngk(irk), irk_loc) = isort_(1:kp%ngk(irk))
     !    intwfnc%cgk(1:kp%ngk(irk),:,:,irk_loc) = cg_c_(1:kp%ngk(irk),:,:)
     ! endif

     irk_loc = peinf%irk_g2l(irk, peinf%inode+1)
     !> Current proc needs this wavefunction
     if (irk_loc .ne. 0) then
        !> ROOT needs this wavefunction
        if (peinf%inode .ne. 0) then
           call MPI_IRECV(cg_c_(1,1,1), kp%ngkmax*pol%ncb*kp%nspin*kp%nspinor, MPI_SCALAR, MPI_ANY_SOURCE, irk, MPI_COMM_WORLD, request_c, mpierr)
        endif
     endif
     !> ROOT sends out this wavefunction to target procs, there could be many target procs!
     if (peinf%inode .eq. 0) then
        do ipes = 2, peinf%npes
           irk_loc_target = peinf%irk_g2l(irk, ipes)
           if (irk_loc_target .ne. 0) then
              call MPI_SEND(cg_c_(1,1,1), kp%ngkmax*pol%ncb*kp%nspin*kp%nspinor, MPI_SCALAR, ipes-1, irk, MPI_COMM_WORLD, mpierr)
           endif
        enddo
     endif
     !> non-ROOT finishes receiving the wavefunction
     if (irk_loc .ne. 0) then
        if (peinf%inode .ne. 0) then
           call MPI_WAIT(request_c, MPI_STATUS_IGNORE, mpierr)
        endif
        intwfnc%ng(irk_loc) = kp%ngk(irk)
        intwfnc%isort(1:kp%ngk(irk), irk_loc) = isort_(1:kp%ngk(irk))
        intwfnc%cgk(1:kp%ngk(irk), :, :, irk_loc) = cg_c_(1:kp%ngk(irk), :, :)
     endif

     !> Distribute valence bands of WFN
     !> Old way
     ! if (pol%nq1 > 0) then
     !    irkq = irk
     !    call MPI_BCAST(cg_v_, kpq%ngkmax*kpq%nspin*kpq%nspinor*pol%nvb, MPI_SCALAR, 0, MPI_COMM_WORLD, mpierr)
     !    irkq_loc = peinf%irkq_g2l(irkq, peinf%inode+1)
     !    if (irkq_loc .ne. 0) then
     !       intwfnvq%ng(irkq_loc) = kpq%ngk(irkq)
     !       intwfnvq%isort(1:kpq%ngk(irkq), irkq_loc) = isort_(1:kpq%ngk(irkq))
     !       intwfnvq%cgk(1:kpq%ngk(irkq),:,:,irkq_loc) = cg_v_(1:kpq%ngk(irkq),:,:)
     !    endif
     ! endif

     if (pol%nq1 > 0) then
        irk_loc_ = INDXG2L(irk, rk_blocksize, 0, 0, peinf%npes)
        ipes_ = INDXG2P(irk, rk_blocksize, 0, 0, peinf%npes) + 1 !> = 1, ..., peinf%npes
        !> Current proc needs this wavefunction
        if (ipes_ .eq. peinf%inode+1) then
           if (peinf%inode .ne. 0) then
              call MPI_IRECV(cg_v_(1,1,1), kp%ngkmax*pol%nvb*kp%nspin*kp%nspinor, MPI_SCALAR, MPI_ANY_SOURCE, irk, MPI_COMM_WORLD, request_v, mpierr)
           endif
        endif
        !> ROOT MPI_SEND to target proc
        if (peinf%inode .eq. 0) then
           if (ipes_-1 .ne. 0) then
              call MPI_SEND(cg_v_(1,1,1), kp%ngkmax*pol%nvb*kp%nspin*kp%nspinor, MPI_SCALAR, ipes_-1, irk, MPI_COMM_WORLD, mpierr)
           endif
        endif
        if (ipes_ .eq. peinf%inode+1) then
           if (peinf%inode .ne. 0) then
              call MPI_WAIT(request_v, MPI_STATUS_IGNORE, mpierr)
           endif
           intwfnv%ng(irk_loc_) = kp%ngk(irk)
           !> [important]
           intwfnv%isort(1:kp%ngk(irk), irk_loc_) = isort_(1:kp%ngk(irk))
           intwfnv%cgk(1:kp%ngk(irk), :, :, irk_loc_) = cg_v_(1:kp%ngk(irk), :, :)
           ! write(*,'(A,I5,A,I5,A,2ES20.12)') " irk = ", irk, " irk_loc_ = ", irk_loc_, " intwfnv%cgk(1,1,1,irk_loc_) = ", intwfnv%cgk(1, 1, 1, irk_loc_)
        endif
     endif !> nq1 > 0
  enddo !> irk

  !> Check norm of intwfnc%cgk
  do irk_loc = 1, peinf%nrk(peinf%inode+1)
     irk = peinf%irk_l2g(irk_loc)
     do ib = 1, pol%ncb
        call checknorm('intwfnc', ib, irk, kp%nspin, intwfnc%cgk(1:kp%ngk(irk), ib, :, irk_loc))
     enddo
  enddo

  if (pol%nq1 > 0) then
     !> Check norm of intwfnv%cgk
     do irk_loc_ = 1, nrk_loc
        irk = INDXL2G( irk_loc_, rk_blocksize, peinf%inode, 0, peinf%npes)
        do ib = 1, pol%nvb
           call checknorm('intwfnv', ib, irk, kp%nspin, intwfnv%cgk(1:kp%ngk(irk), ib, :, irk_loc_))
        enddo
     enddo
  endif

  SAFE_DEALLOCATE(cg)
  SAFE_DEALLOCATE(cg_c_)
  SAFE_DEALLOCATE_P(gvec_kpt%components)
  SAFE_DEALLOCATE(cg_v_)

  !> Write out info about xtal
  if (peinf%inode .eq. 0) then
     write(6,'(/1x,a)') 'Crystal wavefunctions read from file WFN:'
     write(6,'(1x,a,i0)') '- Number of k-points in WFN: ', kp%nrk
     write(6,'(1x,a,i0)') '- Number of k-points in the full BZ of WFN: ', kg%nf
     if (peinf%verb_high) then
        write(6,'(1x,a)') '- K-points:'
        write(6,'(1(2x,3(1x,f10.6)))') kg%r(1:3, 1:kg%nr)
     endif
     call close_file(25)
  endif

  if (pol%nq0 > 0) then
     SAFE_ALLOCATE(cg_v_, (kpq%ngkmax, pol%nvb, kpq%nspin*kpq%nspinor))
     SAFE_ALLOCATE(gvec_kpt%components, (3, kpq%ngkmax))
     SAFE_ALLOCATE(cg, (kpq%ngkmax, kpq%nspin*kpq%nspinor))

     do irkq = 1, kpq%nrk
        call read_binary_gvectors(26, kpq%ngk(irkq), kpq%ngkmax, gvec_kpt%components)
        do ig = 1, kpq%ngk(irkq)
           call findvector(isort_(ig), gvec_kpt%components(:, ig), gvec)
           if (isort_(ig) .eq. 0) call die('could not find gvec', only_root_writes=.true.)
        enddo

        do ib = 1, kpq%mnband
           !> ROOT read wavefunction, bcase=.false. means no broadcast to all proc
           call read_binary_data(26, kpq%ngk(irkq), kpq%ngkmax, kpq%nspin*kpq%nspinor, cg, bcast=.false.)
           if (peinf%inode .eq. 0) then
              call checknorm('WFNmq', ib, irkq, kpq%nspin, cg(1:kpq%ngk(irkq),:))
              do is = 1, kpq%nspin
                 if ( (ib .le. kpq%ifmax(irkq,is)) .and. (ib .ge. (kpq%ifmax(irkq,is) - pol%nvb + 1)) ) then
                    ivb_relative = kpq%ifmax(irkq,is) - ib + 1
                    !> Here only 1:kpq%ngk(irkq) are meaningful, but we copy the entire array for simplicity
                    cg_v_(:, ivb_relative, is:is*kpq%nspinor) = cg(:, is:is*kpq%nspinor)
                 endif
              enddo
           endif
        enddo !> ib

        ! !> Old way
        ! call MPI_BCAST(cg_v_, kpq%ngkmax*kpq%nspin*kpq%nspinor*pol%nvb, MPI_SCALAR, 0, MPI_COMM_WORLD, mpierr)
        ! irkq_loc = peinf%irkq_g2l(irkq, peinf%inode+1)
        ! if (irkq_loc .ne. 0) then
        !    intwfnvq%ng(irkq_loc) = kpq%ngk(irkq)
        !    intwfnvq%isort(1:kpq%ngk(irkq),irkq_loc) = isort_(1:kpq%ngk(irkq))
        !    intwfnvq%cgk(1:kpq%ngk(irkq),:,:,irkq_loc) = cg_v_(1:kpq%ngk(irkq),:,:)
        ! endif

        irkq_loc = peinf%irkq_g2l(irkq, peinf%inode+1)
        !> Current proc needs this wavefunction
        if (irkq_loc .ne. 0) then
           if (peinf%inode .ne. 0) then
              call MPI_IRECV(cg_v_(1,1,1), kpq%ngkmax*pol%nvb*kpq%nspin*kpq%nspinor, MPI_SCALAR, MPI_ANY_SOURCE, irkq, MPI_COMM_WORLD, request_v, mpierr)
           endif
        endif
        !> ROOT sends out this wavefunction to target procs
        if (peinf%inode .eq. 0) then
           do ipes = 2, peinf%npes
              irkq_loc_target = peinf%irkq_g2l(irkq, ipes)
              if (irkq_loc_target .ne. 0) then
                 call MPI_SEND(cg_v_(1,1,1), kpq%ngkmax*pol%nvb*kpq%nspin*kpq%nspinor, MPI_SCALAR, ipes-1, irkq, MPI_COMM_WORLD, mpierr)
              endif
           enddo
        endif
        !> non-ROOT finishes receiving the wavefunction
        if (irkq_loc .ne. 0) then
           if (peinf%inode .ne. 0) then
              call MPI_WAIT(request_v, MPI_STATUS_IGNORE, mpierr)
           endif
           intwfnvq%ng(irkq_loc) = kpq%ngk(irkq)
           intwfnvq%isort(1:kpq%ngk(irkq), irkq_loc) = isort_(1:kpq%ngk(irkq))
           intwfnvq%cgk(1:kpq%ngk(irkq), :, :, irkq_loc) = cg_v_(1:kpq%ngk(irkq), :, :)
           ! write(*,'(A,I5,A,I5,A,2ES20.12)') " irkq = ", irkq, " irkq_loc = ", irkq_loc, " intwfnvq%cgk(1,1,1,irkq_loc) = ", intwfnvq%cgk(1, 1, 1, irkq_loc)
        endif
     enddo !> irkq

     !> Check norm of intwfnvq%cgk
     do irkq_loc = 1, peinf%nrkq(peinf%inode+1)
        irkq = peinf%irkq_l2g(irkq_loc)
        do ib = 1, pol%nvb
           call checknorm('intwfnvq', ib, irkq, kpq%nspin, intwfnvq%cgk(1:kpq%ngk(irkq), ib, :, irkq_loc))
        enddo
     enddo

     if (peinf%inode.eq.0) then
        write(6,'(/1x,a)') 'Crystal wavefunctions read from file WFNmq:'
        write(6,'(1x,a,i0)') '- Number of k-points in WFNmq: ', kgq%nr
        write(6,'(1x,a,i0)') '- Number of k-points in the full BZ of WFNmq: ', kgq%nf
        if (peinf%verb_high) then
           write(6,'(1x,a)') '- K-points:'
           write(6,'(1(2x,3(1x,f10.6)))') kgq%r(1:3,1:kgq%nr)
        endif
        call close_file(26)
     endif

     SAFE_DEALLOCATE_P(gvec_kpt%components)
     SAFE_DEALLOCATE(cg)
     SAFE_DEALLOCATE(cg_v_)
  endif ! nq0 > 0

  SAFE_DEALLOCATE(isort_)

  return
end subroutine input_chi

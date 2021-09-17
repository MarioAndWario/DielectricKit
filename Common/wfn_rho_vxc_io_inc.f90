!=========================================================================
!
! Included from file wfn_rho_vxc_io.f90.
! You are not expected to understand this. --DAS
!
!=========================================================================

#ifdef READ
#define READ_WRITE(x) read ## x
#define INTENT out
#define FLAVOR_INTENT inout
#else
#define READ_WRITE(x) write ## x
#define INTENT in
#define FLAVOR_INTENT in
#endif

#ifdef FORMATTED
#define NAME(x) READ_WRITE(_format ## x)
#define FORMAT , *
#elif defined BINARY
#define NAME(x) READ_WRITE(_binary ## x)
#define FORMAT
#else
#define NAME(x) READ_WRITE(x)
#endif

#define NAMEHEADER NAME(_header)

#undef GVECTORS
#ifdef TEMP_COMPLEX
#define LONGNAME(x) NAME(x ## _complex_data)
#define TEMP_SCALAR complex(DPC)
#define MPI_TEMP_SCALAR MPI_COMPLEX_DPC
#elif defined TEMP_REAL
#define LONGNAME(x) NAME(x ## _real_data)
#define TEMP_SCALAR real(DP)
#define MPI_TEMP_SCALAR MPI_REAL_DP
#else
#define GVECTORS
#define LONGNAME(x) NAME(x ## _gvectors)
#endif
#ifdef GVECTORS
#define ARGUMENTS gvec
#else
#define ARGUMENTS ns, data
#endif

#if defined FORMATTED || defined BINARY
#ifndef TEMP_SCALAR
!=========================================================================
subroutine NAMEHEADER(iunit, sheader, iflavor, ns, ng, ntran, cell_symmetry, &
     nat, nk, nbands, ngkmax, ecutrho, ecutwfc, FFTgrid, kgrid, kshift, celvol, &
     alat, avec, adot, recvol, blat, bvec, bdot, mtrx, tnp, atyp, apos, ngk, &
     kw, kpt, ifmin, ifmax, energies, occupations, nspinor, warn, dont_warn_kgrid, sdate, stime, version)
  integer, intent(in) :: iunit !< unit number
  character(len=3), intent(inout) :: sheader !< file header 'WFN'/'RHO'/'VXC'/'GET' -- last one is to read, and return it
  integer, intent(FLAVOR_INTENT) :: iflavor !< define type. always must be initialized. modified only if -1 on input
  !! -1 = read from file and return it, 0 = as defined by -DCPLX, 1 = real, 2 = complex
  integer, intent(INTENT) :: ns, ng, ntran, cell_symmetry, nat, nk, nbands, ngkmax !< numbers of spins, G-vectors, symmetries,
  !! cell type (0 = cubic, 1 = hexagonal), numbers of atoms, k-points, bands, max(ngk)
  real(DP), intent(INTENT) :: ecutrho, ecutwfc !< charge-density and wave-function cutoffs, in Ry
  integer, intent(INTENT) :: FFTgrid(3), kgrid(3) !< FFT grid size, k-grid size
  real(DP), intent(INTENT) :: kshift(3) !< k-grid offset
  real(DP), intent(INTENT) :: celvol, alat, avec(3, 3), adot(3, 3) !< cell volume, lattice constant,
  !! lattice vectors, metric tensor in real space (in a.u., avec in units of alat)
  real(DP), intent(INTENT) :: recvol, blat, bvec(3, 3), bdot(3, 3) !< cell volume, lattice constant,
  !! lattice vectors, metric tensor in reciprocal space (in a.u., bvec in units of blat)
  integer, intent(INTENT) :: mtrx(3, 3, 48) !< symmetry matrix
  real(DP), intent(INTENT) :: tnp(3, 48) !< fractional translation
  integer, pointer, intent(INTENT) :: atyp(:) !< atomic species
  real(DP), pointer, intent(INTENT) :: apos(:,:) !< atomic positions (in units of alat)
  integer, pointer, intent(INTENT) :: ngk(:) !< number of G-vectors for each k-pt, ngk(nk)
  real(DP), pointer, intent(INTENT) :: kw(:), kpt(:, :) !< k-weight, kw(nk); k-coord, kpt(3, nk) in crystal coords
  integer, pointer, intent(INTENT) :: ifmin(:, :), ifmax(:, :) !< lowest and highest occupied band, ifmin/max(nk, ns)
  real(DP), pointer, intent(INTENT) :: energies(:, :, :) !< energies(nbands, nk, ns) in Ry
  real(DP), pointer, intent(INTENT) :: occupations(:, :, :) !< occupations(nbands, nk, ns) between 0 and 1
  integer, optional, intent(INTENT) :: nspinor !< 2 if doing a two-component spinor calculation; 1 if not present
  logical, optional, intent(in) :: warn !< if false, suppresses warnings about inversion symmetries etc.
  logical, optional, intent(in) :: dont_warn_kgrid !< if true, validity of kgrid read will not be checked.
  !! use for inteqp to allow interpolation onto non-uniform fine grids
  character(len=32), optional, intent(out) :: sdate, stime !< if read, result from file is returned; if write, current is returned
  integer, optional, intent(INTENT) :: version !< version of the header. -1 if no version was/should be specified

  character :: sdate_*32, stime_*32, stitle*32, sflavor*7
  integer :: ii, jj, is, ib, ik, itran, iat, ierr, nspin
  logical :: wfnflag, readwrite_version

#ifdef READ
  logical :: is_get, warn_, warn_kgrid, die_spinors
#else
  character :: adate*11, atime*14
#endif

  PUSH_SUB(NAMEHEADER)

  if (sheader .eq. 'WFN') then
     wfnflag = .true.
  elseif (sheader .eq. 'RHO' .or. sheader .eq. 'VXC' .or. sheader .eq. 'KER') then
     wfnflag = .false.
  elseif (sheader .ne. 'GET') then
     call die("Unknown file header '" + sheader + "' passed to " + TOSTRING(NAMEHEADER) &
          + ": must be 'WFN'/'RHO'/'VXC'/'KER'/'GET'")
  endif

#ifndef READ
  if(sheader .eq. 'GET') call die("Header 'GET' may not be passed to " + TOSTRING(NAMEHEADER) &
       + ": must be 'WFN'/'RHO'/'VXC'")
#endif

  if (peinf%inode .eq. 0) then

     select case(iflavor)
     case(-1)
#ifndef READ
        call die("iflavor = -1 cannot be passed to " + TOSTRING(NAMEHEADER))
#endif
     case(0)
        sflavor = MYFLAVOR
     case(1)
        sflavor = RFLAVOR
     case(2)
        sflavor = CFLAVOR
     case default
        write(sflavor, '(i7)') iflavor
        call die("Illegal value iflavor = " + TRUNC(sflavor) + " passed to " + TOSTRING(NAMEHEADER) &
             + ": must be -1,0,1,2.")
     end select

#ifndef READ
     call date_time(adate, atime)
     sdate_ = adate
     stime_ = atime
     stitle = sheader + "-" + sflavor
#endif

     ! FHJ: we try to READ the version tag if the argument "version" was passed,
     ! and we only WRITE the version if the argument is present and /= -1.
     readwrite_version=.false.
     if (present(version)) then
#ifdef READ
        readwrite_version = .true.
#else
        readwrite_version = version/=-1
#endif
     endif

     ierr = -1
     ! FHJ: try to read/write a first time, including the version tag
     if (readwrite_version) then
        READ_WRITE()(iunit FORMAT, iostat = ierr) version, stitle, sdate_, stime_
     endif

     ! FHJ: if there was no version tag, or if we don`t care about version...
     if (ierr /= 0) then
        ! FHJ: go back one record if we tried reading before
        if (readwrite_version) backspace(iunit, iostat = ierr)
        READ_WRITE()(iunit FORMAT, iostat = ierr) stitle, sdate_, stime_

        ! ! ======
        ! write(*,*) "stitle = ", stitle, "sdate = ", sdate_, "stime_ = ", stime_
        ! ! ------
        
        if(ierr /= 0) then
           call die("Failed " + TOSTRING(READ_WRITE()) + " operation in " + TOSTRING(NAMEHEADER) &
                + " on header record in mode " + sheader + ".")
        endif
#ifdef READ
        if (present(version)) version = -1
#endif
     endif

#ifdef READ
     is_get = .false.
     if(sheader == 'GET') then
        sheader = stitle(1:3)
        is_get = .true.
        wfnflag = (sheader == 'WFN')
     endif

     if(iflavor == -1) then
        sflavor = TRUNC(stitle(5:))
        if(TRUNC(sflavor) .eq. RFLAVOR) then
           iflavor = 1
        else if (TRUNC(sflavor) .eq. CFLAVOR) then
           iflavor = 2
        else
           call die("Read unknown flavor '" + TRUNC(sflavor) + "' in " + &
                TOSTRING(NAMEHEADER) + ": must be '" + RFLAVOR + "'/'" + CFLAVOR + "'")
        endif
     endif
     if(TRUNC(stitle) .ne. TRUNC(sheader + "-" + sflavor)) then
        call die("File header mismatch: read '" + TRUNC(stitle) + "' but expected '" &
             + sheader + "-" + TRUNC(sflavor) + "'")
     endif
#endif

#ifndef READ
     nspin = ns
     if(present(nspinor)) then
        if(nspinor==2) then
           nspin = 4
        endif
     endif
#endif

     ! scalar variables
     if(wfnflag) then
        READ_WRITE()(iunit FORMAT) nspin, ng, ntran, cell_symmetry, nat, ecutrho, nk, nbands, ngkmax, ecutwfc
#if defined READ
        if(nk < 1) call die("nk out of bounds in wfn")
        if(nbands < 1) call die("nbands out of bounds in wfn")
        if(ngkmax < 1) call die("ngkmax out of bounds in wfn")
        if(ecutwfc < -TOL_Zero) call die("ecutwfc out of bounds in wfn")
#endif
     else
        READ_WRITE()(iunit FORMAT) nspin, ng, ntran, cell_symmetry, nat, ecutrho
     endif

#if defined READ
     if(nspin /= 1 .and. nspin /= 2 .and. nspin /= 4) call die("nspin out of bounds")
     if(nspin == 4) then
        if(.not. present(nspinor)) call die("nspinor not passed but read nspin==4")
        nspinor = 2
        ns = 1

!        die_spinors = .true.
        die_spinors = .false.

        if(die_spinors) call die("Cannot use spinor WFN file (nspin = 4).")
     else
        if(present(nspinor)) nspinor = 1
        ns = nspin
     endif
     if(ng < 1) call die("ng out of bounds")
     if(ntran < 1 .or. ntran > 48) call die("ntran out of bounds")
     if(cell_symmetry /= 0 .and. cell_symmetry /= 1) call die("cell_symmetry out of bounds")
     if(nat < 1) call die("nat out of bounds")
     if(ecutrho < -TOL_Zero) call die("ecutrho out of bounds")
#endif

     ! arrays of fixed size
     if(wfnflag) then
        READ_WRITE()(iunit FORMAT) (FFTgrid(ii), ii = 1, 3), (kgrid(ii), ii = 1, 3), (kshift(ii), ii = 1, 3)
#if defined READ
        warn_kgrid = .true.
        if(present(dont_warn_kgrid)) warn_kgrid = .not. dont_warn_kgrid
        if(warn_kgrid) then
           if (any(kgrid(1:3) < 1)) call die("kgrid out of bounds in wfn")
           if (all(abs(kshift(1:3)) < TOL_Zero) .and. product(kgrid(1:3)) < nk) then
              call die("kgrid too small for nk")
           endif
           ! You might think such a condition would hold always but it does not necessarily, since a shifted uniform grid
           ! generally does not include all the points you can get from unfolding it with symmetries. --DAS
           if (product(kgrid(1:3)) > nk * ntran) call die("kgrid too large compared to unfolded nk")
        endif
        if(any(abs(kshift(1:3)) > 1 + TOL_Zero)) call die("kshift out of bounds in wfn")
#endif
     else
        READ_WRITE()(iunit FORMAT) (FFTgrid(ii), ii = 1, 3)
     endif
     READ_WRITE()(iunit FORMAT) celvol, alat, ((avec(jj, ii), jj = 1, 3), ii = 1, 3), ((adot(jj, ii), jj = 1, 3), ii = 1, 3)
     READ_WRITE()(iunit FORMAT) recvol, blat, ((bvec(jj, ii), jj = 1, 3), ii = 1, 3), ((bdot(jj, ii), jj = 1, 3), ii = 1, 3)
#ifdef READ
     mtrx(:,:,:) = 0
     tnp(:,:) = 0d0
#endif
     READ_WRITE()(iunit FORMAT) (((mtrx(ii, jj, itran), ii = 1, 3), jj = 1, 3), itran = 1, ntran)
     READ_WRITE()(iunit FORMAT) ((tnp(jj, itran), jj = 1, 3), itran = 1, ntran)
#if defined READ
     call make_identity_symmetry_first(ntran, mtrx, tnp)
     if(any(FFTgrid(1:3) < 1)) then
        call die("FFTgrid out of bounds")
     endif
     if(product(FFTgrid(1:3)) * 1.05 < ng * 6 / PI_D) then ! consistency of FFT grid with G-vectors, with a 5% tolerance
        write(0,*) 'FFTgrid = ', FFTgrid
        write(0,*) 'ng = ', ng
        if(product(FFTgrid(1:3)) < ng) then
           call die("FFTgrid inconsistent with ng")
        else
           write(0,*) 'WARNING: FFTgrid is suspiciously small. This is ok only if this system is'
           write(0,*) 'extremely anisotropic. Otherwise, use the Visual/gsphere.py utility to double'
           write(0,*) 'check that the G-vectors are compatible with this FFT grid.'
        endif
     endif
     if(celvol < -TOL_Zero) call die("celvol out of bounds")
     if(recvol < -TOL_Zero) call die("recvol out of bounds")
     if(any(abs(mtrx(:,:,1:ntran)) > 1)) call die("symmetry matrix may only contain -1, 0, 1")
#endif
  endif

  if(present(sdate)) sdate = sdate_
  if(present(stime)) stime = stime_
#ifdef MPI
  if(peinf%npes > 1) then
#ifdef READ
     if(present(version)) call MPI_BCAST(version, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
     call MPI_BCAST(iflavor, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
#endif
     if(present(sdate)) call MPI_BCAST(sdate, 32, MPI_CHARACTER, 0, MPI_COMM_WORLD, mpierr)
     if(present(stime)) call MPI_BCAST(stime, 32, MPI_CHARACTER, 0, MPI_COMM_WORLD, mpierr)
  endif
#endif

#ifdef READ
#ifdef MPI
  if (peinf%npes .gt. 1) then
     call MPI_BCAST(is_get, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, mpierr)
     if(is_get) call MPI_BCAST(sheader, 3, MPI_CHARACTER, 0, MPI_COMM_WORLD, mpierr)
     call MPI_BCAST(ns, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
     if(present(nspinor)) call MPI_BCAST(nspinor, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
     call MPI_BCAST(ng, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
     call MPI_BCAST(ntran, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
     call MPI_BCAST(cell_symmetry, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
     call MPI_BCAST(nat, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
     call MPI_BCAST(ecutrho, 1, MPI_REAL_DP, 0, MPI_COMM_WORLD, mpierr)
     call MPI_BCAST(FFTgrid, 3, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
     call MPI_BCAST(celvol, 1, MPI_REAL_DP, 0, MPI_COMM_WORLD, mpierr)
     call MPI_BCAST(alat, 1, MPI_REAL_DP, 0, MPI_COMM_WORLD, mpierr)
     call MPI_BCAST(avec, 9, MPI_REAL_DP, 0, MPI_COMM_WORLD, mpierr)
     call MPI_BCAST(adot, 9, MPI_REAL_DP, 0, MPI_COMM_WORLD, mpierr)
     call MPI_BCAST(recvol, 1, MPI_REAL_DP, 0, MPI_COMM_WORLD, mpierr)
     call MPI_BCAST(blat, 1, MPI_REAL_DP, 0, MPI_COMM_WORLD, mpierr)
     call MPI_BCAST(bvec, 9, MPI_REAL_DP, 0, MPI_COMM_WORLD, mpierr)
     call MPI_BCAST(bdot, 9, MPI_REAL_DP, 0, MPI_COMM_WORLD, mpierr)
     call MPI_BCAST(mtrx(1,1,1), 3*3*48, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
     call MPI_BCAST(tnp, 3*48, MPI_REAL_DP, 0, MPI_COMM_WORLD, mpierr)
     call MPI_BCAST(wfnflag, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, mpierr)
     if (wfnflag) then
        call MPI_BCAST(nk, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        call MPI_BCAST(nbands, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        call MPI_BCAST(ngkmax, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        call MPI_BCAST(ecutwfc, 1, MPI_REAL_DP, 0, MPI_COMM_WORLD, mpierr)
        call MPI_BCAST(kgrid, 3, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        call MPI_BCAST(kshift, 3, MPI_REAL_DP, 0, MPI_COMM_WORLD, mpierr)
     endif
  endif
#endif
  ! MPI

  SAFE_ALLOCATE(atyp, (nat))
  SAFE_ALLOCATE(apos, (3, nat))
  if (wfnflag) then
     SAFE_ALLOCATE(ngk, (nk))
     SAFE_ALLOCATE(kw, (nk))
     SAFE_ALLOCATE(kpt, (3, nk))
     SAFE_ALLOCATE(ifmin, (nk, ns))
     SAFE_ALLOCATE(ifmax, (nk, ns))
     SAFE_ALLOCATE(energies, (nbands, nk, ns))
     SAFE_ALLOCATE(occupations, (nbands, nk, ns))
  endif
#endif
  ! READ

  ! allocatable arrays
  if (peinf%inode .eq. 0) then

     READ_WRITE()(iunit FORMAT) ((apos(ii, iat), ii = 1, 3), atyp(iat), iat = 1, nat)
     if (wfnflag) then
        READ_WRITE()(iunit FORMAT) (ngk(ik), ik = 1, nk)
        READ_WRITE()(iunit FORMAT) (kw(ik), ik = 1, nk)
        READ_WRITE()(iunit FORMAT) ((kpt(ii, ik), ii = 1, 3), ik = 1, nk)
        READ_WRITE()(iunit FORMAT) ((ifmin(ik, is), ik = 1, nk), is = 1, ns)
        READ_WRITE()(iunit FORMAT) ((ifmax(ik, is), ik = 1, nk), is = 1, ns)
        READ_WRITE()(iunit FORMAT) (((energies(ib, ik, is), ib = 1, nbands), ik = 1, nk), is = 1, ns)
        READ_WRITE()(iunit FORMAT) (((occupations(ib, ik, is), ib = 1, nbands), ik = 1, nk), is = 1, ns)

#if defined READ
        if(any(ngk(1:nk) < 1)) then
           call die("ngk out of bounds in wfn")
        endif
        if(maxval(ngk(1:nk)) /= ngkmax) then
           call die("max(ngk(1:nk)) /= ngkmax in wfn")
        endif
        if(ngkmax > ng) then
           call die("ngkmax > ng in wfn")
        endif
        if(any(kw(1:nk) < -TOL_Zero .or. kw(1:nk) > 1d0 + TOL_Zero)) then
           call die("kw out of bounds in wfn")
        endif
        if(abs(sum(kw(1:nk)) - 1d0) > TOL_SMALL) then
           ! write(0,*) 'kweights sum ', sum(kw(1:nk))
           ! write(0,*) kw(1:nk)
           ! call die("kweights do not sum to 1 in wfn")
        endif
        if(any(ifmin(1:nk, 1:ns) < 0 .or. ifmin(1:nk, 1:ns) > ifmax(1:nk, 1:ns))) then
           call die("ifmin out of bounds in wfn")
        endif
        if(any(ifmax(1:nk, 1:ns) < 0 .or. ifmax(1:nk, 1:ns) > nbands)) then
           call die("ifmax out of bounds in wfn")
        endif
        if(any(ifmin(1:nk, 1:ns) == 0 .and. ifmax(1:nk, 1:ns) /= 0)) then
           call die("ifmin can be zero only when ifmax is zero too; out of bounds in wfn")
        endif
        if(any(occupations(1:nbands, 1:nk, 1:ns) < -0.5 .or. occupations(1:nbands, 1:nk, 1:ns) > 1.5)) then
           call die("occupations out of bounds in wfn")
        endif
        if(any(occupations(1:nbands, 1:nk, 1:ns) < -TOL_Zero .or. occupations(1:nbands, 1:nk, 1:ns) > 1 + TOL_Zero)) then
           write(0,'(a)') "WARNING: occupations outside of range [0,1] in wfn"
           ! this is expected for cold smearing and Methfessel-Paxton smearing
        endif
#endif

     endif

  endif

#if defined READ && defined MPI
  if (peinf%npes .gt. 1) then
     call MPI_BCAST(atyp, nat, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
     call MPI_BCAST(apos, 3*nat, MPI_REAL_DP, 0, MPI_COMM_WORLD, mpierr)
     if (wfnflag) then
        call MPI_BCAST(ngk, nk, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        call MPI_BCAST(kw, nk, MPI_REAL_DP, 0, MPI_COMM_WORLD, mpierr)
        call MPI_BCAST(kpt, 3*nk, MPI_REAL_DP, 0, MPI_COMM_WORLD, mpierr)
        call MPI_BCAST(ifmin, ns*nk, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        call MPI_BCAST(ifmax, ns*nk, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        call MPI_BCAST(energies(1,1,1), nbands*ns*nk, MPI_REAL_DP, 0, MPI_COMM_WORLD, mpierr)
        call MPI_BCAST(occupations(1,1,1), nbands*ns*nk, MPI_REAL_DP, 0, MPI_COMM_WORLD, mpierr)
     endif
  endif
#endif

#ifdef READ
  ! converters and utilities should not write a warning for complex with inversion
  warn_ = .not. is_get
  if(present(warn)) warn_ = warn
  call check_inversion(iflavor, ntran, mtrx, ns, warn_, .false., tnp = tnp)
#endif

  POP_SUB(NAMEHEADER)
  return
end subroutine NAMEHEADER

!=========================================================================
!> wrapper routine that uses typedefs types
subroutine NAME(_header_type)(iunit, sheader, iflavor, kp, gvec, syms, crys, warn, dont_warn_kgrid, sdate, stime, version)
  integer, intent(in) :: iunit
  character(len=3), intent(inout) :: sheader
  integer, intent(FLAVOR_INTENT) :: iflavor
  type(kpoints), intent(INTENT) :: kp
  type(gspace), intent(INTENT) :: gvec
  type(symmetry), intent(INTENT) :: syms
  type(crystal), intent(INTENT) :: crys
  logical, optional, intent(in) :: warn !< if false, suppresses warnings about inversion symmetries etc.
  logical, optional, intent(in) :: dont_warn_kgrid !< if true, validity of kgrid will not be checked.
  !! use for inteqp to allow interpolation onto non-uniform fine grids
  character(len=32), optional, intent(out) :: sdate, stime !< if read, result from file is returned; if write, current is returned
  integer, optional, intent(INTENT) :: version !< version of the header. -1 if no version was/should be specified

  PUSH_SUB(NAME(_header_type))

  call NAMEHEADER(iunit, sheader, iflavor, kp%nspin, gvec%ng, syms%ntran, syms%cell_symmetry, crys%nat, &
       kp%nrk, kp%mnband, kp%ngkmax, gvec%ecutrho, kp%ecutwfc, gvec%FFTgrid, kp%kgrid, kp%shift, crys%celvol, &
       crys%alat, crys%avec, crys%adot, crys%recvol, crys%blat, crys%bvec, crys%bdot, syms%mtrx, syms%tnp, &
       crys%atyp, crys%apos, kp%ngk, kp%w, kp%rk, kp%ifmin, kp%ifmax, kp%el, kp%occ, nspinor = kp%nspinor, &
       warn = warn, dont_warn_kgrid = dont_warn_kgrid, sdate = sdate, stime = stime, version = version)

#ifdef READ
  gvec%nFFTgridpts = product(gvec%FFTgrid(1:3))
#endif

  POP_SUB(NAME(_header_type))
  return
end subroutine NAME(_header_type)

#endif
!TEMP_SCALAR

!=========================================================================
!> This routine is turned into 12 variants by preprocessing.
!! If nrecord and ng_record are not provided for read, the routine simply does not tell you what it found.
!! For write, nrecord is given the default 1 if not provided. ng_record defaults to ng if nrecord = 1, but
!! if nrecord > 1, ng_record must be provided to specify how to divide up the g-vectors.
subroutine LONGNAME()(iunit, ng, ng_bound, ARGUMENTS, nrecord, ng_record, bcast, gindex, dont_read)
  integer, intent(in) :: iunit
  integer, intent(in) :: ng       !< used size of array
  integer, intent(in) :: ng_bound !< actual size of array, >= ng
#ifdef GVECTORS
  integer, intent(INTENT) :: gvec(:, :) !< (3, ng_bound)
#else
  integer, intent(in) :: ns
  TEMP_SCALAR, intent(INTENT) :: data(:, :) !< (ng_bound, ns)
#endif
  integer, optional, intent(INTENT) :: nrecord !< data/gvectors will be distributed among this many records
  integer, optional, pointer, intent(INTENT) :: ng_record(:) !< number of gvectors in each record
  !! must be provided for write if nrecord > 1
  logical, optional, intent(in) :: bcast     !< whether to do MPI_Bcast of what is read
  integer, optional, intent(in) :: gindex(:) !< map of order in file to order in gvec
  logical, optional, intent(in) :: dont_read !< if true, records will just be skipped; only for unformatted

  integer :: ig, igg, irecord, ii, nrecord_internal, ng_irecord
  character*100 :: tmpstr
  logical :: bcast_, dont_read_

  PUSH_SUB(LONGNAME())

  if(ng_bound < ng) then
     write(tmpstr,'(3a,i10,a,i10,a)') "In routine ", TOSTRING(LONGNAME()), ", ng_bound = ", ng_bound, " < ", ng, " = ng"
     call die(tmpstr)
  endif

  dont_read_ = .false.
  if(present(dont_read)) then
#ifdef READ
#ifdef FORMATTED
     if(dont_read) call die("Formatted routine " + TOSTRING(LONGNAME()) + " cannot take argument dont_read = .true.")
#endif
     dont_read_ = dont_read
#else
     call die("Write routine " + TOSTRING(LONGNAME()) + " cannot take argument dont_read.")
#endif
  endif

  if(peinf%inode .eq. 0) then
     if(present(gindex) .and. .not. dont_read_) then
        if(ubound(gindex, 1) < ng) &
             call die("In routine " + TOSTRING(LONGNAME()) + ", gindex array is too small")
        if(any(gindex(1:ng) > ng_bound) .or. any(gindex(1:ng) <= 0)) then
           call die("gindex out of bounds in " + TOSTRING(LONGNAME()))
        endif
     endif
  endif

  bcast_ = .not. dont_read_
  ! if not read, no point in broadcasting
  if(present(bcast)) then
#ifdef READ
     bcast_ = bcast
#else
     call die("Write routine " + TOSTRING(LONGNAME()) + " cannot take argument bcast.")
#endif
  endif

  ! if not reading, the size of the array passed is irrelevant
  if(.not. dont_read_) then
#ifdef GVECTORS
     if(ubound(gvec, 1) /= 3)        &
          call die("In routine " + TOSTRING(LONGNAME()) + ", mismatch of dimension 1 for G-vectors array")
     if(ubound(gvec, 2) < ng_bound)  &
          call die("In routine " + TOSTRING(LONGNAME()) + ", ng_bound is larger than dimension 2 for G-vectors array")
#else
     if(ubound(data, 1) /= ng_bound) &
          call die("In routine " + TOSTRING(LONGNAME()) + ", mismatch of dimension 1 for data array")
     if(ubound(data, 2) /= ns)       &
          call die("In routine " + TOSTRING(LONGNAME()) + ", mismatch of dimension 2 for data array")
#endif
  endif

#ifndef READ
  if(.not. present(nrecord)) then
     nrecord_internal = 1
     ng_irecord = ng
  else
     nrecord_internal = nrecord

     !> check validity of information if going to write
     if(nrecord_internal .ne. 1 .and. .not. present(ng_record)) then
        call die("Routine " + TOSTRING(LONGNAME()) + " requires ng_record array if nrecord > 1.")
     endif

     if(present(ng_record)) then
        if(nrecord_internal .ne. sum(ng_record(1:nrecord))) then
           write(tmpstr,'(a, i10, a, i10, a, a)') "sum(ng_record) = ", sum(ng_record), " ! = ", nrecord , &
                " = nrecord in arguments to routine ", TOSTRING(LONGNAME())
           call die(tmpstr)
        endif
     endif
  endif
#endif

  if(peinf%inode .eq. 0) then
     READ_WRITE()(iunit FORMAT) nrecord_internal
     if(nrecord_internal < 0) then
        write(tmpstr,'(3a, i10)') "In routine ", TOSTRING(LONGNAME()), " illegal nrecord ", nrecord_internal
        call die(tmpstr)
     endif
  endif

#ifdef READ
  if(present(nrecord)) then
     nrecord = nrecord_internal
#ifdef MPI
     if(peinf%npes > 1) call MPI_BCAST(nrecord, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
#endif
  endif
  if(present(ng_record)) then
     SAFE_ALLOCATE(ng_record, (nrecord_internal))
  endif
#endif

  if (peinf%inode.eq.0) call timacc(81,1)

  ! FHJ: zero buffers
#ifdef READ
#ifdef GVECTORS
  gvec(:,:) = 0
#else
  data(:,:) = ZERO
#endif
#endif

  if(peinf%inode .eq. 0) then
     ig = 1
     !This is a PGI pragma to force the optimization level of this routine to -O1.
     !With -O2 or higher (including -fast) the line after the pragma causes a segmentation fault.
     !pgi$r opt=1
     do irecord = 1, nrecord_internal
#ifndef READ
        if(present(ng_record)) ng_irecord = ng_record(irecord)
#endif
        READ_WRITE()(iunit FORMAT) ng_irecord
#ifdef READ
        if(present(ng_record)) ng_record(irecord) = ng_irecord
#endif
        if(dont_read_) then
           READ_WRITE()(iunit FORMAT)
        else
#ifdef GVECTORS
           if(present(gindex)) then
              READ_WRITE()(iunit FORMAT) ((gvec(ii, gindex(igg)), ii = 1, 3), igg = ig, ig + ng_irecord - 1)
           else
              READ_WRITE()(iunit FORMAT) ((gvec(ii, igg), ii = 1, 3), igg = ig, ig + ng_irecord - 1)
           endif
#else
           if(present(gindex)) then
              READ_WRITE()(iunit FORMAT) ((data(gindex(igg), ii), igg = ig, ig + ng_irecord - 1), ii = 1, ns)
           else
              READ_WRITE()(iunit FORMAT) ((data(igg, ii), igg = ig, ig + ng_irecord - 1), ii = 1, ns)
           endif
#endif
        endif
        ig = ig + ng_irecord
     enddo
  endif

  if (peinf%inode.eq.0) call timacc(81,2)
  if (peinf%inode.eq.0) call timacc(82,1)

#if defined READ && defined MPI
  if(peinf%npes > 1) then
     if(present(ng_record)) call MPI_BCAST(ng_record, nrecord_internal, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
     if(bcast_) then
#ifdef GVECTORS
        call MPI_BCAST(gvec, 3 * ng, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
        ! same effect, for the actual data, as 3 * ng_bound
#else
        call MPI_BCAST(data, ng_bound * ns, MPI_TEMP_SCALAR, 0, MPI_COMM_WORLD, mpierr)
#endif
     endif
  endif
#endif

  if (peinf%inode.eq.0) call timacc(82,2)

  POP_SUB(LONGNAME())
  return
end subroutine LONGNAME()

#else
! defined FORMATTED || defined BINARY

#ifndef TEMP_SCALAR
!> wrappers allowing selection of formatted/binary by argument
!=========================================================================
subroutine NAMEHEADER(iunit, is_formatted, sheader, iflavor, ns, ng, ntran, cell_symmetry, &
     nat, nk, nbands, ngkmax, ecutrho, ecutwfc, FFTgrid, kgrid, kshift, celvol, alat, avec, adot, &
     recvol, blat, bvec, bdot, mtrx, tnp, atyp, apos, ngk, kw, kpt, ifmin, ifmax, energies, occupations, &
     nspinor, warn, sdate, stime, version)
  integer, intent(in) :: iunit !< unit number
  logical, intent(in) :: is_formatted !< true = formatted, false = binary
  character(len=3), intent(inout) :: sheader !< file header 'WFN'/'RHO'/'VXC'/'GET'
  integer, intent(FLAVOR_INTENT) :: iflavor !< define type. always must be initialized. modified only if -1 on input
  integer, intent(INTENT) :: ns, ng, ntran, cell_symmetry, nat, nk, nbands, ngkmax !< numbers of spins, G-vectors, symmetries,
  real(DP), intent(INTENT) :: ecutrho, ecutwfc !< charge-density and wave-function cutoffs, in Ry
  integer, intent(INTENT) :: FFTgrid(3), kgrid(3) !< FFT grid size, k-grid size
  real(DP), intent(INTENT) :: kshift(3) !< k-grid offset
  real(DP), intent(INTENT) :: celvol, alat, avec(3, 3), adot(3, 3) !< cell volume, lattice constant,
  real(DP), intent(INTENT) :: recvol, blat, bvec(3, 3), bdot(3, 3) !< cell volume, lattice constant,
  integer, intent(INTENT) :: mtrx(3, 3, 48) !< symmetry matrix
  real(DP), intent(INTENT) :: tnp(3, 48) !< fractional translation
  integer, pointer, intent(INTENT) :: atyp(:) !< atomic species
  real(DP), pointer, intent(INTENT) :: apos(:,:) !< atomic positions (in units of alat)
  integer, pointer, intent(INTENT) :: ngk(:) !< number of G-vectors for each k-pt, ngk(nk)
  real(DP), pointer, intent(INTENT) :: kw(:), kpt(:, :) !< k-weight, kw(nk); k-coord, kpt(3, nk) in crystal coords
  integer, pointer, intent(INTENT) :: ifmin(:, :), ifmax(:, :) !< lowest and highest occupied band, ifmin/max(nk, ns)
  real(DP), pointer, intent(INTENT) :: energies(:, :, :) !< energies(nbands, nk, ns) in Ry
  real(DP), pointer, intent(INTENT) :: occupations(:, :, :) !< occupations(nbands, nk, ns) between 0 and 1
  integer, intent(INTENT) :: nspinor !< 2 if doing a two-component spinor calculation; 1 if not present
  logical, optional, intent(in) :: warn !< if false, suppresses warnings about inversion symmetries etc.
  character(len=32), optional, intent(out) :: sdate, stime !< if read, result from file is returned; if write, current is returned
  integer, optional, intent(INTENT) :: version !< version of the header. -1 if no version was/should be specified

  PUSH_SUB(NAMEHEADER)

  if(is_formatted) then
     call NAME(_format_header)(iunit, sheader, iflavor, ns, ng, ntran, cell_symmetry, nat, nk, nbands, ngkmax, &
          ecutrho, ecutwfc, FFTgrid, kgrid, kshift, celvol, alat, avec, adot, recvol, blat, bvec, bdot, mtrx, tnp, &
          atyp, apos, ngk, kw, kpt, ifmin, ifmax, energies, occupations, &
          nspinor = nspinor, warn = warn, sdate = sdate, stime = stime, version=version)
  else
     call NAME(_binary_header)(iunit, sheader, iflavor, ns, ng, ntran, cell_symmetry, nat, nk, nbands, ngkmax, &
          ecutrho, ecutwfc, FFTgrid, kgrid, kshift, celvol, alat, avec, adot, recvol, blat, bvec, bdot, mtrx, tnp, &
          atyp, apos, ngk, kw, kpt, ifmin, ifmax, energies, occupations, &
          nspinor = nspinor, warn = warn, sdate = sdate, stime = stime, version=version)
  endif

  POP_SUB(NAMEHEADER)
  return
end subroutine NAMEHEADER

!=========================================================================
!> wrapper routine that uses typedefs types
subroutine NAME(_header_type)(iunit, is_formatted, sheader, iflavor, kp, gvec, syms, crys, &
     warn, dont_warn_kgrid, sdate, stime, version)
  integer, intent(in) :: iunit
  logical, intent(in) :: is_formatted !< true = formatted, false = binary
  character(len=3), intent(inout) :: sheader
  integer, intent(FLAVOR_INTENT) :: iflavor
  type(kpoints), intent(INTENT) :: kp
  type(gspace), intent(INTENT) :: gvec
  type(symmetry), intent(INTENT) :: syms
  type(crystal), intent(INTENT) :: crys
  logical, optional, intent(in) :: warn !< if false, suppresses warnings about inversion symmetries etc.
  logical, optional, intent(in) :: dont_warn_kgrid !< if true, validity of kgrid will not be checked.
  !! use for inteqp to allow interpolation onto non-uniform fine grids
  character(len=32), optional, intent(out) :: sdate, stime !< if read, result from file is returned; if write, current is returned
  integer, optional, intent(INTENT) :: version !< version of the header. -1 if no version was/should be specified

  PUSH_SUB(NAME(_header_type))

  if(is_formatted) then
     call NAME(_format_header_type)(iunit, sheader, iflavor, kp, gvec, syms, crys, warn = warn, &
          dont_warn_kgrid = dont_warn_kgrid, sdate = sdate, stime = stime, version=version)
  else
     call NAME(_binary_header_type)(iunit, sheader, iflavor, kp, gvec, syms, crys, warn = warn, &
          dont_warn_kgrid = dont_warn_kgrid, sdate = sdate, stime = stime, version=version)
  endif

  POP_SUB(NAME(_header_type))
  return
end subroutine NAME(_header_type)

#endif
! defined TEMP_SCALAR

!=========================================================================
subroutine LONGNAME()(iunit, is_formatted, ng, ng_bound, ARGUMENTS, nrecord, ng_record, bcast, gindex, dont_read)
  integer, intent(in) :: iunit
  logical, intent(in) :: is_formatted !< true = formatted, false = binary
  integer, intent(in) :: ng       !< used size of array
  integer, intent(in) :: ng_bound !< actual size of array, >= ng
#ifdef GVECTORS
  integer, intent(INTENT) :: gvec(:, :) !< (3, ng_bound)
#else
  integer, intent(in) :: ns
  TEMP_SCALAR, intent(INTENT) :: data(:, :) !< (ng_bound, ns)
#endif
  integer, optional, intent(INTENT) :: nrecord !< data/gvectors will be distributed among this many records
  integer, optional, pointer, intent(INTENT) :: ng_record(:) !< number of gvectors in each record
  logical, optional, intent(in) :: bcast     !< whether to do MPI_Bcast of what is read
  integer, optional, intent(in) :: gindex(:) !< map of order in file to order in gvec
  logical, optional, intent(in) :: dont_read !< if true, records will just be skipped; only for unformatted

  PUSH_SUB(LONGNAME())

  ! note: surprisingly, the Fortran standard permits passing optional arguments to other routines
  ! in which those arguments are also optional.
  if(is_formatted) then
     call LONGNAME(_format)(iunit, ng, ng_bound, ARGUMENTS, nrecord = nrecord, &
          ng_record = ng_record, bcast = bcast, gindex = gindex, dont_read = dont_read)
  else
     call LONGNAME(_binary)(iunit, ng, ng_bound, ARGUMENTS, nrecord = nrecord, &
          ng_record = ng_record, bcast = bcast, gindex = gindex, dont_read = dont_read)
  endif

  POP_SUB(LONGNAME())
  return
end subroutine LONGNAME()

#endif
! defined FORMATTED || defined BINARY

! these undefs prevent lots of cpp warnings
#undef READ_WRITE
#undef INTENT
#undef FLAVOR_INTENT
#undef NAME
#undef FORMAT
#undef NAMEHEADER
#undef GVECTORS
#undef LONGNAME
#undef TEMP_SCALAR
#undef MPI_TEMP_SCALAR
#undef ARGUMENTS

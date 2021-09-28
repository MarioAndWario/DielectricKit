!>=========================================================================
!!
!! Module:
!!
!! (1) wfn_rho_vxc_io_m Originally by DAS Last Modified 10/17/2010 (DAS)
!!
!! Routines to read and write wavefunctions, density, and Vxc, and
!! deallocate associated variables. The "type" routines use typedefs.
!! The code is generated through repeated inclusion of a file with
!! different preprocessor definitions each time. You are not expected to
!! understand this. Consult the resulting .p.f file for clarity.
!!
!! Specification for ASCII matrix element files (vxc.dat, x.dat):
!! Matrix elements are in eV and are always written with real and imaginary parts.
!! They may contain any number of k-points in any order.
!! Each k-point block begins with the line:
!! kx, ky, kz [crystal coordinates], ndiag*nspin, noffdiag*nspin
!! There are then ndiag*nspin lines of the form
!! ispin, idiag, Re <idiag|V|idiag>, Im <idiag|V|idiag>
!! There are then noffdiag*nspin lines of the form
!! ispin, ioff1, ioff2, Re <ioff1|V|ioff2>, Im <ioff1|V|ioff2>
!!
!!=========================================================================
module wfn_rho_vxc_io_m
  use global_m
  ! use check_inversion_m
  use sort_m
  use misc_m
  implicit none
  private
  !> For library usage, do not make global_m contents available
  !! to avoid namespace clashes.
  public :: &
       read_binary_header, &
       write_binary_header, &
       read_format_header, &
       write_format_header, &
       read_binary_header_type, &
       write_binary_header_type, &
       read_format_header_type, &
       write_format_header_type, &
       read_binary_gvectors, &
       write_binary_gvectors, &
       read_format_gvectors, &
       write_format_gvectors, &
       read_binary_real_data, &
       write_binary_real_data, &
       read_format_real_data, &
       write_format_real_data, &
       read_binary_complex_data, &
       write_binary_complex_data, &
       read_format_complex_data, &
       write_format_complex_data, &
       read_binary_data, &
       write_binary_data, &
       read_format_data, &
       write_format_data, &
       read_header, &
       write_header, &
       read_header_type, &
       write_header_type, &
       read_gvectors, &
       write_gvectors, &
       read_real_data, &
       write_real_data, &
       read_complex_data, &
       write_complex_data, &
       read_data, &
       write_data, &
       dealloc_header, &
       dealloc_header_type, &
       dealloc_crys, &
       dealloc_kp, &
       check_header, &
       require_version, &
       write_mf_header, &
       read_mf_header, &
       init_mf_header_from_types
  !> These interfaces can be used to avoid explicit ifdef 1 switches
  interface read_binary_data
     module procedure read_binary_real_data, read_binary_complex_data
  end interface read_binary_data
  interface write_binary_data
     module procedure write_binary_real_data, write_binary_complex_data
  end interface write_binary_data
  interface read_format_data
     module procedure read_format_real_data, read_format_complex_data
  end interface read_format_data
  interface write_format_data
     module procedure write_format_real_data, write_format_complex_data
  end interface write_format_data
  interface read_data
     module procedure read_real_data, read_complex_data
  end interface read_data
  interface write_data
     module procedure write_real_data, write_complex_data
  end interface write_data
contains
  !=========================================================================
  !
  ! Included from file wfn_rho_vxc_io.f90.
  ! You are not expected to understand this. --DAS
  !
  !=========================================================================
  !=========================================================================
  subroutine read_format_header(iunit, sheader, iflavor, ns, ng, ntran, cell_symmetry, &
       nat, nk, nbands, ngkmax, ecutrho, ecutwfc, FFTgrid, kgrid, kshift, celvol, &
       alat, avec, adot, recvol, blat, bvec, bdot, mtrx, tnp, atyp, apos, ngk, &
       kw, kpt, ifmin, ifmax, energies, occupations, nspinor, warn, dont_warn_kgrid, sdate, stime, version)
    integer, intent(in) :: iunit !< unit number
    character(len=3), intent(inout) :: sheader !< file header 'WFN'/'RHO'/'VXC'/'GET' -- last one is to read, and return it
    integer, intent(inout) :: iflavor !< define type. always must be initialized. modified only if -1 on input
    !! -1 = read from file and return it, 0 = as defined by -DCPLX, 1 = real, 2 = complex
    integer, intent(out) :: ns, ng, ntran, cell_symmetry, nat, nk, nbands, ngkmax !< numbers of spins, G-vectors, symmetries,
    !! cell type (0 = cubic, 1 = hexagonal), numbers of atoms, k-points, bands, max(ngk)
    real(DP), intent(out) :: ecutrho, ecutwfc !< charge-density and wave-function cutoffs, in Ry
    integer, intent(out) :: FFTgrid(3), kgrid(3) !< FFT grid size, k-grid size
    real(DP), intent(out) :: kshift(3) !< k-grid offset
    real(DP), intent(out) :: celvol, alat, avec(3, 3), adot(3, 3) !< cell volume, lattice constant,
    !! lattice vectors, metric tensor in real space (in a.u., avec in units of alat)
    real(DP), intent(out) :: recvol, blat, bvec(3, 3), bdot(3, 3) !< cell volume, lattice constant,
    !! lattice vectors, metric tensor in reciprocal space (in a.u., bvec in units of blat)
    integer, intent(out) :: mtrx(3, 3, 48) !< symmetry matrix
    real(DP), intent(out) :: tnp(3, 48) !< fractional translation
    integer, pointer, intent(out) :: atyp(:) !< atomic species
    real(DP), pointer, intent(out) :: apos(:,:) !< atomic positions (in units of alat)
    integer, pointer, intent(out) :: ngk(:) !< number of G-vectors for each k-pt, ngk(nk)
    real(DP), pointer, intent(out) :: kw(:), kpt(:, :) !< k-weight, kw(nk); k-coord, kpt(3, nk) in crystal coords
    integer, pointer, intent(out) :: ifmin(:, :), ifmax(:, :) !< lowest and highest occupied band, ifmin/max(nk, ns)
    real(DP), pointer, intent(out) :: energies(:, :, :) !< energies(nbands, nk, ns) in Ry
    real(DP), pointer, intent(out) :: occupations(:, :, :) !< occupations(nbands, nk, ns) between 0 and 1
    integer, optional, intent(out) :: nspinor !< 2 if doing a two-component spinor calculation; 1 if not present
    logical, optional, intent(in) :: warn !< if false, suppresses warnings about inversion symmetries etc.
    logical, optional, intent(in) :: dont_warn_kgrid !< if true, validity of kgrid read will not be checked.
    !! use for inteqp to allow interpolation onto non-uniform fine grids
    character(len=32), optional, intent(out) :: sdate, stime !< if read, result from file is returned; if write, current is returned
    integer, optional, intent(out) :: version !< version of the header. -1 if no version was/should be specified
    character :: sdate_*32, stime_*32, stitle*32, sflavor*7
    integer :: ii, jj, is, ib, ik, itran, iat, ierr, nspin
    logical :: wfnflag, readwrite_version
    logical :: is_get, warn_, warn_kgrid, die_spinors

    if (sheader .eq. 'WFN') then
       wfnflag = .true.
    elseif (sheader .eq. 'RHO' .or. sheader .eq. 'VXC' .or. sheader .eq. 'KER') then
       wfnflag = .false.
    elseif (sheader .ne. 'GET') then
       call die("Unknown file header '" + sheader + "' passed to " + "read_format_header" &
            + ": must be 'WFN'/'RHO'/'VXC'/'KER'/'GET'")
    endif
    if (peinf%inode .eq. 0) then
       select case(iflavor)
       case(-1)
       case(0)
          sflavor = "Complex"
       case(1)
          sflavor = "Real"
       case(2)
          sflavor = "Complex"
       case default
          write(sflavor, '(i7)') iflavor
          call die("Illegal value iflavor = " + TRUNC(sflavor) + " passed to " + "read_format_header" &
               + ": must be -1,0,1,2.")
       end select
       ! FHJ: we try to the version tag if the argument "version" was passed,
       ! and we only WRITE the version if the argument is present and /= -1.
       readwrite_version=.false.
       if (present(version)) then
          readwrite_version = .true.
       endif
       ierr = -1
       ! FHJ: try to read/write a first time, including the version tag
       if (readwrite_version) then
          read(iunit , *, iostat = ierr) version, stitle, sdate_, stime_
       endif
       ! FHJ: if there was no version tag, or if we don`t care about version...
       if (ierr /= 0) then
          ! FHJ: go back one record if we tried reading before
          if (readwrite_version) backspace(iunit, iostat = ierr)
          read(iunit , *, iostat = ierr) stitle, sdate_, stime_
          ! ! ======
          ! write(*,*) "stitle = ", stitle, "sdate = ", sdate_, "stime_ = ", stime_
          ! ! ------
          if(ierr /= 0) then
             call die("Failed " + "read" + " operation in " + "read_format_header" &
                  + " on header record in mode " + sheader + ".")
          endif
          if (present(version)) version = -1
       endif
       is_get = .false.
       if(sheader == 'GET') then
          sheader = stitle(1:3)
          is_get = .true.
          wfnflag = (sheader == 'WFN')
       endif
       if(iflavor == -1) then
          sflavor = TRUNC(stitle(5:))
          if(TRUNC(sflavor) .eq. "Real") then
             iflavor = 1
          else if (TRUNC(sflavor) .eq. "Complex") then
             iflavor = 2
          else
             call die("Read unknown flavor '" + TRUNC(sflavor) + "' in " + &
                  "read_format_header" + ": must be '" + "Real" + "'/'" + "Complex" + "'")
          endif
       endif
       if(TRUNC(stitle) .ne. TRUNC(sheader + "-" + sflavor)) then
          call die("File header mismatch: read '" + TRUNC(stitle) + "' but expected '" &
               + sheader + "-" + TRUNC(sflavor) + "'")
       endif
       ! scalar variables
       if(wfnflag) then
          read(iunit , *) nspin, ng, ntran, cell_symmetry, nat, ecutrho, nk, nbands, ngkmax, ecutwfc
          if(nk < 1) call die("nk out of bounds in wfn")
          if(nbands < 1) call die("nbands out of bounds in wfn")
          if(ngkmax < 1) call die("ngkmax out of bounds in wfn")
          if(ecutwfc < -TOL_Zero) call die("ecutwfc out of bounds in wfn")
       else
          read(iunit , *) nspin, ng, ntran, cell_symmetry, nat, ecutrho
       endif
       if(nspin /= 1 .and. nspin /= 2 .and. nspin /= 4) call die("nspin out of bounds")
       if(nspin == 4) then
          if(.not. present(nspinor)) call die("nspinor not passed but read nspin==4")
          nspinor = 2
          ns = 1
          ! die_spinors = .true.
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
       ! arrays of fixed size
       if(wfnflag) then
          read(iunit , *) (FFTgrid(ii), ii = 1, 3), (kgrid(ii), ii = 1, 3), (kshift(ii), ii = 1, 3)
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
       else
          read(iunit , *) (FFTgrid(ii), ii = 1, 3)
       endif
       read(iunit , *) celvol, alat, ((avec(jj, ii), jj = 1, 3), ii = 1, 3), ((adot(jj, ii), jj = 1, 3), ii = 1, 3)
       read(iunit , *) recvol, blat, ((bvec(jj, ii), jj = 1, 3), ii = 1, 3), ((bdot(jj, ii), jj = 1, 3), ii = 1, 3)
       mtrx(:,:,:) = 0
       tnp(:,:) = 0d0
       read(iunit , *) (((mtrx(ii, jj, itran), ii = 1, 3), jj = 1, 3), itran = 1, ntran)
       read(iunit , *) ((tnp(jj, itran), jj = 1, 3), itran = 1, ntran)
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
    endif
    if(present(sdate)) sdate = sdate_
    if(present(stime)) stime = stime_
    if(peinf%npes > 1) then
       if(present(version)) call MPI_BCAST(version, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
       call MPI_BCAST(iflavor, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
       if(present(sdate)) call MPI_BCAST(sdate, 32, MPI_CHARACTER, 0, MPI_COMM_WORLD, mpierr)
       if(present(stime)) call MPI_BCAST(stime, 32, MPI_CHARACTER, 0, MPI_COMM_WORLD, mpierr)
    endif
    if (peinf%npes .gt. 1) then
       call MPI_BCAST(is_get, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, mpierr)
       if(is_get) call MPI_BCAST(sheader, 3, MPI_CHARACTER, 0, MPI_COMM_WORLD, mpierr)
       call MPI_BCAST(ns, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
       if(present(nspinor)) call MPI_BCAST(nspinor, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
       call MPI_BCAST(ng, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
       call MPI_BCAST(ntran, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
       call MPI_BCAST(cell_symmetry, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
       call MPI_BCAST(nat, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
       call MPI_BCAST(ecutrho, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
       call MPI_BCAST(FFTgrid, 3, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
       call MPI_BCAST(celvol, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
       call MPI_BCAST(alat, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
       call MPI_BCAST(avec, 9, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
       call MPI_BCAST(adot, 9, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
       call MPI_BCAST(recvol, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
       call MPI_BCAST(blat, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
       call MPI_BCAST(bvec, 9, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
       call MPI_BCAST(bdot, 9, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
       call MPI_BCAST(mtrx(1,1,1), 3*3*48, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
       call MPI_BCAST(tnp, 3*48, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
       call MPI_BCAST(wfnflag, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, mpierr)
       if (wfnflag) then
          call MPI_BCAST(nk, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
          call MPI_BCAST(nbands, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
          call MPI_BCAST(ngkmax, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
          call MPI_BCAST(ecutwfc, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
          call MPI_BCAST(kgrid, 3, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
          call MPI_BCAST(kshift, 3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
       endif
    endif
    ! 1
    allocate(atyp (nat))
    allocate(apos (3, nat))
    if (wfnflag) then
       allocate(ngk (nk))
       allocate(kw (nk))
       allocate(kpt (3, nk))
       allocate(ifmin (nk, ns))
       allocate(ifmax (nk, ns))
       allocate(energies (nbands, nk, ns))
       allocate(occupations (nbands, nk, ns))
    endif
    !
    ! allocatable arrays
    if (peinf%inode .eq. 0) then
       read(iunit , *) ((apos(ii, iat), ii = 1, 3), atyp(iat), iat = 1, nat)
       if (wfnflag) then
          read(iunit , *) (ngk(ik), ik = 1, nk)
          read(iunit , *) (kw(ik), ik = 1, nk)
          read(iunit , *) ((kpt(ii, ik), ii = 1, 3), ik = 1, nk)
          read(iunit , *) ((ifmin(ik, is), ik = 1, nk), is = 1, ns)
          read(iunit , *) ((ifmax(ik, is), ik = 1, nk), is = 1, ns)
          read(iunit , *) (((energies(ib, ik, is), ib = 1, nbands), ik = 1, nk), is = 1, ns)
          read(iunit , *) (((occupations(ib, ik, is), ib = 1, nbands), ik = 1, nk), is = 1, ns)
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
       endif
    endif
    if (peinf%npes .gt. 1) then
       call MPI_BCAST(atyp, nat, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
       call MPI_BCAST(apos, 3*nat, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
       if (wfnflag) then
          call MPI_BCAST(ngk, nk, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
          call MPI_BCAST(kw, nk, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
          call MPI_BCAST(kpt, 3*nk, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
          call MPI_BCAST(ifmin, ns*nk, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
          call MPI_BCAST(ifmax, ns*nk, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
          call MPI_BCAST(energies(1,1,1), nbands*ns*nk, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
          call MPI_BCAST(occupations(1,1,1), nbands*ns*nk, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
       endif
    endif
    ! converters and utilities should not write a warning for complex with inversion
    warn_ = .not. is_get
    if(present(warn)) warn_ = warn
    ! call check_inversion(iflavor, ntran, mtrx, ns, warn_, .false., tnp = tnp)

    return
  end subroutine read_format_header
  !=========================================================================
  !> wrapper routine that uses typedefs types
  subroutine read_format_header_type(iunit, sheader, iflavor, kp, gvec, syms, crys, warn, dont_warn_kgrid, sdate, stime, version)
    integer, intent(in) :: iunit
    character(len=3), intent(inout) :: sheader
    integer, intent(inout) :: iflavor
    type(kpoints), intent(out) :: kp
    type(gspace), intent(out) :: gvec
    type(symmetry), intent(out) :: syms
    type(crystal), intent(out) :: crys
    logical, optional, intent(in) :: warn !< if false, suppresses warnings about inversion symmetries etc.
    logical, optional, intent(in) :: dont_warn_kgrid !< if true, validity of kgrid will not be checked.
    !! use for inteqp to allow interpolation onto non-uniform fine grids
    character(len=32), optional, intent(out) :: sdate, stime !< if read, result from file is returned; if write, current is returned
    integer, optional, intent(out) :: version !< version of the header. -1 if no version was/should be specified

    call read_format_header(iunit, sheader, iflavor, kp%nspin, gvec%ng, syms%ntran, syms%cell_symmetry, crys%nat, &
         kp%nrk, kp%mnband, kp%ngkmax, gvec%ecutrho, kp%ecutwfc, gvec%FFTgrid, kp%kgrid, kp%shift, crys%celvol, &
         crys%alat, crys%avec, crys%adot, crys%recvol, crys%blat, crys%bvec, crys%bdot, syms%mtrx, syms%tnp, &
         crys%atyp, crys%apos, kp%ngk, kp%w, kp%rk, kp%ifmin, kp%ifmax, kp%el, kp%occ, nspinor = kp%nspinor, &
         warn = warn, dont_warn_kgrid = dont_warn_kgrid, sdate = sdate, stime = stime, version = version)
    gvec%nFFTgridpts = product(gvec%FFTgrid(1:3))

    return
  end subroutine read_format_header_type
  !TEMP_SCALAR
  !=========================================================================
  !> This routine is turned into 12 variants by preprocessing.
  !! If nrecord and ng_record are not provided for read, the routine simply does not tell you what it found.
  !! For write, nrecord is given the default 1 if not provided. ng_record defaults to ng if nrecord = 1, but
  !! if nrecord > 1, ng_record must be provided to specify how to divide up the g-vectors.
  subroutine read_format_gvectors(iunit, ng, ng_bound, gvec, nrecord, ng_record, bcast, gindex, dont_read)
    integer, intent(in) :: iunit
    integer, intent(in) :: ng !< used size of array
    integer, intent(in) :: ng_bound !< actual size of array, >= ng
    integer, intent(out) :: gvec(:, :) !< (3, ng_bound)
    integer, optional, intent(out) :: nrecord !< data/gvectors will be distributed among this many records
    integer, optional, pointer, intent(out) :: ng_record(:) !< number of gvectors in each record
    !! must be provided for write if nrecord > 1
    logical, optional, intent(in) :: bcast !< whether to do MPI_Bcast of what is read
    integer, optional, intent(in) :: gindex(:) !< map of order in file to order in gvec
    logical, optional, intent(in) :: dont_read !< if true, records will just be skipped; only for unformatted
    integer :: ig, igg, irecord, ii, nrecord_internal, ng_irecord
    character*100 :: tmpstr
    logical :: bcast_, dont_read_

    if(ng_bound < ng) then
       write(tmpstr,'(3a,i10,a,i10,a)') "In routine ", "read_format_gvectors", ", ng_bound = ", ng_bound, " < ", ng, " = ng"
       call die(tmpstr)
    endif
    dont_read_ = .false.
    if(present(dont_read)) then
       if(dont_read) call die("Formatted routine " + "read_format_gvectors" + " cannot take argument dont_read = .true.")
       dont_read_ = dont_read
    endif
    if(peinf%inode .eq. 0) then
       if(present(gindex) .and. .not. dont_read_) then
          if(ubound(gindex, 1) < ng) &
               call die("In routine " + "read_format_gvectors" + ", gindex array is too small")
          if(any(gindex(1:ng) > ng_bound) .or. any(gindex(1:ng) <= 0)) then
             call die("gindex out of bounds in " + "read_format_gvectors")
          endif
       endif
    endif
    bcast_ = .not. dont_read_
    ! if not read, no point in broadcasting
    if(present(bcast)) then
       bcast_ = bcast
    endif
    ! if not reading, the size of the array passed is irrelevant
    if(.not. dont_read_) then
       if(ubound(gvec, 1) /= 3) &
            call die("In routine " + "read_format_gvectors" + ", mismatch of dimension 1 for G-vectors array")
       if(ubound(gvec, 2) < ng_bound) &
            call die("In routine " + "read_format_gvectors" + ", ng_bound is larger than dimension 2 for G-vectors array")
    endif
    if(peinf%inode .eq. 0) then
       read(iunit , *) nrecord_internal
       if(nrecord_internal < 0) then
          write(tmpstr,'(3a, i10)') "In routine ", "read_format_gvectors", " illegal nrecord ", nrecord_internal
          call die(tmpstr)
       endif
    endif
    if(present(nrecord)) then
       nrecord = nrecord_internal
       if(peinf%npes > 1) call MPI_BCAST(nrecord, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
    endif
    if(present(ng_record)) then
       allocate(ng_record (nrecord_internal))
    endif
    if (peinf%inode.eq.0) call timacc(81,1)
    ! FHJ: zero buffers
    gvec(:,:) = 0
    if(peinf%inode .eq. 0) then
       ig = 1
       !This is a PGI pragma to force the optimization level of this routine to -O1.
       !With -O2 or higher (including -fast) the line after the pragma causes a segmentation fault.
       !pgi$r opt=1
       do irecord = 1, nrecord_internal
          read(iunit , *) ng_irecord
          if(present(ng_record)) ng_record(irecord) = ng_irecord
          if(dont_read_) then
             read(iunit , *)
          else
             if(present(gindex)) then
                read(iunit , *) ((gvec(ii, gindex(igg)), ii = 1, 3), igg = ig, ig + ng_irecord - 1)
             else
                read(iunit , *) ((gvec(ii, igg), ii = 1, 3), igg = ig, ig + ng_irecord - 1)
             endif
          endif
          ig = ig + ng_irecord
       enddo
    endif
    if (peinf%inode.eq.0) call timacc(81,2)
    if (peinf%inode.eq.0) call timacc(82,1)
    if(peinf%npes > 1) then
       if(present(ng_record)) call MPI_BCAST(ng_record, nrecord_internal, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
       if(bcast_) then
          call MPI_BCAST(gvec, 3 * ng, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
          ! same effect, for the actual data, as 3 * ng_bound
       endif
    endif
    if (peinf%inode.eq.0) call timacc(82,2)

    return
  end subroutine read_format_gvectors
  ! defined || defined BINARY
  ! these undefs prevent lots of cpp warnings
  !> write
  !=========================================================================
  !
  ! Included from file wfn_rho_vxc_io.f90.
  ! You are not expected to understand this. --DAS
  !
  !=========================================================================
  !=========================================================================
  subroutine write_format_header(iunit, sheader, iflavor, ns, ng, ntran, cell_symmetry, &
       nat, nk, nbands, ngkmax, ecutrho, ecutwfc, FFTgrid, kgrid, kshift, celvol, &
       alat, avec, adot, recvol, blat, bvec, bdot, mtrx, tnp, atyp, apos, ngk, &
       kw, kpt, ifmin, ifmax, energies, occupations, nspinor, warn, dont_warn_kgrid, sdate, stime, version)
    integer, intent(in) :: iunit !< unit number
    character(len=3), intent(inout) :: sheader !< file header 'WFN'/'RHO'/'VXC'/'GET' -- last one is to read, and return it
    integer, intent(in) :: iflavor !< define type. always must be initialized. modified only if -1 on input
    !! -1 = read from file and return it, 0 = as defined by -DCPLX, 1 = real, 2 = complex
    integer, intent(in) :: ns, ng, ntran, cell_symmetry, nat, nk, nbands, ngkmax !< numbers of spins, G-vectors, symmetries,
    !! cell type (0 = cubic, 1 = hexagonal), numbers of atoms, k-points, bands, max(ngk)
    real(DP), intent(in) :: ecutrho, ecutwfc !< charge-density and wave-function cutoffs, in Ry
    integer, intent(in) :: FFTgrid(3), kgrid(3) !< FFT grid size, k-grid size
    real(DP), intent(in) :: kshift(3) !< k-grid offset
    real(DP), intent(in) :: celvol, alat, avec(3, 3), adot(3, 3) !< cell volume, lattice constant,
    !! lattice vectors, metric tensor in real space (in a.u., avec in units of alat)
    real(DP), intent(in) :: recvol, blat, bvec(3, 3), bdot(3, 3) !< cell volume, lattice constant,
    !! lattice vectors, metric tensor in reciprocal space (in a.u., bvec in units of blat)
    integer, intent(in) :: mtrx(3, 3, 48) !< symmetry matrix
    real(DP), intent(in) :: tnp(3, 48) !< fractional translation
    integer, pointer, intent(in) :: atyp(:) !< atomic species
    real(DP), pointer, intent(in) :: apos(:,:) !< atomic positions (in units of alat)
    integer, pointer, intent(in) :: ngk(:) !< number of G-vectors for each k-pt, ngk(nk)
    real(DP), pointer, intent(in) :: kw(:), kpt(:, :) !< k-weight, kw(nk); k-coord, kpt(3, nk) in crystal coords
    integer, pointer, intent(in) :: ifmin(:, :), ifmax(:, :) !< lowest and highest occupied band, ifmin/max(nk, ns)
    real(DP), pointer, intent(in) :: energies(:, :, :) !< energies(nbands, nk, ns) in Ry
    real(DP), pointer, intent(in) :: occupations(:, :, :) !< occupations(nbands, nk, ns) between 0 and 1
    integer, optional, intent(in) :: nspinor !< 2 if doing a two-component spinor calculation; 1 if not present
    logical, optional, intent(in) :: warn !< if false, suppresses warnings about inversion symmetries etc.
    logical, optional, intent(in) :: dont_warn_kgrid !< if true, validity of kgrid read will not be checked.
    !! use for inteqp to allow interpolation onto non-uniform fine grids
    character(len=32), optional, intent(out) :: sdate, stime !< if read, result from file is returned; if write, current is returned
    integer, optional, intent(in) :: version !< version of the header. -1 if no version was/should be specified
    character :: sdate_*32, stime_*32, stitle*32, sflavor*7
    integer :: ii, jj, is, ib, ik, itran, iat, ierr, nspin
    logical :: wfnflag, readwrite_version
    character :: adate*11, atime*14

    if (sheader .eq. 'WFN') then
       wfnflag = .true.
    elseif (sheader .eq. 'RHO' .or. sheader .eq. 'VXC' .or. sheader .eq. 'KER') then
       wfnflag = .false.
    elseif (sheader .ne. 'GET') then
       call die("Unknown file header '" + sheader + "' passed to " + "write_format_header" &
            + ": must be 'WFN'/'RHO'/'VXC'/'KER'/'GET'")
    endif
    if(sheader .eq. 'GET') call die("Header 'GET' may not be passed to " + "write_format_header" &
         + ": must be 'WFN'/'RHO'/'VXC'")
    if (peinf%inode .eq. 0) then
       select case(iflavor)
       case(-1)
          call die("iflavor = -1 cannot be passed to " + "write_format_header")
       case(0)
          sflavor = "Complex"
       case(1)
          sflavor = "Real"
       case(2)
          sflavor = "Complex"
       case default
          write(sflavor, '(i7)') iflavor
          call die("Illegal value iflavor = " + TRUNC(sflavor) + " passed to " + "write_format_header" &
               + ": must be -1,0,1,2.")
       end select
       call date_time(adate, atime)
       sdate_ = adate
       stime_ = atime
       stitle = sheader + "-" + sflavor
       ! FHJ: we try to READ the version tag if the argument "version" was passed,
       ! and we only WRITE the version if the argument is present and /= -1.
       readwrite_version=.false.
       if (present(version)) then
          readwrite_version = version/=-1
       endif
       ierr = -1
       ! FHJ: try to read/write a first time, including the version tag
       if (readwrite_version) then
          write(iunit , *, iostat = ierr) version, stitle, sdate_, stime_
       endif
       ! FHJ: if there was no version tag, or if we don`t care about version...
       if (ierr /= 0) then
          ! FHJ: go back one record if we tried reading before
          if (readwrite_version) backspace(iunit, iostat = ierr)
          write(iunit , *, iostat = ierr) stitle, sdate_, stime_
          if(ierr /= 0) then
             call die("Failed " + "write" + " operation in " + "write_format_header" &
                  + " on header record in mode " + sheader + ".")
          endif
       endif
       nspin = ns
       if(present(nspinor)) then
          if(nspinor==2) then
             nspin = 4
          endif
       endif
       ! scalar variables
       if(wfnflag) then
          write(iunit , *) nspin, ng, ntran, cell_symmetry, nat, ecutrho, nk, nbands, ngkmax, ecutwfc
       else
          write(iunit , *) nspin, ng, ntran, cell_symmetry, nat, ecutrho
       endif
       ! arrays of fixed size
       if(wfnflag) then
          write(iunit , *) (FFTgrid(ii), ii = 1, 3), (kgrid(ii), ii = 1, 3), (kshift(ii), ii = 1, 3)
       else
          write(iunit , *) (FFTgrid(ii), ii = 1, 3)
       endif
       write(iunit , *) celvol, alat, ((avec(jj, ii), jj = 1, 3), ii = 1, 3), ((adot(jj, ii), jj = 1, 3), ii = 1, 3)
       write(iunit , *) recvol, blat, ((bvec(jj, ii), jj = 1, 3), ii = 1, 3), ((bdot(jj, ii), jj = 1, 3), ii = 1, 3)
       write(iunit , *) (((mtrx(ii, jj, itran), ii = 1, 3), jj = 1, 3), itran = 1, ntran)
       write(iunit , *) ((tnp(jj, itran), jj = 1, 3), itran = 1, ntran)
    endif
    if(present(sdate)) sdate = sdate_
    if(present(stime)) stime = stime_
    if(peinf%npes > 1) then
       if(present(sdate)) call MPI_BCAST(sdate, 32, MPI_CHARACTER, 0, MPI_COMM_WORLD, mpierr)
       if(present(stime)) call MPI_BCAST(stime, 32, MPI_CHARACTER, 0, MPI_COMM_WORLD, mpierr)
    endif
    ! READ
    ! allocatable arrays
    if (peinf%inode .eq. 0) then
       write(iunit , *) ((apos(ii, iat), ii = 1, 3), atyp(iat), iat = 1, nat)
       if (wfnflag) then
          write(iunit , *) (ngk(ik), ik = 1, nk)
          write(iunit , *) (kw(ik), ik = 1, nk)
          write(iunit , *) ((kpt(ii, ik), ii = 1, 3), ik = 1, nk)
          write(iunit , *) ((ifmin(ik, is), ik = 1, nk), is = 1, ns)
          write(iunit , *) ((ifmax(ik, is), ik = 1, nk), is = 1, ns)
          write(iunit , *) (((energies(ib, ik, is), ib = 1, nbands), ik = 1, nk), is = 1, ns)
          write(iunit , *) (((occupations(ib, ik, is), ib = 1, nbands), ik = 1, nk), is = 1, ns)
       endif
    endif

    return
  end subroutine write_format_header
  !=========================================================================
  !> wrapper routine that uses typedefs types
  subroutine write_format_header_type(iunit, sheader, iflavor, kp, gvec, syms, crys, warn, dont_warn_kgrid, sdate, stime, version)
    integer, intent(in) :: iunit
    character(len=3), intent(inout) :: sheader
    integer, intent(in) :: iflavor
    type(kpoints), intent(in) :: kp
    type(gspace), intent(in) :: gvec
    type(symmetry), intent(in) :: syms
    type(crystal), intent(in) :: crys
    logical, optional, intent(in) :: warn !< if false, suppresses warnings about inversion symmetries etc.
    logical, optional, intent(in) :: dont_warn_kgrid !< if true, validity of kgrid will not be checked.
    !! use for inteqp to allow interpolation onto non-uniform fine grids
    character(len=32), optional, intent(out) :: sdate, stime !< if read, result from file is returned; if write, current is returned
    integer, optional, intent(in) :: version !< version of the header. -1 if no version was/should be specified

    call write_format_header(iunit, sheader, iflavor, kp%nspin, gvec%ng, syms%ntran, syms%cell_symmetry, crys%nat, &
         kp%nrk, kp%mnband, kp%ngkmax, gvec%ecutrho, kp%ecutwfc, gvec%FFTgrid, kp%kgrid, kp%shift, crys%celvol, &
         crys%alat, crys%avec, crys%adot, crys%recvol, crys%blat, crys%bvec, crys%bdot, syms%mtrx, syms%tnp, &
         crys%atyp, crys%apos, kp%ngk, kp%w, kp%rk, kp%ifmin, kp%ifmax, kp%el, kp%occ, nspinor = kp%nspinor, &
         warn = warn, dont_warn_kgrid = dont_warn_kgrid, sdate = sdate, stime = stime, version = version)

    return
  end subroutine write_format_header_type
  !TEMP_SCALAR
  !=========================================================================
  !> This routine is turned into 12 variants by preprocessing.
  !! If nrecord and ng_record are not provided for read, the routine simply does not tell you what it found.
  !! For write, nrecord is given the default 1 if not provided. ng_record defaults to ng if nrecord = 1, but
  !! if nrecord > 1, ng_record must be provided to specify how to divide up the g-vectors.
  subroutine write_format_gvectors(iunit, ng, ng_bound, gvec, nrecord, ng_record, bcast, gindex, dont_read)
    integer, intent(in) :: iunit
    integer, intent(in) :: ng !< used size of array
    integer, intent(in) :: ng_bound !< actual size of array, >= ng
    integer, intent(in) :: gvec(:, :) !< (3, ng_bound)
    integer, optional, intent(in) :: nrecord !< data/gvectors will be distributed among this many records
    integer, optional, pointer, intent(in) :: ng_record(:) !< number of gvectors in each record
    !! must be provided for write if nrecord > 1
    logical, optional, intent(in) :: bcast !< whether to do MPI_Bcast of what is read
    integer, optional, intent(in) :: gindex(:) !< map of order in file to order in gvec
    logical, optional, intent(in) :: dont_read !< if true, records will just be skipped; only for unformatted
    integer :: ig, igg, irecord, ii, nrecord_internal, ng_irecord
    character*100 :: tmpstr
    logical :: bcast_, dont_read_

    if(ng_bound < ng) then
       write(tmpstr,'(3a,i10,a,i10,a)') "In routine ", "write_format_gvectors", ", ng_bound = ", ng_bound, " < ", ng, " = ng"
       call die(tmpstr)
    endif
    dont_read_ = .false.
    if(present(dont_read)) then
       call die("Write routine " + "write_format_gvectors" + " cannot take argument dont_read.")
    endif
    if(peinf%inode .eq. 0) then
       if(present(gindex) .and. .not. dont_read_) then
          if(ubound(gindex, 1) < ng) &
               call die("In routine " + "write_format_gvectors" + ", gindex array is too small")
          if(any(gindex(1:ng) > ng_bound) .or. any(gindex(1:ng) <= 0)) then
             call die("gindex out of bounds in " + "write_format_gvectors")
          endif
       endif
    endif
    bcast_ = .not. dont_read_
    ! if not read, no point in broadcasting
    if(present(bcast)) then
       call die("Write routine " + "write_format_gvectors" + " cannot take argument bcast.")
    endif
    ! if not reading, the size of the array passed is irrelevant
    if(.not. dont_read_) then
       if(ubound(gvec, 1) /= 3) &
            call die("In routine " + "write_format_gvectors" + ", mismatch of dimension 1 for G-vectors array")
       if(ubound(gvec, 2) < ng_bound) &
            call die("In routine " + "write_format_gvectors" + ", ng_bound is larger than dimension 2 for G-vectors array")
    endif
    if(.not. present(nrecord)) then
       nrecord_internal = 1
       ng_irecord = ng
    else
       nrecord_internal = nrecord
       !> check validity of information if going to write
       if(nrecord_internal .ne. 1 .and. .not. present(ng_record)) then
          call die("Routine " + "write_format_gvectors" + " requires ng_record array if nrecord > 1.")
       endif
       if(present(ng_record)) then
          if(nrecord_internal .ne. sum(ng_record(1:nrecord))) then
             write(tmpstr,'(a, i10, a, i10, a, a)') "sum(ng_record) = ", sum(ng_record), " ! = ", nrecord , &
                  " = nrecord in arguments to routine ", "write_format_gvectors"
             call die(tmpstr)
          endif
       endif
    endif
    if(peinf%inode .eq. 0) then
       write(iunit , *) nrecord_internal
       if(nrecord_internal < 0) then
          write(tmpstr,'(3a, i10)') "In routine ", "write_format_gvectors", " illegal nrecord ", nrecord_internal
          call die(tmpstr)
       endif
    endif
    if (peinf%inode.eq.0) call timacc(81,1)
    ! FHJ: zero buffers
    if(peinf%inode .eq. 0) then
       ig = 1
       !This is a PGI pragma to force the optimization level of this routine to -O1.
       !With -O2 or higher (including -fast) the line after the pragma causes a segmentation fault.
       !pgi$r opt=1
       do irecord = 1, nrecord_internal
          if(present(ng_record)) ng_irecord = ng_record(irecord)
          write(iunit , *) ng_irecord
          if(dont_read_) then
             write(iunit , *)
          else
             if(present(gindex)) then
                write(iunit , *) ((gvec(ii, gindex(igg)), ii = 1, 3), igg = ig, ig + ng_irecord - 1)
             else
                write(iunit , *) ((gvec(ii, igg), ii = 1, 3), igg = ig, ig + ng_irecord - 1)
             endif
          endif
          ig = ig + ng_irecord
       enddo
    endif
    if (peinf%inode.eq.0) call timacc(81,2)
    if (peinf%inode.eq.0) call timacc(82,1)
    if (peinf%inode.eq.0) call timacc(82,2)

    return
  end subroutine write_format_gvectors
  ! defined || defined BINARY
  ! these undefs prevent lots of cpp warnings
  !=========================================================================
  !
  ! Included from file wfn_rho_vxc_io.f90.
  ! You are not expected to understand this. --DAS
  !
  !=========================================================================
  !=========================================================================
  subroutine read_binary_header(iunit, sheader, iflavor, ns, ng, ntran, cell_symmetry, &
       nat, nk, nbands, ngkmax, ecutrho, ecutwfc, FFTgrid, kgrid, kshift, celvol, &
       alat, avec, adot, recvol, blat, bvec, bdot, mtrx, tnp, atyp, apos, ngk, &
       kw, kpt, ifmin, ifmax, energies, occupations, nspinor, warn, dont_warn_kgrid, sdate, stime, version)
    integer, intent(in) :: iunit !< unit number
    character(len=3), intent(inout) :: sheader !< file header 'WFN'/'RHO'/'VXC'/'GET' -- last one is to read, and return it
    integer, intent(inout) :: iflavor !< define type. always must be initialized. modified only if -1 on input
    !! -1 = read from file and return it, 0 = as defined by -DCPLX, 1 = real, 2 = complex
    integer, intent(out) :: ns, ng, ntran, cell_symmetry, nat, nk, nbands, ngkmax !< numbers of spins, G-vectors, symmetries,
    !! cell type (0 = cubic, 1 = hexagonal), numbers of atoms, k-points, bands, max(ngk)
    real(DP), intent(out) :: ecutrho, ecutwfc !< charge-density and wave-function cutoffs, in Ry
    integer, intent(out) :: FFTgrid(3), kgrid(3) !< FFT grid size, k-grid size
    real(DP), intent(out) :: kshift(3) !< k-grid offset
    real(DP), intent(out) :: celvol, alat, avec(3, 3), adot(3, 3) !< cell volume, lattice constant,
    !! lattice vectors, metric tensor in real space (in a.u., avec in units of alat)
    real(DP), intent(out) :: recvol, blat, bvec(3, 3), bdot(3, 3) !< cell volume, lattice constant,
    !! lattice vectors, metric tensor in reciprocal space (in a.u., bvec in units of blat)
    integer, intent(out) :: mtrx(3, 3, 48) !< symmetry matrix
    real(DP), intent(out) :: tnp(3, 48) !< fractional translation
    integer, pointer, intent(out) :: atyp(:) !< atomic species
    real(DP), pointer, intent(out) :: apos(:,:) !< atomic positions (in units of alat)
    integer, pointer, intent(out) :: ngk(:) !< number of G-vectors for each k-pt, ngk(nk)
    real(DP), pointer, intent(out) :: kw(:), kpt(:, :) !< k-weight, kw(nk); k-coord, kpt(3, nk) in crystal coords
    integer, pointer, intent(out) :: ifmin(:, :), ifmax(:, :) !< lowest and highest occupied band, ifmin/max(nk, ns)
    real(DP), pointer, intent(out) :: energies(:, :, :) !< energies(nbands, nk, ns) in Ry
    real(DP), pointer, intent(out) :: occupations(:, :, :) !< occupations(nbands, nk, ns) between 0 and 1
    integer, optional, intent(out) :: nspinor !< 2 if doing a two-component spinor calculation; 1 if not present
    logical, optional, intent(in) :: warn !< if false, suppresses warnings about inversion symmetries etc.
    logical, optional, intent(in) :: dont_warn_kgrid !< if true, validity of kgrid read will not be checked.
    !! use for inteqp to allow interpolation onto non-uniform fine grids
    character(len=32), optional, intent(out) :: sdate, stime !< if read, result from file is returned; if write, current is returned
    integer, optional, intent(out) :: version !< version of the header. -1 if no version was/should be specified
    character :: sdate_*32, stime_*32, stitle*32, sflavor*7
    integer :: ii, jj, is, ib, ik, itran, iat, ierr, nspin
    logical :: wfnflag, readwrite_version
    logical :: is_get, warn_, warn_kgrid, die_spinors

    if (sheader .eq. 'WFN') then
       wfnflag = .true.
    elseif (sheader .eq. 'RHO' .or. sheader .eq. 'VXC' .or. sheader .eq. 'KER') then
       wfnflag = .false.
    elseif (sheader .ne. 'GET') then
       call die("Unknown file header '" + sheader + "' passed to " + "read_binary_header" &
            + ": must be 'WFN'/'RHO'/'VXC'/'KER'/'GET'")
    endif
    if (peinf%inode .eq. 0) then
       select case(iflavor)
       case(-1)
       case(0)
          sflavor = "Complex"
       case(1)
          sflavor = "Real"
       case(2)
          sflavor = "Complex"
       case default
          write(sflavor, '(i7)') iflavor
          call die("Illegal value iflavor = " + TRUNC(sflavor) + " passed to " + "read_binary_header" &
               + ": must be -1,0,1,2.")
       end select
       ! FHJ: we try to the version tag if the argument "version" was passed,
       ! and we only WRITE the version if the argument is present and /= -1.
       readwrite_version=.false.
       if (present(version)) then
          readwrite_version = .true.
       endif
       ierr = -1
       ! FHJ: try to read/write a first time, including the version tag
       if (readwrite_version) then
          read(iunit , iostat = ierr) version, stitle, sdate_, stime_
       endif
       ! FHJ: if there was no version tag, or if we don`t care about version...
       if (ierr /= 0) then
          ! FHJ: go back one record if we tried reading before
          if (readwrite_version) backspace(iunit, iostat = ierr)
          read(iunit , iostat = ierr) stitle, sdate_, stime_
          if(ierr /= 0) then
             call die("Failed " + "read" + " operation in " + "read_binary_header" &
                  + " on header record in mode " + sheader + ".")
          endif
          if (present(version)) version = -1
       endif
       is_get = .false.
       if(sheader == 'GET') then
          sheader = stitle(1:3)
          is_get = .true.
          wfnflag = (sheader == 'WFN')
       endif
       if(iflavor == -1) then
          sflavor = TRUNC(stitle(5:))
          if(TRUNC(sflavor) .eq. "Real") then
             iflavor = 1
          else if (TRUNC(sflavor) .eq. "Complex") then
             iflavor = 2
          else
             call die("Read unknown flavor '" + TRUNC(sflavor) + "' in " + &
                  "read_binary_header" + ": must be '" + "Real" + "'/'" + "Complex" + "'")
          endif
       endif
       if(TRUNC(stitle) .ne. TRUNC(sheader + "-" + sflavor)) then
          call die("File header mismatch: read '" + TRUNC(stitle) + "' but expected '" &
               + sheader + "-" + TRUNC(sflavor) + "'")
       endif
       ! scalar variables
       if(wfnflag) then
          read(iunit ) nspin, ng, ntran, cell_symmetry, nat, ecutrho, nk, nbands, ngkmax, ecutwfc
          if(nk < 1) call die("nk out of bounds in wfn")
          if(nbands < 1) call die("nbands out of bounds in wfn")
          if(ngkmax < 1) call die("ngkmax out of bounds in wfn")
          if(ecutwfc < -TOL_Zero) call die("ecutwfc out of bounds in wfn")
       else
          read(iunit ) nspin, ng, ntran, cell_symmetry, nat, ecutrho
       endif
       if(nspin /= 1 .and. nspin /= 2 .and. nspin /= 4) call die("nspin out of bounds")
       if(nspin == 4) then
          if(.not. present(nspinor)) call die("nspinor not passed but read nspin==4")
          nspinor = 2
          ns = 1
          ! die_spinors = .true.
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
       ! arrays of fixed size
       if(wfnflag) then
          read(iunit ) (FFTgrid(ii), ii = 1, 3), (kgrid(ii), ii = 1, 3), (kshift(ii), ii = 1, 3)
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
       else
          read(iunit ) (FFTgrid(ii), ii = 1, 3)
       endif
       read(iunit ) celvol, alat, ((avec(jj, ii), jj = 1, 3), ii = 1, 3), ((adot(jj, ii), jj = 1, 3), ii = 1, 3)
       read(iunit ) recvol, blat, ((bvec(jj, ii), jj = 1, 3), ii = 1, 3), ((bdot(jj, ii), jj = 1, 3), ii = 1, 3)
       mtrx(:,:,:) = 0
       tnp(:,:) = 0d0
       read(iunit ) (((mtrx(ii, jj, itran), ii = 1, 3), jj = 1, 3), itran = 1, ntran)
       read(iunit ) ((tnp(jj, itran), jj = 1, 3), itran = 1, ntran)
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
    endif
    if(present(sdate)) sdate = sdate_
    if(present(stime)) stime = stime_
    if(peinf%npes > 1) then
       if(present(version)) call MPI_BCAST(version, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
       call MPI_BCAST(iflavor, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
       if(present(sdate)) call MPI_BCAST(sdate, 32, MPI_CHARACTER, 0, MPI_COMM_WORLD, mpierr)
       if(present(stime)) call MPI_BCAST(stime, 32, MPI_CHARACTER, 0, MPI_COMM_WORLD, mpierr)
    endif
    if (peinf%npes .gt. 1) then
       call MPI_BCAST(is_get, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, mpierr)
       if(is_get) call MPI_BCAST(sheader, 3, MPI_CHARACTER, 0, MPI_COMM_WORLD, mpierr)
       call MPI_BCAST(ns, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
       if(present(nspinor)) call MPI_BCAST(nspinor, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
       call MPI_BCAST(ng, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
       call MPI_BCAST(ntran, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
       call MPI_BCAST(cell_symmetry, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
       call MPI_BCAST(nat, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
       call MPI_BCAST(ecutrho, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
       call MPI_BCAST(FFTgrid, 3, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
       call MPI_BCAST(celvol, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
       call MPI_BCAST(alat, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
       call MPI_BCAST(avec, 9, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
       call MPI_BCAST(adot, 9, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
       call MPI_BCAST(recvol, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
       call MPI_BCAST(blat, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
       call MPI_BCAST(bvec, 9, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
       call MPI_BCAST(bdot, 9, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
       call MPI_BCAST(mtrx(1,1,1), 3*3*48, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
       call MPI_BCAST(tnp, 3*48, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
       call MPI_BCAST(wfnflag, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, mpierr)
       if (wfnflag) then
          call MPI_BCAST(nk, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
          call MPI_BCAST(nbands, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
          call MPI_BCAST(ngkmax, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
          call MPI_BCAST(ecutwfc, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
          call MPI_BCAST(kgrid, 3, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
          call MPI_BCAST(kshift, 3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
       endif
    endif
    ! 1
    allocate(atyp (nat))
    allocate(apos (3, nat))
    if (wfnflag) then
       allocate(ngk (nk))
       allocate(kw (nk))
       allocate(kpt (3, nk))
       allocate(ifmin (nk, ns))
       allocate(ifmax (nk, ns))
       allocate(energies (nbands, nk, ns))
       allocate(occupations (nbands, nk, ns))
    endif
    !
    ! allocatable arrays
    if (peinf%inode .eq. 0) then
       read(iunit ) ((apos(ii, iat), ii = 1, 3), atyp(iat), iat = 1, nat)
       if (wfnflag) then
          read(iunit ) (ngk(ik), ik = 1, nk)
          read(iunit ) (kw(ik), ik = 1, nk)
          read(iunit ) ((kpt(ii, ik), ii = 1, 3), ik = 1, nk)
          read(iunit ) ((ifmin(ik, is), ik = 1, nk), is = 1, ns)
          read(iunit ) ((ifmax(ik, is), ik = 1, nk), is = 1, ns)
          read(iunit ) (((energies(ib, ik, is), ib = 1, nbands), ik = 1, nk), is = 1, ns)
          read(iunit ) (((occupations(ib, ik, is), ib = 1, nbands), ik = 1, nk), is = 1, ns)
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
       endif
    endif
    if (peinf%npes .gt. 1) then
       call MPI_BCAST(atyp, nat, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
       call MPI_BCAST(apos, 3*nat, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
       if (wfnflag) then
          call MPI_BCAST(ngk, nk, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
          call MPI_BCAST(kw, nk, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
          call MPI_BCAST(kpt, 3*nk, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
          call MPI_BCAST(ifmin, ns*nk, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
          call MPI_BCAST(ifmax, ns*nk, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
          call MPI_BCAST(energies(1,1,1), nbands*ns*nk, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
          call MPI_BCAST(occupations(1,1,1), nbands*ns*nk, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
       endif
    endif
    ! converters and utilities should not write a warning for complex with inversion
    warn_ = .not. is_get
    if(present(warn)) warn_ = warn
    ! call check_inversion(iflavor, ntran, mtrx, ns, warn_, .false., tnp = tnp)

    return
  end subroutine read_binary_header
  !=========================================================================
  !> wrapper routine that uses typedefs types
  subroutine read_binary_header_type(iunit, sheader, iflavor, kp, gvec, syms, crys, warn, dont_warn_kgrid, sdate, stime, version)
    integer, intent(in) :: iunit
    character(len=3), intent(inout) :: sheader
    integer, intent(inout) :: iflavor
    type(kpoints), intent(out) :: kp
    type(gspace), intent(out) :: gvec
    type(symmetry), intent(out) :: syms
    type(crystal), intent(out) :: crys
    logical, optional, intent(in) :: warn !< if false, suppresses warnings about inversion symmetries etc.
    logical, optional, intent(in) :: dont_warn_kgrid !< if true, validity of kgrid will not be checked.
    !! use for inteqp to allow interpolation onto non-uniform fine grids
    character(len=32), optional, intent(out) :: sdate, stime !< if read, result from file is returned; if write, current is returned
    integer, optional, intent(out) :: version !< version of the header. -1 if no version was/should be specified

    call read_binary_header(iunit, sheader, iflavor, kp%nspin, gvec%ng, syms%ntran, syms%cell_symmetry, crys%nat, &
         kp%nrk, kp%mnband, kp%ngkmax, gvec%ecutrho, kp%ecutwfc, gvec%FFTgrid, kp%kgrid, kp%shift, crys%celvol, &
         crys%alat, crys%avec, crys%adot, crys%recvol, crys%blat, crys%bvec, crys%bdot, syms%mtrx, syms%tnp, &
         crys%atyp, crys%apos, kp%ngk, kp%w, kp%rk, kp%ifmin, kp%ifmax, kp%el, kp%occ, nspinor = kp%nspinor, &
         warn = warn, dont_warn_kgrid = dont_warn_kgrid, sdate = sdate, stime = stime, version = version)
    gvec%nFFTgridpts = product(gvec%FFTgrid(1:3))

    return
  end subroutine read_binary_header_type
  !TEMP_SCALAR
  !=========================================================================
  !> This routine is turned into 12 variants by preprocessing.
  !! If nrecord and ng_record are not provided for read, the routine simply does not tell you what it found.
  !! For write, nrecord is given the default 1 if not provided. ng_record defaults to ng if nrecord = 1, but
  !! if nrecord > 1, ng_record must be provided to specify how to divide up the g-vectors.
  subroutine read_binary_gvectors(iunit, ng, ng_bound, gvec, nrecord, ng_record, bcast, gindex, dont_read)
    integer, intent(in) :: iunit
    integer, intent(in) :: ng !< used size of array
    integer, intent(in) :: ng_bound !< actual size of array, >= ng
    integer, intent(out) :: gvec(:, :) !< (3, ng_bound)
    integer, optional, intent(out) :: nrecord !< data/gvectors will be distributed among this many records
    integer, optional, pointer, intent(out) :: ng_record(:) !< number of gvectors in each record
    !! must be provided for write if nrecord > 1
    logical, optional, intent(in) :: bcast !< whether to do MPI_Bcast of what is read
    integer, optional, intent(in) :: gindex(:) !< map of order in file to order in gvec
    logical, optional, intent(in) :: dont_read !< if true, records will just be skipped; only for unformatted
    integer :: ig, igg, irecord, ii, nrecord_internal, ng_irecord
    character*100 :: tmpstr
    logical :: bcast_, dont_read_

    if(ng_bound < ng) then
       write(tmpstr,'(3a,i10,a,i10,a)') "In routine ", "read_binary_gvectors", ", ng_bound = ", ng_bound, " < ", ng, " = ng"
       call die(tmpstr)
    endif
    dont_read_ = .false.
    if(present(dont_read)) then
       dont_read_ = dont_read
    endif
    if(peinf%inode .eq. 0) then
       if(present(gindex) .and. .not. dont_read_) then
          if(ubound(gindex, 1) < ng) &
               call die("In routine " + "read_binary_gvectors" + ", gindex array is too small")
          if(any(gindex(1:ng) > ng_bound) .or. any(gindex(1:ng) <= 0)) then
             call die("gindex out of bounds in " + "read_binary_gvectors")
          endif
       endif
    endif
    bcast_ = .not. dont_read_
    ! if not read, no point in broadcasting
    if(present(bcast)) then
       bcast_ = bcast
    endif
    ! if not reading, the size of the array passed is irrelevant
    if(.not. dont_read_) then
       if(ubound(gvec, 1) /= 3) &
            call die("In routine " + "read_binary_gvectors" + ", mismatch of dimension 1 for G-vectors array")
       if(ubound(gvec, 2) < ng_bound) &
            call die("In routine " + "read_binary_gvectors" + ", ng_bound is larger than dimension 2 for G-vectors array")
    endif
    if(peinf%inode .eq. 0) then
       read(iunit ) nrecord_internal
       if(nrecord_internal < 0) then
          write(tmpstr,'(3a, i10)') "In routine ", "read_binary_gvectors", " illegal nrecord ", nrecord_internal
          call die(tmpstr)
       endif
    endif
    if(present(nrecord)) then
       nrecord = nrecord_internal
       if(peinf%npes > 1) call MPI_BCAST(nrecord, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
    endif
    if(present(ng_record)) then
       allocate(ng_record (nrecord_internal))
    endif
    if (peinf%inode.eq.0) call timacc(81,1)
    ! FHJ: zero buffers
    gvec(:,:) = 0
    if(peinf%inode .eq. 0) then
       ig = 1
       !This is a PGI pragma to force the optimization level of this routine to -O1.
       !With -O2 or higher (including -fast) the line after the pragma causes a segmentation fault.
       !pgi$r opt=1
       do irecord = 1, nrecord_internal
          read(iunit ) ng_irecord
          if(present(ng_record)) ng_record(irecord) = ng_irecord
          if(dont_read_) then
             read(iunit )
          else
             if(present(gindex)) then
                read(iunit ) ((gvec(ii, gindex(igg)), ii = 1, 3), igg = ig, ig + ng_irecord - 1)
             else
                read(iunit ) ((gvec(ii, igg), ii = 1, 3), igg = ig, ig + ng_irecord - 1)
             endif
          endif
          ig = ig + ng_irecord
       enddo
    endif
    if (peinf%inode.eq.0) call timacc(81,2)
    if (peinf%inode.eq.0) call timacc(82,1)
    if(peinf%npes > 1) then
       if(present(ng_record)) call MPI_BCAST(ng_record, nrecord_internal, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
       if(bcast_) then
          call MPI_BCAST(gvec, 3 * ng, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
          ! same effect, for the actual data, as 3 * ng_bound
       endif
    endif
    if (peinf%inode.eq.0) call timacc(82,2)

    return
  end subroutine read_binary_gvectors
  ! defined FORMATTED || defined
  ! these undefs prevent lots of cpp warnings
  !> write
  !=========================================================================
  !
  ! Included from file wfn_rho_vxc_io.f90.
  ! You are not expected to understand this. --DAS
  !
  !=========================================================================
  !=========================================================================
  subroutine write_binary_header(iunit, sheader, iflavor, ns, ng, ntran, cell_symmetry, &
       nat, nk, nbands, ngkmax, ecutrho, ecutwfc, FFTgrid, kgrid, kshift, celvol, &
       alat, avec, adot, recvol, blat, bvec, bdot, mtrx, tnp, atyp, apos, ngk, &
       kw, kpt, ifmin, ifmax, energies, occupations, nspinor, warn, dont_warn_kgrid, sdate, stime, version)
    integer, intent(in) :: iunit !< unit number
    character(len=3), intent(inout) :: sheader !< file header 'WFN'/'RHO'/'VXC'/'GET' -- last one is to read, and return it
    integer, intent(in) :: iflavor !< define type. always must be initialized. modified only if -1 on input
    !! -1 = read from file and return it, 0 = as defined by -DCPLX, 1 = real, 2 = complex
    integer, intent(in) :: ns, ng, ntran, cell_symmetry, nat, nk, nbands, ngkmax !< numbers of spins, G-vectors, symmetries,
    !! cell type (0 = cubic, 1 = hexagonal), numbers of atoms, k-points, bands, max(ngk)
    real(DP), intent(in) :: ecutrho, ecutwfc !< charge-density and wave-function cutoffs, in Ry
    integer, intent(in) :: FFTgrid(3), kgrid(3) !< FFT grid size, k-grid size
    real(DP), intent(in) :: kshift(3) !< k-grid offset
    real(DP), intent(in) :: celvol, alat, avec(3, 3), adot(3, 3) !< cell volume, lattice constant,
    !! lattice vectors, metric tensor in real space (in a.u., avec in units of alat)
    real(DP), intent(in) :: recvol, blat, bvec(3, 3), bdot(3, 3) !< cell volume, lattice constant,
    !! lattice vectors, metric tensor in reciprocal space (in a.u., bvec in units of blat)
    integer, intent(in) :: mtrx(3, 3, 48) !< symmetry matrix
    real(DP), intent(in) :: tnp(3, 48) !< fractional translation
    integer, pointer, intent(in) :: atyp(:) !< atomic species
    real(DP), pointer, intent(in) :: apos(:,:) !< atomic positions (in units of alat)
    integer, pointer, intent(in) :: ngk(:) !< number of G-vectors for each k-pt, ngk(nk)
    real(DP), pointer, intent(in) :: kw(:), kpt(:, :) !< k-weight, kw(nk); k-coord, kpt(3, nk) in crystal coords
    integer, pointer, intent(in) :: ifmin(:, :), ifmax(:, :) !< lowest and highest occupied band, ifmin/max(nk, ns)
    real(DP), pointer, intent(in) :: energies(:, :, :) !< energies(nbands, nk, ns) in Ry
    real(DP), pointer, intent(in) :: occupations(:, :, :) !< occupations(nbands, nk, ns) between 0 and 1
    integer, optional, intent(in) :: nspinor !< 2 if doing a two-component spinor calculation; 1 if not present
    logical, optional, intent(in) :: warn !< if false, suppresses warnings about inversion symmetries etc.
    logical, optional, intent(in) :: dont_warn_kgrid !< if true, validity of kgrid read will not be checked.
    !! use for inteqp to allow interpolation onto non-uniform fine grids
    character(len=32), optional, intent(out) :: sdate, stime !< if read, result from file is returned; if write, current is returned
    integer, optional, intent(in) :: version !< version of the header. -1 if no version was/should be specified
    character :: sdate_*32, stime_*32, stitle*32, sflavor*7
    integer :: ii, jj, is, ib, ik, itran, iat, ierr, nspin
    logical :: wfnflag, readwrite_version
    character :: adate*11, atime*14

    if (sheader .eq. 'WFN') then
       wfnflag = .true.
    elseif (sheader .eq. 'RHO' .or. sheader .eq. 'VXC' .or. sheader .eq. 'KER') then
       wfnflag = .false.
    elseif (sheader .ne. 'GET') then
       call die("Unknown file header '" + sheader + "' passed to " + "write_binary_header" &
            + ": must be 'WFN'/'RHO'/'VXC'/'KER'/'GET'")
    endif
    if(sheader .eq. 'GET') call die("Header 'GET' may not be passed to " + "write_binary_header" &
         + ": must be 'WFN'/'RHO'/'VXC'")
    if (peinf%inode .eq. 0) then
       select case(iflavor)
       case(-1)
          call die("iflavor = -1 cannot be passed to " + "write_binary_header")
       case(0)
          sflavor = "Complex"
       case(1)
          sflavor = "Real"
       case(2)
          sflavor = "Complex"
       case default
          write(sflavor, '(i7)') iflavor
          call die("Illegal value iflavor = " + TRUNC(sflavor) + " passed to " + "write_binary_header" &
               + ": must be -1,0,1,2.")
       end select
       call date_time(adate, atime)
       sdate_ = adate
       stime_ = atime
       stitle = sheader + "-" + sflavor
       ! FHJ: we try to READ the version tag if the argument "version" was passed,
       ! and we only WRITE the version if the argument is present and /= -1.
       readwrite_version=.false.
       if (present(version)) then
          readwrite_version = version/=-1
       endif
       ierr = -1
       ! FHJ: try to read/write a first time, including the version tag
       if (readwrite_version) then
          write(iunit , iostat = ierr) version, stitle, sdate_, stime_
       endif
       ! FHJ: if there was no version tag, or if we don`t care about version...
       if (ierr /= 0) then
          ! FHJ: go back one record if we tried reading before
          if (readwrite_version) backspace(iunit, iostat = ierr)
          write(iunit , iostat = ierr) stitle, sdate_, stime_
          if(ierr /= 0) then
             call die("Failed " + "write" + " operation in " + "write_binary_header" &
                  + " on header record in mode " + sheader + ".")
          endif
       endif
       nspin = ns
       if(present(nspinor)) then
          if(nspinor==2) then
             nspin = 4
          endif
       endif
       ! scalar variables
       if(wfnflag) then
          write(iunit ) nspin, ng, ntran, cell_symmetry, nat, ecutrho, nk, nbands, ngkmax, ecutwfc
       else
          write(iunit ) nspin, ng, ntran, cell_symmetry, nat, ecutrho
       endif
       ! arrays of fixed size
       if(wfnflag) then
          write(iunit ) (FFTgrid(ii), ii = 1, 3), (kgrid(ii), ii = 1, 3), (kshift(ii), ii = 1, 3)
       else
          write(iunit ) (FFTgrid(ii), ii = 1, 3)
       endif
       write(iunit ) celvol, alat, ((avec(jj, ii), jj = 1, 3), ii = 1, 3), ((adot(jj, ii), jj = 1, 3), ii = 1, 3)
       write(iunit ) recvol, blat, ((bvec(jj, ii), jj = 1, 3), ii = 1, 3), ((bdot(jj, ii), jj = 1, 3), ii = 1, 3)
       write(iunit ) (((mtrx(ii, jj, itran), ii = 1, 3), jj = 1, 3), itran = 1, ntran)
       write(iunit ) ((tnp(jj, itran), jj = 1, 3), itran = 1, ntran)
    endif
    if(present(sdate)) sdate = sdate_
    if(present(stime)) stime = stime_
    if(peinf%npes > 1) then
       if(present(sdate)) call MPI_BCAST(sdate, 32, MPI_CHARACTER, 0, MPI_COMM_WORLD, mpierr)
       if(present(stime)) call MPI_BCAST(stime, 32, MPI_CHARACTER, 0, MPI_COMM_WORLD, mpierr)
    endif
    ! READ
    ! allocatable arrays
    if (peinf%inode .eq. 0) then
       write(iunit ) ((apos(ii, iat), ii = 1, 3), atyp(iat), iat = 1, nat)
       if (wfnflag) then
          write(iunit ) (ngk(ik), ik = 1, nk)
          write(iunit ) (kw(ik), ik = 1, nk)
          write(iunit ) ((kpt(ii, ik), ii = 1, 3), ik = 1, nk)
          write(iunit ) ((ifmin(ik, is), ik = 1, nk), is = 1, ns)
          write(iunit ) ((ifmax(ik, is), ik = 1, nk), is = 1, ns)
          write(iunit ) (((energies(ib, ik, is), ib = 1, nbands), ik = 1, nk), is = 1, ns)
          write(iunit ) (((occupations(ib, ik, is), ib = 1, nbands), ik = 1, nk), is = 1, ns)
       endif
    endif

    return
  end subroutine write_binary_header
  !=========================================================================
  !> wrapper routine that uses typedefs types
  subroutine write_binary_header_type(iunit, sheader, iflavor, kp, gvec, syms, crys, warn, dont_warn_kgrid, sdate, stime, version)
    integer, intent(in) :: iunit
    character(len=3), intent(inout) :: sheader
    integer, intent(in) :: iflavor
    type(kpoints), intent(in) :: kp
    type(gspace), intent(in) :: gvec
    type(symmetry), intent(in) :: syms
    type(crystal), intent(in) :: crys
    logical, optional, intent(in) :: warn !< if false, suppresses warnings about inversion symmetries etc.
    logical, optional, intent(in) :: dont_warn_kgrid !< if true, validity of kgrid will not be checked.
    !! use for inteqp to allow interpolation onto non-uniform fine grids
    character(len=32), optional, intent(out) :: sdate, stime !< if read, result from file is returned; if write, current is returned
    integer, optional, intent(in) :: version !< version of the header. -1 if no version was/should be specified

    call write_binary_header(iunit, sheader, iflavor, kp%nspin, gvec%ng, syms%ntran, syms%cell_symmetry, crys%nat, &
         kp%nrk, kp%mnband, kp%ngkmax, gvec%ecutrho, kp%ecutwfc, gvec%FFTgrid, kp%kgrid, kp%shift, crys%celvol, &
         crys%alat, crys%avec, crys%adot, crys%recvol, crys%blat, crys%bvec, crys%bdot, syms%mtrx, syms%tnp, &
         crys%atyp, crys%apos, kp%ngk, kp%w, kp%rk, kp%ifmin, kp%ifmax, kp%el, kp%occ, nspinor = kp%nspinor, &
         warn = warn, dont_warn_kgrid = dont_warn_kgrid, sdate = sdate, stime = stime, version = version)

    return
  end subroutine write_binary_header_type
  !TEMP_SCALAR
  !=========================================================================
  !> This routine is turned into 12 variants by preprocessing.
  !! If nrecord and ng_record are not provided for read, the routine simply does not tell you what it found.
  !! For write, nrecord is given the default 1 if not provided. ng_record defaults to ng if nrecord = 1, but
  !! if nrecord > 1, ng_record must be provided to specify how to divide up the g-vectors.
  subroutine write_binary_gvectors(iunit, ng, ng_bound, gvec, nrecord, ng_record, bcast, gindex, dont_read)
    integer, intent(in) :: iunit
    integer, intent(in) :: ng !< used size of array
    integer, intent(in) :: ng_bound !< actual size of array, >= ng
    integer, intent(in) :: gvec(:, :) !< (3, ng_bound)
    integer, optional, intent(in) :: nrecord !< data/gvectors will be distributed among this many records
    integer, optional, pointer, intent(in) :: ng_record(:) !< number of gvectors in each record
    !! must be provided for write if nrecord > 1
    logical, optional, intent(in) :: bcast !< whether to do MPI_Bcast of what is read
    integer, optional, intent(in) :: gindex(:) !< map of order in file to order in gvec
    logical, optional, intent(in) :: dont_read !< if true, records will just be skipped; only for unformatted
    integer :: ig, igg, irecord, ii, nrecord_internal, ng_irecord
    character*100 :: tmpstr
    logical :: bcast_, dont_read_

    if(ng_bound < ng) then
       write(tmpstr,'(3a,i10,a,i10,a)') "In routine ", "write_binary_gvectors", ", ng_bound = ", ng_bound, " < ", ng, " = ng"
       call die(tmpstr)
    endif
    dont_read_ = .false.
    if(present(dont_read)) then
       call die("Write routine " + "write_binary_gvectors" + " cannot take argument dont_read.")
    endif
    if(peinf%inode .eq. 0) then
       if(present(gindex) .and. .not. dont_read_) then
          if(ubound(gindex, 1) < ng) &
               call die("In routine " + "write_binary_gvectors" + ", gindex array is too small")
          if(any(gindex(1:ng) > ng_bound) .or. any(gindex(1:ng) <= 0)) then
             call die("gindex out of bounds in " + "write_binary_gvectors")
          endif
       endif
    endif
    bcast_ = .not. dont_read_
    ! if not read, no point in broadcasting
    if(present(bcast)) then
       call die("Write routine " + "write_binary_gvectors" + " cannot take argument bcast.")
    endif
    ! if not reading, the size of the array passed is irrelevant
    if(.not. dont_read_) then
       if(ubound(gvec, 1) /= 3) &
            call die("In routine " + "write_binary_gvectors" + ", mismatch of dimension 1 for G-vectors array")
       if(ubound(gvec, 2) < ng_bound) &
            call die("In routine " + "write_binary_gvectors" + ", ng_bound is larger than dimension 2 for G-vectors array")
    endif
    if(.not. present(nrecord)) then
       nrecord_internal = 1
       ng_irecord = ng
    else
       nrecord_internal = nrecord
       !> check validity of information if going to write
       if(nrecord_internal .ne. 1 .and. .not. present(ng_record)) then
          call die("Routine " + "write_binary_gvectors" + " requires ng_record array if nrecord > 1.")
       endif
       if(present(ng_record)) then
          if(nrecord_internal .ne. sum(ng_record(1:nrecord))) then
             write(tmpstr,'(a, i10, a, i10, a, a)') "sum(ng_record) = ", sum(ng_record), " ! = ", nrecord , &
                  " = nrecord in arguments to routine ", "write_binary_gvectors"
             call die(tmpstr)
          endif
       endif
    endif
    if(peinf%inode .eq. 0) then
       write(iunit ) nrecord_internal
       if(nrecord_internal < 0) then
          write(tmpstr,'(3a, i10)') "In routine ", "write_binary_gvectors", " illegal nrecord ", nrecord_internal
          call die(tmpstr)
       endif
    endif
    if (peinf%inode.eq.0) call timacc(81,1)
    ! FHJ: zero buffers
    if(peinf%inode .eq. 0) then
       ig = 1
       !This is a PGI pragma to force the optimization level of this routine to -O1.
       !With -O2 or higher (including -fast) the line after the pragma causes a segmentation fault.
       !pgi$r opt=1
       do irecord = 1, nrecord_internal
          if(present(ng_record)) ng_irecord = ng_record(irecord)
          write(iunit ) ng_irecord
          if(dont_read_) then
             write(iunit )
          else
             if(present(gindex)) then
                write(iunit ) ((gvec(ii, gindex(igg)), ii = 1, 3), igg = ig, ig + ng_irecord - 1)
             else
                write(iunit ) ((gvec(ii, igg), ii = 1, 3), igg = ig, ig + ng_irecord - 1)
             endif
          endif
          ig = ig + ng_irecord
       enddo
    endif
    if (peinf%inode.eq.0) call timacc(81,2)
    if (peinf%inode.eq.0) call timacc(82,1)
    if (peinf%inode.eq.0) call timacc(82,2)

    return
  end subroutine write_binary_gvectors
  ! defined FORMATTED || defined
  ! these undefs prevent lots of cpp warnings
  !=========================================================================
  !
  ! Included from file wfn_rho_vxc_io.f90.
  ! You are not expected to understand this. --DAS
  !
  !=========================================================================
  ! defined FORMATTED || defined BINARY
  !> wrappers allowing selection of formatted/binary by argument
  !=========================================================================
  subroutine read_header(iunit, is_formatted, sheader, iflavor, ns, ng, ntran, cell_symmetry, &
       nat, nk, nbands, ngkmax, ecutrho, ecutwfc, FFTgrid, kgrid, kshift, celvol, alat, avec, adot, &
       recvol, blat, bvec, bdot, mtrx, tnp, atyp, apos, ngk, kw, kpt, ifmin, ifmax, energies, occupations, &
       nspinor, warn, sdate, stime, version)
    integer, intent(in) :: iunit !< unit number
    logical, intent(in) :: is_formatted !< true = formatted, false = binary
    character(len=3), intent(inout) :: sheader !< file header 'WFN'/'RHO'/'VXC'/'GET'
    integer, intent(inout) :: iflavor !< define type. always must be initialized. modified only if -1 on input
    integer, intent(out) :: ns, ng, ntran, cell_symmetry, nat, nk, nbands, ngkmax !< numbers of spins, G-vectors, symmetries,
    real(DP), intent(out) :: ecutrho, ecutwfc !< charge-density and wave-function cutoffs, in Ry
    integer, intent(out) :: FFTgrid(3), kgrid(3) !< FFT grid size, k-grid size
    real(DP), intent(out) :: kshift(3) !< k-grid offset
    real(DP), intent(out) :: celvol, alat, avec(3, 3), adot(3, 3) !< cell volume, lattice constant,
    real(DP), intent(out) :: recvol, blat, bvec(3, 3), bdot(3, 3) !< cell volume, lattice constant,
    integer, intent(out) :: mtrx(3, 3, 48) !< symmetry matrix
    real(DP), intent(out) :: tnp(3, 48) !< fractional translation
    integer, pointer, intent(out) :: atyp(:) !< atomic species
    real(DP), pointer, intent(out) :: apos(:,:) !< atomic positions (in units of alat)
    integer, pointer, intent(out) :: ngk(:) !< number of G-vectors for each k-pt, ngk(nk)
    real(DP), pointer, intent(out) :: kw(:), kpt(:, :) !< k-weight, kw(nk); k-coord, kpt(3, nk) in crystal coords
    integer, pointer, intent(out) :: ifmin(:, :), ifmax(:, :) !< lowest and highest occupied band, ifmin/max(nk, ns)
    real(DP), pointer, intent(out) :: energies(:, :, :) !< energies(nbands, nk, ns) in Ry
    real(DP), pointer, intent(out) :: occupations(:, :, :) !< occupations(nbands, nk, ns) between 0 and 1
    integer, intent(out) :: nspinor !< 2 if doing a two-component spinor calculation; 1 if not present
    logical, optional, intent(in) :: warn !< if false, suppresses warnings about inversion symmetries etc.
    character(len=32), optional, intent(out) :: sdate, stime !< if read, result from file is returned; if write, current is returned
    integer, optional, intent(out) :: version !< version of the header. -1 if no version was/should be specified

    if(is_formatted) then
       call read_format_header(iunit, sheader, iflavor, ns, ng, ntran, cell_symmetry, nat, nk, nbands, ngkmax, &
            ecutrho, ecutwfc, FFTgrid, kgrid, kshift, celvol, alat, avec, adot, recvol, blat, bvec, bdot, mtrx, tnp, &
            atyp, apos, ngk, kw, kpt, ifmin, ifmax, energies, occupations, &
            nspinor = nspinor, warn = warn, sdate = sdate, stime = stime, version=version)
    else
       call read_binary_header(iunit, sheader, iflavor, ns, ng, ntran, cell_symmetry, nat, nk, nbands, ngkmax, &
            ecutrho, ecutwfc, FFTgrid, kgrid, kshift, celvol, alat, avec, adot, recvol, blat, bvec, bdot, mtrx, tnp, &
            atyp, apos, ngk, kw, kpt, ifmin, ifmax, energies, occupations, &
            nspinor = nspinor, warn = warn, sdate = sdate, stime = stime, version=version)
    endif

    return
  end subroutine read_header
  !=========================================================================
  !> wrapper routine that uses typedefs types
  subroutine read_header_type(iunit, is_formatted, sheader, iflavor, kp, gvec, syms, crys, &
       warn, dont_warn_kgrid, sdate, stime, version)
    integer, intent(in) :: iunit
    logical, intent(in) :: is_formatted !< true = formatted, false = binary
    character(len=3), intent(inout) :: sheader
    integer, intent(inout) :: iflavor
    type(kpoints), intent(out) :: kp
    type(gspace), intent(out) :: gvec
    type(symmetry), intent(out) :: syms
    type(crystal), intent(out) :: crys
    logical, optional, intent(in) :: warn !< if false, suppresses warnings about inversion symmetries etc.
    logical, optional, intent(in) :: dont_warn_kgrid !< if true, validity of kgrid will not be checked.
    !! use for inteqp to allow interpolation onto non-uniform fine grids
    character(len=32), optional, intent(out) :: sdate, stime !< if read, result from file is returned; if write, current is returned
    integer, optional, intent(out) :: version !< version of the header. -1 if no version was/should be specified

    if(is_formatted) then
       call read_format_header_type(iunit, sheader, iflavor, kp, gvec, syms, crys, warn = warn, &
            dont_warn_kgrid = dont_warn_kgrid, sdate = sdate, stime = stime, version=version)
    else
       call read_binary_header_type(iunit, sheader, iflavor, kp, gvec, syms, crys, warn = warn, &
            dont_warn_kgrid = dont_warn_kgrid, sdate = sdate, stime = stime, version=version)
    endif

    return
  end subroutine read_header_type
  ! defined TEMP_SCALAR
  !=========================================================================
  subroutine read_gvectors(iunit, is_formatted, ng, ng_bound, gvec, nrecord, ng_record, bcast, gindex, dont_read)
    integer, intent(in) :: iunit
    logical, intent(in) :: is_formatted !< true = formatted, false = binary
    integer, intent(in) :: ng !< used size of array
    integer, intent(in) :: ng_bound !< actual size of array, >= ng
    integer, intent(out) :: gvec(:, :) !< (3, ng_bound)
    integer, optional, intent(out) :: nrecord !< data/gvectors will be distributed among this many records
    integer, optional, pointer, intent(out) :: ng_record(:) !< number of gvectors in each record
    logical, optional, intent(in) :: bcast !< whether to do MPI_Bcast of what is read
    integer, optional, intent(in) :: gindex(:) !< map of order in file to order in gvec
    logical, optional, intent(in) :: dont_read !< if true, records will just be skipped; only for unformatted

    ! note: surprisingly, the Fortran standard permits passing optional arguments to other routines
    ! in which those arguments are also optional.
    if(is_formatted) then
       call read_format_gvectors(iunit, ng, ng_bound, gvec, nrecord = nrecord, &
            ng_record = ng_record, bcast = bcast, gindex = gindex, dont_read = dont_read)
    else
       call read_binary_gvectors(iunit, ng, ng_bound, gvec, nrecord = nrecord, &
            ng_record = ng_record, bcast = bcast, gindex = gindex, dont_read = dont_read)
    endif

    return
  end subroutine read_gvectors
  ! defined FORMATTED || defined BINARY
  ! these undefs prevent lots of cpp warnings
  !> write
  !=========================================================================
  !
  ! Included from file wfn_rho_vxc_io.f90.
  ! You are not expected to understand this. --DAS
  !
  !=========================================================================
  ! defined FORMATTED || defined BINARY
  !> wrappers allowing selection of formatted/binary by argument
  !=========================================================================
  subroutine write_header(iunit, is_formatted, sheader, iflavor, ns, ng, ntran, cell_symmetry, &
       nat, nk, nbands, ngkmax, ecutrho, ecutwfc, FFTgrid, kgrid, kshift, celvol, alat, avec, adot, &
       recvol, blat, bvec, bdot, mtrx, tnp, atyp, apos, ngk, kw, kpt, ifmin, ifmax, energies, occupations, &
       nspinor, warn, sdate, stime, version)
    integer, intent(in) :: iunit !< unit number
    logical, intent(in) :: is_formatted !< true = formatted, false = binary
    character(len=3), intent(inout) :: sheader !< file header 'WFN'/'RHO'/'VXC'/'GET'
    integer, intent(in) :: iflavor !< define type. always must be initialized. modified only if -1 on input
    integer, intent(in) :: ns, ng, ntran, cell_symmetry, nat, nk, nbands, ngkmax !< numbers of spins, G-vectors, symmetries,
    real(DP), intent(in) :: ecutrho, ecutwfc !< charge-density and wave-function cutoffs, in Ry
    integer, intent(in) :: FFTgrid(3), kgrid(3) !< FFT grid size, k-grid size
    real(DP), intent(in) :: kshift(3) !< k-grid offset
    real(DP), intent(in) :: celvol, alat, avec(3, 3), adot(3, 3) !< cell volume, lattice constant,
    real(DP), intent(in) :: recvol, blat, bvec(3, 3), bdot(3, 3) !< cell volume, lattice constant,
    integer, intent(in) :: mtrx(3, 3, 48) !< symmetry matrix
    real(DP), intent(in) :: tnp(3, 48) !< fractional translation
    integer, pointer, intent(in) :: atyp(:) !< atomic species
    real(DP), pointer, intent(in) :: apos(:,:) !< atomic positions (in units of alat)
    integer, pointer, intent(in) :: ngk(:) !< number of G-vectors for each k-pt, ngk(nk)
    real(DP), pointer, intent(in) :: kw(:), kpt(:, :) !< k-weight, kw(nk); k-coord, kpt(3, nk) in crystal coords
    integer, pointer, intent(in) :: ifmin(:, :), ifmax(:, :) !< lowest and highest occupied band, ifmin/max(nk, ns)
    real(DP), pointer, intent(in) :: energies(:, :, :) !< energies(nbands, nk, ns) in Ry
    real(DP), pointer, intent(in) :: occupations(:, :, :) !< occupations(nbands, nk, ns) between 0 and 1
    integer, intent(in) :: nspinor !< 2 if doing a two-component spinor calculation; 1 if not present
    logical, optional, intent(in) :: warn !< if false, suppresses warnings about inversion symmetries etc.
    character(len=32), optional, intent(out) :: sdate, stime !< if read, result from file is returned; if write, current is returned
    integer, optional, intent(in) :: version !< version of the header. -1 if no version was/should be specified

    if(is_formatted) then
       call write_format_header(iunit, sheader, iflavor, ns, ng, ntran, cell_symmetry, nat, nk, nbands, ngkmax, &
            ecutrho, ecutwfc, FFTgrid, kgrid, kshift, celvol, alat, avec, adot, recvol, blat, bvec, bdot, mtrx, tnp, &
            atyp, apos, ngk, kw, kpt, ifmin, ifmax, energies, occupations, &
            nspinor = nspinor, warn = warn, sdate = sdate, stime = stime, version=version)
    else
       call write_binary_header(iunit, sheader, iflavor, ns, ng, ntran, cell_symmetry, nat, nk, nbands, ngkmax, &
            ecutrho, ecutwfc, FFTgrid, kgrid, kshift, celvol, alat, avec, adot, recvol, blat, bvec, bdot, mtrx, tnp, &
            atyp, apos, ngk, kw, kpt, ifmin, ifmax, energies, occupations, &
            nspinor = nspinor, warn = warn, sdate = sdate, stime = stime, version=version)
    endif

    return
  end subroutine write_header
  !=========================================================================
  !> wrapper routine that uses typedefs types
  subroutine write_header_type(iunit, is_formatted, sheader, iflavor, kp, gvec, syms, crys, &
       warn, dont_warn_kgrid, sdate, stime, version)
    integer, intent(in) :: iunit
    logical, intent(in) :: is_formatted !< true = formatted, false = binary
    character(len=3), intent(inout) :: sheader
    integer, intent(in) :: iflavor
    type(kpoints), intent(in) :: kp
    type(gspace), intent(in) :: gvec
    type(symmetry), intent(in) :: syms
    type(crystal), intent(in) :: crys
    logical, optional, intent(in) :: warn !< if false, suppresses warnings about inversion symmetries etc.
    logical, optional, intent(in) :: dont_warn_kgrid !< if true, validity of kgrid will not be checked.
    !! use for inteqp to allow interpolation onto non-uniform fine grids
    character(len=32), optional, intent(out) :: sdate, stime !< if read, result from file is returned; if write, current is returned
    integer, optional, intent(in) :: version !< version of the header. -1 if no version was/should be specified

    if(is_formatted) then
       call write_format_header_type(iunit, sheader, iflavor, kp, gvec, syms, crys, warn = warn, &
            dont_warn_kgrid = dont_warn_kgrid, sdate = sdate, stime = stime, version=version)
    else
       call write_binary_header_type(iunit, sheader, iflavor, kp, gvec, syms, crys, warn = warn, &
            dont_warn_kgrid = dont_warn_kgrid, sdate = sdate, stime = stime, version=version)
    endif

    return
  end subroutine write_header_type
  ! defined TEMP_SCALAR
  !=========================================================================
  subroutine write_gvectors(iunit, is_formatted, ng, ng_bound, gvec, nrecord, ng_record, bcast, gindex, dont_read)
    integer, intent(in) :: iunit
    logical, intent(in) :: is_formatted !< true = formatted, false = binary
    integer, intent(in) :: ng !< used size of array
    integer, intent(in) :: ng_bound !< actual size of array, >= ng
    integer, intent(in) :: gvec(:, :) !< (3, ng_bound)
    integer, optional, intent(in) :: nrecord !< data/gvectors will be distributed among this many records
    integer, optional, pointer, intent(in) :: ng_record(:) !< number of gvectors in each record
    logical, optional, intent(in) :: bcast !< whether to do MPI_Bcast of what is read
    integer, optional, intent(in) :: gindex(:) !< map of order in file to order in gvec
    logical, optional, intent(in) :: dont_read !< if true, records will just be skipped; only for unformatted

    ! note: surprisingly, the Fortran standard permits passing optional arguments to other routines
    ! in which those arguments are also optional.
    if(is_formatted) then
       call write_format_gvectors(iunit, ng, ng_bound, gvec, nrecord = nrecord, &
            ng_record = ng_record, bcast = bcast, gindex = gindex, dont_read = dont_read)
    else
       call write_binary_gvectors(iunit, ng, ng_bound, gvec, nrecord = nrecord, &
            ng_record = ng_record, bcast = bcast, gindex = gindex, dont_read = dont_read)
    endif

    return
  end subroutine write_gvectors
  ! defined FORMATTED || defined BINARY
  ! these undefs prevent lots of cpp warnings
  !=========================================================================
  !
  ! Included from file wfn_rho_vxc_io.f90.
  ! You are not expected to understand this. --DAS
  !
  !=========================================================================
  !complex(DPC)
  !=========================================================================
  !> This routine is turned into 12 variants by preprocessing.
  !! If nrecord and ng_record are not provided for read, the routine simply does not tell you what it found.
  !! For write, nrecord is given the default 1 if not provided. ng_record defaults to ng if nrecord = 1, but
  !! if nrecord > 1, ng_record must be provided to specify how to divide up the g-vectors.
  subroutine read_format_complex_data(iunit, ng, ng_bound, ns, data, nrecord, ng_record, bcast, gindex, dont_read)
    integer, intent(in) :: iunit
    integer, intent(in) :: ng !< used size of array
    integer, intent(in) :: ng_bound !< actual size of array, >= ng
    integer, intent(in) :: ns
    complex(DPC), intent(out) :: data(:, :) !< (ng_bound, ns)
    integer, optional, intent(out) :: nrecord !< data/gvectors will be distributed among this many records
    integer, optional, pointer, intent(out) :: ng_record(:) !< number of gvectors in each record
    !! must be provided for write if nrecord > 1
    logical, optional, intent(in) :: bcast !< whether to do MPI_Bcast of what is read
    integer, optional, intent(in) :: gindex(:) !< map of order in file to order in gvec
    logical, optional, intent(in) :: dont_read !< if true, records will just be skipped; only for unformatted
    integer :: ig, igg, irecord, ii, nrecord_internal, ng_irecord
    character*100 :: tmpstr
    logical :: bcast_, dont_read_

    if(ng_bound < ng) then
       write(tmpstr,'(3a,i10,a,i10,a)') "In routine ", "read_format_complex_data", ", ng_bound = ", ng_bound, " < ", ng, " = ng"
       call die(tmpstr)
    endif
    dont_read_ = .false.
    if(present(dont_read)) then
       if(dont_read) call die("Formatted routine " + "read_format_complex_data" + " cannot take argument dont_read = .true.")
       dont_read_ = dont_read
    endif
    if(peinf%inode .eq. 0) then
       if(present(gindex) .and. .not. dont_read_) then
          if(ubound(gindex, 1) < ng) &
               call die("In routine " + "read_format_complex_data" + ", gindex array is too small")
          if(any(gindex(1:ng) > ng_bound) .or. any(gindex(1:ng) <= 0)) then
             call die("gindex out of bounds in " + "read_format_complex_data")
          endif
       endif
    endif
    bcast_ = .not. dont_read_
    ! if not read, no point in broadcasting
    if(present(bcast)) then
       bcast_ = bcast
    endif
    ! if not reading, the size of the array passed is irrelevant
    if(.not. dont_read_) then
       if(ubound(data, 1) /= ng_bound) &
            call die("In routine " + "read_format_complex_data" + ", mismatch of dimension 1 for data array")
       if(ubound(data, 2) /= ns) &
            call die("In routine " + "read_format_complex_data" + ", mismatch of dimension 2 for data array")
    endif
    if(peinf%inode .eq. 0) then
       read(iunit , *) nrecord_internal
       if(nrecord_internal < 0) then
          write(tmpstr,'(3a, i10)') "In routine ", "read_format_complex_data", " illegal nrecord ", nrecord_internal
          call die(tmpstr)
       endif
    endif
    if(present(nrecord)) then
       nrecord = nrecord_internal
       if(peinf%npes > 1) call MPI_BCAST(nrecord, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
    endif
    if(present(ng_record)) then
       allocate(ng_record (nrecord_internal))
    endif
    if (peinf%inode.eq.0) call timacc(81,1)
    ! FHJ: zero buffers
    data(:,:) = (0.0d0,0.0d0)
    if(peinf%inode .eq. 0) then
       ig = 1
       !This is a PGI pragma to force the optimization level of this routine to -O1.
       !With -O2 or higher (including -fast) the line after the pragma causes a segmentation fault.
       !pgi$r opt=1
       do irecord = 1, nrecord_internal
          read(iunit , *) ng_irecord
          if(present(ng_record)) ng_record(irecord) = ng_irecord
          if(dont_read_) then
             read(iunit , *)
          else
             if(present(gindex)) then
                read(iunit , *) ((data(gindex(igg), ii), igg = ig, ig + ng_irecord - 1), ii = 1, ns)
             else
                read(iunit , *) ((data(igg, ii), igg = ig, ig + ng_irecord - 1), ii = 1, ns)
             endif
          endif
          ig = ig + ng_irecord
       enddo
    endif
    if (peinf%inode.eq.0) call timacc(81,2)
    if (peinf%inode.eq.0) call timacc(82,1)
    if(peinf%npes > 1) then
       if(present(ng_record)) call MPI_BCAST(ng_record, nrecord_internal, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
       if(bcast_) then
          call MPI_BCAST(data, ng_bound * ns, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, mpierr)
       endif
    endif
    if (peinf%inode.eq.0) call timacc(82,2)

    return
  end subroutine read_format_complex_data
  ! defined || defined BINARY
  ! these undefs prevent lots of cpp warnings
  !> write
  !=========================================================================
  !
  ! Included from file wfn_rho_vxc_io.f90.
  ! You are not expected to understand this. --DAS
  !
  !=========================================================================
  !complex(DPC)
  !=========================================================================
  !> This routine is turned into 12 variants by preprocessing.
  !! If nrecord and ng_record are not provided for read, the routine simply does not tell you what it found.
  !! For write, nrecord is given the default 1 if not provided. ng_record defaults to ng if nrecord = 1, but
  !! if nrecord > 1, ng_record must be provided to specify how to divide up the g-vectors.
  subroutine write_format_complex_data(iunit, ng, ng_bound, ns, data, nrecord, ng_record, bcast, gindex, dont_read)
    integer, intent(in) :: iunit
    integer, intent(in) :: ng !< used size of array
    integer, intent(in) :: ng_bound !< actual size of array, >= ng
    integer, intent(in) :: ns
    complex(DPC), intent(in) :: data(:, :) !< (ng_bound, ns)
    integer, optional, intent(in) :: nrecord !< data/gvectors will be distributed among this many records
    integer, optional, pointer, intent(in) :: ng_record(:) !< number of gvectors in each record
    !! must be provided for write if nrecord > 1
    logical, optional, intent(in) :: bcast !< whether to do MPI_Bcast of what is read
    integer, optional, intent(in) :: gindex(:) !< map of order in file to order in gvec
    logical, optional, intent(in) :: dont_read !< if true, records will just be skipped; only for unformatted
    integer :: ig, igg, irecord, ii, nrecord_internal, ng_irecord
    character*100 :: tmpstr
    logical :: bcast_, dont_read_

    if(ng_bound < ng) then
       write(tmpstr,'(3a,i10,a,i10,a)') "In routine ", "write_format_complex_data", ", ng_bound = ", ng_bound, " < ", ng, " = ng"
       call die(tmpstr)
    endif
    dont_read_ = .false.
    if(present(dont_read)) then
       call die("Write routine " + "write_format_complex_data" + " cannot take argument dont_read.")
    endif
    if(peinf%inode .eq. 0) then
       if(present(gindex) .and. .not. dont_read_) then
          if(ubound(gindex, 1) < ng) &
               call die("In routine " + "write_format_complex_data" + ", gindex array is too small")
          if(any(gindex(1:ng) > ng_bound) .or. any(gindex(1:ng) <= 0)) then
             call die("gindex out of bounds in " + "write_format_complex_data")
          endif
       endif
    endif
    bcast_ = .not. dont_read_
    ! if not read, no point in broadcasting
    if(present(bcast)) then
       call die("Write routine " + "write_format_complex_data" + " cannot take argument bcast.")
    endif
    ! if not reading, the size of the array passed is irrelevant
    if(.not. dont_read_) then
       if(ubound(data, 1) /= ng_bound) &
            call die("In routine " + "write_format_complex_data" + ", mismatch of dimension 1 for data array")
       if(ubound(data, 2) /= ns) &
            call die("In routine " + "write_format_complex_data" + ", mismatch of dimension 2 for data array")
    endif
    if(.not. present(nrecord)) then
       nrecord_internal = 1
       ng_irecord = ng
    else
       nrecord_internal = nrecord
       !> check validity of information if going to write
       if(nrecord_internal .ne. 1 .and. .not. present(ng_record)) then
          call die("Routine " + "write_format_complex_data" + " requires ng_record array if nrecord > 1.")
       endif
       if(present(ng_record)) then
          if(nrecord_internal .ne. sum(ng_record(1:nrecord))) then
             write(tmpstr,'(a, i10, a, i10, a, a)') "sum(ng_record) = ", sum(ng_record), " ! = ", nrecord , &
                  " = nrecord in arguments to routine ", "write_format_complex_data"
             call die(tmpstr)
          endif
       endif
    endif
    if(peinf%inode .eq. 0) then
       write(iunit , *) nrecord_internal
       if(nrecord_internal < 0) then
          write(tmpstr,'(3a, i10)') "In routine ", "write_format_complex_data", " illegal nrecord ", nrecord_internal
          call die(tmpstr)
       endif
    endif
    if (peinf%inode.eq.0) call timacc(81,1)
    ! FHJ: zero buffers
    if(peinf%inode .eq. 0) then
       ig = 1
       !This is a PGI pragma to force the optimization level of this routine to -O1.
       !With -O2 or higher (including -fast) the line after the pragma causes a segmentation fault.
       !pgi$r opt=1
       do irecord = 1, nrecord_internal
          if(present(ng_record)) ng_irecord = ng_record(irecord)
          write(iunit , *) ng_irecord
          if(dont_read_) then
             write(iunit , *)
          else
             if(present(gindex)) then
                write(iunit , *) ((data(gindex(igg), ii), igg = ig, ig + ng_irecord - 1), ii = 1, ns)
             else
                write(iunit , *) ((data(igg, ii), igg = ig, ig + ng_irecord - 1), ii = 1, ns)
             endif
          endif
          ig = ig + ng_irecord
       enddo
    endif
    if (peinf%inode.eq.0) call timacc(81,2)
    if (peinf%inode.eq.0) call timacc(82,1)
    if (peinf%inode.eq.0) call timacc(82,2)

    return
  end subroutine write_format_complex_data
  ! defined || defined BINARY
  ! these undefs prevent lots of cpp warnings
  !=========================================================================
  !
  ! Included from file wfn_rho_vxc_io.f90.
  ! You are not expected to understand this. --DAS
  !
  !=========================================================================
  !complex(DPC)
  !=========================================================================
  !> This routine is turned into 12 variants by preprocessing.
  !! If nrecord and ng_record are not provided for read, the routine simply does not tell you what it found.
  !! For write, nrecord is given the default 1 if not provided. ng_record defaults to ng if nrecord = 1, but
  !! if nrecord > 1, ng_record must be provided to specify how to divide up the g-vectors.
  subroutine read_binary_complex_data(iunit, ng, ng_bound, ns, data, nrecord, ng_record, bcast, gindex, dont_read)
    integer, intent(in) :: iunit
    integer, intent(in) :: ng !< used size of array
    integer, intent(in) :: ng_bound !< actual size of array, >= ng
    integer, intent(in) :: ns
    complex(DPC), intent(out) :: data(:, :) !< (ng_bound, ns)
    integer, optional, intent(out) :: nrecord !< data/gvectors will be distributed among this many records
    integer, optional, pointer, intent(out) :: ng_record(:) !< number of gvectors in each record
    !! must be provided for write if nrecord > 1
    logical, optional, intent(in) :: bcast !< whether to do MPI_Bcast of what is read
    integer, optional, intent(in) :: gindex(:) !< map of order in file to order in gvec
    logical, optional, intent(in) :: dont_read !< if true, records will just be skipped; only for unformatted
    integer :: ig, igg, irecord, ii, nrecord_internal, ng_irecord
    character*100 :: tmpstr
    logical :: bcast_, dont_read_

    if(ng_bound < ng) then
       write(tmpstr,'(3a,i10,a,i10,a)') "In routine ", "read_binary_complex_data", ", ng_bound = ", ng_bound, " < ", ng, " = ng"
       call die(tmpstr)
    endif
    dont_read_ = .false.
    if(present(dont_read)) then
       dont_read_ = dont_read
    endif
    if(peinf%inode .eq. 0) then
       if(present(gindex) .and. .not. dont_read_) then
          if(ubound(gindex, 1) < ng) &
               call die("In routine " + "read_binary_complex_data" + ", gindex array is too small")
          if(any(gindex(1:ng) > ng_bound) .or. any(gindex(1:ng) <= 0)) then
             call die("gindex out of bounds in " + "read_binary_complex_data")
          endif
       endif
    endif
    bcast_ = .not. dont_read_
    ! if not read, no point in broadcasting
    if(present(bcast)) then
       bcast_ = bcast
    endif
    ! if not reading, the size of the array passed is irrelevant
    if(.not. dont_read_) then
       if(ubound(data, 1) /= ng_bound) &
            call die("In routine " + "read_binary_complex_data" + ", mismatch of dimension 1 for data array")
       if(ubound(data, 2) /= ns) &
            call die("In routine " + "read_binary_complex_data" + ", mismatch of dimension 2 for data array")
    endif
    if(peinf%inode .eq. 0) then
       read(iunit ) nrecord_internal
       if(nrecord_internal < 0) then
          write(tmpstr,'(3a, i10)') "In routine ", "read_binary_complex_data", " illegal nrecord ", nrecord_internal
          call die(tmpstr)
       endif
    endif
    if(present(nrecord)) then
       nrecord = nrecord_internal
       if(peinf%npes > 1) call MPI_BCAST(nrecord, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
    endif
    if(present(ng_record)) then
       allocate(ng_record (nrecord_internal))
    endif
    if (peinf%inode.eq.0) call timacc(81,1)
    ! FHJ: zero buffers
    data(:,:) = (0.0d0,0.0d0)
    if(peinf%inode .eq. 0) then
       ig = 1
       !This is a PGI pragma to force the optimization level of this routine to -O1.
       !With -O2 or higher (including -fast) the line after the pragma causes a segmentation fault.
       !pgi$r opt=1
       do irecord = 1, nrecord_internal
          read(iunit ) ng_irecord
          if(present(ng_record)) ng_record(irecord) = ng_irecord
          if(dont_read_) then
             read(iunit )
          else
             if(present(gindex)) then
                read(iunit ) ((data(gindex(igg), ii), igg = ig, ig + ng_irecord - 1), ii = 1, ns)
             else
                read(iunit ) ((data(igg, ii), igg = ig, ig + ng_irecord - 1), ii = 1, ns)
             endif
          endif
          ig = ig + ng_irecord
       enddo
    endif
    if (peinf%inode.eq.0) call timacc(81,2)
    if (peinf%inode.eq.0) call timacc(82,1)
    if(peinf%npes > 1) then
       if(present(ng_record)) call MPI_BCAST(ng_record, nrecord_internal, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
       if(bcast_) then
          call MPI_BCAST(data, ng_bound * ns, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, mpierr)
       endif
    endif
    if (peinf%inode.eq.0) call timacc(82,2)

    return
  end subroutine read_binary_complex_data
  ! defined FORMATTED || defined
  ! these undefs prevent lots of cpp warnings
  !> write
  !=========================================================================
  !
  ! Included from file wfn_rho_vxc_io.f90.
  ! You are not expected to understand this. --DAS
  !
  !=========================================================================
  !complex(DPC)
  !=========================================================================
  !> This routine is turned into 12 variants by preprocessing.
  !! If nrecord and ng_record are not provided for read, the routine simply does not tell you what it found.
  !! For write, nrecord is given the default 1 if not provided. ng_record defaults to ng if nrecord = 1, but
  !! if nrecord > 1, ng_record must be provided to specify how to divide up the g-vectors.
  subroutine write_binary_complex_data(iunit, ng, ng_bound, ns, data, nrecord, ng_record, bcast, gindex, dont_read)
    integer, intent(in) :: iunit
    integer, intent(in) :: ng !< used size of array
    integer, intent(in) :: ng_bound !< actual size of array, >= ng
    integer, intent(in) :: ns
    complex(DPC), intent(in) :: data(:, :) !< (ng_bound, ns)
    integer, optional, intent(in) :: nrecord !< data/gvectors will be distributed among this many records
    integer, optional, pointer, intent(in) :: ng_record(:) !< number of gvectors in each record
    !! must be provided for write if nrecord > 1
    logical, optional, intent(in) :: bcast !< whether to do MPI_Bcast of what is read
    integer, optional, intent(in) :: gindex(:) !< map of order in file to order in gvec
    logical, optional, intent(in) :: dont_read !< if true, records will just be skipped; only for unformatted
    integer :: ig, igg, irecord, ii, nrecord_internal, ng_irecord
    character*100 :: tmpstr
    logical :: bcast_, dont_read_

    if(ng_bound < ng) then
       write(tmpstr,'(3a,i10,a,i10,a)') "In routine ", "write_binary_complex_data", ", ng_bound = ", ng_bound, " < ", ng, " = ng"
       call die(tmpstr)
    endif
    dont_read_ = .false.
    if(present(dont_read)) then
       call die("Write routine " + "write_binary_complex_data" + " cannot take argument dont_read.")
    endif
    if(peinf%inode .eq. 0) then
       if(present(gindex) .and. .not. dont_read_) then
          if(ubound(gindex, 1) < ng) &
               call die("In routine " + "write_binary_complex_data" + ", gindex array is too small")
          if(any(gindex(1:ng) > ng_bound) .or. any(gindex(1:ng) <= 0)) then
             call die("gindex out of bounds in " + "write_binary_complex_data")
          endif
       endif
    endif
    bcast_ = .not. dont_read_
    ! if not read, no point in broadcasting
    if(present(bcast)) then
       call die("Write routine " + "write_binary_complex_data" + " cannot take argument bcast.")
    endif
    ! if not reading, the size of the array passed is irrelevant
    if(.not. dont_read_) then
       if(ubound(data, 1) /= ng_bound) &
            call die("In routine " + "write_binary_complex_data" + ", mismatch of dimension 1 for data array")
       if(ubound(data, 2) /= ns) &
            call die("In routine " + "write_binary_complex_data" + ", mismatch of dimension 2 for data array")
    endif
    if(.not. present(nrecord)) then
       nrecord_internal = 1
       ng_irecord = ng
    else
       nrecord_internal = nrecord
       !> check validity of information if going to write
       if(nrecord_internal .ne. 1 .and. .not. present(ng_record)) then
          call die("Routine " + "write_binary_complex_data" + " requires ng_record array if nrecord > 1.")
       endif
       if(present(ng_record)) then
          if(nrecord_internal .ne. sum(ng_record(1:nrecord))) then
             write(tmpstr,'(a, i10, a, i10, a, a)') "sum(ng_record) = ", sum(ng_record), " ! = ", nrecord , &
                  " = nrecord in arguments to routine ", "write_binary_complex_data"
             call die(tmpstr)
          endif
       endif
    endif
    if(peinf%inode .eq. 0) then
       write(iunit ) nrecord_internal
       if(nrecord_internal < 0) then
          write(tmpstr,'(3a, i10)') "In routine ", "write_binary_complex_data", " illegal nrecord ", nrecord_internal
          call die(tmpstr)
       endif
    endif
    if (peinf%inode.eq.0) call timacc(81,1)
    ! FHJ: zero buffers
    if(peinf%inode .eq. 0) then
       ig = 1
       !This is a PGI pragma to force the optimization level of this routine to -O1.
       !With -O2 or higher (including -fast) the line after the pragma causes a segmentation fault.
       !pgi$r opt=1
       do irecord = 1, nrecord_internal
          if(present(ng_record)) ng_irecord = ng_record(irecord)
          write(iunit ) ng_irecord
          if(dont_read_) then
             write(iunit )
          else
             if(present(gindex)) then
                write(iunit ) ((data(gindex(igg), ii), igg = ig, ig + ng_irecord - 1), ii = 1, ns)
             else
                write(iunit ) ((data(igg, ii), igg = ig, ig + ng_irecord - 1), ii = 1, ns)
             endif
          endif
          ig = ig + ng_irecord
       enddo
    endif
    if (peinf%inode.eq.0) call timacc(81,2)
    if (peinf%inode.eq.0) call timacc(82,1)
    if (peinf%inode.eq.0) call timacc(82,2)

    return
  end subroutine write_binary_complex_data
  ! defined FORMATTED || defined
  ! these undefs prevent lots of cpp warnings
  !=========================================================================
  !
  ! Included from file wfn_rho_vxc_io.f90.
  ! You are not expected to understand this. --DAS
  !
  !=========================================================================
  ! defined FORMATTED || defined BINARY
  ! defined complex(DPC)
  !=========================================================================
  subroutine read_complex_data(iunit, is_formatted, ng, ng_bound, ns, data, nrecord, ng_record, bcast, gindex, dont_read)
    integer, intent(in) :: iunit
    logical, intent(in) :: is_formatted !< true = formatted, false = binary
    integer, intent(in) :: ng !< used size of array
    integer, intent(in) :: ng_bound !< actual size of array, >= ng
    integer, intent(in) :: ns
    complex(DPC), intent(out) :: data(:, :) !< (ng_bound, ns)
    integer, optional, intent(out) :: nrecord !< data/gvectors will be distributed among this many records
    integer, optional, pointer, intent(out) :: ng_record(:) !< number of gvectors in each record
    logical, optional, intent(in) :: bcast !< whether to do MPI_Bcast of what is read
    integer, optional, intent(in) :: gindex(:) !< map of order in file to order in gvec
    logical, optional, intent(in) :: dont_read !< if true, records will just be skipped; only for unformatted

    ! note: surprisingly, the Fortran standard permits passing optional arguments to other routines
    ! in which those arguments are also optional.
    if(is_formatted) then
       call read_format_complex_data(iunit, ng, ng_bound, ns, data, nrecord = nrecord, &
            ng_record = ng_record, bcast = bcast, gindex = gindex, dont_read = dont_read)
    else
       call read_binary_complex_data(iunit, ng, ng_bound, ns, data, nrecord = nrecord, &
            ng_record = ng_record, bcast = bcast, gindex = gindex, dont_read = dont_read)
    endif

    return
  end subroutine read_complex_data
  ! defined FORMATTED || defined BINARY
  ! these undefs prevent lots of cpp warnings
  !> write
  !=========================================================================
  !
  ! Included from file wfn_rho_vxc_io.f90.
  ! You are not expected to understand this. --DAS
  !
  !=========================================================================
  ! defined FORMATTED || defined BINARY
  ! defined complex(DPC)
  !=========================================================================
  subroutine write_complex_data(iunit, is_formatted, ng, ng_bound, ns, data, nrecord, ng_record, bcast, gindex, dont_read)
    integer, intent(in) :: iunit
    logical, intent(in) :: is_formatted !< true = formatted, false = binary
    integer, intent(in) :: ng !< used size of array
    integer, intent(in) :: ng_bound !< actual size of array, >= ng
    integer, intent(in) :: ns
    complex(DPC), intent(in) :: data(:, :) !< (ng_bound, ns)
    integer, optional, intent(in) :: nrecord !< data/gvectors will be distributed among this many records
    integer, optional, pointer, intent(in) :: ng_record(:) !< number of gvectors in each record
    logical, optional, intent(in) :: bcast !< whether to do MPI_Bcast of what is read
    integer, optional, intent(in) :: gindex(:) !< map of order in file to order in gvec
    logical, optional, intent(in) :: dont_read !< if true, records will just be skipped; only for unformatted

    ! note: surprisingly, the Fortran standard permits passing optional arguments to other routines
    ! in which those arguments are also optional.
    if(is_formatted) then
       call write_format_complex_data(iunit, ng, ng_bound, ns, data, nrecord = nrecord, &
            ng_record = ng_record, bcast = bcast, gindex = gindex, dont_read = dont_read)
    else
       call write_binary_complex_data(iunit, ng, ng_bound, ns, data, nrecord = nrecord, &
            ng_record = ng_record, bcast = bcast, gindex = gindex, dont_read = dont_read)
    endif

    return
  end subroutine write_complex_data
  ! defined FORMATTED || defined BINARY
  ! these undefs prevent lots of cpp warnings
  !=========================================================================
  !
  ! Included from file wfn_rho_vxc_io.f90.
  ! You are not expected to understand this. --DAS
  !
  !=========================================================================
  !real(DP)
  !=========================================================================
  !> This routine is turned into 12 variants by preprocessing.
  !! If nrecord and ng_record are not provided for read, the routine simply does not tell you what it found.
  !! For write, nrecord is given the default 1 if not provided. ng_record defaults to ng if nrecord = 1, but
  !! if nrecord > 1, ng_record must be provided to specify how to divide up the g-vectors.
  subroutine read_format_real_data(iunit, ng, ng_bound, ns, data, nrecord, ng_record, bcast, gindex, dont_read)
    integer, intent(in) :: iunit
    integer, intent(in) :: ng !< used size of array
    integer, intent(in) :: ng_bound !< actual size of array, >= ng
    integer, intent(in) :: ns
    real(DP), intent(out) :: data(:, :) !< (ng_bound, ns)
    integer, optional, intent(out) :: nrecord !< data/gvectors will be distributed among this many records
    integer, optional, pointer, intent(out) :: ng_record(:) !< number of gvectors in each record
    !! must be provided for write if nrecord > 1
    logical, optional, intent(in) :: bcast !< whether to do MPI_Bcast of what is read
    integer, optional, intent(in) :: gindex(:) !< map of order in file to order in gvec
    logical, optional, intent(in) :: dont_read !< if true, records will just be skipped; only for unformatted
    integer :: ig, igg, irecord, ii, nrecord_internal, ng_irecord
    character*100 :: tmpstr
    logical :: bcast_, dont_read_

    if(ng_bound < ng) then
       write(tmpstr,'(3a,i10,a,i10,a)') "In routine ", "read_format_real_data", ", ng_bound = ", ng_bound, " < ", ng, " = ng"
       call die(tmpstr)
    endif
    dont_read_ = .false.
    if(present(dont_read)) then
       if(dont_read) call die("Formatted routine " + "read_format_real_data" + " cannot take argument dont_read = .true.")
       dont_read_ = dont_read
    endif
    if(peinf%inode .eq. 0) then
       if(present(gindex) .and. .not. dont_read_) then
          if(ubound(gindex, 1) < ng) &
               call die("In routine " + "read_format_real_data" + ", gindex array is too small")
          if(any(gindex(1:ng) > ng_bound) .or. any(gindex(1:ng) <= 0)) then
             call die("gindex out of bounds in " + "read_format_real_data")
          endif
       endif
    endif
    bcast_ = .not. dont_read_
    ! if not read, no point in broadcasting
    if(present(bcast)) then
       bcast_ = bcast
    endif
    ! if not reading, the size of the array passed is irrelevant
    if(.not. dont_read_) then
       if(ubound(data, 1) /= ng_bound) &
            call die("In routine " + "read_format_real_data" + ", mismatch of dimension 1 for data array")
       if(ubound(data, 2) /= ns) &
            call die("In routine " + "read_format_real_data" + ", mismatch of dimension 2 for data array")
    endif
    if(peinf%inode .eq. 0) then
       read(iunit , *) nrecord_internal
       if(nrecord_internal < 0) then
          write(tmpstr,'(3a, i10)') "In routine ", "read_format_real_data", " illegal nrecord ", nrecord_internal
          call die(tmpstr)
       endif
    endif
    if(present(nrecord)) then
       nrecord = nrecord_internal
       if(peinf%npes > 1) call MPI_BCAST(nrecord, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
    endif
    if(present(ng_record)) then
       allocate(ng_record (nrecord_internal))
    endif
    if (peinf%inode.eq.0) call timacc(81,1)
    ! FHJ: zero buffers
    data(:,:) = (0.0d0,0.0d0)
    if(peinf%inode .eq. 0) then
       ig = 1
       !This is a PGI pragma to force the optimization level of this routine to -O1.
       !With -O2 or higher (including -fast) the line after the pragma causes a segmentation fault.
       !pgi$r opt=1
       do irecord = 1, nrecord_internal
          read(iunit , *) ng_irecord
          if(present(ng_record)) ng_record(irecord) = ng_irecord
          if(dont_read_) then
             read(iunit , *)
          else
             if(present(gindex)) then
                read(iunit , *) ((data(gindex(igg), ii), igg = ig, ig + ng_irecord - 1), ii = 1, ns)
             else
                read(iunit , *) ((data(igg, ii), igg = ig, ig + ng_irecord - 1), ii = 1, ns)
             endif
          endif
          ig = ig + ng_irecord
       enddo
    endif
    if (peinf%inode.eq.0) call timacc(81,2)
    if (peinf%inode.eq.0) call timacc(82,1)
    if(peinf%npes > 1) then
       if(present(ng_record)) call MPI_BCAST(ng_record, nrecord_internal, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
       if(bcast_) then
          call MPI_BCAST(data, ng_bound * ns, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
       endif
    endif
    if (peinf%inode.eq.0) call timacc(82,2)

    return
  end subroutine read_format_real_data
  ! defined || defined BINARY
  ! these undefs prevent lots of cpp warnings
  !> write
  !=========================================================================
  !
  ! Included from file wfn_rho_vxc_io.f90.
  ! You are not expected to understand this. --DAS
  !
  !=========================================================================
  !real(DP)
  !=========================================================================
  !> This routine is turned into 12 variants by preprocessing.
  !! If nrecord and ng_record are not provided for read, the routine simply does not tell you what it found.
  !! For write, nrecord is given the default 1 if not provided. ng_record defaults to ng if nrecord = 1, but
  !! if nrecord > 1, ng_record must be provided to specify how to divide up the g-vectors.
  subroutine write_format_real_data(iunit, ng, ng_bound, ns, data, nrecord, ng_record, bcast, gindex, dont_read)
    integer, intent(in) :: iunit
    integer, intent(in) :: ng !< used size of array
    integer, intent(in) :: ng_bound !< actual size of array, >= ng
    integer, intent(in) :: ns
    real(DP), intent(in) :: data(:, :) !< (ng_bound, ns)
    integer, optional, intent(in) :: nrecord !< data/gvectors will be distributed among this many records
    integer, optional, pointer, intent(in) :: ng_record(:) !< number of gvectors in each record
    !! must be provided for write if nrecord > 1
    logical, optional, intent(in) :: bcast !< whether to do MPI_Bcast of what is read
    integer, optional, intent(in) :: gindex(:) !< map of order in file to order in gvec
    logical, optional, intent(in) :: dont_read !< if true, records will just be skipped; only for unformatted
    integer :: ig, igg, irecord, ii, nrecord_internal, ng_irecord
    character*100 :: tmpstr
    logical :: bcast_, dont_read_

    if(ng_bound < ng) then
       write(tmpstr,'(3a,i10,a,i10,a)') "In routine ", "write_format_real_data", ", ng_bound = ", ng_bound, " < ", ng, " = ng"
       call die(tmpstr)
    endif
    dont_read_ = .false.
    if(present(dont_read)) then
       call die("Write routine " + "write_format_real_data" + " cannot take argument dont_read.")
    endif
    if(peinf%inode .eq. 0) then
       if(present(gindex) .and. .not. dont_read_) then
          if(ubound(gindex, 1) < ng) &
               call die("In routine " + "write_format_real_data" + ", gindex array is too small")
          if(any(gindex(1:ng) > ng_bound) .or. any(gindex(1:ng) <= 0)) then
             call die("gindex out of bounds in " + "write_format_real_data")
          endif
       endif
    endif
    bcast_ = .not. dont_read_
    ! if not read, no point in broadcasting
    if(present(bcast)) then
       call die("Write routine " + "write_format_real_data" + " cannot take argument bcast.")
    endif
    ! if not reading, the size of the array passed is irrelevant
    if(.not. dont_read_) then
       if(ubound(data, 1) /= ng_bound) &
            call die("In routine " + "write_format_real_data" + ", mismatch of dimension 1 for data array")
       if(ubound(data, 2) /= ns) &
            call die("In routine " + "write_format_real_data" + ", mismatch of dimension 2 for data array")
    endif
    if(.not. present(nrecord)) then
       nrecord_internal = 1
       ng_irecord = ng
    else
       nrecord_internal = nrecord
       !> check validity of information if going to write
       if(nrecord_internal .ne. 1 .and. .not. present(ng_record)) then
          call die("Routine " + "write_format_real_data" + " requires ng_record array if nrecord > 1.")
       endif
       if(present(ng_record)) then
          if(nrecord_internal .ne. sum(ng_record(1:nrecord))) then
             write(tmpstr,'(a, i10, a, i10, a, a)') "sum(ng_record) = ", sum(ng_record), " ! = ", nrecord , &
                  " = nrecord in arguments to routine ", "write_format_real_data"
             call die(tmpstr)
          endif
       endif
    endif
    if(peinf%inode .eq. 0) then
       write(iunit , *) nrecord_internal
       if(nrecord_internal < 0) then
          write(tmpstr,'(3a, i10)') "In routine ", "write_format_real_data", " illegal nrecord ", nrecord_internal
          call die(tmpstr)
       endif
    endif
    if (peinf%inode.eq.0) call timacc(81,1)
    ! FHJ: zero buffers
    if(peinf%inode .eq. 0) then
       ig = 1
       !This is a PGI pragma to force the optimization level of this routine to -O1.
       !With -O2 or higher (including -fast) the line after the pragma causes a segmentation fault.
       !pgi$r opt=1
       do irecord = 1, nrecord_internal
          if(present(ng_record)) ng_irecord = ng_record(irecord)
          write(iunit , *) ng_irecord
          if(dont_read_) then
             write(iunit , *)
          else
             if(present(gindex)) then
                write(iunit , *) ((data(gindex(igg), ii), igg = ig, ig + ng_irecord - 1), ii = 1, ns)
             else
                write(iunit , *) ((data(igg, ii), igg = ig, ig + ng_irecord - 1), ii = 1, ns)
             endif
          endif
          ig = ig + ng_irecord
       enddo
    endif
    if (peinf%inode.eq.0) call timacc(81,2)
    if (peinf%inode.eq.0) call timacc(82,1)
    if (peinf%inode.eq.0) call timacc(82,2)

    return
  end subroutine write_format_real_data
  ! defined || defined BINARY
  ! these undefs prevent lots of cpp warnings
  !=========================================================================
  !
  ! Included from file wfn_rho_vxc_io.f90.
  ! You are not expected to understand this. --DAS
  !
  !=========================================================================
  !real(DP)
  !=========================================================================
  !> This routine is turned into 12 variants by preprocessing.
  !! If nrecord and ng_record are not provided for read, the routine simply does not tell you what it found.
  !! For write, nrecord is given the default 1 if not provided. ng_record defaults to ng if nrecord = 1, but
  !! if nrecord > 1, ng_record must be provided to specify how to divide up the g-vectors.
  subroutine read_binary_real_data(iunit, ng, ng_bound, ns, data, nrecord, ng_record, bcast, gindex, dont_read)
    integer, intent(in) :: iunit
    integer, intent(in) :: ng !< used size of array
    integer, intent(in) :: ng_bound !< actual size of array, >= ng
    integer, intent(in) :: ns
    real(DP), intent(out) :: data(:, :) !< (ng_bound, ns)
    integer, optional, intent(out) :: nrecord !< data/gvectors will be distributed among this many records
    integer, optional, pointer, intent(out) :: ng_record(:) !< number of gvectors in each record
    !! must be provided for write if nrecord > 1
    logical, optional, intent(in) :: bcast !< whether to do MPI_Bcast of what is read
    integer, optional, intent(in) :: gindex(:) !< map of order in file to order in gvec
    logical, optional, intent(in) :: dont_read !< if true, records will just be skipped; only for unformatted
    integer :: ig, igg, irecord, ii, nrecord_internal, ng_irecord
    character*100 :: tmpstr
    logical :: bcast_, dont_read_

    if(ng_bound < ng) then
       write(tmpstr,'(3a,i10,a,i10,a)') "In routine ", "read_binary_real_data", ", ng_bound = ", ng_bound, " < ", ng, " = ng"
       call die(tmpstr)
    endif
    dont_read_ = .false.
    if(present(dont_read)) then
       dont_read_ = dont_read
    endif
    if(peinf%inode .eq. 0) then
       if(present(gindex) .and. .not. dont_read_) then
          if(ubound(gindex, 1) < ng) &
               call die("In routine " + "read_binary_real_data" + ", gindex array is too small")
          if(any(gindex(1:ng) > ng_bound) .or. any(gindex(1:ng) <= 0)) then
             call die("gindex out of bounds in " + "read_binary_real_data")
          endif
       endif
    endif
    bcast_ = .not. dont_read_
    ! if not read, no point in broadcasting
    if(present(bcast)) then
       bcast_ = bcast
    endif
    ! if not reading, the size of the array passed is irrelevant
    if(.not. dont_read_) then
       if(ubound(data, 1) /= ng_bound) &
            call die("In routine " + "read_binary_real_data" + ", mismatch of dimension 1 for data array")
       if(ubound(data, 2) /= ns) &
            call die("In routine " + "read_binary_real_data" + ", mismatch of dimension 2 for data array")
    endif
    if(peinf%inode .eq. 0) then
       read(iunit ) nrecord_internal
       if(nrecord_internal < 0) then
          write(tmpstr,'(3a, i10)') "In routine ", "read_binary_real_data", " illegal nrecord ", nrecord_internal
          call die(tmpstr)
       endif
    endif
    if(present(nrecord)) then
       nrecord = nrecord_internal
       if(peinf%npes > 1) call MPI_BCAST(nrecord, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
    endif
    if(present(ng_record)) then
       allocate(ng_record (nrecord_internal))
    endif
    if (peinf%inode.eq.0) call timacc(81,1)
    ! FHJ: zero buffers
    data(:,:) = (0.0d0,0.0d0)
    if(peinf%inode .eq. 0) then
       ig = 1
       !This is a PGI pragma to force the optimization level of this routine to -O1.
       !With -O2 or higher (including -fast) the line after the pragma causes a segmentation fault.
       !pgi$r opt=1
       do irecord = 1, nrecord_internal
          read(iunit ) ng_irecord
          if(present(ng_record)) ng_record(irecord) = ng_irecord
          if(dont_read_) then
             read(iunit )
          else
             if(present(gindex)) then
                read(iunit ) ((data(gindex(igg), ii), igg = ig, ig + ng_irecord - 1), ii = 1, ns)
             else
                read(iunit ) ((data(igg, ii), igg = ig, ig + ng_irecord - 1), ii = 1, ns)
             endif
          endif
          ig = ig + ng_irecord
       enddo
    endif
    if (peinf%inode.eq.0) call timacc(81,2)
    if (peinf%inode.eq.0) call timacc(82,1)
    if(peinf%npes > 1) then
       if(present(ng_record)) call MPI_BCAST(ng_record, nrecord_internal, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
       if(bcast_) then
          call MPI_BCAST(data, ng_bound * ns, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpierr)
       endif
    endif
    if (peinf%inode.eq.0) call timacc(82,2)

    return
  end subroutine read_binary_real_data
  ! defined FORMATTED || defined
  ! these undefs prevent lots of cpp warnings
  !> write
  !=========================================================================
  !
  ! Included from file wfn_rho_vxc_io.f90.
  ! You are not expected to understand this. --DAS
  !
  !=========================================================================
  !real(DP)
  !=========================================================================
  !> This routine is turned into 12 variants by preprocessing.
  !! If nrecord and ng_record are not provided for read, the routine simply does not tell you what it found.
  !! For write, nrecord is given the default 1 if not provided. ng_record defaults to ng if nrecord = 1, but
  !! if nrecord > 1, ng_record must be provided to specify how to divide up the g-vectors.
  subroutine write_binary_real_data(iunit, ng, ng_bound, ns, data, nrecord, ng_record, bcast, gindex, dont_read)
    integer, intent(in) :: iunit
    integer, intent(in) :: ng !< used size of array
    integer, intent(in) :: ng_bound !< actual size of array, >= ng
    integer, intent(in) :: ns
    real(DP), intent(in) :: data(:, :) !< (ng_bound, ns)
    integer, optional, intent(in) :: nrecord !< data/gvectors will be distributed among this many records
    integer, optional, pointer, intent(in) :: ng_record(:) !< number of gvectors in each record
    !! must be provided for write if nrecord > 1
    logical, optional, intent(in) :: bcast !< whether to do MPI_Bcast of what is read
    integer, optional, intent(in) :: gindex(:) !< map of order in file to order in gvec
    logical, optional, intent(in) :: dont_read !< if true, records will just be skipped; only for unformatted
    integer :: ig, igg, irecord, ii, nrecord_internal, ng_irecord
    character*100 :: tmpstr
    logical :: bcast_, dont_read_

    if(ng_bound < ng) then
       write(tmpstr,'(3a,i10,a,i10,a)') "In routine ", "write_binary_real_data", ", ng_bound = ", ng_bound, " < ", ng, " = ng"
       call die(tmpstr)
    endif
    dont_read_ = .false.
    if(present(dont_read)) then
       call die("Write routine " + "write_binary_real_data" + " cannot take argument dont_read.")
    endif
    if(peinf%inode .eq. 0) then
       if(present(gindex) .and. .not. dont_read_) then
          if(ubound(gindex, 1) < ng) &
               call die("In routine " + "write_binary_real_data" + ", gindex array is too small")
          if(any(gindex(1:ng) > ng_bound) .or. any(gindex(1:ng) <= 0)) then
             call die("gindex out of bounds in " + "write_binary_real_data")
          endif
       endif
    endif
    bcast_ = .not. dont_read_
    ! if not read, no point in broadcasting
    if(present(bcast)) then
       call die("Write routine " + "write_binary_real_data" + " cannot take argument bcast.")
    endif
    ! if not reading, the size of the array passed is irrelevant
    if(.not. dont_read_) then
       if(ubound(data, 1) /= ng_bound) &
            call die("In routine " + "write_binary_real_data" + ", mismatch of dimension 1 for data array")
       if(ubound(data, 2) /= ns) &
            call die("In routine " + "write_binary_real_data" + ", mismatch of dimension 2 for data array")
    endif
    if(.not. present(nrecord)) then
       nrecord_internal = 1
       ng_irecord = ng
    else
       nrecord_internal = nrecord
       !> check validity of information if going to write
       if(nrecord_internal .ne. 1 .and. .not. present(ng_record)) then
          call die("Routine " + "write_binary_real_data" + " requires ng_record array if nrecord > 1.")
       endif
       if(present(ng_record)) then
          if(nrecord_internal .ne. sum(ng_record(1:nrecord))) then
             write(tmpstr,'(a, i10, a, i10, a, a)') "sum(ng_record) = ", sum(ng_record), " ! = ", nrecord , &
                  " = nrecord in arguments to routine ", "write_binary_real_data"
             call die(tmpstr)
          endif
       endif
    endif
    if(peinf%inode .eq. 0) then
       write(iunit ) nrecord_internal
       if(nrecord_internal < 0) then
          write(tmpstr,'(3a, i10)') "In routine ", "write_binary_real_data", " illegal nrecord ", nrecord_internal
          call die(tmpstr)
       endif
    endif
    if (peinf%inode.eq.0) call timacc(81,1)
    ! FHJ: zero buffers
    if(peinf%inode .eq. 0) then
       ig = 1
       !This is a PGI pragma to force the optimization level of this routine to -O1.
       !With -O2 or higher (including -fast) the line after the pragma causes a segmentation fault.
       !pgi$r opt=1
       do irecord = 1, nrecord_internal
          if(present(ng_record)) ng_irecord = ng_record(irecord)
          write(iunit ) ng_irecord
          if(dont_read_) then
             write(iunit )
          else
             if(present(gindex)) then
                write(iunit ) ((data(gindex(igg), ii), igg = ig, ig + ng_irecord - 1), ii = 1, ns)
             else
                write(iunit ) ((data(igg, ii), igg = ig, ig + ng_irecord - 1), ii = 1, ns)
             endif
          endif
          ig = ig + ng_irecord
       enddo
    endif
    if (peinf%inode.eq.0) call timacc(81,2)
    if (peinf%inode.eq.0) call timacc(82,1)
    if (peinf%inode.eq.0) call timacc(82,2)

    return
  end subroutine write_binary_real_data
  ! defined FORMATTED || defined
  ! these undefs prevent lots of cpp warnings
  !=========================================================================
  !
  ! Included from file wfn_rho_vxc_io.f90.
  ! You are not expected to understand this. --DAS
  !
  !=========================================================================
  ! defined FORMATTED || defined BINARY
  ! defined real(DP)
  !=========================================================================
  subroutine read_real_data(iunit, is_formatted, ng, ng_bound, ns, data, nrecord, ng_record, bcast, gindex, dont_read)
    integer, intent(in) :: iunit
    logical, intent(in) :: is_formatted !< true = formatted, false = binary
    integer, intent(in) :: ng !< used size of array
    integer, intent(in) :: ng_bound !< actual size of array, >= ng
    integer, intent(in) :: ns
    real(DP), intent(out) :: data(:, :) !< (ng_bound, ns)
    integer, optional, intent(out) :: nrecord !< data/gvectors will be distributed among this many records
    integer, optional, pointer, intent(out) :: ng_record(:) !< number of gvectors in each record
    logical, optional, intent(in) :: bcast !< whether to do MPI_Bcast of what is read
    integer, optional, intent(in) :: gindex(:) !< map of order in file to order in gvec
    logical, optional, intent(in) :: dont_read !< if true, records will just be skipped; only for unformatted

    ! note: surprisingly, the Fortran standard permits passing optional arguments to other routines
    ! in which those arguments are also optional.
    if(is_formatted) then
       call read_format_real_data(iunit, ng, ng_bound, ns, data, nrecord = nrecord, &
            ng_record = ng_record, bcast = bcast, gindex = gindex, dont_read = dont_read)
    else
       call read_binary_real_data(iunit, ng, ng_bound, ns, data, nrecord = nrecord, &
            ng_record = ng_record, bcast = bcast, gindex = gindex, dont_read = dont_read)
    endif

    return
  end subroutine read_real_data
  ! defined FORMATTED || defined BINARY
  ! these undefs prevent lots of cpp warnings
  !> write
  !=========================================================================
  !
  ! Included from file wfn_rho_vxc_io.f90.
  ! You are not expected to understand this. --DAS
  !
  !=========================================================================
  ! defined FORMATTED || defined BINARY
  ! defined real(DP)
  !=========================================================================
  subroutine write_real_data(iunit, is_formatted, ng, ng_bound, ns, data, nrecord, ng_record, bcast, gindex, dont_read)
    integer, intent(in) :: iunit
    logical, intent(in) :: is_formatted !< true = formatted, false = binary
    integer, intent(in) :: ng !< used size of array
    integer, intent(in) :: ng_bound !< actual size of array, >= ng
    integer, intent(in) :: ns
    real(DP), intent(in) :: data(:, :) !< (ng_bound, ns)
    integer, optional, intent(in) :: nrecord !< data/gvectors will be distributed among this many records
    integer, optional, pointer, intent(in) :: ng_record(:) !< number of gvectors in each record
    logical, optional, intent(in) :: bcast !< whether to do MPI_Bcast of what is read
    integer, optional, intent(in) :: gindex(:) !< map of order in file to order in gvec
    logical, optional, intent(in) :: dont_read !< if true, records will just be skipped; only for unformatted

    ! note: surprisingly, the Fortran standard permits passing optional arguments to other routines
    ! in which those arguments are also optional.
    if(is_formatted) then
       call write_format_real_data(iunit, ng, ng_bound, ns, data, nrecord = nrecord, &
            ng_record = ng_record, bcast = bcast, gindex = gindex, dont_read = dont_read)
    else
       call write_binary_real_data(iunit, ng, ng_bound, ns, data, nrecord = nrecord, &
            ng_record = ng_record, bcast = bcast, gindex = gindex, dont_read = dont_read)
    endif

    return
  end subroutine write_real_data
  ! defined FORMATTED || defined BINARY
  ! these undefs prevent lots of cpp warnings

  !=========================================================================
  !> deallocate variables allocated by read_header
  subroutine dealloc_header(sheader, atyp, apos, ngk, kw, kpt, ifmin, ifmax, energies, occupations)
    character(len=3), intent(in) :: sheader
    integer, pointer, intent(inout) :: atyp(:)
    real(DP), pointer, intent(inout) :: apos(:,:)
    integer, pointer, intent(inout) :: ngk(:)
    real(DP), pointer, intent(inout) :: kw(:), kpt(:, :)
    integer, pointer, intent(inout) :: ifmin(:, :), ifmax(:, :)
    real(DP), pointer, intent(inout) :: energies(:, :, :)
    real(DP), pointer, intent(inout) :: occupations(:, :, :)
    logical :: wfnflag

    if (sheader .eq. 'WFN') then
       wfnflag = .true.
    elseif (sheader .eq. 'RHO' .or. sheader .eq. 'VXC') then
       wfnflag = .false.
    else
       call die("unknown file header: '" + sheader + "' (should be 'WFN'/'RHO'/'VXC')")
    endif
    if (associated(atyp)) then;
       deallocate(atyp);
       nullify(atyp);
    endif
    if (associated(apos)) then;
       deallocate(apos);
       nullify(apos);
    endif
    if (wfnflag) then
       if (associated(ngk)) then;
          deallocate(ngk);
          nullify(ngk);
       endif

       if (associated(kw)) then;
          deallocate(kw);
          nullify(kw);
       endif

       if (associated(kpt)) then;
          deallocate(kpt);
          nullify(kpt);
       endif

       if (associated(ifmin)) then;
          deallocate(ifmin);
          nullify(ifmin);
       endif
       if (associated(ifmax)) then;
          deallocate(ifmax);
          nullify(ifmax);
       endif
       if (associated(energies)) then;
          deallocate(energies);
          nullify(energies);
       endif
       if (associated(occupations)) then;
          deallocate(occupations);
          nullify(occupations);
       endif
    endif

    return
  end subroutine dealloc_header
  !=========================================================================
  subroutine dealloc_header_type(sheader, crys, kp)
    character(len=3), intent(in) :: sheader
    type(crystal), intent(inout) :: crys
    type(kpoints), intent(inout) :: kp

    call dealloc_header(sheader, crys%atyp, crys%apos, kp%ngk, kp%w, kp%rk, kp%ifmin, kp%ifmax, kp%el, kp%occ)

    return
  end subroutine dealloc_header_type
  !=========================================================================
  subroutine dealloc_kp(kp)
    type(kpoints), intent(inout) :: kp

    if (associated(kp%ngk)) then;
       deallocate(kp%ngk);
       nullify(kp%ngk);
    endif

    if (associated(kp%w)) then;
       deallocate(kp%w);
       nullify(kp%w);
    endif

    if (associated(kp%rk)) then;
       deallocate(kp%rk);
       nullify(kp%rk);
    endif

    if (associated(kp%ifmin)) then;
       deallocate(kp%ifmin);
       nullify(kp%ifmin);
    endif

    if (associated(kp%ifmax)) then;
       deallocate(kp%ifmax);
       nullify(kp%ifmax);
    endif

    if (associated(kp%el)) then;
       deallocate(kp%el);
       nullify(kp%el);
    endif

    if (associated(kp%occ)) then;
       deallocate(kp%occ);
       nullify(kp%occ);
    endif

    return
  end subroutine dealloc_kp
  !=========================================================================
  subroutine dealloc_crys(crys)
    type(crystal), intent(inout) :: crys

    if (associated(crys%atyp)) then;
       deallocate(crys%atyp);
       nullify(crys%atyp);
    endif

    if (associated(crys%apos)) then;
       deallocate(crys%apos);
       nullify(crys%apos);
    endif

    return
  end subroutine dealloc_crys
  !=========================================================================
  !> detect incompatibility between header info for wfns supposedly describing same system
  subroutine check_header(name, kp, gvec, syms, crys, name2, kp2, gvec2, syms2, crys2, is_wfn, tolerant)
    character(len=*), intent(in) :: name
    type(kpoints), intent(in) :: kp
    type(gspace), intent(in) :: gvec
    type(symmetry), intent(in) :: syms
    type(crystal), intent(in) :: crys
    character(len=*), intent(in) :: name2
    type(kpoints), intent(in) :: kp2
    type(gspace), intent(in) :: gvec2
    type(symmetry), intent(in) :: syms2
    type(crystal), intent(in) :: crys2
    logical, intent(in) :: is_wfn
    !< set to false if RHO or VXC is one of the two being compared to avoid checking uninitialized fields
    logical, optional, intent(in) :: tolerant !< set to true to allow difference in symmetries and atoms
    character*100 :: string
    logical :: tolerant_
    integer :: isym

    tolerant_ = .false.
    if (present(tolerant)) then
       tolerant_ = tolerant
    endif
    string = TRUNC(name) + " vs. " + TRUNC(name2)
    !> kpoints
    if (kp%nspin .ne. kp2%nspin) call die(TRUNC(string) + ": spin mismatch")
    if (kp%nspinor .ne. kp2%nspinor) call die(TRUNC(string) + ": nspinor mismatch")
    if (is_wfn .and. abs(kp%ecutwfc - kp2%ecutwfc) > TOL_Small) call die(TRUNC(string) + ": wfn cutoff mismatch")
    !> gspace
    if (gvec%ng .ne. gvec2%ng) call die(TRUNC(string) + ": total number of G-vectors mismatch")
    if (abs(gvec%ecutrho - gvec2%ecutrho) > TOL_Small) call die(TRUNC(string) + ": charge-density cutoff mismatch")
    if (any(gvec%FFTgrid(1:3) .ne. gvec2%FFTgrid(1:3))) call die(TRUNC(string) + ": FFT grid mismatch")
    if (.not. tolerant_) then
       !> symmetries
       if (syms%ntran .ne. syms2%ntran) call die(TRUNC(string) + ": number of symmetries mismatch")
       if (syms%cell_symmetry .ne. syms2%cell_symmetry) call die(TRUNC(string) + ": type of cell symmetry mismatch")
       if (any(syms%mtrx(1:3, 1:3, 1:syms%ntran) .ne. syms2%mtrx(1:3, 1:3, 1:syms2%ntran))) then
          call die(TRUNC(string) + ": symmetry rotation matrix mismatch")
       endif
       if (any(abs(syms%tnp(1:3, 1:syms%ntran) - syms2%tnp(1:3, 1:syms2%ntran)) > TOL_Small)) then
          write(*,'(A,I5,A,I5)') "syms%ntran = ", syms%ntran, " syms2%ntran = ", syms2%ntran
          do isym = 1, syms%ntran
             write(*,'(6F12.5)') syms%tnp(:,isym), syms2%tnp(:,isym)
          enddo
          call die(TRUNC(string) + ": symmetry fractional translation mismatch")
       endif
       !> atoms
       if (crys%nat .ne. crys2%nat) call die(TRUNC(string) + ": number of atoms mismatch")
       if (any(crys%atyp(1:crys%nat) .ne. crys2%atyp(1:crys2%nat))) then
          write(*,*) "crys%nat WFN1: ", crys%nat, "WFN2: ", crys2%nat
          write(*,*) "crys%atyp WFN1: ", crys%atyp(1:crys%nat), "WFN2: ", crys2%atyp(1:crys2%nat)
          call die(TRUNC(string) + ": atom species mismatch")
       endif
       if (any(abs(crys%alat * crys%apos(1:3, 1:crys%nat) - crys2%alat * crys2%apos(1:3, 1:crys2%nat)) > TOL_Small)) &
            call die(TRUNC(string) + ": atom position mismatch")
    endif
    !> lattice
    if (abs(crys%celvol - crys2%celvol) > TOL_Small) call die(TRUNC(string) + ": cell volume mismatch")
    if (abs(crys%recvol - crys2%recvol) > TOL_Small) call die(TRUNC(string) + ": reciprocal cell volume mismatch")
    if (any(abs(crys%alat * crys%avec(1:3, 1:3) - crys2%alat * crys2%avec(1:3, 1:3)) > TOL_Small)) then
       call die(TRUNC(string) + ": lattice vector mismatch")
    endif
    if (any(abs(crys%blat * crys%bvec(1:3, 1:3) - crys2%blat * crys2%bvec(1:3, 1:3)) > TOL_Small)) then
       call die(TRUNC(string) + ": reciprocal lattice vector mismatch")
    endif
    if (any(abs(crys%adot(1:3, 1:3) - crys2%adot(1:3, 1:3)) > TOL_Small)) then
       call die(TRUNC(string) + ": real-space metric mismatch")
    endif
    if (any(abs(crys%bdot(1:3, 1:3) - crys2%bdot(1:3, 1:3)) > TOL_Small)) then
       call die(TRUNC(string) + ": reciprocal-space metric mismatch")
    endif

    return
  end subroutine check_header
  !> require `version` to be the same as `version_ref`
  subroutine require_version(fname, version, version_ref)
    character(len=*), intent(in) :: fname
    integer, intent(in) :: version
    integer, intent(in) :: version_ref

    if (version/=version_ref) then
       if (peinf%inode==0) then
          write(0,*) 'ERROR: Wrong version for file "',TRUNC(fname),'".'
          write(0,*) 'Expected: ', version_ref
          write(0,*) 'Got: ', version
       endif
       call die('Wrong version for file "'+TRUNC(fname)+'".', only_root_writes=.true.)
    endif

  end subroutine require_version
  !> A high-level wrapper for write_*_header* functions
  subroutine write_mf_header(iunit, mf)
    integer, intent(in) :: iunit
    type(mf_header_t), intent(in) :: mf
    character(len=3) :: sheader
    character(len=16) :: fmt_str
    logical :: is_fmt = .false.

    if (peinf%inode==0) then
       inquire(unit=iunit, form=fmt_str)
       if (TRUNC(fmt_str)=='FORMATTED') then
          is_fmt = .true.
       else if (TRUNC(fmt_str)/='UNFORMATTED') then
          call die('Unknown value for formatted string: '+TRUNC(fmt_str), &
               only_root_writes=.true.)
       endif
    endif
    ! FHJ: No need to bcast is_fmt because only inode==0 reads the file.
    sheader = mf%sheader
    call write_header_type(iunit, is_fmt, sheader, mf%iflavor, &
         mf%kp, mf%gvec, mf%syms, mf%crys, version=mf%version)

  end subroutine write_mf_header
  !> A high-level wrapper for write_*_header* functions
  subroutine read_mf_header(iunit, mf, iflavor, sheader, warn, dont_warn_kgrid)
    integer, intent(in) :: iunit
    type(mf_header_t), intent(out) :: mf
    integer, intent(in), optional :: iflavor
    character(len=3), intent(in), optional :: sheader
    logical, intent(in), optional :: warn
    logical, intent(in), optional :: dont_warn_kgrid
    character(len=16) :: fmt_str
    logical :: is_fmt = .false.

    if (peinf%inode==0) then
       inquire(unit=iunit, form=fmt_str)
       if (TRUNC(fmt_str)=='FORMATTED') then
          is_fmt = .true.
       else if (TRUNC(fmt_str)/='UNFORMATTED') then
          call die('Unknown value for formatted string: '+TRUNC(fmt_str), &
               only_root_writes=.true.)
       endif
    endif
    ! FHJ: No need to bcast is_fmt because only inode==0 reads the file.
    if (present(sheader)) then
       mf%sheader = sheader
    else
       mf%sheader = 'GET'
    endif
    if (present(iflavor)) then
       mf%iflavor = iflavor
    else
       mf%iflavor = -1
    endif
    call read_header_type(iunit, is_fmt, mf%sheader, mf%iflavor, &
         mf%kp, mf%gvec, mf%syms, mf%crys, version=mf%version, sdate=mf%sdate, stime=mf%stime, &
         warn=warn, dont_warn_kgrid=dont_warn_kgrid)

  end subroutine read_mf_header
  !> Routine to initialize the mf_header_t type from a bunch of separated data types
  subroutine init_mf_header_from_types(mf_header, sheader, iflavor, version, kp, gvec, syms, crys)
    type(mf_header_t), intent(out) :: mf_header
    character(len=*), intent(in) :: sheader
    integer, intent(in) :: iflavor
    integer, intent(in) :: version
    type(kpoints), intent(in) :: kp
    type(gspace), intent(in) :: gvec
    type(symmetry), intent(in) :: syms
    type(crystal), intent(in) :: crys

    mf_header%version = version
    mf_header%sheader = sheader
    mf_header%sdate = ''
    mf_header%stime = ''
    mf_header%iflavor = iflavor
    mf_header%kp = kp
    mf_header%gvec = gvec
    mf_header%syms = syms
    mf_header%crys = crys

  end subroutine init_mf_header_from_types
end module wfn_rho_vxc_io_m

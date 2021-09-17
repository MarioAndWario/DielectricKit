#include "f_defs.h"

!>=========================================================================
!!
!! Module:
!!
!! (1) wfn_rho_vxc_io_m     Originally by DAS     Last Modified 10/17/2010 (DAS)
!!
!!     Routines to read and write wavefunctions, density, and Vxc, and
!!     deallocate associated variables. The "type" routines use typedefs.
!!     The code is generated through repeated inclusion of a file with
!!     different preprocessor definitions each time. You are not expected to
!!     understand this. Consult the resulting .p.f file for clarity.
!!
!! Specification for ASCII matrix element files (vxc.dat, x.dat):
!! Matrix elements are in eV and are always written with real and imaginary parts.
!! They may contain any number of k-points in any order.
!! Each k-point block begins with the line:
!!   kx, ky, kz [crystal coordinates], ndiag*nspin, noffdiag*nspin
!! There are then ndiag*nspin lines of the form
!!   ispin, idiag, Re <idiag|V|idiag>, Im <idiag|V|idiag>
!! There are then noffdiag*nspin lines of the form
!!   ispin, ioff1, ioff2, Re <ioff1|V|ioff2>, Im <ioff1|V|ioff2>
!!
!!=========================================================================

module wfn_rho_vxc_io_m

  use global_m
  use check_inversion_m
  use sort_m
  use misc_m

  implicit none

  private
  !> For library usage, do not make global_m contents available
  !! to avoid namespace clashes.

  public ::                        &
       bgw_conf_test,                 &
       read_binary_header,            &
       write_binary_header,           &
       read_format_header,            &
       write_format_header,           &
       read_binary_header_type,       &
       write_binary_header_type,      &
       read_format_header_type,       &
       write_format_header_type,      &
       read_binary_gvectors,          &
       write_binary_gvectors,         &
       read_format_gvectors,          &
       write_format_gvectors,         &
       read_binary_real_data,         &
       write_binary_real_data,        &
       read_format_real_data,         &
       write_format_real_data,        &
       read_binary_complex_data,      &
       write_binary_complex_data,     &
       read_format_complex_data,      &
       write_format_complex_data,     &
       read_binary_data,              &
       write_binary_data,             &
       read_format_data,              &
       write_format_data,             &
       read_header,                   &
       write_header,                  &
       read_header_type,              &
       write_header_type,             &
       read_gvectors,                 &
       write_gvectors,                &
       read_real_data,                &
       write_real_data,               &
       read_complex_data,             &
       write_complex_data,            &
       read_data,                     &
       write_data,                    &
       dealloc_header,                &
       dealloc_header_type,           &
       dealloc_crys,                  &
       dealloc_kp,                    &
       check_header,                  &
       write_matrix_elements,         &
       write_matrix_elements_type,    &
       read_matrix_elements,          &
       read_matrix_elements_type,     &
       require_version,               &
       write_mf_header,               &
       read_mf_header,                &
       init_mf_header_from_types,     &
       prepare_syms



  !> These interfaces can be used to avoid explicit ifdef CPLX switches
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

  interface read_matrix_elements
     module procedure read_matrix_elements_real, read_matrix_elements_cplx
  end interface read_matrix_elements

  interface read_matrix_elements_type
     module procedure read_matrix_elements_type_real, read_matrix_elements_type_cplx
  end interface read_matrix_elements_type

contains

#define FORMATTED
#define READ
#include "wfn_rho_vxc_io_inc.f90"
  !> write
#undef READ
#include "wfn_rho_vxc_io_inc.f90"

#undef FORMATTED
#define BINARY
#define READ
#include "wfn_rho_vxc_io_inc.f90"
  !> write
#undef READ
#include "wfn_rho_vxc_io_inc.f90"

#undef BINARY
#define READ
#include "wfn_rho_vxc_io_inc.f90"
  !> write
#undef READ
#include "wfn_rho_vxc_io_inc.f90"

#define TEMP_COMPLEX
#define FORMATTED
#define READ
#include "wfn_rho_vxc_io_inc.f90"
  !> write
#undef READ
#include "wfn_rho_vxc_io_inc.f90"

#undef FORMATTED
#define BINARY
#define READ
#include "wfn_rho_vxc_io_inc.f90"
  !> write
#undef READ
#include "wfn_rho_vxc_io_inc.f90"

#undef BINARY
#define READ
#include "wfn_rho_vxc_io_inc.f90"
  !> write
#undef READ
#include "wfn_rho_vxc_io_inc.f90"

#undef TEMP_COMPLEX
#define TEMP_REAL
#define FORMATTED
#define READ
#include "wfn_rho_vxc_io_inc.f90"
  !> write
#undef READ
#include "wfn_rho_vxc_io_inc.f90"

#undef FORMATTED
#define BINARY
#define READ
#include "wfn_rho_vxc_io_inc.f90"
  !> write
#undef READ
#include "wfn_rho_vxc_io_inc.f90"

#undef BINARY
#define READ
#include "wfn_rho_vxc_io_inc.f90"
  !> write
#undef READ
#include "wfn_rho_vxc_io_inc.f90"

  !=========================================================================
  !> this routine is purely for use in configure scripts to test module accessibility
  subroutine bgw_conf_test()

    write(6,*) 'Yes, it works.'

  end subroutine bgw_conf_test

  !=========================================================================
  !> this routine is used to write the ASCII files vxc.dat and x.dat, without typedefs types
  subroutine write_matrix_elements(iunit, kk, nspin, ndiag, noffdiag, spin_index, diag, offdiag1, offdiag2, mtxel)
    integer, intent(in) :: iunit !< file unit to write to
    real(DP), intent(in) :: kk(3) !< kpoint to write, in crystal coords
    integer, intent(in) :: nspin !< number of spins to write
    integer, intent(in) :: ndiag !< number of diagonal elements to write
    integer, intent(in) :: noffdiag !< number of offdiagonal elements to write
    integer, intent(in) :: spin_index(:) !< (nspin) mapping of 1:nspin to actual spins. 3 choices:
    !! spin-unpolarized: spin_index(1) = 1
    !! spin-polarized: spin_index(1) = 1, spin_index(2) = 2
    !! spin-polarized, spin 2 only: spin_index(1) = 1, spin_index(2) = 2
    integer, intent(in) :: diag(:) !< (ndiag) mapping of 1:ndiag onto band indices for diagonals
    integer, intent(in) :: offdiag1(:) !< (noffdiag) mapping of 1:noffdiag onto band indices for offdiagonal 1
    integer, intent(in) :: offdiag2(:) !< (noffdiag) mapping of 1:noffdiag onto band indices for offdiagonal 2
    complex(DPC), intent(in) :: mtxel(:,:) !< (ndiag+noffdiag,nspin)

    integer :: ispin, idiag, ioff

    PUSH_SUB(write_matrix_elements)

    if(any(spin_index(1:nspin) < 1 .or. spin_index(1:nspin) > 2)) &
         call die("write_matrix_elements: spin_index out of bounds")
    if(nspin < 1 .or. nspin > 2) call die("write_matrix_elements: nspin out of bounds")
    if(ndiag < 0) call die("write_matrix_elements: ndiag < 0")
    if(noffdiag < 0) call die("write_matrix_elements: noffdiag < 0")

    !> write header
    write(iunit,'(3f15.10,2i12)') kk(1:3), ndiag*nspin, noffdiag*nspin
    !> write diagonal matrix elements
    do idiag=1,ndiag
       do ispin=1,nspin
          write(iunit,'(i3,i10,2g22.10)') spin_index(ispin), diag(idiag), dble(mtxel(idiag,ispin)), IMAG(mtxel(idiag,ispin))
       enddo
    enddo
    !> write off-diagonal matrix elements
    do ioff=1,noffdiag
       do ispin=1,nspin
          write(iunit,'(i3,2i10,2g22.10)') spin_index(ispin), offdiag1(ioff), offdiag2(ioff), &
               dble(mtxel(ioff+ndiag,ispin)), IMAG(mtxel(ioff+ndiag,ispin))
       enddo
    enddo

    POP_SUB(write_matrix_elements)
    return
  end subroutine write_matrix_elements

  !=========================================================================
  !> this routine is used to write the ASCII files vxc.dat and x.dat, with typedefs types
  subroutine write_matrix_elements_type(iunit, kk, sig, mtxel)
    integer, intent(in) :: iunit !< file unit to write to
    real(DP), intent(in) :: kk(3) !< kpoint to write, in crystal coords
    type(siginfo), intent(in) :: sig !< structure containing other needed data
    complex(DPC), intent(in) :: mtxel(:,:) !< (sig%ndiag+sig%noffdiag,sig%nspin)

    PUSH_SUB(write_matrix_elements_type)

    call write_matrix_elements(iunit, kk, sig%nspin, sig%ndiag, sig%noffdiag, &
         sig%spin_index, sig%diag, sig%off1, sig%off2, mtxel)

    POP_SUB(write_matrix_elements_type)
    return
  end subroutine write_matrix_elements_type

  !=========================================================================
  !> this routine is used to read the ASCII files vxc.dat and x.dat, without typedefs types. do not call directly
  ! ======
  ! from read_matrix_elements_type_cplx:
  !   call read_matrix_elements_base(iunit, iostat, kk, nspin, ndiag, noffdiag,spin_index, diag, offdiag1, offdiag2, mtxel_cplx = mtxel)
  ! ------
  !
  subroutine read_matrix_elements_base(iunit, iostat, kk, nspin, ndiag, noffdiag, &
       spin_index, diag, offdiag1, offdiag2, mtxel_real, mtxel_cplx)
    integer, intent(in) :: iunit !< file unit to read from
    integer, intent(out) :: iostat !< status of reading header, to detect end of file
    real(DP), intent(out) :: kk(3) !< kpoint read, in crystal coords
    integer, intent(in) :: nspin !< number of spins to read
    integer, intent(in) :: ndiag !< number of diagonal elements to read
    integer, intent(in) :: noffdiag !< number of offdiagonal elements to read
    integer, intent(in) :: spin_index(:) !< (nspin) mapping of 1:nspin to actual spins. 3 choices:
    !! spin-unpolarized: spin_index(1) = 1
    !! spin-polarized: spin_index(1) = 1, spin_index(2) = 2
    !! spin-polarized, spin 2 only: spin_index(1) = 1, spin_index(2) = 2
    integer, intent(in) :: diag(:) !< (ndiag) mapping of 1:ndiag onto band indices for diagonals
    integer, intent(in) :: offdiag1(:) !< (noffdiag) mapping of 1:noffdiag onto band indices for offdiagonal 1
    integer, intent(in) :: offdiag2(:) !< (noffdiag) mapping of 1:noffdiag onto band indices for offdiagonal 2
    real(DP), optional, intent(out) :: mtxel_real(:,:) !< (ndiag+noffdiag,nspin) this or mtxel_cplx must be given
    complex(DPC), optional, intent(out) :: mtxel_cplx(:,:) !< (ndiag+noffdiag,nspin) this or mtxel_real must be given

    integer :: ispin, idiag, ioff, ime, ndiag_read, noffdiag_read, ispin_read, diag_read, off1_read, off2_read
    real(DP) :: mtxel_re, mtxel_im
    logical, allocatable :: found(:,:)

    PUSH_SUB(read_matrix_elements_base)

    if(any(spin_index(1:nspin) < 1 .or. spin_index(1:nspin) > 2)) &
         call die("read_matrix_elements: spin_index out of bounds")
    if(nspin < 1 .or. nspin > 2) call die("read_matrix_elements: nspin out of bounds")
    if(ndiag < 0) call die("read_matrix_elements: ndiag < 0")
    if(noffdiag < 0) call die("read_matrix_elements: noffdiag < 0")
    if((.not. present(mtxel_real) .and. .not. present(mtxel_cplx)) .or. (present(mtxel_real) .and. present(mtxel_cplx))) &
         call die("read_matrix elements must be called with exactly one of mtxel_real and mtxel_cplx")

    !> read header
    read(iunit,*,iostat = iostat) kk(1:3), ndiag_read, noffdiag_read
    if(iostat /= 0) then !< we are at the end of the file
       POP_SUB(read_matrix_elements)
       return
    endif

    if(ndiag_read < ndiag * nspin) call die("read_matrix_elements: not enough diagonals present")
    if(noffdiag_read < noffdiag * nspin) call die("read_matrix_elements: not enough offdiagonals present")

    SAFE_ALLOCATE(found, (ndiag+noffdiag,nspin))
    found(:,:) = .false.

    !> read diagonal matrix elements
    do ime = 1, ndiag_read
       read(iunit,*) ispin_read, diag_read, mtxel_re, mtxel_im
       if(present(mtxel_real) .and. abs(mtxel_im) > TOL_Zero) call die("cannot have complex matrix elements in real version")
       do ispin = 1, nspin
          do idiag = 1, ndiag
             if (spin_index(ispin) == ispin_read .and. diag(idiag) == diag_read) then
                if(present(mtxel_real)) then
                   mtxel_real(idiag, ispin) = mtxel_re
                else
                   mtxel_cplx(idiag, ispin) = DCMPLX(mtxel_re, mtxel_im)
                endif
                found(idiag, ispin) = .true.
             endif
          enddo
       enddo
    enddo

    if(any(.not. found(1:ndiag, :))) then
       if(peinf%inode == 0) then
          write(0,*) 'missing diagonal matrix elements (band, spin): '
          do ispin = 1, nspin
             do idiag = 1, ndiag
                write(0,*) '(', diag(idiag), ',', spin_index(ispin), ')'
             enddo
          enddo
       endif
       call die("read_matrix_elements: not all needed data present")
    endif

    !> read off-diagonal matrix elements
    do ime = 1, noffdiag_read
       read(iunit,*) ispin_read, off1_read, off2_read, mtxel_re, mtxel_im
       if(present(mtxel_real) .and. abs(mtxel_im) > TOL_Zero) call die("cannot have complex matrix elements in real version")
       do ispin = 1, nspin
          do ioff = 1, noffdiag
             if (spin_index(ispin) == ispin_read .and. offdiag1(ioff) == off1_read .and. offdiag2(ioff) == off2_read) then
                if(present(mtxel_real)) then
                   mtxel_real(ndiag + ioff, ispin) = mtxel_re
                else
                   mtxel_cplx(ndiag + ioff, ispin) = DCMPLX(mtxel_re, mtxel_im)
                endif
                found(ndiag + ioff, ispin) = .true.
             endif
          enddo
       enddo
    enddo

    if(any(.not. found(ndiag + 1:ndiag + noffdiag, :))) then
       if(peinf%inode == 0) then
          write(0,*) 'missing off-diagonal matrix elements (band, band, spin): '
          do ispin = 1, nspin
             do ioff = 1, noffdiag
                write(0,*) '(', offdiag1(ioff), ',', offdiag2(ioff), ',', spin_index(ispin), ')'
             enddo
          enddo
       endif
       call die("read_matrix_elements: not all needed data present")
    endif

    SAFE_DEALLOCATE(found)

    POP_SUB(read_matrix_elements_base)
    return
  end subroutine read_matrix_elements_base

  !=========================================================================
  !> this routine is used to read the ASCII files vxc.dat and x.dat, without typedefs types, with real output
  subroutine read_matrix_elements_real(iunit, iostat, kk, nspin, ndiag, noffdiag, spin_index, diag, offdiag1, offdiag2, mtxel)
    integer, intent(in) :: iunit !< file unit to read from
    integer, intent(out) :: iostat !< status of reading header, to detect end of file
    real(DP), intent(out) :: kk(3) !< kpoint read, in crystal coords
    integer, intent(in) :: nspin !< number of spins to read
    integer, intent(in) :: ndiag !< number of diagonal elements to read
    integer, intent(in) :: noffdiag !< number of offdiagonal elements to read
    integer, intent(in) :: spin_index(:) !< (nspin) mapping of 1:nspin to actual spins. 3 choices:
    !! spin-unpolarized: spin_index(1) = 1
    !! spin-polarized: spin_index(1) = 1, spin_index(2) = 2
    !! spin-polarized, spin 2 only: spin_index(1) = 1, spin_index(2) = 2
    integer, intent(in) :: diag(:) !< (ndiag) mapping of 1:ndiag onto band indices for diagonals
    integer, intent(in) :: offdiag1(:) !< (noffdiag) mapping of 1:noffdiag onto band indices for offdiagonal 1
    integer, intent(in) :: offdiag2(:) !< (noffdiag) mapping of 1:noffdiag onto band indices for offdiagonal 2
    real(DP), intent(out) :: mtxel(:,:) !< (ndiag+noffdiag,nspin)

    PUSH_SUB(read_matrix_elements_real)

    call read_matrix_elements_base(iunit, iostat, kk, nspin, ndiag, noffdiag, &
         spin_index, diag, offdiag1, offdiag2, mtxel_real = mtxel)

    POP_SUB(read_matrix_elements_real)
    return
  end subroutine read_matrix_elements_real

  !=========================================================================
  !> this routine is used to read the ASCII files vxc.dat and x.dat, without typedefs types, with cplx output
  subroutine read_matrix_elements_cplx(iunit, iostat, kk, nspin, ndiag, noffdiag, spin_index, diag, offdiag1, offdiag2, mtxel)
    integer, intent(in) :: iunit !< file unit to read from
    integer, intent(out) :: iostat !< status of reading header, to detect end of file
    real(DP), intent(out) :: kk(3) !< kpoint read, in crystal coords
    integer, intent(in) :: nspin !< number of spins to read
    integer, intent(in) :: ndiag !< number of diagonal elements to read
    integer, intent(in) :: noffdiag !< number of offdiagonal elements to read
    integer, intent(in) :: spin_index(:) !< (nspin) mapping of 1:nspin to actual spins. 3 choices:
    !! spin-unpolarized: spin_index(1) = 1
    !! spin-polarized: spin_index(1) = 1, spin_index(2) = 2
    !! spin-polarized, spin 2 only: spin_index(1) = 1, spin_index(2) = 2
    integer, intent(in) :: diag(:) !< (ndiag) mapping of 1:ndiag onto band indices for diagonals
    integer, intent(in) :: offdiag1(:) !< (noffdiag) mapping of 1:noffdiag onto band indices for offdiagonal 1
    integer, intent(in) :: offdiag2(:) !< (noffdiag) mapping of 1:noffdiag onto band indices for offdiagonal 2
    complex(DPC), intent(out) :: mtxel(:,:) !< (ndiag+noffdiag,nspin)

    PUSH_SUB(read_matrix_elements_cplx)

    call read_matrix_elements_base(iunit, iostat, kk, nspin, ndiag, noffdiag, &
         spin_index, diag, offdiag1, offdiag2, mtxel_cplx = mtxel)

    POP_SUB(read_matrix_elements_cplx)
    return
  end subroutine read_matrix_elements_cplx

  !=========================================================================
  !> this routine is used to read the ASCII files vxc.dat and x.dat, with typedefs types, with real output
  subroutine read_matrix_elements_type_real(iunit, iostat, kk, sig, mtxel)
    integer, intent(in) :: iunit !< file unit to read from
    integer, intent(out) :: iostat !< status of reading header, to detect end of file
    real(DP), intent(out) :: kk(3) !< kpoint read, in crystal coords
    type(siginfo), intent(in) :: sig !< structure containing other needed data
    real(DP), intent(out) :: mtxel(:,:) !< (ndiag+noffdiag,nspin)

    PUSH_SUB(read_matrix_elements_type_real)

    call read_matrix_elements(iunit, iostat, kk, sig%nspin, sig%ndiag, sig%noffdiag, &
         sig%spin_index, sig%diag, sig%off1, sig%off2, mtxel)

    POP_SUB(read_matrix_elements_type_real)
    return
  end subroutine read_matrix_elements_type_real

  !=========================================================================
  !> this routine is used to read the ASCII files vxc.dat and x.dat, with typedefs types, with cplx output
  ! ======
  ! Call from sigma_main.90:
  ! call read_matrix_elements_type(120, ierr, qk, sig, alda)
  ! ------
  ! kk = qk
  ! sig = sig
  ! mtxel = alda
  ! -------
  subroutine read_matrix_elements_type_cplx(iunit, iostat, kk, sig, mtxel)
    integer, intent(in) :: iunit !< file unit to read from
    integer, intent(out) :: iostat !< status of reading header, to detect end of file
    real(DP), intent(out) :: kk(3) !< kpoint read, in crystal coords
    type(siginfo), intent(in) :: sig !< structure containing other needed data
    complex(DPC), intent(out) :: mtxel(:,:) !< (ndiag+noffdiag,nspin)

    PUSH_SUB(read_matrix_elements_type_cplx)

    call read_matrix_elements(iunit, iostat, kk, sig%nspin, sig%ndiag, sig%noffdiag, &
         sig%spin_index, sig%diag, sig%off1, sig%off2, mtxel)

    POP_SUB(read_matrix_elements_type_cplx)
    return
  end subroutine read_matrix_elements_type_cplx

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

    PUSH_SUB(dealloc_header)

    if (sheader .eq. 'WFN') then
       wfnflag = .true.
    elseif (sheader .eq. 'RHO' .or. sheader .eq. 'VXC') then
       wfnflag = .false.
    else
       call die("unknown file header: '" + sheader + "' (should be 'WFN'/'RHO'/'VXC')")
    endif

    SAFE_DEALLOCATE_P(atyp)
    SAFE_DEALLOCATE_P(apos)
    if (wfnflag) then
       SAFE_DEALLOCATE_P(ngk)
       SAFE_DEALLOCATE_P(kw)
       SAFE_DEALLOCATE_P(kpt)
       SAFE_DEALLOCATE_P(ifmin)
       SAFE_DEALLOCATE_P(ifmax)
       SAFE_DEALLOCATE_P(energies)
       SAFE_DEALLOCATE_P(occupations)
    endif

    POP_SUB(dealloc_header)
    return
  end subroutine dealloc_header

  !=========================================================================
  subroutine dealloc_header_type(sheader, crys, kp)
    character(len=3), intent(in) :: sheader
    type(crystal), intent(inout) :: crys
    type(kpoints), intent(inout) :: kp

    PUSH_SUB(dealloc_header_type)

    call dealloc_header(sheader, crys%atyp, crys%apos, kp%ngk, kp%w, kp%rk, kp%ifmin, kp%ifmax, kp%el, kp%occ)

    POP_SUB(dealloc_header_type)
    return
  end subroutine dealloc_header_type

  !=========================================================================
  subroutine dealloc_kp(kp)
    type(kpoints), intent(inout) :: kp

    PUSH_SUB(dealloc_kp)

    SAFE_DEALLOCATE_P(kp%ngk)
    SAFE_DEALLOCATE_P(kp%w)
    SAFE_DEALLOCATE_P(kp%rk)
    SAFE_DEALLOCATE_P(kp%ifmin)
    SAFE_DEALLOCATE_P(kp%ifmax)
    SAFE_DEALLOCATE_P(kp%el)
    SAFE_DEALLOCATE_P(kp%occ)

    POP_SUB(dealloc_kp)
    return
  end subroutine dealloc_kp

  !=========================================================================
  subroutine dealloc_crys(crys)
    type(crystal), intent(inout) :: crys

    PUSH_SUB(dealloc_crys)

    SAFE_DEALLOCATE_P(crys%atyp)
    SAFE_DEALLOCATE_P(crys%apos)

    POP_SUB(dealloc_crys)
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

    PUSH_SUB(check_header)

    tolerant_ = .false.
    if(present(tolerant)) tolerant_ = tolerant
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
    
    POP_SUB(check_header)
    return
  end subroutine check_header

  !> require `version` to be the same as `version_ref`
  subroutine require_version(fname, version, version_ref)
    character(len=*), intent(in) :: fname
    integer, intent(in) :: version
    integer, intent(in) :: version_ref

    PUSH_SUB(require_version)

    if (version/=version_ref) then
       if (peinf%inode==0) then
          write(0,*) 'ERROR: Wrong version for file "',TRUNC(fname),'".'
          write(0,*) 'Expected: ', version_ref
          write(0,*) 'Got: ', version
       endif
       call die('Wrong version for file "'+TRUNC(fname)+'".', only_root_writes=.true.)
    endif

    POP_SUB(require_version)

  end subroutine require_version

  !> A high-level wrapper for write_*_header* functions
  subroutine write_mf_header(iunit, mf)
    integer, intent(in) :: iunit
    type(mf_header_t), intent(in) :: mf

    character(len=3) :: sheader
    character(len=16) :: fmt_str
    logical :: is_fmt = .false.

    PUSH_SUB(write_mf_header)

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

    POP_SUB(write_mf_header)

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

    PUSH_SUB(read_mf_header)

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

    POP_SUB(read_mf_header)

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

    PUSH_SUB(init_mf_header_from_types)

    mf_header%version = version
    mf_header%sheader = sheader
    mf_header%sdate = ''
    mf_header%stime = ''
    mf_header%iflavor = iflavor
    mf_header%kp = kp
    mf_header%gvec = gvec
    mf_header%syms = syms
    mf_header%crys = crys

    POP_SUB(init_mf_header_from_types)

  end subroutine init_mf_header_from_types

  subroutine prepare_syms(crys, syms)

    type(crystal), intent(in) :: crys
    type(symmetry), intent(inout) :: syms
    real(DP), dimension(3,3) :: atTat
    real(DP), dimension(3,3) :: inv_atTat
    real(DP), dimension(3,3) :: inv_at
    integer :: itran
    logical :: flag_inv
    
    PUSH_SUB(prepare_syms)

    ! crys%avec(:,:) = at(:,:)
    ! at(:,:) = [a1, a2, a3]
    atTat = MATMUL(TRANSPOSE(crys%avec),crys%avec)

    ! -------------------------------
    !                               [a1]
    ! ------adot = [a1 a2 a3] \cdot [a2]
    !                               [a3]
    ! -------------------------------
    
    call M33INV(atTat, inv_atTat, flag_inv)
    if (.not. flag_inv) then
       call die("prepare_syms: atTat not invertable", only_root_writes=.true.)
    endif

    call M33INV(crys%avec,inv_at, flag_inv)
    if (.not. flag_inv) then
       call die("prepare_syms: crys%avec not invertable", only_root_writes=.true.)
    endif

    syms%mtrx_reci(:,:,:) = 0.0D0
    syms%mtrx_cart(:,:,:) = 0.0D0

    do itran = 1, syms%ntran
       syms%mtrx_reci(1:3,1:3,itran) = nint(MATMUL(MATMUL(atTat, dble(syms%mtrx(1:3,1:3,itran))), inv_atTat))       
       syms%mtrx_cart(1:3,1:3,itran) = MATMUL(MATMUL(crys%avec, dble(syms%mtrx(1:3,1:3,itran))), inv_at)
    enddo
    
    POP_SUB(prepare_syms)

  end subroutine prepare_syms

end module wfn_rho_vxc_io_m

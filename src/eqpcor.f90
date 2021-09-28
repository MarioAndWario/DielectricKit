#include "f_defs.h"

!================================================================================
!
! Module eqpcor_m
!
! 1. eqpcor()        Originally By gsm       Last Modified 8/4/2015 (FHJ)
!
!    Reads quasiparticle energy corrections from eqp.dat-type files.
!    Such files are made by bin/eqp.py.
!
!    There are some differences in how eqpcor is called from different codes.
!
!    In Epsilon, eqpcor is called from proc 0, so set inode = 0 and npes = 1.
!    In other places, set inode = peinf%inode and npes = peinf%npes.
!
!    Epsilon and BSE require quasiparticle energies in Rydbergs, so set
!    irydflag = 1. Sigma (outer) requires energies in eVs, so set irydflag = 0.
!
!    In Epsilon and Sigma, the quasiparticle energies are returned in
!    array eqp(nbmin:nbmax,1:kp%nrk,1:kp%nspin). In BSE, the valence and
!    conduction energies indexed with respect to the Fermi level are
!    returned in arrays eqpv(:,1:kp%nrk,1:kp%nspin) and eqpc(:,1:kp%nrk,1:kp%nspin).
!
!    In inteqp, set ivalflag = 2 to return the difference Eqp - Edft.
!    In other places, set ivalflag = 0 to return Eqp.
!
!    DO NOT try to use both styles, they will overwrite each other!
!
!================================================================================

module eqpcor_m
  use global_m
  use misc_m
  implicit none
  private
  public :: eqpcor

contains

  !> use ROOT to update kp%el using fn = 'eqp.dat
  !> and then bcast to all procs
  !> kp%el are in units of Ryd

  subroutine eqpcor(crys, kp, ib_min, ib_max, file_eqp)
    type (crystal), intent(in) :: crys
    type(kpoints), intent(inout) :: kp
    integer, intent(in) :: ib_min, ib_max
    character(LEN=*), intent(in) :: file_eqp
    integer :: eof, irk, irk_target, is, ib, ib_eqp, nb_eqp, nk_found, gumk(3)
    !> e_dft is the meanfield energy read from eqp.dat
    !> e_qp is the quasi-particle energy read from eqp.dat
    real(DP) :: e_dft, e_qp, k_eqp(3), delta_k(3)
    character(100) :: errmsg
    ! real(DP), parameter :: TOL_eqp = 1D-5
    real(DP), parameter :: TOL_eqp = 1D-3
    PUSH_SUB(eqpcor)

    if (peinf%inode .eq. 0) then
       write(6,'(1X,A)') "Reading quasiparticle energy corrections from "//trim(file_eqp)
       call open_file(9, file=TRUNC(file_eqp), form='formatted', status='old')

       nk_found = 0
       !> number of kpoint in eqp.dat could be larger than kp%nrk
       do while (nk_found < kp%nrk)
          nb_eqp = 0
          read(9, *, iostat=eof) k_eqp(:), nb_eqp

          !> eof /= 0 means we reach the end of eqp file
          if (eof .ne. 0) then
             call die("Missing kpoints in file"//TRUNC(file_eqp), only_root_writes=.true.)
          endif

          irk_target = 0
          rk_loop: do irk = 1, kp%nrk
             delta_k(:) = k_eqp(:) - kp%rk(:,irk)
             call get_gumk3(crys%bdot, delta_k, gumk)
             delta_k(:) = delta_k(:) - DBLE(gumk(:))
             if (NORM2(delta_k(:)) < TOL_SMALL ) then
                irk_target = irk
                nk_found = nk_found + 1
                exit rk_loop
             endif
          enddo rk_loop

          do ib_eqp = 1, nb_eqp
             read(9, *, iostat=eof) is, ib, e_dft, e_qp
             if (eof .ne. 0) then
                write(errmsg,'(a)') 'Wrong contents of k-point blocks in file '//TRUNC(file_eqp)
                call die(errmsg, only_root_writes=.true.)
             endif
             if (irk_target .eq. 0) cycle
             if ((ib .ge. ib_min) .and. (ib .le. ib_max)) then
                !> eqp = kp%el initially stores DFT energies in units of Ryd
                if (ABS(e_dft - kp%elda(ib, irk_target, is) * RYD) > TOL_eqp) then
                   write(*,'(A,I5,A,I5,A,I5,A,F15.8,A,F15.8)') "ib = ", ib, " irk = ", irk_target, " is = ", is, " e_dft = ", e_dft, " kpela = ", kp%elda(ib, irk_target, is)
                   call die("eqpcor: eqpcor mean-field energy mismatch")
                endif

                kp%el(ib, irk_target, is) = e_qp/RYD
             endif
          enddo !> ib_eqp
       enddo
       call close_file(9)
    endif

#ifdef MPI
    if (peinf%npes > 1) then
       call MPI_Bcast(kp%el(1,1,1), size(kp%el), MPI_REAL_DP, 0, MPI_COMM_WORLD, mpierr)
    endif
#endif

    POP_SUB(eqpcor)
    return
  end subroutine eqpcor
end module eqpcor_m

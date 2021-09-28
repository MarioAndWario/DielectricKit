#include "f_defs.h"

!==============================================================================
!
! Routines()
!
! (1) write_eps()    By Meng Wu (2020)
!
!     Writes out polarizability or dielectric response function in real space
!     in the .xsf format (http://www.xcrysden.org/doc/XSF.html).
!
!     Note that .xsf requires that the first and last indices correspond to the
!     same point, even for a periodic crystal structure.
!
!=============================================================================

module write_eps_m
  use global_m
  use realspace_common_m
  implicit none
  private
  public :: write_eps

contains

  subroutine write_eps(peps, crys, nfft, scfft, ir2)
    type(realspace_t), intent(in) :: peps
    type(crystal), intent(in) :: crys
    integer, intent(in) :: nfft(3)
    complex(DP), intent(in) :: scfft(:,:,:) !< (nw_eff(1),nw_eff(2),nw_eff(3))
    integer, intent(in) :: ir2
    integer :: i1,i2,i3, i1_, i2_, i3_
    integer :: nw(3), nw_eff(3), natom_tot, count, isuper1, isuper2, isuper3, iatom
    character :: filename*100
    real(DP) :: r2_cart_ang(3), shift_cart_ang(3)
    PUSH_SUB(write_eps)

    nw(:) = nfft(:) * peps%nsuper(:)
    nw_eff(:) = (nw(:) + peps%downsample(:) - 1) / peps%downsample(:)
    r2_cart_ang(:) = MATMUL(crys%avec(:,:), peps%r2(:,ir2)) * crys%alat * BOHR ! in angstrom

    !! output xsf file
    !! Use wrapper functions to output
    
    !! real part of epsinv
    if (peps%realpart) then
       write(filename,'(a,i0,a,i0,a)') 'OUT.REAL.ir2_', ir2, '_is_', peps%ispin,'.xsf'
       write(6,'(A)')
       write(6,'(1X,A)') "Writing data set to "//TRUNC(filename)
       call open_file(30,file=TRUNC(filename),form='formatted',status='replace')

       write(30, '(A,3f15.8,A,2f15.8)') "# r2 = ", peps%r2(:,ir2), " freq = ", peps%freq_target
       write(30, '(A)') "CRYSTAL"
       write(30, '(A)') "PRIMVEC"
       write(30, '(3f15.8)') peps%nsuper(1) * crys%avec(1:3,1) * crys%alat * BOHR ! in Angstrom
       write(30, '(3f15.8)') peps%nsuper(2) * crys%avec(1:3,2) * crys%alat * BOHR ! in Angstrom
       write(30, '(3f15.8)') peps%nsuper(3) * crys%avec(1:3,3) * crys%alat * BOHR ! in Angstrom
       write(30, '(A)') "CONVVEC"
       write(30, '(3f15.8)') peps%nsuper(1) * crys%avec(1:3,1) * crys%alat * BOHR ! in Angstrom
       write(30, '(3f15.8)') peps%nsuper(2) * crys%avec(1:3,2) * crys%alat * BOHR ! in Angstrom
       write(30, '(3f15.8)') peps%nsuper(3) * crys%avec(1:3,3) * crys%alat * BOHR ! in Angstrom

       write(30, '(A)') "PRIMCOORD"
       natom_tot=peps%nsuper(1)*peps%nsuper(2)*peps%nsuper(3)*crys%nat
       write(30, '(I5,1X,I5)') natom_tot+1, 1

       !! Print all the atomic positions, considering supercell
       do isuper3 = 1, peps%nsuper(3)
          do isuper2 = 1, peps%nsuper(2)
             do isuper1 = 1, peps%nsuper(1)
                !! crys%apos(1:3,1:crys%nat) in cartesian bohr units
                do iatom = 1, crys%nat
                   shift_cart_ang(1:3) = ( dble(isuper1-1) * crys%avec(1:3,1) * crys%alat + dble(isuper2-1) * crys%avec(1:3,2) * crys%alat + dble(isuper3-1) * crys%avec(1:3,3) * crys%alat ) * BOHR

                   write(30,'(I5,1X,3F15.8)') crys%atyp(iatom), crys%apos(1:3,iatom) * crys%alat * BOHR + shift_cart_ang(1:3) ! in Angstrom
                enddo
             enddo
          enddo
       enddo

       !! Print the position of r2
       write(30, '(I5,1X,3F15.8)') 1, r2_cart_ang(1:3) ! r2

       write(30, '(A)') "BEGIN_BLOCK_DATAGRID_3D"
       write(30, '(A)') "  real_epsinv"
       write(30, '(A)') "  BEGIN_DATAGRID_3D_real_epsinv"

       !! Output number of data-points in each direction (i.e. nx ny nz for 3D grids), note that xsf use general grid (with dim nw_eff(:)+1)
       write(30, '(3I5)') nw_eff(:) + 1
       !! origin of the datagrid
       write(30, '(3f15.8)') 0.0D0, 0.0D0, 0.0D0
       !! spanning vector of the datagrid
       write(30, '(3f15.8)') peps%nsuper(1) * crys%avec(1:3,1) * crys%alat * BOHR ! in Angstrom
       write(30, '(3f15.8)') peps%nsuper(2) * crys%avec(1:3,2) * crys%alat * BOHR ! in Angstrom
       write(30, '(3f15.8)') peps%nsuper(3) * crys%avec(1:3,3) * crys%alat * BOHR ! in Angstrom
       do i3_ = 1, nw_eff(3)+1
          if (i3_ .eq. nw_eff(3)+1) then
             i3 = 1
          else
             i3 = i3_
          endif
          count = 0
          do i2_ = 1, nw_eff(2)+1
             if (i2_ .eq. nw_eff(2)+1) then
                i2 = 1
             else
                i2 = i2_
             endif
             do i1_ = 1, nw_eff(1)+1
                if (i1_ .eq. nw_eff(1)+1) then
                   i1 = 1
                else
                   i1 = i1_
                endif
                write(30, '(es13.5,2X)', advance="no") DBLE(scfft(i1, i2, i3))
                count = count + 1
                if (mod(count,5) .eq. 0) then
                   write(30,'(A)')
                endif
             enddo
          enddo
          if (i3_ .ne. nw_eff(3)+1) then
             write(30,'(A)') NEW_LINE('A')
          endif
       enddo
       write(30, '(A)') "  END_DATAGRID_3D"
       write(30, '(A)') "END_BLOCK_DATAGRID_3D"
       call close_file(30)
    endif

    !! imaginary part of epsinv
    if (peps%imagpart) then
       write(filename,'(a,i0,a,i0,a)') 'OUT.IMAG.ir2_', ir2, '_is_', peps%ispin,'.xsf'
       write(6,'(A)')       
       write(6,'(1X,A)') "Writing data set to "//TRUNC(filename)
       call open_file(30,file=TRUNC(filename),form='formatted',status='replace')

       write(30, '(A,3f15.8,A,2f15.8)') "# r2 = ", peps%r2(:,ir2), " freq = ", peps%freq_target
       write(30, '(A)') "CRYSTAL"
       write(30, '(A)') "PRIMVEC"
       write(30, '(3f15.8)') peps%nsuper(1) * crys%avec(1:3,1) * crys%alat * BOHR ! in Angstrom
       write(30, '(3f15.8)') peps%nsuper(2) * crys%avec(1:3,2) * crys%alat * BOHR ! in Angstrom
       write(30, '(3f15.8)') peps%nsuper(3) * crys%avec(1:3,3) * crys%alat * BOHR ! in Angstrom
       write(30, '(A)') "CONVVEC"
       write(30, '(3f15.8)') peps%nsuper(1) * crys%avec(1:3,1) * crys%alat * BOHR ! in Angstrom
       write(30, '(3f15.8)') peps%nsuper(2) * crys%avec(1:3,2) * crys%alat * BOHR ! in Angstrom
       write(30, '(3f15.8)') peps%nsuper(3) * crys%avec(1:3,3) * crys%alat * BOHR ! in Angstrom

       write(30, '(A)') "PRIMCOORD"
       natom_tot=peps%nsuper(1)*peps%nsuper(2)*peps%nsuper(3)*crys%nat
       write(30, '(I5,1X,I5)') natom_tot+1, 1

       !! Print all the atomic positions, considering supercell
       do isuper3 = 1, peps%nsuper(3)
          do isuper2 = 1, peps%nsuper(2)
             do isuper1 = 1, peps%nsuper(1)
                !! crys%apos(1:3,1:crys%nat) in cartesian bohr units
                do iatom = 1, crys%nat
                   shift_cart_ang(1:3) = ( dble(isuper1-1) * crys%avec(1:3,1) * crys%alat + dble(isuper2-1) * crys%avec(1:3,2) * crys%alat + dble(isuper3-1) * crys%avec(1:3,3) * crys%alat ) * BOHR
                   write(30,'(I5,1X,3F15.8)') crys%atyp(iatom), crys%apos(1:3,iatom) * crys%alat * BOHR + shift_cart_ang(1:3) ! in Angstrom
                enddo
             enddo
          enddo
       enddo

       !! Print the position of r2
       write(30, '(I5,1X,3F15.8)') 1, r2_cart_ang(1:3) ! r2
       write(30, '(A)') "BEGIN_BLOCK_DATAGRID_3D"
       write(30, '(A)') "  imag_epsinv"
       write(30, '(A)') "  BEGIN_DATAGRID_3D_imag_epsinv"
       !! Output number of data-points in each direction (i.e. nx ny nz for 3D grids), note that xsf use general grid (with dim nw_eff(:)+1)
       write(30, '(3I5)') nw_eff(:) + 1
       !! origin of the datagrid
       write(30, '(3f15.8)') 0.0D0, 0.0D0, 0.0D0
       !! spanning vector of the datagrid
       write(30, '(3f15.8)') peps%nsuper(1) * crys%avec(1:3,1) * crys%alat * BOHR ! in Angstrom
       write(30, '(3f15.8)') peps%nsuper(2) * crys%avec(1:3,2) * crys%alat * BOHR ! in Angstrom
       write(30, '(3f15.8)') peps%nsuper(3) * crys%avec(1:3,3) * crys%alat * BOHR ! in Angstrom
       do i3_ = 1, nw_eff(3)+1
          if (i3_ .eq. nw_eff(3)+1) then
             i3 = 1
          else
             i3 = i3_
          endif
          count = 0
          do i2_ = 1, nw_eff(2)+1
             if (i2_ .eq. nw_eff(2)+1) then
                i2 = 1
             else
                i2 = i2_
             endif
             do i1_ = 1, nw_eff(1)+1
                if (i1_ .eq. nw_eff(1)+1) then
                   i1 = 1
                else
                   i1 = i1_
                endif
                write(30, '(es13.5,2X)', advance="no") DIMAG(scfft(i1, i2, i3))
                count = count + 1
                if (mod(count,5) .eq. 0) then
                   write(30,'(A)')
                endif
             enddo
          enddo
          if (i3_ .ne. nw_eff(3)+1) then
             write(30,'(A)') NEW_LINE('A')
          endif
       enddo
       write(30, '(A)') "  END_DATAGRID_3D"
       write(30, '(A)') "END_BLOCK_DATAGRID_3D"
       call close_file(30)
    endif

    !! abs of epsinv
    if (peps%abs) then
       write(filename,'(a,i0,a,i0,a)') 'OUT.ABS.ir2_', ir2, '_is_', peps%ispin,'.xsf'
       write(6,'(A)')       
       write(6,'(1X,A)') "Writing data set to "//TRUNC(filename)
       call open_file(30,file=TRUNC(filename),form='formatted',status='replace')

       write(30, '(A,3f15.8,A,2f15.8)') "# r2 = ", peps%r2(:,ir2), " freq = ", peps%freq_target
       write(30, '(A)') "CRYSTAL"
       write(30, '(A)') "PRIMVEC"
       write(30, '(3f15.8)') peps%nsuper(1) * crys%avec(1:3,1) * crys%alat * BOHR ! in Angstrom
       write(30, '(3f15.8)') peps%nsuper(2) * crys%avec(1:3,2) * crys%alat * BOHR ! in Angstrom
       write(30, '(3f15.8)') peps%nsuper(3) * crys%avec(1:3,3) * crys%alat * BOHR ! in Angstrom
       write(30, '(A)') "CONVVEC"
       write(30, '(3f15.8)') peps%nsuper(1) * crys%avec(1:3,1) * crys%alat * BOHR ! in Angstrom
       write(30, '(3f15.8)') peps%nsuper(2) * crys%avec(1:3,2) * crys%alat * BOHR ! in Angstrom
       write(30, '(3f15.8)') peps%nsuper(3) * crys%avec(1:3,3) * crys%alat * BOHR ! in Angstrom

       write(30, '(A)') "PRIMCOORD"
       natom_tot=peps%nsuper(1)*peps%nsuper(2)*peps%nsuper(3)*crys%nat
       write(30, '(I5,1X,I5)') natom_tot+1, 1

       !! Print all the atomic positions, considering supercell
       do isuper3 = 1, peps%nsuper(3)
          do isuper2 = 1, peps%nsuper(2)
             do isuper1 = 1, peps%nsuper(1)
                !! crys%apos(1:3,1:crys%nat) in cartesian bohr units
                do iatom = 1, crys%nat
                   shift_cart_ang(1:3) = ( dble(isuper1-1) * crys%avec(1:3,1) * crys%alat + dble(isuper2-1) * crys%avec(1:3,2) * crys%alat + dble(isuper3-1) * crys%avec(1:3,3) * crys%alat ) * BOHR

                   write(30,'(I5,1X,3F15.8)') crys%atyp(iatom), crys%apos(1:3,iatom) * crys%alat * BOHR + shift_cart_ang(1:3) ! in Angstrom
                enddo
             enddo
          enddo
       enddo

       !! Print the position of r2
       write(30, '(I5,1X,3F15.8)') 1, r2_cart_ang(1:3) ! r2

       write(30, '(A)') "BEGIN_BLOCK_DATAGRID_3D"
       write(30, '(A)') "  abs_epsinv"
       write(30, '(A)') "  BEGIN_DATAGRID_3D_abs_epsinv"

       !! Output number of data-points in each direction (i.e. nx ny nz for 3D grids), note that xsf use general grid (with dim nw_eff(:)+1)
       write(30, '(3I5)') nw_eff(:) + 1
       !! origin of the datagrid
       write(30, '(3f15.8)') 0.0D0, 0.0D0, 0.0D0
       !! spanning vector of the datagrid
       write(30, '(3f15.8)') peps%nsuper(1) * crys%avec(1:3,1) * crys%alat * BOHR ! in Angstrom
       write(30, '(3f15.8)') peps%nsuper(2) * crys%avec(1:3,2) * crys%alat * BOHR ! in Angstrom
       write(30, '(3f15.8)') peps%nsuper(3) * crys%avec(1:3,3) * crys%alat * BOHR ! in Angstrom
       do i3_ = 1, nw_eff(3)+1
          if (i3_ .eq. nw_eff(3)+1) then
             i3 = 1
          else
             i3 = i3_
          endif
          count = 0
          do i2_ = 1, nw_eff(2)+1
             if (i2_ .eq. nw_eff(2)+1) then
                i2 = 1
             else
                i2 = i2_
             endif
             do i1_ = 1, nw_eff(1)+1
                if (i1_ .eq. nw_eff(1)+1) then
                   i1 = 1
                else
                   i1 = i1_
                endif
                write(30, '(es13.5,2X)', advance="no") ABS(scfft(i1, i2, i3))
                count = count + 1
                if (mod(count,5) .eq. 0) then
                   write(30,'(A)')
                endif
             enddo
          enddo
          if (i3_ .ne. nw_eff(3)+1) then
             write(30,'(A)') NEW_LINE('A')
          endif
       enddo
       write(30, '(A)') "  END_DATAGRID_3D"
       write(30, '(A)') "END_BLOCK_DATAGRID_3D"
       call close_file(30)
    endif

    !! abs2 of epsinv
    if (peps%abs2) then
       write(filename,'(a,i0,a,i0,a)') 'OUT.ABS2.ir2_', ir2, '_is_', peps%ispin,'.xsf'
       write(6,'(A)')       
       write(6,'(1X,A)') "Writing data set to "//TRUNC(filename)
       call open_file(30,file=TRUNC(filename),form='formatted',status='replace')

       write(30, '(A,3f15.8,A,2f15.8)') "# r2 = ", peps%r2(:,ir2), " freq = ", peps%freq_target
       write(30, '(A)') "CRYSTAL"
       write(30, '(A)') "PRIMVEC"
       write(30, '(3f15.8)') peps%nsuper(1) * crys%avec(1:3,1) * crys%alat * BOHR ! in Angstrom
       write(30, '(3f15.8)') peps%nsuper(2) * crys%avec(1:3,2) * crys%alat * BOHR ! in Angstrom
       write(30, '(3f15.8)') peps%nsuper(3) * crys%avec(1:3,3) * crys%alat * BOHR ! in Angstrom
       write(30, '(A)') "CONVVEC"
       write(30, '(3f15.8)') peps%nsuper(1) * crys%avec(1:3,1) * crys%alat * BOHR ! in Angstrom
       write(30, '(3f15.8)') peps%nsuper(2) * crys%avec(1:3,2) * crys%alat * BOHR ! in Angstrom
       write(30, '(3f15.8)') peps%nsuper(3) * crys%avec(1:3,3) * crys%alat * BOHR ! in Angstrom

       write(30, '(A)') "PRIMCOORD"
       natom_tot=peps%nsuper(1)*peps%nsuper(2)*peps%nsuper(3)*crys%nat
       write(30, '(I5,1X,I5)') natom_tot+1, 1

       !! Print all the atomic positions, considering supercell
       do isuper3 = 1, peps%nsuper(3)
          do isuper2 = 1, peps%nsuper(2)
             do isuper1 = 1, peps%nsuper(1)
                !! crys%apos(1:3,1:crys%nat) in cartesian bohr units
                do iatom = 1, crys%nat
                   shift_cart_ang(1:3) = ( dble(isuper1-1) * crys%avec(1:3,1) * crys%alat + dble(isuper2-1) * crys%avec(1:3,2) * crys%alat + dble(isuper3-1) * crys%avec(1:3,3) * crys%alat ) * BOHR
                   write(30,'(I5,1X,3F15.8)') crys%atyp(iatom), crys%apos(1:3,iatom) * crys%alat * BOHR + shift_cart_ang(1:3) ! in Angstrom
                enddo
             enddo
          enddo
       enddo

       !! Print the position of r2
       write(30, '(I5,1X,3F15.8)') 1, r2_cart_ang(1:3) ! r2
       write(30, '(A)') "BEGIN_BLOCK_DATAGRID_3D"
       write(30, '(A)') "  abs2_epsinv"
       write(30, '(A)') "  BEGIN_DATAGRID_3D_abs2_epsinv"

       !! Output number of data-points in each direction (i.e. nx ny nz for 3D grids), note that xsf use general grid (with dim nw_eff(:)+1)
       write(30, '(3I5)') nw_eff(:) + 1
       !! origin of the datagrid
       write(30, '(3f15.8)') 0.0D0, 0.0D0, 0.0D0
       !! spanning vector of the datagrid
       write(30, '(3f15.8)') peps%nsuper(1) * crys%avec(1:3,1) * crys%alat * BOHR ! in Angstrom
       write(30, '(3f15.8)') peps%nsuper(2) * crys%avec(1:3,2) * crys%alat * BOHR ! in Angstrom
       write(30, '(3f15.8)') peps%nsuper(3) * crys%avec(1:3,3) * crys%alat * BOHR ! in Angstrom
       do i3_ = 1, nw_eff(3)+1
          if (i3_ .eq. nw_eff(3)+1) then
             i3 = 1
          else
             i3 = i3_
          endif
          count = 0
          do i2_ = 1, nw_eff(2)+1
             if (i2_ .eq. nw_eff(2)+1) then
                i2 = 1
             else
                i2 = i2_
             endif
             do i1_ = 1, nw_eff(1)+1
                if (i1_ .eq. nw_eff(1)+1) then
                   i1 = 1
                else
                   i1 = i1_
                endif
                write(30, '(es13.5,2X)', advance="no") ABS(scfft(i1, i2, i3))**2
                count = count + 1
                if (mod(count,5) .eq. 0) then
                   write(30,'(A)')
                endif
             enddo
          enddo
          if (i3_ .ne. nw_eff(3)+1) then
             write(30,'(A)') NEW_LINE('A')
          endif
       enddo
       write(30, '(A)') "  END_DATAGRID_3D"
       write(30, '(A)') "END_BLOCK_DATAGRID_3D"
       call close_file(30)
    endif       
    POP_SUB(write_eps)
    return
  end subroutine write_eps

end module write_eps_m

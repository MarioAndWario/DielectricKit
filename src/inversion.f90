!==============================================================================
!
! Routines:
!
! (1) Xinvert_with_scalapack_d() Originally by JRD Last Modified 02/2015 (FHJ)
!
! This routine inverts a matrix which is already distributed in block
! cyclic form with ScaLAPACK.
!
! (2) Xinvert_serial() Originally by JRD Last Modified 02/2015 (FHJ)
!
! Inverts a matrix using LAPACK.
!
!==============================================================================
module inversion_m
  use global_m
  use lapack_m
  use scalapack_m
  implicit none
  private
  public :: &
       dinvert_with_scalapack, &
       zinvert_with_scalapack, &
       dinvert_serial, &
       zinvert_serial
contains

  !---------------------- Use scaLAPACK For Inversion -----------------------------------
  subroutine dinvert_with_scalapack(nmtx, scal, matrix)
    integer, intent(in) :: nmtx
    type (scalapack), intent(in) :: scal
    REAL(DP), intent(inout) :: matrix(scal%npr,scal%npc)
    integer :: info, desca(9), lwork, liwork, ipiv(scal%npr+scal%nbl)
    integer, allocatable :: iwork(:)
    REAL(DP), allocatable :: work(:)

    if (peinf%verb_debug .and. peinf%inode==0) then
       write(6,*)
       write(6,*) 'We are doing inversion with ScaLAPACK.'
       write(6,*)
       write(6,*) 'Peinf: ', peinf%inode,' Scalapack Grid:'
       write(6,*) scal%nprow,scal%npcol,scal%myprow,scal%mypcol
       write(6,*) nmtx,scal%nbl,scal%npr,scal%npc
       write(6,*)
    endif
    call descinit(desca, nmtx, nmtx, scal%nbl, scal%nbl, 0, 0, &
         scal%icntxt, max(1,scal%npr), info)
    if(info < 0) then
       write(0,'(a,i3,a)') 'Argument number ', -info, ' had an illegal value on entry.'
       call die("descinit error for descaA in inversion")
    else if(info > 0) then
       write(0,*) 'info = ', info
       call die("descinit error for descaA in inversion")
    endif
    ! FHJ: LU factorization of the matrix
    call pdgetrf(nmtx, nmtx, matrix, 1, 1, desca, ipiv, info)
    if (info/=0) then
       if (peinf%inode==0) write(0,*) 'ERROR: got info = ', info, ' in p?getrf'
       call die('p?getrf failed')
    endif
    ! FHJ: tringular inversion of LU decomposition
    allocate(work (10))
    allocate(iwork (10))
    call pdgetri(nmtx, matrix, 1, 1, desca, ipiv, work, -1, iwork, -1, info)
    if (info/=0) then
       if (peinf%inode==0) write(0,*) 'ERROR: got info = ', info, ' in p?getri'
       call die('p?getri failed for query mode')
    endif
    lwork = max(1,int(work(1)))
    liwork = max(1,iwork(1))
    if (allocated(work)) then
       deallocate(work)
    endif
    if (allocated(iwork)) then
       deallocate(iwork)
    endif
    allocate(work (lwork))
    allocate(iwork (liwork))
    call pdgetri(nmtx, matrix, 1, 1, desca, ipiv, work, lwork, iwork, liwork, info)
    if (info/=0) then
       if (peinf%inode==0) write(0,*) 'ERROR: got info = ', info, ' in p?getri'
       call die('p?getri failed')
    endif
    if (allocated(iwork)) then
       deallocate(iwork)
    endif
    if (allocated(work)) then
       deallocate(work)
    endif

    return
  end subroutine dinvert_with_scalapack
  !------------------------------------------------------------
  subroutine dinvert_serial(nmtx, matrix)
    integer, intent(in) :: nmtx
    REAL(DP), intent(inout) :: matrix(nmtx,nmtx)
    integer :: info, lwork, ipiv(nmtx)
    REAL(DP), allocatable :: work(:)

    ! FHJ: LU factorization of the matrix
    call dgetrf(nmtx, nmtx, matrix, nmtx, ipiv, info)
    if (info/=0) then
       if (peinf%inode==0) write(0,*) 'ERROR: got info = ', info, ' in ?getrf'
       call die('?getrf failed')
    endif
    ! FHJ: tringular inversion of LU decomposition
    allocate(work (10))
    call dgetri(nmtx, matrix, nmtx, ipiv, work, -1, info)
    if (info/=0) then
       if (peinf%inode==0) write(0,*) 'ERROR: got info = ', info, ' in ?getri'
       call die('?getri failed for query mode')
    endif
    lwork = max(1,int(work(1)))
    if (allocated(work)) then
       deallocate(work)
    endif
    allocate(work (lwork))
    call dgetri(nmtx, matrix, nmtx, ipiv, work, lwork, info)
    if (info/=0) then
       if (peinf%inode==0) write(0,*) 'ERROR: got info = ', info, ' in ?getri'
       call die('?getri failed')
    endif
    if (allocated(work)) then
       deallocate(work)
    endif

    return
  end subroutine dinvert_serial

  !---------------------- Use scaLAPACK For Inversion -----------------------------------
  subroutine zinvert_with_scalapack(nmtx, scal, matrix)
    integer, intent(in) :: nmtx
    type (scalapack), intent(in) :: scal
    COMPLEX(DPC), intent(inout) :: matrix(scal%npr,scal%npc)
    integer :: info, desca(9), lwork, liwork, ipiv(scal%npr+scal%nbl)
    integer, allocatable :: iwork(:)
    COMPLEX(DPC), allocatable :: work(:)

    if (peinf%verb_debug .and. peinf%inode==0) then
       write(6,*)
       write(6,*) 'We are doing inversion with ScaLAPACK.'
       write(6,*)
       write(6,*) 'Peinf: ', peinf%inode,' Scalapack Grid:'
       write(6,*) scal%nprow,scal%npcol,scal%myprow,scal%mypcol
       write(6,*) nmtx,scal%nbl,scal%npr,scal%npc
       write(6,*)
    endif
    call descinit(desca, nmtx, nmtx, scal%nbl, scal%nbl, 0, 0, &
         scal%icntxt, max(1,scal%npr), info)
    if(info < 0) then
       write(0,'(a,i3,a)') 'Argument number ', -info, ' had an illegal value on entry.'
       call die("descinit error for descaA in inversion")
    else if(info > 0) then
       write(0,*) 'info = ', info
       call die("descinit error for descaA in inversion")
    endif
    ! FHJ: LU factorization of the matrix
    call pzgetrf(nmtx, nmtx, matrix, 1, 1, desca, ipiv, info)
    if (info/=0) then
       if (peinf%inode==0) write(0,*) 'ERROR: got info = ', info, ' in p?getrf'
       call die('p?getrf failed')
    endif
    ! FHJ: tringular inversion of LU decomposition
    allocate(work (10))
    allocate(iwork (10))
    call pzgetri(nmtx, matrix, 1, 1, desca, ipiv, work, -1, iwork, -1, info)
    if (info/=0) then
       if (peinf%inode==0) write(0,*) 'ERROR: got info = ', info, ' in p?getri'
       call die('p?getri failed for query mode')
    endif
    lwork = max(1,int(work(1)))
    liwork = max(1,iwork(1))
    if (allocated(work)) then
       deallocate(work)
    endif
    if(allocated(iwork)) then
       deallocate(iwork)
    endif
    allocate(work (lwork))
    allocate(iwork (liwork))
    call pzgetri(nmtx, matrix, 1, 1, desca, ipiv, work, lwork, iwork, liwork, info)
    if (info/=0) then
       if (peinf%inode==0) write(0,*) 'ERROR: got info = ', info, ' in p?getri'
       call die('p?getri failed')
    endif
    if (allocated(iwork)) then
       deallocate(iwork)
    endif
    if(allocated(work)) then
       deallocate(work)
    endif

    return
  end subroutine zinvert_with_scalapack

  !------------------------------------------------------------
  subroutine zinvert_serial(nmtx, matrix)
    integer, intent(in) :: nmtx
    COMPLEX(DPC), intent(inout) :: matrix(nmtx,nmtx)
    integer :: info, lwork, ipiv(nmtx)
    COMPLEX(DPC), allocatable :: work(:)

    ! FHJ: LU factorization of the matrix
    call zgetrf(nmtx, nmtx, matrix, nmtx, ipiv, info)
    if (info/=0) then
       if (peinf%inode==0) write(0,*) 'ERROR: got info = ', info, ' in ?getrf'
       call die('?getrf failed')
    endif
    ! FHJ: tringular inversion of LU decomposition
    allocate(work (10))
    call zgetri(nmtx, matrix, nmtx, ipiv, work, -1, info)
    if (info/=0) then
       if (peinf%inode==0) write(0,*) 'ERROR: got info = ', info, ' in ?getri'
       call die('?getri failed for query mode')
    endif
    lwork = max(1,int(work(1)))
    if (allocated(work)) then
       deallocate(work)
    endif
    allocate(work (lwork))
    call zgetri(nmtx, matrix, nmtx, ipiv, work, lwork, info)
    if (info/=0) then
       if (peinf%inode==0) write(0,*) 'ERROR: got info = ', info, ' in ?getri'
       call die('?getri failed')
    endif
    if (allocated(work)) then
       deallocate(work)
    endif

    return
  end subroutine zinvert_serial
end module inversion_m
!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

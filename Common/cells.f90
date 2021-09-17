#include "f_defs.h"

!===============================================================================
!
! MODULE:
!
!   cells_m                Originally by FHJ     Last Modified 10/24/2011 (FHJ)
!
!>    Implements cell lists to partition the space and permit quick search of
!!    points.
!
!===============================================================================

module cells_m

  use global_m
  implicit none

  private

  type :: cells_t
     integer :: dims !< Number of dimensions
     integer :: npts !< Number of points/particles/etc
     integer :: n_cells(3) !< Number of cells in each particular dimension
     real(DP) :: dmin(3) !< Min coordinate of cell
     real(DP) :: dmax(3) !< Max coordinate of cell
     real(DP) :: shift(3) !<
     real(DP) :: factor(3) !<
     real(DP), pointer :: pts(:,:) !< The points (dims, npts)
     integer, pointer :: head(:,:,:)
     integer, pointer :: list(:)
     logical :: periodic !< Use periodic boundary conditions?
     logical :: active_dim(3) !< Whether a direction is active or degenerate
  end type cells_t

  public :: cells_t, &
       cells_init, &
       cells_free, &
       cells_get_cell_idx, &
       cells_find_exactly

contains


  !> Puts all points pts(3,:) into a cell structure.
  !! If the cell is periodic, all points are moved to the [0, 1) range.
  !! This subroutine automatically detects the number of dimension spanned by
  !! the points, and it might be a good idea for the caller to check that.
  !! eg: this%dims == xct%idimensions
  ! =======
  ! intkernel.f90:
  ! call cells_init(cells_fi, dqs, periodic=.true.)
  ! dqs(:,ik=1,xct%nkpt_fi) = UNIT_RANGE(kg%f(:,ik=1,xct%nkpt_fi) - kg%f(:,1))
  ! -------
  ! 
  ! -------
  subroutine cells_init(this, pts, periodic)
    type(cells_t), intent(out) :: this
    real(DP), intent(in) :: pts(:,:)
    logical, intent(in) :: periodic

    real(DP) :: dim_vol, prop_const
    integer :: jdim
    logical :: should_write

    PUSH_SUB(cells_init)

    should_write = peinf%verb_debug

    if (peinf%inode==0 .and. should_write) &
         write(6,'(/1x,a)') 'Initializing cell structure.'

    this%periodic = periodic
    this%npts = size(pts, dim=2)

    SAFE_ALLOCATE(this%pts, (3, this%npts))

    ! <->
    if (this%periodic) then
       ! ======
       ! f_defs.h
       ! #define UNIT_RANGE(x) (x - floor(x))
       ! ------
       ! FHJ: this is guaranteed to be in [0, 1)
       this%pts(:,:) = UNIT_RANGE(pts(:,:))
    else
       this%pts(:,:) = pts(:,:)
    endif

    ! ======
    ! FHJ: Avoid problems that a very small dimension may cause
    this%active_dim(:) = .false.
    this%dims = 0
    this%dmin(:) = minval(this%pts(:,:), dim=2)
    this%dmax(:) = maxval(this%pts(:,:), dim=2)

    do jdim=1,3
       if ((this%dmax(jdim)-this%dmin(jdim)) > TOL_Small) then
          this%active_dim(jdim) = .true.
          this%dims = this%dims + 1
          ! <->
          if (this%periodic) then
             this%dmin(jdim) = 0d0
             this%dmax(jdim) = 1d0
          endif
       endif
    enddo

    ! FHJ: This is only to avoid division by zero later on
    this%dmin(:) = this%dmin(:) - TOL_Small
    this%dmax(:) = this%dmax(:) + TOL_Small

    ! FHJ: number of cells is proportional to the length of each active dimension
    !      this should work for any number of dimensions. Proof:
    !        d1*d2*... = V
    !        n1 = c*d1*(N)^(1/D)
    !        n1*n2*... = c^(D)*V*N = N => c = (1/V)^(1/D)
    !      so, const = (N/V)^(1/D)
    this%n_cells(:) = 1

    ! <->
    if (this%dims>0) then
       dim_vol=1.0d0
       do jdim=1,3
          if (this%active_dim(jdim)) dim_vol = dim_vol * (this%dmax(jdim)-this%dmin(jdim))
       enddo

       prop_const = (this%npts/dim_vol)**(1.0d0/this%dims)

       do jdim=1,3
          if (this%active_dim(jdim)) then
             this%n_cells(jdim) = nint(prop_const * (this%dmax(jdim)-this%dmin(jdim)))
             this%n_cells(jdim) = min(max(this%n_cells(jdim), 1), this%npts)
          endif
       enddo
    endif

    this%factor(:) = this%n_cells(:) / (this%dmax(:) - this%dmin(:))
    this%shift(:) = (this%dmax(:) - this%dmin(:)) / (2.0d0*this%n_cells(:))

    if (peinf%inode==0 .and. should_write) then
       write(6,'(1x,a,i0,a)') 'Automatically creating a cell structure for ',this%npts,' points'
       write(6,'(1x,a,i0,a)') 'Found ',this%dims,' dimension(s)'
       write(6,'(1x,a,3(1x,i0))') 'Number of cells:', this%n_cells(1),this%n_cells(2),this%n_cells(3)
       write(6,'(1x,a,i0)') 'Total number of cells: ', this%n_cells(1)*this%n_cells(2)*this%n_cells(3)
    endif

    SAFE_ALLOCATE(this%head, (this%n_cells(1),this%n_cells(2),this%n_cells(3)))
    SAFE_ALLOCATE(this%list, (this%npts))

    ! FHJ: Initialize and populate cells
    this%head(:,:,:) = 0
    this%list(:) = 0

    if (peinf%inode==0 .and. should_write) then
       do jdim=1,3
          write(6,801) jdim,this%dmin(jdim),this%dmax(jdim),this%shift(jdim)*2.0d0
801       format(' Cells [',i1,'], dmin= ',f8.5,' dmax= ',f8.5,' length= ',f12.5)
       enddo
    endif

    call cells_populate(this)

    if (peinf%inode==0 .and. should_write) then
       call cells_show_population(this)
       write(6,'(a)') 'Finished initializing cells.'
    endif

    POP_SUB(cells_init)

  end subroutine cells_init


  !> Used internally by cells_init. Puts all points into the cells, update
  !! head and list structures.
  subroutine cells_populate(this)
    type(cells_t), intent(inout) :: this

    integer :: idx_co
    integer :: cell_idx(3)

    PUSH_SUB(cells_populate)
    do idx_co=1,this%npts
       call cells_get_cell_idx(this, this%pts(:,idx_co), cell_idx)
       this%list(idx_co) = this%head(cell_idx(1),cell_idx(2),cell_idx(3))
       this%head(cell_idx(1),cell_idx(2),cell_idx(3)) = idx_co
    enddo
    POP_SUB(cells_populate)

  end subroutine cells_populate


  !> Returns the cell index for point `pt`.
  subroutine cells_get_cell_idx(this, pt, cell_idx)
    type(cells_t), intent(in) :: this
    real(DP), intent(in) :: pt(3)
    integer, intent(out) :: cell_idx(3)

    ! no push/pop, called too frequently
    cell_idx(:) = idint((pt(:)-this%dmin(:)+this%shift(:))*this%factor(:))
    if (this%periodic) then
       cell_idx(:) = modulo(cell_idx(:), this%n_cells(:)) + 1
    else
       cell_idx(:) = min(max(cell_idx(:)+1, 1), this%n_cells(:))
    endif
  end subroutine cells_get_cell_idx


  !> Find the index idx such that this%pts(:, idx) == pt(:)
  !! The algorithm also looks for nearby cells.
  subroutine cells_find_exactly(this, pt_in, idx_found)
    type(cells_t), intent(in) :: this
    real(DP), intent(in) :: pt_in(3)
    integer, intent(out) :: idx_found

    real(DP) :: pt(3), delta(3)
    integer :: cell_idx(3)
    integer :: i1,i2,i3,j1,j2,j3
    integer :: cell_min(3), cell_max(3)

    ! no push/pop, called too frequently
    pt(:) = pt_in(:)
    if (this%periodic) pt(:) = UNIT_RANGE(pt(:))

    ! See if the point is in this central cell. The subroutine should
    ! usually return here
    call cells_get_cell_idx(this, pt, cell_idx)
    call find_exactly_in_cell(cell_idx, idx_found)
    if (idx_found/=0) return

    ! Otherwise look in neighboring cells. The routine will probably fail...
    call cells_get_search_bounds(this, 1, cell_idx, cell_min, cell_max)
    do i1 = cell_min(1), cell_max(1)
       j1 = cells_fix_index(this, i1, 1)
       cell_idx(1) = j1
       do i2 = cell_min(2), cell_max(2)
          j2 = cells_fix_index(this, i2, 2)
          cell_idx(2) = j2
          do i3 = cell_min(3), cell_max(3)
             j3 = cells_fix_index(this, i3, 3)
             cell_idx(3) = j3
             call find_exactly_in_cell(cell_idx, idx_found)
             if (idx_found/=0) then
                if (peinf%inode==0) then
                   write(0,'(a)') &
                        'WARNING: find_exactly: point was not located in central cell.'
                   write(0,'(a)') &
                        'The cell structure is non-optimal.'
                endif
                return
             endif
          enddo
       enddo
    enddo

  contains

    subroutine find_exactly_in_cell(c_idx, i_found)
      integer, intent(in) :: c_idx(3) !< only search in this cell index
      integer, intent(out) :: i_found !< this%pts(:,i_found) == pts

      integer :: idx_pt

      ! no push/pop, called too frequently
      i_found = 0
      idx_pt = this%head(c_idx(1), c_idx(2), c_idx(3))
      do while (idx_pt/=0)
         delta(:) = this%pts(:, idx_pt) - pt(:)
         if (this%periodic) delta(:) = MIN_RANGE(delta(:))
         if (all(dabs(delta)<TOL_SMALL)) then
            i_found = idx_pt
            return
         endif
         idx_pt = this%list(idx_pt)
      enddo

    end subroutine find_exactly_in_cell

  end subroutine cells_find_exactly

  !> Given a central cell and a Manhattan radius, returns the minimum and
  !! maximum indices to be used in do-loops. The min and max indices obeying
  !! the periodicity or finiteness of the system, and never overlap.
  !! Away from borders, the subroutine simply returns:
  !!  cell_dmin(:) = cell_idx(:) - radius
  !!  cell_dmax(:) = cell_idx(:) + radius
  subroutine cells_get_search_bounds(this, radius, center_idx, cell_min, cell_max)
    type(cells_t), intent(in) :: this
    integer, intent(in) :: radius
    integer, intent(in) :: center_idx(3)
    integer, intent(out) :: cell_min(3), cell_max(3)

    integer :: ii

    ! no push/pop, called too often
    do ii=1,3
       if ((2*radius+1) > this%n_cells(ii)) then
          cell_min(ii) = 1
          cell_max(ii) = this%n_cells(ii)
       else
          cell_min(ii) = center_idx(ii) - radius
          cell_max(ii) = center_idx(ii) + radius
       endif
    enddo
    if (.not.this%periodic) then
       cell_min(:) = max(cell_min(:), 1)
       cell_max(:) = min(cell_max(:), this%n_cells(:))
    endif
  end subroutine cells_get_search_bounds

  !> Fixed an index so that it respects the periodic boundary conditions.
  !! No check is made for non-periodic cells, since functions like
  !! cells_get_search_bounds are already bounded by [1, n_cells]
  integer function cells_fix_index(this, idx_in, dim_)
    type(cells_t), intent(in) :: this
    integer, intent(in) :: idx_in
    integer, intent(in) :: dim_

    ! no push/pop since called too frequently
    cells_fix_index = idx_in
    if (this%periodic) then
       cells_fix_index = modulo(idx_in-1, this%n_cells(dim_)) + 1
    endif
  end function cells_fix_index

  !> Prints a summary of the cells population
  subroutine cells_show_population(this, iunit_)
    type(cells_t), intent(in) :: this
    integer, intent(in), optional :: iunit_

    integer :: j1, j2, j3, idx_co, iunit

    PUSH_SUB(cells_show_population)
    iunit = 6
    if (present(iunit_)) iunit = iunit_
    write(iunit,*) 'Cell population analysis:'
    write(iunit,*)
    write(iunit,900) ' x ',' y ',' z ',' members '
    write(iunit,900) '---','---','---','---------'
900 format(2x, 3(a3,1x), a9)
    do j1 = 1, this%n_cells(1)
       do j2 = 1, this%n_cells(2)
          do j3 = 1, this%n_cells(3)
             write(iunit,'(2x,3(i3,1x))',advance='no') j1,j2,j3
             idx_co = this%head(j1,j2,j3)
             do while (idx_co>0)
                write(iunit,'(1x,i5)',advance='no') idx_co
                idx_co = this%list(idx_co)
             enddo
             write(iunit,*)
          enddo
       enddo
    enddo
    write(iunit,*)
    POP_SUB(cells_show_population)

  end subroutine cells_show_population

  !> Deallocates buffers related to a cell structure
  subroutine cells_free(this)
    type(cells_t), intent(inout) :: this

    PUSH_SUB(cells_free)
    SAFE_DEALLOCATE_P(this%head)
    SAFE_DEALLOCATE_P(this%list)
    SAFE_DEALLOCATE_P(this%pts)
    POP_SUB(cells_free)

  end subroutine cells_free

end module cells_m

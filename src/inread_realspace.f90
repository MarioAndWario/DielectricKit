#include "f_defs.h"

!==============================================================================
!
! Routine:
!
! (1) inread
!
!     Read information from input file realspace.inp.
!
!==============================================================================

module inread_m
  use global_m
  use realspace_common_m
  implicit none
  public :: inread
contains

  subroutine inread(peps)
    type (realspace_t), intent(out) :: peps
    character*256 :: blockword,keyword,line,errmsg
    integer :: ii,iostat,jj
    integer :: ir2
    real(DP) :: r2_temp(3,1000), freq_real, freq_imag

    r2_temp(:,:)=0.0D0
    PUSH_SUB(inread)
#ifdef MPI
    ! Non-root nodes should wait for root to read the whole file.
    ! That way, we can be sure root gets a chance to write errors before
    ! any die call is issued by another node. Root calls MPI_Barrier below.
    if(peinf%inode /= 0) call MPI_Barrier(MPI_COMM_WORLD, mpierr)
#endif
    call open_file(8,file='realspace.inp',form='formatted',status='old')

    ! Set default values
    peps%nfq = 0
    peps%filename = "epsmat.h5"
    peps%qgrid(:) = 0
    peps%nsuper(:) = 0
    peps%nr2 = 0
    peps%ispin=1
    peps%unfold = .false.
    peps%downsample(:)=1
    peinf%verbosity = 1
    peps%ecut = 0.0D0
    peps%freq_target = (0.0D0, 0.0D0)
    peps%FFTfactor = 1
    peps%low_comm = .true.
    peps%mid_comm = .false.

    !> Output flags
    peps%realpart = .false.
    peps%imagpart = .false.
    peps%abs = .false.
    peps%abs2 = .false.
    peps%high_resolution = .false.
    peps%epsinvhead = 0.0D0

    ! Never ending loop...
    do while(0 .eq. 0)
       ! Actually the loop ends when the end of the file is reached
       read(8,'(a256)',iostat=iostat) line
       if(iostat < 0) exit
       ! Skip comment lines
       if(len_trim(line).eq.0) cycle
       if(line(1:1).eq.'#') cycle
       ! Determine keyword:
       keyword=line(1:scan(line," ")-1)
       line=adjustl(line(scan(line," ")+1:256))

       if(trim(keyword).eq.'begin') then
          blockword=line(1:scan(line," ")-1)
          ii=0
          do while(trim(line).ne.'end')
             read(8,'(a256)',iostat=iostat) line
             if(iostat /= 0) then
                write(errmsg,'(3a)') 'The end of the file was reached while reading elements of the ', &
                     trim(blockword),' block.'
                call die(errmsg, only_root_writes = .true.)
             endif
             if(trim(line).ne.'end') then
                ii=ii+1
                if(trim(blockword).eq.'r2') then
                   read(line,*) (r2_temp(jj,ii),jj=1,3)
                endif
             end if
          enddo
          if(trim(blockword).eq.'r2') then
             peps%nr2=ii
          endif
       elseif(trim(keyword).eq.'verbosity') then
          read(line,*,err=110) peinf%verbosity
       elseif(trim(keyword).eq.'use_symmetries') then
          peps%unfold = .true.
       elseif(trim(keyword).eq.'filename') then
          read(line,*,err=110) peps%filename
          peps%filename = TRUNC(peps%filename)
       elseif(trim(keyword).eq.'cutoff') then
          read(line,*,err=110) peps%ecut
          ! elseif(trim(keyword).eq.'qgrid') then
          !    read(line,*,err=110) peps%qgrid(1:3)
          !    peps%nfq = PRODUCT(peps%qgrid(:))
       elseif(trim(keyword).eq.'frequency') then
          read(line,*,err=110) freq_real, freq_imag
          peps%freq_target = DCMPLX(freq_real, freq_imag)
          ! elseif(trim(keyword).eq.'supercell_size') then
          !    read(line,*,err=110) (peps%nsuper(ii),ii=1,3)
       elseif(trim(keyword).eq.'real_part') then
          peps%realpart = .TRUE.
       elseif(trim(keyword).eq.'imag_part') then
          peps%imagpart = .TRUE.
       elseif(trim(keyword).eq.'abs') then
          peps%abs = .TRUE.
       elseif(trim(keyword).eq.'abs2') then
          peps%abs2 = .TRUE.
       elseif(trim(keyword).eq.'high_resolution') then
          peps%high_resolution = .TRUE.
       elseif(trim(keyword).eq.'low_comm') then
          peps%low_comm = .TRUE.
          peps%mid_comm = .FALSE.
       elseif(trim(keyword).eq.'mid_comm') then
          peps%low_comm = .FALSE.
          peps%mid_comm = .TRUE.
       elseif(trim(keyword).eq.'downsample') then
          read(line,*,err=110) peps%downsample
       elseif(trim(keyword).eq.'plot_spin') then
          read(line,*,err=110) peps%ispin
       elseif(trim(keyword).eq.'FFTfactor') then
          read(line,*,err=110) peps%FFTfactor
       elseif(trim(keyword).eq.'epsinvhead') then
          read(line,*,err=110) peps%epsinvhead
       else
          write(errmsg,'(3a)') 'Unexpected keyword ', trim(keyword), ' was found in realspace.inp.'
          call die(errmsg, only_root_writes = .true.)
       end if
    enddo

    call close_file(8)

    call peinfo_set_verbosity()

    if (peinf%inode .eq. 0) then
       if (peps%ecut < TOL_ZERO) then
          call die("Must set peps%ecut", only_root_writes = .true.)
       endif
       ! if (peps%nfq .le. 0) then
       !    call die("Number of qpoints <= 0, please check qgrid", only_root_writes = .true.)
       ! endif
       if (peps%nr2 .le. 0) then
          call die("Number of r2 points <= 0, please check r2 list", only_root_writes = .true.)
       endif
       if ( .not. (peps%realpart .or. peps%imagpart .or. peps%abs .or. peps%abs2) ) then
          call die("Must select one output option: real_part, imag_part, abs, abs2", only_root_writes = .true.)
       endif

       if (peps%low_comm) then
          write(6,'(1X,A)') "low_comm: We will store entire epsmat.h5 on each proc."
       endif

       if (peps%mid_comm) then
          write(6,'(1X,A)') "mid_comm: We will store part of epsmat.h5 on each proc."
       endif

       write(6,'(1X,A)') "Dealing with file:"//TRUNC(peps%filename)
       write(6,'(1X)')

       write(6,'(1X,A)') "List of output files:"
       if (peps%realpart) then
          write(6,'(1X,A)') 'Output real part of epsinv.'
       endif
       if (peps%imagpart) then
          write(6,'(1X,A)') 'Output imaginary part of epsinv.'
       endif
       if (peps%abs) then
          write(6,'(1X,A)') 'Output absolute value of epsinv.'
       endif
       if (peps%abs2) then
          write(6,'(1X,A)') 'Output module square of epsinv.'
       endif
       if (peps%high_resolution) then
          write(6,'(1X,A)') 'We will use dense FFT grid (gvec%FFTgrid) to plot high-resolution image.'
       else
          write(6,'(1X,A)') 'We will use self-adaptive FFT grid (using peps cutoff) to plot low-resolution image.'
       endif

       write(6,'(A)')
       if (peps%unfold) then
          write(6,'(1X,A)') 'Unfolding epsinv matrix using symmetry.'
       else
          write(6,'(1X,A)') 'No symmetry is used.'
       endif
       write(6,'(A)')
       ! write(6,'(1X,A,3(i4,1x))') 'Supercell size:', peps%nsuper
       ! write(6,'(1X,A,3(i4,1x),A,I4,A)') 'qgrid :', peps%qgrid(:), " with ", peps%nfq, " qpoints in full BZ"
       write(6,'(1X,A,3(i4,1x))') 'Real-space grid will be downsampled by a factor of:', peps%downsample
       write(6,'(1X,A,3(i4,1x))') 'Real-space FFTgrid will be multiplied by a factor of:', peps%FFTfactor

       write(6,'(1X,A,F10.3,A)') 'G-cutoff for epsinv used in FFT : ', peps%ecut, " Ry"
       write(6,'(A)')
    endif

    if (TRUNC(peps%filename) .eq. "epsmat.h5") then
       if ( ABS(peps%epsinvhead) < 1.0D-10 ) then
          call die("Must set epsinvhead for epsmat.h5.", only_root_writes=.true.)
       else
          if (peinf%inode .eq. 0) then
             write(6,'(1X,A,F17.10)') 'epsinv(G=Gp=0)(q-->0) :', peps%epsinvhead
          endif
       endif
    elseif (TRUNC(peps%filename) .eq. "chimat.h5") then
       peps%epsinvhead = 0.0D0
       if (peinf%inode .eq. 0) then
          write(6,'(1X,A)') 'chi(G=Gp=0)(q-->0) = 0'
       endif
    else
       call die("Invalid filename.", only_root_writes=.true.)
    endif

    if (ANY(peps%downsample<1)) then
       call die('Invalid value for real-space grid subsampling.', only_root_writes = .true.)
    endif
    if(peps%ispin < 0 .or. peps%ispin > 2) call die("plot_spin out of bounds", only_root_writes = .true.)
    SAFE_ALLOCATE(peps%r2, (3, peps%nr2))
    peps%r2(1:3, 1:peps%nr2) = r2_temp(1:3, 1:peps%nr2)

    if (peinf%inode .eq. 0) then
       ! !> Check that all the r2 are within the peps%nsuper(3) range
       ! if ( (any(peps%r2(1,:) > peps%nsuper(1))) .or. (any(peps%r2(2,:) > peps%nsuper(2))) .or. (any(peps%r2(3,:) > peps%nsuper(3))) ) then
       !    call die("r2 beyond range", only_root_writes=.TRUE.)
       ! elseif ( (any(peps%r2(1,:) < 0)) .or. (any(peps%r2(2,:) < 0)) .or. (any(peps%r2(3,:) < 0)) ) then
       !    call die("r2 beyond range", only_root_writes=.TRUE.)
       ! endif
       ! if ( ANY( (peps%nsuper(:) - peps%qgrid(:)) > 0 ) ) then
       !    call die("nsuper cannot exceed qgrid", only_root_writes = .TRUE.)
       ! endif
       ! if ( ANY( (peps%nsuper(:) .le. 0 ) ) ) then
       !    call die("Please check supercell_size", only_root_writes = .TRUE.)
       ! endif
       ! if ( ANY( (peps%qgrid(:) .le. 0 ) ) ) then
       !    call die("Please check qgrid", only_root_writes = .TRUE.)
       ! endif
       write(6,'(1X,A,I5,A)') "We will consider ", peps%nr2, " r2 points (in fractional coordinates): "
       do ir2 = 1, peps%nr2
          write(6,'(1X,"#",I5,3F15.6)') ir2, peps%r2(:,ir2)
       enddo
    endif

#ifdef MPI
    !> root lets the others go after it is done reading (see beginning of function)
    if (peinf%inode .eq. 0) call MPI_Barrier(MPI_COMM_WORLD, mpierr)
#endif
    POP_SUB(inread)
    return
110 write(errmsg,'(3a)') 'Unexpected characters were found while reading the value for the keyword ', trim(keyword), '. '
    call die(errmsg, only_root_writes = .true.)
  end subroutine inread
end module inread_m

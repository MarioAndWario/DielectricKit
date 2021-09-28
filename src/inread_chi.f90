#include "f_defs.h"

module inread_chi_m
  use global_m
  implicit none
  private
  public :: inread_chi
contains
  !> supply optionals for absorption, and do not for inteqp
  subroutine inread_chi(pol)
    type (polarizability), intent(out) :: pol
    character*256 :: blockword, keyword, line, errmsg
    integer :: ii,iostat, i, qflag, qflag_, rq_counter, iq
    logical :: unknown_keyword, found
    real(DP) :: rq_input(3, MAX_KPTS)
    integer, parameter :: unit_input=234
    real(DP) :: tmpFreq, omegapp_max_, tany_max_, y_
    integer :: ifreq_real, ifreq_imag, ifreq, ifreqCounter, nfreq_imag_, ifreq_imag_

    PUSH_SUB(inread_chi)
    qflag = 0
    qflag_ = 0

    call open_file(unit_input, file='chi.inp', form='formatted', status='old')

    !> Every proc will read in absorption.inp file
#ifdef MPI
    ! Non-root nodes should wait for root to read the whole file.
    ! That way, we can be sure root gets a chance to write errors before
    ! any die call is issued by another node. Root calls MPI_Barrier below.
    if (peinf%inode .ne. 0) call MPI_Barrier(MPI_COMM_WORLD, mpierr)
#endif

    !> Set default values of pol
    pol%nq0 = 0
    pol%nq1 = 0
    pol%nq = 0
    pol%nvb = 0
    pol%ncb = 0
    pol%skip_nvb=0
    pol%skip_ncb=0
    pol%ecuts = 0.0d0
    pol%non_uniform = .false.
    pol%subsample = .false.

    !> ------------------
    !> Frequency related variables
    pol%freq_dep=0
    pol%freq_dep_method=2
    pol%nFreq=1
    pol%nfreq_real=1
    pol%nfreq_imag=15
    pol%dInitFreq=0.0d0
    pol%dDeltaFreq=-1d0
    pol%delta_freq_imag=2.0D0
    pol%dFreqCutoff=-1d0
    pol%dBrdning= -1000D0
    pol%nBrdning = 1
    pol%Brdning_stepsize = 0.0D0

    pol%nSFreq=1
    pol%dInitSFreq=0.0d0
    pol%dDeltaSFreq=-1d0
    pol%dSFreqStepIncrease=0d0
    pol%dSFreqCutoff1=-1d0
    pol%dSFreqCutoff2=-1d0
    !> ------------------
    pol%fullConvLog=0
    pol%icutv=0
    pol%truncval(1)=0.0d0
    pol%truncval(2)=0.0d0
    pol%truncval(3)=0.0d0
    ! pol%valueq0=0
    ! pol%iqexactlyzero=0
    pol%ncrit=0
    pol%use_hdf5 = .true.
    pol%efermi_input=0.0d0
    pol%rfermi=.true.
    pol%gcomm=-1
    pol%os_opt_ffts=0
    pol%restart=.false.
    pol%min_fftgrid=.true.
    pol%lin_denominator=0d0
    ! pol%nfreq_group=1
    pol%os_hdf5=.false.
    pol%skip_epsilon=.false.
    pol%skip_chi=.false.
    pol%degeneracy_check_override=.false.
    peinf%npools=0
    pol%eqp_corrections=.false.
    pol%intraband_flag=0
    pol%intraband_overlap_min=0.9d0
    pol%patched_sampling=.false.
    pol%qgrid(:) = 0
    ! variables for subspace truncation method in epsilon
    pol%subspace = .FALSE.
    pol%chi_eigenvalue_cutoff = 1.0d-6
    pol%neig_sub_input = -1
    pol%use_elpa = .FALSE.
    pol%keep_full_eps_static = .TRUE.
    pol%matrix_in_subspace_basis = .FALSE.
    ! variables for nonblocking scheme
    pol%nonblocking_cyclic = .FALSE.
    pol%mtxel_2=.false.
    pol%matrix_type = 2
    pol%timeordered=.true.
    pol%resetrealfreq=.false.
    pol%k_plus_q=.false.
    pol%time_reversal=.true.
    pol%full_axis_frequency=.false.

    !> Never ending loop...
    do while(.true.)
       !> Actually the loop ends when the end of the file is reached
       read(unit_input,'(a256)',iostat=iostat) line
       if(iostat < 0) exit

       !> Skip comment lines
       if (len_trim(line).eq.0 .or. line(1:1).eq.'#') cycle
       !> Determine keyword:
       keyword=line(1:scan(line," ")-1)
       line=adjustl(line(scan(line," ")+1:256))
       unknown_keyword = .false.

       if (trim(keyword).eq.'begin') then
          blockword = line(1:scan(line," ")-1)
          rq_counter = 0
          do while(trim(line).ne.'end')
             read(unit_input,'(a256)',end=105) line
             if (trim(line).ne.'end') then
                rq_counter = rq_counter + 1
                if (rq_counter > MAX_KPTS) then
                   call die("Increase MAX_KPTS.", only_root_writes=.true.)
                endif
                if (trim(blockword).eq.'qpoints') then
                   read(line,*,iostat=iostat) rq_input(:,rq_counter), qflag
                   if (iostat /= 0) then
                      write(errmsg,'(3a)') 'Unexpected characters were found while reading elements of the ', trim(blockword),' block.'
                      call die(errmsg, only_root_writes = .true.)
                   endif
                   if (qflag .ne. 0 .and. qflag .ne. 1) then
                      call die("Invalid qflag", only_root_writes=.true.)
                   endif
                   !> qflag == 0 ==> q0
                   !> qflag == 1 ==> q1
                   !> all q0 must be placed before all q1
                   if (qflag < qflag_) then
                      call die("All q0 must be placed before all q1.", only_root_writes=.true.)
                   else
                      qflag_ = qflag
                   endif
                   if (qflag .eq. 0) then
                      pol%nq0 = pol%nq0 + 1
                      if (NORM2(rq_input(:, rq_counter)) < TOL_Zero) then
                         call die("q0 cannot be exactly zero.", only_root_writes=.true.)
                      endif
                   endif
                else
                   write(errmsg,'(3a)') 'Unexpected blockword ', trim(blockword), ' was found in epsilon.inp.'
                   call die(errmsg, only_root_writes = .true.)
                endif
             endif
          enddo
          if (trim(blockword).eq.'qpoints') then
             pol%nq = rq_counter
          endif
       elseif(trim(keyword).eq.'number_val_bands') then
          read(line,*,err=110) pol%nvb
       elseif(trim(keyword).eq.'number_cond_bands') then
          read(line,*,err=110) pol%ncb
       elseif(trim(keyword).eq.'chi_cutoff') then
          read(line,*,err=110) pol%ecuts
       elseif(trim(keyword).eq.'restart') then
          pol%restart = .true.
       elseif(trim(keyword).eq.'skip_nvb') then
          read(line,*,err=110) pol%skip_nvb
       elseif(trim(keyword).eq.'skip_ncb') then
          read(line,*,err=110) pol%skip_ncb
       elseif(trim(keyword).eq.'qgrid') then
          read(line,*,err=110) pol%qgrid(1:3)
       elseif(trim(keyword).eq.'frequency_dependence') then
          read(line,*,err=110) pol%freq_dep
       elseif(trim(keyword).eq.'frequency_dependence_method') then
          read(line,*,err=110) pol%freq_dep_method
       elseif(trim(keyword).eq.'init_frequency') then
          read(line,*,err=110) pol%dInitFreq
       elseif(trim(keyword).eq.'delta_frequency') then
          read(line,*,err=110) pol%dDeltaFreq
       elseif(trim(keyword).eq.'frequency_cutoff') then
          read(line,*,err=110) pol%dFreqCutoff
       elseif(trim(keyword).eq.'number_imaginary_freqs') then
          read(line,*,err=110) pol%nfreq_imag
       elseif(trim(keyword).eq.'broadening') then
          read(line,*,err=110) pol%dBrdning
       elseif(trim(keyword).eq.'nbroadening') then
          read(line,*,err=110) pol%nBrdning
       elseif(trim(keyword).eq.'broadening_stepsize') then
          read(line,*,err=110) pol%Brdning_stepsize
       elseif(trim(keyword).eq.'delta_freq_imag') then
          read(line,*,err=110) pol%delta_freq_imag
       elseif(trim(keyword).eq.'eqp_corrections') then
          pol%eqp_corrections=.true.
       elseif(trim(keyword).eq.'time-ordered') then
          pol%timeordered=.true.
       elseif(trim(keyword).eq.'retarded') then
          pol%timeordered=.false.
       elseif(trim(keyword).eq.'reset_real_freq') then
          pol%resetrealfreq=.true.
       elseif(trim(keyword).eq.'k_plus_q') then
          pol%k_plus_q=.true.
       elseif(trim(keyword).eq.'k_minus_q') then
          pol%k_plus_q=.false.
       elseif(trim(keyword).eq.'broken_time_reversal') then
          pol%time_reversal=.false.
       elseif(trim(keyword).eq.'full_axis_frequency') then
          pol%full_axis_frequency=.true.
       elseif(trim(keyword).eq.'spherical_truncation') then
          pol%icutv=2
       elseif(trim(keyword).eq.'cell_wire_truncation') then
          pol%icutv=4
       elseif(trim(keyword).eq.'cell_box_truncation') then
          pol%icutv=5
       elseif(trim(keyword).eq.'cell_slab_truncation') then
          pol%icutv=6
       else
          if (.not. found) then
             write(errmsg,'(3a)') 'Unexpected keyword ', trim(keyword), ' was found.'
             call die(errmsg, only_root_writes = .true.)
          endif
       endif
    enddo

    pol%nq1 = pol%nq - pol%nq0
    SAFE_ALLOCATE(pol%qpt, (3, pol%nq))
    pol%qpt(1:3, 1:pol%nq) = rq_input(1:3, 1:pol%nq)

    if (peinf%inode .eq. 0) then
       write(*,'(1X,A)') "q0 vectors:"
       do iq = 1, pol%nq0
          write(*,'(I5,3F15.8)') iq, pol%qpt(:,iq)
       enddo
       write(*,'(1X,A)') "q1 vectors:"
       do iq = pol%nq0+1, pol%nq
          write(*,'(I5,3F15.8)') iq, pol%qpt(:,iq)
       enddo
    endif

    call peinfo_set_verbosity()

    if (pol%ecuts .le. 0) then
       call die("Must set pol%ecuts", only_root_writes=.true.)
    endif

    if (peinf%inode .eq. 0) then
       if ((pol%nvb <= 0) .or. (pol%ncb <=0 )) then
          call die("Invalid number of bands.", only_root_writes=.true.)
       endif
       write(6,'(1X,A,F10.3,A,I5,A,I5)') "We will calculate polarizability with cutoff = ", pol%ecuts, " nvb = ", pol%nvb, " ncb = ", pol%ncb
       if ( pol%skip_nvb .lt. 0 ) then
          call die('Invalid skip_nvb .', only_root_writes=.true.)
       elseif ( pol%skip_ncb .lt. 0 ) then
          call die('Invalid skip_ncb .', only_root_writes=.true.)
       else
          write(6,'(1X,A,I5,A,I5,A)') "We will skip transitions from the top ", pol%skip_nvb, " valence states to the bottom ", pol%skip_ncb, " conduction states"
       endif
       write(6,'(1X,A)')
       if (pol%k_plus_q) then
          write(6,'(1X,A)') "Use |k+q> for valence states."
       else
          write(6,'(1X,A)') "Use |k-q> for valence states."
       endif
       write(6,'(1X,A)')
    endif

    if (pol%freq_dep .ne. 0) then
       if (pol%freq_dep .ne. 2) then
          call die("Only support freq_dep = 2 for FF", only_root_writes = .true.)
       endif
    endif

    !<- FF ->
    if (pol%freq_dep .eq. 2) then
       call initialize_FF()
       !> <- GPP ->
    elseif (pol%freq_dep .eq. 0) then
       if (peinf%inode .eq. 0) then
          write(*,'(1X,A)') "Use GPP."
          write(*,'(1X,A)')
       endif
       if (pol%nBrdning .ne. 1) then
          call die("GPP must use pol%nBrdning = 1.", only_root_writes=.true.)
       endif
       pol%nfreq_imag = 0
       pol%nfreq_real = 1
       pol%nfreq = 1
       SAFE_ALLOCATE(pol%dFreqGrid,(pol%nfreq))
       SAFE_ALLOCATE(pol%dFreqBrd,(pol%nfreq))
       pol%dFreqGrid = 0D0
       pol%dFreqBrd = ZERO
    else
       call die("inread: only FF+CD or GPP are supported.", only_root_writes=.true.)
    endif

#ifdef MPI
    if (peinf%inode .eq. 0) call MPI_Barrier(MPI_COMM_WORLD, mpierr)
#endif

    call close_file(unit_input)

    POP_SUB(inread_chi)
    return

105 write(errmsg,'(3a)') 'The end of the file was reached while reading elements of the ', trim(blockword),' block.'
    call die(errmsg, only_root_writes = .true.)

110 write(errmsg,'(3a)') 'Unexpected characters were found while reading the value for the keyword ', trim(keyword), '. '
    call die(errmsg, only_root_writes = .true.)

  contains

    subroutine initialize_FF()
      if (peinf%inode .eq. 0) then
         write(*,'(1X,A)') "Use full-frequency."
         if (pol%time_reversal) then
            write(*,'(1X,A)') "With time-reversal symmetry"
            if (pol%full_axis_frequency) then
               write(*,'(1X,A)') "Use full-axis frequencies."
            else
               write(*,'(1X,A)') "Use positive-axis frequencies."
            endif
         else
            pol%full_axis_frequency = .true.
            write(*,'(1X,A)') "Without time-reversal symmetry, must use full-axis frequencies (set pol%full_axis_frequency = T)"
         endif
         write(*,'(1X,A)')
      endif

      if (pol%freq_dep_method .ne. 2) then
         call die("inread: only FF+CD is supported.", only_root_writes=.true.)
      endif
      if (pol%dDeltaFreq .le. TOL_Zero) then
         call die("Check pol%dDeltaFreq for FF calculations.", only_root_writes=.true.)
      endif
      if (pol%delta_freq_imag .lt. TOL_Zero) then
         call die("Check pol%delta_freq_imag for FF calculations.", only_root_writes=.true.)
      endif

      !> We now allow negative broadening pol%dBrdning and negative pol%dFreqCutoff1 as input
      !> If pol%dFreqCutoff1 < 0, there will not be any \omega + i \eta near the real-axis
      !> We will only reset pol%dBrdning when they are default value -1000 (default)
      if (ABS(pol%dBrdning+1000d0) < TOL_Zero) pol%dBrdning = 0.1d0
      if (ABS(pol%dFreqCutoff+1.0D0) < TOL_Zero) pol%dFreqCutoff = 200d0
      if (pol%dDeltaFreq < 0d0) pol%dDeltaFreq = pol%dBrdning
      if (pol%nBrdning > 1 .and. pol%dBrdning .le. 0.0D0) then
         call die("When nBrdning > 1, pol%dBrdning cannot <= 0.", only_root_writes=.true.)
      endif
      if (peinf%inode .eq. 0) then
         write(*,'(1X,A,I5,A,F10.3,A,F10.3,A)') "We will consider ", pol%nBrdning, " broadening parameters starting from ", pol%dBrdning, "eV with a stepsize of ", pol%Brdning_stepsize, "eV"
      endif

      if (peinf%inode .eq. 0) then
         if (.not. pol%timeordered) then
            write(6,705)
         else
            write(6,706)
         endif
         write(6,'(1x,a)') 'Parameters for the full-frequency calculation:'
         write(6,'(1x,a,i0)') '- frequency_dependence_method: ', pol%freq_dep_method
         write(6,'(1x,a,f0.6)') '- init_frequency: ', pol%dInitFreq
         write(6,'(1x,a,f0.6)') '- broadening: ', pol%dBrdning
         write(6,'(1x,a,f0.6)') '- delta_frequency: ', pol%dDeltaFreq
         write(6,'(1x,a,f0.6)') '- frequency_cutoff: ', pol%dFreqCutoff
         write(6,'(1x,a,i0)') '- number_real_freqs: ', pol%nfreq_real
         write(6,'(1x,a,i0)') '- number_imaginary_freqs: ', pol%nfreq_imag
         write(6,'(1x,a,i0)') '- number_freqs: ', pol%nfreq
         write(6,'(1x,a,f0.6)') '- delta_freq_imag: ', pol%delta_freq_imag
         write(6,'(1x,a,i0)') '- number_broadenings: ', pol%nBrdning
         write(6,'(1x,a,f0.6)') '- broadening_stepsize: ', pol%Brdning_stepsize
      endif
705   format(1x,'Computing the full frequency-dependent retarded polarizability matrix (Contour-Deformation formulation)',/)
706   format(1x,'Computing the full frequency-dependent time-ordered polarizability matrix (Contour-Deformation formulation)',/)

      tmpFreq = pol%dInitFreq
      ifreqCounter = 0
      do while (tmpFreq .le. pol%dFreqCutoff+TOL_SMALL)
         ifreqCounter = ifreqCounter+1
         !< since pol%dFreqCutoff2 = pol%dFreqCutoff1 (default), we will skip the else section
         !< (default) pol%dDeltaFreq = pol%dBrdning
         tmpFreq = tmpFreq+pol%dDeltaFreq
      enddo
      pol%nfreq_real = iFreqCounter
      pol%nfreq = pol%nfreq_real + pol%nfreq_imag

      ! !> This condition plays nicely with the condition above when nfreq_group .gt. npes, i.e. the conditions ensure
      ! !> the number of processors will always be re-set to a legitimate/sensible number
      ! if (pol%nfreq_group .gt. pol%nFreq) then
      !    call die("pol%nfreq_group > nFreq", only_root_writes=.true.)
      !    write(0,'(/1x,a)') 'WARNING: Number of frequency groups cannot exceed number of frequencies computed'
      !    write(0,'(/,1x,2(a,i5),/)') 'Resetting nfreq_group',pol%nfreq_group,'to number of frequencies computed', pol%nfreq
      !    pol%nfreq_group = peinf%npes
      ! endif
      ! if (pol%subspace) then
      !    call die("Low-rank not supported", only_root_writes=.true.)
      !    if (pol%nfreq_group .lt. peinf%npes) then
      !       ! this is a workaround for the subspace method in order to use all proc at freq=0
      !       pol%nfreq_group = 1
      !    endif
      ! endif

      !> Consider broken time-reversal symmetry
      if (pol%full_axis_frequency) then
         !> [0, 1] --> [-1, 0, 1]
         pol%nfreq_real = pol%nfreq_real * 2 - 1
         pol%nfreq_imag = pol%nfreq_imag * 2 - 1
         pol%nfreq = pol%nfreq_real + pol%nfreq_imag
      endif

      !> real(DP)
      SAFE_ALLOCATE(pol%dFreqGrid,(pol%nFreq))
      !> complex(DPC)
      SAFE_ALLOCATE(pol%dFreqBrd,(pol%nFreq))

      !> Positive-axis frequency
      if (.not. pol%full_axis_frequency) then
         if (ABS(pol%dInitFreq) > TOL_ZERO) then
            if (peinf%inode .eq. 0) then
               write(6,'(1X,A,F12.5)') "[WARNING] You have non-zero pol%dInitFreq = ", pol%dInitFreq
            endif
         endif

         do ifreq_real = 1, pol%nfreq_real
            pol%dFreqGrid(ifreq_real) = DBLE(ifreq_real-1) * pol%dDeltaFreq + pol%dInitFreq
            pol%dFreqBrd(ifreq_real) = pol%dBrdning * (0.0D0,1.0D0)
         enddo

         if (pol%nfreq_imag .gt. 0) then
            if (pol%nfreq_imag .eq. 1) then
               pol%dFreqGrid(pol%nfreq_real + 1) = 0.0D0
               pol%dFreqBrd(pol%nfreq_real + 1) = (0.0D0, 0.0D0)
            else
               omegapp_max_ = pol%delta_freq_imag * ( 2 * pol%nfreq_imag / PI_D )**2
               tany_max_ = tan(DBLE(pol%nfreq_imag-1)/DBLE(pol%nfreq_imag)*PI_D/2.0D0)
               do ifreq_imag = 1, pol%nfreq_imag
                  y_ = DBLE(ifreq_imag-1) / DBLE(pol%nfreq_imag) * PI_D / 2.0D0
                  pol%dFreqGrid(pol%nfreq_real + ifreq_imag) = 0.0D0
                  pol%dFreqBrd(pol%nfreq_real + ifreq_imag) = omegapp_max_ / tany_max_ * tan(y_) * (0.0D0,1.0D0)
               enddo
            endif
         endif
         !> Full-axis frequency
      else
         if (ABS(pol%dInitFreq) > TOL_ZERO) then
            call die("Full-axis only supports frequencies starting from zero.", only_root_writes=.true.)
         endif

         do ifreq_real = 1, pol%nfreq_real
            pol%dFreqGrid(ifreq_real) = DBLE(ifreq_real-1) * pol%dDeltaFreq - DBLE((pol%nfreq_real-1)/2) * pol%dDeltaFreq
            pol%dFreqBrd(ifreq_real) = pol%dBrdning * (0.0D0,1.0D0)
         enddo

         if (pol%nfreq_imag .gt. 0) then
            if (pol%nfreq_imag .eq. 1) then
               pol%dFreqGrid(pol%nfreq_real + 1) = 0.0D0
               pol%dFreqBrd(pol%nfreq_real + 1) = (0.0D0, 0.0D0)
            else
               nfreq_imag_ = (pol%nfreq_imag + 1) / 2
               omegapp_max_ = pol%delta_freq_imag * ( 2 * nfreq_imag_ / PI_D )**2
               tany_max_ = tan(DBLE(nfreq_imag_-1)/DBLE(nfreq_imag_)*PI_D/2.0D0)
               do ifreq_imag = 1, pol%nfreq_imag
                  ifreq_imag_ = ifreq_imag - nfreq_imag_ + 1
                  y_ = DBLE(ifreq_imag_-1) / DBLE(nfreq_imag_) * PI_D / 2.0D0
                  pol%dFreqGrid(pol%nfreq_real + ifreq_imag) = 0.0D0
                  pol%dFreqBrd(pol%nfreq_real + ifreq_imag) = omegapp_max_ / tany_max_ * tan(y_) * (0.0D0,1.0D0)
               enddo
            endif
         endif
      endif

      if (pol%timeordered) then
         !> If pol%resetrealfreq = T, set first real-axis frequency to be exactly 0
         if (pol%resetrealfreq) then
            if (peinf%inode .eq. 0) then
               write(*,'(1X,A)') "Reset real-axis frequency at omega=0 to be exacly 0."
            endif
            do ifreq_real = 1, pol%nfreq_real
               if (ABS(pol%dFreqGrid(ifreq_real)) < TOL_ZERO) then
                  if (peinf%inode .eq. 0) then
                     write(*,*) "ifreq_zero = ", ifreq_real
                  endif
                  pol%dFreqBrd(ifreq_real) = (0.0D0, 0.0D0)
               endif
            enddo
         endif
      endif

      !> output frequencies
      if (peinf%inode .eq. 0) then
         if (pol%timeordered) then
            write(*,'(1X,A,F0.5,A)') "- omegapp_max = ", omegapp_max_, "eV"
            write(*,'(1X,A,F0.5,A)') "- delta_freq_imag = ", pol%delta_freq_imag, " eV"
            write(*,'(1X,A,F0.5,A)') "- pol%dDeltaFreq = ", pol%dDeltaFreq, " eV"
            write(*,'(1X,A,I5)') "- pol%nfreq_imag = ", pol%nfreq_imag
            write(*,'(1X,A)')

            write(*,'(1X,A)') "Time-ordered <=> varepsilon is non-zero for real-axis frequencies"
            if (pol%nfreq_real .gt. 0) then
               write(*,'(1X,A,I10)') "Real-axis frequencies: nfreq_real = ", pol%nfreq_real
               do ifreq = 1, pol%nfreq_real
                  !> pol%dFreqBrd(ifreq) is frequency-dependent varepsilon in this case
                  write(*,'(1X,A,I5,A,F12.5,A,F12.5,A)') "#", ifreq, " omega = (", pol%dFreqGrid(ifreq), " ) eV varepsilon = ", DIMAG(pol%dFreqBrd(ifreq)), " eV"
               enddo
            endif

            if (pol%nfreq_imag .gt. 0) then
               write(*,'(1X,A,I10)') "Imaginary-axis frequencies: nfreq_imag = ", pol%nfreq_imag
               do ifreq = pol%nfreq_real+1, pol%nfreq
                  write(*,'(1X,A,I5,A,"(",F12.5, " + i ", F12.5,")",A)') "#", ifreq, " omega = ", pol%dFreqGrid(ifreq), DIMAG(pol%dFreqBrd(ifreq)), " eV varepsilon = 0 eV"
               enddo
            endif
         else
            write(*,'(1X,A)') "Retarded <==> varepsilon = 0 eV"
            if (pol%nfreq_real .gt. 0) then
               write(*,'(1X,A,I10)') "Real-axis frequencies: nfreq_real = ", pol%nfreq_real
               do ifreq = 1, pol%nfreq_real
                  write(*,'(1X,A,I5,A,"(",F12.5, " + i ", F12.5,")",A)') "#", ifreq, " omega = ", pol%dFreqGrid(ifreq), DIMAG(pol%dFreqBrd(ifreq)), " eV"
               enddo
            endif

            if (pol%nfreq_imag .gt. 0) then
               write(*,'(1X,A,I10)') "Imaginary-axis frequencies: nfreq_imag = ", pol%nfreq_imag
               do ifreq = pol%nfreq_real+1, pol%nfreq
                  write(*,'(1X,A,I5,A,"(",F12.5, " + i ", F12.5,")",A)') "#", ifreq, " omega = ", pol%dFreqGrid(ifreq), DIMAG(pol%dFreqBrd(ifreq)), " eV"
               enddo
            endif
         endif
         write(*,'(1X,A)')
      endif
    end subroutine initialize_FF

  end subroutine inread_chi

end module inread_chi_m

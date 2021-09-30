#include "f_defs.h"

!!====================================================================
!!
!! Program:
!!
!!    Chi     By Meng Wu (2020)
!!
!!    Chi.x calculates RPA polarizability function and stores it in
!!    chimat.h5 for uniform qgrid, and chi0mat.h5 for q --> 0. 
!!
!! Description:
!!
!!    We will parallel over kpoints and calculate the matrix elements
!!    < c fk | e^{i(q+G).r} | v (fk-q) > with fast Fourier transform
!!    These matrix elements are then assembled to calculate polarizability
!!
!!    Crystal point group symmetry is used to generate wavefunction from
!!    the reduced Brillouin zone to the full Brillouin zone.
!!
!!    We calculate each q-point sequentially. For each q-point, we parallel
!!    over kpoints in the full Brillouin zone.
!!
!! Input files:
!!
!!    chi.inp
!!    WFN and/or WFNmq on a reduced Brillouin zone
!!
!! Output files:
!!
!!    chimat.h5 and/or chi0mat.h5
!!
!!
!!====================================================================

program Chi
  use global_m
  use blas_m
  use fullbz_m
  use misc_m
  use inread_chi_m
  use input_utils_m
  use sort_m
  use mtxel_m
  use gmap_m
  use so32su2_m
  use hdf5
  use h5lt
  use hdf5_io_m
  use epswrite_hdf5_m
  use write_matrix_m
  use scalapack_m
  use lapack_m
  implicit none

  type (crystal) :: crys
  type (polarizability) :: pol
  type (symmetry) :: syms, syms_wfn, symsq_wfn
  type (gspace) :: gvec
  type (grid) :: kg, kgq
  type (kpoints) :: kp, kpq
  type (wavefunction) :: wfnc, wfnv
  type (int_wavefunction) :: intwfnc, intwfnv, intwfnvq
  real(DP) :: tsec(2),tmin(2),tmax(2)
  character*16, allocatable :: routnam(:)
  integer, allocatable :: routsrt(:)
  integer :: error, info, ii, ncount
  integer :: id, ig, nmat, nk_loc_max
  real(DP), allocatable :: ekin(:)
  character(LEN=100) :: filename_chi_hdf5
  INTEGER(HID_T) :: file_id
  integer(HID_T) :: dset_id
  INTEGER(HID_T) :: plist_id
  integer(HID_T) :: group_id
  INTEGER(HID_T) :: filespace
  integer(HID_T) :: memspace
  integer, parameter :: rank_mdat=6
  integer(HSIZE_T) :: count_mdat(6), offset_mdat(6)

  !! mdat(ig, is, iv, ic, ifk_loc)
  complex(DPC), allocatable :: mdat(:,:,:,:,:), vrhoc_noeh_aux(:,:,:,:,:)
  complex(DPC), allocatable :: ph(:), chi_noeh(:,:,:), cg_(:,:,:), cg_2d(:,:), cg_2d_(:,:)
  real(DP), allocatable :: chi_noeh_out(:,:,:,:)
  integer, allocatable :: ifk_map(:,:), irk_map(:,:), gumk_map(:,:,:)
  integer, allocatable :: ind(:), data_temp2(:,:), isort_old2gvec(:), isorti_gvec2old(:)
  integer :: size_temp1(1), size_temp2(2), size_temp6(6)
  real(DP) :: q_vector(3), k_vector(3), fkq_vector(3), rkq_vector(3), lambda
  complex(DPC) :: umtrx(2,2), umtrx_transpose(2,2), prefactor, fact_FF, temp
  integer :: ifk_loc, ifk, ispin, gumk(3), irk_need, irk_start, irk_end, gumk_(3)
  integer :: gumk__(3), rk_blocksize, nrk_loc, nrk_loc_max, nfk_loc, nfk_loc_max
  integer :: request, request_, irk_target, ib, is, ispinor, iq, ipes, ifk_, irk_loc
  integer :: iv, ic, ib_c, ib_v, irk_c, irk_v, ib_vbm, iq_offset, nroutnam, ifreq, imat
  logical :: qpt_done, found  

  info = MPI_INFO_NULL
  call peinfo_init()
  call h5open_f(error)

  !! Read chi.inp
  call inread_chi(pol)
  call MPI_BARRIER(MPI_COMM_WORLD, mpierr)

  call timacc(0,0)
  call timacc(1,1)

  !! Read wavefunctions on the fine grid
  call timacc(2,1)

  call input_chi(pol, crys, gvec, kp, kg, syms, syms_wfn, intwfnv, intwfnc, kpq, &
       kgq, symsq_wfn, intwfnvq)

  call timacc(2,2)

  pol%has_advanced = .false.
  !! for complex wfn case, pol%matrix_flavor = 2
  !! for real wfn case, pol%matrix_flavor = 1
  pol%matrix_flavor = SCALARSIZE
  if ((pol%freq_dep .ne. 0) .and. (pol%matrix_flavor .eq. 1)) then
     call die("When FF, we must use COMPLEX version of BGW.", only_root_writes=.true.)
  endif
  pol%nmatrix = 1
  !! Set number of matrices depending on nspin (chimat is spin resolved)
  !! if output chimat.h5, we enable the spin index
  if (pol%matrix_type .eq. 2) then
     pol%nmatrix = pol%nmatrix * kp%nspin
  endif

  if (pol%matrix_flavor .ne. 2) then
     call die("Only support pol%matrix_flavor = 2.", only_root_writes=.true.)
  endif

  if (pol%nmatrix .ne. 1) then
     call die("Only support pol%nmatrix = 1.", only_root_writes=.true.)
  endif

  prefactor = -4.0D0 / ( DBLE(kg%nf) * crys%celvol * DBLE(kp%nspin) * DBLE(kp%nspinor) )
  ib_vbm = MAXVAL(kp%ifmax(:,:))

  call get_wfn_fftgrid(pol, gvec, kp, intwfnc)
  !! Use pol%WFN_FFTgrid as pol%FFTgrid
  pol%FFTgrid(:) = pol%WFN_FFTgrid(:)

  if (peinf%inode .eq. 0) then
     write(6,'(/1X,A)') 'More job parameters:'
     write(6,'(1X,A,i0)') '- Number of valence bands: ', pol%nvb
     write(6,'(1X,A,i0)') '- Number of valence bands to skip: ', pol%skip_nvb
     write(6,'(1X,A,i0)') '- Number of conduction bands: ', pol%ncb
     write(6,'(1X,A,i0)') '- Number of conduction bands to skip: ', pol%skip_ncb
     write(6,'(1X,A,i0)') '- Number of spins: ', kp%nspin
     write(6,'(1X,A,3I10)') "- FFTgrid = ", pol%FFTgrid
     write(6,'()')
  endif

  rk_blocksize = iceil(kp%nrk, peinf%npes)
  nrk_loc      = NUMROC(kp%nrk, rk_blocksize, peinf%inode, 0, peinf%npes)
  nrk_loc_max  = NUMROC(kp%nrk, rk_blocksize,           0, 0, peinf%npes)
  nfk_loc      = peinf%ikt(peinf%inode+1)
  nfk_loc_max  = peinf%ikt(1)

  SAFE_ALLOCATE(isort_old2gvec, (gvec%ng))
  SAFE_ALLOCATE(isorti_gvec2old, (gvec%ng))
  SAFE_ALLOCATE(pol%isrtx, (gvec%ng))
  SAFE_ALLOCATE(pol%isrtxi, (gvec%ng))
  SAFE_ALLOCATE(ekin, (gvec%ng))
  SAFE_ALLOCATE(pol%nmtx_of_q, (pol%nq))
  if (peinf%inode .eq. 0) then
     do iq = 1, pol%nq
        if (iq <= pol%nq0) then
           call kinetic_energies(gvec, crys%bdot, ekin)
        else
           call kinetic_energies(gvec, crys%bdot, ekin, qvec = pol%qpt(:, iq))
        endif
        call sortrx(gvec%ng, ekin, pol%isrtx, gvec = gvec%components)
        pol%nmtx_of_q(iq) = gcutoff(gvec%ng, ekin, pol%isrtx, pol%ecuts)
     enddo

     if (pol%nq0 > 0) then
        filename_chi_hdf5 = 'chi0mat.h5'
        !! broken time-reversal
        if (.not. pol%time_reversal) then
           if (.not. pol%k_plus_q) then
              filename_chi_hdf5 = "chi0mat.mag.kmq.h5"
           else
              filename_chi_hdf5 = "chi0mat.mag.kpq.h5"
           endif
        endif
        call eps_hdf5_setup_part(kp, gvec, syms, crys, pol, TRUNC(filename_chi_hdf5), 1, &
             pol%nq0, restart = pol%restart)
     endif

     if (pol%nq1 > 0) then
        filename_chi_hdf5 = 'chimat.h5'
        !! broken time-reversal
        if (.not. pol%time_reversal) then
           if (.not. pol%k_plus_q) then
              filename_chi_hdf5 = "chimat.mag.kmq.h5"
           else
              filename_chi_hdf5 = "chimat.mag.kpq.h5"
           endif
        endif
        call eps_hdf5_setup_part(kp, gvec, syms, crys, pol, TRUNC(filename_chi_hdf5), pol%nq0+1, &
             pol%nq, restart = pol%restart)
     endif
  endif
  if (peinf%npes > 1) then
     call MPI_BCAST(pol%nmtx_of_q, pol%nq, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
  endif

  !! Main loop over all q-points
  iq_loop: do iq = 1, pol%nq
     q_vector(:) = pol%qpt(:,iq)

     if (peinf%inode .eq. 0) then
        write(6,'(1X,A)') "======================================================================"
        write(6,'(1X,I5,A,I5,5X,A,3F20.10)') iq, " /", pol%nq, " q = ", q_vector(:)
        write(6,'(1X,A)') "======================================================================"
        write(6,'(1X,A)')
        write(6,'(1X,A,F15.8,A)') "|q| = ", SQRT(DOT_PRODUCT(pol%qpt(:,iq), &
             MATMUL(crys%bdot, pol%qpt(:,iq)))), " Bohr^{-1}"
     endif
     !! |q+G|^2
     if (iq <= pol%nq0) then
        iq_offset = iq
        if (peinf%inode .eq. 0) then
           write(6,'(1X,A)') 'This is the special q->0 point.'
           write(6,'(1X,A)')
        endif
        filename_chi_hdf5 = 'chi0mat.h5'
        !! broken time-reversal
        if (.not. pol%time_reversal) then
           if (.not. pol%k_plus_q) then
              filename_chi_hdf5 = "chi0mat.mag.kmq.h5"
           else
              filename_chi_hdf5 = "chi0mat.mag.kpq.h5"
           endif
        endif
        call kinetic_energies(gvec, crys%bdot, ekin)
     else
        iq_offset = iq - pol%nq0
        if (peinf%inode .eq. 0) then
           write(6,'(1X,A)') 'This is a regular non-zero q-point.'
           write(6,'(1X,A)')
        endif
        filename_chi_hdf5 = 'chimat.h5'
        !! broken time-reversal
        if (.not. pol%time_reversal) then
           if (.not. pol%k_plus_q) then
              filename_chi_hdf5 = "chimat.mag.kmq.h5"
           else
              filename_chi_hdf5 = "chimat.mag.kpq.h5"
           endif
        endif
        call kinetic_energies(gvec, crys%bdot, ekin, qvec = pol%qpt(:,iq))
     endif

     if (pol%restart) then
        if (peinf%inode .eq. 0) qpt_done = is_qpt_done(TRUNC(filename_chi_hdf5), iq_offset)
#ifdef MPI
        if (peinf%npes > 1) then
           call MPI_BCAST(qpt_done, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, mpierr)
        endif
#endif
        if (qpt_done) then
           if (peinf%inode==0) then
              write(6,'(/,1X,A,/)') 'This q-point was already calculated: skipping.'
              write(6,'(1X,A)')
           endif
           cycle iq_loop
        endif
     endif

     call sortrx(gvec%ng, ekin, pol%isrtx, gvec = gvec%components)
     pol%isrtxi = 0
     do ig = 1, gvec%ng
        pol%isrtxi(pol%isrtx(ig)) = ig
     enddo

     pol%nmtx = pol%nmtx_of_q(iq)

     if (peinf%inode .eq. 0) then
        call write_gvec_indices_hdf(gvec%ng, pol%isrtx, pol%isrtxi, ekin, iq_offset, TRUNC(filename_chi_hdf5))
        write(6, '(1X,A,i0)') 'Rank of the polarizability matrix (nmtx): ', pol%nmtx
        write(6,'(1X,A)')
     endif

     if (iq > pol%nq0) then
        !! Build global maps for wfnv at this iq
        SAFE_ALLOCATE(irk_map, (peinf%npes, nfk_loc_max))
        SAFE_ALLOCATE(ifk_map, (peinf%npes, nfk_loc_max))
        SAFE_ALLOCATE(gumk_map, (3, peinf%npes, nfk_loc_max))
        irk_map = 0
        ifk_map = 0
        gumk_map = 0
        do ipes = 1, peinf%npes
           do ifk_loc = 1, peinf%ikt(ipes)
              ifk = peinf%ik(ipes, ifk_loc)

              !! Need to modify for broken TRS
              if (pol%k_plus_q) then
                 k_vector(:) = kg%f(:, ifk) + q_vector(:)
              else
                 k_vector(:) = kg%f(:, ifk) - q_vector(:)
              endif
              call get_gumk3(crys%bdot, k_vector, gumk)
              k_vector(:) = k_vector(:) - DBLE(gumk(:))
              !! Find fk -/+ q in kg%f
              found = .false.
              ifk_loop: do ifk_ = 1, kg%nf
                 if (NORM2(k_vector(:) - kg%f(:, ifk_)) < TOL_SMALL) then
                    found = .true.
                    ifk_map(ipes, ifk_loc) = ifk_          !! refers to kg%f
                    irk_map(ipes, ifk_loc) = kg%indr(ifk_) !! refers to kp%rk or kg%r
                    !! R[rkq] + kg0 = \lfloor fk - q \rfloor = fk - q - gumk
                    !! We need fk - q = R[rkq] + kg0 + gumk = \lfloor fk - q \rfloor + gumk
                    ! gumk_map(:, ipes, ifk_loc) = gumk(:) + kg%kg0(:,ifk_)
                    gumk_map(:, ipes, ifk_loc) = gumk(:)
                    exit ifk_loop
                 endif
              enddo ifk_loop
              if (.not. found) then
                 write(6,'(A,I5,A,3F15.7)') "ifk = ", ifk, " fk - q - gumk = ", k_vector(:)
                 call die("Cannot find fk-q.", only_root_writes=.true.)
              endif
           enddo
        enddo
     endif !! iq > pol%nq0

     if (nfk_loc > 0) then
        SAFE_ALLOCATE(mdat, (pol%nmtx, kp%nspin, pol%nvb, pol%ncb, nfk_loc))
     else
        SAFE_ALLOCATE(mdat, (1, 1, 1, 1, 1))
     endif
     SAFE_ALLOCATE(chi_noeh, (pol%nmtx, pol%nmtx, pol%nFreq))
     mdat = ZERO
     chi_noeh = ZERO

     do ifk_loc = 1, nfk_loc_max
        call timacc(3,1)
        !! All procs work together to generate wfnv
        if (iq > pol%nq0) then
           SAFE_ALLOCATE(cg_, (kp%ngkmax, pol%nvb, kp%nspin*kp%nspinor))
           cg_ = ZERO
           isort_old2gvec = 0

           !! If current proc stores some of intwfnv
           if (nrk_loc > 0) then
              irk_start = INDXL2G(      1, rk_blocksize, peinf%inode, 0, peinf%npes)
              irk_end   = INDXL2G(nrk_loc, rk_blocksize, peinf%inode, 0, peinf%npes)
           else
              irk_start = 0
              irk_end = 0
           endif

           !! current proc is not idle
           if (ifk_loc <= nfk_loc) then
              !! if the rk needed by current proc is elsewhere, MPI_IRECV
              !! if the rk needed by current proc is just here, copy intwfnv to wfnv
              irk_need = irk_map(peinf%inode+1, ifk_loc)
              !! Current proc has the irk_need
              if ((irk_need <= irk_end) .and. (irk_need >= irk_start)) then
                 irk_loc = INDXG2L(irk_need, rk_blocksize, 0, 0, peinf%npes)
                 !! intwfnv%cgk(ig, ib, is, irk_loc)
                 cg_(:, :, :) = intwfnv%cgk(:, :, :, irk_loc)
                 isort_old2gvec(:) = intwfnv%isort(:, irk_loc)
              else
                 call MPI_IRECV(cg_(1,1,1), kp%ngkmax*pol%nvb*kp%nspin*kp%nspinor, MPI_SCALAR, MPI_ANY_SOURCE, &
                      irk_need, MPI_COMM_WORLD, request, mpierr)
                 !! tags must be non-negative
                 call MPI_IRECV(isort_old2gvec(1), gvec%ng, MPI_INTEGER, MPI_ANY_SOURCE, kp%nrk+irk_need, &
                      MPI_COMM_WORLD, request_, mpierr)
              endif
           endif

           !! loop over all procs, see if current proc needs to send out epscol to a target proc
           !! All procs take care of the ipes-th procs
           do ipes = 1, peinf%npes
              !! ipes-th proc is not idle
              if (peinf%ik(ipes, ifk_loc) .ne. 0) then
                 irk_target = irk_map(ipes, ifk_loc)
                 !! See if ipes-th proc needs intwfnv from current procs
                 if ((irk_target <= irk_end) .and. (irk_target >= irk_start)) then
                    irk_loc = INDXG2L(irk_target, rk_blocksize, 0, 0, peinf%npes)
                    if (ipes .ne. (peinf%inode+1)) then
                       call MPI_SEND(intwfnv%cgk(1,1,1,irk_loc), kp%ngkmax*pol%nvb*kp%nspin*kp%nspinor, MPI_SCALAR, &
                            ipes-1, irk_target, MPI_COMM_WORLD, mpierr)
                       call MPI_SEND(intwfnv%isort(1,irk_loc), gvec%ng, MPI_INTEGER, ipes-1, kp%nrk+irk_target, &
                            MPI_COMM_WORLD, mpierr)
                    endif
                 endif
              endif
           enddo

           !! Once current procs receive the epsinv_loc_rq it needs, it is free to move on
           if (ifk_loc <= nfk_loc) then
              !! Current proc does not have the irq_need
              if ( .not. ((irk_need <= irk_end) .and. (irk_need >= irk_start)) ) then
                 call MPI_WAIT(request, MPI_STATUS_IGNORE, mpierr)
                 call MPI_WAIT(request_, MPI_STATUS_IGNORE, mpierr)
              endif
              wfnv%nband = pol%nvb
              wfnv%ng = kp%ngk(irk_need)
              wfnv%nspin = kp%nspin
              wfnv%nspinor = kp%nspinor

              ifk_ = ifk_map(peinf%inode+1, ifk_loc)
              gumk_(:) = gumk_map(:, peinf%inode+1, ifk_loc)
              !! fkq_vector could be outside FBZ
              fkq_vector = kg%f(:, ifk_) + DBLE(gumk_(:))
              rkq_vector = kg%r(:, irk_need)

              SAFE_ALLOCATE(wfnv%isort, (gvec%ng))
              !! intwfnv%cgk(ig, ib, is, irk)
              SAFE_ALLOCATE(wfnv%cg, (wfnv%ng, wfnv%nband, wfnv%nspin*wfnv%nspinor))

              call kinetic_energies(gvec, crys%bdot, ekin, qvec = fkq_vector)
              !! wfnv%isort(ig_new) = ig_gvec
              call sortrx(gvec%ng, ekin, wfnv%isort, gvec = gvec%components)

              isorti_gvec2old = 0
              do ig = 1, wfnv%ng
                 isorti_gvec2old(isort_old2gvec(ig)) = ig
              enddo

              !! Rotate cg_ into wfnv%cg
              SAFE_ALLOCATE(ind, (wfnv%ng))
              SAFE_ALLOCATE(ph,  (wfnv%ng))
              ind = 0
              ph  = ZERO
              gumk__(:) = gumk_(:) + kg%kg0(:,ifk_)
              call gmap(gvec, syms, wfnv%ng, kg%itran(ifk_), gumk__, wfnv%isort, isorti_gvec2old, ind, ph)

              !$OMP PARALLEL DO collapse(4)
              do is = 1, wfnv%nspin
                 do ispinor = 1, wfnv%nspinor
                    do ib = 1, wfnv%nband
                       do ig = 1, wfnv%ng
                          if (ind(ig) > 0) then
                             wfnv%cg(ig, ib, is*ispinor) = ph(ig) * cg_(ind(ig), ib, is*ispinor)
                          else
                             wfnv%cg(ig, ib, is*ispinor) = ZERO
                          endif
                       enddo ! ig
                    enddo ! ib
                 enddo ! ispinor
              enddo ! is
              !$OMP END PARALLEL DO

              SAFE_DEALLOCATE(ind)
              SAFE_DEALLOCATE(ph)
#ifdef CPLX
              if (wfnv%nspinor .eq. 2) then
                 SAFE_ALLOCATE(cg_2d,  (wfnv%ng*pol%nvb, wfnv%nspin*wfnv%nspinor))
                 SAFE_ALLOCATE(cg_2d_, (wfnv%ng*pol%nvb, wfnv%nspin*wfnv%nspinor))

                 call so32su2(syms%mtrx_cart(1:3, 1:3, kg%itran(ifk_)), umtrx)

                 umtrx_transpose = TRANSPOSE(umtrx)
                 !$OMP PARALLEL WORKSHARE
                 cg_2d   = RESHAPE(wfnv%cg, (/wfnv%ng*pol%nvb, wfnv%nspin*wfnv%nspinor/))
                 cg_2d_  = MATMUL(cg_2d, umtrx_transpose)
                 wfnv%cg = RESHAPE(cg_2d_, (/wfnv%ng, pol%nvb, wfnv%nspin*wfnv%nspinor/))
                 !$OMP END PARALLEL WORKSHARE

                 SAFE_DEALLOCATE(cg_2d)
                 SAFE_DEALLOCATE(cg_2d_)
              endif
#endif
              !! Checknorm
              do ib = 1, wfnv%nband
                 call checknorm('wfnv', ib, ifk_, wfnv%nspin, wfnv%cg(1:wfnv%ng, ib, :))
              enddo
           endif !! ifk_loc <= nfk_loc
           SAFE_DEALLOCATE(cg_)
        endif !! iq > pol%nq0
        call timacc(3,2)

        !! Local proc generate wfnc and wfnvq
        if (ifk_loc <= nfk_loc) then
           ifk = peinf%ik(peinf%inode+1, ifk_loc)
           call timacc(3,1)
           call genwf_chi(pol, ifk, iq, crys, gvec, kg, kgq, syms_wfn, symsq_wfn, wfnc, wfnv, intwfnc, intwfnvq)
           call timacc(3,2)

           call timacc(4,1)
           do ispin = 1, kp%nspin
              !! Use wfnc and wfnv to calculate M matrix elements for all valence states,
              !! all conduction states, and all Gvectors determined by the epsilon cutoff
              !! Here mdat = crhov_noeh
              !! If k_plus_q = F, mdat = < c, k   | e^{i(q+G).r} | v, k-q >
              !! If k_plus_q = T, mdat = < v, k+q | e^{i(q+G).r} | c, k   >
              call mtxel(gvec, pol, ispin, wfnc, wfnv, mdat_complex=mdat(:,:,:,:,ifk_loc))
           enddo

           SAFE_DEALLOCATE_P(wfnv%cg)
           SAFE_DEALLOCATE_P(wfnv%isort)
           SAFE_DEALLOCATE_P(wfnc%cg)
           SAFE_DEALLOCATE_P(wfnc%isort)

           !! GPP
           if (pol%freq_dep .eq. 0) then
              !$OMP PARALLEL DO collapse(3) PRIVATE(ib_v, ib_c, irk_v, irk_c, lambda)
              do ic = 1, pol%ncb
                 do iv = 1, pol%nvb
                    do is = 1, kp%nspin
                       ib_v  = ib_vbm - iv + 1
                       ib_c = ib_vbm + ic
                       irk_c = kg%indr(ifk)
                       !! In units of Ryd
                       if (iq <= pol%nq0) then
                          irk_v = kgq%indr(pol%indexq(ifk))
                          lambda = kp%el(ib_c, irk_c, is) - kpq%el(ib_v, irk_v, is)
                       else
                          irk_v = kg%indr(ifk_)
                          lambda = kp%el(ib_c, irk_c, is) - kp%el(ib_v, irk_v, is)
                       endif
                       ! !! mdat(ig, is, iv, ic, ifk_loc)
                       mdat(:, is, iv, ic, ifk_loc) = CONJG(mdat(:, is, iv, ic, ifk_loc)) / SQRT(lambda)
                    enddo ! is
                 enddo ! iv
              enddo ! ic
              !$OMP END PARALLEL DO

              !! FF
           else
              !$OMP PARALLEL WORKSHARE
              mdat(:, :, :, :, ifk_loc) = CONJG(mdat(:, :, :, :, ifk_loc))
              !$OMP END PARALLEL WORKSHARE
           endif
           !! After the complex conjugate:
           !! if pol%k_plus_q = F, mdat = < v, k-q | e^{-i(q+G).r} | c, k   >
           !! if pol%k_plus_q = T, mdat = < c, k   | e^{-i(q+G).r} | v, k+q >
           call timacc(4,2)
        else
           ifk = -1
        endif
     enddo !! ifk_loc

     call timacc(5,1)
     !! DOT PRODUCT of mdat to get chi_noeh
     if (nfk_loc > 0) then
        nmat = kp%nspin*pol%nvb*pol%ncb*nfk_loc
        !! GPP
        if (pol%freq_dep .eq. 0) then
           do is = 1, kp%nspin
              if (pol%time_reversal) then
                 !! mdat = vrhoc_noeh
                 call zgemm('n', 'c', pol%nmtx, pol%nmtx, nmat, prefactor, mdat(1,1,1,1,1), pol%nmtx, &
                      mdat(1,1,1,1,1), pol%nmtx, ZERO, chi_noeh(:,:,1), pol%nmtx)
                 !! if use broken-time-reversal formalism, only have contributions
                 !! from either <v(k-q)|e^{-i(q+G).r}|ck> or <ck|e^{-i(q+G).r}|v(k+q) >
                 !! which means the prefactor is only half of the time-reversal case ==> prefactor/2.0D0
              else
                 call zgemm('n', 'c', pol%nmtx, pol%nmtx, nmat, prefactor/2.0D0, mdat(1,1,1,1,1), &
                      pol%nmtx, mdat(1,1,1,1,1), pol%nmtx, ZERO, chi_noeh(:,:,1), pol%nmtx)
              endif
           enddo
           !! FF
        else
           SAFE_ALLOCATE(vrhoc_noeh_aux, (pol%nmtx, kp%nspin, pol%nvb, pol%ncb, nfk_loc))
           !! real-axis frequencies
           !! Here varepsilon = pol%dFreqBrd
           do ifreq = 1, pol%nfreq_real
              vrhoc_noeh_aux = ZERO
              !$OMP PARALLEL DO collapse(4) PRIVATE(ifk, ifk_, ib_v, ib_c, irk_c, irk_v, lambda, fact_FF)
              do ifk_loc = 1, nfk_loc
                 do ic = 1, pol%ncb
                    do iv = 1, pol%nvb
                       do is = 1, kp%nspin
                          ifk = peinf%ik(peinf%inode+1, ifk_loc)
                          if (iq > pol%nq0) then
                             ifk_ = ifk_map(peinf%inode+1, ifk_loc)
                          endif
                          ib_v  = ib_vbm - iv + 1
                          ib_c = ib_vbm + ic
                          irk_c = kg%indr(ifk)
                          !! In units of Ryd
                          if (iq <= pol%nq0) then
                             irk_v = kgq%indr(pol%indexq(ifk))
                             lambda = kp%el(ib_c, irk_c, is) - kpq%el(ib_v, irk_v, is)
                          else
                             irk_v = kg%indr(ifk_)
                             lambda = kp%el(ib_c, irk_c, is) - kp%el(ib_v, irk_v, is)
                          endif
                          !! Retarded
                          if (.not. pol%timeordered) then
                             !! With time-reversal
                             if (pol%time_reversal) then
                                fact_FF = -0.5D0 * ( 1.0D0/(pol%dFreqGrid(ifreq)/ryd + pol%dFreqBrd(ifreq)/ryd - lambda) &
                                     - 1.0D0/(pol%dFreqGrid(ifreq)/ryd + pol%dFreqBrd(ifreq)/ryd + lambda))
                                !! Broken time-reversal
                             else
                                !! Use |k-q>
                                if (.not. pol%k_plus_q) then
                                   fact_FF = -0.5D0/(pol%dFreqGrid(ifreq)/ryd + pol%dFreqBrd(ifreq)/ryd - lambda)
                                else
                                   fact_FF =  0.5D0/(pol%dFreqGrid(ifreq)/ryd + pol%dFreqBrd(ifreq)/ryd + lambda)
                                endif
                             endif
                             !! Time-ordered
                          else
                             !! With time-reversal
                             if (pol%time_reversal) then
                                fact_FF = -0.5D0 * ( 1.0D0/(pol%dFreqGrid(ifreq)/ryd - lambda + pol%dFreqBrd(ifreq)/ryd) &
                                     - 1.0D0/(pol%dFreqGrid(ifreq)/ryd + lambda - pol%dFreqBrd(ifreq)/ryd))
                                !! Broken time-reversal
                             else
                                !! Use |k-q>, the first term in Eqn. (213), Tips in Field Theory
                                if (.not. pol%k_plus_q) then
                                   fact_FF = -0.5D0/(pol%dFreqGrid(ifreq)/ryd - lambda + pol%dFreqBrd(ifreq)/ryd)
                                   !! Use |k+q>, the second term in Eqn. (213), Tips in Field Theory
                                else
                                   fact_FF =  0.5D0/(pol%dFreqGrid(ifreq)/ryd + lambda - pol%dFreqBrd(ifreq)/ryd)
                                endif
                             endif
                          endif
                          vrhoc_noeh_aux(:, is, iv, ic, ifk_loc) = fact_FF * mdat(:, is, iv, ic, ifk_loc)
                       enddo ! is
                    enddo ! iv
                 enddo ! ic
              enddo !! ifk_loc
              !$OMP END PARALLEL DO

              !! Calculate chi_noeh
              do is = 1, kp%nspin
                 call zgemm('n','c', pol%nmtx, pol%nmtx, nmat, prefactor, vrhoc_noeh_aux(1,1,1,1,1), pol%nmtx, mdat(1,1,1,1,1), &
                      pol%nmtx, ZERO, chi_noeh(:,:,ifreq), pol%nmtx)
              enddo
           enddo !! ifreq_real

           !! imaginary-axis frequencies
           !! varepsilon = 0
           do ifreq = pol%nfreq_real+1, pol%nfreq
              vrhoc_noeh_aux = ZERO
              !$OMP PARALLEL DO collapse(4) PRIVATE(ifk, ifk_, ib_v, ib_c, irk_c, irk_v, lambda, fact_FF)
              do ifk_loc = 1, nfk_loc
                 do ic = 1, pol%ncb
                    do iv = 1, pol%nvb
                       do is = 1, kp%nspin
                          ifk = peinf%ik(peinf%inode+1, ifk_loc)
                          if (iq > pol%nq0) then
                             ifk_ = ifk_map(peinf%inode+1, ifk_loc)
                          endif
                          ib_v  = ib_vbm - iv + 1
                          ib_c = ib_vbm + ic
                          irk_c = kg%indr(ifk)
                          !! In units of Ryd
                          if (iq <= pol%nq0) then
                             irk_v = kgq%indr(pol%indexq(ifk))
                             lambda = kp%el(ib_c, irk_c, is) - kpq%el(ib_v, irk_v, is)
                          else
                             irk_v = kg%indr(ifk_)
                             lambda = kp%el(ib_c, irk_c, is) - kp%el(ib_v, irk_v, is)
                          endif
                          !! Retarded
                          !! With time-reversal
                          if (pol%time_reversal) then
                             fact_FF = -0.5D0 * ( 1.0D0/(pol%dFreqGrid(ifreq)/ryd + pol%dFreqBrd(ifreq)/ryd - lambda) &
                                  - 1.0D0/(pol%dFreqGrid(ifreq)/ryd + pol%dFreqBrd(ifreq)/ryd + lambda))
                             !! Broken time-reversal
                          else
                             !! Use |k-q>
                             if (.not. pol%k_plus_q) then
                                fact_FF = -0.5D0/(pol%dFreqGrid(ifreq)/ryd + pol%dFreqBrd(ifreq)/ryd - lambda)
                             else
                                fact_FF =  0.5D0/(pol%dFreqGrid(ifreq)/ryd + pol%dFreqBrd(ifreq)/ryd + lambda)
                             endif
                          endif
                          vrhoc_noeh_aux(:, is, iv, ic, ifk_loc) = fact_FF * mdat(:, is, iv, ic, ifk_loc)
                       enddo ! is
                    enddo ! iv
                 enddo ! ic
              enddo ! ifk_loc
              !$OMP END PARALLEL DO

              !! Calculate chi_noeh
              do is = 1, kp%nspin
                 call zgemm('n','c', pol%nmtx, pol%nmtx, nmat, prefactor, vrhoc_noeh_aux(1,1,1,1,1), pol%nmtx, &
                      mdat(1,1,1,1,1), pol%nmtx, ZERO, chi_noeh(:,:,ifreq), pol%nmtx)
              enddo
           enddo !! ifreq_imag
           SAFE_DEALLOCATE(vrhoc_noeh_aux)
        endif
     endif
     SAFE_DEALLOCATE(mdat)

     if (iq > pol%nq0) then
        SAFE_DEALLOCATE(irk_map)
        SAFE_DEALLOCATE(ifk_map)
        SAFE_DEALLOCATE(gumk_map)
     endif

     call MPI_ALLREDUCE(MPI_IN_PLACE, chi_noeh(1,1,1), pol%nmtx*pol%nmtx*pol%nFreq, MPI_COMPLEX_DPC, &
          MPI_SUM, MPI_COMM_WORLD, mpierr)

     call timacc(5,2)

     call timacc(6,1)
     if (peinf%inode .eq. 0) then
        !! use chi_noeh_out(2, nmtx_max, nmtx_max, nfreq) instead of chi_noeh_out(2, pol%nmtx, pol%nmtx, pol%nfreq)
        !! To save memory, we can write one frequency at a time!
        SAFE_ALLOCATE(chi_noeh_out, (2, pol%nmtx, pol%nmtx, pol%nFreq))
        chi_noeh_out(1,:,:,:) = DBLE(chi_noeh(:,:,:))
        chi_noeh_out(2,:,:,:) = DIMAG(chi_noeh(:,:,:))

        call h5fopen_f(TRUNC(filename_chi_hdf5), H5F_ACC_RDWR_F, file_id, error)
        if (error .ne. 0) then
           call die("HDF5 error", only_root_writes=.true.)
        endif

        !! Write out chi_noeh_out using hyperslab
        call hdf5_write_double_hyperslab(file_id, 'mats/matrix', (/pol%matrix_flavor, pol%nmtx, pol%nmtx, pol%nfreq, &
             pol%nmatrix, 1/), (/ 0, 0, 0, 0, 0, iq_offset-1/), chi_noeh_out, error)
        call h5fclose_f(file_id, error)
        
        call set_qpt_done(TRUNC(filename_chi_hdf5), iq_offset)
        SAFE_DEALLOCATE(chi_noeh_out)

        if (pol%freq_dep .eq. 0) then
           write(6,'(1X,A,3I5,A,3I5,A,2F15.8)') "chi_noeh(  G = ",gvec%components(:,1),", Gp = ",gvec%components(:,1), &
                " ) = ", chi_noeh(pol%isrtxi(1),pol%isrtxi(1),1)
           write(6,'(1X,A,3I5,A,3I5,A,2F15.8)') "chi_noeh(  G = ",gvec%components(:,2),", Gp = ",gvec%components(:,2), &
                " ) = ", chi_noeh(pol%isrtxi(2),pol%isrtxi(2),1)
           write(6,'(1X,A,3I5,A,3I5,A,2F15.8)') "chi_noeh(  G = ",gvec%components(:,1),", Gp = ",gvec%components(:,2), &
                " ) = ", chi_noeh(pol%isrtxi(1),pol%isrtxi(2),1)
        else
           !! Time-ordered
           do ifreq = 1, pol%nfreq_real
              if (pol%timeordered) then
                 write(6,'(A,I5,A,F12.5,A,F12.5,A)') "Frequency #", ifreq, " : ", pol%dFreqGrid(ifreq), " eV varepsilon = ", &
                      DIMAG(pol%dFreqBrd(ifreq))," eV"
              else
                 write(6,'(A,I5,A,"(",F12.5,",",F12.5,")",A)') "Frequency #", ifreq, " : ", &
                      pol%dFreqGrid(ifreq)+pol%dFreqBrd(ifreq), " eV varepsilon = 0 eV"
              endif
              write(6,'(1X,A,3I5,A,3I5,A,2F15.8)') "chi_noeh(  G = ",gvec%components(:,1),", Gp = ",gvec%components(:,1), &
                   " ) = ", chi_noeh(pol%isrtxi(1),pol%isrtxi(1),ifreq)
              write(6,'(1X,A,3I5,A,3I5,A,2F15.8)') "chi_noeh(  G = ",gvec%components(:,2),", Gp = ",gvec%components(:,2), &
                   " ) = ", chi_noeh(pol%isrtxi(2),pol%isrtxi(2),ifreq)
              write(6,'(1X,A,3I5,A,3I5,A,2F15.8)') "chi_noeh(  G = ",gvec%components(:,1),", Gp = ",gvec%components(:,2), &
                   " ) = ", chi_noeh(pol%isrtxi(1),pol%isrtxi(2),ifreq)
           enddo

           !! Imaginary frequencies
           do ifreq = pol%nfreq_real+1, pol%nfreq
              write(6,'(A,I5,A,"(",F12.5,",",F12.5,")",A)') "Frequency #", ifreq, " : ", &
                   pol%dFreqGrid(ifreq)+pol%dFreqBrd(ifreq), " eV varepsilon = 0 eV"
              write(6,'(1X,A,3I5,A,3I5,A,2F15.8)') "chi_noeh(  G = ",gvec%components(:,1),", Gp = ",gvec%components(:,1), &
                   " ) = ", chi_noeh(pol%isrtxi(1),pol%isrtxi(1),ifreq)
              write(6,'(1X,A,3I5,A,3I5,A,2F15.8)') "chi_noeh(  G = ",gvec%components(:,2),", Gp = ",gvec%components(:,2), &
                   " ) = ", chi_noeh(pol%isrtxi(2),pol%isrtxi(2),ifreq)
              write(6,'(1X,A,3I5,A,3I5,A,2F15.8)') "chi_noeh(  G = ",gvec%components(:,1),", Gp = ",gvec%components(:,2), &
                   " ) = ", chi_noeh(pol%isrtxi(1),pol%isrtxi(2),ifreq)
           enddo
        endif
     endif !! ROOT
     call timacc(6,2)
     SAFE_DEALLOCATE(chi_noeh)
  enddo iq_loop !! iq

  SAFE_DEALLOCATE_P(pol%isrtx)
  SAFE_DEALLOCATE_P(pol%isrtxi)
  SAFE_DEALLOCATE_P(pol%nmtx_of_q)
  SAFE_DEALLOCATE_P(gvec%components)
  SAFE_DEALLOCATE_P(gvec%index_vec)
  SAFE_DEALLOCATE_P(pol%qpt)
  SAFE_DEALLOCATE(ekin)
  SAFE_DEALLOCATE(isort_old2gvec)
  SAFE_DEALLOCATE(isorti_gvec2old)

#ifdef MPI
  call MPI_BARRIER(MPI_COMM_WORLD, mpierr)
#endif

  !! Time accounting
  nroutnam = 6
  SAFE_ALLOCATE(routnam, (nroutnam))

  routnam(1)='TOTAL:'
  routnam(2)='INPUT_CHI:'
  routnam(3)='GENWF:'
  routnam(4)='GET_MDAT:'
  routnam(5)='GET_CHI:'
  routnam(6)='OUTPUT_CHI:'

  SAFE_ALLOCATE(routsrt, (nroutnam))
  routsrt=(/ (ii, ii=2,nroutnam), 1 /)

  call timacc(1,2)

  if(peinf%inode.eq.0) then
     write(6,*)
     write(6,9000) 'CPU (s)','WALL (s)','#'
     write(6,*)
  endif

  do ii = 1, nroutnam
     call timacc(routsrt(ii), 3, tsec, ncount)
#ifdef MPI
     call MPI_ALLREDUCE(tsec, tmin, 2, MPI_REAL_DP, MPI_MIN, MPI_COMM_WORLD, mpierr)
     call MPI_ALLREDUCE(tsec, tmax, 2, MPI_REAL_DP, MPI_MAX, MPI_COMM_WORLD, mpierr)
#else
     tmin = tsec
     tmax = tsec
#endif
     if (peinf%inode .eq. 0) then
        write(6,9001) routnam(routsrt(ii)),tmin(1),tmin(2),ncount
        write(6,9002) tsec(1),tsec(2)
        write(6,9003) tmax(1),tmax(2)
     endif
  enddo

9000 format(23x,a13,3x,a13,3x,a8)
9001 format(1X,A16,'(min.)',f13.3,3x,f13.3,3x,i8)
9002 format(   17x,'(PE 0)',f13.3,3x,f13.3)
9003 format(   17x,'(max.)',f13.3,3x,f13.3)

  call h5close_f(error)

#ifdef MPI
  call MPI_FINALIZE(mpierr)
#endif

contains

  !! Find the minimum fftbox that holds all the WFNs
  subroutine get_wfn_fftgrid(pol, gvec, kp, intwfn)
    type(polarizability), intent(inout) :: pol
    type(gspace), intent(in) :: gvec
    type(kpoints), intent(in) :: kp
    type(int_wavefunction), intent(in) :: intwfn
    integer :: ik, wfn_box_min(3), wfn_box_max(3)
    integer :: wfn_box_min_(3), wfn_box_max_(3), total_gvec_grid, total_pol_WFN_grid
    PUSH_SUB(get_wfn_fftgrid)

    wfn_box_min(:) = 0; wfn_box_max(:) = 0
    wfn_box_min_(:) = 0; wfn_box_max_(:) = 0

    ! int_wavefunction in epsilon.x will store sorting indices for all the RBZ kpoints
    do ik = 1, peinf%nrk(peinf%inode+1) ! kp%nrk
       call get_gvecs_bounds(gvec, intwfn%ng(ik), intwfn%isort(:,ik), wfn_box_min, wfn_box_max)
    enddo

    call MPI_ALLREDUCE(wfn_box_max(1), wfn_box_max_(1), 3, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD,mpierr)
    call MPI_ALLREDUCE(wfn_box_min(1), wfn_box_min_(1), 3, MPI_INTEGER, MPI_MIN, MPI_COMM_WORLD,mpierr)

    wfn_box_max(:) = wfn_box_max_(:)
    wfn_box_min(:) = wfn_box_min_(:)
    pol%WFN_FFTgrid(1:3) = wfn_box_max(1:3) - wfn_box_min(1:3) + 1

    !! Use gvec%FFTgrid to correct pol%WFN_FFTgrid, to make sure that the ratio betwen gvec%FFTgrid(:) is preseved
    total_gvec_grid    = SUM(gvec%FFTgrid(:))
    total_pol_WFN_grid = SUM(pol%WFN_FFTgrid(:))
    pol%WFN_FFTgrid(1) = NINT(DBLE(gvec%FFTgrid(1)) / DBLE(Total_gvec_grid) * DBLE(total_pol_WFN_grid))
    pol%WFN_FFTgrid(2) = NINT(DBLE(gvec%FFTgrid(2)) / DBLE(Total_gvec_grid) * DBLE(total_pol_WFN_grid))
    pol%WFN_FFTgrid(3) = NINT(DBLE(gvec%FFTgrid(3)) / DBLE(Total_gvec_grid) * DBLE(total_pol_WFN_grid))

    if (peinf%verb_debug .and. peinf%inode==0) then
       write(6,*) 'WFN min. FFT grid:',pol%WFN_FFTgrid
    endif

    POP_SUB(get_wfn_fftgrid)
    return
  end subroutine get_wfn_fftgrid

  subroutine get_gvecs_bounds(gvec, ng, isort, box_min, box_max)
    type(gspace), intent(in) :: gvec
    integer, intent(in) :: ng
    integer, intent(in) :: isort(:)
    integer, intent(inout) :: box_min(3), box_max(3)
    integer :: ig
    PUSH_SUB(get_gvecs_bounds)

    do ig = 1, ng
       box_min(1:3) = MIN(box_min(1:3),  gvec%components(1:3, isort(ig)))
       box_max(1:3) = MAX(box_max(1:3),  gvec%components(1:3, isort(ig)))
    enddo

    POP_SUB(get_gvecs_bounds)
    return
  end subroutine get_gvecs_bounds

end program Chi

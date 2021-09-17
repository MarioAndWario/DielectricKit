!============================================================================
!
! Included from fftw.F90
!
!============================================================================

!> This routine takes data(1:ndata) and puts it into the FFT box fftbox(:,:,:).
!! The FFT box is zeroed out first, and the data is entered into it.
!! 
!!   ndata -- number of data items in data(:)
!!   data -- the data set, real or complex, depending on ifdef CPLX
!!   ng -- number of g vectors in glist
!!   glist -- a master list of g vectors
!!   gindex(1:ng) -- which g vector (in the master list) the data(1:ndata)
!!                   actually refer to:  so data(j) is for the g-vector
!!                   glist(1:3,gindex(j))
!!   fftbox(:,:,:) -- 3D complex FFT box where the data is put
!!   Nfft(1:3) -- sizes of FFT box Nx,Ny,Nz
subroutine X(put_into_fftbox)(ndata, data, glist, gindex, fftbox, Nfft)
  integer, intent(in) :: ndata
  SCALAR,  intent(in) :: data(:) !< (ndata) this is to avoid creation of array temporary
  integer, intent(in) :: glist(:,:) !< (3, ng)
  integer, intent(in) :: gindex(:) !< (ng)
  integer, intent(in) :: Nfft(:) !< (3)
  complex(DPC), intent(out) :: fftbox(:,:,:) !< (Nfft(1), Nfft(2), Nfft(3))
  
  integer :: j, k, bidx(3)
  
  PUSH_SUB(X(put_into_fftbox))
  
  ! Zero out FFT box and put data into it
  if(peinf%inode.eq.0) call timacc(91,1)

!$OMP PARALLEL
!$OMP DO COLLAPSE(2)
  do j=1,Nfft(3)
    do k=1,Nfft(2)
      fftbox(:,k,j) = (0.0d0,0.0d0)
    enddo
  enddo
!$OMP END DO
!$OMP END PARALLEL

  if(peinf%inode.eq.0) call timacc(91,2)
  if(peinf%inode.eq.0) call timacc(92,1)

!$OMP PARALLEL PRIVATE(bidx,j) SHARED(fftbox, glist, gindex, data)
!$OMP DO
  do j=1,ndata
    call gvec_to_fft_index(glist(:,gindex(j)),bidx,Nfft)
    fftbox(bidx(1),bidx(2),bidx(3)) = data(j)
  end do
!$OMP END DO
!$OMP END PARALLEL

  if(peinf%inode.eq.0) call timacc(92,2)

  POP_SUB(X(put_into_fftbox))

  return
end subroutine X(put_into_fftbox)

!> Does the inverse of the above routine:  takes the data in the
!! fftbox(:,:,:) and puts it into the data(1:ndata) array.  ndata entries
!! are extracted, and the gindex and glist specify which ones to get:
!! data(j) corresponds to the g-vector glist(:,gindex(j)).  The data
!! in fftbox is multiplied by scale before storage into data(:).
!!
!! data(:) is zeroed first and then the data is put into it.
!!
subroutine X(get_from_fftbox)(ndata, data, glist, gindex, fftbox, Nfft, scale)
  integer, intent(in) :: ndata
  SCALAR, intent(out) :: data(:) !< (ndata)
  integer, intent(in) :: glist(:,:) !< (3, ng)
  integer, intent(in) :: gindex(:) !< (ng)
  integer, intent(in) :: Nfft(:)
  complex(DPC), intent(in) :: fftbox(:,:,:) !< (Nfft(1), Nfft(2), Nfft(3))
  real(DP), intent(in) :: scale
  
  integer :: j, bidx(3)
  
  PUSH_SUB(X(get_from_fftbox))
  
  if(peinf%inode.eq.0) call timacc(97,1)

!$OMP PARALLEL PRIVATE(bidx)
!$OMP DO
  do j=1,ndata
    !data(j) = 0.0
    call gvec_to_fft_index(glist(:,gindex(j)),bidx,Nfft)
    data(j) = fftbox(bidx(1),bidx(2),bidx(3))*scale
  end do
!$OMP END DO
!$OMP END PARALLEL
  
  if(peinf%inode.eq.0) call timacc(97,2)

  POP_SUB(X(get_from_fftbox))
  
  return
end subroutine X(get_from_fftbox)

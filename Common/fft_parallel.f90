!===============================================================================
!
! Module fft_parallel_m
!
! 1. fft_set_p()        Originally By gsm        Last Modified 9/20/2010 (gsm)
!
!    Sets up distribution of planes and rods for parallel 3D (serial 2D/1D) FFT
!
! 2. fft_map_s()        Originally By gsm        Last Modified 9/20/2010 (gsm)
!
!    Builds box to sphere map for serial 3D FFT
!
! 3. fft_map_p()        Originally By gsm        Last Modified 9/20/2010 (gsm)
!
!    Builds box to sphere map for parallel 3D (serial 2D/1D) FFT
!
! 4. fft_r2g_s()        Originally By gsm        Last Modified 9/20/2010 (gsm)
!
!    Performs forward (from R-space to G-space) serial 3D FFT
!
! 5. fft_r2g_p()        Originally By gsm        Last Modified 9/20/2010 (gsm)
!
!    Performs forward (from R-space to G-space) parallel 3D (serial 2D/1D) FFT
!
! 6. fft_g2r_s()        Originally By gsm        Last Modified 9/20/2010 (gsm)
!
!    Performs backward (from G-space to R-space) serial 3D FFT
!
! 7. fft_g2r_p()        Originally By gsm        Last Modified 9/20/2010 (gsm)
!
!    Performs backward (from G-space to R-space) parallel 3D (serial 2D/1D) FFT
!
!===============================================================================
!
!  Forward FFT (from R-space to G-space):
!    f(G) = \sum_r f(r) exp(-iGr) = (N / Omega) \int f(r) exp(-iGr) dr
!
!  Backward FFT (from G-space to R-space):
!    f(r) = \sum_G f(G) exp(iGr)
!
!  Omega = unit cell volume, N = FFTgrid(1)*FFTgrid(2)*FFTgrid(3), FFTgrid = FFT grid size
!
!  Periodic Bloch function:
!    Forward FFT: scale = sqrt(Omega) / N
!      u(G) = (1 / sqrt(Omega)) \int u(r) exp(-iGr) dr
!    Backward FFT: scale = 1 / sqrt(Omega)
!      u(r) = (1 / sqrt(Omega)) \sum_G u(G) exp(iGr)
!  \int |u(r)|^2 dr = \sum_G |u(G)|^2 = 1
!
!  Electron density:
!    Forward FFT: scale = Omega / N
!      n(G) = \int n(r) exp(-iGr) dr
!    Backward FFT: scale = 1 / Omega
!      n(r) = (1 / Omega) \sum_G n(G) exp(iGr)
!  \int n(r) dr = n(G = 0) = number of electrons
!
!  Ionic pseudopotential, Hartree potential, exchange-correlation potential:
!    Forward FFT: scale = 1 / N
!      V(G) = (1 / Omega) \int V(r) exp(-iGr) dr
!    Backward FFT: scale = 1
!      V(r) = \sum_G V(G) exp(iGr)
!  (1 / Omega) \int V(r) dr = V(G = 0) = average potential
!
!  Coulomb potential:
!    Forward FFT: scale = Omega / N
!      V_c(G) = \int V_c(r) exp(-iGr) dr = 4 pi e^2 / G^2
!    Backward FFT: scale = 1 / Omega
!      V_c(r) = (1 / Omega) \sum_G V_c(G) exp(iGr) = e^2 / r
!  e^2 = 2 in Rydberg atomic units
!
!===============================================================================

#include "f_defs.h"

module fft_parallel_m

  use global_m
  use fftw_m

  implicit none

  private

  public :: &
    fft_set_p, &
    fft_map_s, &
    fft_map_p, &
    fft_r2g_s, &
    fft_r2g_p, &
    fft_g2r_s, &
    fft_g2r_p

contains

subroutine fft_set_p(FFTgrid, Nplane, Nrod)

  integer, intent(in) :: FFTgrid(3)
  integer, intent(out) :: Nplane
  integer, intent(out) :: Nrod
  
  PUSH_SUB(fft_set_p)
  
  if (mod(FFTgrid(3),peinf%npes) == 0) then
    Nplane = FFTgrid(3)/peinf%npes
  else
    Nplane = FFTgrid(3)/peinf%npes+1
  endif
  if (mod(FFTgrid(1)*FFTgrid(2),peinf%npes) == 0) then
    Nrod = (FFTgrid(1)*FFTgrid(2))/peinf%npes
  else
    Nrod = (FFTgrid(1)*FFTgrid(2))/peinf%npes+1
  endif
  
  POP_SUB(fft_set_p)
  
  return
  
end subroutine fft_set_p

!-------------------------------------------------------------------------------

subroutine fft_map_s(ng,FFTgrid,isort,gvec,inv_indx)
  
  integer, intent(in) :: ng
  integer, intent(in) :: FFTgrid(3)
  integer, intent(in) :: isort(:) !< (ng)
  integer, intent(in) :: gvec(:,:) !< (3,ngtot)
  integer, intent(out) :: inv_indx(:,:,:) !< (FFTgrid(1),FFTgrid(2),FFTgrid(3))

  integer :: ig,id,idx(3)

  PUSH_SUB(fft_map_s)
  
  inv_indx(:,:,:)=0
  do ig=1,ng
    do id=1,3
      idx(id)=gvec(id,isort(ig))+1
      if (idx(id).lt.1) idx(id)=idx(id)+FFTgrid(id)
    enddo
    inv_indx(idx(1),idx(2),idx(3))=ig
  enddo

  POP_SUB(fft_map_s)
  
  return

end subroutine fft_map_s

!-------------------------------------------------------------------------------

subroutine fft_map_p(imod,ng,ng_l,Nrod,FFTgrid,isort,gvec, &
  mrod,msph,crod,csph,drod,dsph,irod,isph)
  
  integer, intent(in) :: imod
  integer, intent(in) :: ng
  integer, intent(in) :: ng_l
  integer, intent(in) :: Nrod
  integer, intent(in) :: FFTgrid(3)
  integer, intent(in) :: isort(:) !< (ng)
  integer, intent(in) :: gvec(:,:) !< (3,ngtot)
  integer, intent(inout) :: mrod
  integer, intent(inout) :: msph
  integer, intent(inout) :: crod(:) !< (peinf%npes)
  integer, intent(inout) :: csph(:) !< (peinf%npes)
  integer, intent(inout) :: drod(:) !< (peinf%npes)
  integer, intent(inout) :: dsph(:) !< (peinf%npes)
  integer, intent(out) :: irod(:,:) !< (2,MAX(1,mrod))
  integer, intent(out) :: isph(:) !< (MAX(1,msph))

  integer :: i,j,k,ig,id,idx(3)
  integer, allocatable :: orod(:)
  integer, allocatable :: osph(:)
  
  PUSH_SUB(fft_map_p)
  
  if (imod.eq.0) then
    
    crod(:)=0
    csph(:)=0
    do ig=1,ng
      do id=1,3
        idx(id)=gvec(id,isort(ig))+1
        if (idx(id).lt.1) idx(id)=idx(id)+FFTgrid(id)
      enddo
      i=(idx(2)-1)*FFTgrid(1)+idx(1)-1
      j=i/Nrod
      k=(ig-1)/ng_l
      if (peinf%inode.eq.j) crod(k+1)=crod(k+1)+1
      if (peinf%inode.eq.k) csph(j+1)=csph(j+1)+1
    enddo
    do i=1,peinf%npes
      drod(i)=sum(crod(1:i-1))
      dsph(i)=sum(csph(1:i-1))
    enddo
    mrod=sum(crod(:))
    msph=sum(csph(:))
    
  elseif (imod.eq.1) then
    
    irod(:,:)=0
    isph(:)=0
    SAFE_ALLOCATE(orod, (peinf%npes))
    SAFE_ALLOCATE(osph, (peinf%npes))
    orod(:)=drod(:)
    osph(:)=dsph(:)
    do ig=1,ng
      do id=1,3
        idx(id)=gvec(id,isort(ig))+1
        if (idx(id).lt.1) idx(id)=idx(id)+FFTgrid(id)
      enddo
      i=(idx(2)-1)*FFTgrid(1)+idx(1)-1
      j=i/Nrod
      k=(ig-1)/ng_l
      if (peinf%inode.eq.j) then
        orod(k+1)=orod(k+1)+1
        irod(1,orod(k+1))=idx(3)
        irod(2,orod(k+1))=i-j*Nrod+1
      endif
      if (peinf%inode.eq.k) then
        osph(j+1)=osph(j+1)+1
        isph(osph(j+1))=ig-k*ng_l
      endif
    enddo
    SAFE_DEALLOCATE(orod)
    SAFE_DEALLOCATE(osph)
    
  endif
  
  POP_SUB(fft_map_p)
  
  return
  
end subroutine fft_map_p

!-------------------------------------------------------------------------------

subroutine fft_r2g_s(FFTgrid,scale,inv_indx,fftbox_3D,gsphere)
  
  integer, intent(in) :: FFTgrid(3)
  real(DP), intent(in) :: scale
  integer, intent(in) :: inv_indx(:,:,:) !< (FFTgrid(1),FFTgrid(2),FFTgrid(3))
  complex(DPC), intent(inout) :: fftbox_3D(:,:,:) !< (FFTgrid(1),FFTgrid(2),FFTgrid(3))
  complex(DPC), intent(out) :: gsphere(:) !< (ng)
  
  integer :: i1,i2,i3,j1,j2,j3,ig
  
  PUSH_SUB(fft_r2g_s)

  call do_FFT(fftbox_3D,FFTgrid,-1)

  gsphere(:)=(0.0d0,0.0d0)
  do j3=-FFTgrid(3)/2,FFTgrid(3)/2-1
    i3=j3+1
    if (j3.lt.0) i3=FFTgrid(3)+i3
    do j2=-FFTgrid(2)/2,FFTgrid(2)/2-1
      i2=j2+1
      if (j2.lt.0) i2=FFTgrid(2)+i2
      do j1=-FFTgrid(1)/2,FFTgrid(1)/2-1
        i1=j1+1
        if (j1.lt.0) i1=FFTgrid(1)+i1
        ig=inv_indx(i1,i2,i3)
        if (ig.eq.0) cycle
        gsphere(ig)=fftbox_3D(i1,i2,i3)
      enddo
    enddo
  enddo
  
  gsphere(:)=gsphere(:)*scale
  
  POP_SUB(fft_r2g_s)
  
  return
  
end subroutine fft_r2g_s

!-------------------------------------------------------------------------------

subroutine fft_r2g_p(FFTgrid,Nplane,Nrod,mrod,msph,crod,csph, &
  drod,dsph,irod,isph,scale,fftbox_2D,fftbox_1D,buffer_2D,buffer_1D, &
  buffer_rod,buffer_sph,gsphere_d)
  
  integer, intent(in) :: FFTgrid(3)
  integer, intent(in) :: Nplane
  integer, intent(in) :: Nrod
  integer, intent(in) :: mrod
  integer, intent(in) :: msph
  integer, intent(in) :: crod(:) !< (peinf%npes)
  integer, intent(in) :: csph(:) !< (peinf%npes)
  integer, intent(in) :: drod(:) !< (peinf%npes)
  integer, intent(in) :: dsph(:) !< (peinf%npes)
  integer, intent(in) :: irod(:,:) !< (2,MAX(1,mrod))
  integer, intent(in) :: isph(:) !< (MAX(1,msph))
  real(DP), intent(in) :: scale
  complex(DPC), intent(inout) :: fftbox_2D(:,:,:) !< (FFTgrid(1),FFTgrid(2),Nplane)
  complex(DPC), intent(out) :: fftbox_1D(:,:) !< (FFTgrid(3),Nrod)
  complex(DPC), intent(out) :: buffer_2D(:,:,:) !< (Nrod,Nplane,peinf%npes)
  complex(DPC), intent(out) :: buffer_1D(:,:,:) !< (Nrod,Nplane,peinf%npes)
  complex(DPC), intent(out) :: buffer_rod(:) !< (MAX(1,mrod))
  complex(DPC), intent(out) :: buffer_sph(:) !< (MAX(1,msph))
  complex(DPC), intent(out) :: gsphere_d(:) !< (ng_l)

  integer :: i,j,k,i1,i2,i3
  integer :: NSize(3)
  complex(DPC), allocatable :: fftbox_temp(:,:,:)
  
  PUSH_SUB(fft_r2g_p)

  NSize(:)=FFTgrid(:)
  NSize(3)=1
  SAFE_ALLOCATE(fftbox_temp,(NSize(1),NSize(2),NSize(3)))
  do i=1,Nplane
    fftbox_temp(:,:,1)=fftbox_2D(:,:,i)
    call do_FFT(fftbox_temp,NSize,-1)
    fftbox_2D(:,:,i)=fftbox_temp(:,:,1)
  enddo
  SAFE_DEALLOCATE(fftbox_temp) 
  
  buffer_2D(:,:,:)=(0.0d0,0.0d0)
  do k=1,Nplane
    do i2=1,FFTgrid(2)
      do i1=1,FFTgrid(1)
        i=(i2-1)*FFTgrid(1)+i1-1
        j=i/Nrod
        buffer_2D(i-j*Nrod+1,k,j+1)=fftbox_2D(i1,i2,k)
      enddo
    enddo
  enddo
#ifdef MPI
  call MPI_Alltoall(buffer_2D(1,1,1),Nrod*Nplane,MPI_COMPLEX_DPC, &
    buffer_1D(1,1,1),Nrod*Nplane,MPI_COMPLEX_DPC,MPI_COMM_WORLD,mpierr)
#else
  buffer_1D(:,:,:)=buffer_2D(:,:,:)
#endif
  do i3=1,FFTgrid(3)
    do i=1,Nrod
      j=(i3-1)/Nplane
      fftbox_1D(i3,i)=buffer_1D(i,i3-j*Nplane,j+1)
    enddo
  enddo

  NSize(:)=1
  NSize(1)=FFTgrid(3)
  SAFE_ALLOCATE(fftbox_temp,(NSize(1),NSize(2),NSize(3)))
  do i=1,Nrod
    fftbox_temp(:,1,1)=fftbox_1D(:,i)
    call do_FFT(fftbox_temp,NSize,-1)
    fftbox_1D(:,i)=fftbox_temp(:,1,1)
  enddo
  SAFE_DEALLOCATE(fftbox_temp)
  
  buffer_rod(:)=(0.0d0,0.0d0)
  do i=1,mrod
    buffer_rod(i)=fftbox_1D(irod(1,i),irod(2,i))
  enddo
#ifdef MPI
  call MPI_Alltoallv(buffer_rod,crod,drod,MPI_COMPLEX_DPC, &
    buffer_sph,csph,dsph,MPI_COMPLEX_DPC,MPI_COMM_WORLD,mpierr)
#else
  buffer_sph(:)=buffer_rod(:)
#endif
  gsphere_d(:)=(0.0d0,0.0d0)
  do i=1,msph
    gsphere_d(isph(i))=buffer_sph(i)
  enddo
  
  gsphere_d(:)=gsphere_d(:)*scale
  
  POP_SUB(fft_r2g_p)
  
  return
end subroutine fft_r2g_p

!-------------------------------------------------------------------------------

subroutine fft_g2r_s(FFTgrid,scale,inv_indx,gsphere,fftbox_3D)

  integer, intent(in) :: FFTgrid(3)
  real(DP), intent(in) :: scale
  integer, intent(in) :: inv_indx(:,:,:) !< (FFTgrid(1),FFTgrid(2),FFTgrid(3))
  complex(DPC), intent(in) :: gsphere(:) ! <ng)
  complex(DPC), intent(out) :: fftbox_3D(:,:,:) !< (FFTgrid(1),FFTgrid(2),FFTgrid(3))
  
  integer :: i1,i2,i3,j1,j2,j3,ig
  
  PUSH_SUB(fft_g2r_s)
  
  fftbox_3D(:,:,:)=(0.0d0,0.0d0)
  do j3=-FFTgrid(3)/2,FFTgrid(3)/2-1
    i3=j3+1
    if (j3.lt.0) i3=FFTgrid(3)+i3
    do j2=-FFTgrid(2)/2,FFTgrid(2)/2-1
      i2=j2+1
      if (j2.lt.0) i2=FFTgrid(2)+i2
      do j1=-FFTgrid(1)/2,FFTgrid(1)/2-1
        i1=j1+1
        if (j1.lt.0) i1=FFTgrid(1)+i1
        ig=inv_indx(i1,i2,i3)
        if (ig.eq.0) cycle
        fftbox_3D(i1,i2,i3)=gsphere(ig)
      enddo
    enddo
  enddo

  call do_FFT(fftbox_3D,FFTgrid,1)
 
  fftbox_3D(:,:,:)=fftbox_3D(:,:,:)*scale
  
  POP_SUB(fft_g2r_s)
  
  return
end subroutine fft_g2r_s

!-------------------------------------------------------------------------------

subroutine fft_g2r_p(FFTgrid,Nplane,Nrod,mrod,msph,crod,csph, &
  drod,dsph,irod,isph,scale,gsphere_d,fftbox_2D,fftbox_1D,buffer_2D, &
  buffer_1D,buffer_rod,buffer_sph)

  integer, intent(in) :: FFTgrid(3)
  integer, intent(in) :: Nplane
  integer, intent(in) :: Nrod
  integer, intent(in) :: mrod
  integer, intent(in) :: msph
  integer, intent(in) :: crod(:) !< (peinf%npes)
  integer, intent(in) :: csph(:) !< (peinf%npes)
  integer, intent(in) :: drod(:) !< (peinf%npes)
  integer, intent(in) :: dsph(:) !< (peinf%npes)
  integer, intent(in) :: irod(:,:) !< (2,MAX(1,mrod))
  integer, intent(in) :: isph(:) !< (MAX(1,msph))
  real(DP), intent(in) :: scale
  complex(DPC), intent(in) :: gsphere_d(:) !< (ng_l)
  complex(DPC), intent(out) :: fftbox_2D(:,:,:) !< (FFTgrid(1),FFTgrid(2),Nplane)
  complex(DPC), intent(out) :: fftbox_1D(:,:) !< (FFTgrid(3),Nrod)
  complex(DPC), intent(out) :: buffer_2D(:,:,:) !< (Nrod,Nplane,peinf%npes)
  complex(DPC), intent(out) :: buffer_1D(:,:,:) !< (Nrod,Nplane,peinf%npes)
  complex(DPC), intent(out) :: buffer_rod(:) !< (MAX(1,mrod))
  complex(DPC), intent(out) :: buffer_sph(:) !< (MAX(1,msph))

  integer :: i,j,k,i1,i2,i3
  integer :: Nsize(3)
  complex(DPC), allocatable :: fftbox_temp(:,:,:)

  PUSH_SUB(fft_g2r_p)

  buffer_sph(:)=(0.0d0,0.0d0)
  do i=1,msph
    buffer_sph(i)=gsphere_d(isph(i))
  enddo
#ifdef MPI
  call MPI_Alltoallv(buffer_sph,csph,dsph,MPI_COMPLEX_DPC, &
    buffer_rod,crod,drod,MPI_COMPLEX_DPC,MPI_COMM_WORLD,mpierr)
#else
  buffer_rod(:)=buffer_sph(:)
#endif
  fftbox_1D(:,:)=(0.0d0,0.0d0)
  do i=1,mrod
    fftbox_1D(irod(1,i),irod(2,i))=buffer_rod(i)
  enddo

  NSize(:)=1
  NSize(1)=FFTgrid(3)
  SAFE_ALLOCATE(fftbox_temp,(NSize(1),NSize(2),NSize(3)))
  do i=1,Nrod
    fftbox_temp(:,1,1)=fftbox_1D(:,i)
    call do_FFT(fftbox_temp,NSize,1)
    fftbox_1D(:,i)=fftbox_temp(:,1,1)
  enddo
  SAFE_DEALLOCATE(fftbox_temp)
  
  buffer_1D(:,:,:)=(0.0d0,0.0d0)
  do i3=1,FFTgrid(3)
    do i=1,Nrod
      j=(i3-1)/Nplane
      buffer_1D(i,i3-j*Nplane,j+1)=fftbox_1D(i3,i)
    enddo
  enddo
#ifdef MPI
  call MPI_Alltoall(buffer_1D(1,1,1),Nrod*Nplane,MPI_COMPLEX_DPC, &
    buffer_2D(1,1,1),Nrod*Nplane,MPI_COMPLEX_DPC,MPI_COMM_WORLD,mpierr)
#else
  buffer_2D(:,:,:)=buffer_1D(:,:,:)
#endif
  do k=1,Nplane
    do i2=1,FFTgrid(2)
      do i1=1,FFTgrid(1)
        i=(i2-1)*FFTgrid(1)+i1-1
        j=i/Nrod
        fftbox_2D(i1,i2,k)=buffer_2D(i-j*Nrod+1,k,j+1)
      enddo
    enddo
  enddo

  NSize(:)=FFTgrid(:)
  NSize(3)=1
  SAFE_ALLOCATE(fftbox_temp,(NSize(1),NSize(2),NSize(3)))
  do i=1,Nplane
    fftbox_temp(:,:,1)=fftbox_2D(:,:,i)
    call do_FFT(fftbox_temp,NSize,1)
    fftbox_2D(:,:,i)=fftbox_temp(:,:,1)
  enddo
  SAFE_DEALLOCATE(fftbox_temp)
  
  fftbox_2D(:,:,:)=fftbox_2D(:,:,:)*scale
  
  POP_SUB(fft_g2r_p)
  
  return
end subroutine fft_g2r_p

end module fft_parallel_m

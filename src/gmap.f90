#include "f_defs.h"

module gmap_m

  use global_m
  use misc_m

  implicit none

  private

  public :: gmap_2 !, gmap

  ! interface gmap
  !    module procedure dgmap, zgmap
  ! end interface gmap

  interface gmap_2
     module procedure dgmap_2, zgmap_2
  end interface gmap_2
  
contains

  ! subroutine dgmap(gvec, syms, ngk, itran, kgq, isortc, isorti, ind, phase, die_outside_sphere)
  !   ! subroutine dgmap(gvec, syms, ngk, itran, kgq, fk, isortc, isorti, ind, phase, die_outside_sphere)
  !   type (gspace), intent(in) :: gvec         !< uses gvec%components(:,isrtc(1:ngk)), and index_vec
  !   type (symmetry), intent(in) :: syms       !< uses syms%mtrx(:,:,itran) & syms%tnp(:,itran)
  !   integer, intent(in) :: ngk                !< number of g-vector entries in a wavefunction
  !   integer, intent(in) :: itran              !< index of transformation
  !   integer, intent(in) :: kgq(3)             !< an umklapp vector (i.e. integer 3-vector)
  !   ! real(DP), intent(in) :: fk(3)
  !   integer, intent(in) :: isortc(:)          !< index array for R(q) (1:ngk)
  !   integer, intent(in) :: isorti(:)          !< inverse index array for q (1:gvec%ng)
  !   integer, intent(out) :: ind(:)            !< indices for the vectors inv(symm(itran))*(g+kgq) (ngk)
  !   real(DP), intent(out) :: phase(:)         !< exp(-i*(g+kgq).dot.syms%tnp(itran)) (ngk)
  !   logical, intent(in) :: die_outside_sphere !< specifies whether to die if G-vectors are falling outside of the sphere

  !   PUSH_SUB(dgmap)

  !   call gmap_base( gvec, syms, ngk, itran, kgq, isortc, isorti, ind, die_outside_sphere, dphase = phase)

  !   POP_SUB(dgmap)
  !   return
  ! end subroutine dgmap

  subroutine dgmap_2(gvec, syms, ngk, itran, kgq, isortc, isorti, ind, phase)
    type (gspace), intent(in) :: gvec         !< uses gvec%components(:,isrtc(1:ngk)), and index_vec
    type (symmetry), intent(in) :: syms       !< uses syms%mtrx(:,:,itran) & syms%tnp(:,itran)
    integer, intent(in) :: ngk                !< number of g-vector entries in a wavefunction
    integer, intent(in) :: itran              !< index of transformation
    integer, intent(in) :: kgq(3)             !< an umklapp vector (i.e. integer 3-vector)
    integer, intent(in) :: isortc(:)          !< index array for R(q) (1:ngk)
    integer, intent(in) :: isorti(:)          !< inverse index array for q (1:gvec%ng)
    integer, intent(out) :: ind(:)            !< indices for the vectors inv(symm(itran))*(g+kgq) (ngk)
    real(DP), intent(out) :: phase(:)         !< exp(-i*(g+kgq).dot.syms%tnp(itran)) (ngk)

    PUSH_SUB(dgmap_2)

    call gmap_base_2( gvec, syms, ngk, itran, kgq, isortc, isorti, ind, dphase = phase)

    POP_SUB(dgmap_2)
    return
  end subroutine dgmap_2
  
  ! subroutine zgmap(gvec, syms, ngk, itran, kgq, isortc, isorti, ind, phase, die_outside_sphere)

  !   type (gspace), intent(in) :: gvec         !< uses gvec%components(:,isrtc(1:ngk)), and index_vec
  !   type (symmetry), intent(in) :: syms       !< uses syms%mtrx(:,:,itran) & syms%tnp(:,itran)
  !   integer, intent(in) :: ngk                !< number of g-vector entries in a wavefunction
  !   integer, intent(in) :: itran              !< index of transformation
  !   integer, intent(in) :: kgq(3)             !< an umklapp vector (i.e. integer 3-vector)
  !   integer, intent(in) :: isortc(:)          !< index array for R(q) (1:ngk)
  !   integer, intent(in) :: isorti(:)          !< inverse index array for q (1:gvec%ng)
  !   integer, intent(out) :: ind(:)            !< indices for the vectors inv(symm(itran))*(g+kgq) (ngk)
  !   complex(DPC), intent(out) :: phase(:)     !< exp(-i*(g+kgq).dot.syms%tnp(itran)) (ngk)
  !   logical, intent(in) :: die_outside_sphere !< specifies whether to die if G-vectors are falling outside of the sphere

  !   PUSH_SUB(zgmap)

  !   call gmap_base(gvec, syms, ngk, itran, kgq, isortc, isorti, ind, die_outside_sphere, zphase = phase)

  !   POP_SUB(zgmap)
  !   return
  ! end subroutine zgmap

  subroutine zgmap_2(gvec, syms, ngk, itran, kgq, isortc, isorti, ind, phase)
    type (gspace), intent(in) :: gvec         !< uses gvec%components(:,isrtc(1:ngk)), and index_vec
    type (symmetry), intent(in) :: syms       !< uses syms%mtrx(:,:,itran) & syms%tnp(:,itran)
    integer, intent(in) :: ngk                !< number of g-vector entries in a wavefunction
    integer, intent(in) :: itran              !< index of transformation
    integer, intent(in) :: kgq(3)             !< an umklapp vector (i.e. integer 3-vector)
    integer, intent(in) :: isortc(:)          !< index array for R(q) (1:ngk)
    integer, intent(in) :: isorti(:)          !< inverse index array for q (1:gvec%ng)
    integer, intent(out) :: ind(:)            !< indices for the vectors inv(symm(itran))*(g+kgq) (ngk)
    complex(DPC), intent(out) :: phase(:)     !< exp(-i*(g+kgq).dot.syms%tnp(itran)) (ngk)

    PUSH_SUB(zgmap_2)

    call gmap_base_2(gvec, syms, ngk, itran, kgq, isortc, isorti, ind, zphase = phase)

    POP_SUB(zgmap_2)
    return
  end subroutine zgmap_2
  
  !==================================================================

!   subroutine gmap_base(gvec, syms, ngk, itran, kgq, isortc, isorti, ind, die_outside_sphere, dphase, zphase)
!     type (gspace), intent(in) :: gvec         !< uses gvec%components(:,isrtc(1:ngk)), and index_vec
!     type (symmetry), intent(in) :: syms       !< uses syms%mtrx(:,:,itran) & syms%tnp(:,itran)
!     integer, intent(in) :: ngk                !< number of g-vector entries in a wavefunction
!     integer, intent(in) :: itran              !< index of transformation
!     integer, intent(in) :: kgq(3)             !< an umklapp vector (i.e. integer 3-vector)
!     integer, intent(in) :: isortc(:)          !< index array for R(q) (1:ngk)
!     integer, intent(in) :: isorti(:)          !< inverse index array for q (1:gvec%ng)
!     integer, intent(out) :: ind(:)            !< indices for the vectors inv(symm(itran))*(g+kgq) (ngk)
!     logical, intent(in) :: die_outside_sphere !< specifies whether to die if G-vectors are falling outside of the sphere
!     real(DP),     optional, intent(out) :: dphase(:) !< exp(-i*(g+kgq).dot.syms%tnp(itran)) (ngk)
!     complex(DPC), optional, intent(out) :: zphase(:) !< exp(-i*(g+kgq).dot.syms%tnp(itran)) (ngk)
!     integer :: ig, kd(3), kadd, kgrad, kgrad1
!     integer :: kg(3), kgr(3), nout, nin
!     integer, dimension(3,3) :: inv_mtrx_reci
!     real(DP) :: fi

!     PUSH_SUB(gmap_base)

!     if(present(dphase) .and. present(zphase)) then
!        call die("gmap: cannot pass both dphase and zphase")
!     else if(.not. present(dphase) .and. .not. present(zphase)) then
!        call die("gmap: must pass either dphase or zphase")
!     endif

!     if(ngk > gvec%ng) call die("gmap: ngk (wfn cutoff) is greater than gvec%ng (rho cutoff)")
!     if(ubound(isorti, 1) < gvec%ng)  call die("gmap: isorti size < gvec%ng")
!     if(any(isorti(1:gvec%ng) > gvec%ng)) call die("gmap: isorti cannot be greater than gvec%ng.")

!     if(ubound(isortc, 1) < ngk)      call die("gmap: isortc size < ngk")
!     if(any(isortc(1:ngk) < 1))       call die("gmap: isortc cannot be less than 1.")
!     if(any(isortc(1:ngk) > gvec%ng)) call die("gmap: isortc cannot be greater than ng.")

!     if(ubound(gvec%index_vec, 1) /= gvec%nFFTgridpts) call die("gmap: gvec%index_vec has wrong size")
!     if(any(gvec%index_vec(1:gvec%nFFTgridpts) < 0)) call die("gmap: index_vec cannot be less than 0")
!     if(any(gvec%index_vec(1:gvec%nFFTgridpts) > gvec%ng)) call die("gmap: index_vec cannot be greater than ng")

!     if(present(dphase)) then
!        if(ubound(dphase, 1) < ngk) call die("gmap: dphase size < ngk")
!     else
!        if(ubound(zphase, 1) < ngk) call die("gmap: zphase size < ngk")
!     endif
!     if(ubound(ind, 1) < ngk) call die("gmap: ind size < ngk")

!     ! Invert rotation matrix that gives rq
!     ! -------------------------------------------------
!     ! itran(if) = itran
!     ! kg%kg0(1:3,if) + syms%mtrx(1:3, 1:3, itran(if))*kg%r(:,indr(if))
!     ! = kg%f(1:3,if)

!     ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!     ! call invert_matrix_int(syms%mtrx(1:3, 1:3, itran), mtrxi(1:3, 1:3))
!     ! first transform mtrx_reci into cartesian coordinates
!     ! and then take inverse: mtrx_reci ==> mtrx_cart ==> inv_mtrx_cart
!     ! and then transform it back into mtrx_reci: inv_mtrx_cart ==> inv_mtrx_reci

!     ! !write(6,*) peinf%inode,'mtrxi: '
!     ! do i = 1, 3
!     !    write(6,*) syms%mtrx_reci(i,:,itran)
!     ! enddo
!     ! ! write(6,*) peinf%inode,'kgq',kgq

!     call invert_matrix_int(syms%mtrx_reci(1:3,1:3,itran), inv_mtrx_reci(1:3,1:3))
!     ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<
!     ! Construct u(alpha) here using syms%mtrx(1:3, 1:3, itran)
!     ! so3_2_su2(crys%avec(1:3,1:3), syms%mtrx(1:3, 1:3, itran),umtrx)
!     ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>
!     ! JRD: Temporary Debugging

!     ! Loop over g-vectors in function of r(q)

!     nout = 0  ! number of waves outside sphere
!     nin = 0  ! number of waves inside sphere

!     ! ngk = wfn%ng = work%ng
!     !     = intwfn%ng(ii) = ngk_g(kg%indr(peinf%ik(peinf%inode+1,ii)))
!     ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!     ! Should use ig_ifk

!     do ig = 1, ngk ! = pol%nmtx (Epsilon)
!        ! kgq = kg%kg0(:,ik2), ik2 = ik = peinf%ik(peinf%inode+1,ikt)
!        ! kg = g(ig) + kgq
!        ! ------------------------------------------
!        ! in genwf.f90:
!        ! call kinetic_energies(gvec, crys%bdot, ekin, qvec = kg%f(:, ik2))
!        ! call sortrx(gvec%ng, ekin, work%isort, gvec = gvec%components)
!        !
! !!!!!! Note that isort is prepared after the rearrangement of gvec%components
!        !
!        ! isortc = work%isort
!        ! ------------------------------------------
!        !< gvec%components(1:3,isortc(ig)) is the Gvector in global G list with
!        !< increasing magnitude of |G(1:3,:)+kg%f(1:3,peinf%ik(peinf%inode+1, ikt))|^2
!        !< here G' = gvec%components(1:3, isortc(ig)) is the G' in C_{G' n [\alpha k]}
!        !< That is the ig-th G' for intwfn with rotated k (irk = indr(ifk))
!        !< And now we want to find out the corresponding index of \alpha^{-1} (G'+G^{0}_{\alpha k}) on the Gvector list for RBZ kpoint (original intwfn)
!        ! ------------------------------------------------
!        !< kg(1:3) = G'+G^{0}_{\alpha k}
!        kg(1:3) = gvec%components(1:3, isortc(ig)) + kgq(1:3)
!        ! T is the itran-th symmetry operation
!        ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!        ! kgr(1:3) = \alpha^{-1}[G'+G^{0}_{\alpha k}]
!        ! kgr(1:3) = MATMUL(mtrxi(1:3, 1:3), kg(1:3))
!        kgr(1:3) = MATMUL(dble(inv_mtrx_reci(1:3, 1:3)), kg(1:3))
!        ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!        ! -----------------------------------
!        !< The following shows, given a Gvector kgr(1:3), how can we
!        !< find its unique index 'kadd', and use this address and the map
!        !< gvec%index_vec(1:) to calculate the index of this kgr(1:3) in
!        !< global Gvector list gvec%components(1:3,1:ngm_g)
!        ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!        ! ***************************************
!        !< We can capsulate this functionality into a function
!        !< call Get_Gindex(kgr(1:3), gvec, Gindex)
!        ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!        ! -----------------------------------
!        kd(1:3) = kgr(1:3) + gvec%FFTgrid(1:3) / 2 + 1

!        if (any(kd(1:3) .lt. 1 .or. kd(1:3) .gt. gvec%FFTgrid(1:3))) then
!           call die('gmap: kd out of bounds')
!        endif

!        !< in Common/input_utils.f90:
!        !< address=((g(1)+FFTgrid(1)/2)*FFTgrid(2)+g(2)+FFTgrid(2)/2)*FFTgrid(3)+g(3)+FFTgrid(3)/2+1
!        !< kadd gives the index on FFTgrid of the Gvector kgr(1:3)
!        !< kadd <==> kgr(1:3) the same thing
!        kadd = ((kd(1) - 1) * gvec%FFTgrid(2) + kd(2) - 1) * gvec%FFTgrid(3) + kd(3)
!        ! in BSE/input.f90:
!        ! call gvec_index(gvec)
!        !
!        ! kgrid1 gives the index of kgr(1:3) in global g_g(1:3,1:gvec%ng = ngm_g) list
!        ! -----------------------------------------------------------
!        ! kgrad1 <==> kadd <==> kgr(1:3) all denote the same Gvector
!        kgrad1 = gvec%index_vec(kadd)

!        if (kgrad1 .lt. 1 .or. kgrad1 .gt. gvec%ng) then
!           write(0,*) 'itran = ', itran, 'ig = ', ig, ', kadd = ', kadd, ', kgrad1 = ', kgrad1
!           call die('gmap: G-vectors falling outside of the charge-density G-sphere')
!        endif

!        !< isorti: G index on global g_g(:,:) list |-----> G index on local CPU
!        ! ----------------------------------------------------
!        !< kgrad gives Gindex of \alpha^{-1}(G'+G^{0}_{\alpha k}) on the Gvector list for RBZ kpoint (original intwfn)
!        !< for the indr(peinf%ik(peinf%inode+1,ikt))-th RBZ kpoint
!        ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!        !< If kgrid1 is really out of sphere around kg%r(1:3,irk=indr(ifk)),
!        !< then isorti will be 0
!        ! -------------------------------------------------------
!        !< in genwf.f90: use old isort(:) to generate isorti
!        !< here work%isort gives the global Gindex for irk=indr(ifk)-th RBZ kpt
!        !<  SAFE_ALLOCATE(isorti, (gvec%ng))
!        !<  isorti(:)=0
!        !< do ii=1,wfn%ng ! = ngk_g(irk)
!        !<    isorti(work%isort(ii))=ii
!        !< enddo
!        ! ******************************************
!        ! Should use kgrad ==> ig_irk
!        kgrad = isorti(kgrad1)

!        ! write(*,"(A,I7,A,I7,A,I7,A,I7,A,I7)") "ngk = ", ngk, "kadd = ", kadd, "kgrad1 = ", kgrad1, "ig = ", ig, "kgrad = ", kgrad
!        ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!        ! --------------------------------------------------
!        !< SIB: if kgr is outside the sphere, then increment out counter,
!        !< set its phase to zero, and its ind() entry to the maximum.
!        ! -------------------------------------------
!        !< ifk = peinf%ik(peinf%inode+1, ikt)
!        !< irk = indr(ifk)
!        ! -------------------------------------------
!        !< if the G list for ikt-th kpoint on current CPU does not contain
!        !< kgrad
!        ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!        ! in genwf.f90:
!        ! isorti(:)=0
!        ! do ii=1,wfn%ng
!        !    isorti(work%isort(ii))=ii
!        ! enddo
!        ! -------------------------------------------
!        ! outside sphere Gvectors center on irk-th RBZ kpoint
!        ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!        ! *******************************************
!        !< According to how isorti is constructed in genwf.f90
!        !< kgrad = 0 : ngk, will never > ngk
!        !> In fact, we do not need this IF,
!        !> Since kgrad is a gvector from old wavefunction, and should not be bound by the number of Gvectors in the new wavefunction
!        if (kgrad .gt. ngk) then
!           nout = nout + 1
!           ! Set ind(ig) = MAX index of Glist for irk-th RBZ kpoint
!           ind(ig) = ngk
!           if(present(zphase)) then
!              zphase(ig) = DCMPLX(0d0, 0d0)
!           else
!              dphase(ig) = 0d0
!           endif
!        else
!           ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!           ! SIB:  Otherwise, record the index of kgr (kgrad) into ind(ig)
!           !< and compute the phase = exp(-i*kg.dot.syms%tnp(:,itran))
!           ! ------------ IMPORTANT -------------
!           !< This means, the ig-th Gvector for the ifk-th FBZ kpt corresponds
!           !< to the kgrad-th Gvector for the irk=indr(ifk)-th RBZ kpt
!           ind(ig) = kgrad
!           ! --------------------------------------------
!           ! kg(1:3) = gvec%components(1:3, isortc(ig)) + kgq(1:3)
!           ! itran(if) = itran
!           ! kg%kg0(1:3,if) + syms%mtrx(1:3, 1:3, itran(if))*kg%r(:,indr(if))
!           ! = kg%f(1:3,if)
!           ! -----------------------------------------------
!           ! <<<<<<<<<<<<<<<<<<<<<<<<<<<
!           !< syms%tnp = ft_minus(1:3,1:ntran) from pw2bgw.f90
!           !< fi = ([\alpha k]+G) \cdot t
!           ! ----------------------------
!           ! ======
!           ! [Original code] :
!           ! fi = dot_product(dble(kg(:)), syms%tnp(:,itran))
!           ! here syms%tnp use * 2pi units
!           ! ======
!           !< B in units 2 \Pi/alat
!           !< A^{T} in units alat
!           ! fi = 2 * PI_D * DOT_PRODUCT(fk(1:3)+gvec%components(1:3, isortc(ig)),syms%tnp(1:3,itran))
!           fi = 2.0D0 * PI_D * DOT_PRODUCT(dble(gvec%components(1:3, isortc(ig))),syms%tnp(1:3,itran))
!           ! >>>>>>>>>>>>>>>>>>>>>>>>>>
!           if(present(zphase)) then
!              ! Exp[-I*fi]
!              zphase(ig) = DCMPLX(cos(fi), -sin(fi))
!           else
!              ! ------------ IMPORTANT -----------
!              !< DAS: The imaginary part can be thrown away because it is always zero
!              !< if we have inversion and time-reversal symmetries, and so can use the real version.
!              !< phase = +/- 1. Otherwise the wavefunction would not be normalized.
!              !< c(G) -> c(G) e^iGt with fractional translation.
!              !< c(G) e^iGt = c(-G)* e^-iGt by time-reversal symmetry
!              !<  = c(G)* e^-iGt by inversion symmetry. c(G) = c(G)* since real.
!              !< Therefore e^iGt = e^-iGt. e^iGt is real, and hence 1 or -1.
!              !< Note there is also a global phase e^ikt, but it is just a convention
!              !< and can be safely ignored here.
!              ! ----------------------------------
!              dphase(ig) = cos(fi)

!              if(abs(abs(dphase(ig)) - 1) .gt. TOL_Small) then
!                 write(0,'(a,i8,a,f12.8,a)') 'phase(', ig, ') = ', dphase(ig), ' != +/- 1'
!                 call die("Illegal non-unity phase in gmap, error in fractional translation.")
!              endif
!              if(abs(sin(fi)) .gt. TOL_Small) then
!                 write(0,'(a,i8,a,f12.8)') 'Im phase(', ig, ') = ', -sin(fi)
!                 call die("Illegal complex phase in gmap, error in fractional translation.")
!              endif
!           endif
!        endif

!     enddo !end loop over g-vectors (ig)

!     if (die_outside_sphere .and. nout .gt. 0) then
!        call die('G-vectors are falling outside of the sphere in gmap')
!     endif

!     !  if (die_outside_sphere .and. nin .gt. 0) then
!     !    call die('G-vectors are falling inside of the sphere in gmap')
!     !  endif

!     POP_SUB(gmap_base)

!     return
!   end subroutine gmap_base

  !> This gmap_base_2 allow ig_old exceed ngk, which is the largest index of ig_new
  subroutine gmap_base_2(gvec, syms, ngk, itran, kgq, isortc, isorti, ind, dphase, zphase)
    type (gspace), intent(in) :: gvec         !< uses gvec%components(:,isrtc(1:ngk)), and index_vec
    type (symmetry), intent(in) :: syms       !< uses syms%mtrx(:,:,itran) & syms%tnp(:,itran)
    integer, intent(in) :: ngk                !< number of g-vector entries in a wavefunction
    integer, intent(in) :: itran              !< index of transformation
    integer, intent(in) :: kgq(3)             !< an umklapp vector (i.e. integer 3-vector)
    integer, intent(in) :: isortc(:)          !< index array for R(q) (1:ngk)
    integer, intent(in) :: isorti(:)          !< inverse index array for q (1:gvec%ng)
    integer, intent(out) :: ind(:)            !< indices for the vectors inv(symm(itran))*(g+kgq) (ngk)
    real(DP),     optional, intent(out) :: dphase(:) !< exp(-i*(g+kgq).dot.syms%tnp(itran)) (ngk)
    complex(DPC), optional, intent(out) :: zphase(:) !< exp(-i*(g+kgq).dot.syms%tnp(itran)) (ngk)
    integer :: ig, kd(3), kadd, kgrad, kgrad1
    integer :: kg(3), kgr(3)
    integer, dimension(3,3) :: inv_mtrx_reci_int
    real(DP), dimension(3,3) :: inv_mtrx_reci, mtrx_reci, diff_inv_mtrx_reci
    logical :: flag_inv
    real(DP) :: fi
    PUSH_SUB(gmap_base_2)

    if(present(dphase) .and. present(zphase)) then
       call die("gmap: cannot pass both dphase and zphase")
    else if(.not. present(dphase) .and. .not. present(zphase)) then
       call die("gmap: must pass either dphase or zphase")
    endif

    if (ngk > gvec%ng) call die("gmap: ngk (wfn cutoff) is greater than gvec%ng (rho cutoff)")
    if (ubound(isorti, 1) < gvec%ng) call die("gmap: isorti size < gvec%ng")
    if (any(isorti(1:gvec%ng) > gvec%ng)) call die("gmap: isorti cannot be greater than gvec%ng.")
    if (ubound(isortc, 1) < ngk)      call die("gmap: isortc size < ngk")
    if (any(isortc(1:ngk) < 1))       call die("gmap: isortc cannot be less than 1.")
    if (any(isortc(1:ngk) > gvec%ng)) call die("gmap: isortc cannot be greater than ng.")
    if (ubound(gvec%index_vec, 1) /= gvec%nFFTgridpts) call die("gmap: gvec%index_vec has wrong size")
    if (any(gvec%index_vec(1:gvec%nFFTgridpts) < 0)) call die("gmap: index_vec cannot be less than 0")
    if (any(gvec%index_vec(1:gvec%nFFTgridpts) > gvec%ng)) call die("gmap: index_vec cannot be greater than ng")
    if(present(dphase)) then
       if(ubound(dphase, 1) < ngk) call die("gmap: dphase size < ngk")
    else
       if(ubound(zphase, 1) < ngk) call die("gmap: zphase size < ngk")
    endif
    if(ubound(ind, 1) < ngk) call die("gmap: ind size < ngk")

    mtrx_reci(1:3,1:3) = dble(syms%mtrx_reci(1:3,1:3,itran))
    call M33INV(mtrx_reci,inv_mtrx_reci,flag_inv)
    if (.not. flag_inv) then
       call die("gmap_2: mtrx_reci not invertible.", only_root_writes=.true.)
    endif

    inv_mtrx_reci_int(:,:) = NINT(inv_mtrx_reci(:,:))
    diff_inv_mtrx_reci(:,:) = DBLE(inv_mtrx_reci_int(:,:)) - inv_mtrx_reci(:,:)
    IF ( NORM2(diff_inv_mtrx_reci) > TOL_Zero ) THEN
       call die("gmap_2: inv_mtrx_reci not integer ", only_root_writes=.true.)
    ENDIF   
    
    !> Loop over g-vectors in new wave function
    do ig = 1, ngk
       kg(1:3) = gvec%components(1:3, isortc(ig)) + kgq(1:3)
       kgr(1:3) = MATMUL(inv_mtrx_reci_int(1:3, 1:3), kg(1:3))      
       kd(1:3) = kgr(1:3) + gvec%FFTgrid(1:3) / 2 + 1
       if (any(kd(1:3) .lt. 1 .or. kd(1:3) .gt. gvec%FFTgrid(1:3))) then
          call die('gmap: kd out of bounds')
       endif

       kadd = ((kd(1) - 1) * gvec%FFTgrid(2) + kd(2) - 1) * gvec%FFTgrid(3) + kd(3)
       !> Common/input_utils.f90:
       !> gvec%index_vec(ig_FFT) = ig_gvec \in [1, gvec%ng]
       kgrad1 = gvec%index_vec(kadd)
       if (kgrad1 .lt. 1 .or. kgrad1 .gt. gvec%ng) then
          write(0,*) 'itran = ', itran, 'ig = ', ig, ', kadd = ', kadd, ', kgrad1 = ', kgrad1
          call die('gmap: G-vectors falling outside of the charge-density G-sphere')
       endif
       !> isorti(ig_gvec) = ig_old
       kgrad = isorti(kgrad1)
       !> kgrad can be 0 if kgrad1 is outside the cutoff of wfn_old!
       ind(ig) = kgrad
       fi = 2.0D0 * PI_D * DOT_PRODUCT(DBLE(gvec%components(1:3, isortc(ig))),syms%tnp(1:3,itran))
       if(present(zphase)) then
          zphase(ig) = DCMPLX(cos(fi), -sin(fi))
       else
          dphase(ig) = cos(fi)
          if(abs(abs(dphase(ig)) - 1.0D0) .gt. TOL_Small) then
             write(0,'(a,i8,a,f12.8,a)') 'phase(', ig, ') = ', dphase(ig), ' != +/- 1'
             call die("Illegal non-unity phase in gmap, error in fractional translation.")
          endif

          if(abs(sin(fi)) .gt. TOL_Small) then
             write(0,'(a,i8,a,f12.8)') 'Im phase(', ig, ') = ', -sin(fi)
             call die("Illegal complex phase in gmap, error in fractional translation.")
          endif         
       endif
    enddo !end loop over g-vectors (ig)

    POP_SUB(gmap_base_2)
    return
  end subroutine gmap_base_2
end module gmap_m

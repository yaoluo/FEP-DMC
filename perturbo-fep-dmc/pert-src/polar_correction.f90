!===============================================================================
! Copyright (C) 2016-2020 Jin-Jian Zhou, Jinsoo Park, I-Te Lu, Marco Bernardi
!
! This program is distributed under the terms of the GNU General Public License.
! See the file `LICENSE' in the root directory of this distribution, or obtain 
! a copy of the License at <https://www.gnu.org/licenses/gpl-3.0.txt>.
!
! Author: jjzhou <jjchou.comphy@gmail.com>
! Comment:
!   subroutines for interpolation of long-range electron-phonon matrix
!   element in polar semiconductor. substrate and reintroduce on el-optimal
!   smooth Bloch states (or wannier gauge).
!
!   Note: the approximation <psi_m,k+q| e^{i(q+G)r} |psi_n,k> ~ delta_nm
!   holds under the wannier gauge, aka the rotated bloch wave-function,
!   which gives the localized wannier functions after Fourier transform.
!  
!   Calculate the long-range el-ph matrix element for given q in phonon 
!   Cartesian coordinate (g_tau,s(q), tau->atomic index, s->x,y,z). 
!   See eq.(12) in PRB 92, 054307. the eph matrix in eigen-mode coordinate 
!   g_mu(q) is simply dot product of eigen-displacement u_mu(q) and g_tau,s(q).
!   see eq.(2) in PRB 92, 054307. also see eq.(4) in PRL 115, 176401.
!   (note, 4pi*epsilon_0 = 1 in ryd atom unit, Volume Omega can't be ignored)
!
!   refers: J.Sjakste PRB 92, 054307 & C.Verdi PRL 115, 176401
!     2D polar correction: T.Sohier PRB 94, 085415.
!  
! Maintenance:
!   jpark: polar correction for 2D system.
!===============================================================================
module polar_correction
   use kinds, only: dp
   implicit none
   private
   !
   complex(dp), parameter :: czero=(0.E0_dp, 0.E0_dp), cone=(1.E0_dp, 0.E0_dp)
   complex(dp), parameter :: conei=(0.E0_dp, 1.E0_dp)
   !
   real(dp),    parameter :: pi = 3.141592653589793238462643383279_dp
   real(dp),    parameter :: tpi = 2.0_dp*pi
   real(dp),    parameter :: fpi = 4.0_dp*pi
   !
   real(dp),    parameter :: e2 = 2.0_dp  !e^2 in Rydberg atomic unit
   
   !pack all the info needed in polar correction
   type, public :: polar_parameter
      !alpha is the Ewald summation parameter (in unit of (tpiba)**2 )
      real(dp) :: alpha
      real(dp) :: gmax
      integer  :: nrx(3)
      ! put 'nat' right after 'nrx(3)' to align memory layout
      ! number of atoms
      integer :: nat  
      integer :: nrx_ph(3)
      ! nrx for phonons (identical to rgd_blk in rigid.f90)
      !
      ! unit cell volume
      real(dp) :: omega, tpiba
      ! dielectric tensor, bg(:,i): i-th recip. lattice vector.
      real(dp) :: epsil(3,3), bg(3,3)
      ! born effective charge
      real(dp), pointer :: bcharge(:,:,:) => null()
      ! atomic position in cartesian coordinate (in the unit of alat)
      real(dp), pointer :: tau(:,:) => null()
      !
      !thickness of the 2d system (unit of alat)
      ! if thickness_2d_ba > 0, then it's 2D system, we do 2d_polar
      real(dp) :: thickness_2d_ba
      !
      !on-site correction to enforce acoustic sum rule for phonon
      real(dp), pointer :: onsite_correction(:,:,:) => null()
   end type

   public set_polar_parameter
   public eph_wan_longrange
   public dyn_mat_longrange
contains

subroutine set_polar_parameter &
      (pdim, nat, omega, tpiba, bg, epsil, bcharge, tau, pdata, thickness_in, alpha_in)
   !use cell_base, only: omega, tpiba, bg
   !use ions_base, only: nat, tau
   implicit none
   ! number of phonon modes, 3*nat
   integer,  intent(in) :: nat, pdim(3) !pdim: specify system dimemsion
   ! 2\pi/alat, unit cell volume, and recip. lattice vector. (see cell_base)
   real(dp), intent(in) :: tpiba, omega, bg(3,3)
   !dielectric constant tensor
   real(dp), intent(in) :: epsil(3,3) 
   ! Born effective charges tensor
   real(dp), intent(in) :: bcharge(3,3,nat) 
   ! atomic positions in cartesian coordiante (unit of alat)
   real(dp), intent(in) :: tau(3, nat)
   ! 
   real(dp), intent(in), optional :: alpha_in
   ! thickness of the 2d system (unit of bohr)
   ! negative value means it is a 3D system, not a 2D material.
   real(dp), intent(in), optional :: thickness_in
   !
   type(polar_parameter), intent(inout) :: pdata

   !local
   integer :: nrx1, nrx2, nrx3, nrx1_ph, nrx2_ph, nrx3_ph
   real(dp) :: alpha, ggmax

   pdata%nat   = nat
   pdata%omega = omega
   pdata%tpiba = tpiba
   pdata%epsil = epsil
   pdata%bg    = bg
   
   if( associated(pdata%bcharge) ) deallocate( pdata%bcharge )
   if( associated(pdata%tau) ) deallocate( pdata%tau )
   allocate( pdata%bcharge(3, 3, nat), pdata%tau(3, nat) )
   pdata%bcharge(1:3, 1:3, 1:nat) = bcharge(1:3, 1:3, 1:nat)
   pdata%tau = tau
   
   alpha = 1.0E0_dp  ! default value
   if( present(alpha_in) ) then
      if( alpha_in < 1.0E-12 ) &
         call errore('set_polar_parameter','alpha too small !!!',1)
      alpha = alpha_in
   endif
   pdata%alpha = alpha
   !
   pdata%thickness_2d_ba = -1.0_dp !default, negative value means 3D system (no 2d polar)
   if( present(thickness_in) ) then
      pdata%thickness_2d_ba = thickness_in * tpiba / tpi
   endif
   !  
   ! exp(-14) ~ 10^-6
   pdata%gmax = 14.0E0_dp
   ggmax = pdata%gmax * 4.0E0_dp * alpha
   !
   ! Estimate of nrx1,nrx2,nrx3 generating all vectors up to G*epsi*G < ggmax
   ! bg is in the unit of tpiba, need to check !
   ! ggmax is the maximum value of (q+G)*epsi*(q+G)
   nrx1 = ceiling(sqrt(ggmax/dot_product(bg(:,1), matmul(epsil, bg(:,1)))))
   nrx2 = ceiling(sqrt(ggmax/dot_product(bg(:,2), matmul(epsil, bg(:,2)))))
   nrx3 = ceiling(sqrt(ggmax/dot_product(bg(:,3), matmul(epsil, bg(:,3)))))
   if(pdim(1) < 2) nrx1 = 0
   if(pdim(2) < 2) nrx2 = 0
   if(pdim(3) < 2) nrx3 = 0
   !
   pdata%nrx(1:3) = (/nrx1, nrx2, nrx3/)
   !
   ! for phonons : follow rgd_blk in rigid.f90
   ! Estimate of nrx1_ph,nrx2_ph,nrx3_ph generating all vectors up to G^2 < geg
   ! bg is in the unit of tpiba, need to check !
   ! ggmax is the maximum value of (q+G)*epsi*(q+G)
   !
   nrx1_ph = ceiling(sqrt(ggmax/dot_product(bg(:,1), bg(:,1))))
   nrx2_ph = ceiling(sqrt(ggmax/dot_product(bg(:,2), bg(:,2))))
   nrx3_ph = ceiling(sqrt(ggmax/dot_product(bg(:,3), bg(:,3))))
   if(pdim(1) < 2) nrx1_ph = 0
   if(pdim(2) < 2) nrx2_ph = 0
   if(pdim(3) < 2) nrx3_ph = 0
   !
   pdata%nrx_ph(1:3) = (/nrx1_ph, nrx2_ph, nrx3_ph/)
   !
   call init_onsite_correction(pdata)
   !
end subroutine set_polar_parameter


subroutine eph_wan_longrange(pol, xq, epmatlr, uf)
   implicit none
   type(polar_parameter), intent(in) :: pol
   real(kind=dp), intent(in) :: xq(3)
   ! the el-ph matrix element in phonon cartesian coordinate
   complex(kind=dp), intent(out) :: epmatlr(pol%nat*3)
   ! if uf present, epmatlr is in phonon eigen-displacement coordinate
   ! uf(:,i), the eigen-mode of i-th phonon mode. (has divided by sqrt(amass))
   complex(kind=dp), intent(in), optional:: uf(pol%nat*3, pol%nat*3)
   !
   integer :: nmd
   complex(kind=dp) :: ephlr(pol%nat*3)
   
   if ( pol%thickness_2d_ba > 0 ) then
      call eph_wan_longrange_2d(pol, xq, epmatlr)
   else
      call eph_wan_longrange_3d(pol, xq, epmatlr)
   endif
   
   nmd = pol%nat * 3
   if(present(uf)) then
      !tranform to phonon eigen-mode coordinate, g_mu = sum_i uf(i,mu)*epmatlr(i)
      call zgemv ('t', nmd, nmd, cone, uf, nmd, epmatlr(1:nmd), 1, czero, ephlr(1:nmd), 1)
      epmatlr(1:nmd) = ephlr(1:nmd)
   endif
end subroutine eph_wan_longrange


subroutine eph_wan_longrange_3d(pol, xq, epmatlr)
   implicit none
   type(polar_parameter), intent(in) :: pol
   ! coordinate of q-points, in crystal coordinate
   real(kind=dp), intent(in) :: xq(3)
   ! the el-ph matrix element in phonon cartesian coordinate
   complex(kind=dp), intent(out) :: epmatlr(pol%nat*3)
   !
   !unpack polar_parameter to the following variables
   integer :: nrx1, nrx2, nrx3, nat, nmodes
   real(dp) :: bg(3,3), epsil(3,3), ggmax, alpha, falph, omega, tpiba
   real(dp) :: tau(3, pol%nat), zeu(3, 3, pol%nat)
   ! local variable
   integer :: m1, m2, m3, mu, ia
   real(kind=dp) :: xqr(3), qeq, qfac, arg
   complex(kind=dp) :: prefac
   
   !unpack pol
   nat = pol%nat;    nrx1 = pol%nrx(1);   nrx2 = pol%nrx(2);   nrx3 = pol%nrx(3)
   nmodes = 3*nat;   zeu = pol%bcharge;   omega = pol%omega;   tpiba = pol%tpiba
   tau = pol%tau;    alpha = pol%alpha;   epsil = pol%epsil;   bg = pol%bg
   !
   falph = 4.0E0_dp * alpha
   ggmax = pol%gmax * falph
   
   !initialize
   epmatlr(:) = czero
   ! if xq is too small, we just return zero, it's a cutoff.
   if(sqrt( dot_product(xq(:),xq(:)) ) < 1.0E-8_dp) return
   ! sum over G
   do m1 = -nrx1,nrx1
   do m2 = -nrx2,nrx2
   do m3 = -nrx3,nrx3
      !xq is crystal coordinate of q-point. xqr is (q+G)
      xqr(1:3) = xq(1:3) + (/m1, m2, m3/)
      ! transform to cartesian coordinate (in unit of tpiba)
      call cryst_to_cart(1, xqr(1:3), bg, 1)
      ! (q+G)*epsilon*(q+G)
      qeq = dot_product(xqr(1:3), matmul(epsil(1:3,1:3), xqr(1:3)))
      !skip when (q+G)*epsilon*(q+G) equals 0, or extremely small. 
      if(qeq < 1.0E-14_dp .or. qeq > ggmax) cycle
      qfac = exp( -qeq/falph ) / qeq
      !loop over atom index, (q+G)*Z_s, s->atomic index, nmodes = 3*nat
      do ia = 1, nat
         ! phase factor,  e^{-i(q+G)*tau_k}
         arg = -tpi*dot_product(xqr, tau(:,ia))
         mu = (ia-1)*3 + 1
         epmatlr(mu:mu+2) = epmatlr(mu:mu+2) + &
            matmul(xqr, zeu(:,:,ia)) * qfac * cmplx(cos(arg),sin(arg), kind=dp)
      enddo
   enddo; enddo; enddo
   ! prefactor: 4pi*i*e^2/Omega, the tpiba is added back here, 
   ! since the dimension is (q+G)/(q+G)^2, we need to divide a tpiba 
   prefac = conei * fpi * e2 / omega / tpiba
   ! the el-ph matrix element in cartesian coordinate, Ryd atom unit.
   epmatlr(1:nmodes) = prefac * epmatlr(1:nmodes)
end subroutine


subroutine eph_wan_longrange_2d(pol, xq, epmatlr)
   implicit none
   type(polar_parameter), intent(in) :: pol
   ! coordinate of q-points, in crystal coordinate
   real(kind=dp), intent(in) :: xq(3)
   ! the el-ph matrix element in phonon cartesian coordinate
   complex(kind=dp), intent(out) :: epmatlr(pol%nat*3)
   !
   !unpack polar_parameter to the following variables
   integer :: nrx1, nrx2, nrx3, nat, nmodes
   real(dp) :: bg(3,3), epsil(3,3), ggmax, alpha, falph, omega, tpiba
   real(dp) :: tau(3, pol%nat), zeu(3, 3, pol%nat)
   ! local variable
   integer :: m1, m2, m3, mu, ia
   real(kind=dp) :: xqr(3), qdotq, qfac, arg
   complex(kind=dp) :: prefac
   ! 2d related variable
   real(kind=dp) :: alat, r, gp2, reff(2,2), epsil_2d(2,2), t, c, epsil_0(2,2)
   integer :: i, j
   !unpack pol
   nat = pol%nat;    nrx1 = pol%nrx(1);   nrx2 = pol%nrx(2);   nrx3 = pol%nrx(3)
   nmodes = 3*nat;   zeu = pol%bcharge;   omega = pol%omega;   tpiba = pol%tpiba
   tau = pol%tau;    alpha = pol%alpha;   epsil = pol%epsil;   bg = pol%bg
   !
   falph = 4.0E0_dp * alpha
   ggmax = pol%gmax * falph
   
   alat=tpi/tpiba
   c=1.0_dp/bg(3,3) ! cell thickness in alat units (in z direction)
   t=pol%thickness_2d_ba ! material thickness in alat units
   epsil_0=0.0_dp;epsil_0(1,1)=1.0_dp;epsil_0(2,2)=1.0_dp !epsilon of vacuum

   !2d effective epsilon (Eq. 10 in PRB 94, 085415)
   epsil_2d = epsil_0 + (epsil(1:2,1:2)-epsil_0) * (c/t) 

   reff=0.0d0
   DO i=1,2
      DO j= 1,2
         reff(i,j)=epsil_2d(i,j)*0.5d0*tpi*t ! (eps_2d)*t/2 in 2pi/a units
      ENDDO
   ENDDO

   !initialize
   epmatlr(:) = czero
   ! if xq is too small, we just return zero, it's a cutoff.
   if(sqrt( dot_product(xq(:),xq(:)) ) < 1.0E-8_dp) return
   ! sum over G
   do m1 = -nrx1,nrx1
   do m2 = -nrx2,nrx2
   do m3 = -nrx3,nrx3
      !xq is crystal coordinate of q-point. xqr is (q+G)
      xqr(1:3) = xq(1:3) + (/m1, m2, m3/)
      ! transform to cartesian coordinate (in unit of tpiba)
      call cryst_to_cart(1, xqr(1:3), bg, 1)
      qdotq = dot_product(xqr(:),xqr(:))
      r=0.0d0
      gp2 = dot_product(xqr(1:2),xqr(1:2))
      IF (gp2>1.0d-8) THEN
         r = dot_product(xqr(1:2), matmul(reff, xqr(1:2)))
         r=r/gp2
      ENDIF
      !skip when (q+G)*epsilon*(q+G) equals 0, or extremely small. 
      if(qdotq < 1.0E-14_dp .or. qdotq > ggmax) cycle
      qfac = exp( -qdotq/falph ) / SQRT(qdotq) / (1.0d0+r*SQRT(qdotq))
      !loop over atom index, (q+G)*Z_s, s->atomic index, nmodes = 3*nat
      do ia = 1, nat
         ! phase factor,  e^{-i(q+G)*tau_k}
         arg = -tpi*dot_product(xqr, tau(:,ia))
         mu = (ia-1)*3 + 1
         epmatlr(mu:mu+2) = epmatlr(mu:mu+2) + &
            matmul(xqr, zeu(:,:,ia)) * qfac * cmplx(cos(arg),sin(arg), kind=dp)
      enddo
   enddo; enddo; enddo
   ! prefactor: 2pi*i*e^2/A, the tpiba is added back here, 
   ! since the dimension is (q+G)/(q+G), we DON'T need to divide a tpiba 
   prefac = tpi * conei * e2 / omega*alat/bg(3,3)
   ! the el-ph matrix element in cartesian coordinate, Ryd atom unit.
   epmatlr(1:nmodes) = prefac * epmatlr(1:nmodes)
   !
end subroutine eph_wan_longrange_2d


subroutine dyn_mat_longrange(pol, xq, dynq)
   implicit none
   type(polar_parameter), intent(in) :: pol
   ! coordinate of q-points, in crystal coordinate
   real(kind=dp), intent(in) :: xq(3)
   ! contribution to the dynamical matrix due to dipole-dipole interactio.
   complex(kind=dp), intent(out) :: dynq(3, 3, pol%nat * (pol%nat+1) / 2)
   !
   integer :: ia, idx
   !
   if ( pol%thickness_2d_ba > 0 ) then
      ! two-dimensional treatment of LO-TO splitting
      call dyn_mat_longrange_2d(pol, xq, dynq)
   else
      call dyn_mat_longrange_3d(pol, xq, dynq)
   endif

   !add onsite correction to the diagonal terms (:,:, ia, ia)
   do ia = 1, pol%nat
      idx = ia * (ia+1) / 2
      dynq(:,:,idx) = dynq(:,:,idx) + pol%onsite_correction(:,:,ia)
   enddo

end subroutine dyn_mat_longrange


subroutine init_onsite_correction(pol)
   implicit none
   type(polar_parameter), intent(inout) :: pol
   !local 
   integer :: ia, ja, n, nat, nelem
   real(dp) :: q_gamma(3), dd0(3,3)
   real(dp), allocatable :: dd(:,:,:,:)
   complex(dp), allocatable :: dyn(:,:,:)

   nat = pol%nat
   nelem = nat * (nat + 1) / 2
   allocate( dyn(3, 3, nelem), dd(3, 3, nat, nat) )

   q_gamma = 0.0_dp
   if ( pol%thickness_2d_ba > 0 ) then
      ! two-dimensional treatment of LO-TO splitting
      call dyn_mat_longrange_2d(pol, q_gamma, dyn)
   else
      call dyn_mat_longrange_3d(pol, q_gamma, dyn)
   endif
      
   if(any( abs(aimag(dyn(:,:,:))) > 1.0E-16_dp )) &
      call errore('polar_correction', 'On-site polar correction is not real!',1)
   
   dd = 0.0_dp
   do ja = 1, nat
   do ia = 1, ja
      n = ja * (ja-1) / 2 + ia
      dd0(:,:) = real(dyn(:,:,n), dp)

      dd(:,:, ia, ja) = dd0
      if(ia .ne. ja) dd(:,:, ja, ia) = transpose(dd0)
   enddo; enddo
   deallocate( dyn )

   if( associated(pol%onsite_correction) ) deallocate( pol%onsite_correction )
   allocate( pol%onsite_correction(3, 3, nat) )

   do ia = 1, nat
      dd0 = 0.0_dp
      do ja = 1, nat
         dd0 = dd0 + dd(:,:, ia, ja)
      enddo
      pol%onsite_correction(:,:,ia) = -dd0
   enddo
   deallocate( dd )
end subroutine 

! similar to rigid.
subroutine dyn_mat_longrange_3d(pol, xq, dynq)
   implicit none
   type(polar_parameter), intent(in) :: pol
   ! coordinate of q-points, in crystal coordinate
   real(kind=dp), intent(in) :: xq(3)
   ! contribution to the dynamical matrix due to dipole-dipole interactio.
   complex(kind=dp), intent(out) :: dynq(3, 3, pol%nat * (pol%nat+1) / 2)
   !
   !unpack pol to the following variables
   integer :: nrx1, nrx2, nrx3, nat, nelem
   real(dp) :: bg(3,3), epsil(3,3), ggmax, alpha, falph, omega, tpiba
   real(dp) :: tau(3, pol%nat), zeu(3, 3, pol%nat)
   !
   ! local variable
   integer :: m1, m2, m3, ia, ja, i, j, n
   real(kind=dp) :: xqr(3), qeq, qfac, arg, fac
   complex(dp) :: dd0(3,3)

   !unpack pol
   nat = pol%nat;    nrx1 = pol%nrx_ph(1);   nrx2 = pol%nrx_ph(2);   nrx3 = pol%nrx_ph(3)
   bg = pol%bg;      omega = pol%omega;   tpiba = pol%tpiba;   nelem = nat*(nat+1)/2
   tau = pol%tau;    alpha = pol%alpha;   epsil = pol%epsil;   zeu = pol%bcharge
   !
   falph = 4.0E0_dp * alpha
   ggmax = pol%gmax * falph

   !initialize
   dynq(:,:,:) = czero
   ! comment out the following line to be consistent with PH/rigid.f90/rgd_blk
   !if(sqrt( dot_product(xq(:),xq(:)) ) < 1.0E-8_dp) return
   ! 4pi*e^2/Omega (3D)
   fac = fpi * e2 / omega 

   do m1 = -nrx1,nrx1
   do m2 = -nrx2,nrx2
   do m3 = -nrx3,nrx3
      !xq is crystal coordinate of q-point. xqr is (q+G)
      xqr(1:3) = xq(1:3) + (/m1, m2, m3/)
      ! transform to cartesian coordinate (in unit of tpiba)
      call cryst_to_cart(1, xqr(1:3), bg, 1)
      ! (q+G)*epsilon*(q+G)
      qeq = dot_product( xqr, matmul(epsil, xqr) )
      !skip when (q+G)*epsilon*(q+G) equals 0, or extremely small.
      if( qeq < 1.0E-14_dp .or. qeq > ggmax ) cycle
      ! (q+G)*Z vector
      !do ia = 1, nat
      !   qzeu(1:3,ia) = matmul(xqr, zeu(1:3,1:3,ia))
      !enddo
      qfac = exp( -qeq/falph ) / qeq

      ! do the summation
      n = 0
      do ja = 1, nat
      do ia = 1, ja
         n = n + 1
         ! phase factor, e^{ i(q+G)*(tau_i - tau_j) }
         arg = tpi*dot_product(xqr, (tau(:,ia)-tau(:,ja)) )
         do j = 1, 3
         do i = 1, 3
            ! collect all the contributions.
            dynq(i, j, n) = dynq(i, j, n) + &
               xqr(i) * xqr(j) * qfac * cmplx(cos(arg), sin(arg), kind=dp)
         enddo; enddo
      enddo; enddo
   enddo; enddo; enddo
  
   !multiply the Born Effective Charge
   do ja = 1, nat
   do ia = 1, ja
      n = ja * (ja-1) / 2 + ia
      dd0(:,:) = dynq(:,:, n)
      dd0 = matmul(dd0, zeu(:,:,ja))
      dynq(:,:, n) = matmul(transpose(zeu(:,:,ia)), dd0)
   enddo; enddo
   
   ! prefactor
   dynq(:,:,:) = dynq(:,:,:) * fac
end subroutine


! similar to rigid in QE/Phonon.
! two-dimensional treatment of LO-TO splitting
! results are identical with rgd_blk
subroutine dyn_mat_longrange_2d(pol, xq, dynq)
   implicit none
   type(polar_parameter), intent(in) :: pol
   ! coordinate of q-points, in crystal coordinate
   real(kind=dp), intent(in) :: xq(3)
   ! contribution to the dynamical matrix due to dipole-dipole interactio.
   complex(kind=dp), intent(out) :: dynq(3, 3, pol%nat * (pol%nat+1) / 2)
   !
   !unpack pol to the following variables
   integer :: nrx1, nrx2, nrx3, nat, nelem
   real(dp) :: bg(3,3), epsil(3,3), ggmax, alpha, falph, omega, tpiba
   real(dp) :: tau(3, pol%nat), zeu(3, 3, pol%nat)
   !
   ! local variable
   integer :: m1, m2, m3, ia, ja, i, j, n
   real(kind=dp) :: xqr(3), qeq, qfac, arg, fac
   complex(dp) :: dd0(3,3)
   !
   ! 2D variables
   real(kind=dp) :: reff(2,2), epsil_2d(2,2), gp2, r, alat

   !unpack pol
   nat = pol%nat;    nrx1 = pol%nrx_ph(1);   nrx2 = pol%nrx_ph(2);   nrx3 = pol%nrx_ph(3)
   bg = pol%bg;      omega = pol%omega;   tpiba = pol%tpiba;   nelem = nat*(nat+1)/2
   tau = pol%tau;    alpha = pol%alpha;   epsil = pol%epsil;   zeu = pol%bcharge
   !
   falph = 4.0E0_dp * alpha
   ggmax = pol%gmax * falph

   !initialize
   dynq(:,:,:) = czero
   ! comment out the following line to be consistent with PH/rigid.f90/rgd_blk
   !if(sqrt( dot_product(xq(:),xq(:)) ) < 1.0E-8_dp) return
  
   ! cell thickness is not used, in order to be consistent with rgd_blk
   alat=tpi/tpiba
   fac = fpi * e2 / omega *0.5d0*alat/bg(3,3)
   reff=0.0d0
   do i=1,2
      do j=1,2
         reff(i,j)=epsil(i,j)*0.5d0*tpi/bg(3,3)! (eps)*c/2 in 2pi/a units
      enddo
   enddo
   do i=1,2
      reff(i,i)=reff(i,i)-0.5d0*tpi/bg(3,3) ! (-1)*c/2 in 2pi/a units
   enddo

   do m1 = -nrx1,nrx1
   do m2 = -nrx2,nrx2
   do m3 = -nrx3,nrx3
      !xq is crystal coordinate of q-point. xqr is (q+G)
      xqr(1:3) = xq(1:3) + (/m1, m2, m3/)
      ! transform to cartesian coordinate (in unit of tpiba)
      call cryst_to_cart(1, xqr(1:3), bg, 1)
         
      qeq = dot_product(xqr(:),xqr(:))
      r=0.0d0
      gp2 = dot_product(xqr(1:2),xqr(1:2))
      if (gp2>1.0d-8) then
         !r=xqr(1)*reff(1,1)*xqr(1)+xqr(1)*reff(1,2)*xqr(2)+xqr(2)*reff(2,1)*xqr(1)+xqr(2)*reff(2,2)*xqr(2)
         r = dot_product( xqr(1:2), matmul(reff, xqr(1:2)) )
         r=r/gp2
      endif

      !skip when (q+G)*epsilon*(q+G) equals 0, or extremely small.
      if( qeq < 1.0E-14_dp .or. qeq > ggmax ) cycle
      ! (q+G)*Z vector
      !do ia = 1, nat
      !   qzeu(1:3,ia) = matmul(xqr, zeu(1:3,1:3,ia))
      !enddo
      qfac = exp( -qeq/falph ) / sqrt(qeq) / (1.0d0+r*sqrt(qeq))

      ! do the summation
      n = 0
      do ja = 1, nat
      do ia = 1, ja
         n = n + 1
         ! phase factor, e^{ i(q+G)*(tau_i - tau_j) }
         arg = tpi*dot_product(xqr, (tau(:,ia)-tau(:,ja)) )
         do j = 1, 3
         do i = 1, 3
            ! collect all the contributions.
            dynq(i, j, n) = dynq(i, j, n) + &
               xqr(i) * xqr(j) * qfac * cmplx(cos(arg), sin(arg), kind=dp)
         enddo; enddo
      enddo; enddo
   enddo; enddo; enddo

   !multiply the Born Effective Charge
   do ja = 1, nat
   do ia = 1, ja
      n = ja * (ja-1) / 2 + ia
      dd0(:,:) = dynq(:,:, n)
      dd0 = matmul(dd0, zeu(:,:,ja))
      dynq(:,:, n) = matmul(transpose(zeu(:,:,ia)), dd0)
   enddo; enddo

   ! prefactor
   dynq(:,:,:) = dynq(:,:,:) * fac
end subroutine

end module polar_correction

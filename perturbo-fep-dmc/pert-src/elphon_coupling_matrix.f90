!===============================================================================
! Copyright (C) 2016-2020 Jin-Jian Zhou, Jinsoo Park, I-Te Lu, Marco Bernardi
!
! This program is distributed under the terms of the GNU General Public License.
! See the file `LICENSE' in the root directory of this distribution, or obtain 
! a copy of the License at <https://www.gnu.org/licenses/gpl-3.0.txt>.
!
! Author: jjzhou <jjchou.comphy@gmail.com>
! Comment:
!
! Maintenance:
!===============================================================================

module elphon_coupling_matrix
   use pert_utils, only: get_exp_ikr
   use wigner_seitz_cell, only: ws_cell
   use polar_correction,  only: set_polar_parameter, eph_wan_longrange
   use pert_const, only: dp, twopi, czero, ci, e2, pi, cone
   use elph_matrix_wannier, only: elph_mat_wann, eph_wannier, set_ws_cell_eph
   use qe_mpi_mod, only: ionode, mp_bcast, inter_pool_comm, ionode_id, ionode, stdout
   implicit none
   public
   
   type(elph_mat_wann), save :: elph

   public :: init_elph_mat_wann, &
             load_elph_mat_wann_ep_hop, &
             release_elph_mat_wann_ep_hop, &
             release_elph_mat_wann, &
             eph_fourier_el_para, &
             eph_fourier_el, &
             eph_fourier_ph, &
             eph_fourier_elph, &
             eph_fourier_phel, & 
             eph_transform, &
             eph_transform_fast, &
             eph_transform_polar

contains

subroutine init_elph_mat_wann(file_id, kdim, qdim, ep, lep_hop)
   use hdf5_utils
   use epwan_hdf5_io, only: read_elph_mat_wann
   use pert_data, only: zstar, epsil, nat, num_wann, at, bg, volume, tpiba, &
                        tau, wannier_center_cryst, polar_alpha, thickness_2d, lpolar
   implicit none
   integer(HID_T), intent(in) :: file_id
   integer, intent(in) :: kdim(3), qdim(3)
   type(elph_mat_wann), intent(out) :: ep
   logical, intent(in), optional :: lep_hop
   !
   logical :: load_ep_hop
   integer :: ia, jw, iw
   real(dp), allocatable :: cryst_tau(:,:) 
   
   call start_clock('init_eph_wan')
   !whether or not to allocate and load the ep_hop data
   load_ep_hop = .true.
   if( present(lep_hop) ) load_ep_hop = lep_hop

   allocate( cryst_tau(3,nat) )
   !transform atomic position to crystal coordinate
   cryst_tau(:,:) = tau(:,:)
   call cryst_to_cart(nat, cryst_tau, bg, -1)
   
   call set_ws_cell_eph &
      (ep, kdim, qdim, nat, num_wann, at, cryst_tau, wannier_center_cryst, load_ep_hop)

   if(load_ep_hop) then
      if(ionode) call read_elph_mat_wann(file_id, ep)
      !bcast
      do ia = 1, ep%na
      do jw = 1, ep%nb
      do iw = 1, ep%nb
         call mp_bcast(ep%epwan(iw,jw,ia)%ep_hop(:,:,:), ionode_id, inter_pool_comm)
      enddo; enddo; enddo
   endif
   
   ep%lpol = lpolar
   ep%l_2d = merge(.true., .false., thickness_2d > 0.0_dp)
   !setup for polar correction
   if(ep%lpol) call set_polar_parameter &
      (qdim, nat, volume, tpiba, bg, epsil, zstar, tau, ep%pol, thickness_2d, polar_alpha)

   deallocate( cryst_tau )
   if(ionode) then 
      write(*,'(A10,5i5)')'size_epw = ',3*ep%na,ep%nb,ep%nb,ep%max_nre,ep%max_nrp
   endif 
   call stop_clock('init_eph_wan')
end subroutine init_elph_mat_wann

subroutine load_elph_mat_wann_ep_hop(file_id, ep)
   use hdf5_utils
   use epwan_hdf5_io, only: read_elph_mat_wann
   implicit none
   integer(HID_T), intent(in) :: file_id
   type(elph_mat_wann), intent(inout) :: ep
   !
   integer :: ia, jw, iw, nre, nrp
   type(eph_wannier), pointer :: ptr
   !
   !allocate space
   do ia = 1, ep%na
   do jw = 1, ep%nb
   do iw = 1, ep%nb
      ptr => ep%epwan(iw, jw, ia)
      !
      nre = ptr%ws_el%nr
      nrp = ptr%ws_ph%nr
      if(.not. associated(ptr%ep_hop)) allocate( ptr%ep_hop(3, nre, nrp) )
   enddo; enddo; enddo
   !load data   
   if(ionode) call read_elph_mat_wann(file_id, ep)
   !bcast
   do ia = 1, ep%na
   do jw = 1, ep%nb
   do iw = 1, ep%nb
      call mp_bcast(ep%epwan(iw,jw,ia)%ep_hop(:,:,:), ionode_id, inter_pool_comm)
   enddo; enddo; enddo
end subroutine


subroutine release_elph_mat_wann_ep_hop(ep)
   implicit none
   type(elph_mat_wann), intent(inout) :: ep
   !
   integer :: ia, jw, iw
   type(eph_wannier), pointer :: ptr

   if( associated(ep%epwan) ) then
      do ia = 1, ep%na
      do jw = 1, ep%nb
      do iw = 1, ep%nb
         ptr => ep%epwan(iw, jw, ia)
         if( associated(ptr%ep_hop) ) deallocate( ptr%ep_hop )
      enddo; enddo; enddo
   endif
end subroutine


subroutine release_elph_mat_wann(ep)
   implicit none
   type(elph_mat_wann), intent(inout) :: ep
   !
   integer :: ia, jw, iw
   type(eph_wannier), pointer :: ptr
   
   if( associated(ep%rvec_set_el) ) deallocate( ep%rvec_set_el )
   if( associated(ep%rvec_set_ph) ) deallocate( ep%rvec_set_ph )
   
   if( associated(ep%epwan) ) then
      do ia = 1, ep%na
      do jw = 1, ep%nb
      do iw = 1, ep%nb
         ptr => ep%epwan(iw, jw, ia)
         if( associated(ptr%ep_hop) ) deallocate( ptr%ep_hop )
         !
         if( associated(ptr%ws_el%ndeg) ) deallocate( ptr%ws_el%ndeg )
         if( associated(ptr%ws_el%rvec) ) deallocate( ptr%ws_el%rvec )
         !
         if( associated(ptr%ws_ph%ndeg) ) deallocate( ptr%ws_ph%ndeg )
         if( associated(ptr%ws_ph%rvec) ) deallocate( ptr%ws_ph%rvec )
      enddo; enddo; enddo
      deallocate( ep%epwan )
   endif
   
   if( associated(ep%pol%bcharge) ) deallocate( ep%pol%bcharge )
   if( associated(ep%pol%tau) ) deallocate( ep%pol%tau )
end subroutine


!sum_{Re} <0,n| dV_{Rp,mod} |Re,m> * e^{ikRe}
subroutine eph_fourier_el_para(ep, kpt, g_kerp)
   implicit none
   type(elph_mat_wann), intent(in) :: ep
   real(dp), intent(in) :: kpt(3)
   complex(dp), intent(out) :: g_kerp(3, ep%max_nrp, ep%nb, ep%nb, ep%na)
   !blas function zdotu: product of two complex vector: sum(a(:) * b(:))
   !NOTE: use zdotc or dot_product to perform sum(conjg(a(:)) * b(:))
   complex(dp), external :: zdotu
   complex(dp) :: cfac_tmp(ep%max_nre)
   complex(dp), allocatable :: exp_ikr(:) 
   !local variables
   integer :: idx, nelem, ia, jw, iw, ir, i
   type(eph_wannier), pointer :: pep
      
   nelem = ep%nb * ep%nb * ep%na
   !
   allocate( exp_ikr( ep%nrvec_e ) )
   exp_ikr = czero
   call get_exp_ikr(kpt, ep%rvec_set_el, exp_ikr)

   g_kerp = czero
!$omp parallel do schedule(guided) default(shared) private(idx, ia, jw, iw, pep, cfac_tmp, ir, i)
   do idx = 0, nelem-1
      !idx = (ia-1) * nb*nb + (jw-1)*nb + (iw-1)
      ia = idx / (ep%nb * ep%nb) + 1
      jw = mod( idx/ep%nb, ep%nb ) + 1
      iw = mod( idx, ep%nb ) + 1
      !
      pep => ep%epwan(iw, jw, ia)
      call copy_exp_ikr(pep%ws_el, exp_ikr, cfac_tmp(1))
      !
      do ir = 1, pep%ws_ph%nr
         do i = 1, 3
            g_kerp(i,ir,iw,jw,ia) = zdotu(pep%ws_el%nr, cfac_tmp(1), 1, pep%ep_hop(i,1,ir), 3)
         enddo
      enddo
   enddo
!$omp end parallel do
   !
   deallocate( exp_ikr )
end subroutine eph_fourier_el_para


!sum_{Re} <0,n| dV_{Rp,mod} |Re,m> * e^{ikRe}
subroutine eph_fourier_el(ep, kpt, g_kerp)
   implicit none
   type(elph_mat_wann), intent(in) :: ep
   real(dp), intent(in) :: kpt(3)
   complex(dp), intent(out) :: g_kerp(3, ep%max_nrp, ep%nb, ep%nb, ep%na)
   !blas function zdotu: product of two complex vector: sum(a(:) * b(:))
   !NOTE: use zdotc or dot_product to perform sum(conjg(a(:)) * b(:))
   complex(dp), external :: zdotu
   complex(dp) :: cfac_tmp(ep%max_nre)
   complex(dp), allocatable :: exp_ikr(:) 
   !local variables
   integer :: ia, jw, iw, ir, i
   type(eph_wannier), pointer :: pep
   
   allocate( exp_ikr( ep%nrvec_e ) )
   exp_ikr = czero
   call get_exp_ikr(kpt, ep%rvec_set_el, exp_ikr)
      
   g_kerp = czero
   do ia = 1, ep%na
   do jw = 1, ep%nb
   do iw = 1, ep%nb
      pep => ep%epwan(iw, jw, ia)
      call copy_exp_ikr(pep%ws_el, exp_ikr, cfac_tmp(1))
      !
      do ir = 1, pep%ws_ph%nr
         do i = 1, 3
            g_kerp(i,ir,iw,jw,ia) = zdotu(pep%ws_el%nr, cfac_tmp(1), 1, pep%ep_hop(i,1,ir), 3)
         enddo
      enddo
      !
   enddo; enddo; enddo
   !
   deallocate( exp_ikr )
end subroutine eph_fourier_el

!sum_{Rp} <0,n| dV_{Rp,mod} |Re,m> * e^{ikRp}
subroutine eph_fourier_ph(ep, kpt, g_rekp)
   implicit none
   type(elph_mat_wann), intent(in) :: ep
   real(dp), intent(in) :: kpt(3)
   complex(dp), intent(out) :: g_rekp(3, ep%max_nre, ep%nb, ep%nb, ep%na)
   !blas function zdotu: product of two complex vector: sum(a(:) * b(:))
   !NOTE: use zdotc or dot_product to perform sum(conjg(a(:)) * b(:))
   complex(dp), external :: zdotu
   complex(dp) :: cfac_tmp(ep%max_nrp), ctemp(ep%max_nrp)
   complex(dp), allocatable :: exp_ikr(:) 
   !local variables
   integer :: ia, jw, iw, ir, i
   type(eph_wannier), pointer :: pep
   
   allocate( exp_ikr( ep%nrvec_p ) )
   exp_ikr = czero
   call get_exp_ikr(kpt, ep%rvec_set_ph, exp_ikr)
      
   g_rekp = czero
   do ia = 1, ep%na
   do jw = 1, ep%nb
   do iw = 1, ep%nb
      pep => ep%epwan(iw, jw, ia)
      call copy_exp_ikr(pep%ws_ph, exp_ikr, cfac_tmp(1))
      !
      do ir = 1, pep%ws_el%nr
         do i = 1, 3
            ctemp(1:pep%ws_ph%nr) = pep%ep_hop(i,ir,1:pep%ws_ph%nr)
            g_rekp(i,ir,iw,jw,ia) = zdotu( pep%ws_ph%nr, cfac_tmp(1), 1, ctemp, 1 )
            !g_rekp(i,ir,iw,jw,ia) = sum(pep%ws_el%nr, cfac_tmp(1), ctemp,  )
         enddo
      enddo
      !
   enddo; enddo; enddo
   !
   deallocate( exp_ikr )
end subroutine eph_fourier_ph

subroutine copy_exp_ikr(ws, exp_ikr_all, exp_ikr_loc)
   use wigner_seitz_cell, only: ws_cell
   implicit none
   type(ws_cell), intent(in) :: ws
   complex(dp), intent(in) :: exp_ikr_all(:)
   complex(dp), intent(out) :: exp_ikr_loc(ws%nr)
   !
   integer :: i, ir

   do i = 1, ws%nr
      ir = ws%rvec(i)
      exp_ikr_loc(i) = exp_ikr_all(ir)
   enddo
end subroutine copy_exp_ikr


!<k+q,n|dV_{q,mod}|k,m> = sum_{Re,Rp} <0,n|dV_{Rp,mod}|Re,m> e^{ikRe} e^{iqRp}
!where |k,m> = sum_{Re} |Re,m> e^{ikRe} is bloch wavefunction in wannier gauge.
! dV_{q,mod} = sum_{Rp} dV_{Rp,mod} e^{iqRp}, mod is displacement cart. coord.
subroutine eph_fourier_elph(ep, xpt, g_kerp, gkq)
   implicit none
   type(elph_mat_wann), intent(in) :: ep
   real(dp), intent(in) :: xpt(3)
   !gmat can be either g_rekp or g_kerp
   complex(dp), intent(in) :: g_kerp(3, ep%max_nrp, ep%nb, ep%nb, ep%na)
   complex(dp), intent(out) :: gkq(ep%nb, ep%nb, 3*ep%na)
   !
   integer :: ia, jw, iw, ir, irp, i, idx
   complex(dp) :: gkq_tmp(3)
   complex(dp), allocatable :: exp_ikr(:)
   type(ws_cell), pointer :: ws_ph

   allocate( exp_ikr( ep%nrvec_p ) )
   exp_ikr = czero
   call get_exp_ikr(xpt, ep%rvec_set_ph, exp_ikr)
   
   do ia = 1, ep%na
   do jw = 1, ep%nb
   do iw = 1, ep%nb
      ws_ph => ep%epwan(iw, jw, ia)%ws_ph
      !
      gkq_tmp = czero
      do ir = 1, ws_ph%nr
         irp = ws_ph%rvec(ir)
         gkq_tmp(1:3) = gkq_tmp(1:3) + exp_ikr(irp) * g_kerp(1:3,ir,iw,jw,ia)
      enddo
      
      do i = 1, 3
         idx = (ia-1)*3 + i
         gkq(iw, jw, idx) = gkq_tmp(i)
      enddo
   enddo; enddo; enddo
   
   deallocate( exp_ikr )
end subroutine eph_fourier_elph

!<k+q,n|dV_{q,mod}|k,m> = sum_{Re,Rp} <0,n|dV_{Rp,mod}|Re,m> e^{ikRe} e^{iqRp}
!where |k,m> = sum_{Re} |Re,m> e^{ikRe} is bloch wavefunction in wannier gauge.
! dV_{q,mod} = sum_{Rp} dV_{Rp,mod} e^{iqRp}, mod is displacement cart. coord.
subroutine eph_fourier_phel(ep, xpt, g_rekp, gkq)
   implicit none
   type(elph_mat_wann), intent(in) :: ep
   real(dp), intent(in) :: xpt(3)
   !gmat can be either g_rekp or g_kerp
   complex(dp), intent(in) :: g_rekp(3, ep%max_nre, ep%nb, ep%nb, ep%na)
   complex(dp), intent(out) :: gkq(ep%nb, ep%nb, 3*ep%na)
   !
   integer :: ia, jw, iw, ir, irp, ire, i, idx
   complex(dp) :: gkq_tmp(3)
   complex(dp), allocatable :: exp_ikr(:)
   type(ws_cell), pointer :: ws_ph, ws_el

   allocate( exp_ikr( ep%nrvec_e ) )
   exp_ikr = czero
   call get_exp_ikr(xpt, ep%rvec_set_el, exp_ikr)
   
   do ia = 1, ep%na
   do jw = 1, ep%nb
   do iw = 1, ep%nb
      !ws_ph => ep%epwan(iw, jw, ia)%ws_ph
      ws_el => ep%epwan(iw, jw, ia)%ws_el
      !
      gkq_tmp = czero
      do ir = 1, ws_el%nr
         ire = ws_el%rvec(ir)
         gkq_tmp(1:3) = gkq_tmp(1:3) + exp_ikr(ire) * g_rekp(1:3,ir,iw,jw,ia)
      enddo
      
      do i = 1, 3
         idx = (ia-1)*3 + i
         gkq(iw, jw, idx) = gkq_tmp(i)
      enddo
   enddo; enddo; enddo
   
   deallocate( exp_ikr )
end subroutine eph_fourier_phel

!transfer from wannier gauge to H eigen-states gauge
!   |i-th eig-state> = sum_{m} |k,m> U_{m,i}(k), 
! where U_{:,i} is the i-th eigen-vector of H^(W) = <k,n| H |k,m>
!.........from displace coord. to phonon eigen-modes
!    dV_{q,mu} = sum_{mod} dV_{q,mod} V_{mod,mu}(q)
! where V_{:,mu} is the mu-th phonon-eigenmode of D(q)
! to summary, the subroutine below does the following transformation. 
!     U_(k+q)^\dagger <k+q,n|dV_{q,mod}|k,m> V(q) U(k) 
subroutine eph_transform(ep, qpt, mq, uk, ukq, gkq)
   implicit none
   type(elph_mat_wann), intent(in) :: ep
   real(dp), intent(in) :: qpt(3) !for polar correction
   !mq(:,i): i-th phonon eigen-displacement
   complex(dp), intent(in) :: mq(3*ep%na, 3*ep%na)
   !uk(:,i) : the i-th eigenvector of H(k)
   complex(dp), intent(in) :: uk(ep%nb, ep%nb), ukq(ep%nb, ep%nb)
   !el-ph matrx elements in the basis of electron- and phonon- eigen-states 
   complex(dp), intent(inout) :: gkq(ep%nb, ep%nb, 3*ep%na)
   !local 
   integer :: i, j
   complex(dp) :: ctmp(3*ep%na), gtmp(ep%nb, ep%nb), ukq_h(ep%nb, ep%nb)
   
   !mq(:,j): the j-th phonon eigen-mode
   !transfer to phonon mode coord.: gkq'(j) = sum_i gkq(i)*mq(i,j)
   do j = 1, ep%nb
      do i = 1, ep%nb
         ctmp = gkq(i, j, :)
         gkq(i, j, :) = matmul(ctmp, mq)
      enddo
   enddo

   !polar correction for polar materials: add the long-range part.
   if(ep%lpol) then
      call eph_wan_longrange(ep%pol, qpt, ctmp, mq)
      do j = 1, 3*ep%na
         do i = 1, ep%nb
            gkq(i,i,j) = gkq(i,i,j) + ctmp(j)
         enddo
      enddo
   endif

   !rotation: transfer from wannier gauge to bloch gauge.
   ! g^(H) = U_kq^\degger g^(W) U_k
   ukq_h = transpose(dconjg(ukq))
   do i = 1, 3*ep%na
      gtmp = matmul(ukq_h, gkq(:,:,i))
      gkq(:,:,i) = matmul(gtmp, uk)
   enddo
end subroutine eph_transform



subroutine eph_transform_fast(ep, uq, uk, ukq, gkq, epmatlr)
   implicit none
   type(elph_mat_wann), intent(in) :: ep
   !uq(:,i): i-th phonon eigen-displacement
   complex(dp), intent(in) :: uq(3*ep%na, 3*ep%na)
   !uk(:,i) : the i-th eigenvector of H(k)
   complex(dp), intent(in) :: uk(ep%nb, ep%nb), ukq(ep%nb, ep%nb)
   !el-ph matrx elements in the basis of electron- and phonon- eigen-states 
   complex(dp), intent(inout) :: gkq(ep%nb, ep%nb, 3*ep%na)
   !long-range e-ph coupling in wannier, phonon eigen-displacement basis
   complex(dp), intent(in), optional :: epmatlr(3*ep%na)
   !local 
   integer :: nb, nm, nb2, i
   complex(dp) :: gtmp(ep%nb, ep%nb, 3*ep%na)
   
   nb = ep%nb;   nb2 = ep%nb*ep%nb;   nm = 3*ep%na
   gtmp = czero
   if( present(epmatlr) ) then
      do i = 1, nb
         gtmp(i,i,:) = epmatlr(:)
      enddo
   endif

   !uq(:,j): the j-th phonon eigen-mode
   !transfer to phonon mode coord and add polar correction: 
   !      gkq'(j) = sum_i gkq(i)*uq(i,j) + gpol(j)
   call zgemm('n','n', nb2, nm, nm, cone, gkq, nb2, uq, nm, cone, gtmp, nb2)
   
   !rotation: transfer from wannier gauge to bloch gauge.
   ! g^(H) = U_kq^\degger g^(W) U_k
   call zgemm('c','n', nb, nb*nm, nb, cone, ukq, nb, gtmp, nb, czero, gkq, nb)
   do i = 1, nm
      gkq(:,:,i) = matmul(gkq(:,:,i), uk)
   enddo
end subroutine eph_transform_fast


!for computing only the polar part of e-ph matrix elements, 
!taking epmatlr from eph_wan_longrange and then perform the transformation:
!    g^(H) = U_kq^\degger g^(W) U_k  
!to get the polar part the e-ph matrix elements.
subroutine eph_transform_polar(ep, epmatlr, uk, ukq, gkq_pol)
   implicit none
   type(elph_mat_wann), intent(in) :: ep
   !the long-range el-ph matrix elements in phonon cartesian coordinate
   complex(dp), intent(in) :: epmatlr(3*ep%na) !epmatlr(3*nat)
   !uk(:,i) : the i-th eigenvector of H(k)
   complex(dp), intent(in) :: uk(ep%nb, ep%nb), ukq(ep%nb, ep%nb)
   !el-ph matrx elements in the basis of electron- and phonon- eigen-states 
   complex(dp), intent(out) :: gkq_pol(ep%nb, ep%nb, 3*ep%na)
   !local 
   integer :: i, nb, nm
   complex(dp) :: gtmp(ep%nb, ep%nb, 3*ep%na)
   
   gkq_pol = czero
   nb = ep%nb;    nm = 3*ep%na
   
   gtmp = czero
   do i = 1, nb
      gtmp(i,i,:) = epmatlr(:)
   enddo

   !g^(H) = U_kq^\degger g^(W) U_k
   call zgemm('c','n', nb, nb*nm, nb, cone, ukq, nb, gtmp, nb, czero, gkq_pol, nb)
   do i = 1, nm
      gkq_pol(:,:,i) = matmul(gkq_pol(:,:,i), uk)
   enddo
end subroutine eph_transform_polar



end module elphon_coupling_matrix

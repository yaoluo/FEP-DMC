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

module elph_matrix_wannier
   use kinds, only: dp
   use polar_correction, only: polar_parameter
   use wigner_seitz_cell, only: ws_cell, vector_set, init_rvec_images, &
      reset_rvec_images, set_wigner_seitz_cell
   implicit none
   private

   type, public :: elph_mat_wann
      integer :: na  !number of atoms
      integer :: nb  !number of bands (or wannier orbitals)
      !
      integer :: max_nre ! max value of epwan(:,:,:)%ep_hop(:)%ws_el%nr
      integer :: max_nrp

      integer :: nrvec_e
      integer :: nrvec_p
      real(dp), pointer :: rvec_set_el(:,:) => null()
      real(dp), pointer :: rvec_set_ph(:,:) => null()
      
      !e-ph matrix elements in wannier basis: <0,n| dV_{Rp,mod} |Re,m>
      type(eph_wannier), pointer :: epwan(:,:,:) => null() ! epwan(nb, nb, na)

      !for polar system
      logical :: lpol
      logical :: l_2d !.true. if polar correction for 2D system
      type(polar_parameter) :: pol
   end type

   type, public :: eph_wannier
      type(ws_cell), pointer :: ws_ph  => null()
      type(ws_cell), pointer :: ws_el  => null()
      complex(dp),   pointer :: ep_hop(:,:,:) => null() ! ep_hop(3, ws_el%nr, ws_ph%nr)
   end type

   public :: set_ws_cell_eph
contains

subroutine set_ws_cell_eph(ep, nk, nq, nat, nwan, at, cryst_tau, wan_center, allocate_ep)
   implicit none
   type(elph_mat_wann), intent(inout) :: ep
   integer, intent(in) :: nk(3), nq(3), nat, nwan
   !cryst_tau: atomic position in crystal coordinate
   real(dp), intent(in) :: at(3,3), cryst_tau(3, nat), wan_center(3, nwan)
   logical, intent(in), optional :: allocate_ep
   !
   logical :: allocate_ephop
   integer :: ia, jw, iw, nre, nrp, max_evec, max_pvec
   type(vector_set) :: rvec_el_images, rvec_ph_images
   type(eph_wannier), pointer :: ptr
   
   if(associated(ep%rvec_set_el) .or. &
      associated(ep%rvec_set_ph) .or. associated(ep%epwan) ) &
      call errore('set_ws_cell_eph','elph_mat_wann object is already initialized',1)
   
   !whether or not pre-allocate space for ep_hop.
   allocate_ephop = .true.
   if( present(allocate_ep) ) allocate_ephop = allocate_ep
   
   ep%na = nat
   ep%nb = nwan
   allocate( ep%epwan(nwan, nwan, nat) )
   !
   call init_rvec_images(nk(1), nk(2), nk(3), at, rvec_el_images)
   call init_rvec_images(nq(1), nq(2), nq(3), at, rvec_ph_images)
   
   max_evec = 0
   max_pvec = 0
   ! ws cell for el-ph coupling
   do ia = 1, ep%na
   do jw = 1, ep%nb
   do iw = 1, ep%nb
      ptr => ep%epwan(iw, jw, ia)
      
      allocate( ptr%ws_el )
      call set_wigner_seitz_cell( nk(1), nk(2), nk(3), at, rvec_el_images, &
         ptr%ws_el, wan_center(1:3,iw), wan_center(1:3,jw) )

      allocate( ptr%ws_ph )
      call set_wigner_seitz_cell( nq(1), nq(2), nq(3), at, rvec_ph_images, &
         ptr%ws_ph, wan_center(1:3,iw), cryst_tau(1:3,ia) )

      nre = ptr%ws_el%nr
      nrp = ptr%ws_ph%nr
      if(max_evec < nre) max_evec = nre
      if(max_pvec < nrp) max_pvec = nrp
      
      !allocate space
      if(allocate_ephop)  allocate( ptr%ep_hop(3, nre, nrp) )
   enddo; enddo; enddo

   ep%max_nre = max_evec
   ep%max_nrp = max_pvec

   call setup_rvec_set_eph(ep, rvec_el_images, rvec_ph_images)

   call reset_rvec_images(rvec_el_images)
   call reset_rvec_images(rvec_ph_images)
end subroutine set_ws_cell_eph


subroutine setup_rvec_set_eph(ep, r_images_el, r_images_ph)
   implicit none
   type(elph_mat_wann), intent(inout) :: ep
   type(vector_set), intent(in) :: r_images_el, r_images_ph
   !
   integer :: ia, jw, iw, ir, ire, irp
   integer, allocatable :: rvec_label_el(:), rvec_label_ph(:)
   type(eph_wannier), pointer :: ptr

   allocate( rvec_label_el( r_images_el%nvec ) )
   allocate( rvec_label_ph( r_images_ph%nvec ) )
   rvec_label_el = 0
   rvec_label_ph = 0

   do ia = 1, ep%na
   do jw = 1, ep%nb
   do iw = 1, ep%nb
      ptr => ep%epwan(iw, jw, ia)

      do ir = 1, ptr%ws_el%nr
         ire = ptr%ws_el%rvec(ir)
         rvec_label_el(ire) = rvec_label_el(ire) + 1
      enddo
      
      do ir = 1, ptr%ws_ph%nr
         irp = ptr%ws_ph%rvec(ir)
         rvec_label_ph(irp) = rvec_label_ph(irp) + 1
      enddo
   enddo; enddo; enddo

   ep%nrvec_e = count( rvec_label_el > 0 )
   ep%nrvec_p = count( rvec_label_ph > 0 )
   
   allocate( ep%rvec_set_el(3, ep%nrvec_e) )
   allocate( ep%rvec_set_ph(3, ep%nrvec_p) )
   ep%rvec_set_el = 0.0E0_dp
   ep%rvec_set_ph = 0.0E0_dp

   ir = 0
   do ire = 1, r_images_el%nvec
      if( rvec_label_el(ire) > 0 ) then
         ir = ir + 1
         ep%rvec_set_el(1:3, ir) = real(r_images_el%vec(1:3, ire), dp)
         rvec_label_el(ire) = ir
      endif
   enddo
   
   ir = 0
   do irp = 1, r_images_ph%nvec
      if( rvec_label_ph(irp) > 0 ) then
         ir = ir + 1
         ep%rvec_set_ph(1:3, ir) = real(r_images_ph%vec(1:3, irp), dp)
         rvec_label_ph(irp) = ir
      endif
   enddo

   !update ws_el and ws_ph
   do ia = 1, ep%na
   do jw = 1, ep%nb
   do iw = 1, ep%nb
      ptr => ep%epwan(iw, jw, ia)

      do ir = 1, ptr%ws_el%nr
         ire = ptr%ws_el%rvec(ir)
         ptr%ws_el%rvec(ir) = rvec_label_el(ire)
      enddo
      
      do ir = 1, ptr%ws_ph%nr
         irp = ptr%ws_ph%rvec(ir)
         ptr%ws_ph%rvec(ir) = rvec_label_ph(irp)
      enddo
   enddo; enddo; enddo

   deallocate(rvec_label_el, rvec_label_ph)
end subroutine setup_rvec_set_eph

end module elph_matrix_wannier

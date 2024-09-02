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

module force_constant
   use kinds, only: dp
   use polar_correction, only: polar_parameter
   use wigner_seitz_cell, only: ws_cell, vector_set, init_rvec_images, &
      reset_rvec_images,  set_wigner_seitz_cell
   implicit none
   private

   type, public :: lattice_ifc
      integer :: na  ! number of atoms per unit cell
      integer :: nm  ! number of modes = 3*na
      ! phi( na*(na+1)/2 )
      type(force_const), pointer :: phi(:) => null()
      !
      integer :: max_nr !  max value of phi(:)%ws_ph%nr
      !
      integer :: nrvec
      real(dp), pointer :: rvec_set(:,:) => null() ! rvec_set(3,nrvec)
      !type(vector_set), pointer :: rvec_set  => null()
      !
      !for polar system
      logical :: lpol ! true if the system is polar
      logical :: l_2d ! true if polar correction for 2D
      type(polar_parameter) :: pol
   end type
   
   type, public :: force_const
      type(ws_cell), pointer :: ws_ph => null()
      complex(dp), pointer :: ifc(:,:,:) => null() ! ifc(3, 3, ws_ph%nr)
   end type
   
   public :: set_ws_cell_ph
contains

subroutine set_ws_cell_ph(ph, qdim, nat, at, cryst_tau)
   implicit none
   type(lattice_ifc), intent(inout) :: ph
   integer, intent(in) :: qdim(3), nat
   ! atomic position in crystal coordinate
   real(dp), intent(in) :: at(3,3), cryst_tau(3,nat) 
   !
   integer :: nelem, m, ia, ja, nvec, max_nvec
   type(vector_set) :: rvec_images

   ! make sure ph is an empty object
   if(associated(ph%phi) .or. associated(ph%rvec_set)) &
      call errore('set_ws_cell_ph','ws_cell object is already initialized',1)
   
   ph%na = nat
   ph%nm = 3*nat
   nelem = nat * (nat + 1) / 2
   !setup all the the vector set for wigner seitz cell
   allocate( ph%phi(nelem) )
   call init_rvec_images( qdim(1), qdim(2), qdim(3), at, rvec_images )

   max_nvec = 0
   m = 0
   do ja = 1, nat
   do ia = 1, ja
      !m = (ja * (ja - 1)) / 2 + ia
      m = m + 1
      allocate( ph%phi(m)%ws_ph )
      !
      call set_wigner_seitz_cell( qdim(1), qdim(2), qdim(3), at, rvec_images, &
         ph%phi(m)%ws_ph, cryst_tau(:,ia), cryst_tau(:,ja) )

      !allocate space and init ham_r
      nvec = ph%phi(m)%ws_ph%nr
      allocate( ph%phi(m)%ifc(3, 3, nvec) )
      !ph%phi(m)%ifc(:,:,:) = cmplx(0.0_dp, 0.0_dp, kind=dp)

      if( max_nvec < nvec ) max_nvec = nvec
   enddo; enddo
   !
   ph%max_nr = max_nvec
   
   call setup_rvec_set_ph(ph, rvec_images)

   call reset_rvec_images( rvec_images )
end subroutine set_ws_cell_ph


subroutine setup_rvec_set_ph(ph, r_images)
   implicit none
   type(lattice_ifc), intent(inout) :: ph
   type(vector_set), intent(in) :: r_images
   !
   integer :: m, ir, irp, nelem
   integer, allocatable :: rvec_label(:)
   type(force_const), pointer :: ptr
   
   nelem = ph%na * (ph%na + 1) / 2
   !set up rvec_set
   allocate( rvec_label( r_images%nvec ) )

   rvec_label = 0
   do m = 1, nelem
      ptr => ph%phi(m)
      
      do ir = 1, ptr%ws_ph%nr
         irp = ptr%ws_ph%rvec(ir)
         rvec_label( irp ) = rvec_label( irp ) + 1
      enddo
   enddo
   
   ph%nrvec = count(rvec_label > 0)
   allocate( ph%rvec_set(3, ph%nrvec) )
   ph%rvec_set = 0.0E0_dp
   
   m = 0
   do irp = 1, r_images%nvec
      !this R is used in the wigner seitz cell
      if( rvec_label(irp) > 0 ) then
         m = m + 1
         ph%rvec_set(1:3, m) = real(r_images%vec(1:3, irp), dp)
         rvec_label(irp) = m
      endif
   enddo
   
   !update ws_ph%rvec
   do m= 1, nelem
      ptr => ph%phi(m)

      do ir = 1, ptr%ws_ph%nr
         irp = ptr%ws_ph%rvec(ir)
         ptr%ws_ph%rvec(ir) = rvec_label(irp)
      enddo
   enddo

   deallocate( rvec_label )
end subroutine setup_rvec_set_ph

end module force_constant

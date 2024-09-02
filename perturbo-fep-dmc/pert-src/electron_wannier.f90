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

module electron_wannier
   use kinds, only: dp
   use wigner_seitz_cell, only: ws_cell, vector_set, init_rvec_images, &
      reset_rvec_images, set_wigner_seitz_cell
   implicit none
   private

   type, public :: electron_wann
      integer :: nb   !number of bands (or wannier functions)
      !
      integer :: max_nr ! max value of ham_r(:)%ws_el%nr
      !
      type(hopping), pointer :: ham_r(:) => null() ! ham_r( nb*(nb+1)/2 )
      !
      integer :: nelem ! nelem = nb*(nb+1)/2, not used, for data alignment
      !
      integer :: nrvec
      ! vector R in crystal coordiante
      real(dp), pointer :: rvec_set(:,:) => null() !rvec_set(3, nrvec)
      ! vector in rvec_set in cartesian coordinate
      real(dp), pointer :: rset_cart(:,:) => null()
   end type electron_wann

   type, public :: hopping
      type(ws_cell), pointer :: ws_el => null()
      complex(dp),   pointer :: hop(:) => null()
   end type

   public :: set_ws_cell_el
contains
   
subroutine set_ws_cell_el(el, nk, nwan, at, wan_center)
   implicit none
   type(electron_wann), intent(inout) :: el
   integer, intent(in) :: nk(3), nwan
   ! wannier center in crystal coordinate
   real(dp), intent(in) :: at(3,3), wan_center(3, nwan)
   !local
   integer :: nelem, m, ib, jb, nvec, max_nvec
   type(vector_set) :: rvec_images

   ! make sure el is an empty object
   if(associated(el%ham_r) .or. associated(el%rvec_set)) &
      call errore('set_ws_cell_el','electron_wann object is already initialized',1)
   
   el%nb = nwan
   nelem = nwan * (nwan + 1) / 2
   el%nelem = nelem
   !setup all the the vector set for wigner seitz cell
   allocate( el%ham_r(nelem) )
   call init_rvec_images( nk(1), nk(2), nk(3), at, rvec_images )

   max_nvec = 0
   m = 0
   do jb = 1, nwan
   do ib = 1, jb
      !m = (jb * (jb - 1)) / 2 + ib
      m = m + 1
      allocate( el%ham_r(m)%ws_el )
      !
      call set_wigner_seitz_cell( nk(1), nk(2), nk(3), at, rvec_images, &
         el%ham_r(m)%ws_el, wan_center(:,ib), wan_center(:,jb) )

      !allocate space and init ham_r
      nvec = el%ham_r(m)%ws_el%nr
      allocate( el%ham_r(m)%hop( nvec ) )
      !el%ham_r(m)%hop(:) = cmplx(0.0_dp, 0.0_dp, kind=dp)

      if( max_nvec < nvec ) max_nvec = nvec
   enddo; enddo
   !
   el%max_nr = max_nvec

   call setup_rvec_set_el(el, rvec_images)

   call reset_rvec_images(rvec_images)
end subroutine set_ws_cell_el


subroutine setup_rvec_set_el(el, r_images)
   implicit none
   type(electron_wann), intent(inout) :: el
   type(vector_set), intent(in) :: r_images
   !
   integer :: m, ir, ire, nelem
   integer, allocatable :: rvec_label(:)
   type(hopping), pointer :: ptr
   
   nelem = el%nb * (el%nb + 1) / 2
   !set up rvec_set
   allocate( rvec_label( r_images%nvec ) )

   rvec_label = 0
   do m = 1, nelem
      ptr => el%ham_r(m)
      
      do ir = 1, ptr%ws_el%nr
         ire = ptr%ws_el%rvec(ir)
         rvec_label( ire ) = rvec_label( ire ) + 1
      enddo
   enddo
   
   el%nrvec = count(rvec_label > 0)
   allocate( el%rvec_set(3, el%nrvec) )
   el%rvec_set = 0.0E0_dp
   
   m = 0
   do ire = 1, r_images%nvec
      !this R is used in the wigner seitz cell
      if( rvec_label(ire) > 0 ) then
         m = m + 1
         el%rvec_set(1:3, m) = real(r_images%vec(1:3, ire), dp)
         rvec_label(ire) = m
      endif
   enddo
   
   !update ws_el%rvec
   do m= 1, nelem
      ptr => el%ham_r(m)

      do ir = 1, ptr%ws_el%nr
         ire = ptr%ws_el%rvec(ir)
         ptr%ws_el%rvec(ir) = rvec_label(ire)
      enddo
   enddo

   deallocate( rvec_label )
end subroutine setup_rvec_set_el

end module electron_wannier

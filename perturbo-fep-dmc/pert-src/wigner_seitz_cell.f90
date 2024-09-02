!===============================================================================
! Copyright (C) 2016-2020 Jin-Jian Zhou, Jinsoo Park, I-Te Lu, Marco Bernardi
!
! This program is distributed under the terms of the GNU General Public License.
! See the file `LICENSE' in the root directory of this distribution, or obtain 
! a copy of the License at <https://www.gnu.org/licenses/gpl-3.0.txt>.
!
! Author: jjzhou <jjchou.comphy@gmail.com>
! Comment:
!  setup the Wigner Seitz supercells for wannier interpolation.
!
! Maintenance:
!===============================================================================

module wigner_seitz_cell
   use kinds, only: dp
   implicit none
   public

   type :: ws_cell
      integer :: nr
      integer :: reserved_ = 0 !not used, for data alignment.
      integer, pointer :: ndeg(:) => null() !ndeg(nr)
      integer, pointer :: rvec(:) => null() !rvec(nr)
   end type

   type :: vector_set
      integer :: nvec
      integer :: reserved_ = 0 !not used, for data alignment.
      integer, pointer :: vec(:,:) => null()  !vec(3,nvec)
      integer, pointer :: idx(:) => null()  ! idx(nr1*nr2*nr3)
      integer, pointer :: nim(:) => null()  ! nim(nr1*nr2*nr3)
      !NOTE: idx and nim are used to access image vectors of (i, j, k)
      ! with ir = (k-1)*nr1*nr2 + (j-1)*nr1 + (i-1) + 1
      !
      ! the images of R=(i,j,k) are stored in vec: 
      !    from index idx(ir) to [idx(ir)+nim(ir)-1]
      ! where nim(ir) stores the number of images of R
      ! and nvec = sum(nim)
   end type

   integer, parameter, private :: ws_search_range = 3
   real(dp), parameter, private :: eps6 = 1.0E-6_dp

   public :: init_rvec_images, reset_rvec_images, set_wigner_seitz_cell, &
      get_length, set_cutoff_large, set_cutoff_small
contains

subroutine init_rvec_images(nr1, nr2, nr3, at, rvec_images, large_cutoff)
   implicit none
   integer, intent(in) :: nr1, nr2, nr3
   real(dp), intent(in) :: at(3,3)
   type(vector_set), intent(inout) :: rvec_images
   logical, intent(in), optional :: large_cutoff
   !local
   logical :: lcut
   real(dp) :: cutoff
   integer :: ws_dim(3), vec(3), rvec(3), nim, i, j, k, mi, mj, mk, tot, ir
   
   ws_dim(1) = nr1
   ws_dim(2) = nr2
   ws_dim(3) = nr3
   
   lcut = .false.
   if(present(large_cutoff)) lcut = large_cutoff
   if( lcut ) then
      cutoff = set_cutoff_large(ws_dim, at)
   else
      cutoff = set_cutoff_small(ws_dim, at)
   endif

   if( associated(rvec_images%idx) ) deallocate( rvec_images%idx )
   if( associated(rvec_images%nim) ) deallocate( rvec_images%nim )
   allocate( rvec_images%idx( nr1*nr2*nr3 ), rvec_images%nim( nr1*nr2*nr3 ) )
   
   tot = 0
   ir  = 0
   do k = 1, nr3
   do j = 1, nr2
   do i = 1, nr1
      ir = ir + 1

      vec(1:3) = (/i-1, j-1, k-1/)
      nim = 0
      do mk = -ws_search_range, ws_search_range
      do mj = -ws_search_range, ws_search_range
      do mi = -ws_search_range, ws_search_range
         rvec(1:3) = (/mi, mj, mk/) * ws_dim(1:3) + vec(1:3)
         !
         if( get_length(real(rvec, kind=dp), at) < cutoff ) nim = nim + 1
      enddo; enddo; enddo
      if(nim < 1) call errore('init_wigner_seitz_cell','nim < 1',1)
      
      ! init idx, nim
      rvec_images%idx(ir) = tot + 1
      rvec_images%nim(ir) = nim
      tot = tot + nim
   enddo; enddo; enddo

   rvec_images%nvec = tot
   if( associated(rvec_images%vec) ) deallocate( rvec_images%vec )
   allocate( rvec_images%vec(3, tot) )

   ! collect all the vectors 
   nim = 0
   do k = 1, nr3
   do j = 1, nr2
   do i = 1, nr1
      vec(1:3) = (/i-1, j-1, k-1/)
      
      do mk = -ws_search_range, ws_search_range
      do mj = -ws_search_range, ws_search_range
      do mi = -ws_search_range, ws_search_range
         rvec(1:3) = (/mi, mj, mk/) * ws_dim(1:3) + vec(1:3)
         
         if( get_length(real(rvec, kind=dp), at) < cutoff ) then
            nim = nim + 1
            rvec_images%vec(:,nim) = rvec
         endif
      enddo; enddo; enddo
   enddo; enddo; enddo
end subroutine init_rvec_images


! <0,a| H_[c] | R,b>
subroutine set_wigner_seitz_cell &
      (nr1, nr2, nr3, at, rvec_images, ws, tau_a, tau_b, ref_c)
   implicit none
   integer,  intent(in) :: nr1, nr2, nr3
   real(dp), intent(in) :: at(3,3)
   type(vector_set), intent(in) :: rvec_images
   type(ws_cell), intent(inout) :: ws
   ! tau_a and tau_b are within unit cell and in crystal coordinate.
   real(dp), intent(in) :: tau_a(3), tau_b(3)
   real(dp), intent(in), optional :: ref_c(3) ! in crystal coordinate
   !local
   integer :: nsize, idx, nim, n, m, tot, ir, tot_r
   real(dp) :: vec_b(3)
   integer, allocatable :: itmp(:), ndeg(:)
   real(dp), allocatable :: dist(:)

   nsize = maxval( rvec_images%nim(:) )
   tot_r = nr1 * nr2 * nr3
   !work array
   allocate( itmp( rvec_images%nvec ), ndeg(tot_r) )
   itmp = 0;  ndeg = 0

!$omp parallel default(shared) private(ir, idx, nim, n, vec_b, dist)
   allocate( dist(nsize) )
!$omp do schedule(guided)
   do ir = 1, tot_r
      ! R = (i, j, k) => 
      !    ir = (k-1)*nr1*nr2 + (j-1)*nr1 + (i-1) + 1
      !
      idx = rvec_images%idx(ir)
      nim = rvec_images%nim(ir)
      do n = 1, nim
         vec_b = real(rvec_images%vec(:, n+idx-1), kind=dp) + tau_b
         dist(n) = get_length(vec_b-tau_a, at)
         ! 
         if(present(ref_c)) dist(n) = dist(n) + get_length(vec_b-ref_c, at)
      enddo
      ! find the shortest one(s) or its equivalences
      dist(1:nim) = dist(1:nim) - minval(dist(1:nim))
      ndeg(ir) = count( dist(1:nim) < eps6 )
      
      do n = 1, nim
         if( dist(n) < eps6 ) itmp( n+idx-1 ) = ndeg(ir)
      enddo
   enddo
!$omp end do
   deallocate( dist )
!$omp end parallel
   !sanity check
   if( any( ndeg < 1 ) ) call errore('set_wigner_seitz_cell','ndeg < 1',1)
   
   tot  = sum( ndeg(:) )
   ws%nr = tot
   if( associated( ws%ndeg ) ) deallocate( ws%ndeg )
   if( associated( ws%rvec ) ) deallocate( ws%rvec )
   allocate( ws%rvec(tot), ws%ndeg(tot) )

   !collect vectors into ws%rvec, only store their index in rvec_images%vec
   m = 0
   do n = 1, rvec_images%nvec
      if( itmp(n) > 0 ) then
         m = m + 1
         ws%rvec(m) = n
         ws%ndeg(m) = itmp(n)
      endif
   enddo
   
   deallocate(itmp, ndeg)
end subroutine set_wigner_seitz_cell


subroutine reset_rvec_images(rvec_images)
   implicit none
   type(vector_set), intent(inout) :: rvec_images
   
   rvec_images%nvec = 0
   if( associated(rvec_images%vec) ) deallocate( rvec_images%vec )
   if( associated(rvec_images%idx) ) deallocate( rvec_images%idx )
   if( associated(rvec_images%nim) ) deallocate( rvec_images%nim )
end subroutine reset_rvec_images


real(dp) pure function set_cutoff_small(rdim, at) result(cutoff)
   implicit none
   integer, intent(in) :: rdim(3)
   real(dp), intent(in) :: at(3,3)
   !local variables
   integer :: i, j, k, ndim(3)
   real(dp) :: vec(3), dist
   
   ndim(1) = rdim(1) / 2 + 1
   ndim(2) = rdim(2) / 2 + 1
   ndim(3) = rdim(3) / 2 + 1

   cutoff = 0.0
   do k = -1, 1, 2
   do j = -1, 1, 2
   do i = -1, 1, 2
      vec(1:3) = real( ndim(1:3) * (/i, j, k/), kind=dp )
      ! length of vec: |vec|
      dist = get_length(vec, at)
      ! set cutoff to the longest |vec|
      if(dist > cutoff) cutoff = dist
   enddo; enddo; enddo
end function set_cutoff_small


real(dp) pure function set_cutoff_large(rdim, at) result(cutoff)
   implicit none
   integer, intent(in) :: rdim(3)
   real(dp), intent(in) :: at(3,3)
   !local variables
   integer :: i, j, k
   real(dp) :: vec(3), dist
      
   cutoff = 0.0
   do k = -1, 1, 2
   do j = -1, 1, 2
   do i = -1, 1, 2
      vec(1:3) = real( rdim(1:3) * (/i, j, k/), kind=dp )
      ! length of vec: |vec|
      dist = get_length(vec, at)
      ! set cutoff to the longest |vec|
      if(dist > cutoff) cutoff = dist
   enddo; enddo; enddo
end function set_cutoff_large


real(dp) pure function get_length(vec, at) result(length)
   implicit none
   real(dp), intent(in) :: vec(3) !R vector in crystal coordinate
   real(dp), intent(in) :: at(3,3) ! lattice vector i:  at(:, i)
   !local
   real(dp) :: cvec(3)

   ! fron crystal coord. to cart. coordinate. call cryst_to_cart(1,vec,at,1)
   cvec = matmul(at, vec)
   ! length of vec: |vec|
   length = norm2( cvec )

end function get_length

end module wigner_seitz_cell

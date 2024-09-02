!===============================================================================
! Copyright (C) 2016-2020 Jin-Jian Zhou, Jinsoo Park, I-Te Lu, Marco Bernardi
!
! This program is distributed under the terms of the GNU General Public License.
! See the file `LICENSE' in the root directory of this distribution, or obtain 
! a copy of the License at <https://www.gnu.org/licenses/gpl-3.0.txt>.
!
! Author: jjzhou <jjchou.comphy@gmail.com>
! Comment:
!  setup up the neighboring k and weights for the k-space derivatives 
!  on a regular gamma centered mp grid. using the finite-difference formula 
!  in wannier90 (see kmesh.f90 herein). 
!  Refer to Mostofi et.al comp. phys. commu. 178 (2008) 685-699 for detail.
!
! Maintenance:
!===============================================================================

module boltz_grid_neighbors
   use kinds, only: dp
   use pert_data, only: bg, at, tpiba
   use boltz_grid, only: grid
   use qe_mpi_mod, only: stdout, ionode, mp_split_pools, mp_sum, inter_pool_comm
   use pert_const, only: bohr2ang
   implicit none
   private
   ! data structure to store bvectors in one shell.
   type :: shell
      real(dp) :: dist  ! the radius of the shell
      real(dp) :: weight
      integer  :: numv ! number of vectors in this shell
      real(dp), pointer :: bvec(:,:)=>null() ! bvectors
   end type
   !some predifined parameters
   real(dp), parameter :: b1_tol = 1.0E-5_dp
   !threshold to distinguish shells
   real(dp), parameter :: shell_tol =  1.0E-3_dp
   
   public :: set_grid_neighbors
contains

subroutine set_grid_neighbors(kgrid)
   use pert_utils,  only: binary_search
   use boltz_utils, only: num2kpt, kpt2num
   implicit none
   type(grid), intent(inout) :: kgrid

   integer :: nbvec, ik_start, ik_stop, ik, i, ikb
   real(dp) :: kpt(3), xkb(3)
   real(dp), allocatable :: bvec(:,:)

   call init_grid_neighbors(kgrid)

   nbvec = kgrid%num_neighbors
   allocate(bvec(3, nbvec), kgrid%neighbors(nbvec, kgrid%nk))
   !find and store neighbors of each kpts on the k-grid.
   kgrid%neighbors(:,:) = 0
   bvec(:,:) = kgrid%bvec(:,:)
   !transform from cartesian coord. to crystal coordinate
   call cryst_to_cart(nbvec, bvec, at, -1)

   !mpi + openmp parallelization
   !distribute kpts over all mpi processors.
   call mp_split_pools(kgrid%nk, ik_start, ik_stop)

!$omp parallel do schedule(guided) default(shared) private(ik, kpt, i, xkb, ikb)
   do ik = ik_start, ik_stop
      kpt(1:3) = num2kpt( kgrid%kpt(ik), kgrid%ndim )
      !loop over all neighbors
      do i = 1, nbvec
         xkb(1:3) = kpt(1:3) + bvec(1:3, i)
         ikb = kpt2num( xkb(1:3), kgrid%ndim )
         if( ikb < 0 ) call errore('set_grid_neighbors', 'neighbor not on grid', 1)
         !if not found, neighbors(i, ik) is 0.
         !if found, neighbors(i,ik) stores this neighbor's position in kgrid%kpt
         kgrid%neighbors(i,ik) = binary_search(kgrid%kpt, ikb)
      enddo
   enddo
!$omp end parallel do
   call mp_sum(kgrid%neighbors, inter_pool_comm)

   deallocate(bvec)
end subroutine set_grid_neighbors

subroutine init_grid_neighbors(kgrid)
   implicit none
   type(grid), intent(inout) :: kgrid

   integer, parameter :: maxshell = 6
   ! define area where to search neighboring k
   integer, parameter :: kcutoff = 5, nkpts=(2*kcutoff+1)**3-1

   ! store the neighboring k in shells.
   type(shell) :: bsh(maxshell)
   ! local variables
   logical :: lpar, lreject, lsat
   integer :: indx(nkpts), ik, i, j, k, ish, nvec
   real(dp) :: xk(3), mindist
   real(dp), allocatable :: kpts(:,:), dist(:)

   allocate( kpts(3, nkpts), dist(nkpts) )
   ! find the near neighbors. 
   ! use gamma point as a reference point, setup a sub-grid around gamma, 
   ! then sort the sub-grid based on their distance to gamma.
   ik = 0
   do i = -kcutoff, kcutoff
   do j = -kcutoff, kcutoff
   do k = -kcutoff, kcutoff
      if( i.eq.0 .and. j.eq.0 .and. k.eq.0) cycle
      ik = ik + 1
      xk(1:3) = (/i, j, k/)
      kpts(1:3, ik) = xk(1:3) / kgrid%ndim(1:3)
   enddo; enddo; enddo
   if(ik .ne. nkpts) call errore('init_grid_neighbors','failed: (ik .ne. nkpts)',1)
   !transform to cartesian coordiantes
   call cryst_to_cart(nkpts, kpts, bg, 1)
   !compute distance to gamma
   do ik = 1, nkpts
      dist(ik) = sqrt( dot_product(kpts(:,ik), kpts(:,ik)) )
   end do
   !sort distances in ascend order using shell sort algorithm
   call shellsort(nkpts, dist, indx)
   mindist = dist(1)
   
   !setup up bvectors and weights (kmesh_shell_automatic in wannier90/kmesh.f90)
   !start from the first shell. take the bvectors from the next shell
   !Reject if parallel to existing b vectors, or result in very small singular value
   ish = 1;  bsh(ish)%dist = dist(1);  nvec = 1
   do ik = 2, nkpts
      if( abs( dist(ik) - bsh(ish)%dist ) < shell_tol*mindist ) then
         !belongs to the same shell.
         nvec = nvec + 1
      else
         !end of the current shell
         !check the current shell is not parallel to previous shells
         lpar = .false.
         if(ish > 1) then
            ! check each vector in current shell
            do i = 1, nvec
               lpar = check_parallel(ish-1, bsh(1:ish-1), kpts(:,indx(ik-i)))
               ! if there is vector parallel with prevous shell, skip this shell
               if(lpar) exit
            enddo
         endif
         if(lpar) then
            ! abandom the current shell. redo the current shell with new kpts
            nvec = 1;  bsh(ish)%dist = dist(ik)
         else
            ! save the current shell.
            bsh(ish)%numv = nvec
            allocate( bsh(ish)%bvec(3, bsh(ish)%numv) )
            ! backward copy kpoints that belong to this shell.
            do i = 1, bsh(ish)%numv
               bsh(ish)%bvec(:, i) = kpts(:, indx(ik-i) )
            enddo
            ! compute weight, 
            ! Test to see if we satisfy B1, if not add another shell and repeat.
            call calc_shell_weight(ish, bsh(1:ish), lsat, lreject)
            ! job done. obtain all weight
            if(lsat) exit
            ! needs more shells to satisfy B1
            if(lreject) then
               ! reject the current shell, redo the current shell with new kpts
               deallocate( bsh(ish)%bvec )
               nvec = 1;   bsh(ish)%dist = dist(ik)
            else
               ! does not satisfy B1, add new shell.
               ish = ish + 1
               ! no more than max_shell shells
               if(ish > maxshell) call errore('init_grid_neighbors', &
                  'Need more shells, please increase maxshell (it is hard-coded)!',1)
               ! initialize the new shell.
               nvec = 1;  bsh(ish)%dist = dist(ik)
            endif
         endif
      endif
   enddo
   if(.not. lsat) call errore('init_grid_neighbors', &
      'B1 not satisfied, please increase kcutoff (it is hard-coded)!', 1)

   !save the results to boltz_data
   nvec = 0
   do i = 1, ish
      nvec = nvec + bsh(i)%numv
   enddo
   kgrid%num_neighbors = nvec
   allocate( kgrid%bvec(3, nvec), kgrid%bweight(nvec) )
   ik = 0
   do i = 1, ish
      do j = 1, bsh(i)%numv
         ik = ik + 1
         !cartesian coordinate, in the unit of tpiba/bohr
         kgrid%bvec(:,ik) = bsh(i)%bvec(:,j)
         !bweight is in the unit of (tpiba/bohr)^-2 (bweight*bvec*bvec -> dimensionless)
         !so bweight*bvec is in the unit of (tpiba/bohr)^-1 (df/dk is in (tpiba/bohr)^-1)
         kgrid%bweight(ik) = bsh(i)%weight
      enddo
   enddo

   !! for test, comparison with wannier90 (in the unit of Ang^-1, Ang**2). Passed.
   !if(ionode) then
   ! do i=1, num_bvec
   !  write(stdout, '(3f12.6,5x,f12.6)') bvec(:,i)*tpiba/bohr2ang, bweight(i)*(bohr2ang**2)/(tpiba**2)
   ! enddo
   !endif
   ! clean up
   deallocate(kpts, dist)
end subroutine init_grid_neighbors

function check_parallel(nshell, ash, kvec)
   implicit none
   integer, intent(in) :: nshell
   type(shell), intent(in) :: ash(nshell)
   real(dp), intent(in) :: kvec(3)
   logical :: check_parallel
   ! local variables
   integer :: is, j
   real(dp) :: bvec(3), delta

   check_parallel = .false.
   do is = 1, nshell
   do j = 1, ash(is)%numv
      bvec(:) = ash(is)%bvec(1:3,j)
      delta = dot_product(kvec, bvec) / &
         sqrt( dot_product(kvec,kvec) * dot_product(bvec,bvec) )
      if(abs( abs(delta)-1.0E0_dp ) < 1.0E-5_dp) then
         check_parallel = .true.
         return
      endif
   enddo; enddo
end function check_parallel

subroutine calc_shell_weight(nshell, ash, lsatb1, lrej)
   implicit none
   integer, intent(in) :: nshell
   type(shell), intent(inout) :: ash(nshell)
   logical, intent(out) :: lsatb1, lrej
   !local variables
   integer :: is, i, j, k, info
   integer, parameter :: lwork = 60
   real(dp) :: work(lwork), bwt(nshell), delta
   real(dp), allocatable :: amat(:,:), umat(:,:), vmat(:,:), smat(:,:), singv(:)
   ! delta_alpha_beta in vector form
   real(dp), parameter :: qtarget(6)=(/1.0_dp,1.0_dp,1.0_dp,0.0_dp,0.0_dp,0.0_dp/)

   ! allocate space for the current number of shells.
   allocate( amat(6, nshell), umat(6, 6), vmat(nshell, nshell) )
   allocate( smat(nshell, 6), singv(nshell) )
   amat=0.0_dp; umat=0.0_dp; vmat=0.0_dp; smat=0.0_dp; singv=0.0_dp
   ! initialize amat
   do is = 1, nshell
      do i = 1, ash(is)%numv
         amat(1,is) = amat(1,is) + ash(is)%bvec(1,i) * ash(is)%bvec(1,i)
         amat(2,is) = amat(2,is) + ash(is)%bvec(2,i) * ash(is)%bvec(2,i)
         amat(3,is) = amat(3,is) + ash(is)%bvec(3,i) * ash(is)%bvec(3,i)
         amat(4,is) = amat(4,is) + ash(is)%bvec(1,i) * ash(is)%bvec(2,i)
         amat(5,is) = amat(5,is) + ash(is)%bvec(2,i) * ash(is)%bvec(3,i)
         amat(6,is) = amat(6,is) + ash(is)%bvec(3,i) * ash(is)%bvec(1,i)
      enddo
   enddo
   ! computing weight using singular value decomposition, see Sec.3.2 in Ref.
   info = 0
   call dgesvd('A','A',6,nshell,amat,6,singv,umat,6,vmat,nshell,work,lwork,info)
   if(info < 0) then
      write(stdout,'(1x,a,1x,I1,1x,a)') &
         'boltz_neighbors_setup: Argument',abs(info),'of dgesvd is incorrect'
      call errore('boltz_neighbors_setup','Problem with Singular Value Decomposition',1)
   else if (info > 0) then
      call errore('boltz_neighbors_setup','Singular Value Decomposition did not converge',1)
   endif
   
   lrej = .false.
   if(any( abs(singv) < 1.0E-5_dp )) then
      if(nshell .eq. 1) then
         call errore('boltz_neighbors_setup',&
            'Singular Value Decomposition has found a very small singular value',1)
      else
         write(stdout,'(1x,a)') &
            '| SVD found small singular value, Rejecting this shell and trying the next   |'
         lrej = .true.
         return
      end if
   end if
   ! D^-1 in eq.27 of Reference.
   do is = 1, nshell
      smat(is, is) = 1.0E0_dp / singv(is)
   enddo
   ! compute the weight
   bwt(1:nshell)=matmul(transpose(vmat),matmul(smat,matmul(transpose(umat),qtarget)))
   ! check if the weight satisfy the B1 equation.
   lsatb1 = .true.
   do i = 1, 3
   do j = i, 3
      delta = 0.0E0_dp
      do is = 1, nshell
      do k = 1, ash(is)%numv
         delta = delta + bwt(is) * ash(is)%bvec(i,k) * ash(is)%bvec(j,k)
      enddo; enddo
      if(i.eq.j .and. abs(delta-1.0E0_dp) > b1_tol ) lsatb1 = .false.
      if(i.ne.j .and. abs(delta) > b1_tol ) lsatb1 = .false.
   enddo; enddo
   if(lsatb1) then
      do is = 1, nshell
         ash(is)%weight = bwt(is)
      enddo
   endif
   ! clean up
   deallocate(amat, umat, vmat, smat, singv)
end subroutine calc_shell_weight

end module boltz_grid_neighbors


subroutine shellsort(num, a, ind)
   use kinds, only: dp
   implicit none
   integer, intent(in) :: num
   real(dp), intent(inout) :: a(num)
   integer, intent(out) :: ind(num)
   ! local variables
   integer  :: i, j, increment, itmp
   real(dp) :: temp
   ! initialize the index array.
   do i = 1, num
      ind(i) = i
   enddo
   ! sort a using shellsort algorithm, and move the ind accordingly.
   ! note on shellsort: if increment = 1, it's basically insertion sort.
   increment = num / 2
   do while (increment > 0)
      do i = increment+1, num
         j = i
         temp = a(i);  itmp = ind(i)
         do while (j >= increment + 1) 
            if( a(j-increment) > temp ) then
               a(j)   = a(j-increment)
               ind(j) = ind(j-increment)
               j = j - increment
            else
               exit
            endif
         enddo
         a(j) = temp; ind(j) = itmp
      enddo
      if(increment == 2) then
         increment = 1
      else
         increment = increment * 5 / 11
      endif
   enddo
end subroutine shellsort

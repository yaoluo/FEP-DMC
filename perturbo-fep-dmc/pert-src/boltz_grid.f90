!===============================================================================
! Copyright (C) 2016-2020 Jin-Jian Zhou, Jinsoo Park, I-Te Lu, Marco Bernardi
!
! This program is distributed under the terms of the GNU General Public License.
! See the file `LICENSE' in the root directory of this distribution, or obtain 
! a copy of the License at <https://www.gnu.org/licenses/gpl-3.0.txt>.
!
! Author: jjzhou <jjchou.comphy@gmail.com>
! Comment:
!   define data used in boltzmann dynamics and transport calculations
!
! Maintenance:
!===============================================================================

module boltz_grid
   use pert_const, only: dp, ryd2mev
   implicit none
   private
 
   type, public :: grid
      !k-point
      integer :: ndim(3) ! size of the grid
      integer :: nk  !number of selected kpoints, includes reducible kpoints.
      integer :: nk_irr !number of irreducible k-points
      !tetrahedra
      !number of tetrahedra, and weight of each tetrahedra
      integer  :: num_tet
      real(dp) :: tweight
      ! tet(2,4,num_tet): tet(1,:,:)->irreducible k; tet(2,:,:)->full grid
      integer, pointer :: tetra(:,:,:)=>null()
      !
      real(dp):: kweight !weight of each k-point, = 1/(kdim(1)*kdim(2)*kdim(3))
      !
      !index representation of kpoints in full grid
      integer, pointer :: kpt(:)=>null() ! kpt(nk)
      !irreducible k-list, irrk(nk_irr)
      ! store index representation of this kpoints
      integer, pointer :: irrk(:)=>null()  !irrk(nk_irr)
      !map the i-th kpt in full grid to its irreducible kpt 
      ! irrk( kpt2ir(i) ): i-th kpoints in full grid -> irr-kpoints.
      integer, pointer :: kpt2ir(:)=>null() ! kpt2ir(nk)
      !
      !band structure and band velocity information
      ! kgrid_enk(numb, nk)
      real(dp), pointer :: enk(:,:)=>null()
      ! E_nk for irreducible k-points in irrk.
      real(dp), pointer :: enk_irr(:,:)=>null()
      ! band velocity, vnk(3,numb,nk)
      real(dp), pointer :: vnk(:,:,:)=>null()
      ! only bands in [bmin, bmax] contribute to transport
      ! 1 < bmin < bmax < num_bands;  numb = bmax - bmin + 1
      integer :: numb, bmin, bmax

      !neighbor k-points, for kspace derivatives
      !using the finite-difference formula in wannier90,
      !Refer to Mostofi et.al comp. phys. commu. 178 (2008) 685-699
      integer :: num_neighbors   ! number of neighbor
      !b vectors in cartesian coordinates, unit tpiba/bohr.
      real(dp), pointer :: bvec(:,:)=>null()
      !weight used to compute df/dk, unit (tpiba/bohr)^-2
      real(dp), pointer :: bweight(:)=>null()
      !neighbors(i,ik): index of i-th neighbor of the ik-points.
      integer,  pointer :: neighbors(:,:)=>null()
   end type grid

   !real(dp), parameter :: eig_tol = 20.0E0_dp/ryd2mev ! 50 meV

   public :: init_boltz_grid
   public :: init_boltz_grid_sym
   public :: boltz_grid_generate
   public :: boltz_grid_load
   public :: sum_equivalent_kpoint
contains

subroutine init_boltz_grid(kg, el, bmin, bmax, load)
   use boltz_utils, only: num2kpt
   use qe_mpi_mod, only: mp_split_pools, mp_sum, inter_pool_comm
   use band_structure, only: electron_wann, solve_band_velocity
   implicit none
   type(electron_wann), intent(in) :: el
   integer, intent(in) :: bmin, bmax
   type(grid), intent(inout) :: kg !kgrid
   !if false, then kg is already initialized, no need to load from files.
   logical, intent(in), optional :: load
   !local variable
   logical :: load_grid
   integer :: ik, kst, kend
   real(dp) :: xk(3), eig(el%nb), vel(3,el%nb)

   !initialize numb
   kg%numb = bmax - bmin + 1
   kg%bmin = bmin
   kg%bmax = bmax
   !
   !read input files, init nk, nk_irr, kweight, irrk_k,
   load_grid = .true.
   if( present(load) ) load_grid = load
   if( load_grid ) call boltz_grid_load(kg)
   !
   !calc. eig and velocity
   allocate(kg%enk(kg%numb, kg%nk), kg%vnk(3, kg%numb, kg%nk))

   kg%enk = 0.0_dp;  kg%vnk = 0.0_dp
   call mp_split_pools(kg%nk, kst, kend)
!$omp parallel do schedule(guided) default(shared) private(ik, xk, eig, vel)
   do ik = kst, kend
      xk = num2kpt(kg%kpt(ik), kg%ndim)
      call solve_band_velocity(el, xk, vel, eigvalue=eig)
      kg%enk(:,ik) = eig(bmin:bmax)
      kg%vnk(:,:,ik) = vel(:, bmin:bmax)
   enddo
!$omp end parallel do
   call mp_sum(kg%enk, inter_pool_comm)
   call mp_sum(kg%vnk, inter_pool_comm)
end subroutine init_boltz_grid

!compute enk and vel for irreducible kpoint only, 
!those of reducible k are determined from symmetry.
subroutine init_boltz_grid_sym(kg, el, bmin, bmax, load)
   use pert_data, only: nsym, symop, at, bg
   use boltz_utils, only: num2kpt, kpt2num
   use qe_mpi_mod, only: mp_split_pools, mp_sum, inter_pool_comm
   use band_structure, only: electron_wann, solve_band_velocity
   implicit none
   type(electron_wann), intent(in) :: el
   !nsym: number of symmetry operation, symop: symmetry operation in crystal coord.
   integer, intent(in) :: bmin, bmax
   type(grid), intent(inout) :: kg !kgrid
   !if false, then kg is already initialized, no need to load from files.
   logical, intent(in), optional :: load
   !local variable
   logical :: load_grid
   integer :: ik, kst, kend, i, j, irk, ib
   real(dp) :: xk(3), eig(el%nb), vel(3,el%nb), symop_inv_cart(3,3,nsym), xkr(3), vabs2
   real(dp), allocatable :: vnk_irr(:,:,:)
   ! Identical matrix
   integer, parameter :: eye(3,3) = reshape((/1, 0, 0, 0, 1, 0, 0, 0, 1/), (/3,3/))

   !initialize numb
   kg%numb = bmax - bmin + 1
   kg%bmin = bmin
   kg%bmax = bmax
   
   !read input files, init nk, nk_irr, kweight, irrk_k, 
   load_grid = .true.
   if( present(load) ) load_grid = load
   if( load_grid ) call boltz_grid_load(kg)
   
   !calc. eig and velocity
   allocate(kg%enk(kg%numb,kg%nk), kg%vnk(3,kg%numb, kg%nk), kg%enk_irr(kg%numb,kg%nk_irr))

   allocate( vnk_irr(3, kg%numb, kg%nk_irr) )
   kg%enk_irr = 0.0_dp;  vnk_irr = 0.0_dp
   
   call mp_split_pools(kg%nk_irr, kst, kend)
!$omp parallel do schedule(guided) default(shared) private(ik, xk, eig, vel)
   do ik = kst, kend
      xk = num2kpt( kg%irrk(ik), kg%ndim )
      call solve_band_velocity(el, xk, vel, eigvalue=eig)
      kg%enk_irr(:,ik) = eig(bmin:bmax)
      vnk_irr(:,:,ik) = vel(:, bmin:bmax)
   enddo
!$omp end parallel do
   call mp_sum(kg%enk_irr, inter_pool_comm)
   call mp_sum(vnk_irr, inter_pool_comm)

   !compute symop in cartesian coord.
   ! k_cart = bg * k, so k = (bg)^-1 * k_cart; and (bg)^-1 = transpose(at)
   ! if S*k = k'; then S * [(bg)^-1 * k_cart] = [(bg)^-1 * k'_cart]
   ! therefore, [bg * S * (bg)^-1] * k_cart = k'_cart
   symop_inv_cart = 0.0_dp
   do i = 1, nsym
      do j = 1, nsym
         if( all(matmul(symop(:,:,i), symop(:,:,j)) .eq. eye) ) then
            symop_inv_cart(:,:,i) = matmul( matmul(bg, symop(:,:,j)), transpose(at) )
            exit
         endif
      enddo
      !throw an error if not initialized
      if(sum( abs(symop_inv_cart(:,:,i)) ) < 1.0E-12_dp) &
         call errore('init_boltz_grid_sym','failed to compute symop_inv_cart',1)
   enddo

   kg%enk = 0.0_dp;  kg%vnk = 0.0_dp
   do ik = 1, kg%nk
      irk = kg%kpt2ir(ik)
      ! E(Sk) = E(k) 
      kg%enk(:, ik) = kg%enk_irr(:, irk)
      !find the symmetry operation that transform irk to ik
      xkr = num2kpt( kg%irrk(irk), kg%ndim )
      do i = 1, nsym
         xk = matmul(symop(1:3, 1:3, i), xkr)
         !if find the operation, compute v(Sk)
         ! v(Sk) = v(k) * S^-1:  v(Sk)_j = \sum_{i} v(k)_i * (S^-1)_ij
         if( kpt2num(xk, kg%ndim) .eq. kg%kpt(ik) ) then
            do ib = 1, kg%numb
               kg%vnk(:, ib, ik) = matmul(vnk_irr(:, ib, irk), symop_inv_cart(:,:,i))
            enddo
            exit
         endif
      enddo
      !check
      vabs2 = dot_product(vnk_irr(:,1,irk), vnk_irr(:,1,irk)) 
      if(abs(dot_product(kg%vnk(:,1,ik), kg%vnk(:,1,ik)) - vabs2) > 1.E-6_dp*vabs2) &
         call errore('init_boltz_grid_sym','failed to unfold velocity', 1)
   enddo
   
   deallocate(vnk_irr)
end subroutine init_boltz_grid_sym


!compute sum_{k} v(k)*v(k) for k corresponding to the same irreducible kpoints
subroutine sum_equivalent_kpoint(kg, vvprod, ndeg)
   use boltz_utils, only: velprod
   implicit none
   type(grid), intent(in) :: kg
   !vvprod(6, numb, kg%nk_irr), ndeg(kg%nk_irr): number of k-points
   real(dp), intent(out) :: vvprod(:,:,:), ndeg(:) 
   !
   integer :: ik, irk, ib, nd(3)
   real(dp) :: vel(3)
   
   nd = shape(vvprod)
   if(nd(1).ne.6 .or. nd(2).ne.kg%numb .or. nd(3).ne.kg%nk_irr .or. size(ndeg).ne.nd(3)) &
      call errore('sum_equivalent_kpoint','inconsistent demension of array vvprod',1)

   vvprod = 0.0_dp
   ndeg   = 0.0_dp
   do ik = 1, kg%nk
      irk = kg%kpt2ir(ik)
      ndeg(irk) = ndeg(irk) + 1.0_dp
      do ib = 1, kg%numb
         vel = kg%vnk(1:3, ib, ik)
         vvprod(1:6, ib, irk) = vvprod(1:6, ib, irk) + velprod(vel, vel)
      enddo
   enddo
end subroutine sum_equivalent_kpoint

#include "boltz_grid_generate.f90"

#include "boltz_grid_load.f90"

end module boltz_grid

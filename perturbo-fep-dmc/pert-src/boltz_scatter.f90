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

module boltz_scatter
   use boltz_grid, only: grid
   use pert_const, only: dp, czero
   use pert_data,  only: epwan_fid, qc_dim, kc_dim
   use pert_param, only: phfreq_cutoff, delta_smear, use_mem, load_scatter_eph
   use band_structure, only: electron_wann, solve_eigenvalue_vector
   use phonon_dispersion, only: lattice_ifc, solve_phonon_modes
   use qe_mpi_mod, only: stdout, mp_barrier, inter_pool_comm, npool, my_pool_id, mp_sum
   use elphon_coupling_matrix, only: elph_mat_wann, release_elph_mat_wann, &
      init_elph_mat_wann, eph_fourier_el_para, eph_fourier_elph, eph_transform_fast
   implicit none
   public
   !user-defined data strucutre for scatering channel. based on kgrid
   type :: t_scatt
      !location of k and kq in kgrid%kpt: xk(1:3) = num2kpt(kgrid%kpt(ik), ndim)
      !NOTE: (k, kq) and (kq, k) are equivalent, only one of them are stored.
      integer :: ik
      integer :: ikq
      !location of q = (k+q) - k in qgrid%list: xq(1:3) = num2kpt(qgrid%list(iq), ndim)
      integer :: iq
      !number of (mkq, nk, mu) terms, mu is phonon modes index
      !only g_m,n,mu of (k, kq) satisfy the following condition are stored. 
      !abs( abs(e_nk - e_mkq) - w_qmu ) < delta_E. (delta_E = e_thr*broadening)
      integer :: nchl
      
      !use logical table to determine nt and select (mkq, nk, mu) terms:
      ! if table_bands(mkq, nk, mu) is true, then
      !   matrix elements:  |g_{ mu, iq}(nk,ik; mkq,ikq)|^2 is stored
      !use any() or count() to check if there is .true. in table_bands and nt
      !
      !index and g2 are stored for selected (mkq, nk, mu). 
      ! map (mkq, nk, mu) to a integer 'index' and store it in bands_index
      ! (mu-1)*numb*numb + (nk-1)*numb + (mkq-1)  + 1 = index
      integer, pointer :: bands_index(:)=>null()
      !eph_g2(nchl), nchl = size( bands_index ).
      real(dp), pointer :: eph_g2(:)=>null()
      !eph_g2(it) : |g_{ mu, iq}(nk,ik; mkq,ikq)|^2
   end type t_scatt

   type, private :: channel
      integer :: ik
      integer :: nkq
      ! 
      integer, pointer :: iq(:) => null()
      integer, pointer :: ikq(:) => null()
      integer, pointer :: nchl(:) => null()
   end type

   type :: phonon_grid
      integer :: ndim(3)  !size of the grid
      !
      integer :: npts  !number of points selected
      integer, pointer :: list(:) => null()
      real(dp), pointer :: freq(:,:) => null()
      !
      real(dp):: weight  ! 1/(ndim(1) * ndim(2) * ndim(3))
   end type

   integer, private, save :: nk_loc
   type(channel), allocatable, private, save :: kq_pair(:)
   !
   type(phonon_grid), save :: qgrid
   !
   integer, save :: num_scat
   type(t_scatt), allocatable, save :: scat(:)
contains

subroutine boltz_scatter_setup(kg, el, ph, qdim)
   implicit none
   type(grid), intent(in) :: kg
   type(electron_wann), intent(in) :: el
   type(lattice_ifc), intent(in) :: ph
   !dimension for the q-grid. q-should be commensurate with the k-grid
   integer, intent(in) :: qdim(3)
   !local 
   integer :: nkpt, i, ik, collect_pools(npool), nq_col(npool)
   real(dp) :: wcut, e_thr
   type(elph_mat_wann) :: elph
   
   call start_clock('scat_setup')
   write(stdout,'(5x,a20)') '>start scatter_setup'
   !init
   nkpt = kg%nk
   wcut  = phfreq_cutoff
   e_thr = delta_smear*3.0_dp  !exp(-3*3) ~ 10^-4
   
   !interleave distribution of kpoints over pools, for better load balance
   nk_loc = nkpt/npool + merge(1, 0, my_pool_id < mod(nkpt, npool))
   allocate( kq_pair(nk_loc) )
   !store the allocated kpts (location in kg%kpt) of the current pool
   i = 0
   do ik = my_pool_id+1, nkpt, npool
      i = i + 1
      kq_pair(i)%ik = ik
   enddo
   !debug
   if(i .ne. nk_loc) call errore('boltz_scatter_setup', 'i .ne. nk_loc', 1)
   !
   qgrid%ndim(1:3) = qdim(1:3)
   qgrid%weight = 1.0_dp / dble(qdim(1)) / dble(qdim(2)) / dble(qdim(3))
   !generate (k,q) pair (initialized kq_pair, qgrid)
   call setup_kq_pair(kg, el, ph, wcut, e_thr, nk_loc, kq_pair, qgrid)
   !
   num_scat = 0
   do i = 1, nk_loc
      num_scat = num_scat + kq_pair(i)%nkq
   enddo
   !
   nq_col = 0
   nq_col( my_pool_id+1 ) = qgrid%npts
   call mp_sum( nq_col, inter_pool_comm )
   !
   collect_pools = 0
   collect_pools( my_pool_id+1 ) = num_scat
   call mp_sum( collect_pools, inter_pool_comm )
   !output info
   write(stdout, '(5x, a)') "No. of q-points and (k,q) pairs on each processor:"
   do i = 1, npool
      write(stdout, '(5x, a, i4, a, i9, a, i12)')  &
         "Proc.", i, ':  #.q ', nq_col(i), ' |  #.(k,q)', collect_pools(i)
   enddo
   !
   ! compute and save g2 to disk
   if(.not. load_scatter_eph) then
      call init_elph_mat_wann(epwan_fid, kc_dim, qc_dim, elph)
      !compute g2 and sove them to temperary files
      call compute_scatter_eph_g2 &
         (kg, el, ph, elph, qgrid, wcut, e_thr, use_mem)
      !release space
      call release_elph_mat_wann(elph)
   endif
   !
   allocate( scat(num_scat) )
   !read g2 from files
   call load_scatter_eph_g2(nk_loc, kq_pair, num_scat, scat)
   deallocate( kq_pair )
   !
   write(stdout,'(5x,a20)') '>scatter_setup done.'
   call stop_clock('scat_setup')
end subroutine boltz_scatter_setup

#include "setup_kq_pair.f90"

#include "compute_scatter_eph_g2.f90"

#include "load_scatter_eph_g2.f90"

end module boltz_scatter

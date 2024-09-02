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

subroutine calc_band_structure()
   use pert_const, only: dp, ryd2ev
   use pert_data,  only: bg, epwan_fid, qc_dim, kc_dim
   use pert_param, only: prefix, fklist, read_H
   use vector_list,only: vlist, load_vector_list
   use qe_mpi_mod, only: ionode, mp_split_pools, mp_sum, inter_pool_comm, npool
   use band_structure, only: electron_wann, init_electron_wann, solve_eigenvalue_vector
   use HamiltonianTable, only : Tabulate_for_DMC, solve_electron_fast_double
   use elphon_coupling_matrix, only: elph_mat_wann, init_elph_mat_wann
   use phonon_dispersion, only: lattice_ifc, init_lattice_ifc
   implicit none
   type(vlist) :: kl
   integer :: ik, kst, kend
   character(len=80) :: fname
   real(dp), allocatable :: enk(:,:), xloc(:)
   !

   type(lattice_ifc)   :: phon
   type(electron_wann) :: elec
   type(elph_mat_wann) :: elph

   call init_electron_wann(epwan_fid, kc_dim, elec)
   call init_lattice_ifc(epwan_fid, qc_dim, phon)
   call init_elph_mat_wann(epwan_fid, kc_dim, qc_dim, elph)

   !
   call load_vector_list(fklist, kl)
   call mp_split_pools(kl%nvec, kst, kend)

   allocate(enk(elec%nb, kl%nvec))
   enk = 0.0_dp
   if(npool > kl%nvec) &
      call errore('calc_bandstructure','too many pools (npool > #.k-points)',1)
   
   if(read_H .eq. .true.) then 
      call Tabulate_for_DMC(elec, phon, elph)
      !call Tabulate_frohlich_only(elec, phon, elph)
      do ik = kst, kend
         call solve_electron_fast_double(kl%vec(:,ik), enk(:,ik))
      enddo
      enk = enk / ryd2ev
      
   else 
      !$omp parallel do schedule(guided) default(shared) private(ik)
      do ik = kst, kend
         call solve_eigenvalue_vector(elec, kl%vec(:,ik), enk(:,ik))
      enddo
      !$omp end parallel do
   endif  
   call mp_sum(enk, inter_pool_comm)
   enk = enk * ryd2ev

   allocate(xloc(kl%nvec))
   call generate_path(kl, bg, xloc)
   fname=trim(prefix)//'.bands'
   if(ionode) call output_bands(fname, elec%nb, kl%nvec, xloc, kl%vec, enk)
end subroutine calc_band_structure

subroutine calc_phonon_spectra()
   use pert_const, only: dp, ryd2mev, ryd2ev
   use pert_data,  only: bg, epwan_fid, qc_dim, kc_dim
   use pert_param, only: fqlist, prefix, read_H
   use vector_list,only: vlist, load_vector_list
   use qe_mpi_mod, only: ionode, mp_split_pools, mp_sum, inter_pool_comm, npool
   use phonon_dispersion, only: lattice_ifc, init_lattice_ifc, solve_phonon_modes
   use band_structure, only: electron_wann, init_electron_wann, solve_eigenvalue_vector
   use HamiltonianTable, only : tabulate_for_DMC, solve_phonon_fast_double 
   use elphon_coupling_matrix, only: elph_mat_wann, init_elph_mat_wann
   !use lattice_dynamical_tdep, only: solve_phonon_modes_tdep
   implicit none
   type(vlist) :: ql
   integer :: qst, qend, iq
   character(len=90) :: fname
   real(dp), allocatable :: wq(:,:), xloc(:)
   !
   type(lattice_ifc)   :: phon
   type(electron_wann) :: elec
   type(elph_mat_wann) :: elph

   call init_electron_wann(epwan_fid, kc_dim, elec)
   call init_lattice_ifc(epwan_fid, qc_dim, phon)
   call init_elph_mat_wann(epwan_fid, kc_dim, qc_dim, elph)


   call load_vector_list(fqlist, ql)
   call mp_split_pools(ql%nvec, qst, qend)
   if(npool > ql%nvec) &
      call errore('calc_band_structure','too many pools (npool > #.q-points)',1)
    
   allocate(wq(3*phon%na, ql%nvec))
   wq = 0.0_dp

   if(read_H) then 
      call tabulate_for_DMC(elec, phon, elph)
      do iq = qst, qend
         call solve_phonon_fast_double(ql%vec(:,iq), wq(:,iq))
      enddo
      wq = wq / ryd2ev
   else 
      !$omp parallel do schedule(guided) default(shared) private(iq)
      do iq = qst, qend
         call solve_phonon_modes(phon, ql%vec(:,iq), wq(:,iq))
      enddo
      !$omp end parallel do 
   endif 

   call mp_sum(wq, inter_pool_comm)
   wq = wq * ryd2mev

   allocate(xloc(ql%nvec))
   call generate_path(ql, bg, xloc)
   fname=trim(prefix)//'.phdisp'
   if(ionode) call output_bands(fname, 3*phon%na, ql%nvec, xloc, ql%vec, wq)
end subroutine calc_phonon_spectra

subroutine generate_path(list, tmat, path)
   use pert_const, only: dp
   use vector_list, only: vlist
   implicit none
   type(vlist), intent(in) :: list
   real(dp), intent(in) :: tmat(3,3) !transformation matrix
   real(dp), intent(out) :: path(list%nvec)
   integer :: nvec, iv
   real(dp) :: dvec(3)
   real(dp), allocatable :: xvec(:,:)

   nvec = list%nvec
   allocate(xvec(3, nvec))

   xvec = list%vec
   call cryst_to_cart(nvec, xvec, tmat, 1) !kpt in cart. coor.
   path(1) = 0.0_dp
   do iv = 2, nvec
      dvec = xvec(:,iv) - xvec(:,iv-1)
      path(iv) = path(iv-1) + sqrt(dot_product(dvec, dvec))
   enddo
end subroutine generate_path

!output band structure or phonon dispersion
subroutine output_bands(fname, nb, nvec, xloc, vectors, bands)
   use pert_const, only: dp
   use pert_utils, only: find_free_unit
   implicit none
   character(len=80), intent(in) :: fname
   integer, intent(in) :: nb, nvec
   real(dp), intent(in) :: xloc(nvec), vectors(3,nvec), bands(nb,nvec)
   integer :: iunit, iv, ib

   iunit = find_free_unit()
   open(iunit, file=trim(fname), form='formatted',status='unknown')
   do ib = 1, nb
      do iv = 1, nvec
         write(iunit,'(1x, f12.7, 2x, 3(1x,f10.5),2x,f16.10)') &
            xloc(iv), vectors(:,iv), bands(ib,iv)
      enddo
      write(iunit,'(a)') " "
   enddo
end subroutine output_bands

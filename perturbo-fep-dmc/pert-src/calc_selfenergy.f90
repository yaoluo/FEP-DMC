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

subroutine calc_selfenergy()
   use pert_const, only: dp
   use qe_mpi_mod, only: ionode
   use pert_data,  only: epwan_fid, kc_dim, qc_dim
   use pert_param, only: fklist, fqlist, band_max, band_min, cum_inner_emin, &
      cum_inner_emax, cum_outer_emin, cum_outer_emax, cum_de, cum_outer_np, &
      ntemper, temper, efermi, prefix, cauchy_scale, sampling, nsamples
   use vector_list, only: load_vector_list, vlist, rand_vector_list
   use cumulant_utility, only: cum_wgrid, init_cum_wgrid, save_selfenergy
   !
   use band_structure, only: electron_wann, init_electron_wann
   use phonon_dispersion, only: lattice_ifc, init_lattice_ifc
   use elphon_coupling_matrix, only: elph_mat_wann, init_elph_mat_wann
   !
   use selfenergy, only: selfenergy_fan
   implicit none
   type(vlist) :: kl, ql
   type(cum_wgrid) :: wg
   integer :: numb, nestep
   real(dp), allocatable :: enk(:,:), se_imag(:,:,:,:)
   !
   type(lattice_ifc)   :: phon
   type(electron_wann) :: elec
   type(elph_mat_wann) :: elph
   !init
   call init_lattice_ifc(epwan_fid, qc_dim, phon)
   call init_electron_wann(epwan_fid, kc_dim, elec)
   call init_elph_mat_wann(epwan_fid, kc_dim, qc_dim, elph)

   if(ionode)write(*,'(3E15.5)') cum_inner_emin,cum_inner_emax,cum_de
   call init_cum_wgrid(cum_outer_emin, cum_outer_emax, cum_inner_emin, &
      cum_inner_emax, cum_de, cum_outer_np, wg)
   nestep = size(wg%w)
   numb = band_max - band_min + 1

   call load_vector_list(fklist, kl)
   if( trim(fqlist) .ne. '' ) then
      call load_vector_list(fqlist, ql)
   else
      call rand_vector_list(ql, nsamples, sampling, cauchy_scale)
   endif
   allocate( se_imag(nestep, ntemper, numb, kl%nvec), enk(numb, kl%nvec) )

   call selfenergy_fan(elph, elec, phon, kl%vec, wg%w, ql%vec, ql%weight,&
      temper, efermi, band_min, enk, se_imag)
   
   !write data to hdf5 files
   if(ionode) then
      call save_selfenergy(trim(prefix)//'_selfenergy.h5', enk, se_imag, wg%w)
   endif

   deallocate(enk, se_imag)
   return
end subroutine calc_selfenergy

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

subroutine electron_imsigma()
   use pert_const, only: dp
   use qe_mpi_mod, only: ionode, stdout
   use pert_output, only: output_imsigma
   use pert_data,  only:  epwan_fid, kc_dim, qc_dim
   use vector_list, only: vlist, load_vector_list, rand_vector_list
   use pert_param, only: fklist, fqlist, ntemper, temper, efermi, band_min, &
      band_max, polar_split, cauchy_scale, sampling, nsamples
   use imaginary_selfenergy, only: calc_selfel_polar, calc_selfel
   !
   use band_structure, only: electron_wann, init_electron_wann
   use phonon_dispersion, only: lattice_ifc, init_lattice_ifc
   use elphon_coupling_matrix, only: elph_mat_wann, init_elph_mat_wann
   implicit none
   integer :: numb, ntpr
   !load k- and q- list
   type(vlist) :: kl, ql
   real(dp), allocatable:: enk(:,:), imsigma(:,:,:,:)
   !
   type(lattice_ifc)   :: phon
   type(electron_wann) :: elec
   type(elph_mat_wann) :: elph

   !check if temperatures have been readin 
   if(ntemper .eq. 0) call errore('electron_imsigma',&
   'ftemper is not specified, missing temperatures and chemical potential',1)
   ntpr = size(temper)
   numb = band_max - band_min + 1
   !load k- and q-grid
   call load_vector_list(fklist, kl)
   if( trim(fqlist) .ne. '' ) then
      call load_vector_list(fqlist, ql)
   else
      call rand_vector_list(ql, nsamples, sampling, cauchy_scale)
   endif
   if(ionode) then
      write(stdout,'(5x,a,i11)') 'total number of k-points:', kl%nvec
      write(stdout,'(5x,a,i11)') 'total number of q-points:', ql%nvec
   endif
   !init
   call init_electron_wann(epwan_fid, kc_dim, elec)
   call init_lattice_ifc(epwan_fid, qc_dim, phon)
   call init_elph_mat_wann(epwan_fid, kc_dim, qc_dim, elph)
   
   allocate(imsigma(numb, phon%nm, ntpr, kl%nvec), enk(numb, kl%nvec))
   
   select case ( polar_split )
      case('polar')
         call calc_selfel_polar(elph, elec, phon, kl%vec, ql%vec, ql%weight, &
            temper, efermi, band_min, enk, imsigma)
      case('rmpol')
         call calc_selfel(elph, elec, phon, kl%vec, ql%vec, ql%weight, &
            temper, efermi, band_min, enk, imsigma, rm_pol=.true.)
      case default
         call calc_selfel(elph, elec, phon, kl%vec, ql%vec, ql%weight, &
            temper, efermi, band_min, enk, imsigma)
   end select
   !output results
   if(ionode) then
      call output_imsigma(.true. , temper, efermi, enk, imsigma)
      call output_imsigma(.false., temper, efermi, enk, imsigma)
   endif
   !release memory
   deallocate(enk, imsigma)
   
   if(ionode) write(stdout,'(5x, a)') '>Im[Sigma] done.'
end subroutine electron_imsigma

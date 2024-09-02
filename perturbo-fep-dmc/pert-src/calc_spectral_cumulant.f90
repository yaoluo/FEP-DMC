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

subroutine calc_spectral_cumulant()
   use hdf5_utils
   use pert_const, only: dp, pi, ryd2ev, kelvin2eV
   use qe_mpi_mod, only: ionode, stdout, mp_split_pools, inter_pool_comm, mp_sum
   use pert_param, only: ntemper, band_min, band_max, prefix, cum_inner_emin, &
      cum_inner_emax, cum_outer_emin, cum_outer_emax, spectral_emin, cum_de, &
      spectral_emax, cum_outer_np, spectral_np, fklist, temper, efermi
   use vector_list, only: load_vector_list, vlist
   use pert_output, only: progressbar_init, set_progress_step, load_imsigma
   use cumulant_utility, only: cum_wgrid, init_cum_wgrid, load_selfenergy, compute_beta
   use cumulant_expansion,only: cumulant_spectral, cum_time, init_cum_time
   implicit none 
   type(vlist) :: kl
   type(cum_wgrid) :: wg
   type(cum_time) :: ctime
   logical :: success
   character(len=80) :: fname, dset_name
   real(dp) :: dw, e_qp
   real(dp), allocatable :: se_imag(:,:,:,:), resgm(:,:,:), enk(:,:), sp_cum(:,:,:,:), beta(:)
   integer(HID_T) :: file_id, group_id
   character(len=6), external :: int_to_char
   integer  :: it, ib, ik, numb, numk, beta_lw, beta_up, spf_lw, spf_up, &
               jst, jend, j, icount, ncount, step, nstep, iter
   
   if(ionode) write(stdout,'(5x,a)') 'Start calc_spectral_cumulant ...'
   if(ntemper .eq. 0) call errore('calc_spectral_cumulant', &
      'Set ftemper to file where temperatures(Fermi levels) are listed!',1)

   !readin k-points
   call load_vector_list(fklist, kl)
   numb = band_max - band_min + 1
   numk = kl%nvec
   
   !generate energy grid
   call init_cum_wgrid(cum_outer_emin, cum_outer_emax, &
      cum_inner_emin, cum_inner_emax, cum_de, cum_outer_np, wg, spectral_np)

   allocate( se_imag(size(wg%w), ntemper, numb, numk), &
             resgm(numb, numk, ntemper), enk(numb, numk) )
   !load Imsigma(w)
   if(ionode) write(stdout,'(5x,a)') '>loading fan selfenergy from file.'
   fname = trim(prefix)//'_selfenergy.h5'
   call load_selfenergy(fname, enk, se_imag, wg%w)
   !load ReSigma(w=0) = ReSigma(Enk)
   call load_imsigma(resgm, success, suffix='resigma')
   if(.not. success) &
      call errore('calc_trans_spectral','reading '// trim(prefix)// '.resigma failed ',1)
   
   !setup w-grid and t-grid for cumulant calculations
   dw = wg%sp_de
   beta_lw = wg%out_lw
   beta_up = wg%out_up
   spf_lw =   floor(spectral_emin / dw)
   spf_up = ceiling(spectral_emax / dw)
   call init_cum_time(beta_lw, beta_up, spf_lw, spf_up, ctime)
   !
   allocate( sp_cum(spf_lw:spf_up, ntemper, numb, numk) )
   sp_cum = 0.0_dp
   !
   write(stdout,'(5x,a)') '>computing spectral function.'

   call mp_split_pools(numk*numb*ntemper, jst, jend, ncount) !mpi parallel
   call set_progress_step(ncount, step, nstep)
   call progressbar_init('spectral_cumulant:')
   
   icount = 0  ! global counter
!$omp parallel default(shared) private(j, ik, ib, it, beta, e_qp, iter)
   allocate( beta(beta_lw:beta_up) )
!$omp do schedule(guided)
   do j = jst, jend
      ! j = (ik-1)*(numb*ntemper) + (ib-1)*ntemper + (it-1) + 1
      ik = (j-1) / (numb*ntemper) + 1
      ib = mod( (j-1)/ntemper, numb) + 1
      it = mod( (j-1), ntemper) + 1
   
      !extract beta(w) = -Im[Sigma](w) / pi, and interpolate it if required.
      call compute_beta(wg, se_imag(:,it,ib,ik), beta)
      !compute spectral function A(w)
      e_qp = resgm(ib, ik, it)
      call cumulant_spectral(beta, beta_lw, sp_cum(:,it,ib,ik), spf_lw, e_qp, dw, ctime)

      !track progress
!$omp atomic update
      icount = icount + 1
      !output progress
      iter = icount
      if( mod(iter, step).eq.0  .or. iter.eq.ncount ) then
!$omp critical (write_progress)
         write(stdout,'(8x, f7.2, a1)') (100.0_dp*iter)/ncount, '%'
!$omp end critical (write_progress)
      endif
   enddo
!$omp end do
   deallocate(beta)
!$omp end parallel
   call mp_sum(sp_cum, inter_pool_comm)
   
   !output results to hdf5 files
   if(ionode) then
      fname = trim(prefix)//'_spectral_cumulant.h5'
      call hdf_open_file(file_id, trim(fname), status='NEW')
      call hdf_write_dataset(file_id, 'kpoints', kl%vec)
      call hdf_write_dataset(file_id, 'temperatures', temper(1:ntemper)*ryd2ev/kelvin2eV)
      call hdf_write_dataset(file_id, 'fermi_levels', efermi(1:ntemper))
      call hdf_write_dataset(file_id, 'band_energy', enk)
      call hdf_write_dataset(file_id, 'w_lower_index', spf_lw)
      call hdf_write_dataset(file_id, 'w_upper_index', spf_up)
      call hdf_write_dataset(file_id, 'wfreq_step_eV', dw*ryd2ev)
      !
      call hdf_create_group(file_id, 'spectral_functions')
      call hdf_open_group(file_id, 'spectral_functions',group_id)
      do ik = 1, numk
         dset_name = "kpt_" // trim(int_to_char(ik))
         call hdf_write_dataset(group_id, trim(dset_name), sp_cum(:,:,:,ik))
      enddo
      call hdf_close_group(group_id)
      call hdf_close_file(file_id)
   endif

   deallocate(se_imag, resgm, sp_cum, enk)
end subroutine calc_spectral_cumulant

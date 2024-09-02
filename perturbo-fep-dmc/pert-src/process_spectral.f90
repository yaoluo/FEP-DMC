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

subroutine process_spectral()
   use pert_data, only: bg
   use pert_const, only: dp, pi, ryd2ev
   use qe_mpi_mod, only: ionode, mp_bcast, ionode_id, inter_pool_comm, mp_sum, mp_split_pools
   use pert_param, only: fklist, band_max, band_min, cum_inner_emin, cum_inner_emax, cum_de, &
      cum_outer_emin, cum_outer_emax, cum_outer_np, spectral_emin, spectral_emax, spectral_np, &
      ntemper, prefix, delta_smear, boltz_emax, boltz_emin
   use vector_list, only: load_vector_list, vlist
   use pert_output, only: load_imsigma
   use pert_utils,  only: find_free_unit
   use cumulant_expansion, only: cumulant_spectral, spectral_qp, cum_time, init_cum_time
   use cumulant_utility, only: cum_wgrid, init_cum_wgrid, compute_beta, load_selfenergy
   implicit none
   type(vlist) :: kl
   type(cum_wgrid) :: wg
   type(cum_time) :: ctime
   logical :: success
   real(dp) :: dw, e_qp
   
   integer :: numb, numk, nestep, iout, ios, reclen, beta_lw, beta_up, spf_lw, spf_up, &
      j, jst, jend, ik, ib, it, i
   real(dp), allocatable :: se_imag(:,:,:,:), enk(:,:), wtmp(:), resgm(:,:,:), &
      sp_cum(:,:,:,:), beta(:,:,:,:), sp_qp(:,:,:,:), enk_tmp(:,:), xloc(:)
   
   !readin k-points
   call load_vector_list(fklist, kl)
   numb = band_max - band_min + 1
   numk = kl%nvec
   !generate energy grid
   call init_cum_wgrid(cum_outer_emin, cum_outer_emax, &
      cum_inner_emin, cum_inner_emax, cum_de, cum_outer_np, wg, spectral_np)
   nestep = size(wg%w)
   
   allocate( se_imag(nestep, ntemper, numb, numk), enk_tmp(numb, numk), &
             resgm(numb, numk, ntemper), wtmp(nestep), enk(numb, numk) )
   se_imag = 0.0_dp;    enk = 0.0_dp;     wtmp = 0.0_dp;  resgm = 0.0_dp
   
   call load_selfenergy(trim(prefix)//'_selfenergy.h5', enk, se_imag, wg%w)

   call load_imsigma(resgm, success, enk_tmp, 'resigma')
   if( any(abs(enk_tmp - enk) > 1.0E-6_dp) ) &
      call errore('process_spectral','inconsistent resigma, selfenergy files',1)
   deallocate(enk_tmp, wtmp)

   if(ionode) write(*,*) "finishing loading data"

   beta_up = wg%out_up;  beta_lw = wg%out_lw;  dw = wg%sp_de;
   spf_lw = int(   floor(spectral_emin / dw) )
   spf_up = int( ceiling(spectral_emax / dw) )

   allocate( beta(  beta_lw:beta_up, ntemper, numb, numk), &
             sp_qp( spf_lw:spf_up, ntemper, numb, numk), &
             sp_cum(spf_lw:spf_up, ntemper, numb, numk) )
   
   !setup time grid for cumulant calculations
   call init_cum_time(beta_lw, beta_up, spf_lw, spf_up, ctime)
   
   if(ionode) write(*,*) "start cumulant_expansion"
   
   !perform cumulant calculation to obtain the spectral functions.
   call mp_split_pools( numk*numb*ntemper, jst, jend )
   sp_cum = 0.0_dp;  sp_qp = 0.0_dp;  beta = 0.0_dp
   !
!$omp parallel do schedule(guided) default(shared) private(j, ik, ib, it, i, e_qp)
   do j = jst, jend
      ! j = (ik-1)*(numb*ntemper) + (ib-1)*ntemper + (it-1) + 1
      ik = (j-1) / (numb*ntemper) + 1
      ib = mod( (j-1)/ntemper, numb) + 1
      it = mod( (j-1), ntemper) + 1
      
      !e_qp = ReSigma(w=0), note: w=0 means E = e_nk
      e_qp = resgm(ib, ik, it)
      call compute_beta(wg, se_imag(:,it,ib,ik), beta(:,it,ib,ik))

      !The cumulant spectral function
      call cumulant_spectral &
         (beta(:,it,ib,ik), beta_lw, sp_cum(:,it,ib,ik), spf_lw, e_qp, dw, ctime)
      !The Quasi-Particle spectral function
      do i = spf_lw, spf_up
         sp_qp(i,it,ib,ik) = spectral_qp(beta(0,it,ib,ik), e_qp, real(i, dp)*dw)
      enddo
      !
      write(*,1000) ik, ib, it, 'Cumulants: ', sum(sp_cum(:,it,ib,ik))*dw, &
         '; FanQuasiP: ', sum(sp_qp(:,it,ib,ik))*dw
      !if(ionode) write(*,1000) ik, ib, it, 'FanQuasiP: ', sum(sp_qp(:,it,ib,ik))*dw
      !   sum(sp_cum(:,it,ib,ik), mask=(sp_cum(:,it,ib,ik)>1.0E-4_dp))*dw
   enddo
!$omp end parallel do
   !
   call mp_sum(sp_cum, inter_pool_comm)
   call mp_sum(sp_qp , inter_pool_comm)
   call mp_sum(beta  , inter_pool_comm)

   !output spectral function: A(nk, w)
   if(ionode) then
      allocate( xloc(numk) )
      call generate_path(kl, bg, xloc)
      
      iout = find_free_unit()
      open(iout, file=trim(prefix)//".beta", form='formatted',status='unknown')
      write(iout,'(a)') "# kloc, enk, w, beta(w)"
      do it = 1, ntemper
      do ik = 1, numk
      do ib = 1, numb
         write(iout,'(a)') ' '
         do i = 1, nestep
            write(iout,1001) xloc(ik), enk(ib,ik)*ryd2ev, wg%w(i)*ryd2ev, &
               -se_imag(i,it,ib,ik)*ryd2ev/pi
         enddo
         write(iout,'(a)') ' '
      enddo; enddo; enddo
      close(iout)
      
      iout = find_free_unit()
      open(iout, file=trim(prefix)//".beta_fine", form='formatted',status='unknown')
      write(iout,'(a)') "# kloc, enk, w, beta(w) [interpolated]"
      do it = 1, ntemper
      do ik = 1, numk
      do ib = 1, numb
         write(iout,'(a)') ' '
         do i = beta_lw, beta_up
            write(iout,1001) xloc(ik), &
               enk(ib,ik)*ryd2ev, real(i,dp)*dw*ryd2ev, beta(i,it,ib,ik)*ryd2ev
         enddo
         write(iout,'(a)') ' '
      enddo; enddo; enddo
      close(iout)
      
      iout = find_free_unit()
      open(iout, file=trim(prefix)//".spectral_cum", form='formatted',status='unknown')
      write(iout,'(a)') "# kloc, enk, w, sp_cum(w), sp_qp(w)"
      do it = 1, ntemper
      do ik = 1, numk
      do ib = 1, numb
         write(iout,'(a)') ' '
         do i = spf_lw, spf_up
            write(iout,1001) xloc(ik), enk(ib,ik)*ryd2ev, real(i,dp)*dw*ryd2ev, &
               sp_cum(i,it,ib,ik)/ryd2ev, sp_qp(i,it,ib,ik)/ryd2ev
         enddo
         write(iout,'(a)') ' '
      enddo; enddo; enddo
      close(iout)
   endif
   return
101 call errore('process_spectral', 'open '//trim(prefix)//'.selfenergy', abs(ios))
102 call errore('process_spectral', 'read '//trim(prefix)//'.selfenergy', 1)
1000 format(1x,'# ik, ib, it:', i6, 1x, i3, 1x, i3, 1x, a, ES23.16, a, ES23.16)
1001 format(1x, 3(f12.6, 2x), 4(E15.8, 1x))
end subroutine process_spectral

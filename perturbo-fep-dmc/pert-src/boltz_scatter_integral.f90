!===============================================================================
! Copyright (C) 2016-2020 Jin-Jian Zhou, Jinsoo Park, I-Te Lu, Marco Bernardi
!
! This program is distributed under the terms of the GNU General Public License.
! See the file `LICENSE' in the root directory of this distribution, or obtain 
! a copy of the License at <https://www.gnu.org/licenses/gpl-3.0.txt>.
!
! Author: jjzhou <jjchou.comphy@gmail.com>
! Comment:
!  solve iterative linearized Boltzmann equation, and real-time dynamics
! 
! Maintenance:
!===============================================================================

module boltz_scatter_integral
   use pert_const, only: dp, pi
   use boltz_grid, only: grid
   use boltz_utils, only: idx2band
   use pert_utils, only: bose, fermi, gauss
   use pert_param, only: delta_smear
   use pert_output, only: stopwatch
   use boltz_scatter, only: num_scat, scat, qgrid
   use qe_mpi_mod, only: stdout, ionode, mp_barrier, mp_sum, inter_pool_comm
   implicit none
   private
   
   public :: rates_scat_int, trans_scat_int, cdyna_scat_int
contains

! compute scattering rate
subroutine rates_scat_int(kg, eptemp, efermi, rates)
   implicit none
   type(grid), intent(in) :: kg
   real(dp), intent(in) :: eptemp, efermi
   real(dp), intent(out) :: rates(:,:) !rates(kg%numb, kg%numk)
   !local variables
   integer :: i, iq, ik, ikq, ns, im, n, m, bands(3)
   real(dp) :: wq, enk, emkq, g2, dt1, dt2, ph_n, fkq, nk2mkq, mkq2nk
   
   call stopwatch('rates_scat_int', 'start')
   rates(:,:) = 0.0_dp
!$omp parallel do schedule(guided) default(shared) private(i, iq, ik, ikq, ns, &
!$omp& im, n, m, wq, enk, emkq, g2, dt1, dt2, ph_n, fkq, mkq2nk, nk2mkq, bands)
   do i = 1, num_scat
      iq = scat(i)%iq;    ik = scat(i)%ik;    ikq = scat(i)%ikq
      ! loop over all the (n, m, mu) pairs in (ik, ikq) of this pool
      do ns = 1, scat(i)%nchl
         ! map index to (mkq, nk, mu)
         bands = idx2band( scat(i)%bands_index(ns), kg%numb )
         m = bands(1);  n = bands(2);  im = bands(3)
         ! eigenvalue
         wq   = qgrid%freq(im, iq)
         enk  = kg%enk(n, ik)
         emkq = kg%enk(m, ikq)
         ! |g_{mu}(nk,mkq)|^2 .eq. |g_{mu}(mkq,nk)|^2
         g2 = scat(i)%eph_g2(ns)
         
         ! add the factor of pi and 2 (2/hbar, hbar is omitted in atomic unit.)
         ! delta(Enk - Emkq + wq) 
         dt1 = qgrid%weight * pi * 2.0_dp * gauss(delta_smear, enk - emkq + wq)
         ! delta(Enk - Emkq -wq)
         dt2 = qgrid%weight * pi * 2.0_dp * gauss(delta_smear, enk - emkq - wq)
         ! phonon occupation
         ! the constrain of 'wq > eps_acustic' is aready done in scat_setup.
         ph_n = bose(eptemp, wq) ! N_q\mu

         ! compute scattering rates
         fkq = fermi(eptemp, emkq-efermi)
         mkq2nk = g2 * (dt1*(ph_n + fkq) + dt2*(ph_n + 1.0_dp - fkq))
! lock memory when updating this variable
!$omp atomic
         rates(n, ik)  = rates(n, ik)  + mkq2nk

         !update mk+q with contribution from k
         fkq = fermi(eptemp, enk-efermi)
         !delta(Ekq-Ek+wq)=delta(Ek-Ekq-wq); delta(Ekq-Ek-wq)=delta(Ek-Ekq+wq)
         nk2mkq = g2 * (dt2*(ph_n + fkq) + dt1*(ph_n + 1.0_dp - fkq))
!$omp atomic
         rates(m, ikq) = rates(m, ikq) + nk2mkq
      enddo  !; enddo; enddo
   enddo ! i loop
!$omp end parallel do
   ! combine contributions from all the pools/processors.
   call mp_sum(rates, inter_pool_comm)
!   call mp_barrier(inter_pool_comm)
   call stopwatch('rates_scat_int', 'stop')
end subroutine rates_scat_int

! iterative process for tranport calculation
subroutine trans_scat_int(kg, eptemp, efermi, mfd, epint)
   implicit none
   type(grid), intent(in) :: kg
   real(dp), intent(in) :: efermi, eptemp, mfd(:,:,:) !mfd(3 or 6, kg%numb, kg%nk)
   real(dp), intent(out) :: epint(:,:,:)  !epint(3 or 6, kg%numb, kg%nk)
   !local variables
   integer :: i, iq, ik, ikq, ns, im, n, m, ic, bands(3), ncomp
   real(dp) :: wq, enk, emkq, g2, dt1, dt2, ph_n, nk2mkq, mkq2nk, fkq
   
   call stopwatch('trans_scat_int','start')
   !sanity check
   ncomp = size(mfd, 1)
   if( ncomp .ne. size(epint, 1) ) &
      call errore('trans_scat_int',"mismatch dimension in input arguments!",1)

   epint(:,:,:) = 0.0_dp
!$omp parallel do schedule(guided) default(shared) private(i, iq, ik, ikq, ns,&
!$omp&  im, n, m, wq, enk, emkq, g2, dt1, dt2, ph_n, fkq, mkq2nk, nk2mkq, ic, bands)
   do i = 1, num_scat
      iq = scat(i)%iq;   ik = scat(i)%ik;   ikq = scat(i)%ikq
      ! loop over all the (n, m, mu) pairs in (ik, ikq) of this pool
      do ns = 1, scat(i)%nchl
         ! map index to (mkq, nk, mu)
         bands = idx2band( scat(i)%bands_index(ns), kg%numb )
         m = bands(1);  n = bands(2);  im = bands(3)
         ! eigenvalue
         wq   = qgrid%freq(im, iq)
         enk  = kg%enk(n, ik)
         emkq = kg%enk(m, ikq)
         ! |g_{mu}(nk,mkq)|^2 .eq. |g_{mu}(mkq,nk)|^2
         g2 = scat(i)%eph_g2(ns)
         
         ! add the factor of pi and 2 (2/hbar, hbar is omitted in atomic unit.)
         ! delta(Enk - Emkq + wq) 
         dt1 = qgrid%weight * pi * 2.0_dp * gauss(delta_smear, enk - emkq + wq)
         ! delta(Enk - Emkq -wq)
         dt2 = qgrid%weight * pi * 2.0_dp * gauss(delta_smear, enk - emkq - wq)
         ! phonon occupation
         ! the constrain of 'wq > eps_acustic' is aready done in scat_setup.
         ph_n = bose(eptemp, wq) ! N_q\mu

         do ic = 1, ncomp
            !update nk with contribution from k+q
            fkq = fermi(eptemp, emkq-efermi)
            mkq2nk = g2 * mfd(ic, m, ikq) * (dt1*(ph_n + fkq) + dt2*(ph_n + 1.0_dp - fkq))
!$omp atomic
            epint(ic, n, ik)  = epint(ic, n, ik)  + mkq2nk
            
            ! update mk+q with contribution from k
            fkq = fermi(eptemp, enk-efermi)
            nk2mkq = g2 * mfd(ic, n, ik)  * (dt2*(ph_n + fkq) + dt1*(ph_n + 1.0_dp - fkq))
!$omp atomic
            epint(ic, m, ikq) = epint(ic, m, ikq) + nk2mkq
         enddo
      enddo
   enddo
!$omp end parallel do
   !combine contributions from all the pools/processors.
   call mp_sum(epint, inter_pool_comm)
   !call mp_barrier(inter_pool_comm)
   call stopwatch('trans_scat_int','stop')
end subroutine trans_scat_int

! Real time dynamics, Refer to eq.17. 18 in Eur. Phys. J. B (2016) 89: 239
! Bernardi, First-principles dynamics of electrons and phonons
subroutine cdyna_scat_int(kg, eptemp, dist, epcol)
   implicit none
   type(grid), intent(in) :: kg
   real(dp), intent(in) :: eptemp, dist(kg%numb, kg%nk)
   ! e-ph collision term computed from given distribution function dist(:,:)
   real(dp), intent(out) :: epcol(kg%numb, kg%nk)
   !local variables
   integer :: i, iq, ik, ikq, ns, im, n, m, bands(3)
   real(dp):: wq, enk, emkq, g2, dt1, dt2, ph_n, fnk, fmkq, fabs, fem, mkq2nk, nk2mkq

   epcol(:,:) = 0.0_dp
!$omp parallel do schedule(guided) default(shared) private(i, iq, ik, ikq, &
!$omp&  ns, im, n, m, wq, bands, enk, emkq, g2, dt1, dt2, ph_n, fnk, fmkq, &
!$omp&  fabs, fem, mkq2nk, nk2mkq)
   do i = 1, num_scat
      iq = scat(i)%iq;   ik = scat(i)%ik;    ikq = scat(i)%ikq
      ! loop over all the (n, m, mu) pairs in (ik, ikq) of this pool
      do ns = 1, scat(i)%nchl
         ! map index to (mkq, nk, mu)
         bands = idx2band( scat(i)%bands_index(ns), kg%numb )
         m = bands(1);  n = bands(2);  im = bands(3)
         ! eigenvalue
         wq   = qgrid%freq(im, iq)
         enk  = kg%enk(n, ik)
         emkq = kg%enk(m, ikq)
         ! |g_{mu}(nk,mkq)|^2 .eq. |g_{mu}(mkq,nk)|^2
         g2 = scat(i)%eph_g2(ns)
         ! add the factor of pi and 2 (2/hbar, hbar is omitted in atomic unit.)
         ! delta(Enk - Emkq + wq) 
         dt1 = qgrid%weight * pi * 2.0_dp * gauss(delta_smear, enk - emkq + wq)
         ! delta(Enk - Emkq -wq)
         dt2 = qgrid%weight * pi * 2.0_dp * gauss(delta_smear, enk - emkq - wq)
         ! phonon occupation
         ph_n = bose(eptemp, wq) ! N_q\mu

         ! electron distributaion
         fnk  = dist(n, ik)
         fmkq = dist(m, ikq)
         ! scattering process:  formula 1
         !fabs = fnk*(one - fmkq)*wgq - fmkq*(one - fnk)*(one + wgq)
         !fem  = fnk*(one - fmkq)*(one + wgq) - fmkq*(one - fnk)*wgq
         ! formula 2
         fabs = fnk*(ph_n + fmkq) - fmkq*(ph_n + 1.0_dp)
         fem  = fnk*(ph_n + 1.0_dp ) - fmkq*(ph_n + fnk)
         ! update nk with contribution from k+q
         mkq2nk = - g2 * (dt1*fabs + dt2*fem)
         !update nk with contribution from k+q
!$omp atomic
         epcol(n, ik) = epcol(n, ik) + mkq2nk
         ! the inverse: update mk+q with contribution from k
         nk2mkq = - mkq2nk
!$omp atomic
         epcol(m, ikq) = epcol(m, ikq) + nk2mkq
      enddo   !; enddo; enddo ! ns loop
   enddo ! i loop
!$omp end parallel do
   ! combine contributions from all the pools/processes.
   call mp_sum(epcol, inter_pool_comm)
!   call mp_barrier(inter_pool_comm)
end subroutine cdyna_scat_int

end module boltz_scatter_integral

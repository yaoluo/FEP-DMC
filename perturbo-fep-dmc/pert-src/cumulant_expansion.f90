!===============================================================================
! Copyright (C) 2016-2020 Jin-Jian Zhou, Jinsoo Park, I-Te Lu, Marco Bernardi
!
! This program is distributed under the terms of the GNU General Public License.
! See the file `LICENSE' in the root directory of this distribution, or obtain
! a copy of the License at <https://www.gnu.org/licenses/gpl-3.0.txt>.
!
! Author: jjzhou <jjchou.comphy@gmail.com>
! Comment:
!   compute spectral function using the cumulant expansion approach.
!   Needs the imaginary part of the Fan-Migdal self-energy as input.
!
! Maintenance:
!===============================================================================

module cumulant_expansion
   use pert_const, only: dp, pi, czero, twopi, cone, ryd2ev
   use qe_mpi_mod, only: mp_split_pools, inter_pool_comm, mp_sum
   use pert_utils, only: find_free_unit
   use pert_param, only: prefix
   implicit none
   private
   
   !for computing C^s(t), and C^s(w)
   type, public :: cum_time
      integer :: ntime  ! the range of t [-ntime, ntime]*dt, total number of t = 2*ntime+1
      !lower and upper boundary for C^s(w): cw_up - cw_lw = 2*ntime
      integer :: cw_lw, cw_up
      !contribution to C^s(t) from the lower extension.
      !cumt_lower(1:ntime),  C^s(t=0) = 0 so t=0 is excluded here.
      complex(dp), pointer :: cumt_lower(:) => null()
   end type cum_time
   
   public :: cumulant_spectral, spectral_qp, init_cum_time
contains

subroutine cumulant_spectral(beta, lb_beta, specfunc, lb_spf, ene_qp, dw, ctime)
   implicit none
   integer, intent(in) :: lb_beta, lb_spf
   ! ene_qp: quasi-particle peak, real part of Fan-Selfenergy at w=0
   real(dp), intent(in) ::  beta(lb_beta:), ene_qp, dw
   real(dp), intent(out) :: specfunc(lb_spf:)
   type(cum_time), intent(in) :: ctime
   !local
   integer :: ub_beta, ub_spf, nt, cw_l, cw_u, i, iout
   real(dp) :: sp_qp, cs_shift, ene_qp_shift
   real(dp),    allocatable :: cum_w(:)
   complex(dp), allocatable :: cs_t(:)
   
   ub_beta = ubound(beta, 1);    ub_spf = ubound(specfunc, 1)
   nt = ctime%ntime;    cw_l = ctime%cw_lw;    cw_u = ctime%cw_up
   !total number of time steps is 2*nt+1; and ub-lb = 2*nt
   allocate( cs_t(0:nt), cum_w(cw_l:cw_u) )

   !if beta(0) is much smaller than dw, then the position of qp peak
   ! will strongly affect the qp spectral (whether or not qp peak sit 
   ! at w-point makes a huge differnt on qp spectral function). 
   ! shift the quasi-particle peak slightly to avoild this.
   ene_qp_shift = (floor(ene_qp/dw) + 0.5_dp) * dw
   
   !compute C^s(t_i), t_i = i*dt, i=0,1,...,nt;  C^s(-t_i) = (C^s(t_i))*
   call cumulant_time(beta, lb_beta, dw, cs_t)
   !add contributions from w beyond the range of lb_beta:ub_beta
   call cumulant_time_extend(beta(0), dw, ctime, cs_t)
   !compute C(t) = exp(C^s(t))
   cs_t(0:nt) = exp( cs_t(0:nt) )
   !shift by real( exp(C^s(nt)) ) 
   cs_shift = real( cs_t(nt) )
   cs_t(0:nt) = cs_t(0:nt) - dcmplx(cs_shift, 0.0_dp)
   !transform to freq. domain and divided by 2pi: C(w)/2pi
   call cumulant_freq_direct(cs_t, cum_w, cw_l, dw)
   !compute convolution: C(w) (x) A^QP(w) 
   call direct_conv_qp(cum_w, cw_l, specfunc, lb_spf, beta(0), ene_qp_shift, dw)
   
   !compute the cumulant spectral function:
   ! A(w) = A^QP(w) * cs_shift + C(w) (x) A^QP(w)
   do i = lb_spf, ub_spf
      sp_qp = spectral_qp( beta(0), ene_qp_shift, real(i,dp)*dw )
      specfunc(i) = specfunc(i) + sp_qp * cs_shift
   enddo

   !debug
   !iout = find_free_unit()
   !open(iout, file=trim(prefix)//".cum", form='formatted',status='unknown')
   !do i = cw_l, cw_u
   !   write(iout,'(1x, f12.6, 2x, (E15.8, 1x))') real(i,dp)*dw*ryd2ev, cum_w(i)/ryd2ev
   !enddo
   !close(iout)
   !!debug
   !iout = find_free_unit()
   !open(iout, file=trim(prefix)//".cum_t", form='formatted',status='unknown')
   !do i = 1, nt
   !   write(iout,'(1x, f12.6, 2x, 2(E15.8, 2x))') real(i,dp), cs_t(i)
   !enddo
   !close(iout)
   deallocate(cum_w, cs_t)
end subroutine cumulant_spectral


subroutine init_cum_time(beta_lb, beta_ub, spf_lb, spf_ub, ctime)
   integer, intent(in) :: beta_lb, beta_ub, spf_lb, spf_ub
   type(cum_time), intent(out) :: ctime
   !
   integer :: cw_l, cw_u, nt, len_t, it, t_st, t_end
   real(dp) :: tt
   !a broad range of C(w) is needed when doing the convolution.
   cw_u = max(beta_ub, 2*spf_ub, -2*spf_lb)
   !determine ntime, the number of time steps sampled for C(t)
   ! T = 2pi/dw; dt = T/ntime; W = 2pi/dt; dw = W/nw; 
   ! therefore, ntime = T/dt = 2pi/dw / (2pi/W) = W/dw = nw
   ! we set t = [-nt:nt]*dt, so total #. of t = 2*nt+1, and dt = 2pi/dw/(2*nt+1)
   ! since t = t_i and t = -t_i are connected, only compute for t = t_i
   cw_l = -cw_u
   nt = cw_u - cw_l
   !the number of t-points is twice as the number of w-points.
   len_t = 2*nt + 1  !total number of t points.
   ctime%cw_lw = cw_l
   ctime%cw_up = cw_u
   ctime%ntime = nt

   if( associated(ctime%cumt_lower) ) deallocate(ctime%cumt_lower)
   allocate( ctime%cumt_lower(nt) )

   ctime%cumt_lower(nt) = czero
   call mp_split_pools(nt, t_st, t_end)
!$omp parallel do schedule(guided) default(shared) private(it, tt)
   do it = t_st, t_end
      tt = real(it, dp) / real(len_t, dp)
      !extend to 3*lb; 3 is hard-coded here, not sure if it's optimal.
      ctime%cumt_lower(it) = cum_time_correction(3*cw_l, beta_lb-1, tt)
   enddo
!$omp end parallel do
   call mp_sum(ctime%cumt_lower, inter_pool_comm)
end subroutine init_cum_time


!compute convolution of two array directly.
! conv = vecl (x) vecr 
subroutine direct_conv(vecl, lb_vecl, vecr, lb_vecr, conv, lb_conv, dw)
   implicit none
   integer, intent(in) :: lb_vecl, lb_vecr, lb_conv
   real(dp), intent(in) :: vecl(lb_vecl:), vecr(lb_vecr:), dw
   real(dp), intent(out) :: conv(lb_conv:)
   !local
   integer :: ub_vecl, ub_vecr, ub_conv, i, j, k
   
   ub_vecl = ubound(vecl, 1)
   ub_vecr = ubound(vecr, 1)
   ub_conv = ubound(conv, 1)
   
   conv = 0.0_dp
   do i = lb_conv, ub_conv
      do j = lb_vecl, ub_vecl
         k = i - j
         if( k < lb_vecr .or. k > ub_vecr ) cycle
         conv(i) = conv(i) + vecl(j) * vecr(k)
      enddo
   enddo
   conv(:) = conv(:) * dw
end subroutine direct_conv


!compute convolution f(w) (x) A^QP(w).
subroutine direct_conv_qp(vecl, lb_vecl, conv, lb_conv, beta0, qp_energy, dw)
   implicit none
   integer, intent(in) :: lb_vecl, lb_conv
   real(dp), intent(in) :: vecl(lb_vecl:), beta0, qp_energy, dw
   real(dp), intent(out) :: conv(lb_conv:)
   !local
   real(dp) :: w
   integer :: ub_vecl, ub_conv, i, j, k
   
   ub_vecl = ubound(vecl, 1)
   ub_conv = ubound(conv, 1)
   
   conv = 0.0_dp
   do i = lb_conv, ub_conv
      do j = lb_vecl, ub_vecl
         w = real(i - j, dp) * dw
         conv(i) = conv(i) + vecl(j) * spectral_qp(beta0, qp_energy, w)
      enddo
   enddo
   conv(:) = conv(:) * dw
end subroutine direct_conv_qp


!compute cumulant C^s in time domian: C^s(t_i); C^s(-t_i) = (C^s(t_i))*
! t_i = i*dt, i=0,1,..,nt; dt = 2*pi/dw/(2*nt+1)
subroutine cumulant_time(beta, lb, dw, cums_t)
   implicit none
   integer, intent(in) :: lb
   real(dp), intent(in) :: beta(lb:), dw  !beta(lb:ub)
   complex(dp), intent(out) :: cums_t(0:) !cums_t(0:nt)
   !local
   integer :: ub, nt, i, it, numt
   real(dp) :: beta_iw, wt, beta0
   
   ub = ubound(beta, 1); 
   nt = ubound(cums_t, 1); 
   if(lb > 0 .or. ub < 0) call errore('cumulant_time','beta should have index 0.', 1)
   
   numt = 2*nt + 1; 
   beta0 = beta(0)
   
   cums_t = czero
   do i = lb, ub
      if( i .eq. 0 ) cycle   !exclude w=0
      !compute (beta(w)-beta(0)) * dw/w^2, where w = i*dw
      beta_iw = (beta(i)-beta0) / (real(i*i, dp) * dw)
      
      !since C^s(t=0) = 0, so skip the first element it=0.
      do it = 1, nt
         !w*t = (i*dw) * (it*dt); dt = 2pi/dw/(2*nt+1); so w*t = i*it*2pi/(2*nt+1)
         wt = real(i*it, dp) / real(numt, dp) * twopi
         !C^s = sum_w (beta(w)-beta(0))*(exp(-iwt)-1)*dw/w^2, C^s(t) is dimensionless.
         cums_t(it) = cums_t(it) + beta_iw * dcmplx(cos(wt)-1.0_dp, -sin(wt))
      enddo
   enddo
end subroutine cumulant_time

! add contribution from w outside lb_beta:ub_beta to cums_t
subroutine cumulant_time_extend(beta0, dw, ctime, cums_t)
   implicit none
   real(dp), intent(in) :: beta0, dw
   type(cum_time), intent(in) :: ctime
   complex(dp), intent(inout) :: cums_t(0:) !cums_t(0:ctime%ntime)
   !local
   integer :: it, nt
   real(dp) :: beta_low

   nt = ctime%ntime
   beta_low = - beta0 / dw
   do it = 1, nt
      !contribution from w below lb
      cums_t(it) = cums_t(it) + beta_low * ctime%cumt_lower(it)
   enddo
end subroutine cumulant_time_extend

! t = tt * T = (it / numt) * T;  T = 2pi/dw
! compute sum_{j=wstart, wend} (e^{-i*w*t}-1)  / (w/dw)^2
complex(dp) pure function cum_time_correction(wstart, wend, tt) result(ct)
   integer,  intent(in) :: wstart, wend !NOTE: wstart < wend < 0.
   real(dp), intent(in) :: tt
   !
   integer :: i
   real(dp) :: wt
   
   ct = czero
   do i = wstart, wend
      !w = i*dw; t = tt*T = tt*2pi/dw;  so w*t = i*tt*2pi
      wt = real(i,dp) * tt * twopi
      ct = ct + dcmplx(cos(wt)-1.0_dp, -sin(wt)) / (real(i,dp)**2)
   enddo
end function cum_time_correction

!compute cumulant C^s in frequency domian: C^s(w_m)/2pi, using direct Fourier trans.
! C^s(w_m)/2pi = 1/2pi sum_{t_i} exp(iw_m*t_i) Ct(t_i)*dt; dt = 2pi/dw/(2*nt+1)
!              = 1/((2*nt+1)*dw) sum_{t_i} exp(iw_m*t_i) Ct(t_i); t_i = i*dt, i=-nt:nt
!Note: C^s(t) is dimensionless, so Ct(t) could be a function of C^s, e.g (C^s)^n
subroutine cumulant_freq_direct(ct, cw, lb, dw)
   implicit none
   integer, intent(in) :: lb
   real(dp), intent(in) :: dw
   complex(dp), intent(in) :: ct(0:)
   real(dp), intent(out) :: cw(lb:)
   !local
   integer :: ub, nt, i, it, numt
   real(dp) :: wt
   
   ub = ubound(cw, 1);  nt = ubound(ct, 1);  numt = 2*nt + 1

   cw(lb:ub) = real( ct(0) )  !contribution from it = 0
   do i = lb, ub
      do it = 1, nt
         wt = real(i*it, dp) / real(numt, dp) * twopi
         !Note that: ct(t_i) = ct(-t_i)*; exp[iw*(-t)]*Ct(-t) = [exp(iwt)*Ct(t)]*
         ! so exp[iw*(-t)]*Ct(-t) + exp(iwt)*Ct(t) = 2*real(exp(iwt)*Ct(t))
         cw(i) = cw(i) + 2.0_dp * real( dcmplx(cos(wt), sin(wt)) * ct(it) )
      enddo
   enddo
   cw(:) = cw(:) / (real(numt, dp)*dw)
end subroutine cumulant_freq_direct

!TODO:
!compute cumulant C^s in frequency domian: C^s(w_m), using fast Fourier trans.
!subroutine cumulant_freq_fft
!
!subroutine cumulant_freq_fft

! calc A^QP(w)
!! We use the real part of Fan-Selfenergy at w=0 for the energy shift
!! since the cauchy principle integral P Int -beta(w)/w dw = Re Sigma(w=0)
real(dp) elemental function spectral_qp(beta0, qp_energy, w) result(specfunc)
   implicit none
   !beta0 = -ImSigma(Enk)/pi, qp_energy: quasi-particle peak shift, Resigma(Enk)
   real(dp), intent(in) :: beta0, qp_energy, w

   specfunc = beta0 / ((w - qp_energy)**2 + (pi*beta0)**2)
end function spectral_qp

!calculate the shift of the quasi-particle peak
! sum_{i.ne.0} (beta(i) - beta(0)) / i = sum_{i.ne.0} ( beta(i) ) / i
! since beta(0) / i and beta(0) / (-1) cancels each other.
! however this summation/integral converge slowly w.r.o the energy range.
! and compute beta(w) in a broaden range is both expensive and problematic.
!------
!! Instead, we use the real part of Fan-Selfenergy at w=0 for the energy shift
!! since the cauchy principle integral P Int -beta(w)/w dw = Re Sigma(w=0)
!------
real(dp) pure function energy_shift(beta, lb) result(ek_c)
   implicit none
   integer, intent(in) :: lb
   real(dp), intent(in) :: beta(lb:)  !beta(lb:ub)
   !local 
   integer :: ub, i
   
   ek_c = 0.0_dp
   !check boundary
   ub = ubound(beta, 1)
   if(lb > 0 .or. ub < 0)  return
   
   do i = lb, ub
      if(i .eq. 0) cycle
      ek_c = ek_c - (beta(i) - beta(0)) / real(i, dp)
   enddo
end function energy_shift

end module cumulant_expansion

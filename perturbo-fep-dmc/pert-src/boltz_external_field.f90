!===============================================================================
! Copyright (C) 2016-2020 Jin-Jian Zhou, Jinsoo Park, I-Te Lu, Marco Bernardi
!
! This program is distributed under the terms of the GNU General Public License.
! See the file `LICENSE' in the root directory of this distribution, or obtain 
! a copy of the License at <https://www.gnu.org/licenses/gpl-3.0.txt>.
!
! Author: jjzhou <jjchou.comphy@gmail.com>
! Comment:
!  compute drift term in the Boltzmann Transport Equation due to external field.
!
! Maintenance:
!===============================================================================

module boltz_external_field
   use kinds, only: dp
   use pert_data, only: tpiba
   use boltz_grid, only: grid
   use qe_mpi_mod, only: mp_split_pools, mp_sum, inter_pool_comm
   implicit none
   private
   ! f = 1 / (e^g + 1), maxg is the max value of |g|
   real(dp), parameter :: maxg = 35.0E0_dp

   public :: calc_drift_efield_plain
   public :: calc_drift_efield  ! drift term due to electric field.
contains

subroutine calc_drift_efield_plain(kgrid, efield, dist, deriv)
   implicit none
   type(grid), intent(in) :: kgrid
   ! electric field e*E in Rydberg atomic unit (Ryd/Bohr).
   real(dp), intent(in)  :: efield(3), dist(kgrid%numb, kgrid%nk)
   real(dp), intent(out) :: deriv(kgrid%numb, kgrid%nk)
   ! local variables
   integer :: ik_start, ik_stop, ik, i, nb, num_bvec
   real(dp), allocatable :: evec(:)
   
   num_bvec = kgrid%num_neighbors
   allocate( evec(num_bvec) )
   
   ! compute E*bvec first.
   ! E*df/dk = E* sum_{b}{w_b * bvec * [f(k+b)-f(b)] }
   !  = sum_{b} { w_b * (E*bvec) * [f(k+b) - f(b)] }
   ! since w_b * bvec is in the unit of (tpiba/bohr)^-1, add the factor 
   ! of tpiba^-1 here, so that the E-field drift term is in atomic unit.
   do i = 1, num_bvec
      evec(i) = dot_product(efield(1:3), kgrid%bvec(1:3,i)) / tpiba
   enddo
   
   ! compute -qE* df/dk = eE* df/dk
   !   = sum_{b}{ w_b * (eE*bvec) * [f(k+b) - f(b)] }
   ! NB. omit 1/hbar here since hbar=1 in Rydberg a.u..
   deriv(:,:) = 0.0E0_dp
   ! mpi + openmp parallelization
   call mp_split_pools(kgrid%nk, ik_start, ik_stop)
!$omp parallel do schedule(guided) default(shared) private(ik, i, nb)
   do ik = ik_start, ik_stop
   do nb = 1, kgrid%numb
      do i = 1, num_bvec
         ! skip neighbors that are not in kgrid_k (assume f(k+b)-f(k) = 0)
         if( kgrid%neighbors(i,ik) .eq. 0 ) cycle  !not sure if it is appropriate.
         ! eE*df(k)/dk = \sum_{b} w_b * (E*bvec) * [f(k+b)-f(k)]
         deriv(nb, ik) = deriv(nb, ik) + kgrid%bweight(i) * evec(i) * &
            ( dist(nb,kgrid%neighbors(i,ik)) - dist(nb,ik) )
      enddo
   enddo
   enddo
!$omp end parallel do
   call mp_sum(deriv, inter_pool_comm)
   deallocate( evec )
end subroutine calc_drift_efield_plain

!  compute the drift term due to electric-field
!  The key quantities is the derivative of f(k): df(k)/dk. 
!  we use a different approach to compute the derivative numerically:
!  define f(k) = 1.0 / (exp(g(k)) + 1), then g(k) = ln( 1/f(k) - 1.0 )
!  so df(k)/dk = f(k)*(1-f(k))*dg(k)/dk.
!  by doing so, we constrain f(k) to be within (0,1).
!  compute dg(k)/dk numerically using finite difference formula:
!    dg(k)/dk = \sum_{b} w_b * bvec * [g(k+b)-g(k)]
!  this numerical derivative formula is exact for linear functions.
subroutine calc_drift_efield(kgrid, efield, dist, deriv)
   implicit none
   type(grid), intent(in) :: kgrid
   ! electric field e*E in Rydberg atomic unit (Ryd/Bohr).
   real(dp), intent(in)  :: efield(3), dist(kgrid%numb, kgrid%nk)
   real(dp), intent(out) :: deriv(kgrid%numb, kgrid%nk)
   ! local variables
   integer :: ik_start, ik_stop, ik, i, nb, num_bvec
   real(dp), allocatable :: evec(:), gdist(:,:)
   
   num_bvec = kgrid%num_neighbors
   allocate( evec(num_bvec), gdist(kgrid%numb, kgrid%nk) )
   ! prepare g( t_i )
   do ik = 1, kgrid%nk
   do nb = 1, kgrid%numb
      gdist(nb, ik) = f2g( dist(nb, ik) )
   enddo; enddo
   ! compute E*bvec first.
   ! E*dg/dk = E* sum_{b}{w_b * bvec * [g(k+b)-g(b)] }
   !  = sum_{b} { w_b * (E*bvec) * [g(k+b) - g(b)] }
   ! since w_b * bvec is in the unit of (tpiba/bohr)^-1, add the factor 
   ! of tpiba^-1 here, so that the E-field drift term is in atomic unit.
   do i = 1, num_bvec
      evec(i) = dot_product(efield(1:3), kgrid%bvec(1:3,i)) / tpiba
   enddo
   
   ! compute -qE* df/dk = eE* -(1-f)*f*dg/dk
   !   = (f-1)*f* sum_{b}{ w_b * (eE*bvec) * [g(k+b) - g(b)] }
   ! NB. omit 1/hbar here since hbar=1 in Rydberg a.u..
   deriv(:,:) = 0.0E0_dp
   ! mpi + openmp parallelization
   call mp_split_pools(kgrid%nk, ik_start, ik_stop)
!$omp parallel do schedule(guided) default(shared) private(ik, i, nb)
   do ik = ik_start, ik_stop
   do nb = 1, kgrid%numb
      do i = 1, num_bvec
         ! skip neighbors that are not in kgrid_k (assume g(k+b)-g(k) = 0)
         if( kgrid%neighbors(i,ik) .eq. 0 ) cycle  !not sure if it is appropriate.
         ! eE*dg(k)/dk = \sum_{b} w_b * (E*bvec) * [g(k+b)-g(k)]
         deriv(nb, ik) = deriv(nb, ik) + kgrid%bweight(i) * evec(i) * &
            ( gdist(nb,kgrid%neighbors(i,ik)) - gdist(nb,ik) )
      enddo
      deriv(nb,ik) = deriv(nb,ik) * dist(nb,ik) * ( dist(nb,ik)-1.0E0_dp )
   enddo
   enddo
!$omp end parallel do
   call mp_sum(deriv, inter_pool_comm)
   deallocate( evec, gdist )
end subroutine calc_drift_efield

!inverse function of f=1/(e^g + 1), input f, return g.
pure function f2g( fval )
   implicit none
   ! fval should be within (0, 1)
   real(dp), intent(in) :: fval 
   real(dp) :: emg, omf, f2g
   
   !fval < e^-maxg (f->0), return maxg; 1.0-fval < e^-maxg (f->1), ruturn -maxg
   ! g = ln( (1/f - 1)) = ln( (1-f)/f ) = ln(1-f) - ln(f)
   ! f should be within (0, 1), slightly lower than 0, or higher than 1 is fine.
   emg = exp(-maxg);   omf = 1.0E0_dp - fval
   f2g = merge( log(omf), -maxg, omf>emg ) - merge( log(fval), -maxg, fval>emg )
end function f2g

!define f = 1/(e^g + 1), input g, return f.
pure function g2f( gval )
   implicit none
   real(dp), intent(in) :: gval
   real(dp) :: x, g2f

   ! -maxg < gval < maxg: f = 1/(e^{gval} + 1)
   ! gval < -maxg : f = 1/(e^{-maxg} + 1) -> 1.0, avoid overflow
   ! gval >  maxg : f = 1/(e^{ maxg} + 1) -> 0.0, avoid overflow
   x = sign( min(abs(gval), maxg), gval)
   g2f = 1.0E0_dp / (1.0E0_dp + exp(x) )
end function g2f

end module boltz_external_field


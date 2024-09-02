!===============================================================================
! Copyright (C) 2016-2020 Jin-Jian Zhou, Jinsoo Park, I-Te Lu, Marco Bernardi
!
! This program is distributed under the terms of the GNU General Public License.
! See the file `LICENSE' in the root directory of this distribution, or obtain 
! a copy of the License at <https://www.gnu.org/licenses/gpl-3.0.txt>.
!
! Author: jjzhou <jjchou.comphy@gmail.com>
! Comment:
!  algorithms to solve the full BTE equation
!
! Maintenance:
!===============================================================================

module boltz_dynamics_solver
   use pert_const, only: dp
   use boltz_grid, only: grid
   use boltz_scatter_integral, only: cdyna_scat_int
   implicit none
   private
   
   public :: euler
   public :: runge_kutta_4th
contains

subroutine euler(kgrid, tstep, cdist, ndist, ept)
   implicit none
   type(grid), intent(in) :: kgrid
   ! time step, temperature for e-ph, and electric field in atomic unit
   real(dp), intent(in)  :: tstep, ept
   ! cdist: f( t_i ); ndist: f( t_(i+1) )
   real(dp), intent(in)  :: cdist(kgrid%numb, kgrid%nk)
   real(dp), intent(out) :: ndist(kgrid%numb, kgrid%nk)
   
   ! compute electron-phonon collision term
   call cdyna_scat_int(kgrid, ept, cdist, ndist)
   !update 
   ndist(:,:) = cdist(:,:) + ndist(:,:) * tstep

end subroutine euler

subroutine runge_kutta_4th(kgrid, tstep, cdist, ndist, df_tmp, dertot, ept)
   implicit none
   type(grid), intent(in) :: kgrid
   ! time step, temperature for e-ph, and electric field in atomic unit
   real(dp), intent(in)  :: tstep, ept
   ! cdist: f( t_i ); ndist: f( t_(i+1) )
   real(dp), intent(in)  :: cdist(kgrid%numb, kgrid%nk)
   real(dp), intent(out) :: ndist(kgrid%numb, kgrid%nk)
   ! workspace, to avoid frequently allocate and deallocate large arrays,
   !   which could hurt efficiency and parallel performance severely.
   real(dp), intent(out) :: df_tmp(kgrid%numb, kgrid%nk), &
                            dertot(kgrid%numb, kgrid%nk)
   ! local variables
   real(dp) :: half_t, t_six, t_three
   
   half_t = tstep*0.5E0_dp;  t_three = tstep/3.E0_dp;  t_six = tstep/6.E0_dp
   !
   df_tmp = 0.0E0_dp; dertot = 0.0E0_dp
   ! compute k1
   ! compute electron-phonon collision term
   call cdyna_scat_int(kgrid, ept, cdist, dertot)
   ! collect contribution from k1
   ndist(:,:) = cdist(:,:) + dertot(:,:) * t_six

   ! compute k2
   df_tmp(:,:) = cdist(:,:) + half_t * dertot(:,:)
   call cdyna_scat_int(kgrid, ept, df_tmp, dertot)
   ! collect contribution from k2
   ndist(:,:) = ndist(:,:) + dertot(:,:) * t_three

   ! compute k3
   df_tmp(:,:) = cdist(:,:) + half_t * dertot(:,:)
   call cdyna_scat_int(kgrid, ept, df_tmp, dertot)
   ! collect contribution from k3
   ndist(:,:) = ndist(:,:) + dertot(:,:) * t_three

   ! compute k4
   df_tmp(:,:) = cdist(:,:) + tstep * dertot(:,:)
   call cdyna_scat_int(kgrid, ept, df_tmp, dertot)
   ! collect contribution from k4
   ndist(:,:) = ndist(:,:) + dertot(:,:) * t_six

end subroutine runge_kutta_4th

end module boltz_dynamics_solver

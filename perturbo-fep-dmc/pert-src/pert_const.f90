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

module pert_const
   use kinds, only: dp
   implicit none
   public
   !double precision
   !integer,  parameter :: dp = kind(1.0d0)
   !to be consistent with QE, set dp to be the same as QE/kinds
   !integer,  parameter :: dp = selected_real_kind(14,200)
   !constants
   real(dp), parameter :: pi = 3.141592653589793238462643383279_dp
   real(dp), parameter :: twopi = 2.0_dp*pi
   real(dp), parameter :: ryd2ev  = 13.605698066_dp
   real(dp), parameter :: ryd2mev = ryd2ev*1.0E3_dp
   real(dp), parameter :: kelvin2eV = 8.6173427909E-05_dp
   real(dp), parameter :: bohr2ang = 0.52917721092_dp
   real(dp), parameter :: e2 = 2.0_dp  !e^2 in Rydberg atomic unit
   real(dp), parameter :: timeunit = 4.8377687E-17_dp !Rydberg atomic unit t0 in s
   real(dp), parameter :: unitcharge = 1.60217733E-19_dp !electron charge e in C
   
   complex(dp), parameter :: czero = (0.0_dp,0.0_dp)
   complex(dp), parameter :: cone  = (1.0_dp,0.0_dp)
   complex(dp), parameter :: ci    = (0.0_dp,1.0_dp)
end module pert_const

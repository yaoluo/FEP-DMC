!===============================================================================
! Copyright (C) 2016-2020 Jin-Jian Zhou, Jinsoo Park, I-Te Lu, Marco Bernardi
!
! This program is distributed under the terms of the GNU General Public License.
! See the file `LICENSE' in the root directory of this distribution, or obtain 
! a copy of the License at <https://www.gnu.org/licenses/gpl-3.0.txt>.
!
! Author: jjzhou <jjchou.comphy@gmail.com>
! Comment:
!  store basic information about the system.
!  most of the variables are copied from quantum esspresso
!  all in Rydberg atomic unit unless indicated explicitly.
!
! Maintenance:
!===============================================================================

module pert_data
   use hdf5_utils
   use pert_const, only: dp
   implicit none
   public
   save
   !! ions_bases in QE
   !number of atoms in the system
   integer :: nat
   integer :: num_wann
   integer :: kc_dim(3)
   integer :: qc_dim(3)
   ! number of types (species)
   !integer :: nsp
   !the type of i-th atom
   !integer,  allocatable :: ityp(:) ! ityp(nat)
   !position of the i-th atom, in cart. coord., unit of alat(same as that in QE)
   real(dp), allocatable :: tau(:,:) ! tau(3,nat)
   !mass of ions, in atomic mass unit
   real(dp), allocatable :: mass(:) ! amass(nat)
   !
   real(dp), allocatable :: wannier_center(:,:)
   real(dp), allocatable :: wannier_center_cryst(:,:)
   !! cell_base in QE
   !volume in au^3
   real(dp) :: volume
   !alat: lattice parameter - often used to scale quantities, or
   !in combination to other parameters/constants to define new units
   real(dp) :: alat
   real(dp) :: tpiba ! = 2*pi/alat
   !direct and reciprocal lattice primitive vectors
   !at(:,i) are the lattice vectors of the simulation cell, a_i,
   !  in alat units: a_i(:) = at(:,i)/alat
   !bg(:,i) are the reciprocal lattice vectors, b_i,
   !  in tpiba=2pi/alat units: b_i(:) = bg(:,i)/tpiba
   !  at^T * bg = I; so bg^-1 = transpose(at); at^-1 = transpose(bg)
   real(dp) :: at(3,3)
   real(dp) :: bg(3,3)
   !symmetry 
   integer :: nsym
   integer :: symop(3,3,48)
   
   ! .true. if it's 2D material, assume at(:,3) is direction that breaks 
   ! periodicity, and at(1:2, 3) should be 0, a(:,1:2) should be in-plane
   logical  :: system_2d  

   !for polar system
   logical  :: lpolar       ! .true. if this is a polar materials
   real(dp) :: polar_alpha  ! Ewald parameter for polar e-ph 
   real(dp) :: loto_alpha   ! Ewald parameter for phonon dispersion
   real(dp) :: thickness_2d !thickness of 2D system, negative for 3D system.
   real(dp) :: epsil(3,3)
   real(dp), allocatable :: zstar(:,:,:) !zstar(3,3,nat)
   !
   integer(HID_T) :: epwan_fid
end module pert_data

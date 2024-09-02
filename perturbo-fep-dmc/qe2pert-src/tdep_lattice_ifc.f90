!===============================================================================
! Copyright (C) 2016-2020 Jin-Jian Zhou, Jinsoo Park, I-Te Lu, Marco Bernardi
!
! This program is distributed under the terms of the GNU General Public License.
! See the file `LICENSE' in the root directory of this distribution, or obtain 
! a copy of the License at <https://www.gnu.org/licenses/gpl-3.0.txt>.
!
! Author: jjzhou <jjchou.comphy@gmail.com>
! Comment:
!   obtain inter-atomic force constants from TDEP
!
! Maintenance:
!===============================================================================

subroutine tdep_lattice_ifc(ph, qdim, nat, at, cryst_tau)
   use kinds, only: dp
   use qe_mpi_mod, only: stdout
   use force_constant, only: lattice_ifc, set_ws_cell_ph
   !
#if defined(__TDEP)
   use lattice_tdep, only: uc_tdep, fc_tdep, init_ifc_tdep
#endif
   !
   implicit none
   type(lattice_ifc), intent(inout) :: ph
   integer,  intent(in) :: qdim(3), nat
   ! at(:,i), lattice vector a_i (scaled by alat)
   ! cryst_tau: atomic position in crystal coordinate
   real(dp), intent(in) :: at(3, 3), cryst_tau(3, nat)

#if defined(__TDEP)
   !setup wigner seitz cell
   call set_ws_cell_ph(ph, qdim, nat, at, cryst_tau)

   write(stdout,'(/5x,a)') "Interface to TDEP: "
   write(stdout,'( 5x,a)') "----------------------------------------------"
   write(stdout,'( 5x,a)') "read infile.ucposcar and infile.forceconstant."
   !read TDEP unit cell
   call uc_tdep%readfromfile('infile.ucposcar')
   !read TDEP forceconstant
   call fc_tdep%readfromfile(uc_tdep, 'infile.forceconstant')

   !init ph with TDEP data
   call init_ifc_tdep(uc_tdep, fc_tdep, ph, at, cryst_tau)
   !
#else
   call errore('tdep_lattice_ifc',"Cannot load TDEP data. Recompile Perturbo with -D__TDEP.",1)
#endif

end subroutine tdep_lattice_ifc

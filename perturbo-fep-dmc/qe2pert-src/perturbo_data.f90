!===============================================================================
! Copyright (C) 2016-2020 Jin-Jian Zhou, Jinsoo Park, I-Te Lu, Marco Bernardi
!
! This program is distributed under the terms of the GNU General Public License.
! See the file `LICENSE' in the root directory of this distribution, or obtain 
! a copy of the License at <https://www.gnu.org/licenses/gpl-3.0.txt>.
!
! Author: jjzhou <jjchou.comphy@gmail.com>
! Comment:
!  save g(R_e, R_p), H(R_e), \Phi(R_p), and other basic info. to a hdf5 file.
!
! Maintenance:
!===============================================================================

module perturbo_data
   use kinds, only: dp
   use ions_base, only: amass, ityp, nat, tau
   use cell_base, only: omega, at, bg, alat
   use hdf5_utils
   implicit none
   public

   logical, save :: lpolar
   ! .true. if the system is a polar insulator/semiconductor
   real(dp), save :: epsil(3,3)
   ! dielectric tensor
   real(dp), pointer, save :: zstar(:,:,:)
   ! born effective charge, zstar(3,3,nat)
   real(dp), save :: loto_alpha !Ewald parameter for dipole-dipole interaction

   public :: save_perturbo_data
contains

subroutine save_perturbo_data()
   use io_global, only: stdout, ionode
   use qe_mpi_mod, only: mp_barrier, inter_pool_comm
   use lattice_data, only: qdim, numq, xq_tot, dynmat, ph_lpolar, ph_zstar, ph_epsil, ph_alpha
   use electronic_data, only: wannier_center_cryst, xk_tot, et_sub, et_corr, rot_wan, num_kpts
   use input_param, only: num_wann, prefix, kdim, num_band, debug, eig_corr, tdep
   !
   use force_constant, only: lattice_ifc
   use electron_wannier, only: electron_wann
   !
   use epwan_hdf5_io, only: write_force_constant, write_electron_wannier
   implicit none
   character(len=80) :: fname
   !atomic position in crystal coordinate
   real(dp) :: cryst_tau(3, nat)
   integer(HID_T), save :: file_id
   integer, external :: find_free_unit
   !
   type(lattice_ifc) :: phon
   type(electron_wann) :: elec
   
   cryst_tau(:,:) = tau(:,:)
   call cryst_to_cart(nat, cryst_tau, bg, -1)

   !forceconstant
   if(tdep) then
      call tdep_lattice_ifc(phon, qdim, nat, at, cryst_tau)
      !
      lpolar = phon%lpol
      epsil  = phon%pol%epsil
      zstar  => phon%pol%bcharge
      loto_alpha = phon%pol%alpha
   else
      call calc_lattice_ifc(phon, qdim, nat, at, cryst_tau, numq, xq_tot, dynmat)
      !
      lpolar = ph_lpolar
      epsil  = ph_epsil
      zstar  => ph_zstar
      loto_alpha = ph_alpha
   endif

   !electron
   if( trim(eig_corr) .eq. '' ) then
      call calc_electron_wann(elec, kdim, num_wann, at, &
         wannier_center_cryst, num_band, num_kpts, xk_tot, et_sub, rot_wan)
   else
      call calc_electron_wann(elec, kdim, num_wann, at, &
         wannier_center_cryst, num_band, num_kpts, xk_tot, et_corr, rot_wan)
   endif
   !
   !write out basic info
   if(ionode) then
      !output to hdf5 file
      fname = trim(prefix)//"_epwan.h5"
      call hdf_open_file(file_id, trim(fname), status='NEW')
      !
      call write_basic_data_hdf5(file_id)
      !
      call write_force_constant(file_id, phon)
      !
      call write_electron_wannier(file_id, elec)
      !close file
      call hdf_close_file(file_id)
      !
      write(stdout, '(5x,a)') "saved electron and phonon data to file."
   endif

   call mp_barrier( inter_pool_comm )
end subroutine save_perturbo_data


subroutine write_basic_data_hdf5(file_id)
   use lattice_data, only: qdim
   use electronic_data, only: wannier_center, wannier_center_cryst, noncolin
   use input_param, only: num_wann, kdim, polar_alpha, thickness_2d, system_2d, tdep
   use symm_base, only: s, nsym
   implicit none
   integer(HID_T), intent(in) :: file_id
   !local
   integer :: i
   integer(HID_T) :: group_id
   real(dp) :: tmp_mass(nat)

   do i = 1, nat
      tmp_mass(i) = amass(ityp(i))
   enddo

   call hdf_create_group(file_id, 'basic_data')
   call hdf_open_group(file_id, 'basic_data', group_id)
   !
   call hdf_write_dataset(group_id, "num_wann", num_wann)
   call hdf_write_dataset(group_id, "nat",       nat)
   call hdf_write_dataset(group_id, "kc_dim",  kdim)
   call hdf_write_dataset(group_id, "qc_dim",  qdim)
   call hdf_write_dataset(group_id, "volume",  omega)
   call hdf_write_dataset(group_id, "epsil",   epsil)
   call hdf_write_dataset(group_id, "alat",     alat)
   call hdf_write_dataset(group_id, "at",         at)
   call hdf_write_dataset(group_id, "bg",         bg)
   call hdf_write_dataset(group_id, "nsym",  nsym)
   call hdf_write_dataset(group_id, "symop", s(:,:,1:nsym))
   ! Ewald parameter for e-ph polar correction only!
   call hdf_write_dataset(group_id, "polar_alpha", polar_alpha)
   call hdf_write_dataset(group_id, "loto_alpha",   loto_alpha)
   call hdf_write_dataset(group_id, 'spinor', merge(1, 0, noncolin))
   call hdf_write_dataset(group_id, 'lpolar', merge(1, 0, lpolar))
   !
   if(tdep) call hdf_write_dataset(group_id, 'tdep', 1)
   
   call hdf_write_dataset(group_id, "system_2d", merge(1, 0, system_2d))
   if(system_2d) call hdf_write_dataset(group_id, "thickness_2d", thickness_2d)

   call hdf_write_dataset(group_id, "mass",  tmp_mass)
   call hdf_write_dataset(group_id, "zstar", zstar(1:3, 1:3, 1:nat) )
   call hdf_write_dataset(group_id, "tau",   tau(:, 1:nat) )
   call hdf_write_dataset(group_id, "wannier_center", wannier_center(:,1:num_wann))
   call hdf_write_dataset(group_id, "wannier_center_cryst", wannier_center_cryst(:,1:num_wann))

   call hdf_close_group(group_id)
end subroutine write_basic_data_hdf5

end module perturbo_data

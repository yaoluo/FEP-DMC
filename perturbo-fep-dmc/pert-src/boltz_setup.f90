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

subroutine boltz_setup()
   use pert_const, only: dp, ryd2ev, bohr2ang
   use pert_data, only: nsym, symop, volume, epwan_fid, qc_dim, kc_dim, system_2d, alat, at
   use qe_mpi_mod, only: ionode, stdout, mp_barrier, inter_pool_comm
   use boltz_grid, only: grid, init_boltz_grid, boltz_grid_generate
   use pert_param, only: ntemper, temper, efermi, ftemper, hole, find_efermi, band_min, &
      band_max, boltz_de, boltz_emax, boltz_emin, spinor, boltz_kdim, doping, debug
   use boltz_trans_mod, only: trans_density_calc
   use boltz_trans_output, only: output_density, output_ftemper, output_dos, output_kgrid
   use band_structure, only: electron_wann, init_electron_wann
   implicit none
   type(grid) :: kg
   integer :: nestep, i
   real(dp) :: emin, emax
   real(dp), allocatable :: ene(:), dos(:), dens(:)
   type(electron_wann) :: elec
   !
   if(ntemper < 1) then
      write(stdout, '(5x, a)') &
         "Warn (setup): ftemper is not specified. skip DOS and N_c calculations."
      return
   endif
   !
   call init_electron_wann(epwan_fid, kc_dim, elec)
   !check if we need to generate k-grid 
   if( check_exist_tet(boltz_kdim, boltz_emin, boltz_emax, band_min, band_max) ) then
      !load from pre-computed prefix_tet.h5 file
      call init_boltz_grid(kg, elec, band_min, band_max)
   else
      write(stdout,'(5x,a)') '>generating grid ...'
      !setup mode: generate kgrid, compute density(doping) or proper fermi-enery.
      call boltz_grid_generate(kg, boltz_kdim, boltz_emin, &
         boltz_emax, band_min, band_max, nsym, symop, elec)
      !
      call init_boltz_grid(kg, elec, band_min, band_max, .false.)
   endif
   call mp_barrier(inter_pool_comm)
   
   ! setup up energy windows for transport calculations
   emin = max(boltz_emin, minval(kg%enk))
   emax = min(boltz_emax, maxval(kg%enk))
   !if set energy range to boltz_emin and boltz_emax if specified
   if(find_efermi) then
      if(boltz_emin > -9000.0_dp) emin = boltz_emin
      if(boltz_emax <  9000.0_dp) emax = boltz_emax
   endif
   write(stdout,'(5x, a, 2(1x,f12.6), /)') 'Energy window (eV):', emin*ryd2ev, emax*ryd2ev
   if(emin > emax) call errore('boltz_setup','illegal energy window',1)
   
   !setup energy grid: emin - boltz_de : emax + boltz_de
   nestep = int( floor((emax-emin)/boltz_de) ) + 3
   allocate(ene(nestep), dos(nestep))
   do i = 1, nestep
      ene(i) = emin + (i-2)*boltz_de
   enddo

   allocate( dens( ntemper ) )
   dens = 0.0E0_dp
   if( find_efermi ) then
      if(system_2d) then
         !for 2D from #./cm^2 to #./bohr^2
         dens(:) = doping(:) * 1.0E-16_dp * (bohr2ang)**2
         !convert to its equivalent 3D counterpart
         dens(:) = dens(:) / (alat * at(3,3))
      else
         !for 3D from #./cm^3 to #./bohr^3
         dens(:) = doping(:) * 1.0E-24_dp * (bohr2ang)**3
      endif
   endif
   !
   !if both setup and find_efermi is true, then compute efermi from dens
   call trans_density_calc &
      (kg, ene, temper, efermi, hole, spinor, volume, dens, find_efermi, dos)

   !that is all for setup mode, output density(doping) info
   if(ionode) then
      call output_dos(ene, dos)
      call output_density(temper, efermi, dens)
      if(debug) call output_kgrid(kg)
      if(find_efermi) call output_ftemper(temper, efermi, dens, ftemper)
   endif
   deallocate(ene, dos, dens)
   !
   return

   contains

   logical function check_exist_tet(kdim, emin, emax, bmin, bmax) result(lexist)
      use pert_const, only: dp
      use pert_param, only: prefix
      use hdf5_utils
      implicit none
      integer, intent(in) :: kdim(3), bmin, bmax
      real(dp), intent(in) :: emin, emax
      !local
      logical :: has_file
      character(len=120) :: fname
      integer(HID_T) :: file_id
      integer :: itmp(3), brange(2)
      real(dp) :: erange(2)
      
      lexist = .false.
      !
      fname = trim(prefix)//"_tet.h5"
      inquire(file=trim(fname), exist=has_file)
      if(.not. has_file) return
      !
      call hdf_open_file(file_id, trim(fname), status='OLD', action='READ')
      call hdf_read_dataset(file_id, 'kgrid_dim', itmp(1:3))
      call hdf_read_dataset(file_id, 'band_window', brange)
      call hdf_read_dataset(file_id, 'energy_window', erange)
      call hdf_close_file(file_id)
      !
      if( any( kdim(1:3) .ne. itmp(1:3) ) ) return
      if( brange(1) .ne. bmin .or. brange(2) .ne. bmax )  return
      if( abs(erange(1)-emin) > 1.0E-12_dp .or. &
          abs(erange(2)-emax) > 1.0E-12_dp ) return
      
      lexist = .true.
   end function check_exist_tet

end subroutine boltz_setup

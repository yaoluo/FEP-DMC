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

subroutine carrier_dynamics_run()
   use pert_const, only: dp
   use qe_mpi_mod, only: ionode, stdout
   use boltz_grid, only: grid, init_boltz_grid
   use pert_data,  only: epwan_fid, qc_dim, kc_dim
   use pert_param, only: ntemper, temper, band_min, band_max, prefix, boltz_qdim, &
      boltz_nstep, time_step, output_nstep, solver
   !
   use boltz_dynamics_mod, only: output_dist
   use pert_output,only: progressbar_init, progressbar
   use boltz_dynamics_solver, only: runge_kutta_4th, euler
   !
   use boltz_scatter, only: boltz_scatter_setup
   use band_structure, only: electron_wann, init_electron_wann
   use phonon_dispersion, only: lattice_ifc, init_lattice_ifc
   use hdf5_utils
   implicit none
   !local variables
   integer  :: i, nstep
   character(len=120) :: fname
   integer(HID_T) :: file_id, group_id
   !distribution function: f_nk(t_n), f_nk(t_n+1), k_nk used in r-k
   real(dp), allocatable:: dist0(:,:), dist1(:,:), dist_t1(:,:), dist_t2(:,:)
   !
   type(grid) :: kg
   type(lattice_ifc)   :: phon
   type(electron_wann) :: elec
   
   !only needs one temperature for phonon occupation N(mod, q)
   if(ntemper > 1) write(stdout,'(5x, a)') &
      "Warn (carrier_dynamics_run): only the first temperature is used."
   !init
   call init_lattice_ifc(epwan_fid, qc_dim, phon)
   call init_electron_wann(epwan_fid, kc_dim, elec)
   !
   !setup k-grid
   call init_boltz_grid(kg, elec, band_min, band_max)
   !setup q-grid, and e-ph scattering channel (a.k.a k-q pair).
   call boltz_scatter_setup(kg, elec, phon, boltz_qdim)
   !
   !allocate space for distribution function
   allocate( dist0(kg%numb, kg%nk), dist1(kg%numb, kg%nk) )
   ! allocate workspace if solver is rk4
   if(trim(solver) .eq. 'rk4') &
      allocate( dist_t1(kg%numb, kg%nk), dist_t2(kg%numb, kg%nk) )
   !
   dist0 = 0.0E0_dp;    dist1 = 0.0E0_dp;
   
   !initialize distribution function: 
   !   restart from previous run or start a new one.
   fname = trim(prefix)//"_cdyna.h5"
   ! get initial distribution
   call cdyna_setup(fname, kg, file_id, group_id, dist1)

   !start dynamics simulation, write the initial step
   if(ionode) call progressbar_init('Carrier Dynamics:')
   nstep = (boltz_nstep / output_nstep) * output_nstep
   !
   do i = 1, nstep
      dist0(:,:) = dist1(:,:)
      !
      !get the f(t_i+1)
      if(trim(solver) .eq. 'rk4') then
         call runge_kutta_4th &
            (kg, time_step, dist0, dist1, dist_t1, dist_t2, temper(1))
      else
         call euler(kg, time_step, dist0, dist1, temper(1))
      endif
      !
      !output f(t_i+1) if needed.
      if( mod(i, output_nstep) .eq. 0 ) then
         if(ionode) then
            call progressbar(i, boltz_nstep)
            call output_dist(kg, group_id, dist1, (i/output_nstep) )
         endif
      endif
   enddo
   !
   !close file
   if(ionode) then
      call hdf_close_group(group_id)
      call hdf_close_file(file_id)
   endif
   ! release work space
   if(trim(solver) .eq. 'rk4') deallocate(dist_t1, dist_t2)
   deallocate(dist0, dist1)
   return

   contains
      !
   subroutine cdyna_setup(filename, kgrid, fid, gid, dist)
      use pert_const, only: dp, timeunit, ryd2ev
      use qe_mpi_mod, only: ionode, mp_barrier, mp_bcast, ionode_id, inter_pool_comm
      use boltz_dynamics_mod, only: output_dist, restart_dist, &
         init_dist_fermi, init_dist_lorentz, init_dist_gaussian
      use pert_param, only: boltz_nstep, output_nstep, time_step, &
         boltz_init_dist, boltz_init_e0, boltz_init_smear, hole
      use boltz_grid, only: grid
      use hdf5_utils
      implicit none
      character(*), intent(in) :: filename
      type(grid), intent(in) :: kgrid
      integer(HID_T), intent(out) :: fid, gid
      real(dp), intent(out) :: dist( kgrid%numb, kgrid%nk )
      !local variables
      logical :: has_file
      real(dp) :: t_step_fs
      character(len=120) :: group_name, msg
      integer :: nrun, current_run
      character(len=6), external :: int_to_char
   
      gid = 0;  fid = 0  ! init
      !restart dynamics simulation from previous run
      if(trim(boltz_init_dist) .eq. 'restart') then
         inquire(file=trim(filename), exist=has_file)
         if(.not. has_file) call errore('cdyna_setup','missing '// trim(filename), 1)
      
         if(ionode) then
            call hdf_open_file(fid, trim(filename), status='OLD', action='READWRITE')
            call hdf_read_dataset(fid, 'num_runs', nrun)
            call restart_dist(kgrid, fid, nrun, dist)
            current_run = nrun + 1
            !
            call hdf_update_dataset(fid, 'num_runs', current_run)
         endif
         call mp_bcast(dist, ionode_id, inter_pool_comm)
      else
         select case (boltz_init_dist)
         case('fermi')
            call init_dist_fermi(kgrid, dist, boltz_init_e0, boltz_init_smear)
         case ('lorentz')
            call init_dist_lorentz(kgrid, dist, boltz_init_e0, boltz_init_smear, hole)
         case('gaussian')
            call init_dist_gaussian(kgrid, dist, boltz_init_e0, boltz_init_smear, hole)
         case default
            msg = "invalid boltz_init_dist, valid options are 'restart', 'fermi', 'lorentz', 'gaussian'."
            call errore('carrier_dynamics_run', trim(msg), 1)
         end select
         
         if(ionode) then
            current_run = 1
            call hdf_open_file(fid, trim(filename), status='NEW')
            call hdf_write_dataset(fid, 'num_runs', current_run)
            !write kg%enk(:,:)
            call hdf_write_dataset(fid, 'band_structure_ryd', kgrid%enk)
            call hdf_write_attribute(fid, 'band_structure_ryd', 'ryd2ev', ryd2ev)
         endif
      endif
      
      !record status of the current run
      if(ionode) then
         group_name = "dynamics_run_" // trim( int_to_char(current_run) )
         call hdf_create_group(fid, trim(group_name))
         call hdf_open_group(fid, trim(group_name), gid)
         !
         t_step_fs = time_step * output_nstep * (timeunit*1.0E15_dp)
         call hdf_write_dataset(gid, 'time_step_fs', t_step_fs)
         call hdf_write_dataset(gid, 'num_steps',  0)
         !
         if(current_run .eq. 1)  call output_dist(kgrid, gid, dist, 0)
      endif
      !
      call mp_barrier(inter_pool_comm)
   end subroutine cdyna_setup
   !
end subroutine carrier_dynamics_run

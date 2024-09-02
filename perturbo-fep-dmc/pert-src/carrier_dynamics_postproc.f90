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

subroutine carrier_dynamics_postproc()
   use pert_const, only: dp, ryd2ev, bohr2ang
   use pert_data,  only: epwan_fid, kc_dim
   use band_structure, only: electron_wann, init_electron_wann
   use qe_mpi_mod, only: ionode, stdout, mp_bcast, ionode_id, inter_pool_comm
   use boltz_grid, only: grid, init_boltz_grid
   use pert_param, only: boltz_de, boltz_emax, boltz_emin, band_min, band_max, prefix, hole
   use boltz_dynamics_mod, only: read_dist, calc_carrier_population
   use pert_output,only: progressbar_init, progressbar
   use pert_utils, only: find_free_unit
   use hdf5_utils
   implicit none
   type(grid) :: kg
   logical :: has_file
   integer :: nestep, i, punit, nrun, tot_nt, it, irun, istep, nstep
   real(dp) :: emin, emax, t, t_step, cnum
   character(len=120) :: fname, group_name, pname, dset_name
   integer(HID_T) :: file_id, group_id, popu_id, popu_gid
   real(dp), allocatable :: ene(:), ene_ev(:), popu(:), dist(:,:), tt(:)
   type(electron_wann) :: elec
   character(len=6), external :: int_to_char
   
   call init_electron_wann(epwan_fid, kc_dim, elec)
   !set up k-grid
   call init_boltz_grid(kg, elec, band_min, band_max)

   !setup up energy windows
   emin = max(boltz_emin, minval(kg%enk))
   emax = min(boltz_emax, maxval(kg%enk))
   !if(ionode) write(stdout,'(2x,a,2(1x,f12.6))') 'Energy window (eV):', emin*ryd2ev, emax*ryd2ev
   if(emin > emax) call errore('transport','illegal energy window',1)
   
   !setup energy grid: emin - boltz_de : emax + boltz_de
   nestep = int( floor((emax-emin)/boltz_de) ) + 3
   allocate(ene(nestep), ene_ev(nestep), popu(nestep), dist(kg%numb, kg%nk))
   do i = 1, nestep
      ene(i) = emin + (i-2)*boltz_de
   enddo
   
   fname = trim(prefix)//"_cdyna.h5"
   pname = trim(prefix)//'_cdyna.dat'
   !
   inquire(file=trim(fname), exist=has_file)
   if(.not. has_file) call errore('carrier_dynamics_postproc','missing '// trim(fname), 1)
   !open files
   if(ionode) then
      call hdf_open_file(file_id, trim(fname), status='OLD', action='READ')
      call hdf_read_dataset(file_id, 'num_runs', nrun)

      !find out total number of snapshots
      it = 1
      do irun = 1, nrun
         group_name = "dynamics_run_" // trim( int_to_char(irun) )
         call hdf_open_group(file_id, trim(group_name), group_id)
         call hdf_read_dataset(group_id, 'num_steps', nstep)
         call hdf_close_group(group_id)
         it = it + nstep 
      enddo
      tot_nt = it
      allocate( tt( 0:tot_nt ) )
      tt = 0.0_dp

      punit = find_free_unit()
      open(punit, file=trim(pname), form='formatted', status='unknown',action='write')
      write(punit, 1001)
         
      !open output hdf5 file
      call hdf_open_file(popu_id, trim(prefix)//"_popu.h5", status='NEW')
      ene_ev(:) = ene(:) * ryd2ev
      call hdf_write_dataset( popu_id, 'energy_grid_ev', ene_ev)
      call hdf_create_group(popu_id, 'energy_distribution')
      call hdf_open_group(popu_id, 'energy_distribution', popu_gid)
   endif
   call mp_bcast(nrun, ionode_id, inter_pool_comm)
   call mp_bcast(tot_nt, ionode_id, inter_pool_comm)
   
   if(ionode) call progressbar_init('Carrier Dynamics - postproc:')

   it = 0
   t = 0.0_dp
   do irun = 1,  nrun
      if(ionode) then
         group_name = "dynamics_run_" // trim( int_to_char(irun) )
         call hdf_open_group(file_id, trim(group_name), group_id)
         call hdf_read_dataset(group_id, 'num_steps', nstep)
         call hdf_read_dataset(group_id, 'time_step_fs', t_step)
      endif
      call mp_bcast(nstep, ionode_id, inter_pool_comm)
      call mp_bcast(t_step, ionode_id, inter_pool_comm)

      do istep = 0, nstep
         if(irun > 1 .and. istep .eq. 0) cycle
         !skip the first snap in the first run as it's the inital f(t=0)
         if(irun > 1 .or.  istep .ne. 0) then
            t = t + t_step
            it = it + 1
         endif

         if(ionode) call read_dist(kg, group_id, dist, istep)
         call mp_bcast(dist, ionode_id, inter_pool_comm)
         call calc_carrier_population(kg, dist, ene, popu, hole)

         !output dist
         if(ionode) then
            cnum  = sum(popu)*boltz_de
            tt(it) = t
            write(punit, 1002) it, t, cnum
            
            dset_name = 'popu_t' // trim( int_to_char(it) )
            call hdf_write_dataset(popu_gid, trim(dset_name), popu)
            !
            call progressbar(it, tot_nt)
         endif
      enddo
      if(ionode) call hdf_close_group(group_id)
   enddo

   if(ionode) then
      call hdf_close_group(popu_gid)
      call hdf_write_dataset(popu_id, 'times_fs', tt)
      call hdf_close_file(popu_id)
      call hdf_close_file(file_id)
      close(punit)
      
      deallocate( tt )
   endif

   deallocate( ene, ene_ev, popu, dist )
   return
1001 format(1x,'#.record         time (fs)       #.carrier/u.c. ')
1002 format(1x, i8, 2x, f16.4, 2x,  es23.16)
end subroutine carrier_dynamics_postproc

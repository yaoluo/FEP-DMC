!===============================================================================
! Copyright (C) 2016-2020 Jin-Jian Zhou, Jinsoo Park, I-Te Lu, Marco Bernardi
!
! This program is distributed under the terms of the GNU General Public License.
! See the file `LICENSE' in the root directory of this distribution, or obtain 
! a copy of the License at <https://www.gnu.org/licenses/gpl-3.0.txt>.
!
! Author: jjzhou <jjchou.comphy@gmail.com>
! Comment:
!    load data from qe2pert output
!
! Maintenance:
!===============================================================================

subroutine load_data_pert(file_id)
   use pert_const, only: dp, twopi, bohr2ang
   use pert_param, only: prefix, band_min, band_max, spinor
   use qe_mpi_mod, only: ionode, mp_bcast, inter_pool_comm, ionode_id, ionode, stdout
   use pert_data,  only:nat, tau, mass, volume, alat, at, bg, zstar, epsil, tpiba, nsym, &
      symop, num_wann, kc_dim, qc_dim, wannier_center, wannier_center_cryst, polar_alpha, &
      thickness_2d, system_2d, loto_alpha, lpolar
   use pert_utils, only: find_free_unit
   use hdf5_utils
   implicit none
   integer(HID_T), intent(in) :: file_id
   !
   integer :: spinor_int, s2d, polar_int
   integer(HID_T) :: group_id

   symop = 0
   if(ionode) then
      call hdf_open_group(file_id, 'basic_data', group_id)
      !
      call hdf_read_dataset(group_id, 'num_wann', num_wann)
      call hdf_read_dataset(group_id, 'nat', nat)
      call hdf_read_dataset(group_id, 'kc_dim', kc_dim)
      call hdf_read_dataset(group_id, 'qc_dim', qc_dim)
      call hdf_read_dataset(group_id, 'volume', volume)
      call hdf_read_dataset(group_id, 'epsil', epsil)
      call hdf_read_dataset(group_id, 'alat', alat)
      call hdf_read_dataset(group_id, 'at', at)
      call hdf_read_dataset(group_id, 'bg', bg)
      call hdf_read_dataset(group_id, 'nsym', nsym)
      call hdf_read_dataset(group_id, 'symop', symop(:,:,1:nsym))

      if( hdf_exists(group_id, 'spinor') ) then
         call hdf_read_dataset(group_id, 'spinor', spinor_int)
         spinor = (spinor_int > 0)
      else
         spinor = .false.
         write(stdout,'(5x,a)') "N.B.: no 'spinor' found, set it to .false."
      endif
      
      if( hdf_exists(group_id, 'polar_alpha') ) then
         call hdf_read_dataset(group_id, 'polar_alpha', polar_alpha)
      else
         polar_alpha = 1.0_dp
         write(stdout,'(5x,a)') "N.B.: no 'polar_alpha' found, set it to 1.0"
      endif
      
      if( hdf_exists(group_id, 'loto_alpha') ) then
         call hdf_read_dataset(group_id, 'loto_alpha', loto_alpha)
      else
         loto_alpha = 1.0_dp  ! default value used in PHonon/rgd_blk
      endif
      
      if( hdf_exists(group_id, 'system_2d') ) then
         call hdf_read_dataset(group_id, 'system_2d', s2d)
         system_2d = (s2d > 0)
         if(system_2d) then
            write(stdout, '(5x,a)') "N.B.: current system is 2D."
            if(any(abs(at(1:2, 3)) > 1.0E-6_dp) .or. any(abs(at(3, 1:2)) > 1.0E-6_dp)) &
               call errore('load_data_pert','2D system requires the third lattice vector along z-dirction',1)
         endif
      else
         system_2d = .false.
         write(stdout,'(5x,a)') "N.B.: no 'system_2d' found, assume 3D system."
      endif
      
      if( system_2d ) then
         if( hdf_exists(group_id, 'thickness_2d') ) then
            call hdf_read_dataset(group_id, 'thickness_2d', thickness_2d)
            if(thickness_2d <= 0.0d0) &
               call errore('load_data_pert','thickness_2d is negative (or too small).',1)
         else
            thickness_2d = 6.0_dp / bohr2ang  !set to 6 angstrom by default
            write(stdout,'(5x,a)') "N.B.: no 'thickness_2d' found, set it to 6 A."
         endif
      else
         thickness_2d = -1.0_dp !set it to a negative value, menas a 3D system.
      endif
   endif
   
   call mp_bcast(spinor, ionode_id, inter_pool_comm)
   call mp_bcast(volume, ionode_id, inter_pool_comm)
   call mp_bcast(alat, ionode_id, inter_pool_comm)
   call mp_bcast(at, ionode_id, inter_pool_comm)
   call mp_bcast(bg, ionode_id, inter_pool_comm)
   call mp_bcast(polar_alpha, ionode_id, inter_pool_comm)
   call mp_bcast(loto_alpha, ionode_id, inter_pool_comm)
   call mp_bcast(system_2d, ionode_id, inter_pool_comm)
   call mp_bcast(thickness_2d, ionode_id, inter_pool_comm)
   call mp_bcast(epsil, ionode_id, inter_pool_comm)
   call mp_bcast(nat, ionode_id, inter_pool_comm)
   call mp_bcast(num_wann, ionode_id, inter_pool_comm)
   call mp_bcast(kc_dim, ionode_id, inter_pool_comm)
   call mp_bcast(qc_dim, ionode_id, inter_pool_comm)
   call mp_bcast(nsym, ionode_id, inter_pool_comm)
   call mp_bcast(symop, ionode_id, inter_pool_comm)

   allocate(mass(nat), zstar(3,3,nat), tau(3,nat))
   allocate( wannier_center(3, num_wann), wannier_center_cryst(3, num_wann) )

   if(ionode) then
      call hdf_read_dataset(group_id, 'mass', mass(1:nat))
      call hdf_read_dataset(group_id, 'zstar', zstar(:,:,1:nat))
      call hdf_read_dataset(group_id, 'tau', tau(:, 1:nat) )
      call hdf_read_dataset(group_id, 'wannier_center', wannier_center(:, 1:num_wann))
      call hdf_read_dataset(group_id, 'wannier_center_cryst', wannier_center_cryst(:,1:num_wann))

      if( hdf_exists(group_id, 'lpolar') ) then
         call hdf_read_dataset(group_id, 'lpolar', polar_int)
         lpolar = (polar_int > 0)
      else
         lpolar = merge(.true., .false., sum(abs(zstar)) > 1.0E-1_dp)
      endif
      !
      call hdf_close_group(group_id)
   endif

   call mp_bcast(mass, ionode_id, inter_pool_comm)
   call mp_bcast(tau, ionode_id, inter_pool_comm)
   call mp_bcast(zstar, ionode_id, inter_pool_comm)
   call mp_bcast(wannier_center, ionode_id, inter_pool_comm)
   call mp_bcast(wannier_center_cryst, ionode_id, inter_pool_comm)
   call mp_bcast(lpolar, ionode_id, inter_pool_comm)
   
   tpiba = twopi/alat
   !band range should not exceed nband
   if(band_min > num_wann) band_min = 1
   if(band_max > num_wann) band_max = num_wann
   if(band_min > num_wann .or. band_max > num_wann) write(stdout, '(5x, a)') &
         "Warning: band_min(or band_max) exceed total number of bands, reset to:"
   write(stdout, '(3x, 2(2x,a,i4))') "band_min: ", band_min, "band_max: ", band_max

end subroutine load_data_pert

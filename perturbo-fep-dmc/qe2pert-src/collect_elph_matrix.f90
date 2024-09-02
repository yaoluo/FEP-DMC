!===============================================================================
! Copyright (C) 2016-2020 Jin-Jian Zhou, Jinsoo Park, I-Te Lu, Marco Bernardi
!
! This program is distributed under the terms of the GNU General Public License.
! See the file `LICENSE' in the root directory of this distribution, or obtain 
! a copy of the License at <https://www.gnu.org/licenses/gpl-3.0.txt>.
!
! Author: jjzhou <jjchou.comphy@gmail.com>
! Comment:
!  collect g(k_c, q_c) from different processes and save to a single HDF5 file.
!
! Maintenance:
!===============================================================================

subroutine collect_elph_matrix()
   use kinds, only: dp
   use ions_base, only: nat 
   use lattice_data, only: numq
   use electronic_data, only: num_kpts, nb_sub
   use qe_mpi_mod, only: mp_barrier, mp_gatherv_cv4d, stdout, ionode, ionode_id, &
      inter_pool_comm, npool, my_pool_id, mp_split_pools, mp_sum
   use input_param, only: prefix, tmp_dir, lwannier
   use buffers, only: open_buffer, get_buffer, close_buffer
   use epwan_hdf5_io, only: write_ephmat_part
   use hdf5_utils
   implicit none
   !
   logical :: exst 
   integer(HID_T) :: file_id
   character(len=120) :: fname
   integer, allocatable :: recvcount(:), displs(:)
   complex(dp), allocatable :: epmat_loc(:,:,:,:), epmat(:,:,:,:)
   integer :: eph_nword, iun_eph, irec, iq_st, iq_end, nq_loc, i, ip, ia, j
   !
   integer, external :: find_free_unit

   if(ionode) write(stdout, '(5x,a)') "collecting electron-phonon matrix..."
   call mp_barrier(inter_pool_comm)

   eph_nword = nb_sub * nb_sub * num_kpts
   !open prefix.ephmat file
   !io_level = 1, read eph matrix from disk
   iun_eph = find_free_unit()
   call open_buffer(iun_eph, 'ephmat', eph_nword, 1, exst)
   !
   call mp_split_pools(numq, iq_st, iq_end, nq_loc)
   !
   !set recvcount and displs for mp_gather
   allocate( recvcount(npool), displs(npool) )
   displs = 0
   recvcount = 0
   recvcount( my_pool_id+1 ) = nq_loc * eph_nword
   call mp_sum(recvcount, inter_pool_comm)
   do i = 2, npool
      displs(i) = sum( recvcount( 1:(i-1) ) )
   enddo 
   !
   !all epmat will be gathered into ionode
   if(ionode) then
      !open hdf5 file
      fname = trim(tmp_dir) // trim(prefix) // "_elph.h5"
      call hdf_open_file(file_id, trim(fname), status='NEW')
      ip = merge(1, -1, lwannier)
      call hdf_write_dataset(file_id, "wannier_gauge", ip)
      !
      allocate( epmat(nb_sub, nb_sub, num_kpts, numq) )
   else
      allocate( epmat(1,1,1,1) )
   endif
   !allocate workspace
   allocate( epmat_loc(nb_sub, nb_sub, num_kpts, nq_loc) )
   !
   !gather elph matrix elements on coarse grid
   do ia = 1, nat
   do j = 1, 3
      ip = (ia - 1) * 3 + j
      ! load data
      !if(ionode) epmat = (0.0_dp, 0.0_dp)
      do i = 1, nq_loc
         irec = (i-1) * 3 * nat + ip
         call get_buffer(epmat_loc(1,1,1,i), eph_nword, iun_eph, irec)
      enddo
      ! gather data into ionode
      call mp_gatherv_cv4d(epmat_loc, epmat, recvcount, displs, ionode_id, inter_pool_comm)
      !
      !output to hdf5 file
      if(ionode) call write_ephmat_part(file_id, ia, j, nb_sub, num_kpts, numq, epmat)
      !call mp_barrier(inter_pool_comm)
   enddo; enddo
   !release space
   deallocate(epmat, epmat_loc, recvcount, displs)
   if(ionode) call hdf_close_file(file_id)
   call mp_barrier(inter_pool_comm)
   
   !if everything goes well, delete these temporary files.
   call close_buffer(iun_eph, 'DELETE')
   return
end subroutine collect_elph_matrix

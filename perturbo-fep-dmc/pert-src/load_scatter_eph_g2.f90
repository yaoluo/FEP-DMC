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

subroutine load_scatter_eph_g2(nc_loc, kq_pair, num_scat, scat)
   use hdf5_utils
   use pert_param, only: prefix, tmp_dir
   implicit none
   integer, intent(in) :: nc_loc, num_scat
   type(channel), intent(inout) :: kq_pair(nc_loc)
   type(t_scatt), intent(out) :: scat(num_scat)
   !
   integer :: n, i, ik, nkq, j, tot_chnl, nchl, it, ib
   integer, allocatable :: ist(:), bnd_idx(:)
   real(dp), allocatable :: eph_g2(:)
   !
   logical :: has_file
   integer(HID_T) :: file_id
   character(len=120) :: fname, dset_name
   character(len=6), external :: int_to_char
   
   call start_clock("load_eph_g2")
   !open hdf5 file
   fname = trim(tmp_dir) // trim(prefix) // "_eph_g2_p" 
   fname = trim(fname) // trim( int_to_char(my_pool_id+1) ) // ".h5"
   ! check if file exist
   inquire(file=trim(fname), exist=has_file)
   if(.not. has_file) call errore('load_scatter_eph_g2', "Missing file: "//trim(fname), 1)
   ! open file
   call hdf_open_file(file_id, trim(fname), status='OLD', action='READ')

   n = 0
   do i = 1, nk_loc
      ik = kq_pair(i)%ik
      nkq = kq_pair(i)%nkq
      
      allocate( ist( nkq ) )
      ist(1) = 0 
      do j = 1, nkq-1
         ist(j+1) = ist(j) + kq_pair(i)%nchl(j)
      enddo

      tot_chnl = sum( kq_pair(i)%nchl(:) )
      allocate( bnd_idx(tot_chnl), eph_g2(tot_chnl) )
      bnd_idx = 0;   eph_g2 = 0.0_dp
      
      !read bnd_idx and eph_g2 from hdf5 file
      dset_name = "bands_index_" // trim( int_to_char(i) )
      call hdf_read_dataset(file_id, trim(dset_name), bnd_idx)
      dset_name = "eph_g2_" // trim( int_to_char(i) )
      call hdf_read_dataset(file_id, trim(dset_name), eph_g2)

      do j = 1, nkq
         n = n + 1
         !
         scat(n)%ik = ik
         scat(n)%iq = kq_pair(i)%iq(j)
         scat(n)%ikq = kq_pair(i)%ikq(j)
         !
         nchl = kq_pair(i)%nchl(j)
         scat(n)%nchl = nchl
         allocate( scat(n)%bands_index(nchl), scat(n)%eph_g2(nchl) )

         it = ist(j)
         do ib = 1, nchl
            it = it + 1
            scat(n)%bands_index(ib) = bnd_idx(it)
            scat(n)%eph_g2(ib) = eph_g2(it)
         enddo
      enddo
      deallocate(ist, bnd_idx, eph_g2)
      !
      ! they are no longer needed.
      deallocate( kq_pair(i)%iq, kq_pair(i)%ikq, kq_pair(i)%nchl )
   enddo
   !
   call hdf_close_file(file_id)
   call mp_barrier( inter_pool_comm )
   call stop_clock("load_eph_g2")
end subroutine load_scatter_eph_g2

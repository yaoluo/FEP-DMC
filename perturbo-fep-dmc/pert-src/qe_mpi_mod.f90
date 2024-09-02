!===============================================================================
! Copyright (C) 2016-2020 Jin-Jian Zhou, Jinsoo Park, I-Te Lu, Marco Bernardi
!
! This program is distributed under the terms of the GNU General Public License.
! See the file `LICENSE' in the root directory of this distribution, or obtain 
! a copy of the License at <https://www.gnu.org/licenses/gpl-3.0.txt>.
!
! Author: jjzhou <jjchou.comphy@gmail.com>
! Comment:
!   mpi wrapper from QE/Modules
!
! Maintenance:
!===============================================================================

module qe_mpi_mod
   use parallel_include
   use kinds, only: dp
   use mp_world,  ONLY : world_comm
   use mp_global, only: mp_startup, mp_global_end
   use mp_pools,  only: npool, my_pool_id, inter_pool_comm
   use io_global, only: stdout, ionode, ionode_id, meta_ionode_id, meta_ionode
   use mp, only: mp_bcast, mp_barrier, mp_sum, mp_stop, mp_end, &
      mp_circular_shift_left, mp_gather
   implicit none
   public :: mp_split_pools, distribute_points, &
      mp_circular_shift, mp_scatter, mp_gatherv_cv4d
   
   !adapted from QE/Module/mp.f90/mp_circular_shift_left
   interface mp_circular_shift
      module procedure &
         mp_circular_shift_i0,  &
         mp_circular_shift_r2d, &
         mp_circular_shift_c2d, &
         mp_circular_shift_c3d
   end interface


   interface mp_scatter
      module procedure mp_scatterv_cv4d
   end interface
   
contains

subroutine mp_split_pools(num, ist, iend, n_pool)
   implicit none
   integer, intent(in) :: num
   integer, intent(out) :: ist, iend
   integer, intent(out), optional :: n_pool
#ifdef __MPI
   !for parallel version
   integer :: base, rest
   base = num / npool
   rest = mod(num, npool)
   ! my_pool_id range from 0, 1, 2, ..., npool-1
   if(my_pool_id < rest) then
      ist = my_pool_id * (base + 1) + 1
      iend = ist + base
   else
      ist = my_pool_id*base + merge(rest, 0, base>0) + 1
      iend = ist + base - 1
   end if
#else
   !for serial version
   ist = 1;  iend = num
#endif
   if(present(n_pool)) n_pool = iend - ist + 1
end subroutine mp_split_pools


subroutine distribute_points(npts, npts_loc, pts_loc, ldist)
   implicit none
   integer, intent(in) :: npts
   logical, intent(in) :: ldist
   integer, intent(out):: npts_loc
   integer, pointer    :: pts_loc(:)
   !
   integer :: i, ik

   if( associated(pts_loc) )  deallocate( pts_loc )
   if(ldist) then
      !interleave distribution of kpoints over pools, for better load balance
      npts_loc = npts/npool + merge(1, 0, my_pool_id < mod(npts, npool))
      allocate( pts_loc( npts_loc ) )
      !
      i = 0
      do ik = my_pool_id+1, npts, npool
         i = i + 1
         pts_loc(i) = ik
      enddo
   else
      !no distribution
      npts_loc = npts
      allocate( pts_loc(npts) )

      do ik = 1, npts
         pts_loc(ik) = ik
      enddo
   endif
end subroutine distribute_points


subroutine mp_circular_shift_i0(buf, itag, gid)
   implicit none
   integer :: buf
   integer, intent(in) :: itag
   integer, intent(in) :: gid

   call mp_circular_shift_left(buf, itag, gid)
end subroutine mp_circular_shift_i0

subroutine mp_circular_shift_r2d(buf, itag, gid)
   implicit none
   real(dp) :: buf(:,:)
   integer, intent(in) :: itag
   integer, intent(in) :: gid

   call mp_circular_shift_left(buf, itag, gid)
end subroutine mp_circular_shift_r2d

subroutine mp_circular_shift_c2d(buf, itag, gid)
   implicit none
   complex(dp) :: buf(:,:)
   integer, intent(in) :: itag
   integer, intent(in) :: gid

   call mp_circular_shift_left(buf, itag, gid)
end subroutine mp_circular_shift_c2d

!adapted from mp_circular_shift_left_c2d
subroutine mp_circular_shift_c3d(buf, itag, gid)
   implicit none
   complex(dp) :: buf(:,:,:)
   integer, intent(in) :: itag
   integer, intent(in) :: gid
   INTEGER :: group, ierr, npe, sour, dest, mype

#if defined (__MPI)
   INTEGER :: istatus( mpi_status_size )
   !
   group = gid
   !
   CALL mpi_comm_size( group, npe, ierr )
   IF (ierr/=0) CALL mp_stop( 8100 )
   CALL mpi_comm_rank( group, mype, ierr )
   IF (ierr/=0) CALL mp_stop( 8101 )
   !
   sour = mype + 1
   IF( sour == npe ) sour = 0
   dest = mype - 1
   IF( dest == -1 ) dest = npe - 1
   !
   CALL MPI_Sendrecv_replace( buf, SIZE(buf), MPI_DOUBLE_COMPLEX, &
        dest, itag, sour, itag, group, istatus, ierr)
   !
   IF (ierr/=0) CALL mp_stop( 8102 )
   !
#else
   ! do nothing
#endif
   RETURN
end subroutine mp_circular_shift_c3d


!based on QE/mp.f90/mp_gatherv_cv
subroutine mp_scatterv_cv4d(alldata, mydata, sendcount, displs, root, gid)
   IMPLICIT NONE
   COMPLEX(DP) :: alldata(:,:,:,:)
   COMPLEX(DP) :: mydata(:,:,:,:)
   INTEGER, INTENT(IN) :: sendcount(:), displs(:), root
   INTEGER, INTENT(IN) :: gid
   !
   INTEGER :: group
   INTEGER :: ierr, npe, myid

#if defined (__MPI)
   group = gid
   CALL mpi_comm_size( group, npe, ierr )
   IF (ierr/=0) CALL mp_stop( 8069 )
   CALL mpi_comm_rank( group, myid, ierr )
   IF (ierr/=0) CALL mp_stop( 8070 )
   !
   IF ( SIZE( sendcount ) < npe .OR. SIZE( displs ) < npe ) CALL mp_stop( 8071 )
   IF ( myid == root ) THEN
      IF ( SIZE( alldata ) < displs( npe ) + sendcount( npe ) ) CALL mp_stop( 8072 )
   END IF
   IF ( SIZE( mydata ) < sendcount( myid + 1 ) ) CALL mp_stop( 8073 )
   !
   CALL MPI_SCATTERV(alldata, sendcount, displs, MPI_DOUBLE_COMPLEX, &
         mydata, sendcount( myid + 1 ), MPI_DOUBLE_COMPLEX, root, group, ierr )
   IF (ierr/=0) CALL mp_stop( 8074 )
#else
   IF ( SIZE( alldata ) < sendcount( 1 ) ) CALL mp_stop( 8075 )
   IF ( SIZE( mydata  ) < sendcount( 1 ) ) CALL mp_stop( 8076 )
   ! by jjzhou
   IF ( SIZE( alldata ) .ne. SIZE( mydata ) ) call mp_stop( 8076 )
   !
   !mydata( 1:sendcount( 1 ) ) = alldata( 1:sendcount( 1 ) )
   mydata(:,:,:,:) = alldata(:,:,:,:)
#endif
   RETURN
end subroutine


SUBROUTINE mp_gatherv_cv4d( mydata, alldata, recvcount, displs, root, gid)
  IMPLICIT NONE
  COMPLEX(DP) :: mydata(:,:,:,:)
  COMPLEX(DP) :: alldata(:,:,:,:)
  INTEGER, INTENT(IN) :: recvcount(:), displs(:), root
  INTEGER, INTENT(IN) :: gid
  INTEGER :: group
  INTEGER :: ierr, npe, myid

#if defined (__MPI)
  group = gid
  CALL mpi_comm_size( group, npe, ierr )
  IF (ierr/=0) CALL mp_stop( 8069 )
  CALL mpi_comm_rank( group, myid, ierr )
  IF (ierr/=0) CALL mp_stop( 8070 )
  !
  IF ( SIZE( recvcount ) < npe .OR. SIZE( displs ) < npe ) CALL mp_stop( 8071 )
  IF ( myid == root ) THEN
     IF ( SIZE( alldata ) < displs( npe ) + recvcount( npe ) ) CALL mp_stop( 8072 )
  END IF
  IF ( SIZE( mydata ) < recvcount( myid + 1 ) ) CALL mp_stop( 8073 )
  !
  CALL MPI_GATHERV( mydata, recvcount( myid + 1 ), MPI_DOUBLE_COMPLEX, &
                   alldata, recvcount, displs, MPI_DOUBLE_COMPLEX, root, group, ierr )
  IF (ierr/=0) CALL mp_stop( 8074 )
#else
  IF ( SIZE( alldata ) < recvcount( 1 ) ) CALL mp_stop( 8075 )
  IF ( SIZE( mydata  ) < recvcount( 1 ) ) CALL mp_stop( 8076 )
   ! by jjzhou
   IF ( SIZE( alldata ) .ne. SIZE( mydata ) ) call mp_stop( 8076 )
  !
  !alldata( 1:recvcount( 1 ) ) = mydata( 1:recvcount( 1 ) )
   alldata(:,:,:,:) = mydata(:,:,:,:)
#endif
  RETURN
END SUBROUTINE mp_gatherv_cv4d
   
end module qe_mpi_mod

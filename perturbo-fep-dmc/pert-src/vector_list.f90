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

module vector_list
   use pert_const, only: dp
   implicit none

   type, public :: vlist
      integer :: nvec
      real(dp), pointer :: vec(:,:)=>null()
      real(dp), pointer :: weight(:)=>null()
   end type vlist

   public :: load_vector_list, rand_vector_list

contains

subroutine load_vector_list(fil, vl)
   use pert_utils, only: find_free_unit
   use qe_mpi_mod, only: stdout, ionode, mp_bcast, ionode_id, inter_pool_comm
   implicit none
   character(len=80), intent(in) :: fil
   type(vlist), intent(out) :: vl
   integer :: iunit, nline, i, ios
   !
   type(vlist) :: vl_tmp

   if (ionode) then
      write(stdout,'(5x,a)') "Reading list from file: "//trim(fil)
      iunit = find_free_unit()
      open(iunit,file=trim(fil),status='old',form='formatted',err=100,iostat=ios)
100   call errore('load_vector_list','opening file '//trim(fil),abs(ios))
      read(iunit, *) nline

      vl_tmp%nvec = nline
      allocate( vl_tmp%vec(3,nline), vl_tmp%weight(nline) )
      do i = 1, nline
         read(iunit,*) vl_tmp%vec(1:3,i), vl_tmp%weight(i)
      enddo
      close(iunit)

      !may need to expand the list
      if( vl_tmp%nvec > 1 .and. vl_tmp%weight(1) > 1.0_dp ) then
         call expand_vector_list(vl_tmp, vl)
         deallocate( vl_tmp%vec, vl_tmp%weight )
      else
         vl%nvec = vl_tmp%nvec
         vl%vec => vl_tmp%vec
         vl%weight => vl_tmp%weight
      endif
   endif
   call mp_bcast(vl%nvec, ionode_id, inter_pool_comm)
   if(.not. associated(vl%vec)) allocate(vl%vec(3,vl%nvec))
   if(.not. associated(vl%weight) ) allocate(vl%weight(vl%nvec))
   call mp_bcast(vl%vec, ionode_id, inter_pool_comm)
   call mp_bcast(vl%weight, ionode_id, inter_pool_comm)
end subroutine load_vector_list

subroutine rand_vector_list(vl, nvec, rand_dist, rand_param)
   use pert_utils, only: random_uniform, random_cauchy
   use qe_mpi_mod, only: stdout, ionode, mp_bcast, ionode_id, inter_pool_comm
   implicit none
   type(vlist), intent(out) :: vl
   integer, intent(in) :: nvec
   character(len=80), intent(in) :: rand_dist
   real(dp), intent(in), optional :: rand_param

   if(ionode) then
      vl%nvec = nvec
      allocate( vl%vec(3,nvec), vl%weight(nvec) )
      if(trim(rand_dist) .eq. 'cauchy') then
         if(present(rand_param)) then
            call random_cauchy(vl%vec, vl%weight, rand_param)
         else
            call random_cauchy(vl%vec, vl%weight)
         endif
         write(stdout,'(5x, a,E13.6)') &
            "Total weight of the Cauchy sampling: ", sum(vl%weight)
      elseif( trim(rand_dist) .eq. 'uniform' ) then
         call random_uniform(vl%vec, vl%weight)
      else
         call errore('rand_vector_list','invalid sampling mode: '//trim(rand_dist),1)
      endif
   endif
   call mp_bcast(vl%nvec, ionode_id, inter_pool_comm)
   if(.not. associated(vl%vec)) allocate(vl%vec(3,vl%nvec))
   if(.not. associated(vl%weight) ) allocate(vl%weight(vl%nvec))
   call mp_bcast(vl%vec, ionode_id, inter_pool_comm)
   call mp_bcast(vl%weight, ionode_id, inter_pool_comm)
end subroutine rand_vector_list

subroutine expand_vector_list(vl_in, vl_out)
   type(vlist), intent(in) :: vl_in
   type(vlist), intent(out) :: vl_out
   !
   integer :: tot_vec, nv, n, i, j
   real(dp) :: st(3), ed(3), dvec(3)

   if(associated(vl_out%vec) .or. associated(vl_out%weight)) &
      call errore('expand_vector_list','vl_out is not empty!',1)
   if( vl_in%nvec < 2 ) &
      call errore('expand_vector_list','cannot generate k-path with only one k-point!',1)

   tot_vec = vl_in%nvec
   do i = 1, vl_in%nvec - 1
      tot_vec = tot_vec + floor( vl_in%weight(i) )
   enddo
   vl_out%nvec = tot_vec
   allocate( vl_out%vec(3, tot_vec), vl_out%weight(tot_vec) )

   !generate k-list
   nv = 0
   do i = 1,  vl_in%nvec - 1
      st(:) = vl_in%vec(:, i)
      ed(:) = vl_in%vec(:, i+1)
      n = 1 + floor( vl_in%weight(i) )
      dvec(:) = (ed(:) - st(:)) / real(n, dp)

      do j = 1, n
         vl_out%vec(:, nv+j) = st(:) + dvec(:) * real(j-1, dp)
      enddo
      nv = nv + n
   enddo
   vl_out%vec(:, tot_vec) = vl_in%vec(:, vl_in%nvec)
   vl_out%weight(:) = 1.0_dp / real(tot_vec, dp)
end subroutine expand_vector_list

end module vector_list

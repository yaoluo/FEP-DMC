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

subroutine setup_kq_pair(kg, el, ph, wcut, e_thr, nk_loc, kq_pair, qg)
   use pert_utils, only: match_table
   use boltz_utils, only: kpts_minus, num2kpt, band2idx, kpts_plus
   implicit none
   type(grid), intent(in) :: kg
   type(electron_wann), intent(in) :: el
   type(lattice_ifc), intent(in) :: ph
   real(dp), intent(in) :: wcut, e_thr
   !
   integer, intent(in) :: nk_loc
   type(channel), intent(inout) :: kq_pair(:)
   !
   type(phonon_grid), intent(inout) :: qg
   !
   !local variables
   integer :: nkpt, nmod, maxq, nkq, nq_tmp, i, ik, ikq, j, qidx, iq, k
   logical :: ltable(kg%numb, kg%numb, ph%nm)
   integer, allocatable :: qcollect(:)
   real(dp) :: xq(3)
   !
   integer,  pointer :: ql_tmp(:), q_check(:)
   real(dp), pointer :: qw_tmp(:,:)
   type(channel), pointer :: kq_tmp(:)

   call start_clock('set_kq_pair')
   !init
   nkpt = kg%nk
   nmod = ph%nm
   
   allocate( kq_tmp(nk_loc) )
   !find iq
   !maxq = qg%ndim(1)*qg%ndim(2)*qg%ndim(3) - 1
   maxq = qg%ndim(3) * (1 +  qg%ndim(2) * (1 + qg%ndim(1)) ) / 2
   allocate(qcollect(0:maxq))
   qcollect = 0
   !(k, k') and (k', k) are connected, only one of them is stored.
   !to balance load, each kpoint is associated with equal number of k+q. 
   !so for each ik in grid%kpt, the (k+q) are choose as ik+1, ik+2, ...., ik+(nkpt-1)/2
   !  if (ik + i) is larger than nkpt, restart from 1. for example, nkpt + 5 -> 5
   !
   !initialize and pre-allocate space
   do i = 1, nk_loc
      ik = kq_pair(i)%ik
      kq_tmp(i)%ik = ik
      !assign (k+q) evenly for each k, for load balance
      nkq = (nkpt-1)/2 + merge(1, 0, ((mod(nkpt-1,2)>0) .and. ik<=nkpt/2))
      kq_tmp(i)%nkq = nkq

      allocate( kq_tmp(i)%iq(nkq) )
      kq_tmp(i)%iq(:) = -1
   enddo
   !
   !collect q-points and store them in qg%list.
!$omp parallel do schedule(guided) default(shared) private(i, ik, nkq, j, ikq, qidx)
   do i = 1, nk_loc
      ik = kq_tmp(i)%ik
      nkq = kq_tmp(i)%nkq

      !loop over (k+q) points assigned to this current k-points
      do j = 1, nkq
         ikq = mod(ik+j-1, nkpt) + 1
         !store (k+q) - k in kq_pair 
         !(since q and -q has the same freqencies, we only store one of them)
         ! if xq = (ikq - iq) is not on the q-grid (qg), then qidx = -1.
         qidx = min( kpts_minus(kg%kpt(ikq), kg%kpt(ik), kg%ndim, qg%ndim), &
               kpts_minus(kg%kpt(ik), kg%kpt(ikq), kg%ndim, qg%ndim) )
         !debug
         !if(qidx > maxq) call errore('setup_kq_pair', 'qidx > maxq', 1)
         !
         !if this xq is on the q-grid qg.
         if(qidx > -1) then
            kq_tmp(i)%iq(j)  = qidx
!$omp atomic
            qcollect(qidx) = qcollect(qidx) + 1
         endif
      enddo 
   enddo
!$omp end parallel do
   
   !collect q-points
   iq = 0
   do qidx = 0, maxq
      if(qcollect(qidx) > 0) then
         iq = iq + 1
         qcollect(qidx) = iq
      endif
   enddo
   nq_tmp = iq
   
   allocate( ql_tmp(nq_tmp) )
   do qidx = 0, maxq
      iq = qcollect(qidx)
      if( iq > 0) ql_tmp( iq ) = qidx
   enddo
   
   !update kq_pair%iq
   do i = 1, nk_loc
      do j = 1, kq_tmp(i)%nkq
         qidx = kq_tmp(i)%iq(j)
         !qidx = -1 if xq is not on the q-grid.
         if(qidx > -1) then
            iq = qcollect(qidx)
            kq_tmp(i)%iq(j) = iq
         endif
      enddo
   enddo
   !release space
   deallocate(qcollect)

   !allocate space
   allocate(qw_tmp(nmod, nq_tmp))
   do i = 1, nk_loc
      nkq = kq_tmp(i)%nkq
      allocate( kq_tmp(i)%ikq(nkq)  )
      allocate( kq_tmp(i)%nchl(nkq) )
      !initialize
      kq_tmp(i)%nchl(:) = 0
   enddo

   !compute phonon freq for qpts
!$omp parallel default(shared) private(iq, xq, i, ik, j, ikq, ltable) 
!$omp do schedule(guided)
   do iq = 1, nq_tmp
      xq  = num2kpt( ql_tmp(iq),  qg%ndim )
      call solve_phonon_modes(ph, xq, qw_tmp(:,iq))
   enddo
!$omp end do

!$omp do schedule(guided)
   do i = 1, nk_loc
      ik = kq_tmp(i)%ik
      !nkq = (nkpt-1)/2 + merge(1, 0, ((mod(nkpt-1,2)>0) .and. ik<=nkpt/2))
      do j = 1, kq_tmp(i)%nkq
         ikq = mod(ik+j-1, nkpt) + 1
         kq_tmp(i)%ikq(j) = ikq
         !
         iq = kq_tmp(i)%iq(j)
         if(iq > 0) then
            call match_table(kg%enk(:,ikq), kg%enk(:,ik), qw_tmp(:,iq), wcut, e_thr, ltable)
            kq_tmp(i)%nchl(j) = count(ltable)
         endif
      enddo
   enddo
!$omp end do
!$omp end parallel

   !update qg%list
   allocate( q_check(nq_tmp) )
   q_check = 0
   !collect (k,q) pair that has at least one valid scattering channel (nchl > 0).
   do i= 1, nk_loc
      nkq = count( kq_tmp(i)%nchl(:) > 0 )
      kq_pair(i)%nkq = nkq
      allocate( kq_pair(i)%iq(nkq), kq_pair(i)%ikq(nkq), kq_pair(i)%nchl(nkq) )
      !
      k = 0
      do j = 1, kq_tmp(i)%nkq
         if( kq_tmp(i)%nchl(j) > 0 ) then
            k = k + 1
            !
            kq_pair(i)%ikq(k) = kq_tmp(i)%ikq(j)
            kq_pair(i)%nchl(k) = kq_tmp(i)%nchl(j)
            iq = kq_tmp(i)%iq(j)
            kq_pair(i)%iq(k) = iq
            !
            q_check(iq) = q_check(iq) + 1
         endif
      enddo
      !debug
      if(k .ne. nkq) call errore('setup_kq_pair', 'k .ne. nkq', 1)
      !
      deallocate(kq_tmp(i)%iq, kq_tmp(i)%ikq, kq_tmp(i)%nchl)
   enddo
   deallocate( kq_tmp )
   
   !if number of q-points is reduced.
   if( count(q_check > 0) .ne. nq_tmp ) then
      k = 0
      do i = 1, nq_tmp
         if(q_check(i) > 0) then
            k = k + 1
            q_check(i) = k
         endif
      enddo
      qg%npts = k
      
      !update qg%list
      allocate( qg%list( qg%npts ), qg%freq(nmod, qg%npts) )
      do i = 1, nq_tmp
         if(q_check(i) > 0) then
            qg%list( q_check(i) ) = ql_tmp(i)
            qg%freq(:, q_check(i)) = qw_tmp(:, i)
         endif
      enddo
      deallocate( ql_tmp, qw_tmp )

      !update kq_pair
      do i = 1, nk_loc
         do j = 1, kq_pair(i)%nkq
            iq = kq_pair(i)%iq(j)
            kq_pair(i)%iq(j) = q_check(iq)
         enddo
      enddo
   else
      qg%npts = nq_tmp
      qg%list => ql_tmp
      qg%freq => qw_tmp
   endif
   deallocate( q_check )

   call stop_clock('set_kq_pair')
end subroutine setup_kq_pair

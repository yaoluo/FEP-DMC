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

subroutine calc_electron_wann(el, kdim, nwan, at, wan_center, nband, numk, xk, eig, rot_wan)
   use kinds, only: dp
   use constants, only: tpi
   use qe_mpi_mod, only: mp_sum, inter_pool_comm, mp_split_pools
   use electron_wannier, only: electron_wann, set_ws_cell_el
   implicit none
   type(electron_wann), intent(out) :: el
   ! number of k-point and bands (may differ from number of wannier functons)
   integer, intent(in) :: kdim(3), nwan, nband, numk
   ! wannier center in crystal coordinate
   real(dp), intent(in) :: at(3,3), wan_center(3, nwan)
   real(dp), intent(in) :: xk(3,numk) ! k-points in crystal coordinate
   real(dp), intent(in) :: eig(nband, numk) ! band energy eig(nband, numk)
   complex(dp), intent(in) :: rot_wan(nband, nwan, numk) ! rot_wan(nband, el%nb, numk)
   !
   integer :: m, ib, jb, ik, n, ir, ire, nk_tmp, ik_st, ik_end, nk_loc
   real(dp) :: fac, rdotk, rvec(3)
   complex(dp) :: cfac, ham_k
   
   nk_tmp = kdim(1) * kdim(2) * kdim(3)
   !sanity check
   if(nk_tmp .ne. numk) &
      call errore('calc_electron_wann','inconsistent kdim and numk',1)
   ! setup wigner seitz cell
   call set_ws_cell_el(el, kdim, nwan, at, wan_center)
   ! for mpi parallel
   call mp_split_pools(numk, ik_st, ik_end, nk_loc)
   
   m = 0
   do jb = 1, nwan
   do ib = 1, jb
      !m = ( jb * (jb-1) ) / 2 + ib
      m = m + 1
      !init
      el%ham_r(m)%hop(:) = cmplx(0.0_dp, 0.0_dp, kind=dp)

      do ik = ik_st, ik_end
         !Preparing H(k)_{ij} = [U^{\dagger} Diag(eig) U]_{ij}
         ham_k = cmplx(0.0_dp, 0.0_dp, kind=dp)
         do n = 1, nband
            ham_k = ham_k + conjg(rot_wan(n,ib,ik)) * eig(n,ik) * rot_wan(n,jb,ik)
         enddo
         !N.B. H(k)_{ji} = conjg( H(k)_{ij} ) since H(k) is herminian

!$omp parallel do schedule(guided) default(shared) private(ir, ire, rvec, rdotk, cfac)
         do ir = 1, el%ham_r(m)%ws_el%nr
            ire = el%ham_r(m)%ws_el%rvec(ir)
            rvec(:) = el%rvec_set(:, ire)
            ! 
            rdotk = tpi * dot_product(xk(:,ik), rvec)
            ! exp(-i k*R) / N_k
            cfac = cmplx(cos(rdotk), -sin(rdotk), kind=dp)
            ! ham_r(ir) = ham_r(ir) + ham_k * factor
            el%ham_r(m)%hop(ir) = el%ham_r(m)%hop(ir) + ham_k * cfac
         enddo
!$omp end parallel do

      enddo
   enddo; enddo
   
   !collect results from different pools
   do m = 1, (nwan*(nwan+1))/2
      call mp_sum(el%ham_r(m)%hop(:), inter_pool_comm)
      ! include the ndeg
      do ir = 1, el%ham_r(m)%ws_el%nr
         fac = 1.0_dp / real( el%ham_r(m)%ws_el%ndeg(ir) * numk, dp)
         el%ham_r(m)%hop(ir) = el%ham_r(m)%hop(ir) * fac
      enddo
   enddo
end subroutine calc_electron_wann


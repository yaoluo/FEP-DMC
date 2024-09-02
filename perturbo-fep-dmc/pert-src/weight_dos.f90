!===============================================================================
! Copyright (C) 2016-2020 Jin-Jian Zhou, Jinsoo Park, I-Te Lu, Marco Bernardi
!
! This program is distributed under the terms of the GNU General Public License.
! See the file `LICENSE' in the root directory of this distribution, or obtain 
! a copy of the License at <https://www.gnu.org/licenses/gpl-3.0.txt>.
!
! Author: jjzhou <jjchou.comphy@gmail.com>
! Comment:
!   evaluate integral weight for the four coner of an tetra.
!
!   the weight are based on appendix B: since delta{e} = d step_func{e} / de
!   and appendix B give the formula to evalute weight for step_func{e_fermi},
!   replace e_fermi with e, the derivative of e is the weihgt we needed.
!   ps: the weights evaluate here are without V_T / V_G.
!
!   refers: Peter Blochl, "Improved tetrahedron method or Brillouin-zone
!   integrations" PRB 49 16 333 (1994),
!
! Maintenance:
!==============================================================================

subroutine weight_dos(e_corner, e_dos, weight)
   use pert_const, only: dp
    implicit none
    real(dp), intent(in) :: e_corner(4)
    real(dp), intent(in) :: e_dos
    real(dp), intent(out) :: weight(4)
    
    integer :: e_index(4)
    real(dp)  :: e(4), etmp, factor
    real(dp)  :: et1, et2, et3, et4, w(4)
    real(dp)  :: c1, c2, c3, dc1, dc2, dc3
    real(dp)  :: e31, e32, e41, e42
    
    integer :: i, j, k, tmp
    
    ! order the e_conear
    ! e1 <= e2 <= e3 <= e4
    do i=1, 4
        e_index(i) = i
        e(i) = e_corner(i)
    end do
        
    do i=1, 3
        k = i
        do j=i+1, 4
            if(e(k) > e(j)) then
                k = j
            end if
        end do
        if (i /= k) then
            !switch e(i) and e(k)
            etmp = e(i)
            e(i) = e(k)
            e(k) = etmp
            tmp = e_index(i)
            e_index(i) = e_index(k)
            e_index(k) = tmp
        end if
    end do
    !write(*,*) e  !write(*,*) e_index
    
    ! evaluate weight of the four corner
    if(e_dos < e(1)) then
        w = 0.0E0_dp
    else if (e_dos < e(2) .or. (e_dos==e(2).and.e(2)==e(3).and.e(3)==e(4).and.e(2)>e(1))) then
        et1 = e_dos - e(1)
        et2 = e(2)  - e(1)
        et3 = e(3)  - e(1)
        et4 = e(4)  - e(1)
        factor = et1*et1 / (et2*et3*et4)
        w(2) = factor * et1  / et2
        w(3) = factor * et1  / et3
        w(4) = factor * et1  / et4
        w(1) = 3.0E0_dp*factor - w(2) - w(3) - w(4)
    else if(e_dos < e(3)) then
        et1 = e_dos - e(1)
        et2 = e_dos - e(2)
        et3 = e(3) - e_dos
        et4 = e(4) - e_dos
        e31 = e(3) - e(1)
        e32 = e(3) - e(2)
        e41 = e(4) - e(1)
        e42 = e(4) - e(2)
        dc1 = 0.5E0_dp*et1 / (e41*e31)
        c1  = 0.5*et1*dc1
        c2 = 0.25E0_dp / (e41*e32*e31)
        dc2 = c2*(et2*et3 + et1*et3 - et1*et2)
        c2 = c2*et1*et2*et3
        c3 = 0.25E0_dp / (e42*e32*e41)
        dc3 = c3*(2.0E0_dp*et2*et4 - et2*et2)
        c3 = c3*(et2*et2*et4)
        w(1) = dc1 + (dc1 + dc2)*et3/e31 + (dc1 + dc2 + dc3)*et4/e41
        w(1) = w(1) - (c1 + c2) / e31 - (c1 + c2 +c3) / e41
        w(2) = dc1 + dc2 + dc3 + (dc2 + dc3)*et3/e32 + dc3*et4/e42
        w(2) = w(2) - (c2 + c3)/e32 - c3/e42
        w(3) = (dc1 + dc2)*et1/e31 + (dc2 + dc3)*et2/e32
        w(3) = w(3) + (c1 + c2)/e31 + (c2 + c3)/e32
        w(4) = (dc1 + dc2 + dc3)*et1/e41 + dc3*et2/e42
        w(4) = w(4) + (c1 + c2 + c3)/e41 + c3/e42
    else if(e_dos < e(4)) then
        et1 = e(4) - e(1)
        et2 = e(4) - e(2)
        et3 = e(4) - e(3)
        et4 = e(4) - e_dos
        factor = et4*et4 / (et1*et2*et3)
        w(1) = factor * et4  / et1
        w(2) = factor * et4  / et2
        w(3) = factor * et4  / et3
        w(4) = 3.0E0_dp*factor - w(1) - w(2) - w(3)
    else
        w = 0.0E0_dp 
    end if
    
    do i=1,4
        weight(e_index(i)) = w(i)
    end do

    return
end subroutine weight_dos

! tetrahedron integration for a given tetrahedron
subroutine tetra_int(numb, et, fnk, nstep, ene, tsum)
   use pert_const,  only: dp
   implicit none
   integer,  intent(in)  :: numb, nstep
   real(dp), intent(in)  :: ene(nstep), et(numb, 4), fnk(numb, 4)  
   real(dp), intent(inout) :: tsum(nstep)
   ! local variables
   real(dp) :: wt(4)
   integer :: ib, ie
   !tsum(:) = 0.E0_dp
   ! contribution of current tetra
   do ie = 1, nstep
      do ib = 1, numb
         call weight_dos(et(ib,:), ene(ie), wt)
         !for test, weight for integrated dos, e.g. number of states
         !call weight(et(ib,:), ene(ie), wt)
         tsum(ie) = tsum(ie) + dot_product(wt(1:4), fnk(ib,1:4))
      enddo
   enddo
end subroutine tetra_int

!program test
!subroutine test_weight_dos
!    implicit none
!
!    real*8 :: e(4) 
!    data e / 0.3d0, 4.0d0, 0.41d0, 0.42d0 /
!    real*8 :: e_d = 0.409d0, temp1 = 0.d0, temp2 = 0.0d0
!    integer :: i
!    real*8 :: weight_1(4), weight_2(4), weight_3(4)
!    write(*,*) e
!
!    call weight_dos(e,e_d,weight_1)
! !   call weight(1.0d0, e, e_d, weight_2)
!!    call T_weight_delta(e,e_d,weight_3)
!    do i = 1, 4
!        temp1 = temp1 + weight_1(i)
!!        temp2 = temp2 + weight_3(i)
!    enddo 
!    write(*,*) temp1
!    write(*,*) weight_1
! !   write(*,*) weight_2
!!    write(*,*) weight_3
!end subroutine test_weight_dos  

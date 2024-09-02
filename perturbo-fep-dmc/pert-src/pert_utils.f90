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

module pert_utils
   use pert_const, only: dp, twopi, czero, pi
   use qe_mpi_mod, only: my_pool_id
   implicit none
   private
   real(dp), parameter :: maxexp = 300.0_dp
   public :: find_free_unit, binary_search, hmat_diag, hmat_upper_diag, match_table
   public :: get_exp_ikr, fermi, mfermi_deriv, abs2, gauss, bose, converged
   public :: random_uniform, random_cauchy, fermi_2nd_deriv, mat_inv3
   public :: spline_coef, spline_interp, spline_interp_interval, deriv_1st, deriv_2nd
contains

!initialize a random seed from the system clock at every run
subroutine init_random_seed()
   integer :: i, n, clock
   integer, allocatable :: seed(:)

   call random_seed(size=n)
   allocate(seed(n))

   call system_clock(count=clock)
   seed = clock + 37 * (/(i-1, i=1, n)/) + my_pool_id
   call random_seed(put=seed)
   
   deallocate(seed)
end subroutine init_random_seed

subroutine random_uniform(pts, weight)
   implicit none
   real(dp), intent(out) :: weight(:), pts(:,:)
   
   call init_random_seed()
   weight = 1.0_dp / real(size(weight), dp)
   call random_number(pts)
end subroutine random_uniform

!generate vectors with random x, y, z compounts from a Cauchy distribution
! i.e. x ~ 1/(pi*eps) * 1/((x/eps)**2 + 1), using inversion method
subroutine random_cauchy(pts, weight, eps_in)
   implicit none
   real(dp), intent(out) :: weight(:), pts(:,:)
   real(dp), intent(in), optional :: eps_in

   integer :: np, i, nv, j
   real(dp) :: eps, pa, pb, wtmp
   real(dp), allocatable :: temp(:)

   nv = size(pts, 1); np = size(pts, 2) !dimension of vector, number of vectors
   allocate( temp(nv) )

   if(present(eps_in)) then
      eps = eps_in
   else
      eps = 0.05_dp   !default value
   endif
   call init_random_seed()
   
   pa = 0.5_dp + atan(-0.5_dp/eps) / pi
   pb = 0.5_dp + atan( 0.5_dp/eps) / pi
   wtmp = 1.0_dp / real(np, dp)
   do i = 1, np
      !uniform random number between (0.0, 1.0)
      call random_number(temp)
      temp(:) = pa + temp(:)*(pb - pa)
      pts(:,i) = eps * tan( (temp(:) - 0.5_dp)*pi )

      !weight is the inverse of probalility density
      temp(:) = ((pts(:,i)/eps)**2 + 1.0_dp) * eps * pi * (pb-pa)

      !weight(i) = temp(1)*temp(2)*temp(3)*wtmp
      weight(i) = wtmp
      do j = 1, nv
         weight(i) = weight(i) * temp(j)
      enddo
   enddo
   
   deallocate(temp)
   !if(abs(sum(weight)-1.0_dp) > 1.0E-2_dp) call errore('random_cauchy',&
   !   'total weight is not 1, need more sampling',-1)
end subroutine random_cauchy



integer function find_free_unit() result(iunit)
   implicit none
   ! local variables
   integer :: test_unit
   logical :: isopen
   !init
   test_unit = 30
   isopen = .true.
   ! find availale iunit 
   do while( isopen )
      test_unit = test_unit + 1
      inquire(test_unit, opened=isopen)
   enddo
   iunit = test_unit
end function find_free_unit


! binary searach in a ordered integer array (list)
integer pure function binary_search(list, val) result(loc)
   implicit none
   integer, intent(in) :: list(:), val
   ! local variables
   integer :: st, ed, mid, dis

   st = 1;  ed = size(list)
   dis = ed - st
   mid = (st + ed) / 2

   do while ( dis > 0 )
      if( val > list(mid) ) then
         st = mid + 1
      elseif( val < list(mid) ) then
         ed = mid - 1
      else
         !find the result: val == list(mid)
         loc = mid
         return
      endif
      dis = ed - st
      mid = (st + ed) / 2
   enddo

   if( mid < 1 .or. mid > size(list) ) then
      loc = 0
   else
      loc = merge(mid,  0, list(mid) .eq. val)
   endif
end function binary_search


subroutine match_table(ekq, ek, wq, wcut, e_thr, ltable)
   implicit none
   real(dp), intent(in):: ekq(:), ek(:), wq(:), wcut, e_thr
   logical, intent(out) :: ltable(:,:,:) !ltable(size(ekq), size(ek), size(wq))
   integer :: i, j, m

   ltable = .false.
   do m = 1, size(wq)
      if(wq(m) < wcut) cycle
      !if e_thr is negtive, then fully relax the energy conservation.
      if(e_thr < 0.0_dp) then
         ltable(:,:,m) = .true.
      else
         do i = 1, size(ek)
         do j = 1, size(ekq)
            ltable(j,i,m) = abs(abs(ek(i)-ekq(j)) - wq(m)) < e_thr
         enddo; enddo
      endif
   enddo
end subroutine match_table

! similar to utility_diagonalize in wannier90 code.
! diagonalize a complex hermitian matrix
subroutine hmat_diag(hmat, ndim, eig, urot)
   implicit none
   integer, intent(in) :: ndim
   complex(dp), intent(in) :: hmat(ndim, ndim)
   real(dp), intent(out) :: eig(ndim)
   complex(dp), intent(out), optional :: urot(ndim,ndim)
   !local variables
   integer :: neig, info, ifail( ndim ), iwork( 5*ndim ), ib, jb
   real(kind=dp) :: rwork( 7*ndim )
   complex(kind=dp) :: hmpack( ndim*(ndim+1)/2 ), cwork(2*ndim), ev(ndim,ndim)
   !hmpack: complex hamiltonian packed (upper triangular part for zhpevx)
   !force hermitian-ization
   do jb = 1, ndim
      do ib = 1, jb
         hmpack( ib + (jb-1)*jb/2 ) = ( hmat(ib,jb) + dconjg(hmat(jb,ib)) )*0.5_dp
      enddo
   enddo
   !diagonalize mat using lacpack
   eig = 0.0_dp; cwork = czero; rwork = 0.0_dp; iwork = 0
   if(present(urot)) then
      call zhpevx ('V', 'A', 'U', ndim, hmpack, 0.0_dp, 0.0_dp, 0, 0, -1.0_dp, &
      neig, eig(1:ndim), ev, ndim, cwork, rwork, iwork, ifail, info)
      urot = ev
   else
      call zhpevx ('N', 'A', 'U', ndim, hmpack, 0.0_dp, 0.0_dp, 0, 0, -1.0_dp, &
      neig, eig(1:ndim), ev, ndim, cwork, rwork, iwork, ifail, info)
   endif

   if(info < 0) &
      call errore('ham_diagonalize','ARGUMENT OF ZHPEVX HAD AN ILLEGAL VALUE',1)
   if(info > 0) &
      call errore('ham_diagonalize','EIGENVECTORS FAILED TO CONVERGE',1)
end subroutine hmat_diag

! similar to utility_diagonalize in wannier90 code.
! diagonalize a complex hermitian matrix
subroutine hmat_upper_diag(hmpack, ndim, eig, urot)
   implicit none
   integer, intent(in) :: ndim
   !hmpack: complex hamiltonian packed (upper triangular part for zhpevx)
   complex(dp), intent(in) :: hmpack(:) ! ndim*(ndim+1)/2
   real(dp), intent(out) :: eig(ndim)
   complex(dp), intent(out), optional :: urot(ndim,ndim)
   !local variables
   integer :: neig, info, ifail( ndim ), iwork( 5*ndim )
   real(kind=dp) :: rwork( 7*ndim )
   complex(kind=dp) :: cwork(2*ndim), ev(ndim,ndim) ! hmpack( ndim*(ndim+1)/2 ), 
   !hmpack: complex hamiltonian packed (upper triangular part for zhpevx)
   !force hermitian-ization
   !do jb = 1, ndim
   !   do ib = 1, jb
   !      hmpack( ib + (jb-1)*jb/2 ) = ( hmat(ib,jb) + dconjg(hmat(jb,ib)) )*0.5_dp
   !   enddo
   !enddo
   if(size(hmpack) .ne. ndim*(ndim+1)/2) &
      call errore('hmat_upper_diag','size of hmpack mismatch with ndim',1)

   !diagonalize mat using lacpack
   eig = 0.0_dp; cwork = czero; rwork = 0.0_dp; iwork = 0
   if(present(urot)) then
      call zhpevx ('V', 'A', 'U', ndim, hmpack, 0.0_dp, 0.0_dp, 0, 0, -1.0_dp, &
      neig, eig(1:ndim), ev, ndim, cwork, rwork, iwork, ifail, info)
      urot = ev
   else
      call zhpevx ('N', 'A', 'U', ndim, hmpack, 0.0_dp, 0.0_dp, 0, 0, -1.0_dp, &
      neig, eig(1:ndim), ev, ndim, cwork, rwork, iwork, ifail, info)
   endif

   if(info < 0) &
      call errore('ham_diagonalize','ARGUMENT OF ZHPEVX HAD AN ILLEGAL VALUE',1)
   if(info > 0) &
      call errore('ham_diagonalize','EIGENVECTORS FAILED TO CONVERGE',1)
end subroutine hmat_upper_diag

!complex(dp) pure function exp_ikr(xpt, vec)
!   implicit none
!   real(dp), intent(in) ::xpt(3), vec(3)
!   real(dp) :: rdotk
!
!   rdotk = dot_product(xpt, vec) * twopi
!   exp_ikr = dcmplx(cos(rdotk), sin(rdotk))
!end function exp_ikr


subroutine get_exp_ikr(xpt, vec, exp_ikr)
   implicit none
   real(dp), intent(in) ::xpt(3), vec(:,:) !vec(3,:)
   complex(dp), intent(out) :: exp_ikr(:)
   ! local variable
   integer :: i
   real(dp) :: rdotk

   do i = 1, size(exp_ikr)
      rdotk = dot_product(xpt, vec(1:3,i)) * twopi
      exp_ikr(i) = dcmplx(cos(rdotk), sin(rdotk))
   enddo
end subroutine get_exp_ikr


real(dp) elemental function fermi(kt, de) result(f)
   implicit none
   real(dp), intent(in) :: kt, de
   ! -maxe < ekt < maxe: f = 1/(e^{ekt} + 1)
   ! ekt < -maxe : f = 1/(e^{-maxe} + 1) -> 1.0
   ! ekt >  maxg : f = 1/(e^{ maxe} + 1) -> 0.0
   f = 1.0_dp / (1.0_dp + exp( sign(min(abs(de/kt), maxexp), de) ))
end function fermi

! compute -df/dE
real(dp) elemental function mfermi_deriv(kt, de) result(df)
   implicit none
   real(dp), intent(in) :: kt, de  ! kt: k*T; 
   !real(dp) :: ekt !ekt: (e-mu)/kT
   !ekt = de/kt
   ! -maxe < ekt < maxe: f = 1/kT * f(ekt)*(1-f(ekt)) = 1/kT * f(ekt) * f(-ekt)
   ! otherwise: 0.0_dp
   !df = merge(exp(ekt)/((exp(ekt)+1.0_dp)**2)/kt, 0.0_dp, abs(ekt) < maxexp)

   df = fermi(kt, de) * fermi(kt, -de) * (1.0_dp / kt)

   !!all the arguments of merge will be evaluated, need to use if .. else .. to avoid.
   !if( abs(ekt) < maxexp ) then
   !   df = exp(ekt)/((exp(ekt)+1.0_dp)**2) / kt
   !else
   !   df = 0.0_dp
   !endif
end function mfermi_deriv

! compute df/dE^2
! df/dE^2 = e^{ekt}/(e^{ekt}+1)^2 * (e{ekt}-1)/(e{ekt}+1) / kT^2
real(dp) elemental function fermi_2nd_deriv(kt, de) result(df)
   implicit none
   real(dp), intent(in) :: kt, de  ! kt: k*T; 
   real(dp) :: ekt, factor !ekt: (e-mu)/kT 
   ekt = de/kt
   
   if( abs(ekt) < maxexp ) then
      factor = (exp(ekt)-1.0E0_dp)/(exp(ekt)+1.0E0_dp) / kt
      df = factor * exp(ekt)/((exp(ekt)+1.0_dp)**2) / kt
   else
      df = 0.0E0_dp
   endif
end function fermi_2nd_deriv


real(dp) elemental function abs2( var ) result(f)
   implicit none
   complex(dp), intent(in) :: var
   f = real(dconjg(var)*var, dp)
end function abs2


real(dp) elemental function gauss(eta, x) result(f)
   implicit none
   real(dp), intent(in) :: eta, x
   real(dp), parameter :: pifactor = 1.0_dp/sqrt(pi)
   !gaussain function f(eta, x) = e^{-(x/eta)**2}/eta/sqrt(pi)
   f = exp( -min((x/eta)**2, 200.0_dp) )*pifactor/eta
end function gauss


real(dp) elemental function bose(kt, wq) result(f)
   implicit none
   real(dp), intent(in) :: kt, wq  ! need to make sure wq > 0
   f = 1.0_dp / (exp( min(wq/kt, maxexp) ) - 1.0_dp)
end function bose

logical pure function converged(prev, current, thr)
   implicit none
   real(dp), intent(in) :: prev(:), current(:), thr
   integer :: i
   real(dp) :: dist, diff
   
   dist = 0.0_dp;  converged = .false.
   do i = 1, size(prev)
      diff = current(i)-prev(i)
      dist = dist + diff*diff
   enddo
   diff = sqrt(dist) / sqrt(dot_product(prev, prev))

   if(diff < thr) converged = .true.
end function converged


!Adapted from spline.f90: Alex G: January 2010
!https://ww2.odu.edu/~agodunov/computing/programs
subroutine spline_coef(x, y, b, c, d)
!======================================================================
!  Calculate the coefficients b(i), c(i), and d(i), i=1,2,...,n
!  for cubic spline interpolation
!  s(x) = y(i) + b(i)*(x-x(i)) + c(i)*(x-x(i))**2 + d(i)*(x-x(i))**3
!----------------------------------------------------------------------
!  input..
!  x = the arrays of data abscissas (in strictly increasing order)
!  y = the arrays of data ordinates
!======================================================================
   implicit none
   real(dp), intent(in)  :: x(:), y(:)  ! x,y,b,c,d have the same size
   real(dp), intent(out) :: b(:), c(:), d(:)  
   !local 
   integer :: n, gap, i, j
   real(dp) :: h

   n = size(x)
   gap = n - 1

   if(n < 2) return
   if(n < 3) then      ! linear interpolation
      b(1) = (y(2)-y(1))/(x(2)-x(1));  c(1) = 0.0_dp;  d(1) = 0.0_dp
      b(2) = b(1);  c(2) = 0.0_dp;  d(2) = 0.0_dp
      return
   end if
    
   !step 1: preparation
   d(1) =  x(2) - x(1)
   c(2) = (y(2) - y(1))/d(1)
   do i = 2, gap
      d(i) = x(i+1) - x(i)
      b(i) = 2.0*(d(i-1) + d(i))
      c(i+1) = (y(i+1) - y(i))/d(i)
      c(i) = c(i+1) - c(i)
   end do
    
   ! step 2: end conditions 
   b(1) = -d(1)
   b(n) = -d(n-1)
   c(1) = 0.0
   c(n) = 0.0
   if(n /= 3) then
     c(1) = c(3)/(x(4)-x(2)) - c(2)/(x(3)-x(1))
     c(n) = c(n-1)/(x(n)-x(n-2)) - c(n-2)/(x(n-1)-x(n-3))
     c(1) = c(1)*d(1)**2/(x(4)-x(1))
     c(n) = -c(n)*d(n-1)**2/(x(n)-x(n-3))
   end if
    
   !step 3: forward elimination 
   do i = 2, n
     h = d(i-1)/b(i-1)
     b(i) = b(i) - h*d(i-1)
     c(i) = c(i) - h*c(i-1)
   end do
    
   !step 4: back substitution
   c(n) = c(n)/b(n)
   do j = 1, gap
     i = n-j
     c(i) = (c(i) - d(i)*c(i+1))/b(i)
   end do
    
   !step 5: compute spline coefficients
   b(n) = (y(n) - y(gap))/d(gap) + d(gap)*(c(gap) + 2.0_dp*c(n))
   do i = 1, gap
     b(i) = (y(i+1) - y(i))/d(i) - d(i)*(c(i+1) + 2.0_dp*c(i))
     d(i) = (c(i+1) - c(i))/d(i)
     c(i) = 3.0_dp*c(i)
   end do
   c(n) = 3.0_dp*c(n)
   d(n) = d(n-1)
end subroutine spline_coef

real(dp) pure function spline_interp(x, y, b, c, d, x_in) result(y_out)
!======================================================================
! evaluates the cubic spline interpolation at point x_in
! y_out = y(i)+b(i)*(x_in-x(i))+c(i)*(x_in-x(i))**2+d(i)*(x_in-x(i))**3
! where  x(i) <= x_in <= x(i+1)
!======================================================================
   implicit none 
   real(dp), intent(in) :: x(:), y(:), b(:), c(:), d(:)
   real(dp), intent(in) :: x_in
   !local
   integer :: n, i, j, k
   real(dp) :: dx
   
   n = size(x)

   if(x_in <= x(1)) then
      y_out = y(1)
      return
   elseif (x_in >= x(n)) then
      y_out = y(n)
      return
   endif

   !binary search for for i, such that x(i) <= x_in <= x(i+1)
   i = 1
   j = n+1
   do while (j > i+1)
      k = (i+j)/2
      if(x_in < x(k)) then
         j=k
      else
         i=k
      end if
   end do
   !evaluate spline interpolation
   dx = x_in - x(i)
   y_out = y(i) + dx*(b(i) + dx*(c(i) + dx*d(i)))
end function spline_interp


subroutine spline_interp_interval(yi, bi, ci, di, dx_in, y_out)
   implicit none
   real(dp), intent(in) :: yi, bi, ci, di, dx_in(:) !dx_in = x_in - xi
   real(dp), intent(out) :: y_out(:) !dx_in and y_out should have the same size

   integer :: n, i
   real(dp) :: dx

   n = size(dx_in)
   do i = 1, n
      !the largest dx should be smaller then x(i+1) - x(i)
      dx = dx_in(i)
      y_out(i) = yi + dx*(bi + dx*(ci + dx*di))
   enddo
end subroutine spline_interp_interval

!compute first derivative of y(x) using five points (O(dx^4))
!see https://en.wikipedia.org/wiki/Finite_difference_coefficient#Central_finite_difference
! and http://web.media.mit.edu/~crtaylor/calculator.html
real(dp) pure function deriv_1st(y, dx) result(dydx)
   ! y(-2), y(-1), y(0), y(1), y(2)
   real(dp), intent(in) :: y(5), dx
   !!coefficients
   real(dp), parameter :: coef(5) = (/ 1.0_dp, -8.0_dp, 0.0_dp, 8.0_dp, -1.0_dp/)
   real(dp), parameter :: denominater = 12.0_dp

   dydx = dot_product(coef, y) / (denominater * dx)
end function deriv_1st

!compute second derivative of y(x) using five points (O(dx^4))
real(dp) pure function deriv_2nd(y, dx) result(dydx2)
   real(dp), intent(in) :: y(5), dx
!   !coefficients
   real(dp), parameter :: coef(5) = (/ -1.0_dp, 16.0_dp, -30.0_dp, 16.0_dp, -1.0_dp/)
   real(dp), parameter :: denominater = 12.0_dp

   dydx2 = dot_product(coef, y) / (denominater * dx * dx)
end function deriv_2nd


!Performs a direct calculation of the inverse of a 3×3 matrix.
pure function mat_inv3(A) result(B)
  real(dp), intent(in) :: A(3,3)   !! Matrix
  real(dp) :: B(3,3)   !! Inverse matrix
  real(dp) :: detinv

  ! Calculate the inverse determinant of the matrix
  detinv = 1/(A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2)&
            - A(1,2)*A(2,1)*A(3,3) + A(1,2)*A(2,3)*A(3,1)&
            + A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1))

  ! Calculate the inverse of the matrix
  B(1,1) = +detinv * (A(2,2)*A(3,3) - A(2,3)*A(3,2))
  B(2,1) = -detinv * (A(2,1)*A(3,3) - A(2,3)*A(3,1))
  B(3,1) = +detinv * (A(2,1)*A(3,2) - A(2,2)*A(3,1))
  B(1,2) = -detinv * (A(1,2)*A(3,3) - A(1,3)*A(3,2))
  B(2,2) = +detinv * (A(1,1)*A(3,3) - A(1,3)*A(3,1))
  B(3,2) = -detinv * (A(1,1)*A(3,2) - A(1,2)*A(3,1))
  B(1,3) = +detinv * (A(1,2)*A(2,3) - A(1,3)*A(2,2))
  B(2,3) = -detinv * (A(1,1)*A(2,3) - A(1,3)*A(2,1))
  B(3,3) = +detinv * (A(1,1)*A(2,2) - A(1,2)*A(2,1))
end function

end module pert_utils

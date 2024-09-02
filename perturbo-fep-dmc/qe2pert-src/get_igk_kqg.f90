!===============================================================================
! Copyright (C) 2016-2020 Jin-Jian Zhou, Jinsoo Park, I-Te Lu, Marco Bernardi
!
! This program is distributed under the terms of the GNU General Public License.
! See the file `LICENSE' in the root directory of this distribution, or obtain 
! a copy of the License at <https://www.gnu.org/licenses/gpl-3.0.txt>.
!
! Author: jjzhou <jjchou.comphy@gmail.com>
! Comment:
!  generate igk for xk+xq. xk+xq might be xkq+G, where xkq is on the regular 
!  k-grid, so we need to map xkq + G1 => xk + xq + G2, a.k.a let G2 = G1 - G, 
!  where G1, G2 represent igk for k=xkq and k=xkq+G, respectively.
!
!  the reason is that we usually only compute \psi_k for k in the First BZ.
!  and use \psi_{k+G} = \psi_{k} to extend \psi to the whole BZ:
!      \psi_{k}  = \sum_{G1}  C_{k}(G1)  exp(i*(k+G1)*r)
!     \psi_{k+G} = \sum_{G2} C_{k+G}(G2) exp(i*(k+G+G2)*r)
!  to make sure \psi_{k+G} = \psi_{k}, for each G1, G2 satisfying 
!  (k+G1) = (k+G+G2), we enforce C_{k}(G1) = C_{k+G}(G2).
!  aka. C_{k}(G1) = C_{k+G}(G1-G), which we used to get \psi_{k+G} from \psi_{k}
!
!  so for k-point: xkq+G, we can use the same evc, just map G1 to G1-G instead.
!
! Maintenance:
!===============================================================================

subroutine get_igk_kqg(npwq, igk_k, xk, xkg, igk_kg)
   use kinds, only: dp
   use gvect, only: ngm, g
   use pwcom, only: npwx
   use cell_base, only: at, bg
   use io_global, only: stdout
   use gvect, only: ngm, g
   use gvecw, only: gcutw, ecutwfc
   
   implicit none
   integer,  intent(in) :: npwq, igk_k(npwq)
   real(dp), intent(in) :: xk(3), xkg(3) ! xk and xk+G in cartesian coordinate
   integer,  intent(out) :: igk_kg(npwq)
   
   real(dp), parameter :: eps6 = 1.0E-6_dp
   !local variables
   integer :: i, n, lw, up, npwqg, i_tmp(3)
   real(dp) :: xg(3), xg_cryst(3), xgg(3), xkk(3), xkgg(3), current_kin, kin, xkg_loc(3)
   integer, allocatable :: igk_tmp(:)
   real(dp), allocatable :: gkin(:)
   
   xkg_loc(:) = xkg(:)
   xg(:) = xkg_loc(:) - xk(:)
   !
   if( norm2( xg ) < eps6 ) then
      ! G = 0
      igk_kg(1:npwq) = igk_k(1:npwq)
      return
   else
      xg_cryst(1:3) = xg(1:3)
      ! transform to crystal coordinate.
      call cryst_to_cart(1, xg_cryst, at, -1)
      
      i_tmp(:) = nint( xg_cryst(:) )
      do i = 1, 3
         if(abs( xg_cryst(i) - i_tmp(i) ) > eps6 ) then
            write(stdout, '(a, 3(3f12.6,2x))') &
               'xk, xkg (cart.),  xkg-xk (cryst.): ', xk, xkg_loc, xg_cryst
            call errore('get_igk_kqg','xkg - xk is not a G,',1)
         else
            !remove the tiny error
            xg_cryst(i) = real(i_tmp(i), dp)
         endif
      enddo
      
      !back to cartesian coordinate
      call cryst_to_cart(1, xg_cryst, bg, 1)
      !update xg and xkg_loc
      xg(:) = xg_cryst(:)
      xkg_loc(:) = xk(:) + xg(:)
   endif
   
   allocate( igk_tmp(npwx), gkin(npwx) )
   
   !N.B. G is sorted based on |(k+G)|^2 (increasing order)
   call gk_sort( xkg_loc, ngm, g, gcutw, npwqg, igk_tmp, gkin )
   if( npwqg .ne. npwq ) call errore('get_igk_kqg','npwqg .ne. npwq', 1)
      
   current_kin = gkin(1)
   lw = 1
   up = npwqg
   do i = lw+1, npwqg
      if(abs(current_kin - gkin(i)) > eps6) then
         up = i - 1
         exit
      endif
   enddo

   igk_kg(1:npwq) = 0
   !
   ! one-to-one mapping: 
   ! if G1, G2 satisfy (k+G1) = (k+G+G2), then C_{k}(G1) = C_{k+G}(G2)
   ! search index of G2 based on |k+G1|^2 = |k+G+G2|^2, 
   ! and igk_tmp stores sorted index fo G2 based on increasing order of |k+G+G2|^2
   do n = 1, npwq
      ! k + G1
      xkk = xk + g(:, igk_k(n))
      kin = dot_product(xkk, xkk)
      ! G2 = G1 - G
      xgg = g(:, igk_k(n)) - xg
         
      ! kin should no smaller than current_kin
      if( kin < (current_kin - eps6) ) &
         call errore('get_igk_kqg','kin < current_kin',1)
      
      ! update current_kin, and [lw, up]
      if( (kin - current_kin) > eps6 ) then
         lw = up + 1
         if(lw > npwqg) call errore('get_igk_kqg','lw > npwqg',1)

         current_kin = gkin(lw)
         up = npwqg
         do i = lw+1, npwqg
            if(abs(current_kin - gkin(i)) > eps6) then
               up = i - 1
               exit
            endif
         enddo
      endif

      !search in [lw, up] for the index of G2
      if( abs(kin - current_kin) < eps6 ) then
         do i = lw, up
            if( norm2( xgg - g(:, igk_tmp(i)) ) < eps6 ) then
               igk_kg(n) = igk_tmp(i)
               exit
            endif
         enddo
      endif
   enddo

   if( any(igk_kg(1:npwq) < 1) ) call errore('get_igk_kqg','get igk_kg failed',1)
   deallocate( igk_tmp, gkin )
end subroutine get_igk_kqg

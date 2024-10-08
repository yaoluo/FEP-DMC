!
! Copyright (C) 2001-2018 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! adapted from PHONON/PH/dwfc.f90
!-------------------------------------------------------------------
SUBROUTINE ph_dwfc (npw_ , igk_ , xk_ , icart_ , func, dfunc)
  !-----------------------------------------------------------------
  !
  ! This routine calculates the first derivative of a wave function 
  ! w.r.t. the position operator r_icart:
  !
  ! dpsi/d(icart) = \sum_G [i(k+G)_icart] psi(G)
  !
  USE kinds,      ONLY : DP
  USE wvfct,      ONLY : npwx
  USE gvect,      ONLY : g
  USE cell_base,  ONLY : tpiba
  !USE klist,      ONLY : xk

  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: npw_
  INTEGER, INTENT(IN) :: igk_(npwx)
  INTEGER, INTENT(IN) :: icart_
  real(dp), intent(in) :: xk_(3)
  COMPLEX(DP), INTENT(IN) :: func(npwx)
  COMPLEX(DP), INTENT(OUT) :: dfunc(npwx)
  !
  ! Local variables
  !
  INTEGER :: ig
  REAL(DP) :: gvec, xk_aux
  !
  !CALL start_clock( 'dwfc' )
  !
  dfunc = (0.d0, 0.d0)
  !
  xk_aux = xk_(icart_) * tpiba
  DO ig = 1, npw_
     gvec = g(icart_, igk_(ig)) * tpiba
     dfunc(ig) = (0.d0,-1.d0) * (gvec + xk_aux) * func(ig)
  ENDDO
  !
  !CALL stop_clock( 'dwfc' )
  !
  RETURN 
  ! 
END SUBROUTINE ph_dwfc
!------------------------------------------------------------------------

!------------------------------------------------------------------------
SUBROUTINE ph_d2wfc (npw_ , igk_ , xk_ , icart_ , jcart_, func, d2func)
  !----------------------------------------------------------------------
  !
  ! This routine calculates the second derivative of a wave function 
  ! w.r.t. the position operator r_icart and r_jcart:
  !
  ! d2 psi/d(icart)d(jcart) = \sum_G [-(k+G)_icart * -(k+G)_gcart] psi(G)
  !
  USE kinds,      ONLY : DP
  USE wvfct,      ONLY : npwx
  USE gvect,      ONLY : g
  USE cell_base,  ONLY : tpiba
  !USE klist,      ONLY : xk
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: npw_
  INTEGER, INTENT(IN) :: igk_(npwx)
  INTEGER, INTENT(IN) :: icart_, jcart_
  real(dp), intent(in) :: xk_(3)
  COMPLEX(DP), INTENT(IN)  :: func(npwx)
  COMPLEX(DP), INTENT(OUT) :: d2func(npwx)
  !
  ! Local variables
  !
  INTEGER :: ig
  REAL(DP) :: gvec, xk_aux, gvec2, xk_aux2
  !
  d2func = (0.d0, 0.d0)
  !
  xk_aux  = xk_(icart_) * tpiba
  xk_aux2 = xk_(jcart_) * tpiba
  !
  DO ig = 1, npw_
     gvec    = g(icart_, igk_(ig)) * tpiba
     gvec2   = g(jcart_, igk_(ig)) * tpiba
     d2func(ig) = - (gvec+xk_aux) * (gvec2+xk_aux2) * func(ig)
  ENDDO
  !
  RETURN
  !
END SUBROUTINE ph_d2wfc
!-----------------------------------------------------------------------------

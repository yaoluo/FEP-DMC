!
! Copyright (C) 2001-2018 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! adapted from QE/LR_Modules/lr_orthoUwfc.f90
!-----------------------------------------------------------------------
!SUBROUTINE lr_orthoUwfc (lflag)
SUBROUTINE ph_lr_orthoUwfc (xk, npw, igk_kk, xkq, npwq, igk_kq, &
      vkb_kk, vkb_kq, wfcatomk, swfcatomk, wfcatomkpq, swfcatomkpq)
  !-----------------------------------------------------------------------
  !
  ! This routine computes and writes to file atomic wavefunctions phi 
  ! at k and k+q. phi(k) is defined in Eq. (A6) in Ref. [1]
  ! [1] Phys. Rev. B 98, 085127 (2018)
  !
  ! lflag=.TRUE.         write phi(k) and phi(k+q) to file (unit iuatwfc)
  ! lflag=.FALSE. DO NOT write phi(k) and phi(k+q) to file (unit iuatwfc)
  !
  ! Independently from the value of the lflag, this routine will always
  ! write S(k)*phi(k) and S(k+q)*phi(k+q) to file with unit iuatswfc
  ! (note that this is not the same unit as iuatwfc).
  !
  ! In the norm-conserving case, S(k)=1 and S(k+q)=1.
  ! Note: here the array wfcU is used as a workspace.
  ! Inspired by PW/src/orthoatwfc.f90
  !
  ! Written by I. Timrov (01.10.2018)
  !
  USE kinds,            ONLY : DP
  !USE io_files,         ONLY : nwordwfcU
  USE basis,            ONLY : natomwfc
  !USE klist,            ONLY : xk, ngk, igk_k
  USE wvfct,            ONLY : npwx
  USE control_flags,    ONLY : gamma_only
  USE uspp,             ONLY : nkb, okvan !vkb, 
  USE becmod,           ONLY : allocate_bec_type, deallocate_bec_type, bec_type !, becp
  !USE buffers,          ONLY : save_buffer
  !USE io_global,        ONLY : stdout
  !USE mp,               ONLY : mp_sum
  !USE mp_pools,         ONLY : intra_pool_comm
  USE noncollin_module, ONLY : noncolin, npol
  !USE qpoint,           ONLY : nksq, ikks, ikqs
  USE control_lr,       ONLY : lgamma
  !USE units_lr,         ONLY : iuatwfc, iuatswfc
  USE ldaU,             ONLY : U_projection, nwfcU, & !, wfcU, copy_U_wfc
                           Hubbard_l, is_hubbard, offsetU, oatwfc
  !
  IMPLICIT NONE
  !
  real(dp), intent(in) :: xk(3), xkq(3)
  integer, intent(in) :: npw, npwq, igk_kk(npwx), igk_kq(npwx)
  complex(dp), intent(in) :: vkb_kk(npwx, nkb), vkb_kq(npwx, nkb)
  !nwfcU is equal to natomwfc, check 'init_lda_plus_u'
  complex(dp), intent(out) :: wfcatomk(npwx*npol, nwfcU), wfcatomkpq(npwx*npol, nwfcU), &
     swfcatomk(npwx*npol, nwfcU), swfcatomkpq(npwx*npol, nwfcU)
  !
  !LOGICAL, INTENT(IN) :: lflag
  !
  !INTEGER :: ik,  & ! dummy k point index
  !           ikk, & ! k point
  !           ikq, & ! k+q point
  !           npw, & ! number of plane waves at k
  !           npwq   ! number of plane waves at k+q
  LOGICAL :: orthogonalize_wfc, normalize_only, lflag
  COMPLEX(DP), ALLOCATABLE :: wfcatom(:,:), &  ! atomic wfc
                              swfcatom(:,:)    ! S * atomic wfc
  !
  type(bec_type) :: becp
  !
  !CALL start_clock ('lr_orthoUwfc')
  !
  lflag = .true.
  !jjz: drop support for gamma only
  IF (gamma_only) CALL errore('ph_lr_orthoUwfc', &
      'Gamma-only calculation for this case not supported', 1)
  !
  IF (U_projection=="atomic") THEN
     orthogonalize_wfc = .FALSE.
     normalize_only = .FALSE.
  !   WRITE( stdout, '(/5x,"Atomic wfc used for the DFT+U projector are NOT orthogonalized")')
  ELSEIF (U_projection=="ortho-atomic") THEN
     orthogonalize_wfc = .TRUE.
     normalize_only = .FALSE.
     !WRITE( stdout, '(/5x,"Atomic wfc used for the DFT+U projector are orthogonalized")')
     !IF (gamma_only) CALL errore('lr_orthoUwfc', &
     !     'Gamma-only calculation for this case not implemented', 1 )
  ELSEIF (U_projection=="norm-atomic") THEN
     orthogonalize_wfc = .TRUE.
     normalize_only = .TRUE.
     !WRITE( stdout, '(/5x,"Atomic wfc used for the DFT+U projector are normalized but NOT orthogonalized")')
     !IF (gamma_only) CALL errore('lr_orthoUwfc', &
     !     'Gamma-only calculation for this case not implemented', 1 )
  ELSE
     !WRITE(stdout,*) "U_projection_type =", U_projection
     CALL errore ("hp_sphi"," This U_projection_type is not valid",1)
  ENDIF
  !
  ALLOCATE (wfcatom(npwx*npol,natomwfc))
  ALLOCATE (swfcatom(npwx*npol,natomwfc))
  !
  IF (okvan) CALL allocate_bec_type (nkb,natomwfc,becp)
  !
  !DO ik = 1, nksq
     !
  !   ikk  = ikks(ik)
  !   ikq  = ikqs(ik)
  !   npw  = ngk(ikk)
  !   npwq = ngk(ikq)
     !
     wfcatom  = (0.d0, 0.d0)
     swfcatom = (0.d0, 0.d0)
     !
     ! Determine the atomic orbital at k : phi(k)
     !
     IF (noncolin) THEN
       CALL pw_atomic_wfc_nc_updown (xk, npw, igk_kk, wfcatom)
     ELSE
       CALL pw_atomic_wfc (xk, npw, igk_kk, wfcatom)
     ENDIF
     !
     ! Compute S(k)*phi(k) (phi means the atomic orbital)
     !
     CALL s_phi (npw, vkb_kk, becp, wfcatom, swfcatom)
     !
     ! Orthonormalize or normalize the atomic orbitals (if needed) 
     !
     IF (orthogonalize_wfc) &
        CALL ortho_swfc (npw, normalize_only, natomwfc, wfcatom, swfcatom, lflag)
     !
     ! If lflag=.TRUE. copy the result from (orthonormalized) wfcatom 
     ! (which uses the offset oatwfc) to wfcU (which uses the offset offsetU),
     ! and then write wfcU = phi(k) to file with unit iuatwfc.
     !
     !IF (lflag) THEN
     !   wfcU = (0.d0, 0.d0)
        wfcatomk = (0.d0, 0.d0)
        CALL pw_copy_U_wfc (wfcatom, wfcatomk, noncolin)
     !   CALL save_buffer (wfcU, nwordwfcU, iuatwfc, ikk)
     !ENDIF
     !
     ! Copy the result from (orthonormalized) swfcatom 
     ! (which uses the offset oatwfc) to wfcU (which uses the offset offsetU),
     ! and then write wfcU = S(k)*phi(k) to file with unit iuatswfc.
     !
     !wfcU = (0.d0, 0.d0)
     swfcatomk = (0.d0, 0.d0)
     CALL pw_copy_U_wfc (swfcatom, swfcatomk, noncolin)
     !CALL save_buffer (wfcU, nwordwfcU, iuatswfc, ikk)
     !
     !IF (.NOT.lgamma) THEN
        !
        wfcatom  = (0.d0, 0.d0)
        swfcatom = (0.d0, 0.d0)
        !
        ! Determine the atomic orbital at k+q : phi(k+q)
        !
        IF (noncolin) THEN
           CALL pw_atomic_wfc_nc_updown (xkq, npwq, igk_kq, wfcatom)
        ELSE 
           CALL pw_atomic_wfc (xkq, npwq, igk_kq, wfcatom)
        ENDIF
        !
        ! Compute S(k+q)*phi(k+q) 
        !
        CALL s_phi (npwq, vkb_kq, becp, wfcatom, swfcatom)
        !
        ! Orthonormalize or normalize the atomic orbitals (if needed)
        !
        IF (orthogonalize_wfc) &
           CALL ortho_swfc (npwq, normalize_only, natomwfc, wfcatom, swfcatom, lflag)
        !
        ! If lflag=.TRUE. copy the result from (orthonormalized) wfcatom 
        ! (which uses the offset oatwfc) to wfcU (which uses the offset offsetU),
        ! and then write wfcU = phi(k) to file with unit iuatwfc.
        !
        !IF (lflag) THEN
        !   wfcU = (0.d0, 0.d0)
           wfcatomkpq = (0.d0, 0.d0)
           CALL pw_copy_U_wfc (wfcatom, wfcatomkpq, noncolin)
        !   CALL save_buffer (wfcU, nwordwfcU, iuatwfc, ikq)
        !ENDIF
        !
        ! Copy the result from (orthonormalized) swfcatom 
        ! (which uses the offset oatwfc) to wfcU (which uses the offset offsetU),
        ! and then write wfcU = S(k+q)*phi(k+q) to file with unit iuatswfc.
        !
        !wfcU = (0.d0, 0.d0)
        swfcatomkpq = (0.d0, 0.d0)
        CALL pw_copy_U_wfc (swfcatom, swfcatomkpq, noncolin)
        !CALL save_buffer (wfcU, nwordwfcU, iuatswfc, ikq)
        !
     !ENDIF
     !
  !ENDDO  
  !
  DEALLOCATE (wfcatom)
  DEALLOCATE (swfcatom)
  ! 
  IF (okvan) CALL deallocate_bec_type (becp)
  !
  CALL stop_clock ('lr_orthoUwfc')
  !
  RETURN
  !
CONTAINS
  !
SUBROUTINE s_phi (npw_, vkb_, becp_, wfc, swfc)
  !-----------------------------------------------------------------------
  !
  ! NCPP: swfc = wfc
  ! USPP: swfc = S * wfc
  !
  USE kinds,          ONLY : DP
  USE becmod,         ONLY : bec_type, calbec
  !
  IMPLICIT NONE
  INTEGER,     INTENT(IN)  :: npw_
  !INTEGER,     INTENT(IN)  :: ik_
  COMPLEX(DP), INTENT(IN)  :: wfc  (npwx, natomwfc)
  COMPLEX(DP), INTENT(OUT) :: swfc (npwx, natomwfc)
  !
  complex(dp), intent(in) :: vkb_(npwx, nkb)
  type(bec_type), INTENT(inout) :: becp_
  !
  ! NCPP case
  !
  IF ( nkb == 0 .OR. .NOT. okvan ) THEN
     swfc = wfc
     RETURN
  ENDIF
  !
  ! USPP case
  !
  ! Compute beta functions vkb at ik_
  !
  !CALL init_us_2 (npw_, igk_k(1,ik_), xk(1,ik_), vkb)
  !
  ! Compute the product of beta functions vkb
  ! with the functions wfc : becp = <vkb|wfc>
  !
  CALL calbec (npw_, vkb_, wfc, becp_)
  !
  ! Calculate S*|wfc> = |wfc> + \sum qq * |vkb> * becp 
  !  
  CALL pw_s_psi (npwx, npw_, natomwfc, vkb_, becp_, wfc, swfc)
  !
  RETURN
  !
END SUBROUTINE s_phi

!-------------------------------------------------------------
! adapted from PW/src/ldaU.f90/copy_U_wfc
!
SUBROUTINE pw_copy_U_wfc( swfcatom_, wfcU_, noncolin_ )
  USE ions_base,    ONLY : nat, ntyp => nsp, ityp
  !
  !  Copy (orthogonalized) atomic wavefunctions "swfcatom"
  !  having a Hubbard U correction to array "wfcU"
  !
  IMPLICIT NONE
  COMPLEX(KIND=DP), INTENT(IN) :: swfcatom_(:,:)
  LOGICAL, INTENT(IN), OPTIONAL :: noncolin_
  ! nwfcU is the same as natomwfc
  COMPLEX(KIND=DP), INTENT(out) :: wfcU_(npwx*npol, nwfcU) 
  !
  LOGICAL :: twice
  INTEGER :: na, nt, m1, m2

  IF ( PRESENT(noncolin_) ) THEN
     twice = noncolin_
  ELSE
     twice = .FALSE.
  ENDIF
  DO na=1,nat
     nt = ityp(na)
     if ( is_hubbard(nt) ) THEN
        m1 = 1
        m2 = 2*hubbard_l(nt)+1
        IF ( twice ) m2 = 2*m2
        wfcU_(:,offsetU(na)+m1:offsetU(na)+m2) = swfcatom_(:,oatwfc(na)+m1:oatwfc(na)+m2)
     ENDIF
  ENDDO
  !
END SUBROUTINE pw_copy_U_wfc
  
END SUBROUTINE ph_lr_orthoUwfc

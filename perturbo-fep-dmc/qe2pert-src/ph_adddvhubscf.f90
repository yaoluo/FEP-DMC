!                                         
! Copyright (C) 2001-2018 Quantum ESPRESSO
! This file is distributed under the terms
! GNU General Public License. See the file
! in the root directory of the present dis
! or http://www.gnu.org/copyleft/gpl.txt .
!
! adapted from QE/Phonon/PH/adddvhubscf.90
!----------------------------------------------------------------------------------------
SUBROUTINE ph_adddvhubscf(ikk, npwq, swfcatomk, swfcatomkpq, dnsscf, dvhubscf)
  !--------------------------------------------------------------------------------------
  ! 
  ! DFPT+U 
  ! This routine calculates the SCF derivative of the Hubbard potential times psi.
  ! is  = current_spin
  ! isi = opposite of the current_spin
  !
  ! |Delta V_SCF_(k+q,is) psi(ibnd,k,is)> = 
  !           - \sum_(I,m1,m2) Hubbard_U(I) * dnsscf(m1,m2,is,I,imode) * 
  !                             |S\phi(I,k+q,m1)><S\phi(I,k,m2)|psi(ibnd,k,is)>
  !
  ! Addition of the J0 terms:
  !
  ! + \sum_(I,m1,m2) Hubbard_J0(I) * dnsscf(m1,m2,isi,I,imode) * 
  !         |S\phi(I,k+q,m1)><S\phi(I,k,m2)|psi(ibnd,k,is)>
  !
  !
  ! Written by A. Floris
  ! Modified by I. Timrov (01.10.2018)
  !
  USE kinds,         ONLY : DP
  !USE io_files,      ONLY : nwordwfcU
  USE ions_base,     ONLY : nat, ityp, ntyp => nsp
  USE ldaU,          ONLY : Hubbard_lmax, Hubbard_l, offsetU, is_hubbard, &
                            Hubbard_J0, nwfcU
  !USE ldaU_ph,       ONLY : dnsscf, swfcatomk, swfcatomkpq, proj1, effU
  !USE wvfct,         ONLY : npwx, nbnd 
  !USE control_lr,    ONLY : lgamma, nbnd_occ
  USE lsda_mod,      ONLY : lsda, nspin !, isk
  !USE wavefunctions, ONLY : evc
  !USE eqv,           ONLY : dvpsi  
  !USE mp_bands,      ONLY : intra_bgrp_comm
  !USE mp,            ONLY : mp_sum
  !USE klist,         ONLY : ngk
  !USE buffers,       ONLY : get_buffer
  !USE qpoint,        ONLY : xq, ikks, ikqs
  !USE units_lr,      ONLY : iuatswfc
  !
  use pwcom, only: npwx
  use dvhub_mod, only: ldim, effU
  use electronic_data, only: ngk_tot, nb_sub, spin_updn, evc_sub
  !
  IMPLICIT NONE
  !
  !INTEGER, INTENT(IN) :: ik, ipert  
  ! the k point under consideration
  ! the index of perturbation
  !
  integer, intent(in) :: ikk, npwq
  ! atomic wfc's at k
  complex(DP), intent(in) :: swfcatomk(npwx, nwfcU),     & ! S * atomic wfc at k
                             swfcatomkpq(npwx, nwfcU),   & ! S * atomic wfc at k+q
                             dnsscf(ldim, ldim, nspin, nat)       ! SCF derivative of ns (one perturbation)
  complex(dp), intent(out) :: dvhubscf(npwx, nb_sub)
  !  
  ! Local variables
  !
  !INTEGER :: npw, npwq, ikk, ikq, op_spin
  INTEGER :: npw, op_spin, current_spin
  INTEGER ::  i, j, k, nt, l, ih, n, ig, ihubst, ihubst1, ihubst2, &
              nah, m, m1, m2, ibnd
  COMPLEX(DP), ALLOCATABLE :: proj1(:,:), dvqi(:,:) !dvhubscf(:,:), 
  COMPLEX(DP), EXTERNAL :: ZDOTC
  !
  !CALL start_clock ( 'adddvhubscf' )
  !
  !
  ALLOCATE (proj1(nb_sub, nwfcU))
  !ALLOCATE (dvhubscf(npwx,nbnd))
  ALLOCATE (dvqi(npwx, nb_sub))
  !  
  !ldim = 2 * Hubbard_lmax + 1
  !
  ! At each iteration dvhubscf is set to zero, 
  ! because it is recomputed through dnsscf 
  !
  dvhubscf(1:npwx, 1:nb_sub) = (0.d0, 0.d0) 
  !
  !ikk = ikks(ik)
  !ikq = ikqs(ik)
  !npw = ngk(ikk)
  !npwq= ngk(ikq)
  !
  npw = ngk_tot(ikk)
  !
  current_spin = spin_updn !isk (ikk)
  IF (lsda) THEN
     IF (current_spin==1) THEN   
        op_spin = 2
     ELSE
        op_spin = 1
     END IF
  ELSE        
     op_spin = 1
  ENDIF
  !
  ! Read S * atomic wfc's at k and k+q on unit iuatswfc
  !  
  !CALL get_buffer (swfcatomk, nwordwfcU, iuatswfc, ikk)
  !IF (.NOT.lgamma) CALL get_buffer (swfcatomkpq, nwordwfcU, iuatswfc, ikq)
  !
  DO nah = 1, nat
     !
     nt = ityp(nah)
     !
     IF (is_hubbard(nt)) THEN
        !   
        DO m = 1, 2*Hubbard_l(nt)+1
           !
           ihubst = offsetU(nah) + m   ! I m index
           !
           ! Calculate proj1(ibnd,ihubst) = < S_{k}\phi_(k,I,m)| psi(inbd,k) >
           !
           DO ibnd = 1, nb_sub
              proj1(ibnd,ihubst) = ZDOTC (npw, swfcatomk(:,ihubst), 1, evc_sub(:,ibnd,ikk), 1)
           ENDDO
           !
        ENDDO
        !
     ENDIF
     !
  ENDDO
  ! 
  !CALL mp_sum (proj1, intra_bgrp_comm)
  !
  DO nah = 1, nat
     !
     nt = ityp(nah)
     !
     dvqi = (0.d0, 0.d0)
     !
     IF (is_hubbard(nt)) THEN
        !
        DO m1 = 1, 2*Hubbard_l(nt)+1
           !
           ihubst1 = offsetU(nah) + m1
           !
           DO m2 = 1, 2*Hubbard_l(nt)+1
              ! 
              ihubst2 = offsetU(nah) + m2
              ! 
              DO ibnd = 1, nb_sub
                 DO ig = 1, npwq
                    dvqi(ig,ibnd) = dvqi(ig,ibnd) - effU(nt) * &
                                    !dnsscf(m1,m2,current_spin,nah,ipert) * &
                                    dnsscf(m1,m2,current_spin,nah) * &
                                    swfcatomkpq(ig,ihubst1) * proj1(ibnd,ihubst2)
                 ENDDO
              ENDDO
              !
           ENDDO
           !
        ENDDO
        !
        DO ibnd = 1, nb_sub
           DO ig = 1, npwq
              dvhubscf(ig,ibnd) = dvhubscf(ig,ibnd) + dvqi(ig,ibnd)
           ENDDO
        ENDDO
        !
     ENDIF
     !
     !  Hubbard_J0
     !
     dvqi = (0.d0, 0.d0)
     !
     IF (Hubbard_J0(nt).NE.0.d0) THEN
        !
        DO m1 = 1, 2*Hubbard_l(nt)+1
           !
           ihubst1 = offsetU(nah) + m1
           !
           DO m2 = 1, 2*Hubbard_l(nt)+1
              !
              ihubst2 = offsetU(nah) + m2
              ! 
              DO ibnd = 1, nb_sub
                 DO ig = 1, npwq
                    ! Notice that sign change here
                    dvqi(ig,ibnd) = dvqi(ig,ibnd) + Hubbard_J0(nt) * &
                                    !dnsscf(m1,m2,op_spin,nah,ipert) * &
                                    dnsscf(m1,m2,op_spin,nah) * &
                                    swfcatomkpq(ig,ihubst1) * proj1(ibnd,ihubst2)
                 ENDDO
              ENDDO
              !
           ENDDO
           !
        ENDDO
        !
        DO ibnd = 1, nb_sub
           DO ig = 1, npwq
              dvhubscf(ig,ibnd) = dvhubscf(ig,ibnd) + dvqi(ig,ibnd)
           ENDDO
        ENDDO
        !
     ENDIF
     !
  ENDDO
  !  
  ! Add the result to dvpsi
  !
  !DO ibnd = 1, nb_sub
  !   DO ig = 1, npwq
  !      dvpsi(ig,ibnd) = dvpsi(ig,ibnd) + dvhubscf(ig,ibnd)  
  !   ENDDO
  !ENDDO
  !
  DEALLOCATE (proj1)
  !DEALLOCATE (dvhubscf)
  DEALLOCATE (dvqi)
  !
  !CALL stop_clock ( 'adddvhubscf' )
  !  
  RETURN
  ! 
END SUBROUTINE ph_adddvhubscf

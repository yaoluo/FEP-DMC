!                                         
! Copyright (C) 2001-2018 Quantum ESPRESSO
! This file is distributed under the terms
! GNU General Public License. See the file
! in the root directory of the present dis
! or http://www.gnu.org/copyleft/gpl.txt .
!
! adapted from QE/PHonon/PH/dvqhub_barepsi_us.f90
!-----------------------------------------------------------------------
!SUBROUTINE ph_dvqhub_barepsi_us (ik, uact)
SUBROUTINE ph_dvqhub_barepsi_us(ikk, xk, xkq, npwq, igk_kq, vkb_, vkbkpq, &
      wfcatomk, swfcatomk, wfcatomkpq, swfcatomkpq, dnsbare, icart, na, dvqhbar)
  !-----------------------------------------------------------------------
  !
  ! DFPT+U 
  ! This routines calculates the BARE derivative of the Hubbard potential times psi.
  ! is  = current_spin
  ! isi = opposite of the current_spin 
  !
  ! |Delta V_BARE_(k+q,is) psi(ibnd,k,is)> = 
  !     + \sum_(I,m1,m2) Hubbard_U(I) * [0.5\delta_(m1,m2) - ns(m1,m2,is,I)] *
  !                         { |dqsphi(imode,I,k+q,m1)><S\phi(I,k,m2)|psi(ibnd,k,is)> + 
  !                           |S\phi(I,k+q,m1)><dmqsphi(imode,I,k,m2)|psi(ibnd,k,is)> } 
  !     - \sum_(I,m1,m2) Hubbard_U(I) * dnsbare(m1,m2,is,I,imode) *
  !                           |S\phi(I,k+q,m1)><S\phi(I,k,m2)|psi(ibnd,k,is)>
  !
  ! Addition of the J0 terms:
  !
  !     + \sum_(I,m1,m2) Hubbard_J0(I) * ns(m1,m2,isi,I) *
  !                         { |dqsphi(imode,I,k+q,m1)><S\phi(I,k,m2)|psi(ibnd,k,is)> + 
  !                           |S\phi(I,k+q,m1)><dmqsphi(imode,I,k,m2)|psi(ibnd,k,is)> } 
  !     + \sum_(I,m1,m2) Hubbard_J0(I) * dnsbare(m1,m2,isi,I,imode) *
  !                           |S\phi(I,k+q,m1)><S\phi(I,k,m2)|psi(ibnd,k,is)>
  !
  ! Important: in this routine vkb is a beta function at k+q, and vkb_ is beta at k.
  ! This is done so because vkb is calculated at k+q in solve_linter (i.e. before calling
  ! this routine), so we keep the same attribution here.
  ! 
  ! Written by A. Floris
  ! Modified by I. Timrov (01.10.2018)
  !                 
  USE kinds,         ONLY : DP
  !USE io_global,     ONLY : stdout, ionode
  !USE io_files,      ONLY : nwordwfcU
  USE ions_base,     ONLY : nat, ityp, ntyp => nsp
  !USE klist,         ONLY : xk, ngk, igk_k
  USE ldaU,          ONLY : U_projection, Hubbard_l, is_hubbard, Hubbard_J0, offsetU, nwfcU
  !USE ldaU_ph,       ONLY : wfcatomk, wfcatomkpq, swfcatomk, swfcatomkpq, dwfcatomkpq, &
  !                          sdwfcatomk, sdwfcatomkpq, dvkb, vkbkpq, dvkbkpq, &
  !                          proj1, proj2, dnsbare, effU 
  !USE wvfct,         ONLY : npwx, nbnd
  USE uspp,          ONLY : nkb !, vkb
  !USE qpoint,        ONLY : nksq, ikks, ikqs
  !USE control_lr,    ONLY : lgamma, ofsbeta
  !USE units_lr,      ONLY : iuatwfc, iuatswfc
  USE uspp_param,    ONLY : nh
  USE lsda_mod,      ONLY : lsda, nspin !,isk , current_spin
  !USE wavefunctions, ONLY : evc
  !USE eqv,           ONLY : dvpsi
  USE scf,           ONLY : rho
  !USE mp_bands,      ONLY : intra_bgrp_comm       
  !USE mp,            ONLY : mp_sum 
  !USE buffers,       ONLY : get_buffer
  !
  !
  USE becmod,           ONLY : allocate_bec_type, deallocate_bec_type, bec_type
  !
  use pwcom, only: npwx
  use dvhub_mod, only: ldim, ofsbeta, effU
  use electronic_data, only: spin_updn, ngk_tot, igk_k_tot, nb_sub, evc_sub
  
  IMPLICIT NONE
  !
  !INTEGER, INTENT(IN) :: ik
  ! the k point under consideration
  !COMPLEX(DP), INTENT(IN) :: uact(3*nat)
  ! the pattern of displacements
  !
  integer, intent(in) :: ikk, npwq, icart, na
  integer, intent(in) :: igk_kq(npwx)
  real(dp), intent(in) :: xk(3), xkq(3) !xk and xk+xq in cartesian coordinate
  ! vkb_: the beta function at k, vkbkpq: the beta function at k+q.
  complex(dp), intent(in) :: vkb_(npwx,nkb), vkbkpq(npwx,nkb)
  !! atomic wfc at k and k+q
  complex(dp), intent(in) :: wfcatomk(npwx,nwfcU), wfcatomkpq(npwx,nwfcU)
  !! S * atomic wfc at k and k+q
  complex(dp), intent(in) :: swfcatomk(npwx,nwfcU), swfcatomkpq(npwx,nwfcU)
  !! bare derivative of ns for displacement of na along icart-direction
  complex(dp), intent(in) :: dnsbare(ldim,ldim,nspin,nat)
  !
  complex(dp), intent(out) :: dvqhbar(npwx,nb_sub)  !dvqhbar for na-atom along icart-direction
  !
  ! Local variables
  !
  type(bec_type) :: becp_tmp   !jjzhou, workspace
  integer :: current_spin
  INTEGER :: i, j, k, counter, nt, l, ih, n, mu, ig, &
             ihubst, ihubst1, ihubst2, nah, m, m1, m2, ibnd, op_spin, npw, ibeta
  complex(dp) :: aux1, aux2, aux3, aux4, aux5
  COMPLEX(DP), ALLOCATABLE :: dqsphi(:,:), dmqsphi(:,:), dvqi(:,:), dwfcatom_(:), &
     proj1(:,:), proj2(:,:), dvkb(:,:), dvkbkpq(:,:), sdwfcatomk(:,:), sdwfcatomkpq(:,:)

  COMPLEX(DP), EXTERNAL :: ZDOTC
  !
  allocate( dvkb(npwx, nkb) )
  allocate( dvkbkpq(npwx, nkb) )
  allocate( sdwfcatomk(npwx,nwfcU) )
  allocate( sdwfcatomkpq(npwx,nwfcU) )
  !
  ALLOCATE (proj1(nb_sub,nwfcU))
  ALLOCATE (proj2(nb_sub,nwfcU))
  !ALLOCATE (aux1(npwx))
  !ALLOCATE (aux2(npwx))
  !ALLOCATE (aux3(npwx))
  !ALLOCATE (aux4(npwx)) 
  !ALLOCATE (aux5(npwx))
  ALLOCATE (dqsphi(npwx,nwfcU))
  ALLOCATE (dmqsphi(npwx,nwfcU))
  ALLOCATE (dvqi(npwx,nb_sub))
  !ALLOCATE (dvqhbar(npwx,nbnd,3,nat))
  !ALLOCATE (vkb_(npwx,nkb))
  ALLOCATE (dwfcatom_(npwx))
  !
  call allocate_bec_type(nkb, 1, becp_tmp)
  !
  dvkb         = (0.0d0, 0.0d0)
  dvkbkpq      = (0.0d0, 0.0d0)
  sdwfcatomk   = (0.0d0, 0.0d0)
  sdwfcatomkpq = (0.0d0, 0.0d0)
  dqsphi       = (0.0d0, 0.0d0)
  dmqsphi      = (0.0d0, 0.0d0)
  dvqi         = (0.0d0, 0.0d0)
  dwfcatom_    = (0.0d0, 0.0d0)
  !
  proj1 = (0.d0, 0.d0)
  proj2 = (0.d0, 0.d0)
  !
  !ikk = ikks(ik)
  !ikq = ikqs(ik)
  !npw = ngk(ikk)
  !npwq= ngk(ikq)
  !
  npw = ngk_tot(ikk)
  !
  current_spin = spin_updn !isk(ikk)
  IF (lsda) THEN 
     IF (current_spin==1) THEN
        op_spin = 2
     ELSE
        op_spin = 1
     ENDIF
  ELSE        
     op_spin = 1
  ENDIF
  !
  ! Compute the beta function at k and put the result in vkb_
  !
  !IF (.NOT.lgamma) THEN   
  !   CALL pw_init_us_2 (npw, igk_k_tot(1,ikk), xk, vkb_)
  !ELSE
  !   vkb_ = vkb
  !ENDIF
  !
  ! The beta function at k+q. Let us put it in the proper array vkbkpq,
  ! because in the following of this routine the array vkb will be 
  ! used as a workspace in the routine swfc. 
  !
  !vkbkpq = vkb
  !
  ! Calculate the derivatives of beta functions 
  ! d^{icart}beta at k and k+q for all the bands and for 
  ! the 3 cartesian directions
  !
  !#commented out icart and na loop by jjz,
  !  icart and na are specified in the input arguments.
  !#commented out 'lgamma' by jjz, 'lgamma' was used to save some calculations
  ! if xq = Gamma, here I choose to do some extra-work, but treat all q equal. 
  !
  !DO icart = 1, 3
  !   DO na = 1, nat
        nt = ityp(na) 
        DO ih = 1, nh(nt)
           !
           ibeta = ofsbeta(na) + ih
           !
           CALL ph_dwfc (npw, igk_k_tot(1,ikk), xk, icart, &
                      vkb_(:,ibeta), dvkb(:,ibeta))
           !IF (.NOT.lgamma) &
           CALL ph_dwfc (npwq, igk_kq(1), xkq, icart, &
                      vkbkpq(:,ibeta), dvkbkpq(:,ibeta))
           !
        ENDDO
  !   ENDDO
  !ENDDO
  !
  ! Read \phi at k and k+q from file (unit iuatwfc)
  ! 
  !CALL get_buffer (wfcatomk, nwordwfcU, iuatwfc, ikk)
  !IF (.NOT.lgamma) CALL get_buffer (wfcatomkpq, nwordwfcU, iuatwfc, ikq)
  !
  ! Read S*\phi at k and k+q from file (unit iuatswfc)
  !
  !CALL get_buffer (swfcatomk, nwordwfcU, iuatswfc, ikk)
  !IF (.NOT.lgamma) CALL get_buffer (swfcatomkpq, nwordwfcU, iuatswfc, ikq)
  ! 
  dvqhbar(1:npwx, 1:nb_sub) = (0.d0, 0.d0)
  !
  ! commented out icart and na loop by jjz,
  !  icart and na are specified in the input arguments.
  !
  !DO na = 1, nat
     !
  !   DO icart = 1, 3 
        !  
        dqsphi  = (0.d0, 0.d0)   
        dmqsphi = (0.d0, 0.d0)
        !
        DO nah = 1, nat
           !
           nt = ityp(nah)
           !
           ! For Hubbard_U - Hubbard_J0
           ! 
           IF (is_hubbard(nt)) THEN
              !
              DO m = 1, 2*Hubbard_l(nt)+1
                 !
                 ihubst = offsetU(nah) + m   ! I m index
                 !
                 IF (nah==na) THEN
                    !
                    ! Calculate | d_icart\phi_(k,I,m)) >
                    !
                    CALL ph_dwfc(npw, igk_k_tot(1,ikk), xk, icart, &
                               wfcatomk(:,ihubst), dwfcatom_) 
                    !
                    ! Apply the S matrix: | S d_^(I,icart)\phi_(k,I,m) >
                    !
                    CALL ph_swfc (npw, 1, vkb_, becp_tmp, dwfcatom_, sdwfcatomk(:,ihubst))
                    ! 
                    !IF (.NOT.lgamma) THEN
                       !
                       ! Calculate |d_icart\phi_(k+q,I,m)) >
                       !
                       CALL ph_dwfc(npwq, igk_kq(1), xkq, icart, &
                                  wfcatomkpq(:,ihubst), dwfcatom_)
                       ! 
                       ! Calculate | S d_^(I,icart)\phi_(k+q,I,m) >
                       !
                       CALL ph_swfc (npwq, 1, vkbkpq, becp_tmp, dwfcatom_, sdwfcatomkpq(:,ihubst))
                       ! 
                    !ENDIF
                    ! 
                 ENDIF
                 !
                 ! Calculate |\Delta_q(S_k \phi_(k,I,m)) >  
                 ! and |\Delta_{-q}(S_{k+q} \phi_(k+q,I,m)) >
                 !
                 CALL ph_delta_sphi (npw, npwq, na, icart, nah, ihubst, wfcatomk, wfcatomkpq,  &
                                  sdwfcatomk, sdwfcatomkpq, vkb_, vkbkpq, dvkb(:,:), &
                                  dvkbkpq(:,:), dqsphi, dmqsphi, 1)  
                 !
                 ! Calculate:
                 ! proj1 (ihubst,ibnd) = < S_{k}\phi_(k,I,m)| psi(ibnd,k) >
                 ! proj2 (ihubst,ibnd) = < \Delta_{-q}(S_{k+q} \phi_(k+q,I,m)) | psi(ibnd,k) > 
                 !
                 DO ibnd = 1, nb_sub
                    proj1(ibnd,ihubst) = ZDOTC (npw, swfcatomk(:,ihubst), 1, evc_sub(:,ibnd,ikk), 1)
                    proj2(ibnd,ihubst) = ZDOTC (npw, dmqsphi(:,ihubst), 1, evc_sub(:,ibnd,ikk), 1)
                 ENDDO
                 !
              ENDDO ! m
              !
           ENDIF
           ! 
        ENDDO ! nah
        !
        deallocate( dvkb )
        deallocate( dvkbkpq )
        deallocate( sdwfcatomk )
        deallocate( sdwfcatomkpq )
        deallocate( dwfcatom_ )
        !
        ! commented out by jjz, perturbo does not support bgrp.
        !
        !CALL mp_sum (proj1, intra_bgrp_comm) 
        !CALL mp_sum (proj2, intra_bgrp_comm)
        !
        DO nah = 1, nat
           !
           nt = ityp(nah)
           !
           dvqi = (0.d0, 0.d0) 
           !
           IF (is_hubbard(nt)) THEN
              !
              DO ibnd = 1, nb_sub
                 !
                 DO m1 = 1, 2*Hubbard_l(nt)+1
                    !
                    ihubst1 = offsetU(nah) + m1 
                    !
                    DO ig = 1, npwq
                       !
                       aux1 = dqsphi(ig,ihubst1) * proj1(ibnd,ihubst1) 
                       !
                       aux3 = swfcatomkpq(ig,ihubst1) * proj2(ibnd,ihubst1)  
                       !
                       dvqi(ig,ibnd) = dvqi(ig,ibnd) + 0.5d0 * (aux1 + aux3)
                       !
                    ENDDO
                    !
                    DO m2 = 1, 2*Hubbard_l(nt)+1
                       !
                       ihubst2 = offsetU(nah) + m2
                       !
                       DO ig = 1, npwq 
                          !                         
                          aux2 = dqsphi(ig,ihubst1) * rho%ns(m1,m2,current_spin,nah) &
                                     * proj1(ibnd, ihubst2)
                          aux4 = swfcatomkpq(ig,ihubst1) * rho%ns(m1,m2,current_spin,nah) &
                                     * proj2(ibnd, ihubst2)
                          aux5 = swfcatomkpq(ig,ihubst1) &  !* dnsbare(m1,m2,current_spin,nah,icart,na) &
                                     * dnsbare(m1,m2,current_spin,nah) &
                                     * proj1(ibnd, ihubst2)
                          !
                          dvqi(ig, ibnd) = dvqi(ig, ibnd) - aux2 - aux4 - aux5
                          ! 
                       ENDDO
                       ! 
                    ENDDO ! m2
                    !
                 ENDDO ! m1
                 !
              ENDDO ! ibnd
              !
              ! effU = Hubbard_U - Hubbard_J0
              !
              dvqi = dvqi * effU(nt)
              ! 
              DO ig = 1, npwq
                 !dvqhbar(ig,:,icart,na) = dvqhbar(ig,:,icart,na) + dvqi(ig,:)
                 dvqhbar(ig,1:nb_sub) = dvqhbar(ig,1:nb_sub) + dvqi(ig,1:nb_sub)
              ENDDO
              ! 
           ENDIF
           !
           ! Hubbard_J0 \= 0 case 
           !
           dvqi = (0.d0, 0.d0) 
           ! 
           IF (Hubbard_J0(nt).NE.0.d0) THEN
              !              
              DO ibnd = 1, nb_sub
                 ! 
                 DO m1 = 1, 2*Hubbard_l(nt)+1
                    ! 
                    ihubst1 = offsetU(nah) + m1 
                    !
                    ! No diagonal term for J0
                    !                          
                    DO m2 = 1, 2*Hubbard_l(nt)+1
                       ! 
                       ihubst2 = offsetU(nah) + m2
                       ! 
                       DO ig = 1, npwq                          
                          aux2 = dqsphi(ig, ihubst1) * rho%ns(m1,m2,op_spin,nah) &
                                     * proj1(ibnd, ihubst2)
                          aux4 = swfcatomkpq(ig,ihubst1) * rho%ns(m1,m2,op_spin,nah) &
                                     * proj2(ibnd, ihubst2)
                          aux5 = swfcatomkpq(ig,ihubst1) &  !* dnsbare (m1,m2,op_spin,nah,icart,na) &
                                     * dnsbare (m1,m2,op_spin,nah) &
                                     * proj1 (ibnd, ihubst2)
                          !
                          ! Note the sign change w.r.t. the case above
                          !
                          dvqi(ig, ibnd) = dvqi(ig, ibnd) + aux2 + aux4 + aux5
                          ! 
                       ENDDO
                       !
                    ENDDO ! m2
                    !
                 ENDDO ! m1
                 !
              ENDDO ! ibnd
              !
              dvqi = dvqi * Hubbard_J0(nt)
              !   
              DO ig = 1, npwq
                 !dvqhbar(ig,:,icart,na) = dvqhbar(ig,:,icart,na) + dvqi(ig,:)
                 dvqhbar(ig,1:nb_sub) = dvqhbar(ig,1:nb_sub) + dvqi(ig,1:nb_sub)
              ENDDO
              !
           ENDIF
           !
        ENDDO ! nah
        !
  !   ENDDO ! icart
  !   ! 
  !ENDDO ! na
  !
  !!commented out by jjz, output  dvqhbar in cartesian coordinate
  !
  !! Compute the displacement along the pattern uact (for each band).
  !! The result is stored in aux1.
  !
  !DO ibnd = 1, nbnd
  !   !   
  !   aux1 = (0.d0, 0.d0)
  !   !
  !   DO na = 1, nat
  !      mu = 3 * (na - 1)
  !      ! Here is the basis transformation from cartesian to pattern uact
  !      DO icart = 1, 3 
  !         DO ig = 1, npwq
  !            aux1(ig) = aux1(ig) + dvqhbar(ig,ibnd,icart,na) * uact(mu+icart) 
  !         ENDDO
  !      ENDDO
  !   ENDDO
  !   !
  !   ! Add the result to dvpsi
  !   !
  !   DO ig = 1, npwq
  !      dvpsi(ig,ibnd) = dvpsi(ig,ibnd) + aux1(ig)
  !   ENDDO
  !   !
  !ENDDO 
  !
  DEALLOCATE (proj1)
  DEALLOCATE (proj2)  
  !DEALLOCATE (aux1) 
  !DEALLOCATE (aux2)
  !DEALLOCATE (aux3) 
  !DEALLOCATE (aux4)
  !DEALLOCATE (aux5)
  DEALLOCATE (dqsphi)
  DEALLOCATE (dmqsphi)
  DEALLOCATE (dvqi)
  !DEALLOCATE (dvqhbar)
  !DEALLOCATE (vkb_)
  !DEALLOCATE (dwfcatom_) 
  call deallocate_bec_type( becp_tmp )
  !
  RETURN
  !
END SUBROUTINE ph_dvqhub_barepsi_us

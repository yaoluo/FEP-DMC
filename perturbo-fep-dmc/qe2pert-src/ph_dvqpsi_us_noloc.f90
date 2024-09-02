!
! Copyright (C) 2001-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! adapted from QE/PHonon/PH/dvqpsi_us_only.f90
!----------------------------------------------------------------------
subroutine ph_dvqpsi_us_noloc (ikk, xkq, npwq, igk_kq, uact, vkb, dvpsi_noloc)
  !----------------------------------------------------------------------
  !
  ! This routine calculates dV_bare/dtau * psi for one perturbation
  ! with a given q. The displacements are described by a vector uact.
  ! The result is stored in dvpsi. The routine is called for each k point
  ! and for each pattern u. It computes simultaneously all the bands.
  ! This routine implements Eq. B29 of PRB 64, 235118 (2001).
  ! Only the contribution of the nonlocal potential is calculated here.
  !
  ! jjzhou & ite.: remove dependence on module
  ! It computes the 1st and 3rd terms in Eq.B29 of PRB 64, 235118(2001)
  !
  USE kinds, only : dp
  USE cell_base, ONLY : tpiba, bg
  USE gvect,     ONLY : g
  USE ions_base, ONLY : nat, ityp, ntyp => nsp
  USE lsda_mod,  ONLY : lsda, isk, nspin
  USE spin_orb,  ONLY : lspinorb
  USE noncollin_module, ONLY : noncolin, npol
  USE uspp, ONLY: okvan, nkb
  USE uspp_param, ONLY: nh, nhm
   
  use pwcom, only: npwx
  use electronic_data, only: num_band, et_sub, becp1_tot, alphap_tot, spin_updn
  use elph_matrix, only: int1, int2, int1_nc, int2_so


  !USE wvfct,     ONLY : nbnd, npwx, et
  !USE klist,     ONLY : xk, ngk, igk_k
  !USE phus,      ONLY : int1, int1_nc, int2, int2_so, alphap
  !USE lrus,       ONLY : becp1
  !USE qpoint,     ONLY : ikks, ikqs
  !USE eqv,        ONLY : dvpsi

  implicit none
  !
  !   The dummy variables
  !
  integer, intent(in) :: ikk, npwq
  integer, intent(in) :: igk_kq(npwx)
  real(dp), intent(in) :: xkq(3)  !xk+xq in cartesian coordinate
  ! the point k
  ! the point k+q
  complex(DP), intent(in) :: uact (3*nat), vkb(npwx, nkb)
  ! input: the pattern of displacements, vkb: beta(k+q+G)
  !
  complex(dp), intent(out) :: dvpsi_noloc(npwx*npol, num_band)
  !   And the local variables
  !

  integer :: na, nb, mu, nu, ig, igg, nt, ibnd, ijkb0, &
       ikb, jkb, ih, jh, ipol, is, js, ijs, current_spin
  ! counter on atoms
  ! counter on modes
  ! counter on G vectors
  ! auxiliary counter on G vectors
  ! counter on atomic types
  ! counter on bands
  ! auxiliary variable for counting
  ! counter on becp functions
  ! counter on becp functions
  ! counter on n index
  ! counter on m index
  ! counter on polarizations

  real(DP), parameter :: eps = 1.d-12
   
  complex(DP), allocatable :: ps1 (:,:), ps2 (:,:,:), aux (:), deff_nc(:,:,:,:)
  real(DP), allocatable :: deff(:,:,:)
  complex(DP), allocatable :: ps1_nc (:,:,:), ps2_nc (:,:,:,:)
  ! work space

  logical :: ok

!!$omp master
!  call start_clock ('ph_dvqpsi_us_noloc')
!!$omp end master

  if (noncolin) then
     allocate (ps1_nc(nkb , npol, num_band))
     allocate (ps2_nc(nkb , npol, num_band , 3))
     allocate (deff_nc(nhm, nhm, nat, nspin))
  else
     allocate (ps1 ( nkb , num_band))
     allocate (ps2 ( nkb , num_band , 3))
     allocate (deff(nhm, nhm, nat))
  end if
  allocate (aux ( npwx))

  dvpsi_noloc = (0.d0, 0.d0)

  !ikk = ikks(ik)
  !ikq = ikqs(ik)

  !if (lsda) current_spin = isk (ikk)
  current_spin = spin_updn  !@jjzhou

  !
  !   we first compute the coefficients of the vectors
  !
  if (noncolin) then
     ps1_nc(:,:,:)   = (0.d0, 0.d0)
     ps2_nc(:,:,:,:) = (0.d0, 0.d0)
  else
     ps1(:,:)   = (0.d0, 0.d0)
     ps2(:,:,:) = (0.d0, 0.d0)
  end if
  do ibnd = 1, num_band
     IF (noncolin) THEN
        CALL compute_deff_nc(deff_nc, et_sub(ibnd,ikk))
     ELSE
        CALL compute_deff(deff,et_sub(ibnd,ikk))
     ENDIF
     ijkb0 = 0
     do nt = 1, ntyp
        do na = 1, nat
           if (ityp (na) .eq.nt) then
              mu = 3 * (na - 1)
              do ih = 1, nh (nt)
                 ikb = ijkb0 + ih
                 do jh = 1, nh (nt)
                    jkb = ijkb0 + jh
                    do ipol = 1, 3
                       if ( abs (uact (mu + 1) ) + &
                            abs (uact (mu + 2) ) + &
                            abs (uact (mu + 3) ) > eps) then
                          IF (noncolin) THEN
                             ijs=0
                             DO is=1,npol
                                DO js=1,npol
                                   ijs=ijs+1
                                   ps1_nc(ikb,is,ibnd)=ps1_nc(ikb,is,ibnd) +  &
                                      deff_nc(ih,jh,na,ijs) * &
                                      alphap_tot(ipol, ikk)%nc(jkb,js,ibnd)* &
                                       uact(mu + ipol)
                                   ps2_nc(ikb,is,ibnd,ipol)=               &
                                          ps2_nc(ikb,is,ibnd,ipol)+        &
                                          deff_nc(ih,jh,na,ijs) *          &
                                          becp1_tot(ikk)%nc(jkb,js,ibnd) *      &
                                          (0.d0,-1.d0) * uact(mu+ipol) * tpiba
                                END DO
                             END DO
                          ELSE
                             ps1 (ikb, ibnd) = ps1 (ikb, ibnd) +      &
                                        deff(ih, jh, na) *            &
                                alphap_tot(ipol, ikk)%k(jkb, ibnd) * uact (mu + ipol)
                             ps2 (ikb, ibnd, ipol) = ps2 (ikb, ibnd, ipol) +&
                                  deff(ih,jh,na)*becp1_tot(ikk)%k (jkb, ibnd) *  &
                                  (0.0_DP,-1.0_DP) * uact (mu + ipol) * tpiba
                          ENDIF
                          IF (okvan) THEN
                             IF (noncolin) THEN
                                ijs=0
                                DO is=1,npol
                                   DO js=1,npol
                                      ijs=ijs+1
                                      ps1_nc(ikb,is,ibnd)=ps1_nc(ikb,is,ibnd)+ &
                                         int1_nc(ih,jh,ipol,na,ijs) *     &
                                         becp1_tot(ikk)%nc(jkb,js,ibnd)*uact(mu+ipol)
                                   END DO
                                END DO
                             ELSE
                                ps1 (ikb, ibnd) = ps1 (ikb, ibnd) + &
                                  (int1 (ih, jh, ipol,na, current_spin) * &
                                  becp1_tot(ikk)%k (jkb, ibnd) ) * uact (mu +ipol)
                             END IF
                          END IF
                       END IF  ! uact>0
                       if (okvan) then
                          do nb = 1, nat
                             nu = 3 * (nb - 1)
                             IF (noncolin) THEN
                                IF (lspinorb) THEN
                                   ijs=0
                                   DO is=1,npol
                                      DO js=1,npol
                                         ijs=ijs+1
                                         ps1_nc(ikb,is,ibnd)= &
                                                   ps1_nc(ikb,is,ibnd)+ &
                                         int2_so(ih,jh,ipol,nb,na,ijs)* &
                                          becp1_tot(ikk)%nc(jkb,js,ibnd)*uact(nu+ipol)
                                      END DO
                                   END DO
                                ELSE
                                   DO is=1,npol
                                      ps1_nc(ikb,is,ibnd)=ps1_nc(ikb,is,ibnd)+ &
                                         int2(ih,jh,ipol,nb,na) * &
                                         becp1_tot(ikk)%nc(jkb,is,ibnd)*uact(nu+ipol)
                                   END DO
                                END IF
                             ELSE
                                ps1 (ikb, ibnd) = ps1 (ikb, ibnd) + &
                                    (int2 (ih, jh, ipol, nb, na) * &
                                     becp1_tot(ikk)%k (jkb, ibnd) ) * uact (nu + ipol)
                             END IF
                          enddo
                       endif  ! okvan
                    enddo ! ipol
                 enddo ! jh
              enddo ! ih
              ijkb0 = ijkb0 + nh (nt)
           endif
        enddo  ! na
     enddo ! nt
  enddo ! nbnd
  !
  !      This term is proportional to beta(k+q+G)
  ! 
  if (nkb.gt.0) then
     if (noncolin) then
        call zgemm ('N', 'N', npwq, num_band*npol, nkb, &
         (1.d0, 0.d0), vkb, npwx, ps1_nc, nkb, (1.d0, 0.d0) , dvpsi_noloc, npwx)
     else
        call zgemm ('N', 'N', npwq, num_band, nkb, &
         (1.d0, 0.d0) , vkb, npwx, ps1, nkb, (1.d0, 0.d0) , dvpsi_noloc, npwx)
     end if
  end if
  !
  !      This term is proportional to (k+q+G)_\alpha*beta(k+q+G)
  !
  !
  do ikb = 1, nkb
     do ipol = 1, 3
        ok = .false.
        IF (noncolin) THEN
           do ibnd = 1, num_band
              ok = ok.or.(abs (ps2_nc (ikb, 1, ibnd, ipol) ).gt.eps).or. &
                         (abs (ps2_nc (ikb, 2, ibnd, ipol) ).gt.eps)
           end do
        ELSE
           do ibnd = 1, num_band
              ok = ok.or. (abs (ps2 (ikb, ibnd, ipol) ) .gt.eps)
           enddo
        ENDIF
        if (ok) then
           do ig = 1, npwq
              igg = igk_kq (ig)
              aux (ig) =  vkb(ig, ikb) * (xkq(ipol) + g(ipol, igg) )
           enddo
           do ibnd = 1, num_band
              IF (noncolin) THEN
                 call zaxpy(npwq,ps2_nc(ikb,1,ibnd,ipol),aux,1,dvpsi_noloc(1,ibnd),1)
                 call zaxpy(npwq,ps2_nc(ikb,2,ibnd,ipol),aux,1, &
                                                         dvpsi_noloc(1+npwx,ibnd),1)
              ELSE
                 call zaxpy (npwq, ps2(ikb,ibnd,ipol), aux, 1, dvpsi_noloc(1,ibnd), 1)
              END IF
           enddo
        endif
     enddo

  enddo
  deallocate (aux)
  IF (noncolin) THEN
     deallocate (ps2_nc)
     deallocate (ps1_nc)
     deallocate (deff_nc)
  ELSE
     deallocate (ps2)
     deallocate (ps1)
     deallocate (deff)
  END IF

!!$omp master
!  call stop_clock ('ph_dvqpsi_us_noloc')
!!$omp end master

  return
end subroutine ph_dvqpsi_us_noloc

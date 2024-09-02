!
! Copyright (C) 2001-2018 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! adapted from QE/LR_modules/adddvscf.90
!----------------------------------------------------------------------
subroutine ph_adddvscf(npwq, ipert, becp1_ikk, vkb_ikq, dvpsi)
  !----------------------------------------------------------------------
  !
  !     This routine computes the contribution of the selfconsistent
  !     change of the potential to the known part of the linear
  !     system.
  !     It implements the second term in Eq. B30 of PRB 64, 235118 (2001).
  !
  USE kinds,      ONLY : DP
  USE uspp_param, ONLY : upf, nh
  USE uspp,       ONLY : okvan, nkb
! modules from pwcom
  USE lsda_mod,   ONLY : lsda, isk
  USE ions_base,  ONLY : ntyp => nsp, nat, ityp
  USE noncollin_module, ONLY : noncolin, npol
  use becmod, only: bec_type
   
  use elph_matrix, only: int3, int3_nc
  use electronic_data, only: num_band, npwx, spin_updn
  !
!  USE klist,      ONLY : ngk
!  USE wvfct,      ONLY : nbnd, npwx
! modules from phcom
!  USE lrus,       ONLY : int3, int3_nc, becp1
!  USE qpoint,     ONLY : ikks, ikqs
!  USE eqv,        ONLY : dvpsi

  implicit none
  !
  !   The dummy variables
  !
  integer, intent(in) :: npwq, ipert  !ikk, ikq, 
  ! number of the plane-waves at point k+q
  ! input: the perturbation
  type(bec_type), intent(in) :: becp1_ikk
  ! vkb_ikq: vkb(k+q+G) at (k+q)
  complex(dp), intent(in) :: vkb_ikq(npwx, nkb)
  complex(dp), intent(out) :: dvpsi(npwx*npol, num_band)
  !
  !   And the local variables
  !
  integer :: na, nt, ibnd, ih, jh, ijkb0, ikb, jkb, is, js, ijs, current_spin
  ! counter on atoms
  ! counter on atomic types
  ! counter on bands
  ! counter on beta functions
  ! counter on beta functions
  ! auxiliary variable for indexing
  ! counter on the k points
  complex(DP) :: sum, sum_nc(npol)
  ! auxiliary variable

  !always initialize it
  dvpsi = (0.d0, 0.d0)
  !
  if (.not.okvan) return
  
!!$omp master
!  call start_clock ('ph_adddvscf')
!!$omp end master 
  !
  !ikk  = ikks(ik)
  !ikq  = ikqs(ik)
  !npwq = ngk(ikq)
  !
  !if (lsda) current_spin = isk(ikk)
  current_spin = spin_updn
  !
  ijkb0 = 0
  do nt = 1, ntyp
     if (upf(nt)%tvanp  ) then
        do na = 1, nat
           if (ityp (na) .eq.nt) then
              !
              !   we multiply the integral for the becp term and the beta_n
              !
              do ibnd = 1, num_band
                 do ih = 1, nh (nt)
                    ikb = ijkb0 + ih
                    IF (noncolin) THEN
                       sum_nc = (0.d0, 0.d0)
                    ELSE
                       sum = (0.d0, 0.d0)
                    END IF
                    do jh = 1, nh (nt)
                       jkb = ijkb0 + jh
                       IF (noncolin) THEN
                          ijs=0
                          do is=1,npol
                             do js=1,npol
                                ijs=ijs+1
                                sum_nc(is)=sum_nc(is)+               &
                                     int3_nc(ih,jh,na,ijs,ipert)*    &
                                     becp1_ikk%nc(jkb, js, ibnd)
                             enddo
                          enddo
                       ELSE
                          sum = sum + int3 (ih, jh, na, current_spin, ipert)*&
                                   becp1_ikk%k(jkb, ibnd)
                       END IF
                    enddo
                    IF (noncolin) THEN
                       call zaxpy(npwq,sum_nc(1),vkb_ikq(1,ikb),1,dvpsi(1,ibnd),1)
                       call zaxpy(npwq,sum_nc(2),vkb_ikq(1,ikb),1, &
                                                 dvpsi(1+npwx,ibnd),1)
                    ELSE
                       call zaxpy(npwq,sum,vkb_ikq(1,ikb),1,dvpsi(1,ibnd),1)
                    END IF
                 enddo
              enddo
              ijkb0 = ijkb0 + nh (nt)
           endif
        enddo
     else
        do na = 1, nat
           if (ityp (na) .eq.nt) ijkb0 = ijkb0 + nh (nt)
        enddo
     endif
  enddo
  !
!!$omp master
!  call stop_clock ('ph_adddvscf')
!!$omp end master
  !
  return
  !
end subroutine ph_adddvscf

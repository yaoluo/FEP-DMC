!
! Copyright (C) 2006 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! adapted from PHONON/PH/transform_int_nc
!
SUBROUTINE ph_transform_int1_nc(int1, int1_nc, na, iflag)
!----------------------------------------------------------------------------
!
! This routine multiply int1 by the identity and the Pauli
! matrices and saves it in int1_nc.
!
USE kinds,                ONLY : DP
USE ions_base,            ONLY : nat, ityp
USE uspp_param,           ONLY : nh, nhm
USE spin_orb,             ONLY : domag
USE noncollin_module,     ONLY : nspin_mag
USE lsda_mod,             ONLY : nspin
!
!USE phus,                 ONLY : int1_nc
!
IMPLICIT NONE

INTEGER, intent(in) :: na, iflag
COMPLEX(DP), intent(in) :: int1(nhm,nhm,3,nat,nspin_mag)
!
! NOTE: only update int1_nc(:,:,:,na,:)
COMPLEX(DP), intent(inout) :: int1_nc(nhm, nhm, 3, nat, nspin)
!
! ... local variables
!
INTEGER :: ih, jh, ipol, np

np=ityp(na)
DO ih = 1, nh(np)
   DO jh = 1, nh(np)
      DO ipol=1,3
         IF (iflag==0) THEN
            IF (domag) THEN
               int1_nc(ih,jh,ipol,na,1)= &
                                 int1(ih,jh,ipol,na,1)+int1(ih,jh,ipol,na,4)
               int1_nc(ih,jh,ipol,na,2)=                                       &
                  int1(ih,jh,ipol,na,2) - (0.d0, 1.d0) * int1(ih,jh,ipol,na,3)
               int1_nc(ih,jh,ipol,na,3)=                                       &
                  int1(ih,jh,ipol,na,2) + (0.d0, 1.d0) * int1(ih,jh,ipol,na,3)
               int1_nc(ih,jh,ipol,na,4)=                                       &
                  int1(ih,jh,ipol,na,1) - int1(ih,jh,ipol,na,4)
            ELSE
               int1_nc(ih,jh,ipol,na,1)=int1(ih,jh,ipol,na,1)
               int1_nc(ih,jh,ipol,na,4)=int1(ih,jh,ipol,na,1)
            END IF
         ELSE
            IF (domag) THEN
               int1_nc(ih,jh,ipol,na,1)= &
                             CONJG(int1(ih,jh,ipol,na,1)+int1(ih,jh,ipol,na,4))
               int1_nc(ih,jh,ipol,na,2)=CONJG(int1(ih,jh,ipol,na,2)) - &
                           (0.d0, 1.d0)*CONJG(int1(ih,jh,ipol,na,3))
               int1_nc(ih,jh,ipol,na,3)=CONJG(int1(ih,jh,ipol,na,2)) + &
                           (0.d0, 1.d0)*CONJG(int1(ih,jh,ipol,na,3))
               int1_nc(ih,jh,ipol,na,4)=                               &
                  CONJG(int1(ih,jh,ipol,na,1) - int1(ih,jh,ipol,na,4))
            ELSE
               int1_nc(ih,jh,ipol,na,1)=CONJG(int1(ih,jh,ipol,na,1))
               int1_nc(ih,jh,ipol,na,4)=CONJG(int1(ih,jh,ipol,na,1))
            END IF
         END IF
      END DO
   END DO
END DO

RETURN
END SUBROUTINE ph_transform_int1_nc


! adapted from PHONON/PH/transform_int_nc
!
!----------------------------------------------------------------------------
SUBROUTINE ph_transform_int2_nc(int2, int2_so, nb, iflag)
!----------------------------------------------------------------------------
!
! This routine sets int2_so for the atomic species which do not
! have a spin-orbit pseudopotential
!
USE kinds,                ONLY : DP
USE ions_base,            ONLY : nat, ityp
USE uspp_param,           ONLY : nh, nhm
USE lsda_mod,             ONLY : nspin
!
!USE phus,                 ONLY : int2_so
!
IMPLICIT NONE
INTEGER, intent(in) :: nb, iflag
COMPLEX(DP), intent(in) :: int2(nhm,nhm,3,nat,nat)
!
! NOTE: only update int2_so(:,:,:,:,nb,1) and int2_so(:,:,:,:,nb,4)
COMPLEX(DP), intent(inout) :: int2_so(nhm, nhm, 3, nat, nat, nspin)
!
! ... local variables
!
INTEGER :: ih, jh, np, na, ipol

np=ityp(nb)
DO ih = 1, nh(np)
   DO jh = 1, nh(np)
      DO na=1,nat
         DO ipol=1,3
            IF (iflag==0) THEN
               int2_so(ih,jh,ipol,na,nb,1)=int2(ih,jh,ipol,na,nb)
               int2_so(ih,jh,ipol,na,nb,4)=int2(ih,jh,ipol,na,nb)
            ELSE
               int2_so(ih,jh,ipol,na,nb,1)=CONJG(int2(ih,jh,ipol,na,nb))
               int2_so(ih,jh,ipol,na,nb,4)=CONJG(int2(ih,jh,ipol,na,nb))
            END IF
         END DO
      END DO
   END DO
END DO

RETURN
END SUBROUTINE ph_transform_int2_nc


! adapted from PHONON/PH/transform_int_so
!
!----------------------------------------------------------------------------
SUBROUTINE ph_transform_int1_so(int1, int1_nc, na, iflag)
!----------------------------------------------------------------------------
!
! This routine multiply int1 by the identity and the Pauli
! matrices, rotate it as appropriate for the spin-orbit case
! and saves it in int1_nc.
!
USE kinds,                ONLY : DP
USE ions_base,            ONLY : nat, ityp
USE uspp_param,           ONLY : nh, nhm
USE noncollin_module,     ONLY : npol, nspin_mag
USE spin_orb,             ONLY : fcoef, domag
USE lsda_mod,             ONLY : nspin
!
!USE phus,                 ONLY : int1_nc
!
IMPLICIT NONE

INTEGER, intent(in) :: na, iflag
COMPLEX(DP), intent(in) :: int1(nhm,nhm,3,nat,nspin_mag)
!
! Note: only update int1_nc(:,:,:,na,:)
COMPLEX(DP), intent(inout) :: int1_nc(nhm, nhm, 3, nat, nspin)
!
! ... local variables
!
INTEGER :: ih, jh, lh, kh, ipol, np, is1, is2, ijs
COMPLEX(DP) :: fact(4)
LOGICAL :: same_lj

np=ityp(na)
DO ih = 1, nh(np)
   DO kh = 1, nh(np)
      IF (same_lj(kh,ih,np)) THEN
         DO jh = 1, nh(np)
            DO lh= 1, nh(np)
               IF (same_lj(lh,jh,np)) THEN
                  DO ipol=1,3
                     ijs=0
                     DO is1=1,npol
                        DO is2=1,npol
                           ijs=ijs+1
                           IF (iflag==0) THEN
                              fact(1)=int1(kh,lh,ipol,na,1)
                           ELSE
                              fact(1)=CONJG(int1(kh,lh,ipol,na,1))
                           ENDIF
                           int1_nc(ih,jh,ipol,na,ijs)=                       &
                               int1_nc(ih,jh,ipol,na,ijs) +                  &
                               fact(1)*                       &
                             (fcoef(ih,kh,is1,1,np)*fcoef(lh,jh,1,is2,np)  + &
                             fcoef(ih,kh,is1,2,np)*fcoef(lh,jh,2,is2,np)   )
                           IF (domag) THEN
                              IF (iflag==0) THEN
                                 fact(2)=int1 (kh,lh,ipol,na,2)
                                 fact(3)=int1 (kh,lh,ipol,na,3)
                                 fact(4)=int1 (kh,lh,ipol,na,4)
                              ELSE
                                 fact(2)=CONJG(int1 (kh,lh,ipol,na,2))
                                 fact(3)=CONJG(int1 (kh,lh,ipol,na,3))
                                 fact(4)=CONJG(int1 (kh,lh,ipol,na,4))
                              ENDIF
                              int1_nc(ih,jh,ipol,na,ijs)=                     &
                                 int1_nc(ih,jh,ipol,na,ijs) +                 &
                                 fact(2)*                       &
                                (fcoef(ih,kh,is1,1,np)*fcoef(lh,jh,2,is2,np)+ &
                                 fcoef(ih,kh,is1,2,np)*fcoef(lh,jh,1,is2,np))+&
                                 (0.D0,-1.D0) * fact(3)*        &
                                (fcoef(ih,kh,is1,1,np)*fcoef(lh,jh,2,is2,np)- &
                                 fcoef(ih,kh,is1,2,np)*fcoef(lh,jh,1,is2,np))+&
                                 fact(4)*                      &
                                (fcoef(ih,kh,is1,1,np)*fcoef(lh,jh,1,is2,np)- &
                                 fcoef(ih,kh,is1,2,np)*fcoef(lh,jh,2,is2,np))
                           END IF
                        END DO
                     END DO
                  END DO
               END IF
            END DO
         END DO
      END IF
   END DO
END DO
       !
RETURN
END SUBROUTINE ph_transform_int1_so


! adapted from PHONON/PH/transform_int_so
!
!----------------------------------------------------------------------------
SUBROUTINE ph_transform_int2_so(int2, int2_so, nb, iflag)
!----------------------------------------------------------------------------
!
! This routine rotates int2 as appropriate for the spin-orbit case
! and saves it in int2_so.
!
USE kinds,                ONLY : DP
USE ions_base,            ONLY : nat, ityp
USE uspp_param,           ONLY : nh, nhm
USE noncollin_module,     ONLY : npol
USE spin_orb,             ONLY : fcoef
USE lsda_mod,             ONLY : nspin
!
!USE phus,                 ONLY : int2_so
!
IMPLICIT NONE
INTEGER, intent(in) :: nb, iflag
COMPLEX(DP), intent(in) :: int2(nhm,nhm,3,nat,nat)
!
! Note: only update int2_so(:,:,:,:,nb,:)
COMPLEX(DP), intent(inout):: int2_so(nhm, nhm, 3, nat, nat, nspin)
!
! ... local variables
!
INTEGER :: ih, jh, lh, kh, ijs, np, is1, is2, na, ipol
COMPLEX(DP) :: fact
LOGICAL :: same_lj

np=ityp(nb)
DO ih = 1, nh(np)
   DO kh = 1, nh(np)
      IF (same_lj(kh,ih,np)) THEN
         DO jh = 1, nh(np)
            DO lh= 1, nh(np)
               IF (same_lj(lh,jh,np)) THEN
                  DO na=1,nat
                     DO ipol=1,3
                        IF (iflag==0) THEN
                           fact=int2(kh,lh,ipol,na,nb)
                        ELSE
                           fact=CONJG(int2(kh,lh,ipol,na,nb))
                        ENDIF
                        ijs=0
                        DO is1=1,npol
                           DO is2=1,npol
                              ijs=ijs+1
                              int2_so(ih,jh,ipol,na,nb,ijs)= &
                              int2_so(ih,jh,ipol,na,nb,ijs)+ &
                              fact* &
                            (fcoef(ih,kh,is1,1,np)*fcoef(lh,jh,1,is2,np) + &
                             fcoef(ih,kh,is1,2,np)*fcoef(lh,jh,2,is2,np)  )
                           END DO
                        END DO
                     END DO
                  END DO
               END IF
            END DO
         END DO
      END IF
   END DO
END DO
       !
RETURN
END SUBROUTINE ph_transform_int2_so

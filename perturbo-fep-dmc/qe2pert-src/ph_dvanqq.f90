!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! adapted from QE/PHonon/PH/dvanqq
!----------------------------------------------------------------------
subroutine ph_dvanqq(xq, eigqts, vlocq)
  !----------------------------------------------------------------------
  !
  ! This routine calculates four integrals of the Q functions and
  ! its derivatives with V_loc and V_eff which are used
  ! to compute term dV_bare/dtau * psi  in addusdvqpsi and in addusdynmat.
  ! The result is stored in int1,int2,int4,int5. The routine is called
  ! only once. int4 and int5 are deallocated after use in addusdynmat.
  ! int1 -> Eq. B20 of Ref.[1]
  ! int2 -> Eq. B21 of Ref.[1]
  ! int4 -> Eq. B23 of Ref.[1]
  ! int5 -> Eq. B24 of Ref.[1]
  !
  ! [1] PRB 64, 235118 (2001).
  !
  ! jjzhou: adapted from PHonon/PH/dvanqq,
  ! remove dependence on module phus, qpoint, control_ph
  ! modified to compute only: int1, int2, int1_nc, int2_so. (originally in phus)
  !
  USE kinds, only : DP
  USE constants, only: tpi
  USE cell_base, ONLY : omega, tpiba2, tpiba
  USE ions_base, ONLY : nat, ityp, tau, ntyp => nsp
  USE fft_base,   ONLY: dfftp
  USE fft_interfaces, ONLY: fwfft
  use gvect, only : ngm, gg, g, mill, eigts1, eigts2, eigts3, gcutm
  use spin_orb, only : lspinorb
  use scf, only : v, vltot
  use noncollin_module, ONLY : noncolin, nspin_mag
  USE uspp, ONLY: okvan, ijtoh
  USE uspp_param, ONLY: upf, lmaxq, nh, nhm, nbetam, lmaxkb
  USE lsda_mod,      ONLY : nspin
  !USE control_lr, ONLY : lgamma

  USE Coul_cut_2D, ONLY: do_cutoff_2D 
  USE Coul_cut_2D_ph, ONLY: lr_Vlocq

  USE us, only: nqxq, qrad, dq
  USE cellmd, only: cell_factor
  USE elph_matrix, only: int1, int2, int1_nc, int2_so
  
  ! no support for task group, only pools level parallelization is supported!
  !USE mp_bands, ONLY: intra_bgrp_comm
  !USE mp,        ONLY: mp_sum
  
  !USE eqv,        ONLY : vlocq
  !USE phus, ONLY : int1, int2, int4, int4_nc, int5, int5_so
  !USE control_ph, ONLY : rec_code_read
  !-------------------------
  ! setqmod from  LR_Modules, compute (q+G) and |q+G|^2
  ! ylmr2   from  Modules, Real spherical harmonics ylm(G) up to l=lmax
  ! fwfft   from  FFTXlib
  ! qvan2   from  PW, computes the fourier transform of the Q functions

  implicit none
  
  !input: xq_: q coordinates in cartesian coordiante
  real(dp), intent(in) :: xq(3)
  complex(dp), intent(in) :: eigqts(nat)
  ! the local potential at q+G
  real(dp), intent(in) :: vlocq(ngm, ntyp)  ! ngm, ntyp)
  !
  ! And the local variables
  logical :: lgamma
  integer :: nt, na, nb, ig, nta, ntb, ir, ih, jh, ijh, ipol, jpol, is, np
  ! counters
  real(dp) :: qnorm
  real(DP), allocatable :: qmod (:), qmodg (:), qpg (:,:), &
       ylmkq (:,:), ylmk0 (:,:)
  ! the modulus of q+G
  ! the modulus of G
  ! the  q+G vectors
  ! the spherical harmonics

  complex(DP) :: fact, fact1, zdotc
  complex(DP), allocatable :: aux1 (:), aux2 (:),&
       aux5 (:), veff (:,:), sk(:)
  ! work space
  complex(DP), allocatable, target :: qgm(:)
  ! the augmentation function at G
  complex(DP), pointer :: qgmq (:)
  ! the augmentation function at q+G

  if (.not.okvan) return
  !check if xq is Gamma point.
  lgamma = ( sqrt(dot_product(xq,xq)) < 1.0E-12_dp )

  ! reallocate the Fourier coefficients of augmentation functions, qrad, which
  ! depends on q coordinates, for int2 [see (B21), PRB 64, 235118, 2001]
  !   qrad is used in qvan2 to compute the augmentation function.
  qnorm = sqrt(xq(1)**2 + xq(2)**2 + xq(3)**2)
  nqxq = int( ( (sqrt(gcutm) + qnorm) / dq + 4) * cell_factor )
  lmaxq = 2*lmaxkb+1
  if (allocated(qrad)) deallocate(qrad)
  IF (lmaxq > 0) allocate( qrad(nqxq, nbetam*(nbetam+1)/2, lmaxq, ntyp) )
  !initialize qrad in the init_us_1 subroutine
  call init_us_1()
  

  !call start_clock ('ph_dvanqq')
  int1(:,:,:,:,:) = (0.d0, 0.d0)
  int2(:,:,:,:,:) = (0.d0, 0.d0)
   
  allocate (sk  (  ngm))
  allocate (aux1(  ngm))
  allocate (aux2(  ngm))
  allocate (aux5(  ngm))
  allocate (qmodg( ngm))
  allocate (ylmk0( ngm , lmaxq * lmaxq))
  allocate (qgm  ( ngm))
  
  if (.not.lgamma) then
     allocate (ylmkq(ngm , lmaxq * lmaxq))
     allocate (qmod( ngm))
     allocate (qgmq( ngm))
  else
     qgmq =>qgm
  endif
  !
  !     compute spherical harmonics
  !
  call ylmr2 (lmaxq * lmaxq, ngm, g, gg, ylmk0)
  do ig = 1, ngm
     qmodg (ig) = sqrt (gg (ig) )
  enddo
  if (.not.lgamma) then
     allocate (qpg (3, ngm))
     call setqmod (ngm, xq, g, qmod, qpg)
     call ylmr2 (lmaxq * lmaxq, ngm, qpg, qmod, ylmkq)
     deallocate (qpg)
     do ig = 1, ngm
        qmod (ig) = sqrt (qmod (ig) )
     enddo
  endif
  !
  !   we start by computing the FT of the effective potential
  !
  allocate (veff ( dfftp%nnr , nspin_mag))
  do is = 1, nspin_mag
     if (nspin_mag.ne.4.or.is==1) then
        do ir = 1, dfftp%nnr
           veff (ir, is) = CMPLX(vltot (ir) + v%of_r (ir, is), 0.d0,kind=DP)
        enddo
     else
        do ir = 1, dfftp%nnr
           veff (ir, is) = CMPLX(v%of_r (ir, is), 0.d0,kind=DP)
        enddo
     endif
     CALL fwfft ('Rho', veff (:, is), dfftp)
  enddo
  !
  !     We compute here four of the five integrals needed in the phonon
  !
  fact1 = CMPLX(0.d0, - tpiba * omega,kind=DP)
  !
  do ntb = 1, ntyp
     if (upf(ntb)%tvanp ) then
        do ih = 1, nh (ntb)
           do jh = ih, nh (ntb)
              ijh = ijtoh(ih,jh,ntb) 
              !
              !    compute the augmentation function
              !
              call qvan2 (ngm, ih, jh, ntb, qmodg, qgm, ylmk0)
              !
              if (.not.lgamma) call qvan2 (ngm, ih, jh, ntb, qmod, qgmq, ylmkq)
              !
              !     NB: for this integral the moving atom and the atom of Q
              !     do not necessarily coincide
              !
              do nb = 1, nat
                 if (ityp (nb) == ntb) then
                    do ig = 1, ngm
                       aux1 (ig) = qgmq (ig) * eigts1 (mill(1,ig), nb) &
                                             * eigts2 (mill(2,ig), nb) &
                                             * eigts3 (mill(3,ig), nb)
                    enddo
                    do na = 1, nat
                       fact = eigqts (na) * CONJG(eigqts (nb) )
                       !
                       !    nb is the atom of the augmentation function
                       !
                       nta = ityp (na)
                       !  
                       IF (do_cutoff_2D) THEN
                          do ig=1, ngm
                             sk(ig)=(vlocq(ig,nta)+lr_Vlocq (ig, nta)) &
                                               * eigts1(mill(1,ig), na) &
                                               * eigts2(mill(2,ig), na) &
                                               * eigts3(mill(3,ig), na)
                          enddo
                       ELSE
                          do ig=1, ngm
                             sk(ig)=vlocq(ig,nta) * eigts1(mill(1,ig), na) &
                                                  * eigts2(mill(2,ig), na) &
                                                  * eigts3(mill(3,ig), na)
                           enddo
                       ENDIF
                       ! 
                       do ipol = 1, 3
                          do ig=1, ngm
                            aux5(ig)= sk(ig) * (g (ipol, ig) + xq (ipol) )
                          enddo
                          int2 (ih, jh, ipol, na, nb) = fact * fact1 * &
                                zdotc (ngm, aux1, 1, aux5, 1)
                       enddo
                    enddo
                    if (.not.lgamma) then
                       do ig = 1, ngm
                          aux1 (ig) = qgm (ig) * eigts1 (mill(1,ig), nb) &
                                               * eigts2 (mill(2,ig), nb) &
                                               * eigts3 (mill(3,ig), nb)
                       enddo
                    endif
                    do is = 1, nspin_mag
                       do ipol = 1, 3
                          do ig = 1, ngm
                             aux2 (ig) = veff (dfftp%nl (ig), is) * g (ipol, ig)
                          enddo
                          int1 (ih, jh, ipol, nb, is) = - fact1 * &
                               zdotc (ngm, aux1, 1, aux2, 1)
                       enddo
                    enddo
                 endif
              enddo
           enddo
        enddo
        do ih = 1, nh (ntb)
           do jh = ih + 1, nh (ntb)
              !
              !    We use the symmetry properties of the integral factor
              !
              do nb = 1, nat
                 if (ityp (nb) == ntb) then
                    do ipol = 1, 3
                       do is = 1, nspin_mag
                          int1(jh,ih,ipol,nb,is) = int1(ih,jh,ipol,nb,is)
                       enddo
                       do na = 1, nat
                          int2(jh,ih,ipol,na,nb) = int2(ih,jh,ipol,na,nb)
                       enddo
                    enddo
                 endif
              enddo
           enddo
        enddo
     endif
  enddo

  !call mp_sum(  int1, intra_bgrp_comm )
  !call mp_sum(  int2, intra_bgrp_comm )

  IF (noncolin) THEN
     !CALL set_int12_nc(0)
     int1_nc=(0.d0,0.d0)
     IF (lspinorb) int2_so=(0.d0,0.d0)
     DO np = 1, ntyp
        IF ( upf(np)%tvanp ) THEN
           DO na = 1, nat
              IF (ityp(na)==np) THEN
                 IF (upf(np)%has_so) THEN
                    CALL ph_transform_int1_so(int1, int1_nc, na, 0)
                    CALL ph_transform_int2_so(int2, int2_so, na, 0)
                 ELSE
                    CALL ph_transform_int1_nc(int1, int1_nc, na, 0)
                    IF (lspinorb) CALL ph_transform_int2_nc(int2, int2_so, na, 0)
                 END IF
              END IF
           END DO
        END IF
     END DO
  END IF


  deallocate (veff)
  if (.not.lgamma) then
     deallocate(qgmq)
     deallocate (qmod)
     deallocate (ylmkq)
  endif
  deallocate (qgm)
  deallocate (ylmk0)
  deallocate (qmodg)
  deallocate (aux5)
  deallocate (aux2)
  deallocate (aux1)
  deallocate (sk)
  !call stop_clock ('ph_dvanqq')
  return
end subroutine ph_dvanqq

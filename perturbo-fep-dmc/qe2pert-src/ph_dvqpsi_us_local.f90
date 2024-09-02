!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! adapted from QE/PHonon/PH/dvqpsi_us.f90
!----------------------------------------------------------------------
!
subroutine ph_dvqpsi_us_local &
      (ikk, xq, eigqts, npwq, igk_kq, uact, dmuxc, vlocq, dvpsi_loc, addnlcc)
  !----------------------------------------------------------------------
  !
  ! This routine calculates dV_bare/dtau * psi for one perturbation
  ! with a given q. The displacements are described by a vector u.
  ! The result is stored in dvpsi. The routine is called for each k point
  ! and for each pattern u. It computes simultaneously all the bands.
  ! It implements Eq. B29 of PRB 64, 235118 (2001). The contribution
  ! of the local pseudopotential is calculated here, that of the nonlocal
  ! pseudopotential in dvqpsi_us_only. (ph_dvqpsi_us_noloc)
  !
  ! jjzhou & ite: remove dependence on module wavefunctions, klist, eqv, nlcc_ph
  ! It computes the second term in Eq.B29 of PRB 64, 235118(2001)
  !
  USE kinds, only : DP
  USE funct,     ONLY : dft_is_gradient, dft_is_nonlocc
  USE ions_base, ONLY : nat, ityp, tau, ntyp => nsp
  USE cell_base, ONLY : tpiba, at
  USE fft_base,  ONLY : dfftp, dffts
  USE fft_interfaces, ONLY: fwfft, invfft
  USE gvect,     ONLY : eigts1, eigts2, eigts3, mill, g, ngm
  USE gvecs,     ONLY : ngms, doublegrid
  USE lsda_mod,  ONLY : nspin, lsda, isk
  USE scf,       ONLY : rho, rho_core
  USE noncollin_module, ONLY : nspin_gga, nspin_mag, npol
  use uspp_param,ONLY : upf
  USE uspp,       ONLY : nlcc_any
  USE gc_lr,      ONLY: grho, dvxc_rr,  dvxc_sr,  dvxc_ss, dvxc_s

  USE Coul_cut_2D, ONLY: do_cutoff_2D  
  USE Coul_cut_2D_ph, ONLY : cutoff_localq
  !
  use pwcom, only: npwx
  use elph_matrix, only: drc
  use electronic_data, only: ngk_tot, igk_k_tot, evc_sub, nb_sub, spin_updn
  !
  !USE nlcc_ph,    ONLY : drc
  !USE wavefunctions,  ONLY: evc
  !USE klist,      ONLY : ngk, igk_k
  !USE qpoint,     ONLY : xq, eigqts, ikqs, ikks
  !USE wvfct,     ONLY : nbnd, npwx
  !USE eqv,        ONLY : dvpsi, dmuxc, vlocq
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(in) :: ikk, npwq
  integer, intent(in) :: igk_kq(npwx)
  ! igk_kq: G in (k+q+G)

  real(dp), INTENT(in) :: xq(3)  ! q in cartesian coordinate
  !! input: the k point
  COMPLEX(DP), INTENT(in) :: uact (3*nat), eigqts(nat)
  !! input: the pattern of displacements
  ! eigqts(nat) : phase associate with q, exp(-i xq*tau)
  !
  ! the derivative of the xc potential w.r.t rho:  d v_xc / d n 
  ! (see Eq.134 in Baroni, RMP 73 515 (2001))
  REAL (DP), INTENT(in) :: dmuxc(:,:,:)  ! nrxx, nspin, nspin)
  !
  ! the local potential at q+G
  REAL (DP), INTENT(in) :: vlocq(ngm, ntyp)  ! vlocq(ngm, ntyp)
  
  LOGICAL, INTENT(in)   :: addnlcc
  ! contain the rhoc (without structure fac) for all atomic types
  ! The variable is needed for non-linear core correction
  !  
  COMPLEX (DP), INTENT(out) :: dvpsi_loc(npwx*npol, nb_sub) !(1:npwx*npol,1:nb_sub)
  ! the dv_bare/dtau * psi, local part
  !!
  !   And the local variables
  !
  INTEGER ::  na  
  !! counter on atoms
  INTEGER :: mu
  !! counter on modes
  INTEGER :: npw
  !! Number of pw
  !INTEGER :: ikk
  !! the point k
  !INTEGER :: npwq
  !! Number of q
  !INTEGER :: ikq
  !! k-q index
  !!
  INTEGER :: ig
  !! counter on G vectors
  INTEGER :: nt
  !! the type of atom
  INTEGER :: ibnd
  !! counter on bands
  INTEGER :: ir 
  !! counter on real mesh
  INTEGER :: is
  !! 
  INTEGER :: ip
  !!
  !
  complex(DP) :: gtau, gu, fact, u1, u2, u3, gu0
  complex(DP) , allocatable :: aux (:,:)
  complex(DP) , allocatable :: aux1 (:), aux2 (:)
  complex(DP) , pointer :: auxs (:)
  COMPLEX(DP), ALLOCATABLE :: drhoc(:,:)
  !
!!$omp master
!  call start_clock ('ph_dvqpsi_us_local')
!!$omp end master

  if (nlcc_any.and.addnlcc) then
     allocate (drhoc( dfftp%nnr, nspin))
     allocate (aux( dfftp%nnr, nspin))
     allocate (auxs(dffts%nnr))
  endif
  allocate (aux1(dffts%nnr))
  allocate (aux2(dffts%nnr))
  !
  !    We start by computing the contribution of the local potential.
  !    The computation of the derivative of the local potential is done in
  !    reciprocal space while the product with the wavefunction is done in
  !    real space
  !
  !ikk = ikks(ik)
  !ikq = ikqs(ik)
  !npw = ngk(ikk)
  !npwq= ngk(ikq)
  !ikk = ik
  !xq: q in cartesian coordinate
  !xq_cryt(1:3) = xq(1:3)
  !transfer to crytal coordinate
  !call cryst_to_cart(1, xq_cryt, at, -1)
  !ikq = ktable( kpt2num( xk_tot(:,ik)+xq_cryt, kdim) )
  npw = ngk_tot(ikk)
  !npwq= ngk_tot(ikq)
  
  ! 
  dvpsi_loc(:,:) = (0.d0, 0.d0)
  aux1(:) = (0.d0, 0.d0)
  do na = 1, nat
     fact = tpiba * (0.d0, -1.d0) * eigqts (na)
     mu = 3 * (na - 1)
     if (abs (uact (mu + 1) ) + abs (uact (mu + 2) ) + abs (uact (mu + &
          3) ) .gt.1.0d-12) then
        nt = ityp (na)
        u1 = uact (mu + 1)
        u2 = uact (mu + 2)
        u3 = uact (mu + 3)
        gu0 = xq (1) * u1 + xq (2) * u2 + xq (3) * u3
        do ig = 1, ngms
           gtau = eigts1 (mill(1,ig), na) * eigts2 (mill(2,ig), na) * &
                  eigts3 (mill(3,ig), na)
           gu = gu0 + g (1, ig) * u1 + g (2, ig) * u2 + g (3, ig) * u3
           aux1 (dffts%nl (ig) ) = aux1 (dffts%nl (ig) ) + vlocq (ig, nt) * gu * &
                fact * gtau
        enddo
        IF (do_cutoff_2D) then  
           call cutoff_localq( aux1, fact, u1, u2, u3, gu0, nt, na) 
        ENDIF
        !
     endif
  enddo
  !
  ! add NLCC when present
  !
  if (nlcc_any.and.addnlcc) then
     drhoc(:,:) = (0.d0, 0.d0)
     aux(:,:) = (0.0_dp, 0.0_dp)
     do na = 1,nat
        fact = tpiba*(0.d0,-1.d0)*eigqts(na)
        mu = 3*(na-1)
        if (abs(uact(mu+1))+abs(uact(mu+2))  &
                        +abs(uact(mu+3)).gt.1.0d-12) then
           nt=ityp(na)
           u1 = uact(mu+1)
           u2 = uact(mu+2)
           u3 = uact(mu+3)
           gu0 = xq(1)*u1 +xq(2)*u2+xq(3)*u3
           if (upf(nt)%nlcc) then
              do ig = 1,ngm
                 gtau = eigts1(mill(1,ig),na)*   &
                        eigts2(mill(2,ig),na)*   &
                        eigts3(mill(3,ig),na)
                 gu = gu0+g(1,ig)*u1+g(2,ig)*u2+g(3,ig)*u3
                 drhoc(dfftp%nl(ig),1)=drhoc(dfftp%nl(ig),1)+drc(ig,nt)*gu*fact*gtau
              enddo
           endif
        endif
     enddo
     CALL invfft ('Rho', drhoc(:,1), dfftp)
     if (.not.lsda) then
        ! compute \partial Vxc / \partial u
        !(see Eq.134 in Baroni, RMP 73 515 (2001))
        do ir=1,dfftp%nnr
           aux(ir,1) = drhoc(ir,1) * dmuxc(ir,1,1)
        end do
     else
        is = spin_updn ! isk(ikk), jjzhou
        do ir=1,dfftp%nnr
           drhoc(ir,1) = 0.5d0 * drhoc(ir,1)
           drhoc(ir,2) = drhoc(ir,1)
           aux(ir,1) = drhoc(ir,1) * ( dmuxc(ir,is,1) + &
                                       dmuxc(ir,is,2) )
        enddo
     endif

     rho%of_r(:,1) = rho%of_r(:,1) + rho_core(:)

     IF ( dft_is_gradient() ) CALL dgradcorr (dfftp, rho%of_r, grho, dvxc_rr, &
                    dvxc_sr, dvxc_ss, dvxc_s, xq, drhoc, nspin, nspin_gga, g, aux)       

     IF (dft_is_nonlocc()) CALL dnonloccorr(rho%of_r, drhoc, xq, aux)

     rho%of_r(:,1) = rho%of_r(:,1) - rho_core(:)

     CALL fwfft ('Rho', aux(:,1), dfftp)
! 
!  This is needed also when the smooth and the thick grids coincide to
!  cut the potential at the cut-off
!
     auxs(:) = (0.d0, 0.d0)
     do ig=1,ngms
        auxs(dffts%nl(ig)) = aux(dfftp%nl(ig),1)
     enddo
     aux1(:) = aux1(:) + auxs(:)
  endif
  !
  ! Now we compute dV_loc/dtau in real space
  !
  CALL invfft ('Rho', aux1, dffts)
  do ibnd = 1, nb_sub
     do ip=1,npol
        aux2(:) = (0.d0, 0.d0)
        if (ip==1) then
           do ig = 1, npw
              aux2 (dffts%nl (igk_k_tot (ig,ikk) ) ) = evc_sub (ig, ibnd, ikk)
           enddo
        else
           do ig = 1, npw
              aux2 (dffts%nl (igk_k_tot (ig,ikk) ) ) = evc_sub (ig+npwx, ibnd, ikk)
           enddo
        end if
        !
        !  This wavefunction is computed in real space
        !
        CALL invfft ('Wave', aux2, dffts)
        do ir = 1, dffts%nnr
           aux2 (ir) = aux2 (ir) * aux1 (ir)
        enddo
        !
        ! and finally dV_loc/dtau * psi is transformed in reciprocal space
        !
        CALL fwfft ('Wave', aux2, dffts)
        if (ip==1) then
           do ig = 1, npwq
              dvpsi_loc(ig, ibnd) = aux2 (dffts%nl ( igk_kq(ig) ) )
           enddo
        else
           do ig = 1, npwq
              dvpsi_loc(ig+npwx, ibnd) = aux2 (dffts%nl ( igk_kq(ig) ) )
           enddo
        end if
     enddo
  enddo
  !
  deallocate (aux2)
  deallocate (aux1)
  if (nlcc_any.and.addnlcc) then
     deallocate (drhoc)
     deallocate (aux)
     deallocate (auxs)
  endif
  !
  !   We add the contribution of the nonlocal potential in the US form
  !   First a term similar to the KB case.
  !   Then a term due to the change of the D coefficients.
  !
  !call dvqpsi_us_only (ikk, uact)
  !J.Z: replace with ph_dvqpsi_us_noloc

!!$omp master
!  call stop_clock ('ph_dvqpsi_us_local')
!!$omp end  master

  return
end subroutine ph_dvqpsi_us_local

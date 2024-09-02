!===============================================================================
! Copyright (C) 2016-2020 Jin-Jian Zhou, Jinsoo Park, I-Te Lu, Marco Bernardi
!
! This program is distributed under the terms of the GNU General Public License.
! See the file `LICENSE' in the root directory of this distribution, or obtain 
! a copy of the License at <https://www.gnu.org/licenses/gpl-3.0.txt>.
!
! Author: jjzhou <jjchou.comphy@gmail.com>; 
! Comment:
!   For L(S)DA+U support. 
!   relevant subroutines in Phonon/PH/: 
!
!===============================================================================

module dvhub_mod
   use kinds,      only: dp
   use parameters, only: ntypx
   use io_global,  only: stdout
   use cell_base,  only: at, bg
   use lsda_mod,   only: nspin
   use ions_base,  only: tau, nat, ityp, ntyp => nsp
   use ldaU, only: Hubbard_U, Hubbard_J0, Hubbard_lmax, is_hubbard, Hubbard_l
   use symm_base,  only: nsym, s, sr, invs, irt, d1, d2, d3
   implicit none
   public

   integer, save :: ldim

   real(dp), save :: effU(ntypx)
   ! !effU(ntypx), same as effU in module ldaU_ph
   !effective Hubbard parameter: effU = Hubbard_U - Hubbard_J0
   
   !for each atom gives the offset of beta functions, moved here from control_lr
   integer, allocatable, save :: ofsbeta(:) ! ofsbeta(nat)


   public :: setup_hubbard, obtain_dnshub
contains

!DFPT+U, following Phonon/PH/phq_setup.f90, setup effU, ofsbeta, and d1, d2, d3
subroutine setup_hubbard()
   use uspp_param,    only: nh
   implicit none
   integer :: nah, na, nt, jkb2, ih
   
   ldim = 2 * Hubbard_lmax + 1

   !Define effU, copied from phq_setup.f90 
   effU = 0.d0
   do nah = 1, nat
      nt = ityp(nah)
      ! For U only calculations
      effU(nt) = Hubbard_U(nt)
      ! When there is also Hubbard_J0/=0
      if( Hubbard_J0(nt).ne.0.d0 ) effU(nt) = Hubbard_U(nt) - Hubbard_J0(nt)
   enddo
   !
   ! Initialize d1, d2, d3 to rotate the spherical harmonics
   call d_matrix (d1, d2, d3)
   
   ! Calculate the offset of beta functions for all atoms. 
   ! call setup_offset_beta() !unfold this subroutine below
   !
   ! Calculate the offset of beta functions for each atom na.
   ! Ordering: first all betas for atoms of type 1,
   ! then all betas for atoms of type 2, and so on.
   !   
   allocate( ofsbeta(nat) )
   ofsbeta = 0
   ! 
   jkb2 = 0
   do nt = 1, ntyp
      do na = 1, nat 
         if( ityp(na).eq.nt ) then
            ofsbeta(na) = jkb2
            do ih = 1, nh(nt)
               jkb2 = jkb2 + 1
            enddo
         endif
      enddo
   enddo

end subroutine setup_hubbard


subroutine obtain_dnshub(suffix, current_irrq, dns_irrq, iq, dns_iq, lnew_irrq)
   use lattice_data, only: xq_tot, identity_op, q2minus_q, pattern, xq_irr, q2irq_symop
   implicit none
   logical, intent(in) :: lnew_irrq ! if .true., init dns_irrq (reading from file)
   character(len=*), intent(in) :: suffix
   integer, intent(in) :: current_irrq
   complex(dp), intent(inout) :: dns_irrq(ldim, ldim, nspin, nat, 3, nat)
   !
   integer, intent(in) :: iq
   complex(dp), intent(out) :: dns_iq(ldim, ldim, nspin, nat, 3, nat)
   !local
   integer :: isym
   real(dp) :: xq_cart(3), xq_irr_cart(3)
      
   ! xq_tot(:,iq): the iq-th point in crystal coordinate
   xq_cart = xq_tot(:, iq)
   !transfer xq to cartesian coordinate
   call cryst_to_cart(1, xq_cart, bg, 1)

   !if lnew_irrq is .true., init dns_irrq by reading from file
   ! dns_irq in cartesian coordinates
   if(lnew_irrq) call load_dns(trim(suffix), current_irrq, pattern(:,:,current_irrq), dns_irrq)
   
   ! xq_irr(:,irrq): the irrq-th irreducible q-point in crystal coordinate
   xq_irr_cart = xq_irr(:, current_irrq)
   ! transfer to cartesian coordinate
   call cryst_to_cart(1, xq_irr_cart, bg, 1)
   
   ! sr(:,:,isym) is the symmetry operation S that Sq = q_irr
   isym = q2irq_symop(iq)
   if( (isym .eq. identity_op) .and. (.not. q2minus_q(iq)) ) then
      !write(stdout,'(a15, i3, 2x, 3f12.8)') "irreducible q: ", current_irrq, xq_irr_cart(:)
      dns_iq(:,:,:,:,:,:) = dns_irrq(:,:,:,:,:,:)
      !
   else
      !write(stdout,'(3f12.8,2x,3f12.8,2x,i2,1x,l)') xq_cart, xq_irr_cart(:), isym, q2minus_q(iq)
      !dvscf in cartesian coordinates for xq, rotated from dvscf_irrq
      call get_dns_star(xq_cart, xq_irr_cart, isym, q2minus_q(iq), dns_irrq, dns_iq)
   endif
end subroutine


!---------------------------------------------------------------------------
! read dnsscf, follow Phonon/PH/elphon.f90 and dnsq_scf.f90
! or read dnsbare (patch dnsq_bare to output it in symmetrized pattern)
! 
! Note that in 'Phonon/PH/dnsq_bare.f90', the (unsymmetrized) dnsbare 
!  in cartesian coordinate is written to file. Symmetrizing it here is quite
!  tedious, so I patched dnsq_bare.f90 to output symmetrized dnsbare in 
!  pattern coordinate, the same way dnsscf is written out, so we can process 
!  dnsbare and dnsscf in a consistent and simpler way. 
!--------------------------------------------------------------------------
subroutine load_dns(suffix, irrq, u, dns)
   use input_param, only: prefix, phdir
   implicit none
   integer, intent(in) :: irrq
   character(len=*), intent(in) :: suffix  ! 'scf' or 'bare'
   complex(dp), intent(in) :: u(3*nat, 3*nat)
   complex(dp), intent(out) :: dns(ldim, ldim, nspin, nat, 3, nat)
   !local
   integer :: iunit, na_icart, imode, nah, nt, is, m1, m2, na, icart
   character(len=120) :: ctmp, fname
   complex(dp), allocatable :: dns_pattern(:,:,:,:,:)
   !functions:
   integer, external :: find_free_unit
   CHARACTER(LEN=6), EXTERNAL :: int_to_char

   allocate( dns_pattern(ldim, ldim, nspin, nat, 3*nat) )
   dns_pattern = (0.d0, 0.d0)
   !
   ctmp = '.dns' // trim(suffix) // '_q' // trim(int_to_char(irrq))
   fname = trim(phdir) // trim(prefix) // trim(ctmp)
   !
   iunit =find_free_unit()
   call open_dns_file(iunit, fname)
   ! 
   ! N.B.: 
   ! dnsscf includes dnsorth in case of USPP, check Phonon/PH/dnsq_scf.90
   ! dnsscf in pattern is already symmetrized, check Phonon/PH/dnsq_scf.90
   !
   ! dnsbare outputed in original dnsq_bare.f90 is not symmetrized and in cartesian cooridnate
   ! To avoid tedious symmetrization here, I patched dnsq_bare.f90 to output 
   ! symmetrized dnsbare in pattern, the same why dnsscf is outputed.
   ! so we can read dnsscf and dnsbare using the same process.
   read(iunit, *) dns_pattern
   close(iunit, status='keep')
   !

   !transfer dnsscf_pattern from pattern to cartesian coordinates
   ! adapted from PHonon/PH/elphon.f90/elphel_read_dnsscf_check()
   !
   dns = (0.d0, 0.d0)
   !
   do na = 1, nat
   do icart = 1, 3
      na_icart = (na-1)*3 + icart
      DO imode = 1, 3*nat
         DO nah = 1, nat
            nt = ityp(nah)
            IF (is_hubbard(nt)) THEN
               DO is = 1, nspin
               DO m1 = 1, 2*Hubbard_l(nt) + 1
               DO m2 = 1, 2*Hubbard_l(nt) + 1
                  !
                  dns(m1, m2, is, nah, icart, na) = dns(m1, m2, is, nah, icart, na) + &
                         dns_pattern (m1, m2, is, nah, imode) * CONJG(u(na_icart, imode))
                  !
               ENDDO;  ENDDO; ENDDO
            ENDIF
         ENDDO
      ENDDO
   ENDDO
   ENDDO

   deallocate( dns_pattern )
end subroutine load_dns


subroutine get_dns_star(xq, xq_irr, isym_inv, q2mq, dns_irr, dns_xq)
   use constants,        only: tpi
   implicit none
   ! xq     : q vector, in cartesian coordinate
   ! xq_irr : q_irr vector, in cartesian coordinate
   real(dp),intent(in) :: xq(3), xq_irr(3)
   ! index of symmetry operation that rotate q -> qirr, e.g. sr(:,:,isym_inv)
   !  is the symmetry operation S^-1 that S^-1 q = q_irr (so q = S q_irr)
   integer, intent(in) :: isym_inv
   ! if q2mq is ture, then S^-1 (-q) = q_irr, 
   ! q and q_irr is connect via S^-1 and time reversal symmetry.
   logical, intent(in) :: q2mq
   !
   !dns matrix in the cartesian coordinates for the irreducible q
   complex(dp), intent(in) :: dns_irr(ldim, ldim, nspin, nat, 3, nat)
   !dns matrix in the cartesian coordinates for the reducible q
   complex(dp), intent(out) :: dns_xq(ldim, ldim, nspin, nat, 3, nat)
   !local
   integer :: na, nb, ipol, j, nt, nar, m1, m2, m0, m00
   ! auxiliary xq\cdot\tau and \xq_s\cdot\tau
   real(dp) :: xq_tau, xqirr_tau, xq_tmp(3), xq_rot(3)
   character(len=6), external :: int_to_char
   !
   complex(dp) :: phase_xqirr, phase_xq
   complex(dp), allocatable :: dns_tmp(:,:,:,:,:,:)
   !
   ! pointer to maxtrix d1, d2, or d3 for rotation of spherical harmonics 
   ! (dl for l=1,2..), check PW/symm_base.f90 or PHonon/PH/sym_dns.f90
   real(dp), target :: d0(1,1)
   real(dp), pointer :: dd(:,:)
   
   !
   if(isym_inv < 1 .or. isym_inv > 48) call errore('get_dns_star', &
      'illegal symmetry operation: '//int_to_char(isym_inv), 1)

   xq_rot = merge(-xq, xq, q2mq)
   !check symmetry operation: S^-1 q = q_irr;
   ! both xq_rot, xq_irr and sr are in cartesian coordinate.
   xq_tmp = matmul(sr(:,:,isym_inv), xq_rot) - xq_irr
   !
   ! Note: we choose the reducible q list so that all the q and their 
   ! corresponding q_irr are connected directly via S^-1.
   ! e.g. S^-1 q = q_irr + G is not used. check qe2pert/lattice_data.f90/find_symmetry_index.
   !
   if(sqrt(dot_product(xq_tmp, xq_tmp)) > 1.0E-5_dp) then
      write(*, '(9(f8.4,1x))') matmul(sr(:,:,isym_inv), xq_rot), xq_irr, xq_tmp
      call errore('dvscf','S^-1 xq is not equal to xq_irr',1)
   endif
   !
   ! Transform to crystalline coordinates (necessary in order to apply s)
   allocate( dns_tmp(ldim, ldim, nspin, nat, 3, nat) )
   dns_tmp = cmplx(0.0_dp, 0.0_dp, kind=dp)
   !
   do na = 1, nat
      do ipol = 1, 3
         do j = 1, 3
            dns_tmp(:,:,:,:,ipol,na) = dns_tmp(:,:,:,:,ipol,na) + &
                                    dns_irr(:,:,:,:,j,na) * at(j,ipol)
         enddo
      enddo
   enddo

   ! take away the phase due to the q-point
   do na = 1, nat
      do nb = 1, nat
         nt = ityp(nb)

         if( is_hubbard(nt) ) then
            xqirr_tau = tpi * dot_product( xq_irr(:), tau(:,na)-tau(:,nb) )
            phase_xqirr= cmplx( cos(xqirr_tau), sin(xqirr_tau), kind=dp )
            !
            dns_tmp(:,:,:,nb,:,na) = phase_xqirr * dns_tmp(:,:,:,nb,:,na)
         endif
      enddo
   enddo
   !
   d0(1,1) = 1.0_dp
   !
   !similar to ruotaijk in dvscf_mod, e.g. dv(S^-1 r) in the case of dvscf
   ! check Phonon/PH/sym_dns.f90, and also compare it with symdvscf.f90 
   ! 
   dns_xq = cmplx(0.0_dp, 0.0_dp, kind=dp)
   !
   !N.B.: for non-collinear case (with SOC), spin and spatial spaces are coupled, 
   ! so one need to do the rotation in both spin and spatial space at the same time. 
   ! [see PW/new_ns.f90/new_ns_nc and PW/plus_u_full.f90/comp_dspinldau]
   ! 
   ! Here We only support DFPT+U with collinear (spinless or LSDA), in which 
   ! spin space spatial space are decoupled, so we only need to do spatial rotation.
   ! [comparing subroutine 'new_ns' to 'new_ns_nc' in PW/new_ns.f90]
   !
   do na = 1, nat
      nt = ityp(na)
      if( is_hubbard(nt) ) then
         !set rotation matrix
         select case( Hubbard_l(nt) )
         case (0)
            dd => d0
         case (1)
            dd => d1(:,:,isym_inv)
         case (2)
            dd => d2(:,:,isym_inv)
         case (3)
            dd => d3(:,:,isym_inv)
         case default
            call errore('get_dns_star','angular momentum not implemented',1)
         end select
         !
         !perform the rotation
         do m1 = 1, 2 * Hubbard_l(nt) + 1
         do m2 = 1, 2 * Hubbard_l(nt) + 1
            !irt is the list of rotated of each atom: nar = [S^-1] na
            nar = irt(isym_inv, na)
            do m0  = 1, 2 * Hubbard_l(nt) + 1
            do m00 = 1, 2 * Hubbard_l(nt) + 1
               dns_xq(m1, m2,:, na,:,:) = dns_xq(m1, m2,:, na,:,:) + &
                  dd(m0, m1) * dd(m00, m2) * dns_tmp(m0, m00,:, nar,:,:)
            enddo
            enddo
         enddo
         enddo
      endif
   enddo

   ! Compute dns at xq
   dns_tmp = cmplx(0.0_dp, 0.0_dp, kind=dp)
   do na = 1, nat
      nar = irt(isym_inv, na)
      
      do ipol = 1, 3
         !N.B. s(a, b, isym_inv) ==> [S^-1]_ba in C13 of PRB97,235146
         do j = 1, 3
            dns_tmp(:,:,:,:,ipol, na) = dns_tmp(:,:,:,:,ipol, na) + &
               s(ipol, j, isym_inv) * dns_xq(:,:,:,:,j, nar)
         enddo
      enddo
   enddo
   !
   !add back the phase factor for the new q-point
   do na = 1, nat
      do nb = 1, nat
         nt = ityp(nb)

         if( is_hubbard(nt) ) then
            xq_tau = - tpi * dot_product( xq_rot(:), tau(:,na)-tau(:,nb) )
            phase_xq= cmplx( cos(xq_tau), sin(xq_tau), kind=dp )
            !
            dns_tmp(:,:,:,nb,:,na) = phase_xq * dns_tmp(:,:,:,nb,:,na)
         endif
      enddo
   enddo
   !
   ! back to cartesian coordinates
   dns_xq = cmplx(0.0_dp, 0.0_dp, kind=dp)
   !
   do na = 1, nat
      do ipol = 1, 3
      do j = 1, 3
         dns_xq(:,:,:,:,ipol, na) = dns_xq(:,:,:,:,ipol, na) + &
            dns_tmp(:,:,:,:,j, na) * bg(ipol, j)
      enddo
      enddo
   enddo

   !the occupation matrix ns is real for collinear case (spinless or LSDA).
   ! [see PW/new_ns.f90/new_ns and PW/v_of_rho.f90/v_hubbard]
   ! so similar to dvscf in dvscf_mod, we can use the following relation.
   ! [see Eq.(42) in PRB 98,085127 (2018)]
   !
   !dns_xq(q) = conjg(dns_xq(-q))
   if(q2mq) dns_xq = conjg( dns_xq )
   !
   !N.B.: for non-collinear case, ns is a complex matrix.
   ! [see PW/new_ns.f90/new_ns_nc or PW/v_of_rho.f90/v_hubbard_nc]
   ! hence the above can not be applied. 
   ! (DFPT+U with non-colliner is not supported yet)

   deallocate(dns_tmp)
end subroutine get_dns_star


subroutine open_dns_file(iunit, fname)
   implicit none
   integer, intent(in) :: iunit
   character(len=*), intent(in) :: fname
   !logical, intent(out) :: exst
   !
   integer :: ios
   logical :: opnd, exst
   !
   if(iunit < 0) call errore('open_dns_file', 'wrong unit', 1)
   !
   inquire(unit=iunit, opened=opnd)
   if(opnd) call errore('open_dns_file', "can't open a connected unit", abs(iunit))
   !
   inquire(file=trim(fname), exist=exst)
   if(.not.exst) call errore('open_dns_file',trim(fname)//' does not exist !',1)
   !
   !   Open the file
   !
   open(unit=iunit, file=trim(fname), form='formatted', status='unknown', &
         action='read', iostat=ios)
   if(ios /= 0) call errore('open_dns_file','error opening '//trim(fname), iunit)
   return
end subroutine open_dns_file

end module dvhub_mod

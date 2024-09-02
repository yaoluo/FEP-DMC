!===============================================================================
! Copyright (C) 2016-2020 Jin-Jian Zhou, Jinsoo Park, I-Te Lu, Marco Bernardi
!
! This program is distributed under the terms of the GNU General Public License.
! See the file `LICENSE' in the root directory of this distribution, or obtain 
! a copy of the License at <https://www.gnu.org/licenses/gpl-3.0.txt>.
!
! Author: jpark <jinsoop412@gmail.com>
! Comment:
!  Based on  PH/q2r.f90 
!            PH/ph_restart.f90
!            PH/q2qstar.f90
!  Reads the relevant phonon datas from a previous phonon run. 
!
! Maintenance:
!===============================================================================

module lattice_data
   use iotk_module
   use kinds,       only : dp
   use constants,   only : amu_ry
   use ions_base,   only : nat, ntyp=>nsp, ityp, amass, tau
   use input_param, only : prefix, phdir, asr, system_2d
   use io_global,   only : ionode, ionode_id, stdout
   use mp,          only : mp_bcast, mp_barrier
   use mp_pools,    only : inter_pool_comm
   !
   implicit none
   public
   !
   integer, save :: numq
   ! total number of full q-points
   integer, save :: nq_irr
   ! total number of irreducible q-points
   integer, save :: qdim(3)
   ! the dimension of the coarse q-grid
   !
   real(dp), allocatable, save :: xq_irr(:,:)
   ! crystal coordinates of irreducible q-points
   ! xq_irr(3, nq_irr)
   real(dp), allocatable, save :: xq_tot(:,:)
   ! crystal ccoordinates of full q-points
   ! xq_tot(3, numq)
   integer, allocatable, save :: iq2irrq(:) 
   ! location of the correspoding irr. q for each q in xq_tot(:,:).
   ! xq_tot(:,i) belongs to the star of the irreducible q-ponit xq_irr(:, iq2irrq(i))
   ! iq2irrq(numq)
   integer, allocatable, save :: q2irq_symop(:)  
   logical, allocatable, save :: q2minus_q(:)
   ! q2irq_symop(iq): symmetry operation i that connect iq to its irreduicble q. 
   ! if q2minus_q(iq) is true, then q2irq_symop(iq) + (q->-q) connect iq to irreduicble q.
   ! only store its index in sr(:,:,:) matrix in module PW/symm_base
   ! q2irq_symop(numq), q2minus_q(numq)
   integer, save :: identity_op=1 
   ! identity operation index, if q2irq_symop(i) .eq. identity_op. this is always 1.
   ! then xq_tot(:,i) is an irreducible q-point
   !
   complex(dp), allocatable, save :: dynmat(:,:,:,:,:) 
   ! dynamical matrices produced by the phonon code
   ! dynmat(numq,3,3,nat,nat) 
   complex(dp), allocatable, save :: pattern(:,:,:) 
   ! patterns of each irreducible q-points produced by the phonon code
   ! pattern(3*nat,3*nat,nq_irr)
   !
   integer, allocatable, save :: num_q_star(:)  
   ! Number of q-points in the star 
   ! num_q_star(nq_irr) 
   integer, allocatable, save :: irr_iq(:)
   ! number of irreducible representation per q point
   ! irr_iq(nq_irr)
   integer, allocatable, save :: npert_irr_iq(:,:) 
   ! the number of perturbations(modes) for each q and irreducible representation
   ! npert_irr_iq(3*nat,nq_irr) 
   integer, allocatable, save :: nsymq_iq(:)
   ! the order of the small group of q for each q (the number of symmetry of the small group)
   ! nsymq_iq(nq_irr)
   !
   logical, save :: ph_lpolar
   ! .true. if the system is a polar insulator/semiconductor
   real(dp), save :: ph_epsil(3,3)
   ! dielectric tensor
   real(dp), allocatable, target, save :: ph_zstar(:,:,:)
   ! born effective charge, ph_zstar(3,3,nat)
   real(dp), parameter :: ph_alpha = 1.0E0_dp !Ewald parameter, set to the value used in rgd_blk
   !
   integer, allocatable, save :: count_q2m1(:)
   integer, allocatable, save :: count_q2m2(:)
   integer, allocatable, save :: count_q2m3(:)
   ! map of count_q to m(3)
   ! used in dynmat->ifc transformation

   public :: init_lattice_data 
contains

subroutine init_lattice_data()
   use cell_base, only : at
   implicit none
   !
   character (len=256) :: phonon_path
   !
   phonon_path=trim(phdir)
   !
   ! Read dynamical matrices and set the values of 
   ! numq, nq_irr, qdim, xq_irr, xq_tot, iq2irrq, dynmat
   call read_dynmat(phonon_path)
   !
   ! Set symmetry matrices and set the values of
   ! q2irq_symop 
   call set_symmetry()
   !
   ! Read patterns and set the values of 
   ! pattern, numq_q_star, irr_iq, npert_irr_iq, nsymq_iq
   call read_patterns(phonon_path)
   !
   !convert both xq_tot and xq_irr to crystal coordinates
   call cryst_to_cart (numq,   xq_tot(:,:), at, -1)
   call cryst_to_cart (nq_irr, xq_irr(:,:), at, -1)
   !
   write(stdout,*) 'end lattice'
   !
   call mp_barrier(inter_pool_comm)
end subroutine init_lattice_data

SUBROUTINE read_dynmat(phonon_path)
   !
   ! adapted from q2r.f90 v6.4.1
   !     reads force constant matrices C(q) produced by the phonon code
   !     for a grid of q-points, calculates the corresponding set of
   !     interatomic force constants (IFC), C(R)
   !
   use cell_base, only: omega, at, bg
   USE io_dyn_mat, ONLY : read_dyn_mat_param, read_dyn_mat_header, &
      read_dyn_mat, read_dyn_mat_tail !, write_dyn_mat_header, write_ifc
   USE rigid, ONLY: rgd_blk
   !
   IMPLICIT NONE
   !
   INTEGER,  PARAMETER :: ntypx = 10
   REAL(DP), PARAMETER :: eps=1.D-5, eps12=1.d-12
   INTEGER  :: nr1, nr2, nr3, nr(3)
   !     dimensions of the FFT grid formed by the q-point grid
   !
   CHARACTER(len=20)  :: crystal
   CHARACTER(len=3)   :: atm(ntypx)
   CHARACTER(LEN=6), EXTERNAL :: int_to_char
   CHARACTER(len=256) :: fildyn, filin, filj, filf, phonon_path
   !
   LOGICAL :: lq, lrigid, lrigid1, xmldyn
   INTEGER :: m1, m2, m3, m(3), l1, l2, l3, i, j, j1, j2, na1, na2, ipol, nn
   INTEGER :: nq, ntyp_, nat_, iq, icar, nfile, ifile, nqs, nq_log
   INTEGER :: na, nt, count_q
   !
   INTEGER :: gid, ibrav_, ierr, nspin_mag_, ios, iunit
   !
   INTEGER, ALLOCATABLE ::  nc(:,:,:), ityp_(:)
   REAL(DP), ALLOCATABLE :: m_loc_(:,:), tau_(:,:)
   COMPLEX(DP), ALLOCATABLE :: phiq(:,:,:,:,:) !, phid(:,:,:,:,:)
   !
   REAL(DP) :: celldm_(6), at_(3,3), bg_(3,3)  
   REAL(DP) :: q(3,48), omega_, xq, amass_(ntypx) !, resi
   ! q(3,48) : list of q vectors in the star (cartesian coordinates)
   !
   integer, external :: find_free_unit
   !
   logical :: loto_2d
   !
   fildyn = trim(phonon_path) // trim(prefix) // '.dyn'
   !
   ! only .xml extension for dynmat is supported
   !
   CALL mp_bcast(fildyn, ionode_id, inter_pool_comm)

   IF (ionode) THEN
      iunit = find_free_unit()
      OPEN (unit=iunit, file=TRIM(fildyn)//'0', status='old', form='formatted', iostat=ierr)
      IF ( ierr /= 0 ) THEN
         ! *.dyn0 file is required.
         CALL errore ('read_dynmat','dyn0 file missing',1)
      ELSE
         WRITE (stdout,'(/,4x," reading grid info from file ",a)') TRIM(fildyn)//'0'
         READ (iunit, *) nr1, nr2, nr3
         READ (iunit, *) nfile
         CLOSE (unit=iunit, status='keep')
      END IF
   ENDIF
   !
   CALL mp_bcast(nr1, ionode_id, inter_pool_comm)
   CALL mp_bcast(nr2, ionode_id, inter_pool_comm)
   CALL mp_bcast(nr3, ionode_id, inter_pool_comm)
   CALL mp_bcast(nfile, ionode_id, inter_pool_comm)
   !
   IF (nr1 < 1 .OR. nr1 > 1024) CALL errore ('read_dynmat',' nr1 wrong or missing',1)
   IF (nr2 < 1 .OR. nr2 > 1024) CALL errore ('read_dynmat',' nr2 wrong or missing',1)
   IF (nr3 < 1 .OR. nr2 > 1024) CALL errore ('read_dynmat',' nr3 wrong or missing',1)
   IF (nfile < 1 .OR. nfile > 1024) &
   CALL errore ('read_dynmat','too few or too many file',MAX(1,nfile))
   !
   ! set and allocate the global variables 
   numq = nr1*nr2*nr3
   nq_irr = nfile
   qdim(1) = nr1
   qdim(2) = nr2
   qdim(3) = nr3
   ph_lpolar  = .false.

   allocate( xq_irr(3,nq_irr) )
   allocate( xq_tot(3,numq) )
   allocate( iq2irrq(numq) )
   allocate( q2irq_symop(numq) )
   allocate( q2minus_q(numq) )
   allocate( nsymq_iq(nq_irr) )
   allocate( ph_zstar(3,3,nat) )
   !
   allocate( dynmat(numq,3,3,nat,nat) )
   allocate( pattern(3*nat,3*nat,nq_irr) )
   allocate( irr_iq(nq_irr) )
   allocate( num_q_star(nq_irr) )
   allocate( npert_irr_iq(3*nat,nq_irr) )
   !
   allocate( count_q2m1(numq) )
   allocate( count_q2m2(numq) )
   allocate( count_q2m3(numq) )
   !
   ! copy nrX -> nr(X)
   !
   nr(1) = nr1
   nr(2) = nr2
   nr(3) = nr3
   !
   ! Dynamical matrix (analytical part)
   !
   ntyp_ = ntypx ! avoids spurious out-of-bound errors
   !
   ! work array
   allocate( m_loc_(3,nat) )
   allocate( tau_(3,nat) )
   allocate( ityp_(nat) )
   allocate( nc(nr1,nr2,nr3) )
   nc = 0
   !
   ! Force constants in reciprocal space read from file
   !
   count_q=0
   !
   DO ifile=1,nfile
      filin = TRIM(fildyn) // TRIM( int_to_char( ifile ) )
      WRITE (stdout,*) ' reading force constants from file ',TRIM(filin) // '.xml'
      !
      CALL read_dyn_mat_param(filin, ntyp_, nat_)
      if(ntyp_ .ne. ntyp) &
         call errore('read_dynmat', 'mismatch ntyp in '//trim(filin), 1)
      if(nat_ .ne. nat) &
         call errore('read_dynmat', 'mismatch nat in '//trim(filin), 1)

      IF (ifile==1) THEN
         CALL read_dyn_mat_header(ntyp, nat, ibrav_, nspin_mag_, &
            celldm_, at_, bg_, omega_, atm, amass_, tau_, ityp_, &
            m_loc_, nqs, lrigid, ph_epsil, ph_zstar )
      ELSE
         CALL read_dyn_mat_header(ntyp, nat, ibrav_, nspin_mag_, &
            celldm_, at_, bg_, omega_, atm, amass_, tau_, ityp_, m_loc_, nqs)
      ENDIF

      ! sanity check
      if( any(abs(ityp(1:nat) - ityp_(1:nat)) > 0) ) &
         call errore('read_dynmat', 'mismatch ityp in '//trim(filin), 1)
      if( any(abs(tau(:,1:nat) - tau_(:,1:nat)) > eps) ) &
         call errore('read_dynmat', 'mismatch tau in '//trim(filin), 1)
      if( any(abs(amass(1:ntyp) - amass_(1:ntyp)) > eps) ) &
         call errore('read_dynmat', 'mismatch amass in '//trim(filin), 1)
      if( abs(omega - omega_) > eps )  &
         call errore('read_dynmat', 'mismatch omega in '//trim(filin), 1)
      
      allocate( phiq(3,3,nat,nat,nqs) )
      DO iq=1,nqs
         CALL read_dyn_mat(nat, iq, q(:,iq), phiq(:,:,:,:,iq))
      ENDDO
      CALL read_dyn_mat_tail(nat)

      IF (ifile == 1) THEN
         ! it must be allocated here because nat is read from file
         !ALLOCATE (phid(nr1*nr2*nr3,3,3,nat,nat) )
         !
         lrigid1=lrigid

         CALL latgen(ibrav_,celldm_,at_(1,1),at_(1,2),at_(1,3),omega_)
         at_ = at_ / celldm_(1)  !  bring at_ in units of alat

         CALL volume(celldm_(1),at_(1,1),at_(1,2),at_(1,3),omega_)
         CALL recips(at_(1,1),at_(1,2),at_(1,3),bg_(1,1),bg_(1,2),bg_(1,3))

         !sanity check
         if( any(abs(at - at_) > eps) ) &
            call errore('read_dynmat', 'mismatch at in '//trim(filin), 1)
         if( any(abs(bg - bg_) > eps) ) &
            call errore('read_dynmat', 'mismatch bg in '//trim(filin), 1)
         
         if(lrigid .AND. (asr.NE.'no')) then
            CALL set_zasr ( asr, nr1,nr2,nr3, nat, ibrav_, tau, ph_zstar)
         endif
      END IF

      IF (lrigid.AND..NOT.lrigid1) CALL errore('read_dynmat', &
       & 'file with dyn.mat. at q=0 should be first of the list',ifile)
      !
      WRITE (stdout,*) ' nqs= ',nqs
      !
      ! set the number of q in the star
      num_q_star(ifile)=nqs
      !
      DO nq = 1,nqs
         WRITE(stdout,'(a,3f12.8)') ' q= ',(q(i,nq),i=1,3)
         count_q = count_q+1

         ! set global variable : xq_tot to xq_irr map
         iq2irrq(count_q)=ifile

         lq = .TRUE.
         DO ipol=1,3
            xq = 0.0d0
            DO icar=1,3
               xq = xq + at(icar,ipol) * q(icar,nq) * nr(ipol)
            END DO

            lq = lq .AND. (ABS(NINT(xq) - xq) .LT. eps)
            iq = NINT(xq)
            !
            m(ipol)= MOD(iq,nr(ipol)) + 1
            IF (m(ipol) .LT. 1) m(ipol) = m(ipol) + nr(ipol)
         END DO

         ! set global variable 
         ! q is in cartesian coordinates
         xq_tot(:,count_q)=q(:,nq)
         if(nq.eq.1) then 
            xq_irr(:,ifile)=xq_tot(:,count_q)
         endif

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !write(stdout,'(3f12.8)') xq_tot(:,count_q)
         !if(nq.eq.1) then 
         !   write(stdout,'(3f12.8)') xq_irr(:,ifile)
         !endif

         IF (.NOT.lq) CALL errore('read_dynmat','q not allowed',1)

         IF(nc(m(1),m(2),m(3)).EQ.0) THEN
            nc(m(1),m(2),m(3))=1
         
            ph_lpolar = merge(.true., .false., sum(abs(ph_zstar)) > 1.0E-1_dp)
            loto_2d = merge(.true., .false., (sum(abs(ph_zstar)) > 1.0E-1_dp).and.system_2d)

            IF (lrigid .and. ph_lpolar) THEN
               CALL rgd_blk (nr1,nr2,nr3,nat,phiq(1,1,1,1,nq),q(1,nq), &
                  tau,ph_epsil,ph_zstar,bg_,omega,celldm_(1), loto_2d,-1.d0) ! 2D added celldm_ and flag
            END IF
            ! set the dynmat values
            CALL trasl_perturbo ( dynmat, phiq, nq, nr1,nr2,nr3, nat, count_q)
            !set up  count_q2m1, count_q2m2, count_q2m3
            count_q2m1(count_q)=m(1)
            count_q2m2(count_q)=m(2)
            count_q2m3(count_q)=m(3)
         ELSE
            WRITE (stdout,'(3i4)') (m(i),i=1,3)
            CALL errore('read_dynmat',' nc already filled: wrong q grid or wrong nr',1)
         END IF
      END DO 
      DEALLOCATE(phiq)
   END DO 
   !
   ! Check grid dimension
   !
   nq_log = SUM (nc)
   IF ( (nq_log == nr1*nr2*nr3).and.(nq_log == count_q) ) THEN
      WRITE (stdout,'(/5x,a,i4)') ' q-space grid ok, #points = ',nq_log
   ELSE
      CALL errore('read_dynmat',' missing q-point(s)!',1)
   END IF
   !
   ! convert unit of ion mass to rydberg atomic unit
   amass = amass * amu_ry
   
   deallocate(nc)
   deallocate(tau_)
   deallocate(ityp_)
   deallocate(m_loc_)
   !
   !
END SUBROUTINE read_dynmat

SUBROUTINE trasl_perturbo( dynmat, phiq, nq, nr1, nr2, nr3, nat, count_q)
   IMPLICIT NONE
   INTEGER, intent(in) ::  nr1, nr2, nr3, nat, nq, count_q
   COMPLEX(DP), intent(in) :: phiq(3,3,nat,nat,48)
   COMPLEX(DP), intent(out) :: dynmat(nr1*nr2*nr3,3,3,nat,nat)
   !
   INTEGER :: j1,j2,  na1, na2
   !
   DO j1=1,3
      DO j2=1,3
         DO na1=1,nat
            DO na2=1,nat
               dynmat(count_q,j1,j2,na1,na2) = &
               0.5d0 * (      phiq(j1,j2,na1,na2,nq) +  &
                CONJG(phiq(j2,j1,na2,na1,nq)))
            END DO
         END DO
      END DO
   END DO
   !
   RETURN
END SUBROUTINE trasl_perturbo

subroutine set_zasr ( zasr, nr1,nr2,nr3, nat, ibrav, tau, zeu)
   !
   ! Impose ASR - refined version by Nicolas Mounet
   !
   USE io_global, ONLY : stdout
   implicit none
   character(len=10) :: zasr
   integer ibrav,nr1,nr2,nr3,nr,m,p,k,l,q,r
   integer n,i,j,n1,n2,n3,na,nb,nat,axis,i1,j1,na1
   !
   real(DP) sum, zeu(3,3,nat)
   real(DP) tau(3,nat), zeu_new(3,3,nat)
   !
   real(DP) zeu_u(6*3,3,3,nat)
   ! These are the "vectors" associated with the sum rules on effective charges
   !
   integer zeu_less(6*3),nzeu_less,izeu_less
   ! indices of vectors zeu_u that are not independent to the preceding ones,
   ! nzeu_less = number of such vectors, izeu_less = temporary parameter
   !
   real(DP) zeu_w(3,3,nat), zeu_x(3,3,nat),scal,norm2
   ! temporary vectors and parameters

   ! Initialization.
   ! n is the number of sum rules to be considered (if zasr.ne.'simple')
   ! and 'axis' is the rotation axis in the case of a 1D system
   ! (i.e. the rotation axis is (Ox) if axis='1', (Oy) if axis='2'
   ! and (Oz) if axis='3')
   !
   if((zasr.ne.'simple').and.(zasr.ne.'crystal').and.(zasr.ne.'one-dim') &
    .and.(zasr.ne.'zero-dim')) then
      call errore('read_dynmat','invalid Acoustic Sum Rulei for Z*:' // zasr, 1)
   endif
   if(zasr.eq.'crystal') n=3
   if(zasr.eq.'one-dim') then
      ! the direction of periodicity is the rotation axis
      ! It will work only if the crystal axis considered is one of
      ! the cartesian axis (typically, ibrav=1, 6 or 8, or 4 along the
      ! z-direction)
      if (nr1*nr2*nr3.eq.1) axis=3
      if ((nr1.ne.1).and.(nr2*nr3.eq.1)) axis=1
      if ((nr2.ne.1).and.(nr1*nr3.eq.1)) axis=2
      if ((nr3.ne.1).and.(nr1*nr2.eq.1)) axis=3
      if (((nr1.ne.1).and.(nr2.ne.1)).or.((nr2.ne.1).and. &
        (nr3.ne.1)).or.((nr1.ne.1).and.(nr3.ne.1))) then
         call errore('read_dynmat','too many directions of &
           &   periodicity in 1D system',axis)
      endif
      if ((ibrav.ne.1).and.(ibrav.ne.6).and.(ibrav.ne.8).and. &
        ((ibrav.ne.4).or.(axis.ne.3)) ) then
         write(stdout,*) 'zasr: rotational axis may be wrong'
      endif
      write(stdout,'("zasr rotation axis in 1D system= ",I4)') axis
      n=4
   endif
   if(zasr.eq.'zero-dim') n=6

   ! Acoustic Sum Rule on effective charges
   !
   if(zasr.eq.'simple') then
      do i=1,3
         do j=1,3
            sum=0.0d0
            do na=1,nat
               sum = sum + zeu(i,j,na)
            end do
            do na=1,nat
               zeu(i,j,na) = zeu(i,j,na) - sum/nat
            end do
         end do
      end do
   else
      ! generating the vectors of the orthogonal of the subspace to project
      ! the effective charges matrix on
      !
      zeu_u(:,:,:,:)=0.0d0
      do i=1,3
         do j=1,3
            do na=1,nat
               zeu_new(i,j,na)=zeu(i,j,na)
            enddo
         enddo
      enddo
      !
      p=0
      do i=1,3
         do j=1,3
            ! These are the 3*3 vectors associated with the
            ! translational acoustic sum rules
            p=p+1
            zeu_u(p,i,j,:)=1.0d0
            !
         enddo
      enddo
      !
      if (n.eq.4) then
         do i=1,3
            ! These are the 3 vectors associated with the
            ! single rotational sum rule (1D system)
            p=p+1
            do na=1,nat
               zeu_u(p,i,MOD(axis,3)+1,na)=-tau(MOD(axis+1,3)+1,na)
               zeu_u(p,i,MOD(axis+1,3)+1,na)=tau(MOD(axis,3)+1,na)
            enddo
            !
         enddo
      endif
      !
      if (n.eq.6) then
         do i=1,3
            do j=1,3
               ! These are the 3*3 vectors associated with the
               ! three rotational sum rules (0D system - typ. molecule)
               p=p+1
               do na=1,nat
                  zeu_u(p,i,MOD(j,3)+1,na)=-tau(MOD(j+1,3)+1,na)
                  zeu_u(p,i,MOD(j+1,3)+1,na)=tau(MOD(j,3)+1,na)
               enddo
               !
            enddo
         enddo
      endif
      !
      ! Gram-Schmidt orthonormalization of the set of vectors created.
      !
      nzeu_less=0
      do k=1,p
         zeu_w(:,:,:)=zeu_u(k,:,:,:)
         zeu_x(:,:,:)=zeu_u(k,:,:,:)
         do q=1,k-1
            r=1
            do izeu_less=1,nzeu_less
               if (zeu_less(izeu_less).eq.q) r=0
            enddo
            if (r.ne.0) then
               call sp_zeu(zeu_x,zeu_u(q,:,:,:),nat,scal)
               zeu_w(:,:,:) = zeu_w(:,:,:) - scal* zeu_u(q,:,:,:)
            endif
         enddo
         call sp_zeu(zeu_w,zeu_w,nat,norm2)
         if (norm2.gt.1.0d-16) then
            zeu_u(k,:,:,:) = zeu_w(:,:,:) / DSQRT(norm2)
         else
            nzeu_less=nzeu_less+1
            zeu_less(nzeu_less)=k
         endif
      enddo
      !
      ! Projection of the effective charge "vector" on the orthogonal of the
      ! subspace of the vectors verifying the sum rules
      !
      zeu_w(:,:,:)=0.0d0
      do k=1,p
         r=1
         do izeu_less=1,nzeu_less
            if (zeu_less(izeu_less).eq.k) r=0
         enddo
         if (r.ne.0) then
            zeu_x(:,:,:)=zeu_u(k,:,:,:)
            call sp_zeu(zeu_x,zeu_new,nat,scal)
            zeu_w(:,:,:) = zeu_w(:,:,:) + scal*zeu_u(k,:,:,:)
         endif
      enddo
      !
   ! Final substraction of the former projection to the initial zeu, to get
      ! the new "projected" zeu
      !
      zeu_new(:,:,:)=zeu_new(:,:,:) - zeu_w(:,:,:)
      call sp_zeu(zeu_w,zeu_w,nat,norm2)
      write(stdout,'("Norm of the difference between old and new effective ", &
       &  "charges: " , F25.20)') SQRT(norm2)
      !
      ! Check projection
      !
      !write(6,'("Check projection of zeu")')
      !do k=1,p
      !  zeu_x(:,:,:)=zeu_u(k,:,:,:)
      !  call sp_zeu(zeu_x,zeu_new,nat,scal)
      !  if (DABS(scal).gt.1d-10) write(6,'("k= ",I8," zeu_new|zeu_u(k)= ",F15.10)') k,scal
      !enddo
      !
      do i=1,3
         do j=1,3
            do na=1,nat
               zeu(i,j,na)=zeu_new(i,j,na)
            enddo
         enddo
      enddo
   endif
   !
   !
   return
end subroutine set_zasr

subroutine set_symmetry()
   !  
   !  This subroutine sets the symmetry matrices s and initializes the modes
   !  of all irreducible representations for all q points. 
   !  adapted from q2qstar.f90 
   !  
   use ions_base,     only : tau, nat
   use io_global,     only : stdout
   use symm_base,     only : nsym, sr, set_sym_bl, find_sym, nrot, invsym
   use noncollin_module, only : m_loc
   !
   implicit none
   !
   integer ::   iq
   ! counters on full q-index
   !
   ! setup bravais lattice symmetry
   !
   call set_sym_bl ( )
   write(stdout, '(5x,a,i3)') "Symmetries of bravais lattice: ", nrot
   !
   ! setup crystal symmetry
   if(.not.allocated(m_loc))  then
      allocate(m_loc(3,nat))
      m_loc = 0._dp
   endif
   call find_sym ( nat, tau, ityp, .false., m_loc )
   write(stdout, '(5x,a,i3)') "Symmetries of crystal:         ", nsym
   !
   !find the symmetry linking xq_irr and xq_tot
   do iq = 1, numq
      call find_symmetry_index(xq_tot(:,iq), xq_irr(:,iq2irrq(iq)), &
         sr, nsym, invsym, q2irq_symop(iq), q2minus_q(iq))
   enddo
   !test
   !write(stdout,*) 'q2irq_symop test'
   !do iq = 1, numq
   !   write(stdout, *) 'total index, irr index, q2irq_symop:', iq, iq2irrq(iq), q2irq_symop(iq)
   !enddo
end subroutine set_symmetry

subroutine find_symmetry_index( xq, xq_irr, sr, nsym, invsym, isym, q2mq)
   !
   ! find the symmetry index q2irq_symop
   ! which links xq_irr and xq_tot (in cartesian coord)
   !
   use cell_base, only : at
   use io_global,     only : stdout, ionode
   !
   implicit none
   !
   real(dp), intent(in):: xq(3), xq_irr(3), sr(3,3,48)
   integer, intent(in) :: nsym
   logical, intent(in) :: invsym  ! .true. if the system has inversion symmetry
   logical ,intent(out) :: q2mq
   integer, intent(out) :: isym
   !
   real(dp) :: sxq(3)
   integer :: i
   real(dp), parameter :: eps=1.d-5
   !
   !
   do isym = 1, nsym
      do i = 1, 3
         sxq(i) = dble(sr(i,1,isym))*xq(1) + dble(sr(i,2,isym))*xq(2) + dble(sr(i,3,isym))*xq(3)
      enddo
      !
      ! if sxq .eq. xq_irr then return index to q2irq_symop
      !
      if (norm2(sxq-xq_irr) .lt. eps) then
         q2mq = .false.
         return
      !if inversion symmetry is absent in [1,nsym], connect q->-q+G using time reversal.
      elseif( (.not. invsym) .and. (norm2(-sxq-xq_irr) .lt. eps)) then
         q2mq = .true.
         return
      endif
   enddo

   !
   ! if the code reaches here then xq and xq_irr are not related by symmetry
   write(stdout,*) 'xq:', xq
   write(stdout,*) 'xq_irr:', xq_irr
   call errore ('find_symmetry_index',' Unable to find the symmetry operation',1)
end subroutine find_symmetry_index

subroutine read_patterns(phonon_path)
   !
   ! read patterns
   ! adapted from PH/ph_restart.f90
   !
   use ions_base,     only : nat
   use modes,         only : u, npert, nirr, name_rap_mode, num_rap_mode
   use io_global,     only : stdout, ionode
   use lr_symm_base, only : nsymq

   implicit none

   integer ::  irr, iq, iunpun
   ! counters
   integer :: ierr
   character (len=256) :: phonon_path, dirname, filename
   character (len=256), external :: trimcheck
   character(len=6) :: int_to_char
   logical :: exst

   allocate (u ( 3 * nat, 3 * nat))
   allocate (name_rap_mode( 3 * nat))
   allocate (num_rap_mode( 3 * nat))
   allocate (npert ( 3 * nat))


   do iq=1, nq_irr
      ! Read the patterns from files
      ierr=0
      u=0.d0
      if ( ionode )  call iotk_free_unit( iunpun, ierr )
      call mp_bcast (ierr, ionode_id, inter_pool_comm)
      call errore('read_patterns', 'no free units to write or read', ierr)
      !
      dirname = trimcheck ( trim( phonon_path ) // trim( prefix ) // '.phsave/' )
      filename= trim( dirname ) // 'patterns.' // trim(int_to_char(iq)) // '.xml'      
      !
      if (ionode ) then
         ierr=0
         INQUIRE( FILE=TRIM(filename), EXIST=exst )
         if (.not.exst) then
            call mp_bcast (exst, ionode_id, inter_pool_comm)
            !
            ! If the file does not exist and we must read from it, we return with
            ! an error message.
            !
            call errore('read_patterns', 'cannot open file for reading', 1)
         end if
         call iotk_open_read( iunpun, file = trim( filename ), binary = .false., IERR = ierr )
      endif
      call mp_bcast (ierr, ionode_id, inter_pool_comm)
      call errore('read_patterns', 'problems opening pattern file', ierr)
      !
      call read_modes( iq, iunpun, ierr)
      call errore('read_patterns', 'problems reading u',ierr)
      !
      IF (ionode) call iotk_close_read(iunpun)
      !
      ! set the global variables
      ! the order of the small group of q for each q (the number of symmetry of the small group)
      nsymq_iq(iq) = nsymq
      irr_iq(iq) = nirr
      do irr=1, nirr
         npert_irr_iq(irr,iq)=npert(irr)
      enddo
      pattern(:,:,iq)=u
      !write(stdout,*) iq, 'nsymq_iq, irr_iq, sum(npert)', nsymq_iq(iq), irr_iq(iq), sum(npert_irr_iq(1:irr_iq(iq),iq))
      !do irr = 1, nirr
      !   write(stdout,*) iq, irr, 'npert_irr_iq ', npert_irr_iq(irr, iq)
      !enddo
      !write(stdout,*) pattern(:,:,iq)
   enddo
   !
   deallocate (u)
   deallocate (npert)
   deallocate (num_rap_mode)
   deallocate (name_rap_mode)
   !
   return
end subroutine read_patterns

subroutine read_modes( current_iq, iunpun, ierr )
   !
   !   Adapted from PH/ph_restart.f90
   !   This routine reads the displacement patterns.
   !
   use modes, only : nirr, npert, u, name_rap_mode, num_rap_mode
   use control_ph, only : trans, zeu

   use lr_symm_base, only : nsymq, minus_q

   implicit none

   integer,          intent(in) :: current_iq, iunpun
   integer,          intent(out) :: ierr
   integer :: imode0, imode, irr, ipert, iq, iu

   logical :: exst
   !
   ierr=0
   if (ionode) then
      call iotk_scan_begin( iunpun, "IRREPS_INFO" )
      !
      call iotk_scan_dat(iunpun,"QPOINT_NUMBER",iq)
   endif
   call mp_bcast( iq,  ionode_id, inter_pool_comm )
   if (iq /= current_iq) call errore('read_modes', &
    'problems with current_iq', 1 )

   if (ionode) then

      call iotk_scan_dat(iunpun, "QPOINT_GROUP_RANK", nsymq)
      call iotk_scan_dat(iunpun, "MINUS_Q_SYM", minus_q)
      call iotk_scan_dat(iunpun, "NUMBER_IRR_REP", nirr)
      imode0=0
      do irr=1,nirr
         call iotk_scan_begin( iunpun, "REPRESENTION"// &
           TRIM( iotk_index( irr ) ) )
         call iotk_scan_dat(iunpun,"NUMBER_OF_PERTURBATIONS", npert(irr))
         do ipert=1,npert(irr)
            imode=imode0+ipert
            call iotk_scan_begin( iunpun, "PERTURBATION"// &
              TRIM( iotk_index( ipert ) ) )
            call iotk_scan_dat(iunpun,"SYMMETRY_TYPE_CODE", &
             num_rap_mode(imode))
            call iotk_scan_dat(iunpun,"SYMMETRY_TYPE", name_rap_mode(imode))
            call iotk_scan_dat(iunpun,"DISPLACEMENT_PATTERN",u(:,imode))
            call iotk_scan_end( iunpun, "PERTURBATION"// &
              TRIM( iotk_index( ipert ) ) )
         enddo
         imode0=imode0+npert(irr)
         call iotk_scan_end( iunpun, "REPRESENTION"// &
           TRIM( iotk_index( irr ) ) )
      enddo
      !
      call iotk_scan_end( iunpun, "IRREPS_INFO" )
      !
   endif

   call mp_bcast( nirr,  ionode_id, inter_pool_comm )
   call mp_bcast( npert,  ionode_id, inter_pool_comm )
   call mp_bcast( nsymq,  ionode_id, inter_pool_comm )
   call mp_bcast( minus_q,  ionode_id, inter_pool_comm )
   call mp_bcast( u,  ionode_id, inter_pool_comm )
   call mp_bcast( name_rap_mode,  ionode_id, inter_pool_comm )
   call mp_bcast( num_rap_mode,  ionode_id, inter_pool_comm )
   return
end subroutine read_modes

end module lattice_data

!===============================================================================
! Copyright (C) 2016-2020 Jin-Jian Zhou, Jinsoo Park, I-Te Lu, Marco Bernardi
!
! This program is distributed under the terms of the GNU General Public License.
! See the file `LICENSE' in the root directory of this distribution, or obtain 
! a copy of the License at <https://www.gnu.org/licenses/gpl-3.0.txt>.
!
! Author: jpark <jinsoop412@gmail.com>; jjzhou <jjchou.comphy@gmail.com>; 
! Comment:
!  adapted from QE/PHonon/PH/dfile_star
!  Formula for rotation of dvscf: Eq. C13 in Luis paper, PRB 97, 235146(2018)
!
!  N.B. sr in symm_base module are symmetry operation matrices.
!  sr is S operation in cartesian cooridnate, sr(i,j) = S_ij, 
!   and the same sr can be applied to both real and recip. space.
!
!  If if we transfer S from cartesian (c) coord. to crystal coord. (r).
!   then the matrix are different for real space and recip. space:
!  in  real sace,  x_c = A x_r, (A -> at in cell_base module) 
!  in recip space. k_c = B k_r  (B -> bg in cell_base module), and A^-1 = B^T
!  so,  
!  S in crystal coord. in real space is:  B^T * sr * A  (T means transpose)
!    S (A x_s) = A x_s'  -->  (A^-1 S A) x_s = x_s'
!   
!  S in crystal coord. in recip. space is:  A^T * sr * B
!    S (B k_s) = B k_s'  -->  (B^-1 S B) k_s = k_s'
!
!  the variable s in symm_base is: 
!    s(i,j) = (B^T * sr * A)_ji  (S in cryst. coord. in real space)
!  quivalently, 
!    s^T = (B^T * sr * A)  ==> s = A^T * sr^T * B. = A^T (sr^-1) *B
!  since sr is rotation operation (in cart.), sr^T = sr^-1
!  so, s(i,j) = (A^T (sr^-1) * B)_ij (S^-1 in cryst, coord. in recip. space)
!  in summary, s has different meaning when applied to real. and recip. space.
!
! Maintenance:
!===============================================================================

module dvscf_mod
   use kinds,     only: dp
   use io_global, only: stdout
   use cell_base, only: at, bg
   use fft_base,  only: dfftp
   use ions_base, only: nat, tau
   use noncollin_module, only: nspin_mag
   use symm_base, only: nsym, s, sr, invs, irt, ft
   implicit none
   private

   public :: obtain_dvscf

contains

! obtain self-consistent part: Eq.B30 of PRB 64 235118
subroutine obtain_dvscf(current_irrq, dvscf_irrq, iq, dvscf_iq, lnew_irrq)
   use lattice_data, only: xq_tot, identity_op, q2minus_q, pattern, xq_irr, q2irq_symop
   implicit none
   logical, intent(in) :: lnew_irrq !if .true., init dvscf_irrq (reading from file)
   integer, intent(in) :: current_irrq
   ! dvscf_irrq : dvscf of the irredicible q (xq_irr) in cartesian coordinate
   complex(dp), intent(inout) :: dvscf_irrq(dfftp%nnr, nspin_mag, 3*nat)
   !
   integer, intent(in) :: iq
   ! dvscf_iq : dvscf of xq in cartesian coordinate
   complex(dp), intent(out) :: dvscf_iq(dfftp%nnr, nspin_mag, 3*nat)
   !
   !local
   integer :: isym
   real(dp) :: xq_cart(3), xq_irr_cart(3)
   !real(dp) :: rtau(3,48,nat)

   ! xq_tot(:,iq): the iq-th point in crystal coordinate
   xq_cart = xq_tot(:,iq)
   !transfer xq_cart to cartesian coordinate
   call cryst_to_cart(1, xq_cart, bg, 1)
   !
   !call sgam_lr (at, bg, nsym, s, irt, tau, rtau, nat)
   !
   !if lnew_irrq is .true., init dvscf_irrq by reading from file
   ! dvscf_irrq in cartesian coordinates
   if(lnew_irrq) call load_dvscf(current_irrq, pattern(:,:,current_irrq), dvscf_irrq)
   !write(stdout,'(a, 2i5, 3f12.8)') "Finish loading dvscf:", iq, irrq, xq
   
   ! xq_irr(:,irrq): the irrq-th irreducible q-point in crystal coordinate
   xq_irr_cart = xq_irr(:, current_irrq)
   ! transfer to cartesian coordinate
   call cryst_to_cart(1, xq_irr_cart, bg, 1)
   
   ! sr(:,:,isym) is the symmetry operation S that Sq = q_irr
   isym = q2irq_symop(iq)
   if( (isym .eq. identity_op) .and. (.not. q2minus_q(iq)) ) then
      write(stdout,'(a15, i3, 2x, 3f12.8)') "irreducible q: ", current_irrq, xq_irr_cart(:)
      dvscf_iq(:,:,:) = dvscf_irrq(:,:,:)
      !
      !uact(:, :) = pattern(:, :, current_irrq)
   else
      !write(stdout,'(3f12.8,2x,3f12.8,2x,i2,1x,l)') xq_cart, xq_irr_cart(:), isym, q2minus_q(iq)
      !dvscf in cartesian coordinates for xq, rotated from dvscf_irrq
      call get_dvscf_star(xq_cart, xq_irr_cart, isym, q2minus_q(iq), dvscf_irrq, dvscf_iq)

      !!rotate displacement pattern
      !isym_inv = invs(isym)
      !call rotate_mod( pattern(:,:,current_irrq), uact, sr(:,:,isym_inv), &
      !      irt, rtau, xq_irr_cart, nat, isym)
      !! u(q) = conjg( u(-q) )
      !if(q2minus_q(iq))  uact = conjg(uact)
   endif

   !!transfer to pattern coordinate
   !nrow = dfftp%nnr * nspin_mag
   !call zgemm('N','N', nrow, nmodes, nmodes, cone, &
   !         dvscf_tmp(1,1,1), nrow, uact, nmodes, czero, dvscf(1,1,1), nrow)
   
end subroutine obtain_dvscf


subroutine load_dvscf(irrq, u, dvscf)
   use input_param, only: prefix, phdir
   implicit none
   integer,     intent(in) :: irrq
   complex(dp), intent(in) :: u(3*nat, 3*nat)
   complex(dp), intent(out):: dvscf(dfftp%nnr, nspin_mag, 3*nat)
    
   logical :: exst
   integer :: rec_len, iunit, imod, nmodes, nrow
   character(len=256) :: fname
   complex(dp), allocatable :: dvscf_tmp(:,:,:)

   complex(dp), parameter :: czero = (0.0_dp, 0.0_dp)
   complex(dp), parameter :: cone  = (1.0_dp, 0.0_dp)
   ! functions:
   integer, external :: find_free_unit
   CHARACTER(LEN=6), EXTERNAL :: int_to_char

   fname = trim(phdir) // trim(prefix) // '.dvscf_q' // trim(int_to_char(irrq))
 !  when the file is read, it must be locked - invoke mpi barrier later
   inquire(file=trim(fname), exist = exst)
   if(.not. exst) call errore('load_dvscf',trim(fname)//' does not exist!',1)
   !
   nmodes  = 3 * nat
   rec_len = 2 * dfftp%nnr * nspin_mag
   !
   ! Open the dvscf file for reading
   iunit = find_free_unit()
   call open_dvscf(iunit, fname, rec_len, nmodes, exst)
   !
   allocate( dvscf_tmp(dfftp%nnr, nspin_mag, nmodes) )
   dvscf_tmp = czero

   do imod = 1, nmodes
      !read in dvscf for all the irreps
      call davcio(dvscf_tmp(:,:,imod), rec_len, iunit, imod, -1)
   end do
   close(iunit)
   
   !Transform from the basis of the patterns to cartesian basis
   ! dvscf(:,:,i) = dvscf(:,:,i) + conjg(u(i,j))*dvscf_tmp(:,:,j) 
   !        for i,j loop over [1, nmodes]
   ! a more efficient way to do the transformation
   nrow = dfftp%nnr * nspin_mag
   call zgemm('N','C', nrow, nmodes, nmodes, cone, &
               dvscf_tmp, nrow, u, nmodes, czero, dvscf, nrow)

   deallocate( dvscf_tmp )
end subroutine load_dvscf


!-----------------------------------------------------------------------
!   rotate dvscf_xqirr to get dvscf_xq
!-----------------------------------------------------------------------
subroutine get_dvscf_star(xq, xq_irr, isym_inv, q2mq, dvscf_irr, dvscf_xq)
   use constants,        only: tpi
   ! #jjzhou: only support nspin_mag = 1 or 2;  nspin_mag = 4 is not supported
   !  nspin_mag = 1: spinless,  or non-colliner without magnetic (domag=.false.)
   !  nspin_mag = 2: LSDA spin-polarized calculations.
   ! (Note: under this constrain, nspin_mag is equivalent to nspin from lsda_mod)
   !
   implicit none
   ! input variables:
   ! xq     : q vector, in cartesian coordinate
   ! xq_irr : q_irr vector, in cartesian coordinate
   real(dp),intent(in) :: xq(3), xq_irr(3)
   ! index of symmetry operation that rotate q -> qirr, e.g. sr(:,:,isym_inv)
   !  is the symmetry operation S^-1 that S^-1 q = q_irr (so q = S q_irr)
   integer, intent(in) :: isym_inv
   ! if q2mq is ture, then S^-1 (-q) = q_irr, 
   ! q and q_irr is connect via S^-1 and time reversal symmetry.
   logical, intent(in) :: q2mq
   ! dvscf_irr : dvscf of the irredicible q (xq_irr) in cartesian coordinate
   complex(dp), intent(in) :: dvscf_irr(dfftp%nnr, nspin_mag, 3*nat)
   ! dvscf_irr : dvscf of xq in cartesian coordinate
   complex(dp), intent(out) :: dvscf_xq(dfftp%nnr, nspin_mag, 3*nat)
   !
   ! local variables
   integer :: na, i, j, nmodes, imode0, ijks, nrpt
   integer :: is, k, n, nn, ri, rj, rk, ipol, index0, nar, ftau(3)
   ! auxiliary xq\cdot\tau and \xq_s\cdot\tau
   real(dp) :: xq_tau, xqirr_tau, xq_tmp(3), xq_rot(3)
   character(len=6), external :: int_to_char
   !
   complex(dp) :: phase_xqirr
   complex(dp), allocatable :: dvscf_tmp(:,:,:), phase_xq(:)

   nmodes = 3*nat
   !
   if(isym_inv < 1 .or. isym_inv > 48) call errore('get_dvscf_star', &
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

   !NOTE: dfftp%nnr = dfftp%nr1x * dfftp%nr2x * dfftp%nr3x 
   ! since we enforce nproc = 1 for each pool
   allocate( dvscf_tmp(dfftp%nnr, nspin_mag, nmodes) )
   !
   ! Transform to crystalline coordinates (necessary in order to apply s)
   dvscf_tmp =cmplx(0._dp, 0._dp, kind=dp)
   do i = 1, nat
      na = (i-1) * 3
      do j = 1, 3
         dvscf_tmp(:,:,na+j) = dvscf_irr(:,:,na+1)*at(1,j) + &
                               dvscf_irr(:,:,na+2)*at(2,j) + &
                               dvscf_irr(:,:,na+3)*at(3,j)
      enddo
   enddo
   !
   ! take away the phase due to the q-point
   do i = 1, nat
     !
     xqirr_tau = tpi * dot_product(xq_irr(:), tau(:,i))
     phase_xqirr= cmplx( cos(xqirr_tau), sin(xqirr_tau), kind=dp)
     !
     do ipol=1,3
        imode0 = (i-1)*3 + ipol
        dvscf_tmp(:,:,imode0) = phase_xqirr * dvscf_tmp(:,:,imode0)
     enddo
   enddo
   !
   ! Now rotate the dvscf
   allocate( phase_xq(nat) )
   !
   do i = 1, nat  
      xq_tau = - tpi * dot_product(xq_rot(:), tau(:,i))
      phase_xq(i) = cmplx(cos(xq_tau), sin(xq_tau), kind=dp)
   enddo
   
   ftau(1) = nint( ft(1,isym_inv) * dfftp%nr1 )
   ftau(2) = nint( ft(2,isym_inv) * dfftp%nr2 )
   ftau(3) = nint( ft(3,isym_inv) * dfftp%nr3 )
   !
   !do is = 1, nspin_mag
   !  kloop : do k = 1, dfftp%nr3
   !    jloop : do j = 1, dfftp%nr2
   !      iloop : do i = 1, dfftp%nr1

   dvscf_xq = cmplx(0._dp, 0._dp, kind=dp)
   nrpt = dfftp%nr1 * dfftp%nr2 * dfftp%nr3

!$omp parallel do schedule(guided) default(shared) private(ijks, is, n, &
!$omp&  i, j, k, ri, rj, rk, nn, na, nar, index0, ipol, imode0)
   do ijks = 1, nspin_mag * nrpt
      is = (ijks-1) / nrpt + 1
      !
      n = mod( (ijks-1), nrpt ) + 1
      !n  = (i-1)  + (j-1)*dfftp%nr1  + (k-1)*dfftp%nr2*dfftp%nr1  + 1
      i = mod(n-1, dfftp%nr1) + 1
      j = mod( (n-1)/dfftp%nr1, dfftp%nr2) + 1
      k = (n-1) / (dfftp%nr1 * dfftp%nr2) + 1
           
      ! ruotaijk find the rotated of i,j,k with the inverse of S
      ! a.k.a ruotaijk perform [S^-1] (i, j, k) = (ri, rj, rk)
      call ruotaijk( s(:, :, isym_inv), ftau, i, j, k, &
                    dfftp%nr1, dfftp%nr2, dfftp%nr3, ri, rj, rk)
      !
      nn = (ri-1) + (rj-1)*dfftp%nr1 + (rk-1)*dfftp%nr2*dfftp%nr1 + 1
      !
      do na = 1, nat
        nar = irt(isym_inv, na)
        !irt is the list of rotated of each atom: nar = [S^-1] na
        index0 = (nar-1) * 3
        !
        do ipol = 1, 3
            imode0 = (na-1)*3 + ipol
            !N.B. s(a, b, isym_inv) ==> [S^-1]_ba in C13 of PRB97,235146
            dvscf_xq(n, is, imode0) = dvscf_xq(n, is, imode0) + &
                ( s(ipol, 1, isym_inv) * dvscf_tmp(nn, is, index0+1) + &
                  s(ipol, 2, isym_inv) * dvscf_tmp(nn, is, index0+2) + &
                  s(ipol, 3, isym_inv) * dvscf_tmp(nn, is, index0+3) )
        enddo
      enddo
   enddo
!$omp end parallel do
   !
   ! Add back the phase factor for the new q-point
   !
   dvscf_tmp=cmplx(0._dp, 0._dp, kind=dp)
   do na = 1, nat
     do ipol = 1, 3
       imode0 = (na-1)*3 + ipol
       dvscf_tmp(:,:,imode0 ) = dvscf_xq(:,:,imode0) * phase_xq(na)
     enddo
   enddo
   !
   ! Back to cartesian coordinates
   !
   dvscf_xq = cmplx(0._dp, 0._dp, kind=dp)
   do i = 1, nat
     imode0 = (i-1)*3
     do j = 1, 3
       dvscf_xq(:,:,imode0+j) = dvscf_tmp(:,:,imode0+1)*bg(j,1) + &
             dvscf_tmp(:,:,imode0+2)*bg(j,2) + dvscf_tmp(:,:,imode0+3)*bg(j,3)
     enddo
   enddo
   
   !since dv/du(q) = \sum_{R} exp(-iq*R) dv/du(R) and dv/du(R) is a real function, 
   ! so dv/du(-q) = conjg( dv/du(q) ) [see Eq.(25) in the perturbo paper]
   ! 
   !dvscf(q) = conjg(dvscf(-q))
   if(q2mq) dvscf_xq = conjg( dvscf_xq )
   !
   !N.B.: the above is not based on time reversal symmetry, 
   ! so in case of LSDA without T-reversal symmetry, it still holds
   !  for spin-up and spin-down channels separately.
   
   deallocate(dvscf_tmp, phase_xq)
   !
end subroutine get_dvscf_star


! Modified version of subroutine diropn in QE/Module/io_file.f90
subroutine open_dvscf(iunit, fname, reclen, tot_rec, exst)
   implicit none
   integer, intent(in) :: iunit, reclen, tot_rec
   !input: unit of the file to open
   !input: length of the records
   !input: total number of records
   character(len=*), intent(in) :: fname
   logical, intent(out) :: exst
   ! output: if true the file exists

   real(dp):: dummy
   integer*8 :: unf_recl, file_size
   ! double precision to prevent integer overflow
   integer :: ios, direct_io_factor
   logical :: opnd
   !
   if (iunit < 0) call errore ('open_dvscf', 'wrong unit', 1)
   inquire( unit = iunit, opened = opnd )
   if (opnd) call errore ('open_dvscf', "can't open a connected unit", abs(iunit))

   inquire (file =trim(fname), exist = exst)
   if(.not. exst) call errore('open_dvscf',trim(fname)//' does not exist!',1)
   
   if (reclen == -1) call errore('open_dvscf','negative record length',1)

   ! the  record length in direct-access I/O is given by the number of
   ! real*8 words times direct_io_factor (may depend on the compiler)
   !
   INQUIRE (IOLENGTH=direct_io_factor) dummy
   unf_recl = direct_io_factor * int(reclen, kind=kind(unf_recl))
   if (unf_recl <= 0) call errore ('open_dvscf', 'wrong record length', 3)
   !
   inquire(file=trim(fname), size=file_size)
   if( file_size .ne. tot_rec*unf_recl ) then
      write(stdout,'(5x,a)') "Warning (open_dvscf): " //trim(fname)// &
         " file size does not match. (turn on -assume byterecl if use intel compiler)"
   endif

   ios = 0
   open (iunit, file=trim(adjustl(fname)), iostat=ios, form='unformatted', &
        status = 'unknown', access = 'direct', recl = unf_recl, action='read')
   if (ios /= 0) call errore ('open_dvscf', 'error opening '//trim(fname), iunit)
   return
end subroutine open_dvscf

end module dvscf_mod

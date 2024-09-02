!===============================================================================
! Copyright (C) 2016-2020 Jin-Jian Zhou, Jinsoo Park, I-Te Lu, Marco Bernardi
!
! This program is distributed under the terms of the GNU General Public License.
! See the file `LICENSE' in the root directory of this distribution, or obtain 
! a copy of the License at <https://www.gnu.org/licenses/gpl-3.0.txt>.
!
! Author: jjzhou <jjchou.comphy@gmail.com>
! Comment:
!
! Maintenance:
!===============================================================================

module phonon_dispersion
   use pert_const, only: dp, czero, twopi, ci, cone, e2, pi
   use pert_utils, only: hmat_upper_diag, get_exp_ikr
   use wigner_seitz_cell, only: ws_cell
   use polar_correction, only: set_polar_parameter, dyn_mat_longrange
   use force_constant, only: lattice_ifc, force_const, set_ws_cell_ph
   use qe_mpi_mod, only: ionode, mp_bcast, inter_pool_comm, ionode_id, ionode, stdout
   implicit none
   public

   type(lattice_ifc), save :: phon
   
   public :: init_lattice_ifc !, init_lattice_ifc_tdep
   public :: solve_phonon_modes
contains

subroutine init_lattice_ifc(file_id, qdim, ph)
   use hdf5_utils
   use epwan_hdf5_io, only: read_force_constant
   use pert_data,  only: zstar, epsil, nat, at, bg, &
                           tau, volume, tpiba, thickness_2d, loto_alpha, lpolar
   implicit none
   integer(HID_T), intent(in) :: file_id
   integer, intent(in) :: qdim(3)
   type(lattice_ifc), intent(out) :: ph
   !
   integer :: i, nelem
   real(dp), allocatable :: cryst_tau(:,:)

   allocate( cryst_tau(3,nat) )
   !transform atomic position to crystal coordinate
   cryst_tau(:,:) = tau(:,:)
   call cryst_to_cart(nat, cryst_tau, bg, -1)

   call set_ws_cell_ph(ph, qdim, nat, at, cryst_tau)
   if(ionode) call read_force_constant(file_id, ph)
   !bcast
   nelem = ph%na * (ph%na + 1) / 2
   do i = 1, nelem
      call mp_bcast(ph%phi(i)%ifc(:,:,:), ionode_id, inter_pool_comm)
   enddo
   
   ph%lpol = lpolar
   ph%l_2d = merge(.true., .false., thickness_2d > 0.0_dp)
   !setup for polar correction.
   !N.B. to be consistent with rigid_bulk, alpha=1.0 is enforced
   if(ph%lpol) call set_polar_parameter(qdim, nat, volume, tpiba, bg, &
                        epsil, zstar, tau, ph%pol, thickness_2d, loto_alpha)
   
   deallocate( cryst_tau )
end subroutine init_lattice_ifc


subroutine solve_phonon_modes(ph, xqt, pheig, phmode)
   use pert_data,  only: mass
   implicit none
   type(lattice_ifc), intent(in) :: ph
   real(dp), intent(in) :: xqt(3)
   real(dp), intent(out) :: pheig(:) !pheig(ph%nm)
   complex(dp), intent(out), optional :: phmode(:,:) !phmode(ph%nm, ph%nm)
   !local variables
   integer :: i, j , ia, ja, ii, jj, m, n, ir, nmodes, nelem, irp
   real(dp) :: mfactor, w2( 3*ph%na ), tmp( 3*ph%na )
   complex(dp), allocatable :: dmat_lr(:,:,:), dmat(:,:,:), ev(:,:), dyn_upper(:), exp_ikr(:)
   type(ws_cell), pointer :: ws
   
   nmodes = 3 * ph%na
   nelem = ph%na * (ph%na + 1) / 2
   !
   allocate( dmat(3,3,nelem), dyn_upper( (nmodes*(nmodes+1))/2 ), exp_ikr(ph%nrvec) )
   !
   call get_exp_ikr(xqt, ph%rvec_set, exp_ikr)
   dmat = czero
   do m = 1, nelem
      ws => ph%phi(m)%ws_ph

      do ir = 1, ws%nr
         irp = ws%rvec(ir)
         ! the 1/ws%ndeg(ir) is already included in ifc(:) !/ real(ws%ndeg(ir), dp)
         !dmat(:,:, m) = dmat(:,:, m) + exp(iqR) * ph%phi(m)%ifc(:,:, ir)
         call zaxpy(9, exp_ikr(irp), ph%phi(m)%ifc(1,1,ir), 1, dmat(1,1,m), 1)
      enddo
   enddo
   
   !apply polar correction if needed.
   if( ph%lpol ) then
      allocate( dmat_lr(3, 3, nelem) )
      call dyn_mat_longrange(ph%pol, xqt, dmat_lr)
      dmat = dmat + dmat_lr
      deallocate(dmat_lr)
   endif
   
   dyn_upper = czero
   !create upper trianglar of the dynamical matrix
   do jj = 1, nmodes
   do ii = 1, jj
      !index of the upper triangluatr dynamical matrix
      n = (jj * (jj-1)) / 2 + ii
      
      !atomic index
      ia = (ii - 1) / 3 + 1
      ja = (jj - 1) / 3 + 1
      !cartesian directions
      i = mod(ii-1, 3) + 1
      j = mod(jj-1, 3) + 1
      
      m = ( ja * (ja-1) ) / 2 + ia
      mfactor = 1.0_dp / sqrt( mass(ia) * mass(ja) )

      if(ia .ne. ja) then
         dyn_upper(n) = dmat(i,j,m) * mfactor
      else
         !enforce herminicity
         dyn_upper(n) = (dmat(i,j,m) + conjg(dmat(j,i,m))) * 0.5_dp * mfactor
      endif
   enddo; enddo
   
   if( present(phmode) ) then
      allocate( ev(nmodes, nmodes) )
      call hmat_upper_diag(dyn_upper, nmodes, w2, ev)
      !return eigen-displacement if requested.
      do ii = 1, nmodes
         ia = (ii - 1) / 3 + 1
         phmode(ii, 1:nmodes) = ev(ii, 1:nmodes) / sqrt( mass(ia) )
      enddo
      deallocate( ev )
   else
      call hmat_upper_diag(dyn_upper, nmodes, w2)
   endif
   !compute phonon frequencies: omega = sqrt(w2)
   tmp = sqrt( abs(w2) )
   ! if w2(i) < 0, return negative frequency.
   pheig(1:nmodes) = merge(tmp, -tmp, w2 > 0.0_dp)
   
   deallocate(dmat, dyn_upper)
end subroutine solve_phonon_modes



end module phonon_dispersion

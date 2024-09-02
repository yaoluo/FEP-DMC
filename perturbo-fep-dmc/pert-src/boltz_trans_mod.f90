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

module boltz_trans_mod
   use pert_const, only: dp
   use pert_output, only: stopwatch
   use boltz_grid, only: grid
   use pert_utils, only: find_free_unit, fermi, mfermi_deriv, mat_inv3
   use qe_mpi_mod, only: mp_split_pools, mp_sum, mp_barrier, inter_pool_comm
   private
   public :: trans_dos_calc !compute density of states
   public :: trans_tdf_calc
   public :: trans_density_calc
   public :: trans_cond_calc
   public :: extract_trans_coeff
contains

! compute dos using tetrahedra integration, spin is taken into account here
subroutine trans_dos_calc(kg, ene, dos)
   implicit none
   type(grid), intent(in) :: kg
   real(dp), intent(in) :: ene(:) 
   real(dp), intent(out) :: dos(:)
   !local variables
   integer :: nstep, nb, it, ic, ie, it_st, it_end
   real(dp) :: fnk(kg%numb,4), et(kg%numb,4)
   real(dp), allocatable :: dos_t(:)
   ! get the size of input array
   nstep = size(ene);   nb = kg%numb
   dos(1:nstep) = 0.0_dp;  fnk = 1.0_dp

   !mpi parallel: distribute tetras over pools
   call mp_split_pools(kg%num_tet, it_st, it_end)
!$omp parallel default(shared) private(it, ic, et, dos_t, ie)
   allocate( dos_t(nstep) );  dos_t = 0.0_dp
!$omp do schedule(guided)
   do it = it_st, it_end
      ! prepare data for the four vertex.
      do ic = 1, 4
         ! tet(2, ic, it): idx of kgrid. full klist
         et(:,ic) = kg%enk(:, kg%tetra(2,ic,it))
      enddo
      !compute the contribution from current tetra, and add it to dos_t
      !NOTE: dos_t is inout. initialization is required for dos_t
      call tetra_int(nb, et, fnk, nstep, ene, dos_t)
   enddo
!$omp end do nowait
   !collect results from each thread
   do ie = 1, nstep
!$omp atomic
      dos(ie) = dos(ie) + dos_t(ie)
   enddo
   deallocate(dos_t)
!$omp end parallel

   !DOS in the units of 'states/Ryd/unit.cell.'
   !tweight = 1/tot_tet_fbz = V_T/V_G, (spin degeneracy is not included here)
   dos = dos * kg%tweight
   ! collect results from different pools
   call mp_sum(dos, inter_pool_comm)
   !call mp_barrier(inter_pool_comm)
end subroutine trans_dos_calc

!compute transport distribution function (TDF)
! TDF_ij(E) = spinor/(N*V) sum_{nk} vnk * vnk * tau * delta(E-e_nk)
! the unit of velocity is Ryd * Bohr * alat, vnk = 1/hbar*dE/dk.
! However, here we does not include spinor and volume V. Be careful!
subroutine trans_tdf_calc(kg, ene, mfd, tdf)
   use boltz_utils, only: velprod
   implicit none
   type(grid), intent(in) :: kg
   !mfd: mean free displacement = vnk*tau, mfd(3 or 6,numb,kg%nk)
   real(dp), intent(in) :: ene(:), mfd(:,:,:)
   real(dp), intent(out) :: tdf(:,:) !tdf(6 or 12,:)
   !local variables
   integer :: nb, nstep, it, ic, ik, ib, ia, ie, it_st, it_end, ncomp
   real(dp) :: et(kg%numb, 4)
   real(dp), allocatable :: tdf_t(:,:), fnk(:,:,:) !fnk(kg%numb, 4, 6 or 12)

   nb = kg%numb;  nstep = size(ene);   ncomp = size(tdf, 1)
   ! array check
   if( size(mfd, 1) .ne. (ncomp/2) ) &
      call errore('trans_tdf_calc', 'mismatch dimension in input arguments!',1)
   if( size(mfd, 2) .ne. nb ) &
      call errore('trans_tdf_calc','array dimension mismatch', 1)

   tdf(:,:) = 0.0_dp
   !mpi parallel: distribute tetras over pools
   call mp_split_pools(kg%num_tet, it_st, it_end)
!$omp parallel default(shared) private(it, ic, et, ik, ib, fnk, ia, tdf_t, ie)
   allocate( tdf_t(nstep, ncomp), fnk(kg%numb, 4, ncomp) );  tdf_t = 0.0_dp
!$omp do schedule(guided)
   do it = it_st, it_end
      ! prepare data for the four vertex.
      do ic = 1, 4
         ! tet(2, ic, it): idx of kgrid. full klist
         ik = kg%tetra(2,ic,it)
         et(:,ic) = kg%enk(:,ik)
         ! get vnk_a * mfd_b
         do ib = 1, nb
            call velprod2( kg%vnk(1:3,ib,ik), mfd(:,ib,ik), fnk(ib,ic,:) )
         enddo
      enddo
      do ia = 1, ncomp
         call tetra_int(nb, et, fnk(:,:,ia), nstep, ene, tdf_t(:,ia))
      enddo
   enddo
!$omp end do nowait
   !collect results from each thread
   do ie = 1, nstep
      do ia = 1, ncomp
!$omp atomic
         tdf(ia, ie) = tdf(ia, ie) + tdf_t(ie, ia)
      enddo
   enddo
   deallocate(tdf_t, fnk)
!$omp end parallel

   !NOTE: volume and spinor is not included here.!!!
   !tweight = 1/tot_tet_fbz = V_T/V_G,
   tdf = tdf * kg%tweight
   ! collect results from different pools
   call mp_sum(tdf, inter_pool_comm)
   !call mp_barrier(inter_pool_comm)
   return
   !
   contains
      subroutine velprod2(v1, vv2, prod2)
         implicit none
         real(dp), intent(in) :: v1(3)
         real(dp), intent(in) :: vv2(:)
         real(dp), intent(out) :: prod2(:)
         !
         integer :: nsize, i

         nsize = size(vv2)
         do i = 0, nsize-1, 3
            prod2( (i*2+1):(i*2+6) ) = velprod(v1, vv2( (i+1):(i+3) ))
         enddo
      end subroutine velprod2
   ! 
end subroutine trans_tdf_calc

!compute carrier concentration: number of carriers per unit cell
!  or determine chemical potential at given concentration.
subroutine trans_density_calc(kg, ene, tmpr, ef, hole, ncol, vol, dens, find_ef, dos)
   implicit none
   type(grid), intent(in) :: kg
   !hole: true -> hole carrier; false -> electron carrier
   logical, intent(in) :: hole, ncol  !ncol: non-colliner
   real(dp), intent(in) :: ene(:), tmpr(:), vol !tmpr: temperature, volume
   !dens: density in Rydberg atomic unit: #.carrier / bohr^3
   real(dp), intent(inout) :: ef(:), dens(:)
   logical, intent(in), optional :: find_ef !if true, calc ef from input density
   real(dp), intent(out), optional :: dos(:)  ! density of states
   ! local variables
   logical :: less
   integer :: i, nestep, ntmpr
   real(dp) :: mid_dens, ist, iend, dist, mid, estep
   real(dp), allocatable :: dos_tmp(:), occup(:)

   call stopwatch('trans_density_calc','start')
   nestep = size(ene)
   allocate( dos_tmp(nestep), occup(nestep) )
   ntmpr = size(tmpr)
   estep = (ene(nestep)-ene(1)) / (nestep-1.0_dp) ! energy step
   ! some check
   if( (size(dens) .ne. ntmpr) .or. (size(ef) .ne. ntmpr)) &
      call errore('trans_density_calc','array dimension mismatch', 1)
   if( (ene(nestep)-ene(1)) < 1.0E-4_dp ) &
      call errore('trans_density_calc','energy range is too small', 1)
   ! compute density of states
   call trans_dos_calc(kg, ene, dos_tmp)
   !account for spin-degeneracy
   dos_tmp = dos_tmp * merge(1.0_dp, 2.0_dp, ncol)
   !  
   if( present(dos) ) then
      if( size(dos) .ne. nestep ) &
         call errore('trans_density_calc',"argument 'dos': mismatch dimension.", 1)
      !DOS in the units of 'states/Ryd/unit.cell.'
      dos(1:nestep) = dos_tmp(1:nestep)
   endif
   ! find ef with density close to dens
   if( present(find_ef) .and. find_ef ) then
      do i = 1, ntmpr
         ! skip if dens is unavailable
         if( dens(i) < 1.0E-15_dp ) cycle
         ! do binary search
         ist = ene(1);  iend = ene(nestep);  dist = iend - ist;
         do while ( dist > 0.5E-4_dp )
            mid = (ist + iend)*0.5_dp
            occup = fermi(tmpr(i), merge(mid-ene, ene-mid, hole))
            mid_dens = dot_product(dos_tmp, occup) * estep / vol
            !if density matches, exit the loop
            if( abs(mid_dens-dens(i))/dens(i) < 1.0E-5_dp ) exit
            !which direction to search
            less = .not. hole .and. (dens(i) > mid_dens)
            less = less .or. (hole .and. (dens(i) < mid_dens) )
            if(less) then
               ist = mid
            else
               iend = mid
            endif
            dist = iend - ist
         enddo
         ef(i) = mid
      enddo
   endif
   ! compute carrier density: int_dos(E)*f(E)*dE 
   do i = 1, ntmpr
      occup = fermi(tmpr(i), merge(ef(i)-ene, ene-ef(i), hole))
      ! carrier density: #. carrier / bohr^3
      dens(i) = dot_product(dos_tmp, occup) * estep / vol
   enddo
   deallocate(dos_tmp, occup)
   !
   call stopwatch('trans_density_calc','stop')
end subroutine trans_density_calc

!calculate conductivity and mobility using tetrahedron integration. 
subroutine trans_cond_calc(kg, ene, tmpr, ef, ncol, vol, mfd, cond, tdf)
   implicit none
   type(grid), intent(in) :: kg
   logical,  intent(in) :: ncol
   real(dp), intent(in) :: tmpr, ef, ene(:), mfd(:,:,:), vol
   real(dp), intent(out) :: cond(6)
   real(dp), intent(out), optional :: tdf(:,:)
   ! local variables
   integer :: i, nestep, ncomp
   real(dp) :: estep
   real(dp), allocatable :: tdf_tmp(:,:)

   call stopwatch('trans_cond_calc','start')
   nestep = size(ene)
   ncomp = size(mfd, 1) * 2
   allocate( tdf_tmp(ncomp, nestep) )
   estep = (ene(nestep)-ene(1)) / (nestep-1.0_dp)
   ! compute TDF_ij(E) = spinor/(N*V) sum_{nk} vnk*vnk*tau*delta(E-e_nk)
   ! the unit of velocity is Ryd*Bohr*alat, vnk = 1/hbar*dE/dk.
   ! NOTE: tdf computed below does not include spinor and volume V !!!
   call trans_tdf_calc(kg, ene, mfd, tdf_tmp)
   !account for spin-degeneracy and volume
   tdf_tmp = tdf_tmp * merge(1.0_dp, 2.0_dp, ncol) / vol
   !
   if( present(tdf) ) then
      if( size(tdf, 1) .ne. ncomp .or. size(tdf, 2) .ne. nestep) &
         call errore('boltz_trans_mod',"argument 'tdf': mismatch dimension.", 1)
      !N.B.: tdf is in Rydberg atomic unit (omitted a factor of (alat)^2).
      tdf(:,:) = tdf_tmp(:,:)
   endif
   !
   !do the integration over energy to compute conductivity.
   do i = 1, 6
      !sigma = \sum_{nk} F_nk * vnk * -df/dE
      tdf_tmp(i,:) = tdf_tmp(i,:) * mfermi_deriv(tmpr, (ene-ef))

      !sigma = q^2*nspin/(N*V) sum_{nk} v_nk*v_nk*tau_nk* f *(1-f)/kT
      ! 1. q^2 is omitted in the calculation of cond(i), 
      ! 2. all the quantities are in Rydberg atomic unit.
      ! 3. (alat)^2 should be applied to cond(i) to get the actual sigma
      cond(i) = sum(tdf_tmp(i,:)) * estep
   enddo
   deallocate(tdf_tmp)
   !
   call stopwatch('trans_cond_calc','stop')
end subroutine trans_cond_calc


! compute seebeck and thermal conductivity
subroutine extract_trans_coeff(tmpr, cond, cond_seebeck, kk_coef, seebeck, t_cond, alpha)
   implicit none
   real(dp), intent(in) :: tmpr, cond(6), cond_seebeck(6), kk_coef(6)
   real(dp), intent(out) :: seebeck(6), t_cond(6)
   real(dp), intent(in), optional :: alpha(6)
   ! local variables
   real(dp) :: cs_tmp(3,3), c_tmp(3,3), s_tmp(3,3), c_aux(3,3)
   
   ! convert rank-2 tensor to matrix
   call tensor2matrix(cond, c_tmp)
   call tensor2matrix(cond_seebeck, cs_tmp)
   
   ! seebeck = sigma^-1 * [sigma*seebeck]
   c_aux = mat_inv3( c_tmp )
   s_tmp = matmul(c_aux, cs_tmp)
   !back to tensor
   call matrix2tensor(s_tmp, seebeck)

   !therm_cond = K - T * [sigma*seebeck] * seebeck
   ! if alpha is present, then replace "T * [sigma*seebeck]" with alpha
   if(present(alpha)) then
      call tensor2matrix(alpha, cs_tmp)
   else
      cs_tmp = tmpr * cs_tmp
   endif
   c_aux = matmul(cs_tmp, s_tmp)
   !
   call matrix2tensor(c_aux, t_cond)
   t_cond(1:6) = kk_coef(1:6) - t_cond(1:6)
   !
end subroutine 


! convert rank-2 tensor to matrix
subroutine tensor2matrix(tensor, matrix)
   implicit none
   real(dp), intent(in) :: tensor(6)
   real(dp), intent(out) :: matrix(3,3)
   !
   integer :: i, j, n

   matrix = 0.0_dp
   n = 0
   do j = 1, 3
   do i = 1, j
      !(1,1), (1,2), (2,2), (1,3), (2,3), (3,3)
      n = n + 1
      matrix(i,j) = tensor(n)
      !
      if(i .ne. j) matrix(j,i) = tensor(n)
   enddo; enddo
end subroutine


! convert 2D matrix to rank-2 tensor
subroutine matrix2tensor(matrix, tensor)
   implicit none
   real(dp), intent(in) :: matrix(3,3)
   real(dp), intent(out) :: tensor(6)
   !
   integer :: i, j, n
   real(dp) :: rtmp

   tensor = 0.0_dp
   n = 0
   do j = 1, 3
   do i = 1, j
      !(1,1), (1,2), (2,2), (1,3), (2,3), (3,3)
      n = n + 1
      tensor(n) = (matrix(i,j) + matrix(j,i)) * 0.5_dp
      !sanity check
      !rtmp = ( abs(matrix(i,j)) + abs(matrix(j,i)) ) * 0.5_dp
      !if( (i.ne.j) .and. rtmp>1.0E-16_dp .and. abs(matrix(i,j)-matrix(j,i)) > 1.0E-6_dp*rtmp ) &
      !   call errore('matrix2tensor','only symmetric matrix can be converted to tensor.',1)
   enddo; enddo
end subroutine

end module boltz_trans_mod

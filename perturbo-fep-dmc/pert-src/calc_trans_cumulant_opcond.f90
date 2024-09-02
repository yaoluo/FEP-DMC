!===============================================================================
! Copyright (C) 2016-2020 Jin-Jian Zhou, Jinsoo Park, I-Te Lu, Marco Bernardi
!
! This program is distributed under the terms of the GNU General Public License.
! See the file `LICENSE' in the root directory of this distribution, or obtain
! a copy of the License at <https://www.gnu.org/licenses/gpl-3.0.txt>.
!
! Author: jjzhou <jjchou.comphy@gmail.com>
! Comment:
!    compute optical conductivity using spectral function from cumulant approach
!  sigma_ab(E) = q^2*hbar * nspin*pi / (nk*vol) &
!  sum_{nk} v_a(nk)*v_b(nk)*{ Int [(f(w)-f(w+E))/E] * A_nk(w) * A_nk(w+E) * dw }
!
! Maintenance:
!===============================================================================

subroutine calc_trans_cumulant_opcond()
   use pert_const, only: dp, ryd2ev, pi, unitcharge, timeunit, bohr2ang, ryd2mev
   use boltz_grid, only: grid, init_boltz_grid_sym, boltz_grid_load, sum_equivalent_kpoint
   use pert_utils, only: find_free_unit, mfermi_deriv, fermi
   use pert_data,  only: alat, volume, epwan_fid, kc_dim
   use band_structure, only: electron_wann, init_electron_wann
   use qe_mpi_mod, only: ionode, stdout, inter_pool_comm, mp_split_pools, mp_sum
   use pert_param, only: ntemper, band_min, band_max, prefix, cum_inner_emin, &
      cum_inner_emax, cum_outer_emin, cum_outer_emax, cum_de, cum_outer_np, &
      spectral_emin, spectral_emax, spectral_np, boltz_emax, spinor, efermi, temper
   use cumulant_utility, only: spectral_cut_tail, cum_wgrid, init_cum_wgrid, load_spectral
   implicit none
   type(grid) :: kg
   type(cum_wgrid) :: wg
   real(dp), parameter :: sp_cut = 1.0E-16_dp  !unit: Ryd^-1;
   real(dp), parameter :: weight_cut = 6.0E-4_dp  !unit: Ryd^-1;
   integer  :: it, ib, ik, numb, numk, spf_lw, spf_up, i, j, ie, &
      jst, jend, iout, nestep, cut_lw
   real(dp) :: dw, ene, prefactor, cond_unit, tmp, de, mdfde
   real(dp), allocatable :: sp_cum(:,:,:,:), vvprod(:,:,:), &
      kwgt(:), op_cond(:,:,:), sptmp(:), enk(:,:)
   type(electron_wann) :: elec
   
   if(ionode) write(stdout,'(2x,a)') '>start loading data.  '
   call init_electron_wann(epwan_fid, kc_dim, elec)
   ! the energy range for optical conductivity is [0, boltz_emax]
   if(boltz_emax < 0.0) call errore('calc_trans_cumulant_opcond','boltz_max should be > 0',1)
   !setup k-grid
   call boltz_grid_load(kg)
   call init_boltz_grid_sym(kg, elec, band_min, band_max)
   numb = kg%numb;   numk = kg%nk_irr
   
   !generate energy grid for spectral function
   if(ionode)write(*,'(3E15.5)') cum_inner_emin,cum_inner_emax,cum_de
   call init_cum_wgrid(cum_outer_emin, cum_outer_emax, &
      cum_inner_emin, cum_inner_emax, cum_de, cum_outer_np, wg, spectral_np)
   dw = wg%sp_de
   spf_lw =   floor(spectral_emin / dw)
   spf_up = ceiling(spectral_emax / dw)

   !load spectral data
   allocate( sp_cum(spf_lw:spf_up, ntemper, numb, numk), enk(numb, numk) )
   call load_spectral(trim(prefix)//'_spectral_cumulant.h5', enk, sp_cum)
      
   if(ionode) write(stdout,'(2x,a)') '>start computing transport.  '
   !energy grid for optical conductivty, the energy step is the same as that for cumulant 
   nestep = nint( boltz_emax/dw ) + 1
   allocate( vvprod(6, numb, numk), kwgt(numk), op_cond(6, nestep, ntemper) )
   !collect the equivalent kpoints of each irreducible k-point.
   call sum_equivalent_kpoint(kg, vvprod, kwgt)
   !TODO: check consistence between enk and kg%enk_irr
   deallocate(enk, kwgt)

   op_cond = 0.0_dp
   call mp_split_pools( numk * numb * ntemper, jst, jend)
!$omp parallel default(shared) private(j, ik, ib, it, i, cut_lw, ie, tmp, ene, mdfde, sptmp)
   allocate(sptmp(spf_lw:spf_up))
!$omp do schedule(guided)
   do j = jst, jend
      ! j = (ik-1)*(numb*ntemper) + (ib-1)*ntemper + (it-1) + 1
      ik = (j-1) / (numb*ntemper) + 1
      ib = mod( (j-1)/ntemper, numb) + 1
      it = mod( (j-1), ntemper) + 1

      !remove very small value and negative value in spectral function (< sp_cut)
      sptmp(spf_lw:spf_up) = 0.0E0_dp
      sptmp = merge(sp_cum(:,it,ib,ik), sptmp, sp_cum(:,it,ib,ik) > sp_cut)
      !remove long tail of spectral function with total spectral weight < weight_cut
      cut_lw = spectral_cut_tail(sptmp, spf_lw, dw, weight_cut)
      ene = kg%enk_irr(ib,ik) - efermi(it)
      !for E = 0
      tmp = 0.0_dp
      do i = cut_lw, spf_up
         !(-df/dw)
         mdfde = mfermi_deriv( temper(it), ene + real(i,dp)*dw )
         tmp = tmp + sptmp(i) * sptmp(i) * mdfde
      enddo
      do i = 1, 6
         mdfde = tmp * dw * vvprod(i, ib, ik)
!$omp atomic
         op_cond(i, 1, it) = op_cond(i, 1, it) + mdfde
      enddo
      
      !for E > 0
      do ie = 1, nestep-1
         ! E = ie*dw
         tmp = 0.0_dp
         do i = cut_lw, spf_up
            if( (ie+i) > spf_up ) cycle  !skip if w+E exceed boundary
            mdfde = (fermi(temper(it), ene + i*dw) - &
                     fermi(temper(it), ene + (i+ie)*dw)) / (ie*dw)
            tmp = tmp + sptmp(i) * sptmp(i+ie) * mdfde
         enddo
         do i = 1, 6
            mdfde = tmp * dw * vvprod(i, ib, ik)
!$omp atomic
            op_cond(i, ie+1, it) = op_cond(i, ie+1, it) + mdfde
         enddo
      enddo
   enddo
!$omp end do
   deallocate(sptmp)
!$omp end parallel
   !collect results from different pools
   call mp_sum(op_cond, inter_pool_comm)
   
   !add prefactors
   prefactor = pi * kg%kweight * merge(1.0_dp, 2.0_dp, spinor) / volume
   !convert cond to (ohm*m)^-1  (ohm = V/A = V*s/C; hbar = Ryd*t0)
   cond_unit = alat*alat*unitcharge*1.0E10_dp/(bohr2ang*ryd2ev*timeunit)
   op_cond = op_cond * prefactor * cond_unit

   !output the optical conductivity
   if(ionode) then
      iout = find_free_unit()
      open(iout, file=trim(prefix)//'.opcond',form='formatted',status='unknown')
      write(iout,'(a,/)') "# it   w(eV)   optical_cond (1/Ohm/m) xx xy yy xz yz zz"
      do it = 1, ntemper
         write(iout, '(a)') ' '
         do ie = 1, nestep
            write(iout,'(i3, 2x, 7(ES15.8,2x))') it, (ie-1)*dw*ryd2ev, (op_cond(i,ie,it), i=1,6)
         enddo
         write(iout, '(a)') ' '
      enddo
      close(iout)
   endif
   deallocate(sp_cum, vvprod, op_cond)
   return
end subroutine calc_trans_cumulant_opcond

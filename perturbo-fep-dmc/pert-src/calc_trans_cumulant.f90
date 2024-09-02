!===============================================================================
! Copyright (C) 2016-2020 Jin-Jian Zhou, Jinsoo Park, I-Te Lu, Marco Bernardi
!
! This program is distributed under the terms of the GNU General Public License.
! See the file `LICENSE' in the root directory of this distribution, or obtain
! a copy of the License at <https://www.gnu.org/licenses/gpl-3.0.txt>.
!
! Author: jjzhou <jjchou.comphy@gmail.com>
! Comment:
!    transport calculations based on spectral functions fron cumulant method.
!
! Maintenance:
!===============================================================================

subroutine calc_trans_cumulant()
   use pert_const, only: dp, pi
   use pert_data,  only: volume, epwan_fid, kc_dim, alat
   use pert_utils, only: mfermi_deriv, fermi
   use band_structure, only: electron_wann, init_electron_wann
   use boltz_grid, only: grid, init_boltz_grid_sym, boltz_grid_load, sum_equivalent_kpoint
   use qe_mpi_mod, only: ionode, stdout, inter_pool_comm, mp_split_pools, mp_sum
   use pert_param, only: ntemper, efermi, temper, band_min, band_max, prefix, &
      cum_inner_emin, cum_inner_emax, cum_outer_emin, cum_outer_emax, cum_de, &
      cum_outer_np, debug, spectral_emin, spectral_emax, spectral_np, spinor, &
      boltz_emax, boltz_emin, boltz_de, hole
   use cumulant_utility, only: cum_wgrid, init_cum_wgrid, spectral_cut_tail, &
      ouput_integrated_dos_tdf, load_spectral
   use pert_output, only: progressbar_init, set_progress_step
   use boltz_trans_output, only: output_trans_coef
   use boltz_trans_mod, only: extract_trans_coeff
   implicit none
   type(grid) :: kg
   type(cum_wgrid) :: wg
   real(dp), parameter :: sp_cut = 1.0E-16_dp  !unit: Ryd^-1;
   real(dp), parameter :: weight_cut = 5.0E-4_dp  !unit: Ryd^-1;
   integer  :: it, ib, ik, numb, numk, spf_lw, spf_up, i, j, &
      jst, jend, iout, nestep, cut_lw, ie, nn, icount, ncount, iter, step, nstep
   real(dp) :: dw, ene, prefactor, cond_unit, tmp, de, mdfde, cumsum, tmp_cs, tmp_kk, fac
   real(dp), allocatable :: sp_cum(:,:,:,:), vvprod(:,:,:), kwgt(:), &
      cond(:,:), dens(:), intdos(:,:), inttdf(:,:), sptmp(:), enk(:,:), &
      egrid(:), cond_seebeck(:,:), seebeck(:,:), kk_coef(:,:), therm_cond(:,:)
   type(electron_wann) :: elec
   
   if(ionode) write(stdout,'(2x,a)') '>start loading data.  '
   call init_electron_wann(epwan_fid, kc_dim, elec)

   !setup k-grid
   call boltz_grid_load(kg)
   call init_boltz_grid_sym(kg, elec, band_min, band_max)
   numb = kg%numb;   numk = kg%nk_irr
   
   !generate energy grid for spectral function
   call init_cum_wgrid(cum_outer_emin, cum_outer_emax, &
      cum_inner_emin, cum_inner_emax, cum_de, cum_outer_np, wg, spectral_np)
   dw = wg%sp_de
   spf_lw =   floor(spectral_emin / dw)
   spf_up = ceiling(spectral_emax / dw)

   !energy grid for integrated DOS and TDF
   nestep = nint( (boltz_emax-boltz_emin)/boltz_de ) + 3
   allocate( egrid(nestep) )
   do i = 1, nestep
      egrid(i) = boltz_emin + (i-2)*boltz_de
   enddo

   !load spectral data
   allocate( sp_cum(spf_lw:spf_up, ntemper, numb, numk), enk(numb, numk) )
   call load_spectral(trim(prefix)//'_spectral_cumulant.h5', enk, sp_cum)
   !TODO: check consistence between enk and kg%enk_irr
   deallocate(enk)
      
   if(ionode) write(stdout,'(2x,a)') '>start computing transport.  '
   !energy grid for optical conductivty, the energy step is the same as that for cumulant 
   allocate(cond(6,ntemper), dens(ntemper), cond_seebeck(6,ntemper), seebeck(6,ntemper))
   allocate( vvprod(6, numb, numk), kwgt(numk), kk_coef(6, ntemper), therm_cond(6,ntemper) )
   !collect the equivalent kpoints of each irreducible k-point.
   call sum_equivalent_kpoint(kg, vvprod, kwgt)
   
   if(debug) then
      allocate( inttdf(nestep,ntemper), intdos(nestep,ntemper) )
      intdos = 0.0_dp
      inttdf = 0.0_dp
   endif
   cond = 0.0_dp
   dens = 0.0_dp
   kk_coef = 0.0_dp
   cond_seebeck = 0.0_dp
   call mp_split_pools(numk * numb * ntemper, jst, jend, ncount)
   call set_progress_step(ncount, step, nstep)
   call progressbar_init('trans_cumulant:')
   
   icount = 0
!$omp parallel default(shared) private(j, ik, ib, it, i, cut_lw, ie, &
!$omp&  tmp, tmp_cs, tmp_kk, ene, mdfde, sptmp, cumsum, iter)
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

      !compute conductivity
      ! cond: Int (-df/dw) A(w)^2 * v*v dw
      ! cond_seebeck: Int (-df/dw) A(w)^2 * v*v * (w-\mu) dw / T
      ! kk_coef: Int (-df/dw) A(w)^2 * v*v * (w-\mu)^2 dw / T
      tmp = 0.0_dp
      tmp_cs = 0.0_dp
      tmp_kk = 0.0_dp
      do i = cut_lw, spf_up
         ! (-df/dw) A(w)^2
         mdfde = mfermi_deriv(temper(it), ene+real(i,dp)*dw) * sptmp(i)*sptmp(i)
         !
         tmp = tmp + mdfde
         tmp_cs = tmp_cs + mdfde * ( ene + real(i,dp)*dw )
         tmp_kk = tmp_kk + mdfde * ( ene + real(i,dp)*dw )**2
      enddo
      !
      do i = 1, 6
         !for cond
         mdfde = tmp * dw * vvprod(i, ib, ik)
!$omp atomic update
         cond(i,it) = cond(i,it) + mdfde

         !for seebeck
         mdfde = tmp_cs * dw * vvprod(i, ib, ik) / temper(it)
!$omp atomic update
         cond_seebeck(i,it) = cond_seebeck(i,it) + mdfde

         !for therm_cond
         mdfde = tmp_kk * dw * vvprod(i, ib, ik) / temper(it)
!$omp atomic update
         kk_coef(i,it) = kk_coef(i,it) + mdfde
      enddo
      
      !using Int f(w)*A_nk(w) dw
      !do i = cut_lw, spf_up
      !   !f(w)*A_nk(w); for hole: using (1-f(w)) instead. [Note: 1-f(E) = f(-E)]
      !   mdfde = (ene + real(i,dp)*dw) * merge(-1.0_dp, 1.0_dp, hole)
      !   tmp = tmp + sptmp(i) * fermi( temper(it), mdfde )
      !enddo

      !compute density
      tmp = 0.0_dp
      !using Int -df(w)/dw * N_nk(w) dw, N(w) is cumulative sums of A
      cumsum = 0.0E0_dp
      do i = cut_lw, spf_up
         mdfde = mfermi_deriv( temper(it), ene + real(i,dp)*dw )
         !cumulative sum of A(w):  cumsum(w) = Int_{-inf}^{w} A(w') dw'
         cumsum = cumsum + sptmp(i) * dw
         tmp = tmp + cumsum * mdfde
      enddo

      mdfde = tmp * dw * kwgt(ik)
!$omp atomic
      dens(it) = dens(it) + mdfde
   
      !compute integrated DOS and TDF if debug mode is on.
      if(debug) then
         !integrated DOS: N(E) = sum_{nk} int_{w<E} A_nk(w) * dw
         do i = 1, nestep
            nn = min( nint((egrid(i) - kg%enk_irr(ib,ik)) / dw), spf_up )
            if( nn < cut_lw ) cycle
            tmp = sum(sptmp(cut_lw:nn)) * dw * kwgt(ik)
!$omp atomic
            intdos(i,it) = intdos(i,it) + tmp
         enddo
         !integrated TDF: sum_{nk} v_a(nk)*v_b(nk)*{ Int_{w<E} (|A_nk(w)|^2) * dw }
         do i = cut_lw, spf_up
            ! A(w)*A(w)*dw
            sptmp(i) = sptmp(i) * sptmp(i) * dw
         enddo
         do i = 1, nestep
            nn = min( nint((egrid(i) - kg%enk_irr(ib,ik)) / dw), spf_up )
            if( nn < cut_lw ) cycle
            ! only compute xx component
            tmp = sum( sptmp(cut_lw:nn) ) * vvprod(1,ib,ik)
!$omp atomic
            inttdf(i,it) = inttdf(i,it) + tmp
         enddo
      endif

      !track progress
!$omp atomic update
      icount = icount + 1
      !output progress
      iter = icount
      if( mod(iter, step).eq.0  .or. iter.eq.ncount ) then
!$omp critical (write_progress_cum_trans)
         write(stdout,'(8x, f7.2, a1)') (100.0_dp*iter)/ncount, '%'
!$omp end critical (write_progress_cum_trans)
      endif
      !
   end do
!$omp end do
   deallocate( sptmp )
!$omp end parallel
   !collect results from all the pools
   call mp_sum(cond, inter_pool_comm)
   call mp_sum(dens, inter_pool_comm)
   call mp_sum(kk_coef, inter_pool_comm)
   call mp_sum(cond_seebeck, inter_pool_comm)
   if(debug) then
      call mp_sum(intdos, inter_pool_comm)
      call mp_sum(inttdf, inter_pool_comm)
      ! #.states / U.C.
      intdos = intdos * kg%kweight
      inttdf = inttdf * kg%kweight
      !output dos & tdf
      if(ionode) call ouput_integrated_dos_tdf &
         (trim(prefix)//'.int_dos_tdf', egrid, intdos, inttdf, temper, efermi)
      deallocate(intdos, inttdf)
   endif

   !include the scale factor (alat)^2 for band velocity
   fac = pi * kg%kweight * merge(1.0_dp, 2.0_dp, spinor) * alat * alat / volume
   !omitted the factor q^2, other than this factor, cond is in Rydberg atomic unit
   cond = cond * fac
   !omitted the factor q*K_B, other than this, cond_seebeck in Rydberg atomic unit
   cond_seebeck = cond_seebeck * fac
   kk_coef = kk_coef * fac
   !carrier density: #. carrier / bohr^3
   dens = dens * kg%kweight * merge(1.0_dp, 2.0_dp, spinor) / volume
   
   seebeck = 0.0_dp
   therm_cond = 0.0_dp
   do it = 1, ntemper
      ! seebeck = sigma^-1 * [sigma*seebeck]
      ! 1, K_B/q is omitted, q is electron charge and K_b is boltzmann constant 
      ! 2, the vlue of seebeck here is dimensionless.
      call extract_trans_coeff &
         (temper(it), cond(:,it), cond_seebeck(:,it), kk_coef(:,it), seebeck(:,it), therm_cond(:,it))
   enddo
   call output_trans_coef(temper, efermi, dens, cond, seebeck, therm_cond)

   deallocate(sp_cum, vvprod, dens, cond, egrid, kwgt, cond_seebeck, seebeck, kk_coef, therm_cond)
end subroutine calc_trans_cumulant

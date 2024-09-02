!===============================================================================
! Copyright (C) 2016-2020 Jin-Jian Zhou, Jinsoo Park, I-Te Lu, Marco Bernardi
!
! This program is distributed under the terms of the GNU General Public License.
! See the file `LICENSE' in the root directory of this distribution, or obtain 
! a copy of the License at <https://www.gnu.org/licenses/gpl-3.0.txt>.
!
! Author: jjzhou <jjchou.comphy@gmail.com>
! Comment:
!  compute the electron mean free path, requires prefix.imsigma as input.
!
! Maintenance:
!===============================================================================

subroutine calc_electron_mfp()
   use pert_const, only: dp, timeunit, ryd2ev, ryd2mev, bohr2ang
   use pert_data,  only: alat, bg, epwan_fid, kc_dim, num_wann
   use pert_output,only: load_imsigma
   use pert_utils, only: find_free_unit
   use vector_list,only: vlist, load_vector_list
   use pert_param, only: prefix, fklist, band_min, band_max, ntemper, temper, efermi
   use qe_mpi_mod, only: ionode, mp_split_pools, mp_sum, inter_pool_comm, npool, stdout
   use band_structure, only: electron_wann, init_electron_wann, solve_band_velocity
   implicit none
   type(vlist) :: kl
   logical :: success
   integer :: numb, ik, kst, kend, uout, it, ib
   real(dp) :: eig(num_wann), v(3,num_wann), vel_abs
   real(dp), allocatable :: enk(:,:), vel(:,:,:), imsigma(:,:,:), enk_t(:,:)
   !
   type(electron_wann) :: elec
   !
   call init_electron_wann(epwan_fid, kc_dim, elec)
   
   numb = band_max - band_min + 1
   call load_vector_list(fklist, kl)
   if(npool > kl%nvec) &
      call errore('calc_electron_mfp','too many pools (npool > #.k-points)',1)
   allocate(enk(numb, kl%nvec), vel(3, numb, kl%nvec), enk_t(numb, kl%nvec))

   enk = 0.0_dp;  vel = 0.0_dp
   call mp_split_pools(kl%nvec, kst, kend)
!$omp parallel do schedule(guided) default(shared) private(ik, eig, v)
   do ik = kst, kend
      call solve_band_velocity(elec, kl%vec(1:3,ik), v, eigvalue=eig)
      enk(:, ik) = eig( band_min:band_max )
      !add the factor 'alat', and convert from Ryd*Bohr/hbar to Ryd*angstrom/hbar
      ! timeunit t0=hbar/Ryd, so the unit is equal to Angstrom/t0
      vel(1:3,1:numb,ik) = v(1:3,band_min:band_max)*alat*bohr2ang
   enddo
!$omp end parallel do
   call mp_sum(enk, inter_pool_comm)
   call mp_sum(vel, inter_pool_comm)

   !load imsigma
   allocate(imsigma(numb, kl%nvec, ntemper))
   call load_imsigma(imsigma, success, enk_t)
   if(.not. success) &
      call errore('calc_electron_mfp','Fail to load imsigma',1)
   if(any(abs(enk - enk_t) > 1.0E-6_dp)) &
      call errore('calc_electron_mfp','inconsistent band energy in'// trim(prefix) //'.imsigma',1)
   if( any(imsigma < 1.0E-16_dp) ) &
      write(stdout,'(5x,a)') 'Warning (calc_electron_mfp): negative imsigma (or less than 1.0E-16)'
      
   !convert imsigma to relaxation time in unit of t0
   ! scattering rate Gamma = 2*Im(Se)/hbar, tau = 1 / Gamma
   imsigma = 0.5_dp / imsigma

   !output band velocity and mean free path MFP = v_nk * tau_nk
   call cryst_to_cart(kl%nvec, kl%vec, bg, 1) !kpt in cart. coord. (in the unit tpiba)
   enk = enk * ryd2ev

if(ionode) then
   uout = find_free_unit()
   open(unit=uout, file=trim(prefix)//'.vel', status='unknown', form='formatted')
   write(uout,'(22x, a)') repeat('#',54)
   write(uout,'(22x, "#                    Band velocity                   #")')
   write(uout,'(22x, a)') repeat('#',54)
   write(uout,'(a)') '#   ik  ibnd   E(ibnd)(eV)     k.coord. (cart. alat)      '&
                        //'          vel-dir                  |vel| (m/s) '
   do ik = 1, kl%nvec
   do ib = 1, numb
      vel_abs = sqrt(dot_product(vel(:,ib,ik), vel(:,ib,ik)))
      write(uout,'(i6, 2x, i4, 2x, f12.6, 2x, 3(f8.5,1x), 4x, 3(f8.5,1x), 1x, es23.16)') &
      ik, ib, enk(ib,ik), kl%vec(:,ik), vel(:,ib,ik)/vel_abs, vel_abs/timeunit*1.0E-10_dp
   enddo
   enddo
   close(uout)

   ! output results
   uout = find_free_unit()
   open(unit=uout, file=trim(prefix)//'.mfp', status='unknown', form='formatted')
   write(uout,'(8x, "#================================================================#")') 
   write(uout,'(8x, "#        Electron Mean Free Path (tau_nk * |v_nk|, in nm)        #")')
   write(uout,'(8x, "#================================================================#")') 
   write(uout, '(8x, a, 10x, a, i7,3x, a10, i4,3x, a6, i4)' ) &
      '#', 'NO.k: ', kl%nvec, 'NO.bands: ', numb, 'NO.T: ', ntemper


   do it = 1, ntemper
      write(uout, '(a)') '#########'
      write(uout, '(a,f10.5,a,f10.5,a)' ) '# Temperature(T)= ',temper(it)*ryd2mev, &
         ' meV;  Chem.Pot.(mu)= ', efermi(it)*ryd2ev, ' eV'
      write(uout,1001) 
      write(uout,'(a)') &
         '# it     ik  ibnd   E(ibnd)(eV)   Relaxation time(in fs)           MFP (in nm) '
      write(uout,1001) 

      do ik = 1, kl%nvec
      do ib = 1, numb
         vel_abs = sqrt(dot_product(vel(:,ib,ik), vel(:,ib,ik)))
         write(uout,'(i3, 2x, i6, 2x, i4, 2x, f12.6, 2x, ES23.16, 2x, ES23.16)') &
         it, ik, ib, enk(ib,ik), imsigma(ib,ik,it)*timeunit*1.0E15_dp, &
         vel_abs*imsigma(ib,ik,it)*0.1_dp
      enddo
      enddo
   enddo
   close(uout)
endif
   deallocate(enk, enk_t, vel, imsigma)
   return
1001 format('#--------------------------------------------------------------------------------')
end subroutine calc_electron_mfp

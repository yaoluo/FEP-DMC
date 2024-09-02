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

module boltz_trans_output
   use pert_const, only: dp, ryd2ev, bohr2ang, kelvin2ev, ryd2mev, unitcharge, timeunit
   use pert_data,  only: tpiba, bg, system_2d, alat, at, volume
   use pert_utils, only: find_free_unit, mfermi_deriv
   use pert_param, only: prefix, boltz_kdim
   use boltz_utils,only: num2kpt
   use boltz_grid, only: grid
   implicit none
   private
   public :: output_kgrid,   output_dos,   output_tdf, output_trans_coef
   public :: output_density, output_rates, output_mobility, output_ftemper
contains

!output kgrid
subroutine output_kgrid(kg)
   implicit none
   type(grid), intent(in) :: kg
   ! local variables
   integer :: ik, uout
   real(dp) :: xk(3), xkc(3)

   uout = find_free_unit()
   open(uout, file=trim(prefix)//'_fullgrid.kpt', status='unknown', form='formatted')
   write(uout, '(1x, i8, 2x, a)')  kg%nk, &
      'Col.1:index;  (2, 3, 4): kpt in crystal; (5, 6, 7): kpt in cartestian (2pi/a).'
   do ik = 1, kg%nk
      xk = num2kpt( kg%kpt(ik), boltz_kdim)
      xkc(:) = xk(:)
      call cryst_to_cart(1, xkc, bg, 1)
      write(uout,'(5x,i8,1x,3(f10.6,2x),2x,3(f10.6,2x))') ik, xk(1:3), xkc(1:3)
   enddo
   close(uout)
end subroutine output_kgrid

subroutine output_dos(energy, dos)
   implicit none
   real(dp), intent(in) :: energy(:), dos(:)
   ! local variables
   integer :: uout, i
   integer, external :: find_free_unit

   uout = find_free_unit()
   open(unit=uout, file=trim(prefix)//'.dos',status='unknown',form='formatted')
   write(uout,'(2x,"#   E (eV)     DOS (#.states/eV/u.c.)")')
   do i = 1, size(energy)
      write(uout,'(2x, E14.6, 3x, E20.10)') energy(i)*ryd2ev, dos(i)/ryd2ev
   enddo
   close(uout)
end subroutine output_dos

subroutine output_rates(kg, tempe, efermi, rates, append)
   implicit none
   type(grid), intent(in) :: kg
   logical, intent(in) :: append ! write in append mode
   real(dp), intent(in) :: tempe, efermi, rates(:,:) ! rates(kg%numb, kg%nk)
   ! local variables
   integer :: uout, ik, ib
   logical :: created

   uout = find_free_unit()
   inquire(file=trim(prefix)//'.rates', exist=created)
   if(append .and. created) then
      open(uout, file=trim(prefix)//'.rates',status='old',action='write',position='append')
   else
      open(uout, file=trim(prefix)//'.rates',status='unknown',form='formatted')
      write(uout,'(1x,"#    Scattering rates computed in transport mode        #")')
      write(uout,'(1x,"#    (\Gamma in meV and N.B. \Gamma = 2 ImSigma)        #")')
   endif
   write(uout, '(a)')
   write(uout,'(1x, a, f9.4, a, f10.6/)') &  
      '#  Temperature: ', tempe*ryd2ev/kelvin2ev, '  Chemical Potential: ', efermi*ryd2ev
   do ik = 1, kg%nk
      do ib = 1, kg%numb
         write(uout,'(1x, i10, 3x, i4.4, 3x, f12.6, 3x, E22.10)') &
            ik, ib, kg%enk(ib,ik)*ryd2ev, rates(ib,ik)*ryd2mev
      enddo
   enddo
   close(uout)
end subroutine output_rates

subroutine output_tdf(energy, tdf, tempe, efermi, append)
   implicit none
   real(dp), intent(in) :: energy(:), tdf(:,:), tempe, efermi   !tdf(6,:)
   logical, intent(in), optional :: append
   ! local variables
   logical :: created
   integer :: uout, i, j
   
   uout = find_free_unit()
   inquire(file=trim(prefix)//'.tdf', exist=created)
   !
   if( present(append) ) created = append .and. created
   if(created) then
      open(uout, file=trim(prefix)//'.tdf',status='old',action='write',position='append')
   else
      open(uout, file=trim(prefix)//'.tdf',status='replace',form='formatted')
      write(uout,'(1x,"#  E(eV)    (-df/dE) (a.u.)    TDF(E)_(xx xy yy xz yz zz) (a.u.)   #")')
   endif
   write(uout,'(/1x, a, f9.4, a, f10.6/)') &
      '# Temperature: ', tempe*ryd2ev/kelvin2ev, '  Chemical Potential: ', efermi*ryd2ev
   do i = 1, size(energy)
      write(uout,'(f12.6, 3x, ES23.16, 2x, 6(E14.6,1x))') energy(i)*ryd2ev, &
         mfermi_deriv(tempe, (energy(i)-efermi)),  (tdf(j,i), j=1,6)
   enddo
   close(uout)
end subroutine output_tdf

subroutine output_density(temper, efermi, density)
   implicit none
   real(dp), intent(in) :: temper(:), efermi(:), density(:)
   character(len=7) :: label
   integer :: uout, it
   real(dp) :: dens_unit

   label = merge("(cm^-2)", "(cm^-3)", system_2d)
   !convert from #./bohr to #./Ang^-3 (3D) or #./cm^2 (2D)
   dens_unit = (1.0E0_dp/bohr2ang) * 1.0E8_dp
   dens_unit = merge(dens_unit**2 * alat * at(3,3), dens_unit**3, system_2d)

   ! output carrer concentration
   uout = find_free_unit()
   open(uout, file=trim(prefix)//'.doping', status='unknown', form='formatted')
   write(uout,'(1x,a,a)') "#", repeat('=',79)
   write(uout,'(1x,a)') "# Temperature(K), Chemical Potential(eV), Carrier concentration" // label
   do it = 1, size(temper)
      write(uout,'(4x, f7.2, 7x, f15.10, 9x, E15.7)')  &
         temper(it)*ryd2ev/kelvin2ev, efermi(it)*ryd2ev, density(it)*dens_unit
   enddo
   close(uout)
end subroutine output_density


subroutine output_ftemper(temper, efermi, density, fname)
   implicit none
   real(dp), intent(in) :: temper(:), efermi(:), density(:)
   character(len=*), intent(in) :: fname
   integer :: uout, it
   real(dp) :: dens_unit
   !convert from #./bohr to #./Ang^-3 (3D) or #./cm^2 (2D)
   dens_unit = (1.0E0_dp/bohr2ang) * 1.0E8_dp
   dens_unit = merge(dens_unit**2 * alat * at(3,3), dens_unit**3, system_2d)
   ! output carrer concentration
   uout = find_free_unit()
   open(uout, file=trim(fname), status='unknown', form='formatted')
   write(uout,'(i5, 5x, a)')  size(temper), 'F'
   do it = 1, size(temper)
      write(uout,'(1x, f7.2, 7x, f15.10, 9x, E15.7)')  &
         temper(it)*ryd2ev/kelvin2ev, efermi(it)*ryd2ev, density(it)*dens_unit
   enddo
   close(uout)
end subroutine


subroutine output_mobility(tempe, efermi, dens, cond, niter_in)
   implicit none
   ! cond(6, max_step, ntempe)
   real(dp), intent(in) :: tempe(:), efermi(:), dens(:), cond(:,:,:) !cond(6,:,:)
   integer, intent(in), optional :: niter_in(:)
   integer :: uout, it, ia, i, max_step, ntmp, last
   integer, allocatable :: niter(:)
   real(dp) :: cond_unit, mobi_unit, dens_unit
   character(len=120) :: ctmp, label
   character(len=3) :: sub(6) = (/'_xx', '_xy', '_yy', '_xz', '_yz', '_zz'/)
         
   label = merge("n_c (cm^-2)", "n_c (cm^-3)", system_2d)
   !convert from #./bohr to #./Ang^-3 (3D) or #./cm^2 (2D)
   dens_unit = (1.0E0_dp/bohr2ang) * 1.0E8_dp
   dens_unit = merge(dens_unit**2 * alat * at(3,3), dens_unit**3, system_2d)

   ! unit conversion of cond. check below.
   ! \sigma = q^2*nspin/(N*V) sum_{nk} v_nk*v_nk*tau_nk* f *(1-f)/kT
   !  (N.B.: the factor alat^2 should have already included in cond at this point)
   !
   ! 1. q^2 is omitted in the calculation of cond(i), 
   ! 2. all the quantities are in Rydberg atomic unit.
   !--------------------------------------------------------------------------
   !cond_unit = q^2* (a0)^-3* (Ryd/(hbar*a0^-1))^2 * (hbar/Ryd) / Ryd
   !          = e^2 * (hbar*a0)^-1 = e^2 *(Ryd*t0*a0)^-1
   !convert cond to (ohm*m)^-1  (ohm = V/A = V*s/C; hbar = Ryd*t0)
   cond_unit = unitcharge*1.0E10_dp/(bohr2ang*ryd2ev*timeunit)
   !for 2D, we output conductance
   if(system_2d) cond_unit = (cond_unit * bohr2ang * 1.0E-10_dp) * alat * at(3,3)

   !convert from Rydberg atomic unit to cm^2/V/s
   !mobility mu = sigma / (n*q); (n:(a0)^-3); mu = sigma * a0^3 / q 
   !mobi_unit = q* (a0)^2/(Ryd*t0)
   mobi_unit = bohr2ang*bohr2ang/(timeunit*1.0E16_dp*ryd2ev)
   
   ntmp = size(tempe)
   max_step = size(cond, 2)
   !
   allocate( niter( ntmp ) )
   niter = 1
   if(present( niter_in )) then
      if(size(niter) .ne. ntmp) call errore('output_mobility', 'niter dimension mismatch', 1)
      niter(:) = niter_in(:)
   else
      do it = 1, ntmp
         do i = 1, max_step
            if(sum(abs(cond(:,i,it))) < 1.0E-16_dp) exit
            niter(it) = i
         enddo
      enddo
   endif
   !
   ! array chcek
   if(size(dens).ne.ntmp .or. size(cond,3).ne.ntmp) &
      call errore('output_mobility','array dimension mismatch',1)
   uout = find_free_unit()
   open(uout, file=trim(prefix)//'.cond', status='unknown', form='formatted')

   write(uout,'(/10x,"#==========================================================#" )')
if(system_2d) then
   write(uout,'( 10x,"#                   Conductance (1/Ohm)                    #" )')
   write(uout,'( 10x,"#----------------------( 2D system )-----------------------#"/)')
else
   write(uout,'( 10x,"#                  Conductivity (1/Ohm/m)                  #" )')
   write(uout,'( 10x,"#----------------------------------------------------------#"/)')
endif
   
   write(ctmp, '(6(5x,a8,2x))') ('sigma'//sub(ia), ia = 1, 6)
   write(uout, '(a8, 3x,a7, 3x,a11, 1x,a90)') '#  T (K)', 'E_f(eV)', trim(label), ctmp
   do it = 1, ntmp
      last = niter(it)
      write(uout,'(f8.2, 1x, f9.5, 1x, E13.5, 1x, 6(1x, E14.6))') &
         tempe(it)*ryd2ev/kelvin2ev, efermi(it)*ryd2ev, dens(it)*dens_unit, &
         (cond(ia, last, it)*cond_unit, ia=1,6)
      !for iterative approach
      if(max_step > 1) then
         write(uout,'( 10x,"#--------------------iterative process---------------------#" )')
         write(uout,'(2x,a6,2x,a90)') "#iter.", ctmp
         do i = 1, niter(it)
            write(uout,'(2x,a1,i4,3x,6(1x, E14.6))') '#', i, (cond(ia,i,it)*cond_unit, ia=1,6)
         enddo
         write(uout,'( 10x,"#----------------------------------------------------------#"/)')
      endif
   enddo

   write(uout,'(/a)') 
   write(uout,'(/10x,"#==========================================================#" )')
   write(uout,'( 10x,"#                    Mobility (cm^2/V/s)                   #" )')
   write(uout,'( 10x,"#--------------------(for semiconductor)-------------------#"/)')

   write(ctmp, '(6(6x,a5,4x))') ('mu'//sub(ia), ia = 1, 6)
   write(uout, '(a8, 3x,a7, 3x,a11, 1x,a90)') '#  T (K)', 'E_f(eV)', trim(label), ctmp
   do it = 1, ntmp
      last = niter(it)
      write(uout,'(f8.2, 1x, f9.5, 1x, E13.5, 1x, 6(1x, E14.6))') &
         tempe(it)*ryd2ev/kelvin2ev, efermi(it)*ryd2ev, dens(it)*dens_unit, &
         (cond(ia, last, it)/dens(it)*mobi_unit, ia=1,6)
   enddo
   deallocate(niter)
   close(uout)
end subroutine output_mobility


subroutine output_trans_coef(tempe, efermi, dens, cond, seebeck, t_cond)
   implicit none
   real(dp), intent(in) :: tempe(:), efermi(:), dens(:), cond(:,:), seebeck(:,:), t_cond(:,:)
   !
   character(len=120) :: ctmp, label
   character(len=3) :: sub(6) = (/'_xx', '_xy', '_yy', '_xz', '_yz', '_zz'/)
   integer  :: uout, it, ia, ntmp
   real(dp) :: dens_unit, seebeck_unit, cond_unit, mobi_unit, t_cond_unit
   ntmp = size(tempe)

   ! (k_b/q) = (- k_b/e) = - k_b*(K) / (e*K) = - (kelvin2ev * eV) / (e* K) = -kelvin2ev * (V/K)
   seebeck_unit = kelvin2ev * 1.0E3_dp  ! convert to mV/K
   
   label = merge("n_c (cm^-2)", "n_c (cm^-3)", system_2d)
   !convert from #./bohr to #./Ang^-3 (3D) or #./cm^2 (2D)
   dens_unit = (1.0E0_dp/bohr2ang) * 1.0E8_dp
   dens_unit = merge(dens_unit**2 * alat * at(3,3), dens_unit**3, system_2d)

   !see comments in boltz_trans_output.f90/output_mobility for more detail.
   !cond_unit = q^2* (a0)^-3* (Ryd/(hbar*a0^-1))^2 * (hbar/Ryd) / Ryd
   !          = e^2 * (hbar*a0)^-1 = e^2 *(Ryd*t0*a0)^-1
   !convert cond to (ohm*m)^-1  (ohm = V/A = V*s/C; hbar = Ryd*t0)
   cond_unit = unitcharge*1.0E10_dp/(bohr2ang*ryd2ev*timeunit)
   if(system_2d) cond_unit = (cond_unit * bohr2ang * 1.0E-10_dp) * alat * at(3,3)

   !thermal conductivity (1/T * Ryd^2 *(cond / e^2) )
   ! \kappa = (1/T)* Ryd^2 * (Ryd*t0*a0)^-1 = (k_b/Ryd) * Ryd^2 * (Ryd*t0*a0)^-1
   !        = k_b / (t0 * a0) = (k_b*K) / (t0 * a0 * K) = kelvin2ev * (e*V/t0) / a0 / K
   ! convert to W/meter/K, W = J/s = C*V / s
   t_cond_unit = kelvin2ev * (unitcharge/timeunit) * 1.0E10_dp / bohr2ang

   !convert from Rydberg atomic unit to cm^2/V/s
   !mobility mu = sigma / (n*q); (n:(a0)^-3); mu = sigma * a0^3 / q 
   !mobi_unit = q* (a0)^2/(Ryd*t0)
   mobi_unit = bohr2ang*bohr2ang/(timeunit*1.0E16_dp*ryd2ev)
   
   uout = find_free_unit()
   open(uout, file=trim(prefix)//'.trans_coef', status='unknown', form='formatted')

   write(uout,'(/10x,"#==========================================================#" )')
if(system_2d) then
   write(uout,'( 10x,"#                   Conductance (1/Ohm)                    #" )')
   write(uout,'( 10x,"#----------------------( 2D system )-----------------------#"/)')
else
   write(uout,'( 10x,"#                  Conductivity (1/Ohm/m)                  #" )')
   write(uout,'( 10x,"#----------------------------------------------------------#"/)')
endif

   write(ctmp, '(6(5x,a8,2x))') ('sigma'//sub(ia), ia = 1, 6)
   write(uout, '(a8, 3x,a7, 3x,a11, 1x,a90)') '#  T (K)', 'E_f(eV)', trim(label), ctmp
   do it = 1, ntmp
      write(uout,'(f8.2, 1x, f9.5, 1x, E13.5, 1x, 6(1x, E14.6))') &
         tempe(it)*ryd2ev/kelvin2ev, efermi(it)*ryd2ev, dens(it)*dens_unit, &
         (cond(ia,it)*cond_unit, ia=1,6)
   enddo

   write(uout,'(/a)') 
   write(uout,'(/10x,"#==========================================================#" )')
   write(uout,'( 10x,"#                    Mobility (cm^2/V/s)                   #" )')
   write(uout,'( 10x,"#--------------------(for semiconductor)-------------------#"/)')

   write(ctmp, '(6(6x,a5,4x))') ('mu'//sub(ia), ia = 1, 6)
   write(uout, '(a8, 3x,a7, 3x,a11, 1x,a90)') '#  T (K)', 'E_f(eV)', trim(label), ctmp
   do it = 1, ntmp
      write(uout,'(f8.2, 1x, f9.5, 1x, E13.5, 1x, 6(1x, E14.6))') &
         tempe(it)*ryd2ev/kelvin2ev, efermi(it)*ryd2ev, dens(it)*dens_unit, &
         (cond(ia,it)/dens(it)*mobi_unit, ia=1,6)
   enddo

   write(uout,'(/a)') 
   write(uout,'(/10x,"#==========================================================#" )')
   write(uout,'( 10x,"#                Seebeck coefficient (mV/K)                #" )')
   write(uout,'( 10x,"#----------------------------------------------------------#"/)')
   !
   write(ctmp, '(6(7x,a4,4x))') ('S'//sub(ia), ia = 1, 6)
   write(uout, '(a8, 3x,a7, 3x,a11, 1x,a90)') '#  T (K)', 'E_f(eV)', trim(label), ctmp
   do it = 1, ntmp
      write(uout,'(f8.2, 1x, f9.5, 1x, E13.5, 1x, 6(1x, E14.6))') &
         tempe(it)*ryd2ev/kelvin2ev, efermi(it)*ryd2ev, dens(it)*dens_unit, &
         (seebeck(ia,it)*seebeck_unit, ia=1,6)
   enddo
   

   write(uout,'(/a)') 
   write(uout,'(/10x,"#==========================================================#" )')
   write(uout,'( 10x,"#                Thermal conductivity (W/m/K)              #" )')
   write(uout,'( 10x,"#------------------(Electronic contribution)---------------#"/)')
   write(uout,'( 10x,"#------------!! Not well tested, use with caution !!-------#"/)')
   !
   write(ctmp, '(6(5x,a8,2x))') ('kappa'//sub(ia), ia = 1, 6)
   write(uout, '(a8, 3x,a7, 3x,a11, 1x,a90)') '#  T (K)', 'E_f(eV)', trim(label), ctmp
   do it = 1, ntmp
      write(uout,'(f8.2, 1x, f9.5, 1x, E13.5, 1x, 6(1x, E14.6))') &
         tempe(it)*ryd2ev/kelvin2ev, efermi(it)*ryd2ev, dens(it)*dens_unit, &
         (t_cond(ia,it)*t_cond_unit, ia=1,6)
   enddo
end subroutine output_trans_coef

end module boltz_trans_output

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
subroutine transport()
   use pert_const, only: dp, ryd2ev
   use pert_utils, only: converged
   use pert_output,only: load_imsigma
   use pert_data,  only: volume, epwan_fid, qc_dim, kc_dim, alat
   use boltz_grid, only: grid, init_boltz_grid, boltz_grid_load
   use qe_mpi_mod, only: ionode, stdout, mp_barrier, inter_pool_comm
   use pert_param, only: ntemper, temper, efermi, ftemper, hole, doping, debug, &
      band_min, band_max, boltz_de, boltz_emax, boltz_emin, spinor, boltz_kdim, &
      trans_thr, prefix, boltz_nstep, boltz_qdim, full_ite
   use boltz_trans_output, only: output_mobility, output_rates, output_tdf
   use boltz_trans_mod, only: trans_cond_calc, trans_density_calc
   use boltz_scatter_integral, only: rates_scat_int, trans_scat_int
   use boltz_scatter, only: boltz_scatter_setup
   use HDF5_utils
   !
   use band_structure, only: electron_wann, init_electron_wann
   use phonon_dispersion, only: lattice_ifc, init_lattice_ifc
   implicit none
   type(grid) :: kg
   logical :: read_rate
   integer(HID_T) :: file_id
   integer :: nestep, i, it, ik, ib, io_nstep, niter(ntemper), ncomp
   real(dp) :: emin, emax, ef, tmpr, enk
   real(dp), allocatable :: ene(:), rates(:,:,:), cond(:,:,:), dos(:), &
      mfd0(:,:,:), mfd1(:,:,:), mfd_tmp(:,:,:), imsgm(:,:,:), tdf(:,:)
   !
   type(lattice_ifc)   :: phon
   type(electron_wann) :: elec
   
   if(ntemper .eq. 0) call errore('tranpsort', &
      'ftemper is not specified, missing temperatures and chemical potential',1)
   !
   call init_electron_wann(epwan_fid, kc_dim, elec)
   call init_boltz_grid(kg, elec, band_min, band_max)
   
   ! setup up energy windows for transport calculations
   emin = max(boltz_emin, minval(kg%enk))
   emax = min(boltz_emax, maxval(kg%enk))
   !write(stdout,'(5x,a,2(1x,f12.6),/)') 'Energy window (eV):', emin*ryd2ev, emax*ryd2ev
   if(emin > emax) call errore('transport','illegal energy window',1)
   !
   call mp_barrier(inter_pool_comm)
   write(stdout,'(5x,a)') '>finish init grid.'
   
   !setup energy grid: emin - boltz_de : emax + boltz_de
   nestep = int( floor((emax-emin)/boltz_de) ) + 3
   allocate( ene(nestep), dos(nestep) )
   do i = 1, nestep
      ene(i) = emin + (i-2)*boltz_de
   enddo
   call trans_density_calc(kg, ene, temper, efermi, hole, spinor, volume, doping, dos=dos)
   
   ! setup scattering if iterative solution or computing rates on the fly (boltz_nstep = 1)
   if( boltz_nstep > 0 ) then
      call init_lattice_ifc(epwan_fid, qc_dim, phon)
      call boltz_scatter_setup(kg, elec, phon, boltz_qdim)
   endif
   
   allocate( rates(kg%numb, kg%nk, ntemper), imsgm(kg%numb, kg%nk_irr, ntemper) )
   !read imsigma files: imaginary part of the selfenergy in meV.
   call load_imsigma(imsgm, read_rate)
   if( read_rate .and. boltz_nstep .ne. 1 ) then
      write(stdout,'(6x,a)') "- Read scattering rate from file: " // trim(prefix) // ".imsigma"
      !
      do it = 1, ntemper
      do ik = 1, kg%nk
      do ib = 1, kg%numb
         !Scattering rates in Rydberg atomic unit
         rates(ib,ik,it) = 2.0_dp*imsgm(ib, kg%kpt2ir(ik), it)
      enddo; enddo; enddo
   elseif (.not. read_rate .and. boltz_nstep .eq. 0) then
      call errore('tranpsort',"Can not find "//trim(prefix)//".imsigma",1)
   else
      write(stdout,'(6x,a)') "- Compute scattering rates on the fly."
      do it = 1, ntemper
         call rates_scat_int(kg, temper(it), efermi(it), rates(:,:,it))
      enddo
      !output rates
      if(debug .and. ionode) then
         do it = 1, ntemper
            call output_rates(kg, temper(it), efermi(it), rates(:,:,it), it>1)
         enddo
      endif
   endif
   deallocate( imsgm )
   if(any(rates(:,:,:) < 1.0E-9_dp)) write(stdout,'(5x, a)') &
      "Warn (transport): scattering rates less than 10^-9 a.u."
   !
   write(stdout,'(5x,a)') '>finish init rates'

   !! now all data we needed are ready. start to do the actual work.
   ! allocate work space first.
   io_nstep = max(1, boltz_nstep)
   niter(:) = 1
   ncomp = merge(6, 3, full_ite)  !if full_ite is true, both mfd for E and T field are computed
   allocate( mfd0(ncomp, kg%numb, kg%nk), cond(ncomp*2, io_nstep, ntemper), tdf(ncomp*2, nestep))
   if(boltz_nstep > 0) allocate(mfd1(ncomp,kg%numb,kg%nk), mfd_tmp(ncomp,kg%numb,kg%nk))
   !open prefix_tdf.h5 files
   if(ionode) call hdf_open_file(file_id, trim(prefix)//"_tdf.h5", status='NEW')
   
   write(stdout,'(5x,a)') '>start transport'
   write(stdout,'(5x,a)') '>progress: (current T / total T)'
   cond = 0.0_dp
   do it = 1, ntemper
      ef = efermi(it);  tmpr = temper(it)
      !compute Fnk = tau*vnk
      do ik = 1, kg%nk
      do ib = 1, kg%numb
         ! F for uniform E-field:  \tau_nk * v_nk
         mfd0(1:3,ib,ik) = kg%vnk(:,ib,ik) / rates(ib, ik, it)
         ! F for T gradienet:  \tau_nk * v_nk * (enk - \mu)
         if(full_ite) mfd0(4:6, ib,ik) = mfd0(1:3, ib,ik) * (kg%enk(ib,ik) - ef)
      enddo; enddo
      !compute conductivity and mobility under RTA, 
      ! or starting point for iter solution
      call trans_cond_calc(kg, ene, tmpr, ef, spinor, volume, mfd0, cond(:,1,it), tdf)
      !add the omitted factor, now tdf is fully in rydberg atomic unit
      tdf = tdf * alat * alat
      if(ionode) call write_tdf_hdf5(file_id, it, 1, nestep, tdf)

      !for iterative solution
      if(boltz_nstep > 0)  mfd1(:,:,:) = mfd0(:,:,:)
      !
      !compute correction from iterative solution
      do i = 2, io_nstep
         niter(it) = i
         !computed the integration terms (without the factor of tau0)
         call trans_scat_int(kg, temper(it), efermi(it), mfd1, mfd_tmp)
         !compute the new mean free displacement(mfd):
         do ik = 1, kg%nk
         do ib = 1, kg%numb
            mfd1(:,ib,ik) = mfd0(:,ib,ik) + mfd_tmp(:,ib,ik)/rates(ib,ik,it)
         enddo; enddo
         !compute conductivity from updated mfd1
         call trans_cond_calc(kg, ene, tmpr, ef, spinor, volume, mfd1, cond(:,i,it), tdf)
         tdf = tdf * alat * alat   
         !save tdf to hdf5 file
         if(ionode) call write_tdf_hdf5(file_id, it, i, nestep, tdf)
         !exit if converged (wthin relative error within trans_thr)
         if(converged(cond(:,i-1,it), cond(:,i,it), trans_thr)) exit
      enddo
      !output the converged tdf in text format
      if(ionode) call output_tdf(ene, tdf, tmpr, ef, it>1)
      !output progress
      write(stdout, '(5x, i4, a, i4)') it, " /", ntemper
      flush(stdout)
   enddo
   write(stdout,'(5x, a)') '>finish transport.'
   
   if(ionode) then
      call write_basic_hdf5(file_id, hole, niter, nestep, ene, dos)
      call hdf_close_file(file_id)
      !add the omitted factor alat**2 to cond.
      cond = cond * alat * alat
      !output conductivity and mobility
      call output_mobility(temper, efermi, doping, cond, niter)
      !compute other transport coefficient
      call trans_postproc()
   endif
   deallocate(ene, rates, mfd0, cond, tdf, dos)
   if(allocated(mfd1))  deallocate(mfd1)
   if(allocated(mfd_tmp)) deallocate(mfd_tmp)
   !
   return

   contains

      subroutine write_tdf_hdf5(fid, tid, iter, nstep, tdf)
         use pert_const, only: dp
         use HDF5_utils
         implicit none
         integer(HID_T), intent(in) :: fid
         integer, intent(in) :: tid, iter, nstep
         real(dp), intent(in) :: tdf(:,:)
         !
         character(len=120) :: dset_name
         character(len=6), external :: int_to_char

         dset_name = 'tdf_t' // trim(int_to_char(it)) // '_i' // trim(int_to_char(iter))
         call hdf_write_dataset(fid, trim(dset_name), tdf(:,1:nstep))
      end subroutine write_tdf_hdf5

      subroutine write_basic_hdf5(fid, is_hole, nit, nstep, energy, dos_out)
         use pert_param, only: ntemper, temper, efermi
         use pert_const, only: dp
         use HDF5_utils
         implicit none
         integer(HID_T), intent(in) :: fid
         logical, intent(in) :: is_hole
         integer, intent(in) :: nit(:), nstep
         real(dp), intent(in) :: energy(:), dos_out(:)
         !
         integer :: hole_i

         call hdf_write_dataset(fid, 'temperature', temper(1:ntemper))
         call hdf_write_attribute(fid, 'temperature', 'nt', ntemper)
         
         call hdf_write_dataset(fid, 'efermi', efermi(1:ntemper))
         call hdf_write_dataset(fid, 'num_iter', nit(1:ntemper) )

         call hdf_write_dataset(fid, 'energy_grid', energy(1:nstep))
         call hdf_write_attribute(fid, 'energy_grid', 'ne', nstep)

         call hdf_write_dataset(fid, 'dos', dos_out(1:nstep))
         hole_i = merge(1, 0, is_hole)
         call hdf_write_attribute(fid, 'dos', 'hole', hole_i)
      end subroutine write_basic_hdf5
      !
end subroutine transport

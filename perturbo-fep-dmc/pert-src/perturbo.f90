!===============================================================================
! Copyright (C) 2016-2020 Jin-Jian Zhou, Jinsoo Park, I-Te Lu, Marco Bernardi
!
! This program is distributed under the terms of the GNU General Public License.
! See the file `LICENSE' in the root directory of this distribution, or obtain 
! a copy of the License at <https://www.gnu.org/licenses/gpl-3.0.txt>.
!
! Author: jjzhou <jjchou.comphy@gmail.com>
! Comment:
!  main entrance for the perturbo code.
!
! Maintenance:
!===============================================================================

program perturbo
   use qe_mpi_mod, only: mp_startup, mp_barrier, mp_global_end, &
         inter_pool_comm, stdout, ionode, ionode_id
   use pert_param,  only: init_input_param, calc_mode, prefix
   use pert_data,   only: epwan_fid
   use hdf5_utils
   implicit none
   logical :: has_file
   
   call mp_startup()
   call environment_setup("PERTURBO", .true.)
   !readin contral parameter
   call init_input_param()
   write(stdout,'(5x,a)') "calc_mode = '" // trim(calc_mode) // "'"
   write(stdout,'(5x,a)') repeat("-",50)
   
   !open hdf file and load data
   epwan_fid = 0  ! only meaningful to ionode
   if(ionode) then
      inquire(file=trim(prefix)//"_epwan.h5", exist=has_file)
      if(.not. has_file) call errore('perturbo','missing '// trim(prefix)// "_epwan.h5 !", 1)
      !
      call hdf_open_file(epwan_fid, trim(prefix)//"_epwan.h5", status='OLD', action='READ')
   endif
   call load_data_pert(epwan_fid)
   
   select case (calc_mode)
      case ('bands')
      !compute wannier interpolated band structure
         call calc_band_structure()
      case ('phdisp')
      !compute phonon dispersion
         call calc_phonon_spectra()
      case ('ephmat')
      ! compute |g| (both with and without 1/sqrt(w) factor)
         call calc_ephmat()
      case ('imsigma')
      !compute on-shell electron self-energy (only imaginary part): Im\Sigma(e_nk)
         call electron_imsigma()
      case ('meanfp')
      !compute electron mean free path (require prefix.imsigma as input)
         call calc_electron_mfp()
      case ('setup')
      !setup a uniform k-grid for BTE dynamics or transport
         call boltz_setup()
      case ('trans')
      !transport calculation, compute mobility using RTA or iterative
         call transport()
      case ('trans-pp')
      !compute more transport coefficients starting from the pre-calculated TDF
         if(ionode) call trans_postproc()
      case ('dynamics-run')
      !carrier dynamics, real-time dynamics by time-steping BTE. (experimental)
         call carrier_dynamics_run()
      case ('dynamics-pp')
      !postprocessing carrier dynamics simulation data. (experimental)
         call carrier_dynamics_postproc()
      case ('spectral-se')
         ! compute off-shell electron self-energy (only imaginary part): Im\Sigma(w)
         call calc_selfenergy()
      case ('spectral-re')
         ! compute real part of the on-shell electron self-enrgy: Re\Sigma(e_nk)
         call electron_resigma()
      case ('spectral-pp')
         ! for debugging, compute spectral function using retarded cumulant approach.
         call process_spectral()
      case ('spectral-cum')
         !compute electron spectral function using retarded cumulant approach
         call calc_spectral_cumulant()
      case ('spectral-trans')
         !transport calculation from spectral function using Kubo formula
         call calc_trans_cumulant()
      case ('spectral-trans-opcond')
         !compute optical conductivity
         call calc_trans_cumulant_opcond()
      case ('diagmc-Gt-Matrix')
         !single k green function
         call diagmc_gt_matrix()
      case ('diagmc-EZ')
         !polaron formation energy 
         call diagmc_EZ_matrix()
      case ('diagmc-JJ')
         !current-current correlation function 
         call diagmc_JJ()
      case ('svd-elph')
         ! SVD e-ph Hamiltonian 
         call cal_eph_svd()
      case ('tabulate-H')
         ! tabulate e-ph Hamiltonian to disk for fine K/Q grid 
         call tabulate_H()
      case default
         write(stdout,'(1x,a,/,5x,a)') "Error: illegal calc_mode!"
   end select

   call mp_barrier(inter_pool_comm)
   if(ionode) call hdf_close_file(epwan_fid)
   
   call environment_stop("PERTURBO")
   call mp_global_end()
end program perturbo

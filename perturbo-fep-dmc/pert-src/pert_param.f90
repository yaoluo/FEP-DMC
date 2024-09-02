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
module pert_param
   use pert_const, only: dp, ryd2mev, ryd2ev, kelvin2eV, bohr2ang, timeunit
   implicit none
   save
   character(len=80) :: prefix
   character(len=80) :: tmp_dir
   type, public :: temperature
      integer :: nt
      real(dp), pointer :: kt(:) => null()
      real(dp), pointer :: ef(:) => null()
   end type
     ! jjzhou
   character(len=80) :: calc_mode
   character(len=80) :: fklist
   character(len=80) :: fqlist
   character(len=80) :: ftemper !file store temperature, chemical potential, doping level
   !polar_split
   !'polar',   only polar part; 
   !'rmpol', whole - polar (remaider part).
   !others,  the whole, no split;
   character(len=80) :: polar_split

!   logical :: lpolar
   logical :: find_efermi
   logical :: spinor ! if spinor is .true., then spin degeneracy will is added in DOS.
   logical :: hole   !.false.: electron, .true.: hole
   logical :: full_ite  ! if .true. solve BTE with both E- and T-field iteratively.
   !logical :: tdep   ! .true.: phonon from TDEP
   logical :: use_mem
   logical :: load_scatter_eph

!   logical :: boltz_restart ! .true. : enable restart in carrier dynamics simulations
   logical :: debug   ! .true. : debug mode
   integer :: boltz_kdim(3)
   integer :: boltz_qdim(3)  ! dimension for q-grid
   ! max number of steps in carrier dynamics or iterative steps in transport calc. 
   ! if boltz_nstep = 0, RTA calculations, otherwise iterative solutions
   integer :: boltz_nstep
   ! band_min, band_max: define bands to be computed.
   integer :: band_min, band_max !, num_bands
   ! boltz_emin, boltz_emax, and boltz_kdim are energy window and kgrid for tranpsort calculation.
   ! boltz_de, in mev energy step in tetrahedron integration
   ! trans_thr,  threshold of iterative procedure.
   !    (converged if relative error of two consecutive steps is within trans_thr)
   real(dp) :: boltz_emin, boltz_emax, boltz_de, trans_thr
   ! parameter for the gaussian broadening of the delta function
   real(dp) :: delta_smear
   ! phonon energy cutoff, phonon with energy smaller than phfreq_cutoff will be excluded.
   real(dp) :: phfreq_cutoff

   integer :: ntemper
   real(dp), allocatable :: temper(:)
   real(dp), allocatable :: doping(:)
   real(dp), allocatable :: efermi(:)
   type(temperature) :: tempers

   !for carrier dynamics
   real(dp) :: boltz_efield(3)  ! External Electric Field in V/cm
   real(dp) :: time_step  !time step of carrier dynamcis,in femto-second
   integer  :: output_nstep  ! output results for every 'output_nstep'
   character(len=80) :: solver !options: 'euler', 'rk4'
   character(len=80) :: boltz_init_dist !options: 'restart', 'lorentz', 'fermi', 'gaussian'
   !used if boltz_init_dist is lorentz, fermi or gaussian:
   !lorentz : f(k) = 0.1 / (((Ek-e0)/smear)^2 + 1.0 )
   !gaussian: f(k) = 0.1 * exp( -((Ek-e0)/smear)^2 )
   !fermi   : f(k) = 1.0 / (exp((Ek-e0)/smear) + 1.0)
   real(dp) :: boltz_init_e0 !eV
   real(dp) :: boltz_init_smear !meV
   !character(len=80) :: boltz_init_fname !restart from the last step of 'fname'

   !for two phonon scattering rate calculation
   !real(dp) :: smear_2ph
   character(len=80) :: sampling !options: 'uniform', 'cauchy'
   real(dp) :: cauchy_scale ! scale parameter gamma for cauchy distribution
   integer :: nsamples ! number of samples to generate
   !and cumulant approach
   ! inner window: energy_step = cum_de;  outer window: energy_step = cum_de*cum_outer_np
   real(dp) :: cum_inner_emin, cum_inner_emax, cum_outer_emin, cum_outer_emax, cum_de
   ! output spectral energy_setp: cum_de / spectral_np
   real(dp) :: spectral_emin, spectral_emax
   integer  :: cum_outer_np, spectral_np
   !for cumulant trans
   !logical  :: cum_velocity  !if ture, turn on band velocity correction

   !eph_svd 
   CHARACTER(len=200) :: svd_dir 
   logical :: read_svd 
   logical :: apply_svd 
   logical :: read_formf 
   integer :: nsvd 
   integer :: nk_svd,nq_svd
   logical :: onRAM 
   !diagmc 
   real(dp) :: FakeRatio 
   logical :: sample_gt 
   integer :: ib_sample 
   logical  :: Holstein = .false. , spinHolstein = .false. 
   real(dp) :: alpha_Holstein = 1.0, flip_strength = 0.1, w0_Holstein = 0.5, gauss_frolich = 0.5 
   integer :: dim_Holstein  
   logical :: LatticeFrohlich = .false.
   logical :: zeroTMC=.false. !zero temperature diagMC
   real(dp) :: alpha_frohlich = 1.0 
   integer :: Ntau 
   real(dp) :: tauMin,tauMax 
   integer :: Nbin, Nsample_MC 
   integer :: lowest_order ! order below this is not counted in green function
   logical :: exact_gkq = .false. 
   logical :: exact_electron = .false.
   logical :: print_sign = .false.
   integer :: DMC_Method = 1   !0,1; 0 = use full gkq for sampling, 1 = using g_q for sampling 
   real(dp) :: g_offset = 0.0_dp
   !Htable 
   real(dp) :: max_Re = 0.2d0 !cutoff for Re
   logical :: read_H = .false.
   logical :: frohlich_only = .false. 
   logical :: drop_frohlich = .false.



contains
subroutine init_input_param()
   use io_files, only: check_tempdir
   use pert_utils, only: find_free_unit
   use qe_mpi_mod, only: meta_ionode, meta_ionode_id, world_comm, mp_bcast, stdout
   implicit none
   logical :: ltmp1, ltmp2
   integer :: iunit, i, ios
   character(len=120) :: ctmp, msg
   !
   CHARACTER(LEN=256), EXTERNAL :: trimcheck

   namelist /perturbo/ prefix, calc_mode, fklist, fqlist, ftemper, hole, full_ite, debug, boltz_kdim, boltz_qdim, &
      band_min, band_max, boltz_emax, boltz_emin, tmp_dir, use_mem, boltz_de, trans_thr, &
      delta_smear, phfreq_cutoff, boltz_nstep, polar_split, load_scatter_eph, solver, &
      boltz_efield, time_step, boltz_init_dist, boltz_init_smear, boltz_init_e0, output_nstep, &
      sampling, cauchy_scale, nsamples, cum_inner_emax, cum_inner_emin, cum_outer_emax, &
      cum_outer_emin, spectral_emin, spectral_emax, cum_de, spectral_np, cum_outer_np, & 
      svd_dir, read_svd, apply_svd, read_formf, nsvd,nk_svd,nq_svd,onRAM,FakeRatio,Holstein,w0_Holstein,alpha_Holstein, &
      LatticeFrohlich, alpha_frohlich,Ntau,tauMin,tauMax,Nbin,lowest_order,dim_Holstein,gauss_frolich,Nsample_MC, &
      spinHolstein,flip_strength,zeroTMC,read_H,frohlich_only,exact_gkq,print_sign,drop_frohlich,exact_electron,max_Re, &
      DMC_Method, g_offset,sample_gt,ib_sample
   
   !set default value
   prefix       = ''
   tmp_dir      = '.'
   calc_mode    = '' !default
   polar_split  = '' !default, on split
   fqlist       = ''
   fklist       = ''
   ftemper      = ''
   debug        = .false.
   !spinor      = .false. ! deprecated, will read from hdf5 file instead.
   !tdep         = .false.
   hole         = .false. ! default, electron carrier
   full_ite     = .false. !
   use_mem      = .true.  ! only used in boltz_scatter
   load_scatter_eph = .false.  ! read eph_g2 from files

   boltz_kdim   =  1      
   boltz_qdim   =  0
   band_min     =  1    !the default will be 1
   band_max     =  9999999  !a very big number, will be set to nband
   boltz_emin   = -9999.0_dp !eV, 
   boltz_emax   =  9999.0_dp !eV, 
   boltz_nstep  =  0
   boltz_de     =  1.0_dp   !meV, 
   trans_thr    =  2.E-3_dp !0.2% as default
   delta_smear  =  10.0_dp  !meV
   phfreq_cutoff= 1.0_dp    !meV
   !for dynamics
   boltz_efield = 0.0_dp !turn off E-field by default
   time_step    = 1.0_dp !fs
   solver       = 'rk4' !use rk4 by default
   boltz_init_dist = ''  !no default value, should be specified.
   boltz_init_e0   = -9999.0_dp !eV, should be specified
   boltz_init_smear  = 20.0_dp  !meV, 20 by default
   output_nstep  = 1
   
   sampling = 'uniform' ! other options 'cauchy'
   cauchy_scale = 0.05
   nsamples = 100000 ! 0.1 million by default
   !energy range for spectral calculations
   spectral_emin = -0.1_dp  !eV
   spectral_emax =  0.1_dp  !eV
   cum_inner_emin = -0.1_dp
   cum_inner_emax =  0.1_dp
   cum_outer_emin = -0.1_dp
   cum_outer_emax =  0.1_dp
   cum_de   =  0.1_dp  !meV
   cum_outer_np = 1 !
   spectral_np = 1  !

   !svd 
   svd_dir = './gkq'
   read_svd = .true. 
   apply_svd = .false. 
   read_formf = .true.
   nsvd = 10 
   nk_svd = 50; nq_svd = 50 
   onRAM = .true.
   !diagmc 
   sample_gt = .false. 
   ib_sample = 1 
   FakeRatio = 1e-2
   Holstein = .false.; dim_Holstein = 3  
   spinHolstein = .false.; flip_strength = 0.1d0 
   LatticeFrohlich = .false. 
   zeroTMC = .false.
   gauss_frolich = 1.d0
   alpha_frohlich = 1.d0
   Ntau = 1
   tauMin = 50 
   tauMax = 200
   Nbin = 500
   Nsample_MC = 1000 
   lowest_order = 1 !=1 for not excluding any diagrams, =2 excludes Fan-Migdal. 
   g_offset = 0
   DMC_Method = 0 
   !Htable 
   max_Re = 0.1d0 !default for Re = 0 only 
   read_H = .false.
   !readin parameters
   if(meta_ionode) then
      call input_from_file()

      read(5, perturbo, err=100, iostat=ios)
100   call errore('init_input_para','reading input namelist',abs(ios))
      
      tmp_dir = trimcheck(tmp_dir)
      !do some check
      !band_min and band_max should be no less than 1
      if(band_min > band_max .or. band_min < 1 .or. band_max < 1) then
         msg = "both band_min and band_max should > 0, and band_min > band_max"
         call errore('init_input_para', trim(msg), 1)
      endif
      
      fqlist     = adjustl(fqlist)
      fklist     = adjustl(fklist)
      ftemper   = adjustl(ftemper)
      calc_mode = adjustl(calc_mode)
      polar_split = adjustl(polar_split)

      if(any(boltz_kdim(1:3) < 1)) &
         call errore('init_input_para','illegal boltz_kdim',1)

      if( all(boltz_qdim(1:3) .eq. 0) ) then
         !by default, boltz_qdim = boltz_kdim
         boltz_qdim(:) = boltz_kdim(:)
      elseif( any(boltz_kdim(1:3) < 1) ) then
         call errore('init_input_para','boltz_qdim should all be positive!',1)
      elseif( any( mod(boltz_kdim(:), boltz_qdim(:)) .ne. 0 ) ) then
         call errore('init_input_para','boltz_qdim is incommensurate with boltz_kdim',1)
      endif

      if(boltz_emin>boltz_emax .or. boltz_de<1.0E-3_dp ) call errore &
         ('init_input_para','illegal boltz_emax, boltz_emin or boltz_de.',1)
      
      boltz_emin = boltz_emin/ryd2ev
      boltz_emax = boltz_emax/ryd2ev
      if(boltz_nstep < 0) boltz_nstep = 0
      ! from mev to Ryd
      boltz_de = boltz_de/ryd2mev
      phfreq_cutoff = phfreq_cutoff/ryd2mev 
      delta_smear   = delta_smear/ryd2mev
      if(trans_thr < 1.0E-16) &
         call errore('init_input_param','trans_thr is too small or negative')

      !open temperature file
      iunit = find_free_unit()
      find_efermi  = .false.
      if(trim(prefix)  .eq. '') prefix = 'pert'
      if(trim(ftemper) .eq. '') then
         ntemper = 0
      else
         open(iunit,file=trim(ftemper),form='formatted',err=101,iostat=ios)
101      call errore('init_input_para','opening file '//trim(ftemper),ios)
         read(iunit,*,iostat=ios) ntemper, ctmp
         if(ios .ne. 0) call errore('init_input_param', &
            'reading ntemper in file '//trim(ftemper), 1)
         ctmp = trim( adjustl(ctmp) )
         call lower_case( ctmp )
         if( ctmp(1:1) .eq. 't' ) then
            find_efermi = .true.
         else if( ctmp(1:1) .eq. 'f' ) then
            find_efermi = .false.
         else
            find_efermi = .false.
            write(stdout,'(5x, a)') &
               "Warning: illegal mode in " // trim(ftemper) // ". Set to 'F' "
         endif
         !
         if(ntemper < 1) call errore('init_input_param', &
            '#. temperatures < 1 in '// trim(ftemper), 1)
         !
         allocate(temper(ntemper), doping(ntemper), efermi(ntemper))
         ! temper in K;  efermi in eV;   doping in cm^-3
         temper = 0.0_dp; efermi = 0.0_dp; doping = 0.0_dp
         do i = 1, ntemper
            read(iunit,*,iostat=ios) temper(i), efermi(i), doping(i)
            if(ios .ne. 0) call errore('init_input_para', &
               'reading temper in file '//trim(ftemper), 1)
         enddo
         close(iunit)
         temper = temper*kelvin2eV/ryd2ev
         efermi = efermi / ryd2ev
         ! do the conversion in a later stage since we don't know it's 2D or 3D system.
         !for 3D from #./cm^3 to #./bohr^3, for 2D from #./cm^2 to #./bohr^2
         !doping = doping*1.0E-24_dp*(bohr2ang)**3 
      endif

      !for dynamics
      if(time_step < 0.0_dp) call errore('init_input_param','negative time step',1)
      !convert to Rydberg atomic unit
      boltz_init_e0 = boltz_init_e0 / ryd2ev
      boltz_init_smear = boltz_init_smear / ryd2mev
      time_step = time_step / (timeunit*1.0E15_dp)
      !from e*V/cm to Rydberg atomic unit (e*E)
      ! convert to Rydberg atomic unit: eE, eV/cm -> E_Rydberg / Bohr
      boltz_efield(1:3)  = boltz_efield(1:3)*bohr2ang/ryd2ev*1.0E-8_dp
      boltz_init_dist = trim(adjustl(boltz_init_dist))
      solver = trim(adjustl(solver))
      if(trim(solver) .ne. 'euler' .and. solver .ne. 'rk4') &
         call errore('init_input_param',"solver should be 'euler' or 'rk4'.", 1)

      sampling = adjustl(sampling)
      if(trim(sampling) .eq. '')  sampling = 'uniform'

      !for cumulant
      if(cum_inner_emin > 0.0_dp .or. cum_inner_emax < 0.0_dp) &
         call errore('init_input_param','(spectral_emin, spectral_emax) should enclose 0.0',1)
      if(cum_outer_emin > cum_inner_emin .or. cum_outer_emax < cum_inner_emax) &
         call errore('init_input_param','outer window should enclose inner window',1)
      spectral_emin = spectral_emin / ryd2ev
      spectral_emax = spectral_emax / ryd2ev
      cum_inner_emin =  cum_inner_emin / ryd2ev 
      cum_inner_emax =  cum_inner_emax / ryd2ev
      cum_outer_emin =  cum_outer_emin / ryd2ev
      cum_outer_emax =  cum_outer_emax / ryd2ev
      ! de in meV
      cum_de   = cum_de  / ryd2mev
      if(spectral_np < 1)  spectral_np = 1
      if(cum_outer_np < 1)  cum_outer_np = 1
   endif
   
   !broadcast 
   call mp_bcast(prefix, meta_ionode_id, world_comm)
   call mp_bcast(tmp_dir, meta_ionode_id, world_comm)
   call mp_bcast(calc_mode, meta_ionode_id, world_comm)
   call mp_bcast(fklist, meta_ionode_id, world_comm)
   call mp_bcast(fqlist, meta_ionode_id, world_comm)
   call mp_bcast(ftemper, meta_ionode_id, world_comm)
   call mp_bcast(polar_split, meta_ionode_id, world_comm)
   !
   call mp_bcast(find_efermi, meta_ionode_id, world_comm)
   call mp_bcast(hole, meta_ionode_id, world_comm)
   call mp_bcast(full_ite, meta_ionode_id, world_comm)
   call mp_bcast(debug, meta_ionode_id, world_comm)
   !call mp_bcast(spinor, meta_ionode_id, world_comm)
   !call mp_bcast(tdep, meta_ionode_id, world_comm)
   call mp_bcast(use_mem, meta_ionode_id, world_comm)
   call mp_bcast(load_scatter_eph, meta_ionode_id, world_comm)
   !
   call mp_bcast(boltz_kdim, meta_ionode_id, world_comm)
   call mp_bcast(boltz_qdim, meta_ionode_id, world_comm)
   call mp_bcast(boltz_nstep, meta_ionode_id, world_comm)
   call mp_bcast(band_min, meta_ionode_id, world_comm)
   call mp_bcast(band_max, meta_ionode_id, world_comm)
   !
   call mp_bcast(boltz_emax, meta_ionode_id, world_comm)
   call mp_bcast(boltz_emin, meta_ionode_id, world_comm)
   call mp_bcast(boltz_de, meta_ionode_id, world_comm)
   call mp_bcast(trans_thr, meta_ionode_id, world_comm)
   call mp_bcast(delta_smear, meta_ionode_id, world_comm)
   call mp_bcast(phfreq_cutoff, meta_ionode_id, world_comm)
   !
   call mp_bcast(ntemper, meta_ionode_id, world_comm)
   tempers%nt = ntemper
   if(ntemper > 0) then
      if(.not. allocated(temper)) allocate(temper(ntemper))
      if(.not. allocated(efermi)) allocate(efermi(ntemper))
      if(.not. allocated(doping)) allocate(doping(ntemper))
      call mp_bcast(temper, meta_ionode_id, world_comm)
      call mp_bcast(doping, meta_ionode_id, world_comm)
      call mp_bcast(efermi, meta_ionode_id, world_comm)

      if(associated(tempers%kt)) deallocate(tempers%kt)
      if(associated(tempers%ef)) deallocate(tempers%ef)
      allocate( tempers%kt(ntemper), tempers%ef(ntemper) )
      tempers%kt(:) = temper(:)
      tempers%ef(:) = efermi(:)
   endif
   call check_tempdir(tmp_dir, ltmp1, ltmp2)
   !
   !for dynamics
   call mp_bcast(boltz_efield, meta_ionode_id, world_comm)
   call mp_bcast(time_step, meta_ionode_id, world_comm)
   call mp_bcast(solver, meta_ionode_id, world_comm)
   call mp_bcast(boltz_init_dist, meta_ionode_id, world_comm)
   call mp_bcast(boltz_init_e0, meta_ionode_id, world_comm)
   call mp_bcast(boltz_init_smear, meta_ionode_id, world_comm)
   call mp_bcast(output_nstep, meta_ionode_id, world_comm)

   call mp_bcast(sampling,      meta_ionode_id, world_comm)
   call mp_bcast(cauchy_scale,  meta_ionode_id, world_comm)
   call mp_bcast(nsamples,      meta_ionode_id, world_comm)
   !for cumulant
   call mp_bcast(spectral_emin, meta_ionode_id, world_comm)
   call mp_bcast(spectral_emax, meta_ionode_id, world_comm)
   call mp_bcast(cum_inner_emin, meta_ionode_id, world_comm)
   call mp_bcast(cum_inner_emax, meta_ionode_id, world_comm)
   call mp_bcast(cum_outer_emin, meta_ionode_id, world_comm)
   call mp_bcast(cum_outer_emax, meta_ionode_id, world_comm)
   call mp_bcast(cum_de  , meta_ionode_id, world_comm)
   call mp_bcast(spectral_np,  meta_ionode_id, world_comm)
   call mp_bcast(cum_outer_np,  meta_ionode_id, world_comm)
   !for eph svd 
   call mp_bcast(svd_dir, meta_ionode_id, world_comm)
   call mp_bcast(read_svd, meta_ionode_id, world_comm)
   call mp_bcast(apply_svd, meta_ionode_id, world_comm)
   call mp_bcast(read_formf, meta_ionode_id, world_comm)
   call mp_bcast(nsvd, meta_ionode_id, world_comm)
   call mp_bcast(nk_svd, meta_ionode_id, world_comm)
   call mp_bcast(nq_svd, meta_ionode_id, world_comm)
   call mp_bcast(onRAM, meta_ionode_id, world_comm)
   call mp_bcast(read_H, meta_ionode_id, world_comm)
   !for diagmc 
   call mp_bcast(FakeRatio, meta_ionode_id, world_comm)
   call mp_bcast(Holstein, meta_ionode_id, world_comm)
   call mp_bcast(alpha_Holstein, meta_ionode_id, world_comm) 
   call mp_bcast(w0_Holstein, meta_ionode_id, world_comm)
   call mp_bcast(dim_Holstein, meta_ionode_id, world_comm)
   call mp_bcast(LatticeFrohlich, meta_ionode_id, world_comm)
   call mp_bcast(gauss_frolich, meta_ionode_id, world_comm)
   call mp_bcast(alpha_frohlich, meta_ionode_id, world_comm)
   call mp_bcast(Ntau, meta_ionode_id, world_comm)
   call mp_bcast(tauMin, meta_ionode_id, world_comm)
   call mp_bcast(tauMax, meta_ionode_id, world_comm)
   call mp_bcast(Nbin, meta_ionode_id, world_comm)
   call mp_bcast(lowest_order, meta_ionode_id, world_comm)
   call mp_bcast(Nsample_MC, meta_ionode_id, world_comm)
   call mp_bcast(spinHolstein, meta_ionode_id, world_comm)
   call mp_bcast(flip_strength, meta_ionode_id, world_comm)
   call mp_bcast(frohlich_only, meta_ionode_id, world_comm)
   call mp_bcast(exact_gkq, meta_ionode_id, world_comm)
   call mp_bcast(print_sign, meta_ionode_id, world_comm)
   call mp_bcast(drop_frohlich, meta_ionode_id, world_comm)
   call mp_bcast(exact_electron, meta_ionode_id, world_comm)
   call mp_bcast(max_Re, meta_ionode_id, world_comm)
   call mp_bcast(DMC_Method, meta_ionode_id, world_comm)
   call mp_bcast(g_offset, meta_ionode_id, world_comm)
end subroutine init_input_param


!convert string to lower case
subroutine lower_case(string)
   character(len=*) :: string
   integer  :: i, ic, nlen

   nlen = len(string)
   do i = 1, nlen
      ic = ichar( string(i:i) )
      if( ic >= 65 .and. ic < 90 ) string(i:i) = achar(ic+32)
   end do
end subroutine lower_case

end module pert_param

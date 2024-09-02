!===============================================================================
! Copyright (C) 2016-2020 Jin-Jian Zhou, Jinsoo Park, I-Te Lu, Marco Bernardi
!
! This program is distributed under the terms of the GNU General Public License.
! See the file `LICENSE' in the root directory of this distribution, or obtain 
! a copy of the License at <https://www.gnu.org/licenses/gpl-3.0.txt>.
!
! Author: jjzhou <jjchou.comphy@gmail.com>
! Comment:
!  read parameters from the input file
!
! Maintenance:
!===============================================================================

module input_param
   use kinds, only: dp
   use constants, only: angstrom_au, RYTOEV
   use io_files,  only: prefix, tmp_dir
   use io_global, only: stdout, stdin, ionode, ionode_id
   use mp_world, only: world_comm, root
   use mp, only: mp_bcast, mp_sum
   implicit none
   public
   !dimension of the regular k-grid.
   integer, save :: kdim(3)
   character(len=256), save :: phdir
   !
   ! the following parameters are optional.
   ! num_band: number of bands in the selected range, [dft_band_min, dft_band_max]
   integer, save :: dft_band_min, dft_band_max
   !'dis_win_min' set in wannier90
   real(dp), save :: dis_win_min
   ! num_wann: number of wannier functions, num_band = dft_band_max - dft_band_min + 1
   integer, save :: num_wann, num_band
   ! using external prefix.eig
   character(len=256), save :: eig_corr
   ! Ewald summation parameter for polar correction of e-ph matrix elements.
   ! N.B. this is for e-ph polar correction only,
   !  to be consistent with rigid_bulk, enforce alpha=1.0 for lattice dynamical.
   real(dp), save :: polar_alpha
   !
   ! if debug is true, then stop after g(k,q) (does not compute g(Re, Rp))
   logical, save :: debug
   ! whether or not to perform wannier transformation when computing g(k,q)
   logical, save :: lwannier
   ! use pre-computed g(k,q)
   logical, save :: load_ephmat
   !
   ! asr           (character) indicates the type of Acoustic Sum Rule imposed
   !               - 'no': no Acoustic Sum Rules imposed
   !               - 'simple':  previous implementation of the asr used
   !                  (3 translational asr imposed by correction of
   !                  the diagonal elements of the force constants matrix)
   !               - 'crystal': 3 translational asr imposed by optimized
   !                  correction of the force constants (projection).
   character(len=10), save :: asr
   ! system_2d = .true. if it is 2d material, assume at(:,3) is direction that 
   ! breaks periodicity,and at(1:2, 3) should be 0, at(:, 1:2) should be in-plane
   logical, save :: system_2d
   ! thickness of the 2d system, only meanful if the system is 2D
   real(dp), save :: thickness_2d
   !
   !for LSDA spin polarized calculation, 'up', 'down', default 'none'
   character(len=4), save :: spin_component

   !for TDEP support
   logical, save :: tdep

   public :: read_input_param
contains

subroutine read_input_param()
! read input from stdin
   implicit none
   ! local variable
   logical :: has_file
   CHARACTER(LEN=256) :: outdir, fname
   integer :: ios, nk1, nk2, nk3

   CHARACTER(LEN=256), EXTERNAL :: trimcheck

   namelist / qe2pert / prefix, outdir, phdir, dft_band_min, dft_band_max, num_wann, system_2d, &
      nk1, nk2, nk3, debug, lwannier, load_ephmat, eig_corr, polar_alpha, asr, thickness_2d, dis_win_min, &
      spin_component, tdep
   
   !default value
   debug = .false.
   lwannier = .true.
   load_ephmat = .false.
   system_2d = .false.
   asr = 'crystal'
   !
   eig_corr = ''
   polar_alpha = 1.0_dp
   dft_band_min = 1
   dft_band_max = 10000  ! set to a large number
   num_wann = 1
   dis_win_min = -9999.0  !eV
   !
   thickness_2d = 6.0d0  !in Angstrom units
   !
   spin_component = 'none' ! by default: on spin or non-collinear,  'up': LSDA-up channel, 'dn'
   tdep = .false.
   !
   ! master read
   if (ionode) then
      ! read input file
      call input_from_file()
      read(stdin, nml=qe2pert, err=100, iostat=ios)

      kdim(1) = nk1
      kdim(2) = nk2
      kdim(3) = nk3
   end if   

   ! broadcast
   call mp_bcast(prefix, ionode_id, world_comm)
   call mp_bcast(outdir, ionode_id, world_comm)
   call mp_bcast(phdir, ionode_id, world_comm)
   
   call mp_bcast(kdim, ionode_id, world_comm)
   call mp_bcast(dft_band_min, ionode_id, world_comm)
   call mp_bcast(dft_band_max, ionode_id, world_comm)
   call mp_bcast(dis_win_min, ionode_id, world_comm)
   call mp_bcast(num_wann, ionode_id, world_comm)
   call mp_bcast(debug,  ionode_id, world_comm)
   call mp_bcast(lwannier, ionode_id, world_comm)
   call mp_bcast(load_ephmat, ionode_id, world_comm)
   call mp_bcast(system_2d, ionode_id, world_comm)
   call mp_bcast(asr, ionode_id, world_comm)
   call mp_bcast(thickness_2d, ionode_id, world_comm)
   call mp_bcast(spin_component, ionode_id, world_comm)
   call mp_bcast(tdep, ionode_id, world_comm)

   call mp_bcast(eig_corr, ionode_id, world_comm)
   call mp_bcast(polar_alpha, ionode_id, world_comm)

   tmp_dir = trimcheck(outdir)
   phdir = trimcheck(phdir)
   eig_corr = trim(adjustl(eig_corr))
   ! covert to rydberg atomic unit
   thickness_2d = thickness_2d * angstrom_au
   dis_win_min = dis_win_min / RYTOEV
   
   !sanity check
   if( polar_alpha < 1.0E-12 ) &
      call errore('input_param','polar_alpha too small !!!',1)
   
   if( system_2d .and. thickness_2d <= 0.0d0 ) then
      call errore('input_param','thickness_2d is negative (or too small) !!!',1)
   elseif (.not. system_2d) then
      thickness_2d = -6.0_dp  ! negative value means 3D system
   endif

   if( dft_band_min < 1 ) &
      call errore('input_param','dft_band_min shuld be larger than 0',1)
   
   if( load_ephmat ) then
      fname = trim(tmp_dir) // trim(prefix) // "_elph.h5"
      inquire(file=trim(fname), exist=has_file)
      if(.not. has_file) call errore('input_param', "Missing file: "//trim(fname), 1)
   endif
   
   if((asr.ne.'simple').and.(asr.ne.'crystal').and.(asr.ne.'no')) then
      call errore('input_param','invalid Acoustic Sum Rule :' // asr, 1)
   endif

#if defined(__TDEP)
   if(tdep .and. system_2d) write(stdout,'(2x,a)') &
      "Warning: TDEP for 2D polar materials is not tested, use with caution!!"
#else
   if(tdep) call errore('input_param',"tdep=.true. requires Perturbo compiled with -D__TDEP",1)
#endif

   return
   100 call errore('init_input_parameters','reading input namelist',abs(ios))
end subroutine read_input_param

end module input_param

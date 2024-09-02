!===============================================================================
! Copyright (C) 2016-2020 Jin-Jian Zhou, Jinsoo Park, I-Te Lu, Marco Bernardi
!
! This program is distributed under the terms of the GNU General Public License.
! See the file `LICENSE' in the root directory of this distribution, or obtain 
! a copy of the License at <https://www.gnu.org/licenses/gpl-3.0.txt>.
!
! Author: ilu <louis161789@gmail.com>; jjzhou <jjchou.comphy@gmail.com>; 
! Comment:
!  
! Maintenance:
!===============================================================================

module electronic_data
   use kinds, only: dp
   use io_files, only: prefix, iunpun, tmp_dir, postfix
   use input_param, only: dft_band_min, dft_band_max, kdim, num_wann, debug, &
      num_band, lwannier, eig_corr, load_ephmat, system_2d, dis_win_min, spin_component
   use pwcom, only: nkstot, nbnd, npwx, xk, ngk, igk_k, nks, et
   use gvect, only: ngm, g
   use gvecw, only: gcutw, ecutwfc
   use uspp,  only: nkb
   use cell_base, only: at, tpiba, bg
   use becmod, only : allocate_bec_type, bec_type, calbec, deallocate_bec_type
   use noncollin_module, only: noncolin, npol, nspin_mag
   !
   use mp, only: mp_bcast, mp_sum
   use mp_pools, only: inter_pool_comm
   use io_global, only: stdout, ionode, ionode_id
   implicit none
   public
   
   ! total number of kpoints:
   !   num_kpts = nkstot in case of nspin=1 (not lsda); 
   !   num_kpts = nkstot/2 in case of nspin = 2 (lsda)
   integer, save :: num_kpts
   ! for LSDA only,  'up' -> spin_updn = 1;  'down' -> spin_updn = 2
   ! for spinless or non-collinear case:  spin_updn = 1
   integer, save :: spin_updn
   ! coordinates of the k-points, xk_tot(3, num_kpts), in crystal coordinate
   ! NOTE: xk in pwcom is in cartesian coordinate (in the unit of 2piba)
   real(dp), allocatable, save :: xk_tot(:,:)
   
   ! nkg_tot, igk_k_tot: collected nkg and igk_k over all the pools (xk_tot)
   ! ktable: index repr. of kpt -> location in kpts, (hence ngk_tot, igk_k_tot)
   !  exmaple: ktable(  kpt2num(xk_tot(1:3,i), kdim)  ) = i
   integer, allocatable, save :: ngk_tot(:), igk_k_tot(:,:), ktable(:)

   ! subset of energy levels and wavefunctions, [bndmin, bndmax]
   real(dp), allocatable, save :: et_sub(:,:), et_corr(:,:)
   ! rot_wan(num_band, num_wann, num_kpts)
   ! unitary matrix that transfers evc from bloch to wannier gauge
   complex(dp), allocatable, save :: rot_wan(:,:,:)
   !
   integer, save :: nb_sub   ! nb_sub = merge(num_wann, num_band, lwannier)
   ! evc_sub(npol*npwx, nb_sub, num_kpts)
   complex(dp), allocatable, save :: evc_sub(:,:,:)
   ! wannier_center(3, num_wann) 
   ! Wannier function centers in cartesian coordinates (in the unit of alat)
   ! wannier_center_cryst(3, num_wann) in crystal coordinate
   real(dp), allocatable, save :: wannier_center(:,:), wannier_center_cryst(:,:)

   ! beta function in the pseudos and its derivative, for k in xk_tot
   ! NOTE: rot_wan does not apply to becp1%nc, 
   type(bec_type), allocatable, save :: becp1_tot(:)
   type(bec_type), allocatable, save :: alphap_tot(:,:)

   
   public :: init_electronic_data, kpt2num, num2kpt, get_kpoint_index, deallocate_evc_bec
contains

subroutine init_electronic_data()
   use lsda_mod, only: nspin, lsda
   !use buffers, only: open_buffer, get_buffer, close_buffer
   implicit none
   real(dp) :: xk_tmp(3)
   integer :: ik, ik_g, kidx
   integer, external :: global_kpoint_index
   
   !read pwscf data, init variables in pwcom, noncolin_module, ... 
   call pw_read_file()
   !
   call printout_info()
   !
   ! if LSDA, either spin_component = 'up'  or 'down' are required
   select case(spin_component)
   case ('up')
      !LSDA, k-points: 1 to nkstot/2 are for spin-up
      spin_updn = 1
      num_kpts = nkstot / 2
   case ('down')
      !LSDA, k-points: nkstot/2 + 1 to nkstot are for spin-dw
      spin_updn = 2
      num_kpts = nkstot / 2
   case default
      !spin-unpolorized and non-collinear,  kpts: 1 to nkstot
      spin_updn = 1
      num_kpts = nkstot
   end select
   !
   !constrain for 2D system
   if(system_2d .and. (any(abs(at(1:2,3))>1.0E-6_dp) .or. any(abs(at(3,1:2))>1.0E-6_dp))) &
      call errore('perturbo','2D system requires the third lattice vector along z-dirction',1)
   !
   if(dft_band_max > nbnd) &
      write(stdout,'(a,i5)') "Warning: dft_band_max exceed nbnd, reset to nbnd: ", nbnd
   dft_band_max = min(dft_band_max, nbnd)

   num_band = dft_band_max - dft_band_min + 1
   ! if lwannier is true: nb_sub = num_wann; else nb_sub = num_band
   nb_sub = merge(num_wann, num_band, lwannier)

   if(num_wann > num_band) call errore('init_electronic_data', &
      'number of WFs is larger than number of bands within [dft_band_min, dft_band_max].',1)
   if( num_kpts .ne. kdim(1)*kdim(2)*kdim(3) ) &
      call errore('init_electronic_data',"nk1, nk2, nk3 does not match with nscf calc.",1)
   
   !allocate space for global arrays
   allocate( et_sub(num_band, num_kpts) )
   allocate( xk_tot(3,num_kpts), ngk_tot(num_kpts), igk_k_tot(npwx, num_kpts) )
   et_sub = 0.0_dp;  xk_tot = 0.0_dp;  ngk_tot = 0;   igk_k_tot = 0;
   
   !collect kpoint info from all the pools
   do ik = 1, nks
      ik_g = global_kpoint_index(nkstot, ik)
      !
      !in case of LSDA and spin component is 'down', we need to subtract nkstot/2
      ik_g = get_kpoint_index(ik_g, num_kpts, spin_updn)
      !skip if not needed, such as spin-up states when spin component is 'down'
      if(ik_g < 0) cycle
      
      xk_tmp = xk(:, ik)
      !transfer xk from cartesian to crystal coordinate
      call cryst_to_cart(1, xk_tmp, at, -1)
      xk_tot(1:3, ik_g) = xk_tmp(1:3)

      ngk_tot(ik_g)  = ngk(ik)
      igk_k_tot(:, ik_g) = igk_k(:, ik)

      !init selected band energies
      et_sub(:, ik_g) = et(dft_band_min:dft_band_max, ik)
   enddo
   call mp_sum(xk_tot, inter_pool_comm)
   call mp_sum(ngk_tot, inter_pool_comm)
   call mp_sum(igk_k_tot, inter_pool_comm)
   call mp_sum(et_sub, inter_pool_comm)
      
   !create mapping, k coordinate -> location in xk_tot
   ! ktable( kpt2num(xk_tot(1:3), kdim) ) -> location of xk_tot(1:3) in xk_tot
   allocate( ktable(0:num_kpts-1) )
   ktable = 0
   do ik = 1, num_kpts
     kidx = kpt2num(xk_tot(:,ik), kdim)
     if(kidx < 0) call errore('init_electronic_data','kpoint is not in the k-grid',1)
     ktable(kidx) = ik
   enddo
   ! check if it's one-to-one mapping
   if(any(ktable < 1)) call errore('init_electronic_data',"missing kpoints in the k-grid",1)
   
   !skip rot_wan if in debug mode, and lwannier=.false. and load_ephmat=.false.
   if(.not. (debug .and. (.not. lwannier) .and. (.not. load_ephmat))) then
      !read rot_wan
      allocate( rot_wan(num_band, num_wann, num_kpts) )
      rot_wan = cmplx(0.0_dp, 0.0_dp, kind=dp)
      if(ionode) call read_u_matrix_wann(rot_wan, et_sub)
      call mp_bcast(rot_wan, ionode_id, inter_pool_comm)

      !read wannier_center
      allocate( wannier_center(3, num_wann ), wannier_center_cryst(3, num_wann) )
      wannier_center = 0.0_dp
      if(ionode) call read_xyz_wann(wannier_center)
      call mp_bcast(wannier_center, ionode_id, inter_pool_comm)
      !in crystal coordinate
      wannier_center_cryst(:,:) = wannier_center(:,:)
      call cryst_to_cart(num_wann, wannier_center_cryst, bg, -1)
   
      !readin eig_corr, if present
      if(trim(eig_corr) .ne. '') then
         allocate( et_corr(num_band, num_kpts) )
         if(ionode) write(stdout,'(a)') "read external eig file: "//trim(eig_corr)
         if(ionode) call read_eig_corr(eig_corr, num_band, num_kpts, et_corr)
         call mp_bcast(et_corr, ionode_id, inter_pool_comm)
      endif
   endif
   
   !if loading pre-computed ephmat, then no need to initialize wavefunctions
   if(.not. load_ephmat) call init_wavefunctions()
   return
end subroutine init_electronic_data

subroutine printout_info()
   use paw_variables,   only : okpaw
   use ldaU,            only : lda_plus_u
   use lsda_mod,        only : nspin, lsda
   use mp_bands,        only : nbgrp
   use uspp,            only : okvan
   use spin_orb,        only : lspinorb
   implicit none
   !
   write(stdout, '(/5x, "QE2PERT  ")') 
   write(stdout, '( 5x, "---------------------------")') 

   !PERTURBO does not support PAW fow now.
   if(okpaw)        call errore('qe2pert','PAW is not supported yet',1)
   if(nbgrp .ne. 1) call errore('qe2pert','band groups parallel (bgrp) not supported',1)
   
   if(lda_plus_u) then
      write(stdout, '(5x, "DFPT + U calculation")')
      if(noncolin) call errore('qe2pert', 'DFPT+U with non-collinear not supported',1)
   endif
   !
   select case(spin_component)
   case ('up')
      if(.not.lsda) call errore('qe2pert',"spin_component='up' only works with LSDA",1)
      write(stdout, '(5x, "Spin polarized: spin-up")')
   case ('down')
      if(.not.lsda) call errore('qe2pert',"spin_component='down' only works with LSDA",1)
      write(stdout, '(5x, "Spin polarized: spin-down")')
   case default
      if(lsda) call errore('qe2pert',"need spin_component='up' or 'down' for LSDA",1)
      if(noncolin) then
         if(lspinorb) then
            write(stdout, '(5x, "Non-collinear with SOC")')
         else
            write(stdout, '(5x, "Non-collinear without SOC")')
         endif
      else
         write(stdout, '(5x, "Spin unpolorized")')
      endif
   end select
   
   if(okvan) then
      write(stdout, '(5x, "USPP pseudopotentials")')
   else
      write(stdout, '(5x, "NC pseudopotentials")')
   endif
   write(stdout, '(5x, "---------------------------"/)') 
end subroutine printout_info

subroutine init_wavefunctions()
   USE pw_restart_new ! only: read_collected_wfc
   implicit none
   complex(dp), parameter :: cone  = (1.0_dp, 0.0_dp)
   complex(dp), parameter :: czero = (0.0_dp, 0.0_dp)
   !
   integer :: ik, ik_g, kidx, npw, ibnd, ig, ipol, npp
   complex(dp), allocatable :: evc(:,:), aux1(:,:), vkb(:,:)
   integer, external :: global_kpoint_index
   CHARACTER( LEN=256 )  :: dirname
   !
   dirname = TRIM( tmp_dir ) // TRIM( prefix ) // postfix
   !
   npp = npol*npwx
   !allocate space for wave-function
   allocate( evc_sub(npol*npwx, nb_sub, num_kpts) )
   !initialize
   evc_sub = cmplx(0.0_dp, 0.0_dp, kind=dp)

   !allocate space for all k, and initialize bec%nc, bec%r, bec%k with czero
   allocate( becp1_tot(num_kpts) )
   allocate( alphap_tot(3, num_kpts) )
   do ik = 1, num_kpts
      ! refer to allocate_phq.f90 for allocation
      call allocate_bec_type( nkb, num_band, becp1_tot(ik) )

      do ipol = 1, 3
         !refer to allcoate_phq.f90 for allocation
         call allocate_bec_type( nkb, num_band, alphap_tot(ipol, ik) )
      enddo
   enddo

   !
!$omp parallel default(shared) private(ik, ik_g, npw, vkb, aux1, evc, ipol, ibnd, ig)
   allocate( vkb(npwx, nkb), aux1(npp, num_band), evc(npp, nbnd) )
!$omp do schedule(guided)
   do ik = 1, nks
      ik_g = global_kpoint_index(nkstot, ik)
      !
      !in case of LSDA and spin component is 'down', we need to subtract nkstot/2
      ik_g = get_kpoint_index(ik_g, num_kpts, spin_updn)
      if(ik_g < 0) cycle
      
      npw = ngk_tot(ik_g)

      !readin Kohn-Sham wavefunctions
!$omp critical (read_wfc)
#if defined(__QE64)
      !work with QE 6.4.x
      call pw_read_collected_to_evc(dirname, iunpun, ik, ik_g, spin_updn, npw, evc)
#else
      !work with QE 6.5.x
      call read_collected_wfc(dirname, ik, evc)
#endif
!$omp end critical (read_wfc)
      
      !rotate to wannier gauge if required.
      if(lwannier) then
         !evc_sub(:,:) = matmul( evc_tmp(:,:), rot_wan(:,:) )
         call zgemm('n','n', npp, num_wann, num_band, cone, evc(1, dft_band_min), &
                  npp, rot_wan(:,:,ik_g), num_band, czero, evc_sub(:,:,ik_g), npp)
      else
         !in this case: nb_sub is equal to num_band
         evc_sub(:,:, ik_g) = evc(:, dft_band_min:dft_band_max)
      endif
      !
      !init beta function related quantities, adapted from PHonon/PH/phq_init.f90
      !
      ! NOTE: current_spin and isk only matter in lsda calc., in which the number 
      ! of k-points is doubled: the first half for spin up (current_spin=1), 
      ! the second half for spin down. (current_spin=2)
      ! for spinless and non-collinear calc. current_spin (=1) does not matter. 
      !
      !comment out, since neither init_us_2 nor calbec depends on current_spin
      !IF ( lsda ) current_spin = isk( ik )
      !
      ! ... d) The functions vkb(k+G)
      !
      call pw_init_us_2( npw, igk_k_tot(1,ik_g), xk(1,ik), vkb )
      !
      ! ... e) we compute the becp terms which are used in the rest of
      ! ...    the code
      !! calbec is defined in QE/Modules/becmod.f90
      call calbec( npw, vkb, evc(:, dft_band_min:dft_band_max), becp1_tot(ik_g) )
      !
      ! ... e') we compute the derivative of the becp term with respect to an
      !         atomic displacement
      !
      do ipol = 1, 3
       
         aux1=(0.d0,0.d0)
         do ibnd = 1, num_band
            do ig = 1, npw
               aux1(ig,ibnd) = evc(ig, dft_band_min+ibnd-1) * tpiba * ( 0.D0, 1.D0 ) * &
                               ( xk(ipol,ik) + g(ipol,igk_k_tot(ig,ik_g)) )
            end do
            if (noncolin) THEN
               do ig = 1, npw
                  aux1(ig+npwx,ibnd) = evc(ig+npwx, dft_band_min+ibnd-1) * tpiba*(0.D0,1.D0)*&
                            ( xk(ipol,ik) + g(ipol,igk_k_tot(ig,ik_g)) )
               end do
            end if
         end do 
         call calbec (npw, vkb, aux1, alphap_tot(ipol,ik_g) )
      end do 
   enddo
!$omp end do
   deallocate( vkb, aux1, evc)
!$omp end parallel
   
   !collect results
   do ik = 1, num_kpts 
      call mp_sum(evc_sub(:,:,ik), inter_pool_comm)
      call mp_sum_bec(becp1_tot(ik))
      do ipol = 1, 3
         call mp_sum_bec(alphap_tot(ipol, ik))
      enddo
   enddo
end subroutine init_wavefunctions


subroutine mp_sum_bec( bec )
   use mp, only: mp_sum
   use mp_pools, only: inter_pool_comm
   implicit none
   type(bec_type), intent(inout) :: bec

   if( allocated(bec%r)  ) call mp_sum(bec%r, inter_pool_comm)
   if( allocated(bec%nc) ) call mp_sum(bec%nc, inter_pool_comm)
   if( allocated(bec%k)  ) call mp_sum(bec%k, inter_pool_comm)
end subroutine mp_sum_bec


subroutine deallocate_evc_bec
   implicit none
   integer :: ik, ipol

   deallocate( evc_sub )
   
   do ik = 1, num_kpts
      call deallocate_bec_type( becp1_tot(ik) )

      do ipol = 1, 3
         call deallocate_bec_type( alphap_tot(ipol, ik) )
      enddo
   enddo
   
   deallocate( becp1_tot )
   deallocate( alphap_tot )
end subroutine deallocate_evc_bec


subroutine read_u_matrix_wann(u_mat_tot, eigval_opt)
   implicit none
   complex(dp), intent(out) :: u_mat_tot(num_band, num_wann, num_kpts)
   real(dp), intent(in) :: eigval_opt(num_band, num_kpts)
   !local variables
   character(len=256) :: fname
   integer :: iunit, ios, i, j, ik, nwan_tmp, numb_tmp, nkpts, kidx, imin
   real(dp), allocatable :: xk_wan(:,:)
   complex(dp), allocatable :: u_mat(:,:,:), u_dis_mat(:,:,:), u_dis_tmp(:,:)

   integer, external :: find_free_unit

   fname = trim(prefix)//'_u.mat'
   iunit = find_free_unit()
   open(iunit, file=trim(fname), status='old', action='read', err=100, iostat=ios)

   read(iunit, *)
   read(iunit, *) nkpts, nwan_tmp

   if(nkpts .ne. num_kpts) &
      call errore('read_u_matrix_wann','nkpts from w90 differs from that of pwscf',1)
   if(nwan_tmp .ne. num_wann) &
      call errore('read_u_matrix_wann','number of WFs in w90 differs from num_wann',1)
      
   allocate( xk_wan( 3, num_kpts ) )
   allocate( u_mat( num_wann, num_wann, num_kpts ) )
   xk_wan = 0.0E0_dp
   u_mat  = cmplx(0.0_dp, 0.0_dp, kind=dp)

   do ik = 1, num_kpts
      read(iunit, *)
      read(iunit, *) xk_wan(1:3, ik)
      read(iunit, "(f15.10,sp,f15.10)") ((u_mat(i,j,ik), i=1, num_wann), j=1, num_wann)
   enddo
   close(iunit)
   
   !for disentanglement
   if(num_band > num_wann) then
      fname = trim(prefix)//'_u_dis.mat'
      iunit = find_free_unit()
      open(iunit, file=trim(fname), status='old', action='read', err=101, iostat=ios)

      read(iunit,*)
      read(iunit,*) nkpts, nwan_tmp, numb_tmp

      if(numb_tmp .ne. num_band) &
         call errore('read_u_matrix_wann','nband in w90 differs from number of selected bands',1)

      allocate( u_dis_mat(num_band, num_wann, num_kpts), u_dis_tmp(num_band, num_wann) )
      u_dis_mat = cmplx(0.0_dp, 0.0_dp, kind=dp)

      do ik = 1, num_kpts
         ! skip the first two line
         u_dis_tmp = cmplx(0.0_dp, 0.0_dp, kind=dp)
         read(iunit,*)
         read(iunit,*) !xk_wan(1:3,ik)
         read(iunit,"(f15.10,sp,f15.10)") ((u_dis_tmp(i,j), i=1, num_band), j=1, num_wann)

         ! from W90/disentangle.f90/dis_window, define the u_dis.mat
         !! Note - in windows eigval_opt are shifted, so the lowest ones go
         !! from nfirstwin(nkp) to nfirstwin(nkp)+ndimwin(nkp)-1, and above
         !! they are set to zero.
         ! find the index of lowest band inside outer window at nkp-th
         imin = 0
         do i = 1, num_band
            if( (eigval_opt(i, ik) .ge. dis_win_min) ) then
               imin = i
               exit
            endif
         enddo
         if(imin < 1) call errore('read_u_matrix_wann', &
            'dis_win_min is too large, dis_win_min should be the same one used in W90!', 1)
         !
         u_dis_mat(imin:num_band, :, ik) = u_dis_tmp(1:(num_band-imin+1), :)
      end do
      close(iunit)
   endif
   
   u_mat_tot = cmplx(0.0_dp, 0.0_dp, kind=dp)
   do ik = 1, num_kpts
      kidx = kpt2num( xk_wan(:, ik), kdim )
      if(kidx < 0) call errore("read_u_matrix_wann","k-grid from w90 differ from nscf", 1)
      kidx = ktable( kidx )

      if( num_band > num_wann ) then
         u_mat_tot(:,:, kidx) = matmul( u_dis_mat(:,:,ik), u_mat(:,:,ik) )
      else
         u_mat_tot(:,:, kidx) = u_mat(:,:,ik)
      endif
   enddo

   deallocate( u_mat, xk_wan )
   if(allocated( u_dis_mat )) deallocate( u_dis_mat )
   if(allocated( u_dis_tmp )) deallocate( u_dis_tmp )
   return
100 call errore('read_u_matrix_wann', 'open '//trim(fname), abs(ios))
101 call errore('read_u_matrix_wann', 'open '//trim(fname), abs(ios))
end subroutine read_u_matrix_wann


subroutine read_xyz_wann(xyz_coord_wann)
   use ions_base, only : nat   
   use constants, only : bohr_radius_angs
   use cell_base, only : alat
 
   implicit none
   real(dp), intent(out) :: xyz_coord_wann(3, num_wann)
   !local variables
   character(len=256) :: fname
   character(len=1) :: dummy_character
   integer :: iunit, ios, num_rows, iw, i

   integer, external :: find_free_unit

   fname = trim(prefix)//'_centres.xyz'
   iunit = find_free_unit()
   open(iunit, file=trim(fname), status='old', action='read', err=100, iostat=ios)
   
   !read number of rows
   read(iunit, *) num_rows 
   if(num_rows .ne. (num_wann + nat)) &
      call errore('read_xyz_wann','number of rows from .xyz differs from num_wann + nat',1)
   
   !skip second line (comment line)
   read(iunit, *) 
   xyz_coord_wann = 0.0_dp
  
   !read the coordinates of the wannier centers and save to xyz_coord_wann
   do iw = 1, num_wann
      read(iunit, *) dummy_character, (xyz_coord_wann(i, iw), i=1, 3)
   enddo
   close(iunit)
  
   !convert the coordinates into atomic units
   xyz_coord_wann(:,:)=xyz_coord_wann(:,:)/bohr_radius_angs

   !convert the coordinates into alat units
   xyz_coord_wann(:,:)=xyz_coord_wann(:,:)/alat
   
   return
100 call errore('read_xyz_wann', 'open '//trim(fname), abs(ios))
end subroutine read_xyz_wann

!
subroutine read_eig_corr(file_eig, nb, numk, et)
   use constants, only: rytoev
   implicit none
   character(len=*), intent(in) :: file_eig
   integer, intent(in) :: nb, numk
   real(dp), intent(out) :: et(nb, numk)
   !local
   real(dp) :: rtmp
   integer :: iunit, ib, ik, itmp1, itmp2
   integer, external :: find_free_unit
   
   iunit = find_free_unit()
   open(iunit, file=file_eig, status='old',form='formatted',action='read')
   
   et = 0.E0_dp
   do ik = 1, numk
   do ib = 1, nb
      read(iunit, *) itmp1, itmp2, rtmp
      !sanity check
      if( (itmp1 .ne. ib) .or. (itmp2 .ne. ik) ) &
         call errore('read_eig_corr','format error in '//trim(file_eig),1)
      ! convert to rydberg atomic unit
      et(ib, ik) = rtmp / rytoev
   enddo; enddo
end subroutine read_eig_corr

! map num to (kx, ky, kz)
pure function num2kpt(num, nk) result(kpt)
   implicit none
   integer, intent(in) :: num
   integer, intent(in) :: nk(3)
   integer :: ipt(3), i
   real(dp) :: kpt(3)
   ! num should be non-negative integer, and 0 <= num <= nk1*nk2*nk3-1
   ! i=0, nk1-1; j=0, nk2-1; k=0, nk3-1;
   ! num =  k + j*nk3 + i*nk2*nk3; map num to (i,j,k)
   ! k = mod(num, nk(3)), and convert to real type
   ipt(3) = mod(num, nk(3))
   ! j = mod(num / nk(3), nk(2))
   ipt(2) = mod(num/nk(3), nk(2))
   ! i =     num / (nk(2)*nk(3))
   ipt(1) = num/(nk(2)*nk(3))

   kpt(:) = real(ipt(:), dp) / real(nk(:), dp)
   !fold to the Gamma-centered FBZ, in crystal coordinate
   do i = 1, 3
      ipt(i) = merge(ipt(i)-nk(i), ipt(i), nint(kpt(i)) > 0)
   enddo
   !avoiding subtraction to prevent numerical noise.
   kpt(:) = real(ipt(:), dp) / real(nk(:), dp)
end function num2kpt

! map (kx, ky, kz) to num, if num < 0, then (kx, ky, kz) is not in the list
pure function kpt2num(kpt, nk) result(num)
   implicit none
   real(dp), intent(in) :: kpt(3)
   integer, intent(in) :: nk(3)
   integer :: num
   real(dp), parameter :: eps = 1.0E-5_dp
   ! local variables
   integer :: r(3)  !, i 
   real(dp) :: xkr(3), dis(3)
   ! init num to the default value 
   num = -1
   ! fold to the Gamma-centered FBZ, in crystal coordinate
   ! and check if kpt(i)*nk(i) is a integer or not.
   xkr(:) = (kpt(:) - nint(kpt(:))) * nk(:)
   dis(:) =  xkr(:) - nint(xkr(:))
   ! return -1 if (kx, ky, kz) is not in the k-mesh.
   if( sqrt(dot_product(dis, dis)) > eps ) return
   ! ri = 0...nki-1; r(1)->i; r(2)->j; r(3)->k
   r(:) = mod( nint(xkr(:)+2*nk(:)), nk(:) )
   ! num = k + j*nk3 + i*nk2*nk3 + 1
   num = r(3) + r(2)*nk(3) + r(1)*nk(2)*nk(3)
end function kpt2num

! in case of LSDA, kpoints are doubled, with the second half of kpoints for spin-down
! use this function to get the real kpoint index, e.g. subtract nkstot/2 for spin-down
pure function get_kpoint_index(ik_g, nkpts, ispin) result(ik)
   implicit none
   ! ispin should be either 1 or 2, 2 for spin-down
   integer, intent(in) :: ik_g, nkpts, ispin
   integer :: ik
   !local
   integer :: ii
   
   !for ispin=1 and ik_g <= nkpts, (LSDA-up, unpolorized and non-collinear)
   !for ispin=2 and ik_g >  nkpts  (LSDA-dn) (ii=2 for both case)
   ! otherwise return a negative value (ii = 1 or 3)
   ii = ispin + 1 - (ik_g-1)/nkpts
   ik = ( mod(ik_g-1, nkpts) + 1 ) * (-1)**ii
end function get_kpoint_index

end module electronic_data

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

subroutine save_elph_mat_wann()
   use kinds, only: dp
   use ions_base, only: nat, tau
   use cell_base, only: at, bg, omega, tpiba
   use lattice_data, only: numq, xq_tot, qdim
   use perturbo_data, only: lpolar, zstar, epsil
   use input_param, only: num_wann, kdim, polar_alpha, prefix, thickness_2d
   use electronic_data, only: xk_tot, num_kpts, wannier_center_cryst
   use qe_mpi_mod, only: mp_sum, mp_barrier, mp_bcast, mp_scatter, my_pool_id, &
      inter_pool_comm, npool, stdout, ionode, ionode_id, mp_split_pools
   use epwan_hdf5_io, only: write_elph_mat_wann_part, read_ephmat_part
   use elph_matrix_wannier, only: elph_mat_wann, eph_wannier, set_ws_cell_eph
   use polar_correction, only: set_polar_parameter, eph_wan_longrange
   use hdf5_utils
   implicit none
   logical :: wan_gauge
   character(len=120) :: fname
   real(dp) :: cryst_tau(3, nat), kqfac
   integer(HID_T) :: elph_fid, epwan_fid, epwan_gid
   integer :: numb, eph_nword, iq_st, iq_end, nq_loc, ndeg_e, ndeg_q
   integer :: iq, i, ia, ic, iw, jw, ire, irp
   integer, allocatable :: sendcount(:), displs(:)
   complex(dp) :: ctmp(3*nat)
   complex(dp), allocatable :: epmat(:,:,:,:), epmat_loc(:,:,:,:,:), &
                               ep_tmp(:,:,:,:,:), eplr(:,:,:)
   !
   type(elph_mat_wann) :: ep
   type(eph_wannier), pointer :: pep
   !
   write(stdout, '(5x,a)') "computing electron-phonon matrix in wannier basis..."
   call mp_barrier( inter_pool_comm )

   ! initialize
   kqfac = 1.0E0_dp / real(num_kpts * numq, dp)
   ep%lpol = lpolar
   !transfer atomic position to crystal coordinate
   cryst_tau(:,:) = tau(:,:)
   call cryst_to_cart(nat, cryst_tau, bg, -1)
   !setup wigner seitz cell (does not allocate space for ep_hop)
   call set_ws_cell_eph &
      (ep, kdim, qdim, nat, num_wann, at, cryst_tau, wannier_center_cryst, .false.)
  
   !open files
   if(ionode) then
      !open epwan.h5 file
      fname = trim(prefix)//"_epwan.h5"
      call hdf_open_file(epwan_fid, trim(fname), status='OLD', action='WRITE')
      call hdf_create_group(epwan_fid, 'eph_matrix_wannier')
      call hdf_open_group(epwan_fid, "eph_matrix_wannier", epwan_gid)
   
      call open_elph_file(elph_fid, wan_gauge, numb)
   endif
   call mp_bcast(wan_gauge, ionode_id, inter_pool_comm)
   call mp_bcast(numb, ionode_id, inter_pool_comm)

   eph_nword = numb * numb * num_kpts
   !set up mpi parallelization
   call mp_split_pools(numq, iq_st, iq_end, nq_loc)
   !set sendcount and displs for mp_scatter
   allocate( sendcount(npool), displs(npool) )
   displs = 0
   sendcount = 0
   sendcount( my_pool_id+1 ) = nq_loc * eph_nword
   call mp_sum(sendcount, inter_pool_comm)
   do i = 2, npool
      displs(i) = sum( sendcount( 1:(i-1) ) )
   enddo

   !setup polar parameter
   if(ep%lpol) then
      call set_polar_parameter &
         (qdim, nat, omega, tpiba, bg, epsil, zstar, tau, ep%pol, thickness_2d, polar_alpha)
      !
      allocate( eplr(3, nat, nq_loc) )
      eplr = (0.0_dp, 0.0_dp)
!$omp parallel do schedule(guided) default(shared) private(i, iq, ctmp)
      do i = 1, nq_loc
         iq = iq_st + i - 1
         ! calculate long-range part directly in cartesian coordinate
         call eph_wan_longrange(ep%pol, xq_tot(:,iq), ctmp)
         eplr(:, :, i) = reshape(ctmp, (/3, nat/) )
      enddo
!$omp end parallel do
   endif
   
   !the main loop
   do ia = 1, nat
      !workspace
      allocate( epmat_loc(numb, numb, num_kpts, nq_loc, 3) )
      
      !read and distribute data
      if(ionode) then
         allocate( epmat(numb, numb, num_kpts, numq) )
      else
         allocate( epmat(1,1,1,1) )
      endif
      !load epmat data
      do ic = 1, 3
         if(ionode) call read_ephmat_part(elph_fid, ia, ic, numb, num_kpts, numq, epmat)
         call mp_scatter(epmat, epmat_loc(:,:,:,:, ic), &
                     sendcount, displs, ionode_id, inter_pool_comm)
      enddo
      deallocate( epmat )
      
      !re-organize elph data
      allocate( ep_tmp(3, num_kpts, nq_loc, num_wann, num_wann) )
      call rearrange_ephmat(numb, iq_st, nq_loc, wan_gauge, epmat_loc, ep_tmp)
      deallocate( epmat_loc )
      
      !subtract long-range part in the wannier gauge
      if(ep%lpol) then
         do ic = 1, 3
         do iw = 1, num_wann
         do iq = 1, nq_loc
            ep_tmp(ic, :, iq,iw,iw) = ep_tmp(ic, :, iq,iw,iw) - eplr(ic, ia, iq)
         enddo; enddo; enddo
      endif
      
      !start to the main task
      do jw = 1, num_wann
      do iw = 1, num_wann
         pep => ep%epwan(iw, jw, ia)
         
         !allocate space for ep_hop here (skipped in set_ws_cell_eph)
         allocate( pep%ep_hop(3, pep%ws_el%nr, pep%ws_ph%nr) )
         pep%ep_hop = cmplx(0.0_dp, 0.0_dp, kind=dp)
         
         !perform the Fourier transformation
         do i = 1, nq_loc
            iq = iq_st + i - 1
            call eph_bloch2wann(xq_tot(:,iq), xk_tot, num_kpts, ep_tmp(:,:,i,iw,jw), &
               ep%nrvec_e, ep%rvec_set_el, ep%nrvec_p, ep%rvec_set_ph, pep)
         enddo
         !collect result over different pools
         call mp_sum(pep%ep_hop, inter_pool_comm)
         
         !add the prefactor and weight
!$omp parallel do schedule(guided) default(shared) private(irp, ndeg_q, ire, ndeg_e)
         do irp = 1, pep%ws_ph%nr
            ndeg_q = pep%ws_ph%ndeg(irp)
            do ire = 1, pep%ws_el%nr
               ndeg_e = pep%ws_el%ndeg(ire)
               !include ndege  
               pep%ep_hop(:,ire,irp) = pep%ep_hop(:,ire,irp) / real(ndeg_e*ndeg_q, dp)
            enddo; 
            !add 1/(Nk * Nq) factor
            pep%ep_hop(:,:,irp) = pep%ep_hop(:,:,irp) * kqfac
         enddo
!$omp end parallel do
         
         !output to hdf5 file
         if(ionode) call write_elph_mat_wann_part(epwan_gid, iw, jw, ia, pep)
         !
         deallocate( pep%ep_hop )
      enddo; enddo
      !
      deallocate(ep_tmp)
   enddo

   !close hdf file
   if(ionode) then
      call hdf_close_group(epwan_gid)
      !
      call hdf_close_file(epwan_fid)
      call hdf_close_file(elph_fid)
   endif

   deallocate(sendcount, displs)
end subroutine save_elph_mat_wann


subroutine open_elph_file(file_id, wannier_gauge, nb)
   use hdf5_utils
   use lattice_data, only: numq
   use electronic_data, only: num_kpts
   use input_param, only: num_band, num_wann, tmp_dir, prefix
   implicit none
   integer, intent(out) :: nb
   integer(HID_T), intent(out) :: file_id
   logical, intent(out) :: wannier_gauge
   !
   integer :: wan_int, ndim(4)
   character(len=120) :: fname
   
   !open hdf5 file: _elph.h5
   fname = trim(tmp_dir) // trim(prefix) // "_elph.h5"
   call hdf_open_file(file_id, trim(fname), status='OLD', action='READ')
   call hdf_read_dataset(file_id, 'wannier_gauge', wan_int)
   wannier_gauge = merge(.true., .false., wan_int > 0)
   ndim = 0
   call hdf_get_dims(file_id, "elph_1_1_r", ndim)
   if(any(ndim < 1)) call errore("open_elph_file","wrong rank in " // trim(fname), 1)

   !sanity check
   if( num_kpts .ne. ndim(3) .or. numq .ne. ndim(4) ) &
      call errore("open_elph_file","mismatch dimension in elph data",1)
   
   nb = ndim(1)
   if(wannier_gauge) then
      if( ndim(1) .ne. num_wann ) &
      call errore("open_elph_file","mismatch num_wann in elph data",1)
   else
      if( ndim(1) .ne. num_band ) &
      call errore("open_elph_file","mismatch num_band in elph data",1)
   endif
end subroutine open_elph_file


subroutine rearrange_ephmat(nb, iq_st, nq_loc, wannier_gauge, epmat_in, epmat_out)
   use kinds, only: dp
   use input_param, only: kdim
   use lattice_data, only: xq_tot
   use electronic_data, only: xk_tot, num_kpts, ktable, kpt2num, rot_wan, num_wann
   implicit none
   integer, intent(in) :: nb, iq_st, nq_loc
   logical, intent(in) :: wannier_gauge
   complex(dp), intent(in) :: epmat_in(nb, nb, num_kpts, nq_loc, 3)
   complex(dp), intent(out) :: epmat_out(3, num_kpts, nq_loc, num_wann, num_wann)
   !local
   real(dp) :: xkqg(3)
   integer :: ic, jw, iw, ik, iq_g, iq, ikq, idx
   complex(dp) :: ctmp(nb, num_wann), ctmp_wan(num_wann, num_wann), czero, cone
   
   cone  = (1.0_dp, 0.0_dp)
   czero = (0.0_dp, 0.0_dp)

   epmat_out = (0.0_dp, 0.0_dp)
   
   if(wannier_gauge) then
      if(nb .ne. num_wann) call errore('rearrange_ephmat','num_band .ne. num_wann',1)
      do ic = 1, 3
      do jw = 1, num_wann
      do iw = 1, num_wann
         epmat_out(ic, :, :, iw, jw) = epmat_in(iw, jw, :, :, ic)
      enddo; enddo; enddo
   else
   !
   !transform to wannier gauge
!$omp parallel default(shared) private(idx, ik, iq, ic, iq_g, xkqg, ikq, ctmp, ctmp_wan)
!$omp do schedule(guided)
      do idx = 0, num_kpts*nq_loc*3 - 1
         !idx = (ic-1)*num_kpts*nq_loc + (iq-1)*num_kpts + (ik-1)
         ik = mod(idx, num_kpts) + 1
         iq = mod(idx/num_kpts, nq_loc) + 1
         ic = idx / (num_kpts*nq_loc) + 1
         
         iq_g = iq_st + iq - 1
         xkqg = xk_tot(:,ik) + xq_tot(:,iq_g)
         ikq = ktable( kpt2num(xkqg, kdim) )
         
         ctmp = czero
         !g_w(:,:) = matmul( ukq(:,:)^T,  matmul(g_b(:,:), uk(:,:)) )
         ! step 1: matmul(g_b(:,:), uk(:,:))
         call zgemm('n','n', nb, num_wann, nb, cone, &
            epmat_in(1,1, ik,iq,ic), nb, rot_wan(1,1, ik), nb, czero, ctmp, nb)
         ! step 2: matmul( ukq(:,:)^T, ctemp )
         call zgemm('c','n', num_wann, num_wann, nb, &
            cone, rot_wan(1,1,ikq), nb, ctmp, nb, czero, ctmp_wan, num_wann)
         
         epmat_out(ic,ik,iq, :,:) = ctmp_wan(:,:)
      enddo
!$omp end do
!$omp end parallel
   !
   endif
end subroutine


subroutine eph_bloch2wann &
      (xq, xk_all, numk, epmat, nrvec_e, rvec_set_el, nrvec_p, rvec_set_ph, epwan)
   use kinds, only: dp
   use constants, only: tpi
   use wigner_seitz_cell, only: ws_cell, vector_set
   use elph_matrix_wannier, only: eph_wannier
   implicit none
   integer, intent(in) :: numk, nrvec_e, nrvec_p
   real(dp), intent(in) :: xq(3), xk_all(3, numk)  ! in crystal coordinate
   real(dp), intent(in) :: rvec_set_el(3, nrvec_e), rvec_set_ph(3, nrvec_p)
   complex(dp), intent(in) :: epmat(3, numk)
   type(eph_wannier), intent(inout) :: epwan
   !
   integer :: ire, irp, ik
   real(dp) :: rdotk, vec(3)
   complex(dp) :: cfac
   complex(dp), allocatable :: epmatw(:,:)
   type(ws_cell), pointer :: ws_el

   ws_el => epwan%ws_el
   !work array
   allocate( epmatw(3, ws_el%nr) )
   epmatw = cmplx(0.0_dp, 0.0_dp, kind=dp)

!$omp parallel default(shared) private(ire, vec, ik, rdotk, cfac, irp)
!
!$omp do schedule(guided)
   do ire = 1, ws_el%nr
      vec = rvec_set_el(:, ws_el%rvec(ire))
      
      do ik = 1, numk
         rdotk = tpi * dot_product(xk_all(:,ik), vec)
         ! exp(-i k*Re)
         cfac = cmplx(cos(rdotk), -sin(rdotk), kind=dp)
         !
         epmatw(1:3, ire)  = epmatw(1:3, ire) + cfac * epmat(1:3, ik)
      enddo
   enddo
!$omp end do

!$omp do schedule(guided)
   do irp = 1, epwan%ws_ph%nr
      vec = rvec_set_ph(:, epwan%ws_ph%rvec(irp) )
      !
      rdotk = tpi * dot_product(xq, vec)
      ! exp(-i k*Rp)
      cfac = cmplx(cos(rdotk), -sin(rdotk), kind=dp)
      !ep_hop(:,:, irp) = ep_hop(:, :, irp) + epmatw(:, :, iq) * factor
      call zaxpy(3*ws_el%nr, cfac, epmatw(1,1), 1, epwan%ep_hop(1,1,irp), 1)
   enddo
!$omp end do
!
!$omp end parallel

   deallocate( epmatw )
end subroutine eph_bloch2wann

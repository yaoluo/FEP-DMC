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

subroutine compute_scatter_eph_g2(kg, el, ph, elph, qg, wcut, e_thr, use_mem)
   use pert_param, only: prefix, tmp_dir
   use pert_utils, only: abs2, match_table
   use polar_correction,  only: eph_wan_longrange
   use boltz_utils, only: kpts_minus, num2kpt, band2idx, kpts_plus
   use pert_output, only: progressbar_init, progressbar
   use hdf5_utils
   implicit none
   type(grid), intent(in) :: kg
   type(phonon_grid), intent(in) :: qg
   type(electron_wann), intent(in) :: el
   type(lattice_ifc), intent(in) :: ph
   type(elph_mat_wann), intent(in) :: elph
   !
   logical, intent(in) :: use_mem
   real(dp), intent(in) :: wcut, e_thr
   !
   integer :: nmod, qidx, numq
   integer :: i, j, ik, iq, ikq, it, im, n, m, nkq, tot_chnl, bands(3)
   real(dp) :: xq(3), xk(3), xkq(3), wq(ph%nm), ek(el%nb), ekq(el%nb)
   complex(dp) :: uqt(ph%nm, ph%nm), uk(el%nb, el%nb), ukq(el%nb, el%nb), gpol(ph%nm)
   logical :: ltable(kg%numb, kg%numb, ph%nm)
   !
   integer, allocatable :: ist(:), bnd_idx(:)
   real(dp), allocatable :: eph_g2(:), g2(:,:,:)
   complex(dp), allocatable :: uq(:,:,:), gp(:,:), g_kerp(:,:,:,:,:), gkq(:,:,:)
   !
   integer(HID_T) :: file_id
   character(len=120) :: fname, dset_name, msg
   character(len=6), external :: int_to_char
   
   call start_clock('compute_g2')
   !sanity check
   if(elph%nb .ne. el%nb .or. elph%na .ne. ph%na) &
      call errore('compute_scatter_eph_g2','array dimension mismatch',1)
   nmod = ph%nm
   numq = qg%npts

   if( use_mem ) allocate( uq(nmod, nmod, numq) )
   if(elph%lpol) allocate( gp(nmod, numq) )
   if(elph%lpol .or. use_mem) then
!$omp parallel do schedule(guided) default(shared) private(iq, xq, wq, uqt)
      do iq = 1, numq
         xq = num2kpt(qg%list(iq), qg%ndim)
         call solve_phonon_modes(ph, xq, wq, uqt)
         if(elph%lpol) call eph_wan_longrange(elph%pol, xq, gp(:, iq), uqt)
         if(use_mem)  uq(:,:, iq) = uqt(:,:)
      enddo
!$omp end parallel do
   endif

   !open hdf5 file
   fname = trim(tmp_dir) // trim(prefix) // "_eph_g2_p" 
   fname = trim(fname) // trim( int_to_char(my_pool_id+1) ) // ".h5"
   call hdf_open_file(file_id, trim(fname), status='NEW')

   ! allocate work space
   allocate( g_kerp(3, elph%max_nrp, elph%nb, elph%nb, elph%na) )
   ! print out progress info
   call progressbar_init('Computing EPhMatrix:')
   ! main loop
   do i = 1, nk_loc
      ik = kq_pair(i)%ik
      nkq = kq_pair(i)%nkq
      
      !kpts in crystal coordinate
      xk = num2kpt( kg%kpt(ik), kg%ndim )
      call solve_eigenvalue_vector(el, xk, ek, uk)
      call eph_fourier_el_para(elph, xk, g_kerp)
      !
      allocate( ist(nkq) )
      ist(1) = 0
      do j = 1, nkq-1
         ist(j+1) = ist(j) + kq_pair(i)%nchl(j)
      enddo
      !allocate space
      tot_chnl = sum(kq_pair(i)%nchl(:))
      allocate( bnd_idx(tot_chnl), eph_g2(tot_chnl) )
      bnd_idx = 0;   eph_g2 = 0.0_dp

!$omp parallel default(shared) private(j, iq, ikq, xkq, ekq, ukq, xq, wq, uqt, &
!$omp& qidx, gpol, gkq, g2, ltable, it, im, n, m, bands)
      allocate( gkq(el%nb, el%nb, nmod), g2(el%nb, el%nb, nmod) )
!$omp do schedule(guided)
      do j = 1, nkq
         iq = kq_pair(i)%iq(j)
         ikq = kq_pair(i)%ikq(j)
         !
         xkq = num2kpt(kg%kpt(ikq), kg%ndim)
         call solve_eigenvalue_vector(el, xkq, ekq, ukq)
         !
         xq = num2kpt( qg%list(iq), qg%ndim )
         if(use_mem) then
            uqt(:,:) = uq(:,:, iq)
         else
            call solve_phonon_modes(ph, xq, wq, uqt)
         endif
         !
         gpol = (0.0_dp,0.0_dp)
         !NOTE: qg%list(iq) could be (k+q)-k or k-(k+q)
         !  Here we need to check which case it is
         qidx= kpts_minus( kg%kpt(ikq), kg%kpt(ik), kg%ndim, qg%ndim)
         ! iq -> (k+q) - k
         if( qidx .eq. qg%list(iq) ) then
            !uqt(:,:) = uq(:,:, iq)
            if(elph%lpol) gpol(:) = gp(:,iq)
         ! iq -> k - (k+q) .e.g -q
         else if( kpts_plus(qidx, qg%list(iq), qg%ndim, qg%ndim) .eq. 0 ) then
            !NOTE: it's important to use xq = -xq; alternatively we can 
            ! use xq'=num2kpt(qidx), but xq' might be equal to xq + G 
            ! (reciprocal lattice vectors). a small numerical noise due 
            ! to G might break the relation u(q)* = u(-q); 
            xq  = -xq
            !u(q)* = u(-q);  g^L(q)* = g^L(-q)
            uqt(:,:) = dconjg( uqt(:,:) )
            if(elph%lpol) gpol(:) = dconjg( gp(:,iq) )
         else
            !write(stdout,*) kpts_plus(qidx, qg%list(iq), qg%ndim, qg%ndim), qg%list(iq), qidx
            write(msg,'(2(3f12.6,2x))') num2kpt(qidx, qg%ndim), xq
            call errore('compute_scatter_eph_g2','failed: ' // trim(msg), 1)
         endif
         
         !get e-ph matrix elements in wannier gauge and cartesian coord.
         call eph_fourier_elph(elph, xq, g_kerp, gkq)
         call eph_transform_fast(elph, uqt, uk, ukq, gkq, gpol)
         !compute |g|^2
         g2 = abs2( gkq )
         !
         call match_table(kg%enk(:,ikq), kg%enk(:,ik), qg%freq(:,iq), wcut, e_thr, ltable)
         !sanity check
         if( kq_pair(i)%nchl(j) .ne. count(ltable) ) &
            call errore('compute_scatter_eph_g2','kq_pair(i)%nchl(j) .ne. count(ltable)',1)
         
         it = ist(j)
         do im = 1, nmod
         do n = 1, kg%numb
         do m = 1, kg%numb
            if( ltable(m,n,im) ) then
               it = it + 1
               bands = (/m, n, im/)
               bnd_idx(it) = band2idx(bands, kg%numb)
               eph_g2(it) = g2(m+kg%bmin-1, n+kg%bmin-1, im) / (2.0_dp*qg%freq(im,iq))
            endif
         enddo; enddo; enddo
      enddo
!$omp end do
      deallocate(gkq, g2)
!$omp end parallel
      !
      ! output to hdf5 files
      dset_name = "bands_index_" // trim( int_to_char(i) )
      call hdf_write_dataset(file_id, trim(dset_name), bnd_idx)
      dset_name = "eph_g2_" // trim( int_to_char(i) )
      call hdf_write_dataset(file_id, trim(dset_name), eph_g2)
      ! show progress
      call progressbar(i, nk_loc)
      !
      deallocate(bnd_idx, eph_g2, ist)
   enddo
   call hdf_close_file(file_id)
   
   deallocate(g_kerp)
   if( use_mem ) deallocate( uq )
   if(elph%lpol) deallocate( gp )
   
   call mp_barrier(inter_pool_comm)
   call stop_clock('compute_g2')
end subroutine compute_scatter_eph_g2

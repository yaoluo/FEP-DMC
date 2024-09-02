!===============================================================================
! Copyright (C) 2016-2020 Jin-Jian Zhou, Jinsoo Park, I-Te Lu, Marco Bernardi
!
! This program is distributed under the terms of the GNU General Public License.
! See the file `LICENSE' in the root directory of this distribution, or obtain
! a copy of the License at <https://www.gnu.org/licenses/gpl-3.0.txt>.
!
! Author: jjzhou <jjchou.comphy@gmail.com>
! Comment:
!  Calculate the real part of the electronic self-energy
!
!  NOTE: We use openmp parallel for qset.
!  NOTE: Segment fault may occur if openmp private space runs out.
!    we allocated memory for large arrays in each threads. however, the 
!    private stack space each thread has is limited (defined by OMP_STACKSIZE)
!    default value is implememention dependent and is usually quite samll. 
!    Segmentation fault may appear if run of of space. We can specify the size 
!    of the private space by 'export OMP_STACKSIZE=xxK (or xxM, xxG)'
!    or use openmp runtime subroutines: kmp_{set, get}_stacksize_s().
!
! Maintenance:
!===============================================================================

module real_selfenergy
   use pert_const, only: dp, pi, czero, cone
   use qe_mpi_mod, only: ionode, stdout, mp_sum, inter_pool_comm, npool, &
      mp_split_pools, distribute_points
   use pert_param, only: delta_smear, phfreq_cutoff, prefix
   use pert_utils, only: abs2, bose, fermi, gauss, match_table
   use pert_output,only: progressbar_init, progressbar, set_progress_step, stopwatch
   use polar_correction, only: eph_wan_longrange
   use band_structure, only: electron_wann, solve_eigenvalue_vector
   use phonon_dispersion, only: lattice_ifc, solve_phonon_modes
   use elphon_coupling_matrix, only: elph_mat_wann, eph_fourier_el, eph_fourier_elph, &
      eph_transform_fast, eph_transform_polar
   implicit none

   public :: real_selfel_polar, real_selfel
contains

subroutine real_selfel_polar&
      (elph, el, ph, kpts, qset, qwt, tmpr, ef, bmin, enk, resgm)
   implicit none
   type(electron_wann), intent(in) :: el
   type(lattice_ifc),   intent(in) :: ph
   type(elph_mat_wann), intent(in) :: elph
   !tmpr: temperature array; ef: chemical potential array; their size is the same.
   !qwt: q-points (qset) weight. kpts: k-points in the current pool. kpts(3,:)
   real(dp), intent(in) :: kpts(:,:), qset(:,:), qwt(:), tmpr(:), ef(:)
   !compute resigma for bands start from bmin.
   integer, intent(in) :: bmin 
   !resgm(nband, nmod, ntmpr, nkpt), enk(nband, nkpt)
   real(dp), intent(out) :: enk(:,:), resgm(:,:,:,:)

   ! local variables
   integer :: nq_pool, qst, qend, step, nstep, icount, iter
   integer :: numk, numq, ntmpr, numb, bmax, ik, iq, it, im, ib, jb, i, ndim(4)
   real(dp) :: xk(3), xq(3), xkq(3), e_thr, wcut, ph_n, el_fkq, dt1, dt2, scat

   !array
   real(dp) :: ek(el%nb), ekq(el%nb), wq(ph%nm)
   complex(dp) :: gpol(ph%nm), ukq(el%nb, el%nb), uq(ph%nm, ph%nm)
   logical,  allocatable :: ltable(:,:,:)
   real(dp), allocatable :: g2(:,:,:)
   complex(dp), allocatable :: uk(:,:,:), gkq(:,:,:)
   
   call start_clock('real_selfel_polar')
   !initialize 
   enk = 0.0_dp;  resgm = 0.0_dp;   ndim = shape(resgm)
   !#. of selected bands, #.temperature, #. of kpts, #. of q-points
   numb = ndim(1);  ntmpr = ndim(3);  numk = ndim(4);  numq = size(qset,2) 
   bmax = bmin + numb -1
   wcut  = phfreq_cutoff
   e_thr = -10.0_dp ! negative value indicates removal of the threshold.
   
   !array size check
   if(size(tmpr).ne.ntmpr .or. size(ef).ne.ntmpr .or. size(kpts,2).ne.numk &
      .or. size(enk,1).ne.numb .or. size(enk,2) .ne. numk) &
      call errore('real_selfel_polar','arguments dimension mismatch.',1)
   !simple array sanity check
   if(bmax > elph%nb) &
      call errore('real_selfel_polar','bmin and resgm mismatch',1)
   if(elph%nb .ne. el%nb .or. elph%na .ne. ph%na) &
      call errore('real_selfel_polar','array dimension mismatch',1)
   if(.not. elph%lpol) &
      call errore('real_selfel_polar','not a polar material',1)
   if(npool > numq) &
      call errore('real_selfel_polar','too many pools (npool > #. k-points)',1)
   
   !prepare ek and uk
   allocate( uk(el%nb, el%nb, numk) )
!$omp parallel do schedule(guided) default(shared) private(ik, xk, ek)
   do ik = 1, numk
      xk = kpts(:,ik)
      call solve_eigenvalue_vector(el, xk, ek, uk(:,:,ik))
      enk(:,ik) = ek(bmin:bmax)
   enddo
!$omp end parallel do
   
   !distribute q among pools
   call mp_split_pools(numq, qst, qend, nq_pool)
   !for tracking progress
   call set_progress_step(nq_pool, step, nstep)
   if(ionode) call progressbar_init('Im[Sigma] (pol):')

   icount = 0
!$omp parallel default(shared) private(iq, xq, wq, uq, gpol, i, ik, xkq, iter, &
!$omp& ekq, ukq, ltable, gkq, g2, im, ib, jb, dt1, dt2, it, ph_n, el_fkq, scat)
   allocate( ltable(el%nb, numb, ph%nm) )
   allocate( gkq(el%nb, el%nb, ph%nm), g2(el%nb, el%nb, ph%nm) )
!$omp do schedule(guided)
   do iq = qst, qend
      xq = qset(:,iq)
      !get phonon energy and modes
      call solve_phonon_modes(ph, xq, wq, uq)
      !for polar correction
      call eph_wan_longrange(elph%pol, xq, gpol, uq)

      do ik = 1, numk
         xkq = kpts(:,ik) + xq
         call solve_eigenvalue_vector(el, xkq, ekq, ukq)
         !check if there is any scattering channel perserve energy conservation
         call match_table(ekq, enk(:,ik), wq, wcut, e_thr, ltable)
         !skip this k if there is no available scattering channel
         if( .not. any(ltable) ) cycle
         call eph_transform_polar(elph, gpol, uk(:,:,ik), ukq, gkq)
         ! |g|^2
         g2 = abs2(gkq)

         do im = 1, ph%nm
         ! compute resgm for bands between bmin:(bmin+numb-1)
         do  i = 1, numb
         do jb = 1, el%nb
            if(.not. ltable(jb,i,im)) cycle
            dt1 = real_smear( delta_smear, enk(i,ik) - ekq(jb) + wq(im) )
            dt2 = real_smear( delta_smear, enk(i,ik) - ekq(jb) - wq(im) )
            ib = i + bmin - 1
         
            do it = 1, ntmpr
               ph_n = bose(tmpr(it), wq(im)) ! N_q\mu
               el_fkq = fermi( tmpr(it), ekq(jb)-ef(it) )
               scat = (el_fkq + ph_n)*dt1 + (1.0_dp - el_fkq + ph_n)*dt2
               !NOTE: the order of jb and ib are important! Tr[g^\dagger g]
               scat = scat * qwt(iq) * g2(jb, ib, im) / (2.0_dp*wq(im))
!!NB: only one thread can update this shared variable at a time. negligible overhead.
!$omp atomic update
               resgm(i, im, it, ik) = resgm(i, im, it, ik) + scat
            enddo
         enddo; enddo; enddo
      enddo

      !track progress
!$omp atomic update
      icount = icount + 1
      iter = icount
      if( mod(iter, step).eq.0 .or. iter.eq.nq_pool ) then
!$omp critical (real_selfel_polar_progress)
         write(stdout,'(8x, f7.2, a1)') (100.0_dp*iter)/nq_pool, '%'
!$omp end critical (real_selfel_polar_progress)
      endif

   enddo
!$omp end do
   deallocate(gkq, g2, ltable)
!$omp end parallel
   deallocate(uk)
   call mp_sum(resgm, inter_pool_comm)

   call stop_clock('real_selfel_polar')
end subroutine real_selfel_polar


subroutine real_selfel &
      (elph, el, ph, kpts, qset, qwt, tmpr, ef, bmin, enk, resgm, rm_pol)
   !use eph_holstein, only : cal_gkq_holstein
   implicit none
   type(electron_wann), intent(in) :: el
   type(lattice_ifc), intent(in) :: ph
   type(elph_mat_wann), intent(in) :: elph
   !tmpr: temperature array; ef: chemical potential array; their size is the same.
   !qwt: q-points (qset) weight. kpts: k-points in the current pool. kpts(3,:)
   real(dp), intent(in) :: kpts(:,:), qset(:,:), qwt(:), tmpr(:), ef(:)
   !compute resigma for bands start from bmin.
   integer, intent(in) :: bmin 
   !resgm(nband, nmod, ntmpr, nkpt), enk(nband, nkpt)
   real(dp), intent(out) :: enk(:,:), resgm(:,:,:,:) 
   logical, intent(in), optional :: rm_pol !if true, (g^2-g_pol^2) is used.
   ! local variables
   integer :: mem_max, step_max, nk_pool, qst_pool, qend_pool, nq_pool, nstep_k, ikc
   integer :: step_q, nstep_q, istep_k, kst, kend, step_k, n, icount, iter
   integer :: numk, numq, ntmpr, numb, ik, iq, it, im, ib, jb, i, bmax, ndim(4), neph
   real(dp) :: xk(3), xq(3), xkq(3), wcut, e_thr, ph_n, el_fkq, dt1, dt2, scat
   
   real(dp) :: ek(el%nb), wq(ph%nm), ekq(el%nb)
   complex(dp) :: gpol(ph%nm),gpolc(ph%nm), ukq(el%nb, el%nb), uq(ph%nm, ph%nm),uqc(ph%nm, ph%nm)
   logical,  allocatable :: ltable(:,:,:)
   real(dp), allocatable :: g2(:,:,:)
   complex(dp), allocatable :: uk(:,:,:), gkq(:,:,:),gkqc(:,:,:), g_kerp(:,:,:,:,:,:),gkq_q(:,:,:)
   integer, pointer :: ik_loc(:)


   !call init_eph_holstein(el,ph,elph)
   call start_clock('real_selfel')
   enk = 0.0_dp;  resgm = 0.0_dp;  ndim = shape(resgm)
   !#. of selected bands, #.temperature, #. of kpts, number of q-points
   numb = ndim(1);  ntmpr = ndim(3);  numk = ndim(4);  numq = size(qset,2)
   bmax = bmin + numb -1
   !array sanity check
   if(size(tmpr).ne.ntmpr .or. size(ef).ne.ntmpr .or. size(kpts,2).ne.numk &
      .or. size(enk,1).ne.numb .or. size(enk,2) .ne. numk) &
      call errore('real_selfel','arguments dimension mismatch.',1)
   if(bmax > elph%nb) &
      call errore('real_selfel','bmin and resgm mismatch',1)
   if(elph%nb .ne. el%nb .or. elph%na .ne. ph%na) &
      call errore('real_selfel','array dimension mismatch',1)
   if(present(rm_pol) .and. rm_pol .and. (.not. elph%lpol)) &
      call errore('real_selfel',"polar_split='rmpol' is for polar material",1)

   e_thr = -10.0_dp ! negative value indicates removal of the threshold.
   wcut  = phfreq_cutoff
   neph = el%nb * el%nb * ph%na * 3
   !max memory per process for g_kerp, hard-coded here, not sure what is the optimal value
   mem_max  = 1024*1024*1024/2 !8GB
   step_max = mem_max / ( neph * elph%max_nrp )

   ik_loc => null()
   if(numk > 2*step_max*npool) then
      !interleave distribution of kpoints over pools, for better load balance
      call distribute_points(numk, nk_pool, ik_loc, .true.)
      !call mp_split_pools(numk, kst_pool, kend_pool, nk_pool)
      qst_pool = 1;     qend_pool = numq;    nq_pool = numq
   else
      call distribute_points(numk, nk_pool, ik_loc, .false.)
      !kst_pool = 1;     kend_pool = numk;    nk_pool = numk
      !distrbite q among pools
      call mp_split_pools(numq, qst_pool, qend_pool, nq_pool)
   endif
   !
   call set_progress_step(nk_pool, step_k, nstep_k, step_max)
   if(nstep_k >= 8) then
      !tracking progess on k only
      step_q = nq_pool;    nstep_q = 1
   else
      !track step on k and q combined
      call set_progress_step(nq_pool, step_q, nstep_q)
   endif
   !output debug info
   !write(stdout,'(5x, 3(a,i8,1x))') &
   !   'step_k:', step_k, 'nstep_k:', nstep_k, 'step_max:', step_max
   
   !allocate work sapce
   allocate( uk(el%nb, el%nb, step_k), &
      g_kerp(3, elph%max_nrp, elph%nb, elph%nb, elph%na, step_k) )

   if(ionode) call progressbar_init('Re[Sigma]:')
   icount = 0
   do istep_k = 1, nstep_k
      kst = (istep_k-1)*step_k + 1
      kend = min(kst+step_k-1, nk_pool)
      
      uk = czero;  g_kerp = czero
!$omp parallel do schedule(guided) default(shared) private(n, ikc, ik, xk, ek)
      do n = kst, kend
         ikc = n - kst + 1
         !electronic wavefunction at ik
         ik = ik_loc(n)
         xk = kpts(:, ik)
         call solve_eigenvalue_vector(el, xk, ek, uk(:,:,ikc))
         enk(:,ik) = ek(bmin:bmax)
         call eph_fourier_el(elph, xk, g_kerp(:,:,:,:,:,ikc))
      enddo
!$omp end parallel do
      
      !
!$omp parallel default(shared) private(iq, xq, wq, uq, gpol, i, n, ik, ikc, xkq, &
!$omp& ekq, ukq, ltable, gkq, g2, im, ib, jb, dt1, dt2, it, ph_n, el_fkq, scat, iter,uqc,gkqc,gpolc,gkq_q)
      allocate( ltable(el%nb, numb, ph%nm) )
      allocate( gkq(el%nb, el%nb, ph%nm), gkq_q(el%nb, el%nb, ph%nm),gkqc(el%nb, el%nb, ph%nm), g2(el%nb, el%nb, ph%nm) )
!$omp do schedule(guided)
      do iq = qst_pool, qend_pool
         xq = qset(1:3, iq)
         !get phonon frequcies and eigen-displacement at iq
         call solve_phonon_modes(ph, xq, wq, uq)
         !get e-ph matrix elements in wannier gauge and cartesian coord.
         gpol = czero; 
         !gpolc = czero
         !uqc = conjg(uq)
         if(elph%lpol) call eph_wan_longrange(elph%pol, xq, gpol, uq)
         !if(elph%lpol) call eph_wan_longrange(elph%pol, -xq, gpolc, uqc)
         !call cal_gkq_holstein(xq, gkq_q )
         do n = kst, kend
            ikc = n - kst + 1
            !electronic wavefunction at ik
            ik = ik_loc(n)
            xkq = kpts(:,ik) + xq
            !get electronic wavefunction at ikq
            call solve_eigenvalue_vector(el, xkq, ekq, ukq)
            !check if there is any scattering channel perserve energy conservation
            !call match_table(ekq, enk(:,ik), wq, wcut, e_thr, ltable)
            !if(.not. any(ltable)) cycle
            !
            call eph_fourier_elph(elph, xq, g_kerp(:,:,:,:,:,ikc), gkq)
            !gkq = gkq_q
            !transfor to phonon mode and bloch gauge.
            call eph_transform_fast(elph, uq, uk(:,:,ikc), ukq, gkq, gpol)
            !compute |g|^2
            g2 = abs2( gkq )
            if(present(rm_pol) .and. rm_pol .and. elph%lpol) then
               call eph_transform_polar(elph, gpol, uk(:,:,ikc), ukq, gkq)
               g2 = g2 - abs2(gkq)
            endif

            do im = 1, ph%nm
            if(wq(im)<phfreq_cutoff) cycle 
            ! compute resgm for bands between bmin:(bmin+numb-1)
            do  i = 1, numb
            do jb = 1, el%nb
            !do jb = 1, 6
            !do jb = 8,12
               !if(.not. ltable(jb,i,im)) cycle
               dt1 = real_smear( delta_smear, enk(i,ik) - ekq(jb) + wq(im) )
               dt2 = real_smear( delta_smear, enk(i,ik) - ekq(jb) - wq(im) )
               ib = i + bmin - 1
            
               do it = 1, ntmpr
                  !ph_n = bose(tmpr(it), wq(im)) ! N_q\mu
                  el_fkq = fermi( tmpr(it), ekq(jb)-ef(it) )
                  !scat = (el_fkq + ph_n)*dt1 + (1.0_dp - el_fkq + ph_n)*dt2
                  scat = (el_fkq)*dt1 + (1.0_dp - el_fkq)*dt2
                  !scat = 1.0_dp * dt2  !zero temperature approx 
                  !NOTE: the order of jb and ib are important! Tr[g^\dagger g]
                  scat = scat * qwt(iq) * g2(jb, ib, im) / (2.0_dp*wq(im))
!!NB: only one thread can update this shared variable at a time. negligible overhead.
!$omp atomic update
                  resgm(i, im, it, ik) = resgm(i, im, it, ik) + scat
               enddo
            enddo; enddo; enddo
         enddo
      
         !track progress
!$omp atomic update
         icount = icount + 1
         iter = icount
         if( mod(iter, step_q*nstep_k).eq.0 .or. iter.eq.(nq_pool*nstep_k) ) then
!$omp critical (real_selfel_progress)
            write(stdout,'(8x, f7.2, a1)') (100.0_dp*iter)/(nq_pool*nstep_k), '%'
!$omp end critical (real_selfel_progress)
         endif

      enddo
!$omp end do
      deallocate(gkq, g2, ltable)
!$omp end parallel
   enddo
   deallocate(uk, g_kerp)
   call mp_sum(resgm, inter_pool_comm)
   if(numk .ne. nk_pool) call mp_sum(enk, inter_pool_comm)
   if( associated(ik_loc) ) deallocate( ik_loc )

   call stop_clock('real_selfel')
end subroutine real_selfel

real(dp) pure function real_smear(eta, x) result(f)
   implicit none
   real(dp), intent(in) :: eta, x
   f = x / (x*x + eta*eta)
end function real_smear

end module real_selfenergy

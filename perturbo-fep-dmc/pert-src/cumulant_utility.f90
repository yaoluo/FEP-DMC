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

module cumulant_utility
   use pert_const, only: dp, pi, ryd2ev, kelvin2ev
   use pert_utils, only: find_free_unit
   use qe_mpi_mod, only: ionode, ionode_id, inter_pool_comm, mp_bcast
   use hdf5_utils
   implicit none
   private
   
   !\omega grid for cumulant calculation
   !*to improve efficiency, double grids are used:
   !  the inner energy window with a small energy step:  de
   !  the outer energy window with a large energy step: np*de
   !*both the inner and outer window enclose w=0 point, 
   !  and outer window enclose the inner window
   type, public :: cum_wgrid
      real(dp) :: de  ! energy step 
      integer  :: np  ! np*de is the step of outer window
      integer  :: in_lw, in_up  ! inner window: in_lw*de:in_up*de
      integer  :: in_lw_idx, in_up_idx  !make sure up_idx >= lw_idx+2
      !outer window: lower*(de*np):upper*(de*np)
      integer  :: lower, upper
      real(dp), pointer ::  w(:)=>null()
      !expand to fine grid for spectral calculation
      real(dp) :: sp_de
      integer  :: sp_np  !sp_de = de / sp_np
      integer  :: out_lw, out_up
   end type cum_wgrid

   public :: init_cum_wgrid, compute_beta, save_selfenergy, load_selfenergy, &
      load_spectral, spectral_cut_tail, ouput_integrated_dos_tdf
contains

subroutine init_cum_wgrid &
   (outer_emin, outer_emax, inner_emin, inner_emax, de, np, wg, npf)
   implicit none
   real(dp), intent(in) :: inner_emin, inner_emax, outer_emin, outer_emax, de
   integer,  intent(in) :: np
   integer,  intent(in), optional :: npf
   type(cum_wgrid), intent(out) :: wg
   !local
   integer :: up, lw, n_co, n_ci, n_fi, ntot, i, ic
   real(dp) :: large_step

   if(inner_emin > 0.0_dp .or. inner_emax < 0.0_dp) &
      call errore('init_cum_wgrid','(emin, emax) should enclose 0.0',1)
   if( outer_emin > inner_emin .or. outer_emax < inner_emax ) &
      call errore('init_cum_wgrid','outer energy window should enclose the inner window',1)
   if(np < 0 .or. (present(npf) .and. npf < 0)) &
      call errore('init_cum_wgrid','negative np or npf',1)
   
   wg%np = np
   wg%de = de
   large_step = real(np, dp) * de
   wg%upper = int( ceiling(outer_emax/ large_step) )
   wg%lower = int(   floor(outer_emin/ large_step) )
   
   !define fine energy grid for spectral function 
   if(present(npf)) then
      wg%sp_np = npf
   else
      wg%sp_np = 1
   endif
   wg%sp_de = wg%de / real(wg%sp_np, dp)
   wg%out_up = wg%upper * wg%np * wg%sp_np
   wg%out_lw = wg%lower * wg%np * wg%sp_np
   
   !calculate wg%w grid used in self-energy calculation.
   up = int( ceiling(inner_emax/large_step) )
   lw = int(   floor(inner_emin/large_step) )
   
   wg%in_up = up * np
   wg%in_lw = lw * np
   !number of points within outer window with step de*np
   n_co = wg%upper - wg%lower + 1
   !number of points within inner window with step de*np
   n_ci = up - lw + 1
   !number of points within inner window with step de
   n_fi = wg%in_up - wg%in_lw + 1
   !total number of points
   ntot = n_co - n_ci + n_fi
   allocate( wg%w( ntot ) )

   ic = 0     
   do i = wg%lower, lw-1
      ic = ic + 1
      wg%w(ic) = real(i*np, dp) * de
   enddo
   
   wg%in_lw_idx = ic + 1
   do i = wg%in_lw, wg%in_up
      ic = ic + 1
      wg%w(ic) = real(i, dp) * de
   enddo
   wg%in_up_idx = ic

   do i = up+1, wg%upper
      ic = ic + 1
      wg%w(ic) = real(i*np, dp) * de
   enddo
   
   !sanity check
   if( wg%in_up_idx < wg%in_lw_idx+2 ) call errore('init_cum_wgrid','in_up_idx < in_lw_idx+2',1)
   if(ic .ne. ntot) call errore('init_cum_wgrid','ic .ne. ntot',1)
end subroutine init_cum_wgrid


subroutine save_selfenergy(fname, enk, se_imag, wgrid)
   implicit none
   character(*), intent(in) :: fname
   !se_real(nestep, ntemper, numb, numk)
   real(dp), intent(in) :: enk(:,:), se_imag(:,:,:,:), wgrid(:)
   !
   integer(HID_T) :: file_id

   call hdf_open_file(file_id, trim(fname), status='NEW')
   
   call hdf_write_dataset(file_id, 'enk', enk)
   call hdf_write_dataset(file_id, 'wgrid', wgrid)
   call hdf_write_dataset(file_id, 'selfenergy_imag', se_imag)

   call hdf_close_file(file_id)
end subroutine save_selfenergy


!read in fan selfenergy
subroutine load_selfenergy(fname, enk, se_imag, wgrid)
   implicit none
   character(*), intent(in) :: fname
   !se_real(nestep, ntemper, numb, numk)
   real(dp), intent(out) :: enk(:,:), se_imag(:,:,:,:)
   real(dp), intent(in), optional :: wgrid(:)
   !local
   logical :: has_file
   integer(HID_T) :: file_id
   integer :: nestep, numb, nkpt
   real(dp), allocatable :: wtmp(:)

   numb = size(enk, 1)
   nkpt = size(enk, 2)
   nestep = size(se_imag, 1)
   write(*,'(A20,2i5)')'ne = ',nestep,size(wgrid)
   if(size(se_imag,3).ne.numb .or. size(se_imag,4).ne.nkpt) &
      call errore('load_selfenergy','mismatch enk and se_imag',1)

   se_imag = 0.0_dp;  enk = 0.0_dp
   if(ionode) then
      inquire(file=trim(fname), exist=has_file)
      if(.not. has_file) call errore('load_selfenergy','missing '// trim(fname), 1)
      !
      allocate( wtmp(nestep) )
      call hdf_open_file(file_id, trim(fname), status='OLD', action='READ')
      call hdf_read_dataset(file_id, 'enk', enk)
      call hdf_read_dataset(file_id, 'wgrid', wtmp)
      call hdf_read_dataset(file_id, 'selfenergy_imag',se_imag)
      call hdf_close_file(file_id)
      
      if(present(wgrid)) then
         if( size(wgrid).ne.nestep .or. any(abs(wgrid - wtmp) > 1.0E-8_dp) ) &
            call errore('load_selfenergy','inconsistence freq. grid.',1)
      endif
      deallocate(wtmp)
   endif
   call mp_bcast(enk, ionode_id, inter_pool_comm)
   call mp_bcast(se_imag, ionode_id, inter_pool_comm)
   return
end subroutine load_selfenergy


subroutine load_spectral(fname, enk, spectral)
   implicit none
   character(*), intent(in) :: fname
   !se_real(:, ntemper, numb, numk)
   real(dp), intent(out) :: enk(:,:), spectral(:,:,:,:)
   character(len=6), external :: int_to_char
   !local
   logical :: has_file
   integer :: ik, numk
   character(len=80) :: dset_name
   integer(HID_T) :: file_id, group_id
   
   numk = size(spectral, 4)
   if( size(enk, 2) .ne. numk) &
      call errore('load_spectral','inconsistent enk, spectral',1)

   enk = 0.0_dp;   spectral = 0.0_dp
   if(ionode) then
      inquire(file=trim(fname), exist=has_file)
      if(.not. has_file) call errore('load_spectral','missing '// trim(fname), 1)
      !
      call hdf_open_file(file_id, trim(fname), status='OLD', action='READ')
      call hdf_read_dataset(file_id, 'band_energy', enk)
      !
      call hdf_open_group(file_id, 'spectral_functions', group_id)
      do ik = 1, numk
         dset_name = "kpt_" // trim(int_to_char(ik))
         call hdf_read_dataset(group_id, trim(dset_name), spectral(:,:,:,ik))
      enddo
      call hdf_close_group(group_id)
      call hdf_close_file(file_id)
   endif

   call mp_bcast(enk, ionode_id, inter_pool_comm)
   do ik = 1, numk
      call mp_bcast(spectral(:,:,:,ik), ionode_id, inter_pool_comm)
   enddo
end subroutine load_spectral


!compute beta(w) from the ImSigma(w): beta(w) = -ImSigma(w)/pi, interpolation if needed.
subroutine compute_beta(wg, se_imag, beta)
   use pert_utils, only: spline_coef, spline_interp_interval
   implicit none
   type(cum_wgrid), intent(in) :: wg
   real(dp), intent(in) :: se_imag(:)
   real(dp), intent(out) :: beta(wg%out_lw:wg%out_up)
   !lacal
   integer :: nw, np, n, i, j, nlen, npf
   real(dp), allocatable :: dx(:), y_out(:), egrid(:), y(:), b(:), c(:), d(:)
   
   npf = wg%sp_np
   !no need for interpolation in inner window
   if(npf .eq. 1) then
      n = wg%in_lw_idx  !start index of the inner window
      do i = wg%in_lw, wg%in_up
         beta(i) = -se_imag(n)/pi
         n = n + 1
      enddo
   else
      !inner window interpolation
      nlen = wg%in_up_idx - wg%in_lw_idx + 1  !number of points within inner windows
      allocate(egrid(nlen), y(nlen), b(nlen), c(nlen), d(nlen))
      b = 0.0_dp;  c = 0.0_dp;  d = 0.0_dp
      n = 0
      do i = wg%in_lw_idx, wg%in_up_idx
         n = n + 1
         egrid(n) = wg%w(i)
         y(n) = -se_imag(i)/pi   ! beta(w) = - ImSigma(w) / pi
      enddo
      !compute coeffieient for interpolation
      call spline_coef(egrid, y, b, c, d)
      !w-grid for interpolation
      allocate( dx(npf), y_out(npf) )
      do i = 1, npf
         dx(i) = real(i-1, dp) * wg%sp_de
      enddo
      !do interpolation
      n = wg%in_lw * npf
      do i = 1, nlen-1
         call spline_interp_interval(y(i), b(i), c(i), d(i), dx, y_out)
         beta(n:(n+npf-1)) = y_out(1:npf)
         n = n + npf
      enddo
      beta(n) = y(nlen)  !the end point.
      deallocate(dx, y_out, egrid, y, b, c, d)
   endif

   !if inner window overlap with outer window, then no need for interpolation.
   if(wg%in_lw .eq. wg%out_lw  .and.  wg%in_up .eq. wg%out_up) return
   
   !for outer window interpolation
   np = wg%np * wg%sp_np
   !dx, y_out for interpolation
   allocate( dx(np), y_out(np) )
   do i = 1, np
      dx(i) = real(i-1, dp) * wg%sp_de
   enddo
   
   !interpolate between wg%out_lw and wg%in_lw
   if(wg%in_lw_idx > 1) then
      nlen = wg%in_lw_idx + 2   !add extra 2 data points for interpolation.
      allocate(egrid(nlen), y(nlen), b(nlen), c(nlen), d(nlen))
      b = 0.0_dp;  c = 0.0_dp;  d = 0.0_dp
      do i = 1, nlen
         egrid(i) = wg%w(i)
         y(i) = -se_imag(i)/pi   ! beta(w) = - ImSigma(w) / pi
      enddo
      !compute coeffieient for interpolation
      call spline_coef(egrid, y, b, c, d)
      !do interpolation
      n = wg%out_lw
      do i = 1, nlen-3
         call spline_interp_interval(y(i), b(i), c(i), d(i), dx, y_out)
         beta(n:(n+np-1)) = y_out(1:np)
         n = n + np
      enddo
      deallocate(egrid, y, b, c, d)
   endif

   nw = size(wg%w)
   !interpolate between wg%in_up and wg%out_up
   if( wg%in_up_idx < nw) then
      nlen = nw - (wg%in_up_idx-2) + 1 !add extra 2 data points for interpolation.
      allocate(egrid(nlen), y(nlen), b(nlen), c(nlen), d(nlen))
      b = 0.0_dp;  c = 0.0_dp;  d = 0.0_dp
      !mapping beta(w) = beta(-w); and interpolating beta(-w)
      n = nw
      do i = 1, nlen
         egrid(i) = -wg%w(n)
         y(i) = -se_imag(n)/pi
         n = n - 1
      enddo
      !compute coeffieient for interpolation
      call spline_coef(egrid, y, b, c, d)
      !do interpolation
      n = wg%out_up
      do i = 1, nlen-3
         call spline_interp_interval(y(i), b(i), c(i), d(i), dx, y_out)
         do j = 0, np-1
            beta(n-j) = y_out(j+1)
         enddo
         n = n - np
      enddo
      deallocate(egrid, y, b, c, d)
   endif

   !beta(w) should be > 0, set it to 0 if it's very small
   do i = wg%out_lw, wg%out_up
      beta(i) = merge( beta(i), 0.0E0_dp, beta(i) > 1.0E-40_dp )
   enddo
end subroutine compute_beta


integer pure function spectral_cut_tail(specfunc, lb_spf, dw, weight_cut) result(lb_cut)
   implicit none
   integer,  intent(in) :: lb_spf
   real(dp), intent(in) :: specfunc(lb_spf:), weight_cut, dw
   !local
   real(dp) :: weight
   integer  :: ub_spf, i

   ub_spf = ubound(specfunc, 1)
   weight = 0.0_dp
   do i = lb_spf, ub_spf
      weight = weight + specfunc(i) * dw
      if( weight > weight_cut ) exit
   enddo
   lb_cut = i
end function spectral_cut_tail


subroutine ouput_integrated_dos_tdf(fname, energy, idos, itdf, tmpr, ef)
   use pert_data, only: volume, alat
   implicit none
   character(*), intent(in) :: fname
   real(dp), intent(in) :: energy(:), idos(:,:), itdf(:,:), tmpr(:), ef(:)
   !local 
   real(dp) :: factor 
   integer :: nestep, ntemp, uout, i, it

   nestep = size(energy)
   ntemp = size(tmpr)
   if( ntemp.ne.size(ef) .or. ntemp.ne.size(idos,2) .or. ntemp.ne.size(itdf,2) &
      .or. nestep.ne.size(idos,1) .or. nestep.ne.size(itdf,1) ) &
      call errore('ouput_dos_tdf','inconsistent argument dimension', 1)
   
   factor = alat*alat/volume
   uout = find_free_unit()
   open(unit=uout, file=trim(fname), status='unknown', form='formatted')
   write(uout,'(2x,"#   E (eV)    Int_DOS (#.states/u.c.)   Int_TDF (a.u.)",/)')
   do it = 1, ntemp
      write(uout,'(/, a, f9.4, a, f10.6, /)') &
      '#  Temperature: ', tmpr(it)*ryd2ev/kelvin2ev, '  Fermi Level: ', ef(it)*ryd2ev
      do i = 1, nestep
         write(uout,'(f12.6, 3x, 6(E14.6,1x))') energy(i)*ryd2ev, idos(i,it), itdf(i,it)*factor
      enddo 
   enddo
   close(uout) 
end subroutine ouput_integrated_dos_tdf

!subroutine cum_trans_spectral_fake(kg, spf, lb_spf, dw)
!   use pert_output, only: load_imsigma
!   implicit none
!   type(grid), intent(in) :: kg
!   integer, intent(in) :: lb_spf  !lower boundary of the leading dimension of spf.
!   real(dp), intent(in) :: dw
!   !spf(lowerb:upper, ntemper, numb, numk)
!   real(dp), intent(out) :: spf(lb_spf:, :, :, :)
!   !local variables
!   logical :: success
!   integer :: ntemp, numk, numb, ik, ib, it, i, ndim(4), ub_spf
!   real(dp) :: de, sgm
!   real(dp), allocatable :: imsgm(:,:,:)
!
!   ndim = shape(spf);  ntemp = ndim(2);   numb = ndim(3);   numk = ndim(4)
!   ub_spf = ndim(1) + lb_spf - 1  ! nestep = ub_spf - lb_spf + 1
!
!   allocate( imsgm(numb, numk, ntemp) )
!   call load_imsigma(imsgm, success)
!   if(.not. success) call errore('cum_trans_spectral_fake','fail to load_imsigma',1)
!   if( any(imsgm < 1.0E-16_dp) ) &
!      call errore('cum_trans_spectral_fake','Warning: negative imsigma (or less than 1.0E-16)',-1)
!   
!   do ik = 1, numk
!   do ib = 1, numb
!   do it = 1, ntemp
!      sgm = imsgm(ib, ik, it)  ! Im(Sigma)(E_nk)
!      do i = lb_spf, ub_spf
!         de = real(i, dp)*dw
!         spf(i, it, ib, ik) = sgm / pi / (de**2 + sgm**2)
!      enddo
!   enddo; enddo; enddo
!
!   deallocate(imsgm)
!end subroutine cum_trans_spectral_fake
!
!
!!compute cumulant spectral for k on the mesh
!subroutine cum_trans_spectral_qp(beta, spf, lb, e_qp, dw)
!   use cumulant_expansion, only: spectral_qp
!   implicit none
!   integer,  intent(in) :: lb  !lower boundary of lead dim of array beta and spf.
!   !beta(w) = -Im[Sigma](w+Enk) / pi; since ImSigma is negative, -Imsigma = |ImSigma|
!   !beta(lb:ub, ntemper, numb, numk), e_qp = se_real(1-lb, it, ib, ik), Re[Imsigma](w=Enk)
!   real(dp), intent(in) :: beta(lb:, :, :, :), e_qp(:,:,:), dw !e_qp(ntemper, numb, numk)
!   real(dp), intent(out) :: spf(lb:, :, :, :)  !spf(lb:ub, ntemper, numb, numk)
!   !local
!   integer :: nestep, ntemp, numb, numk, ub, j, jst, jend, ik, ib, it, i
!   
!   nestep = size(beta, 1);    ntemp = size(beta, 2); 
!   numb   = size(beta, 3);    numk  = size(beta, 4);
!   if(size(beta) .ne. size(spf) .or. size(e_qp) .ne. (ntemp*numb*numk) ) &
!      call errore('cum_trans_spectral_qp','size mismatch: beta, spf',1)
!   
!   ub = nestep + lb - 1  !nestep = ub - lb + 1
!   
!   spf = 0.0_dp
!   call mp_split_pools( ntemp*numb*numk, jst, jend ) ! mpi parallel 
!   do j = jst, jend
!      !j = (ik-1)*(numb*ntemp) + (ib-1)*ntemp + (it-1) + 1
!      ik = (j-1) / (numb*ntemp) + 1
!      ib = mod( (j-1)/ntemp, numb) + 1
!      it = mod( (j-1), ntemp ) + 1
!
!      do i = lb, ub
!         spf(i,it,ib,ik) = spectral_qp(beta(0,it,ib,ik), e_qp(it,ib,ik), real(i,dp)*dw)
!      enddo
!      !check weight
!      if( abs(sum(spf(:,it,ib,ik))*dw - 1.0_dp) > 0.1_dp ) &
!         call errore('cum_trans_spectral_qp', 'spectral weight deviate from 1', 0)
!   enddo
!   call mp_sum(spf, inter_pool_comm)
!end subroutine cum_trans_spectral_qp


end module cumulant_utility

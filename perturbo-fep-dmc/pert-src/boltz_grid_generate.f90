!===============================================================================
! Copyright (C) 2016-2020 Jin-Jian Zhou, Jinsoo Park, I-Te Lu, Marco Bernardi
!
! This program is distributed under the terms of the GNU General Public License.
! See the file `LICENSE' in the root directory of this distribution, or obtain 
! a copy of the License at <https://www.gnu.org/licenses/gpl-3.0.txt>.
!
! Author: jjzhou <jjchou.comphy@gmail.com>
! Comment:
!  - at a given kmesh defined by kdim, select and output tetrahedra that
!    contributed to BZ integration. 
!  - collect vertex of all selected tetrahedra, but only output the 
!    irreducible into tet.kpt file.
!  - only levels between [bmin, bmax] are considered.  only resluts (dos) 
!    inside [emin, emax] are wanted.
!
! Maintenance:
!===============================================================================

subroutine boltz_grid_generate(kg, kdim, emin, emax, bmin, bmax, nsym, symop, el)
   use pert_const, only: dp
   use pert_param, only: prefix
   use pert_utils, only: find_free_unit 
   use boltz_utils, only: num2kpt, kpt2num, inside_win
   use band_structure, only: electron_wann
   use qe_mpi_mod,only: ionode, ionode_id, mp_split_pools, &
      mp_bcast, mp_barrier, mp_sum, inter_pool_comm ,stdout
   implicit none
   type(grid), intent(inout) :: kg
   ! energy range
   real(dp), intent(in) :: emin, emax
   ! dimension of the k-grid; band index range; symmetry operations.
   integer, intent(in) :: kdim(3), bmin, bmax, nsym, symop(3,3,nsym)
   type(electron_wann), intent(in), optional :: el
   ! local variables
   integer :: nkr, ik, num_irk, ns, i, kdx, ist, iend, ntet, nk
   integer, allocatable:: keq(:), irrk(:), irrk_e(:), irrk_t(:)
   real(dp) :: xkg(3), xkr(3)

   call start_clock('Grid_gen')
   nkr = kdim(1)*kdim(2)*kdim(3)
   ! nkr should be smaller than 2147483647, otherwise integer overflow occurs
   if (nkr <= 0) call errore('boltz_grid_generate', &
      'illegal kdim (or integer overflow)', 1)
   kg%ndim = kdim
   kg%kweight = 1.0_dp/dble(kdim(1))/dble(kdim(2))/dble(kdim(3))
   kg%tweight = kg%kweight/6.0_dp

   !keq might be a very large array, only allocate it on root process.
if (ionode) then
   !
   allocate( keq(0:nkr-1) )
   do ik = 0, nkr-1
      keq(ik) = ik
   enddo
   ! find irreducible k-points
   i = 0 !number of irreducible k-points
   do ik = 0, nkr-1
      !skip reducible k-points
      if( keq(ik) .ne. ik ) cycle
      i = i + 1
      !get crystal coordinate of k-points
      xkg(1:3) = num2kpt(ik, kdim)
      
      ! loop over point-group opterations, time-reversal is not included.
      do ns = 1, nsym
         xkr(1:3) = matmul(symop(1:3,1:3,ns), xkg(1:3))
         !map back to index coordinate.
         kdx = kpt2num(xkr, kdim)
         ! kdx < 0 : xkr is not in the grid
         ! from all the equivalent k-poins, choose the one with the smallest kdx
         if ( kdx > ik ) then
            keq( kdx ) = ik
         !sanity check
         elseif ( kdx >= 0 .and. kdx < ik) then
            call errore('boltz_grid_generate','Error in finding irreducible k', 1)
         endif
      enddo
   enddo
   num_irk = i
   !collect irreducible k-points and check consistence
   allocate( irrk(num_irk) )
   i = 0
   do ik = 0, nkr-1
      ! this is a irreducible kpoint
      if(keq(ik) .eq. ik) then 
         i = i + 1
         !from ik, we can get the crystal coord of this kpoints
         irrk(i) = ik
         !for irreducible k, point to its position in irrk.
         keq(ik) = i
      else
         !for reducible k, point to the position of their corresponidng irr-k
         ! the implicit assumption here is that, for reducible k, 'ik' must be
         ! larger than the 'ik' of their corresponing irreducbile kpoints.
         keq(ik) = keq( keq(ik) )
      endif
      ! Up to now, irrk(i) stores the i-th irreducible points.(in index coord)
      ! and, keq(ik) point the corresponding irreducible k in irrk.
   enddo
endif
   call mp_bcast(num_irk, ionode_id, inter_pool_comm)
   if(.not. allocated(irrk) ) allocate( irrk(num_irk) )
   call mp_bcast(irrk, ionode_id, inter_pool_comm)
   
   allocate( irrk_e(num_irk) )
   if( present(el) ) then
      irrk_e(:) = 0
      !this part is time consuming, split num_irk over pools (using mpi+openmp).
      call mp_split_pools(num_irk, ist, iend)
!! openmp parallel.
!$omp parallel do schedule(guided) default(shared) private(i, xkg)
      do i = ist, iend
         xkg(1:3) = num2kpt(irrk(i), kdim)
         ! check whether this irr-k has energy levels inside [emin, emax]
         ! only bands in [bmin, bmax] are counted.
         irrk_e(i) = inside_win(el, xkg, emin, emax, bmin, bmax)
         ! >= 0: no level inside [emin, emax], all are lower than emin or larger
         !  than emax. or some lower than emin and the others larger than emax.
         !   the value is the highest band that have energy lower than emin.
         ! = -2: have levels inside [emin, emax]
         ! = -4 : error, illegal result.
         if(irrk_e(i) .eq. -4) call &
            errore('boltz_grid_generate','Error in finding bands in [emin, emax]',1)
      enddo
!$omp end parallel do
      call mp_sum(irrk_e, inter_pool_comm)
   else
      !select all the kpoints
      irrk_e(:) = -2
   endif

if (ionode) then
   call count_tetra(kdim, keq, num_irk, irrk_e, ntet)
   kg%num_tet = ntet
   allocate( kg%tetra(2, 4, ntet) )
   ! irrk_t is used to reorder irreducible-k according to the order of 
   ! their appearance in tetrahedron. 0 if not selected.
   allocate( irrk_t(num_irk) )
   call collect_tetra(kdim, keq, num_irk, irrk_e, ntet, nk, irrk_t, kg%tetra)
   kg%nk = nk
   kg%nk_irr = count(irrk_t > 0)
   allocate(kg%irrk( kg%nk_irr ))
   !collect selected irreducible kpts 
   kg%irrk = -1
   do i = 1, num_irk
      if(irrk_t(i) > 0) then
         ! irrk(i) : the index coordinate of i-th irreducible-kpoint
         ! irrk_t(i) : the new location of the i-th irreducible-kpoints
         kg%irrk( irrk_t(i) ) = irrk(i)
      endif
   enddo
   if(any( kg%irrk < 0 )) call errore('boltz_grid_generate', 'illegal irrk.', 1)
   deallocate(irrk_t, keq)
   
   !recover the full reducible grid 
   allocate( kg%kpt( kg%nk ), kg%kpt2ir( kg%nk ) )
   ! collect all the reducible k-points into kg%kpt, and their corresponding irreducible k.
   ! and also update thet tetra, so that
   !  tetra(1, ic, it) -> index of irreducible k-list: kg%irrk
   !  tetra(2, ic, it) -> index of   reducible k-list: kg%kpt
   call boltz_grid_setup(kg)

   !output to hdf5 file
   call output_grid(prefix, kg, emin, emax, bmin, bmax)
   !some message to stdout
   write(stdout,'(/5x,a,i11 )') 'number of  tetrahedra selected:', kg%num_tet
   write(stdout,'( 5x,a,i11 )') 'number of   reducible k-points:', kg%nk
   write(stdout,'( 5x,a,i11/)') 'number of irreducible k-points:', kg%nk_irr
endif
   !release space
   deallocate(irrk, irrk_e)

   call mp_bcast(kg%num_tet, ionode_id, inter_pool_comm)   
   call mp_bcast(kg%nk, ionode_id, inter_pool_comm)   
   call mp_bcast(kg%nk_irr, ionode_id, inter_pool_comm)
   !
   if( .not. associated(kg%tetra) ) allocate( kg%tetra(2, 4, kg%num_tet) )
   if( .not. associated(kg%irrk) ) allocate( kg%irrk( kg%nk_irr ) )
   call mp_bcast(kg%tetra, ionode_id, inter_pool_comm)   
   call mp_bcast(kg%irrk, ionode_id, inter_pool_comm)
   !
   if( .not. associated(kg%kpt) ) allocate( kg%kpt( kg%nk ) )
   if( .not. associated(kg%kpt2ir) ) allocate( kg%kpt2ir( kg%nk ) )
   call mp_bcast(kg%kpt, ionode_id, inter_pool_comm)   
   call mp_bcast(kg%kpt2ir, ionode_id, inter_pool_comm)
   !
   call mp_barrier(inter_pool_comm)
   !
   call stop_clock('Grid_gen')
end subroutine boltz_grid_generate


subroutine count_tetra(kdim, kidx, num_irk, irrk_e, numt)
   implicit none
   !kgrid dimension: nk1 * nk2 * nk3
   integer, intent(in) :: kdim(3)
   ! if kidx(i) = i: the i-th kpts is irreducible kpt; 
   !      otherwise kidx(i) is the irreducible kpts equivalent to the i-th kpts
   integer, intent(in) :: kidx( 0:(kdim(1)*kdim(2)*kdim(3)-1) )
   !number of irreducible kpts
   integer, intent(in) :: num_irk
   !label the status of irreducible kpts:
   ! >= 0: no level fall into [emin, emax], the value is the highest band 
   !         that have energy lower than emin.
   ! = -2: have levels fall into [emin, emax]
   integer, intent(in) :: irrk_e(num_irk)
   !number of tetrahedra selected based on the energy of vertex
   integer, intent(out) :: numt
   !
   integer :: nt, nk1, nk2, nk3, i, j, k, corn(6,4), it, ic, ckin(4)
   
   nt = 0
   nk1 = kdim(1); nk2 = kdim(2); nk3 = kdim(3)
   !careful: i,j,k convention should be in consistent with num2kpt
   do i = 0, nk1-1
   do j = 0, nk2-1
   do k = 0, nk3-1
      call get_tetra(nk1, nk2, nk3, i, j, k, corn)
      do it = 1, 6
         do ic = 1, 4
            ! whether the corner have level inside emin, emax
            ckin(ic) = irrk_e( kidx( corn(it,ic) ) )
         enddo 
      ! select tetrahedra that contributes to dos in [emin, emax]. 
      ! if: 1. any vertex of the tetrahedron has energy fall into [emin, emax]
      !     2. eig(n,vertex) have energies cross the [emin, emax] (e.g. ene in 
      !     [emin, emax] can cut this tetrahedron, important for small interval). 
      ! then the tetrahedron and all its corners are selected.
         if(any(ckin(1:4) < -1) .or. any((ckin(1:4)-ckin(1)).ne.0)) then
            nt = nt + 1
         endif
      enddo 
   enddo; enddo; enddo
   numt = nt
end subroutine


subroutine collect_tetra(kdim, kidx, num_irk, irrk_e, numt, totk, irrk_t, tetra)
   implicit none
   integer, intent(in) :: kdim(3), num_irk, irrk_e(num_irk), numt
   integer, intent(inout) :: kidx( 0:(kdim(1)*kdim(2)*kdim(3)-1) ) 
   integer, intent(out) :: totk, irrk_t(num_irk), tetra(:,:,:)
   !
   integer :: i, j, k, it, ic, map_ir(4), irr_pk, all_pk
   integer :: corn(6,4), ckin(4), nt, nk1, nk2, nk3
   !init irrk_t
   irrk_t = 0
   tetra = 0

   nt = 0;  irr_pk = 0;  all_pk = 0
   nk1 = kdim(1); nk2 = kdim(2); nk3 = kdim(3)
   !careful: i,j,k convention should be in consistent with num2kpt
   do i = 0, nk1-1
   do j = 0, nk2-1
   do k = 0, nk3-1
      call get_tetra(nk1, nk2, nk3, i, j, k, corn)
      do it = 1, 6
         do ic = 1, 4
            ! map all k-list to irr_list, 1 ~ num_irk
            map_ir(ic) =  mod( kidx( corn(it,ic) ),  2*num_irk )
            ! whether the corner have level inside emin, emax
            ckin(ic) = irrk_e( map_ir(ic) )
         enddo 
         if(any(ckin(1:4) < -1) .or. any((ckin(1:4)-ckin(1)).ne.0)) then
            nt = nt + 1
            !sanity check
            if(nt > numt) call errore('collect_tetra','nt > numt', 1)
            do ic = 1, 4
               ! re-label the irreducible-k. skip the corner that already labeled.
               if( irrk_t( map_ir(ic) ) .eq. 0 ) then
                  irr_pk = irr_pk + 1
                  irrk_t( map_ir(ic) ) = irr_pk
               endif
               ! label the selected kpoints in full kgrid. includes reducible-k
               if( kidx( corn(it, ic) ) <= num_irk ) then
                  all_pk = all_pk + 1
                  ! if kidx(ik) > num_irk, then this k-points have already been counted.
                  kidx( corn(it, ic) ) = kidx( corn(it, ic) ) + 2*num_irk
               endif
            enddo
            !collect selected tetrahedron.
            do ic = 1, 4
               tetra(1,ic,nt) = irrk_t( map_ir(ic) )
               tetra(2,ic,nt) = corn(it,ic)
            enddo
         endif
      enddo 
      enddo; enddo; enddo
      totk = all_pk
end subroutine


subroutine get_tetra(nk1, nk2, nk3, i, j, k, tetra)
   implicit none
   integer, intent(in) :: nk1, nk2, nk3, i, j, k
   integer, intent(out) :: tetra(6,4)
   !
   integer :: ip1, jp1, kp1, n1, n2, n3, n4, n5, n6, n7, n8
   
   !careful: i,j,k convention should be in consistent with num2kpt
   ! n1-n8 are the indices of k-point 1-8 forming a cube
   ip1 = mod( i+1, nk1)
   jp1 = mod( j+1, nk2)
   kp1 = mod( k+1, nk3)
   n1 =   k +   j*nk3 +   i*nk2*nk3
   n2 =   k +   j*nk3 + ip1*nk2*nk3
   n3 =   k + jp1*nk3 +   i*nk2*nk3
   n4 =   k + jp1*nk3 + ip1*nk2*nk3
   n5 = kp1 +   j*nk3 +   i*nk2*nk3
   n6 = kp1 +   j*nk3 + ip1*nk2*nk3
   n7 = kp1 + jp1*nk3 +   i*nk2*nk3
   n8 = kp1 + jp1*nk3 + ip1*nk2*nk3
   !along 1-8, this is a better choice when i, j, k are symmetric
   ! 2, 4, 1, 8;
   tetra(1,1) = n2; tetra(1,2) = n4; tetra(1,3) = n1; tetra(1,4) = n8
   ! 4, 1, 3, 8;
   tetra(2,1) = n4; tetra(2,2) = n1; tetra(2,3) = n3; tetra(2,4) = n8
   ! 2, 1, 6, 8;
   tetra(3,1) = n2; tetra(3,2) = n1; tetra(3,3) = n6; tetra(3,4) = n8
   ! 1, 3, 8, 7;
   tetra(4,1) = n1; tetra(4,2) = n3; tetra(4,3) = n8; tetra(4,4) = n7
   ! 1, 8, 5, 7;
   tetra(5,1) = n1; tetra(5,2) = n8; tetra(5,3) = n5; tetra(5,4) = n7
   ! 1, 6, 8, 5;
   tetra(6,1) = n1; tetra(6,2) = n6; tetra(6,3) = n8; tetra(6,4) = n5
end subroutine get_tetra


subroutine output_grid(prefix, kg, emin, emax, bmin, bmax)
   use pert_const, only: dp
   use boltz_utils,only: num2kpt
   use pert_utils, only: find_free_unit
   use hdf5_utils
   implicit none
   character(len=80) :: prefix
   type(grid), intent(in) :: kg
   integer, intent(in) :: bmin, bmax
   real(dp), intent(in) :: emin, emax
   !
   integer :: iunit, ir, i, brange(2)
   real(dp) :: xkg(3), erange(2)
   integer(HID_T) :: file_id
   real(dp), allocatable :: kpts(:,:)
   
   !output irreducible k-list
   iunit = find_free_unit()
   open(iunit, file=trim(prefix)//'_tet.kpt',status='replace',form='formatted')
   write(iunit,'(1x,i8,a, 3(2x,i4), 2x,a,i8, 2x,a,i8)') kg%nk_irr, " crystal", &
      (kg%ndim(i), i=1,3), ' #.tetra ', kg%num_tet, ' #.tot.k ', kg%nk
   do ir = 1, kg%nk_irr
      xkg(1:3) = num2kpt( kg%irrk(ir), kg%ndim)
      write(iunit,'(3(f12.8,2x), f16.12)') xkg(1:3), 1.d0/real(kg%nk_irr, dp)
   enddo
   close(iunit)

   allocate( kpts(3, kg%nk) )
   do ir = 1, kg%nk
      kpts(:,ir) = num2kpt(kg%kpt(ir), kg%ndim)
   enddo
   
   erange = (/emin, emax/)
   brange = (/bmin, bmax/)
   !ouput tetra 
   call hdf_open_file(file_id, trim(prefix)//'_tet.h5', status='NEW')
   !
   call hdf_write_dataset(file_id, 'kgrid_dim', kg%ndim(1:3))
   call hdf_write_dataset(file_id, 'num_kpts',  kg%nk)
   call hdf_write_dataset(file_id, 'num_irr_kpts', kg%nk_irr)
   call hdf_write_dataset(file_id, 'num_tetra', kg%num_tet)
   call hdf_write_dataset(file_id, 'band_window', brange)
   call hdf_write_dataset(file_id, 'energy_window', erange)
   !
   call hdf_write_dataset(file_id, 'kpts_irr',  kg%irrk( 1:kg%nk_irr ) )
   call hdf_write_dataset(file_id, 'tetra', kg%tetra(1:2, 1:4, 1:kg%num_tet) )
   call hdf_write_dataset(file_id, 'kpts_all',  kg%kpt( 1:kg%nk ) )
   call hdf_write_dataset(file_id, 'kpt2ir', kg%kpt2ir(1:kg%nk ) )
   
   call hdf_write_dataset(file_id, 'kpts_all_crys_coord', kpts)
   call hdf_close_file(file_id)

   deallocate( kpts )
end subroutine output_grid

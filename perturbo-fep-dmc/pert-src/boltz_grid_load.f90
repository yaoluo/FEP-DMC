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
! Read input files: tet.kpt, tet.tetra
subroutine boltz_grid_load(kg)
   use pert_const, only: dp
   use pert_utils, only: find_free_unit
   use boltz_utils, only: kpt2num
   use pert_param, only: prefix, boltz_kdim
   use qe_mpi_mod, only: stdout, ionode, ionode_id, mp_bcast, inter_pool_comm
   use hdf5_utils
   implicit none
   type(grid), intent(inout) :: kg
   ! local variables
   integer :: kunit, tunit, ios, ik, i, irec, ikpt !, i, ib, ierr, sunit
   integer :: ir_nk, itmp(3), num_t, tot_k, ttmp(8), nkr, tmp !, nept, 
   character(len=120) :: fname, fname1, ctmp1, ctmp2, ctmp3
   real(dp) :: ktmp(3) !, wtmp, rtmp
   logical :: has_file
   integer(HID_T) :: file_id
   
   !sanity check
   if( associated(kg%tetra) .and. associated(kg%irrk) ) then
      write(stdout, '(5x,a)') "Warn: kgrid is already loaded, skip the duplicate loading."
      return
   endif

if (ionode) then
   write(stdout,'(5x,a)') ">loading kpt, tetra from " // trim(prefix) // '_tet.h5'
   !
   fname = trim(prefix)//"_tet.h5"
   inquire(file=trim(fname), exist=has_file)
   if(.not. has_file) call errore('boltz_grid_load','missing '// trim(fname), 1)
   !
   call hdf_open_file(file_id, trim(fname), status='OLD', action='READ')
   call hdf_read_dataset(file_id, 'kgrid_dim', itmp(1:3))
   call hdf_read_dataset(file_id, 'num_kpts',  tot_k)
   call hdf_read_dataset(file_id, 'num_irr_kpts', ir_nk)
   call hdf_read_dataset(file_id, 'num_tetra', num_t)
   
   if(any( (boltz_kdim(1:3)-itmp(1:3)) .ne. 0) ) &
      call errore('boltz_grid_load','boltz_kdim does not match '// trim(fname),1)
   
   kg%nk  = tot_k
   kg%nk_irr  = ir_nk
   kg%num_tet = num_t
   kg%ndim(1:3) = boltz_kdim(1:3)
   allocate( kg%irrk( kg%nk_irr ) )
   allocate( kg%tetra(2, 4, kg%num_tet) )
   allocate( kg%kpt( kg%nk ) )
   allocate( kg%kpt2ir( kg%nk ) )
   !
   call hdf_read_dataset(file_id, 'kpts_irr', kg%irrk)
   call hdf_read_dataset(file_id, 'tetra', kg%tetra)
   call hdf_read_dataset(file_id, 'kpts_all', kg%kpt)
   call hdf_read_dataset(file_id, 'kpt2ir', kg%kpt2ir)
   call hdf_close_file(file_id)

   !if prefix_tet.kpt exist, check consistency
   fname1 = trim(prefix)//"_tet.kpt"
   inquire(file=trim(fname1), exist=has_file)
   if(has_file) then
      !read irreducible k-points  
      kunit = find_free_unit()
      open(kunit,file=trim(fname1),status='old',form='formatted',err=100,iostat=ios)
      read(kunit, *) ir_nk, ctmp1, itmp(1:3), ctmp2, num_t, ctmp3, tot_k

      if((kg%nk_irr .ne. ir_nk).or.(kg%num_tet .ne. num_t).or.(kg%nk .ne. tot_k) &
         .or. any(boltz_kdim-itmp .ne. 0)) &
      call errore('boltz_grid_load', 'inconsistent ' // trim(fname) // ', ' // trim(fname1), 1)

      do ik = 1, kg%nk_irr
         read(kunit,*) ktmp(1:3)
         ikpt = kpt2num( ktmp(1:3), boltz_kdim )
         if( kg%irrk(ik) .ne. ikpt ) call errore('boltz_grid_load', &
            'inconsistent kpoints in ' // trim(fname) // ' and ' // trim(fname1), 1)
      enddo
      
      close(kunit)
   endif
   
   !some message to stdout
   write(stdout,'(/5x,a,i11 )') 'number of  tetrahedra selected:', kg%num_tet
   write(stdout,'( 5x,a,i11 )') 'number of   reducible k-points:', kg%nk
   write(stdout,'( 5x,a,i11/)') 'number of irreducible k-points:', kg%nk_irr
endif

   call mp_bcast( kg%nk,  ionode_id, inter_pool_comm)
   call mp_bcast( kg%nk_irr,  ionode_id, inter_pool_comm)
   call mp_bcast( kg%num_tet, ionode_id, inter_pool_comm)
   call mp_bcast( kg%ndim, ionode_id, inter_pool_comm)
   !
   if(.not. associated(kg%irrk) )  allocate( kg%irrk( kg%nk_irr ) )
   if(.not. associated(kg%tetra))  allocate( kg%tetra(2, 4, kg%num_tet) )
   call mp_bcast(kg%irrk, ionode_id, inter_pool_comm)
   call mp_bcast(kg%tetra, ionode_id, inter_pool_comm)
   !
   if( .not. associated(kg%kpt) ) allocate( kg%kpt( kg%nk ) )
   if( .not. associated(kg%kpt2ir) ) allocate( kg%kpt2ir( kg%nk ) )
   call mp_bcast(kg%kpt, ionode_id, inter_pool_comm)   
   call mp_bcast(kg%kpt2ir, ionode_id, inter_pool_comm)

   !weight of the k-points
   kg%kweight = 1.0_dp / dble(kg%ndim(1)) / dble(kg%ndim(2)) / dble(kg%ndim(3))
   kg%tweight = kg%kweight / 6.0_dp
   
   return
100 call errore('boltz_grid_load','opening file '//trim(fname1),abs(ios))
end subroutine boltz_grid_load


!! read irreducible k-list and tetra, and recover the full grid
subroutine boltz_grid_setup(kg)
   implicit none
   type(grid), intent(inout) :: kg
   ! local variable
   integer, parameter :: nsym = 48 ! max number of symmetry operations
   integer :: it, ic, ir, n, i, ik, ierr
   integer, allocatable :: fktmp(:,:), kcol(:)
   
   !initial work array
   allocate( fktmp(nsym, kg%nk_irr), kcol(kg%nk),  stat=ierr)
   if(ierr /= 0) call errore('boltz_grid_setup','allocating space failed',1)

   fktmp(:,:) = -1
   ! collect kpoint corresponding to the same irreducible k into fktmp(:,ir)
   do it = 1, kg%num_tet
   do ic = 1, 4
      ! ir: index of irreducible-k list; n: kpt coordinate, num2kpt(n, kdim)
      ir = kg%tetra(1, ic, it)
      n  = kg%tetra(2, ic ,it)
      ! scan all the element of fktmp(:,ir)
      do i = 1, nsym
         ! fktmp(:,ir) store all the possible kpoints (n) appears in tetra
         ! that equivalent to irreducible k-list ir.
         ! fktmp(i,ir) == 0: means this position is empty.
         ! scan all the non-empty elements in fktmp(:,ir), if found the same n, 
         ! which means the k-points already appears before, then move to next. 
         if( n .eq. fktmp(i, ir) ) then
            kg%tetra(2, ic, it) = i
            exit
         ! if kpoint n never appears before, add to empty position of fktmp(:,ir)
         elseif(fktmp(i,ir) < 0) then
            fktmp(i,ir) = n
            kg%tetra(2,ic,it) = i
            exit
         ! if fktmp(i,ir) .ne. -1 and .ne. n, then continue.
         ! but after i reach nsym, fktmp still .ne. -1 and .ne. n, raise error.
         elseif(i .eq. nsym) then
            call errore('boltz_grid_setup','error in setup full k-grid', 1)
         endif
      enddo
   enddo; enddo
   
   !check and collect all the kpt in kcol
   kcol(:) = -1;  ik = 0
   do ir = 1, kg%nk_irr
   do i = 1, nsym
      if(fktmp(i,ir) < 0) cycle
      ik = ik + 1
      ! sanity check
      if(ik > kg%nk) call errore('boltz_grid_setup',"ik > kg%nk",1)
      kcol(ik) = fktmp(i, ir)
   enddo; enddo
   ! sanity check
   if(ik .ne. kg%nk) call errore('boltz_grid_setup', "ik .ne. kg%nk.",1)
   
   !sort and store in kg%kpt
   !initialize kg%kpt with kpoints stored in fktmp in ascending order
   kg%kpt2ir(:) = 0;  kg%kpt(:) = -1
   do ir = 1, kg%nk_irr
   do i = 1, nsym
      if(fktmp(i,ir) < 0) cycle
      ! find position of this kpt in the ordered kgrid.
      ic = 1
      do ik = 1, kg%nk
         if( fktmp(i,ir) > kcol(ik) ) ic = ic + 1
      enddo
      ! if this kpt is larger then n kpts, its position is n+1.
      kg%kpt( ic ) = fktmp(i, ir)
      ! this k-points equivalent to irreducible kpoints-ir
      kg%kpt2ir( ic ) = ir
      ! reassign fktmp(i,ir) to its position in kgrid
      fktmp(i,ir) = ic
   enddo; enddo
   !debug & test
   if( any(kg%kpt(:) < 0) .or. any(kg%kpt2ir(:) .eq. 0) ) &
      call errore('boltz_grid_setup','kgrid sorting failed.',1)

   ! update tetra, afterward:
   !  tetra(1, ic, it) -> index of irreducible k-list: kg%irrk
   !  tetra(2, ic, it) -> index of   reducible k-list: kg%kpt
   do it = 1, kg%num_tet
   do ic = 1, 4
      ir = kg%tetra(1, ic, it)
      i  = kg%tetra(2, ic ,it)
      kg%tetra(2, ic, it) = fktmp(i, ir)
   enddo; enddo
   !release memory
   deallocate(fktmp, kcol)
end subroutine boltz_grid_setup

!first principle diagmc is very expansive because of the interpolation of gkq. ryd2ev
!Our solution is to tabulate the quantity in a fine mesh before actual calculation starts. 
!We tabulate : 
!  1. electron : ek, uk 
!  2. phonon : wq, uq
!  3. e-ph : gq, gkq using svd technical
! interpolation of electron and electron phonon matrix 
module HamiltonianTable 
   use pert_const, only : dp,twopi
   use pert_utils, only: abs2,get_exp_ikr,bose, find_free_unit, fermi
   use vector_list,only: vlist, rand_vector_list,load_vector_list
   use pert_data,  only: mass, epwan_fid, kc_dim, qc_dim, alat
   use pert_param, only: fqlist, fklist, prefix, band_min, band_max, phfreq_cutoff, ntemper, temper,&
            cauchy_scale, sampling, nsamples,Ryd2eV,Ryd2meV,svd_dir,read_svd,read_svd,Nsvd,nk_svd,nq_svd,onRAM,read_H, &
            frohlich_only, g_offset 
   use qe_mpi_mod, only: ionode, stdout, mp_sum, inter_pool_comm, npool, &
            mp_split_pools, distribute_points,ionode_id, mp_bcast, my_pool_id,mp_barrier,meta_ionode_id,world_comm
   use band_structure, only: electron_wann, init_electron_wann, solve_eigenvalue_vector,solve_band_velocity
   use phonon_dispersion, only: lattice_ifc, init_lattice_ifc, solve_phonon_modes
   use elphon_coupling_matrix, only: eph_transform_fast,elph_mat_wann,eph_wan_longrange
   use eph_svd, only : nbnd,nmod
   
   logical, parameter :: frohlich = .true.
   logical :: complete = .false.
   real(dp),allocatable :: kgrid(:,:)
   integer, allocatable :: ikgrid(:,:)
   integer :: Nkxyz(3), Nktot
   integer, allocatable :: IndexK(:,:,:)   !(nx,ny,nz)
   !electron 
   real(dp),ALLOCATABLE :: enk_grid(:,:), vnk_grid(:,:,:)  
   complex(dp),allocatable :: uk_grid(:,:,:), Hk_grid(:,:,:)
   !phonon 
   real(dp),ALLOCATABLE :: wq_grid(:,:)
   complex(dp),allocatable :: uq_grid(:,:,:), Dq_grid(:,:,:)  
   !gkq 
   complex(dp),allocatable :: Uk_gf(:,:,:,:)   ! (nb, nb, nsvd, nk)
   complex(dp),allocatable :: Vq_gf(:,:,:,:,:) ! (nb, nb, nmod, nsvd, nk)
   complex(dp), allocatable :: gq_f(:,:,:,:)   ! (nb, nb, nmod, nk )
   complex(dp), allocatable :: gq_sample_f(:,:,:,:)   ! (nb, nb, nmod, nk ), q-dependent part : \sum_{Rp} g(Re=0,Rp) e^{iqRp} 
   !for importance sampling of phonon momentum 
   real(dp), allocatable :: grid_weight(:)


   contains 
   
   !cubic mesh, this 3D grid is for both k & q 
   subroutine setup_kgrid()
      implicit none 
      integer :: ix,iy,iz, icount   
      real(dp) :: kp(3)
      character(len=80) :: fqlist_frohlich
      Nkxyz = nk_svd
      Nktot = Nkxyz(1) * Nkxyz(2) * Nkxyz(3)

      if(ionode)write(stdout,'(A20,6i5)')'k/q mesh = ',Nkxyz
      if(.not.allocated(kgrid)) then 
         ALLOCATE( kgrid(3,Nktot) )
         ALLOCATE( ikgrid(3,Nktot) )
         ALLOCATE( IndexK(Nkxyz(1),Nkxyz(2),Nkxyz(3)))
      endif 

      icount = 0 
      do ix = 1,Nkxyz(1); do iy = 1,Nkxyz(2); do iz = 1,Nkxyz(3)
         icount = icount + 1 
         kp = (/ ix-1,iy-1,iz-1 /)  ; kp = kp / Nkxyz
         ikgrid(:,icount) = (/ ix-1, iy-1, iz-1 /) !k1+k2 -> k3, and ik3 index can be found is IndexK, k3 found in kgrid
         kgrid(:,icount) = kp 
         IndexK(ix,iy,iz) = icount 
      enddo; enddo; enddo 
   end subroutine

   !---------------------------------------------------!
   !------- wrapper for read out e-ph Hamltonian ------!
   !---------------------------------------------------!

   !momentum in integer number  
   subroutine solve_electron_fast(i_kpt,enk,uk,vk)
      implicit none 
      integer, INTENT(INOUT) :: i_kpt(3) 
      real(dp), INTENT(OUT) :: enk(:)
      complex(dp), INTENT(OUT), optional :: uk(:,:)
      real(dp), INTENT(OUT), optional :: vk(:,:)
      
      integer :: ik, it 

      ik = find_ik(i_kpt)
      enk = enk_grid(:,ik)
      if(present(uk))uk = uk_grid(:,:,ik)
      if(present(vk))vk = vnk_grid(:,:,ik)
   end subroutine

   subroutine solve_electron_velocity(i_kpt,vk)
      implicit none 
      integer, INTENT(INOUT) :: i_kpt(3) 
      real(dp), INTENT(OUT) :: vk(:,:)      
      integer :: ik, it 

      ik = find_ik(i_kpt)
      vk = vnk_grid(:,:,ik)
      vk = vk * ryd2ev * alat

   end subroutine

   subroutine solve_phonon_fast(i_qpt,wq,uq)
      implicit none 
      integer, INTENT(INOUT) :: i_qpt(3) 
      real(dp), INTENT(OUT) :: wq(:)
      complex(dp), INTENT(OUT), optional :: uq(:,:)
      integer :: iq, it 

      iq = find_ik(i_qpt)
      wq = wq_grid(:,iq)
      if(present(uq)) then 
         stop 'Hamiltonian table does not keep phonon vector to save memory'
         uq = uq_grid(:,:,iq)
      endif 
   end subroutine

   subroutine solve_gkq_fast(i_kpt, i_qpt, gkq )
      !wannier gauge still
      use elphon_coupling_matrix, only : elph_mat_wann,eph_transform
      implicit none 
      integer, INTENT(INOUT) :: i_kpt(3), i_qpt(3) 
      complex(dp),INTENT(OUT) :: gkq(:,:,:)
      !work 
      integer :: ik,ikq,iq,inq,ig, im,isvd, ikqpt(3), inqpt(3), ib,jb
      real(dp) :: kqt(3), kcop(3)
      complex(dp) :: gkqc(nbnd,nbnd,nmod)
      complex(dp), dimension(nbnd,nbnd,nmod,nsvd) ::vf1, vf2
      complex(dp), dimension(nbnd,nbnd,nsvd) :: uf1, uf2 
      
     
      iq = find_ik(i_qpt)
      inqpt = -i_qpt 
      inq = find_ik(inqpt)
      
      if(frohlich_only) then 
         gkq = gq_f(:,:,:,iq)
      else 
         gkq = gq_sample_f(:,:,:,iq) + gq_f(:,:,:,iq) 
      endif 
      
      do im = 1,nmod 
         if(wq_grid(im, iq) > phfreq_cutoff) then  
            gkq(:,:,im) = gkq(:,:,im) / sqrt(2.d0*wq_grid(im, iq))
            do ib = 1,nbnd 
               gkq(ib,ib,im) = gkq(ib,ib,im) + g_offset
            enddo 
         else 
            gkq(:,:,im) = 0.d0 
         endif 
      enddo 
      return
   end subroutine

   subroutine solve_gkq_full_fast(i_kpt, i_qpt, gkq )
      !wannier gauge still
      use elphon_coupling_matrix, only : elph_mat_wann,eph_transform
      implicit none 
      integer, INTENT(INOUT) :: i_kpt(3), i_qpt(3) 
      complex(dp),INTENT(OUT) :: gkq(:,:,:)
      !work 
      integer :: ik,ikq,iq,inq,ig, im,isvd, ikqpt(3), inqpt(3), ib,jb
      real(dp) :: kqt(3), kcop(3)
      complex(dp) :: gkqc(nbnd,nbnd,nmod)
      complex(dp), dimension(nbnd,nbnd,nmod,nsvd) ::vf1, vf2
      complex(dp), dimension(nbnd,nbnd,nsvd) :: uf1, uf2 
      
      !gkq = 0.d0; gkqc = 0.d0 
      ik = find_ik(i_kpt)
      iq = find_ik(i_qpt)

      gkq = gq_f(:,:,:,iq)
      if(frohlich_only) then 
         do im = 1,nmod 
            if(wq_grid(im, iq) > phfreq_cutoff) then  
               gkq(:,:,im) = gkq(:,:,im) / sqrt(2.d0*wq_grid(im, iq))
            else 
               gkq(:,:,im) = 0.d0 
            endif 
         enddo 
         return 
      endif  

      do isvd = 1,nsvd 
         do im = 1,nmod
            gkq(:,:,im) = gkq(:,:,im) + Uk_gf(:,:,isvd,ik) * Vq_gf(:,:,im,isvd,iq)
         enddo 
      enddo
    
      ikqpt = i_kpt+i_qpt
      ikq = find_ik(ikqpt)
      inqpt = -i_qpt
      inq = find_ik(inqpt)
      gkqc = gq_f(:,:,:,inq)
      do isvd = 1,nsvd 
         do im = 1,nmod
            gkqc(:,:,im) = gkqc(:,:,im) + Uk_gf(:,:,isvd,ikq) * Vq_gf(:,:,im,isvd,inq)
         enddo 
      enddo

      do im = 1,nmod 
         if(wq_grid(im, iq) > phfreq_cutoff) then  
            gkq(:,:,im) = 0.5d0 * (gkq(:,:,im) + transpose(dconjg(gkqc(:,:,im)))) / sqrt(2.d0*wq_grid(im, iq))
         else 
            gkq(:,:,im) = 0.d0 
         endif 
      enddo 

       
   end subroutine

   !momentum in real number 
   subroutine solve_electron_fast_double(kpt,enk,uk)
      implicit none 
      real(dp), INTENT(INOUT) :: kpt(3) 
      real(dp), INTENT(OUT) :: enk(:)
      complex(dp), INTENT(OUT), optional :: uk(:,:)
      integer :: ikpt(3)
      ikpt = int(kpt * Nkxyz)
      if(present(uk)) then 
         call solve_electron_fast(ikpt, enk, uk)
      else 
         call solve_electron_fast(ikpt, enk)
      endif 
   end subroutine

   subroutine solve_phonon_fast_double(qpt,wq,uq)
      implicit none 
      real(dp), INTENT(INOUT) :: qpt(3) 
      integer :: i_qpt(3) 
      real(dp), INTENT(OUT) :: wq(:)
      complex(dp), INTENT(OUT), optional :: uq(:,:)
      integer :: iq, it 
      
      i_qpt = int(qpt*Nkxyz) 
      if(present(uq)) then 
         call solve_phonon_fast(i_qpt, wq, uq)
      else 
         call solve_phonon_fast(i_qpt, wq)
      endif 
   end subroutine

   subroutine solve_gkq_fast_double(kpt, qpt, gkq)
      !wannier gauge still
      use elphon_coupling_matrix, only : elph_mat_wann,eph_transform
      implicit none 
      real(dp), INTENT(IN) :: kpt(3), qpt(3) 
      complex(dp),INTENT(OUT) :: gkq(:,:,:)
      !work 
      integer :: ikpt(3), iqpt(3) 

      ikpt = int(Nkxyz * kpt)
      iqpt = int(Nkxyz * qpt)
      !call solve_gkq_fast(ikpt, iqpt, gkq)
      call solve_gkq_full_fast(ikpt, iqpt, gkq)
   end subroutine

   subroutine Read_gkq_bebug(ep,i_kpt,i_qpt,gkq)
      use elphon_coupling_matrix, only : elph_mat_wann,eph_transform
      implicit none 
      type(elph_mat_wann),INTENT(IN) :: ep
      integer,INTENT(INOUT) :: i_kpt(3),i_qpt(3) 
      complex(dp),INTENT(OUT) :: gkq(:,:,:)
      !work 
      integer :: ik,iq,ig 
      ik = find_ik(i_kpt)
      iq = find_ik(i_qpt)

      gkq = 0.d0 
      return 
   end subroutine

   subroutine random_qfromlist_omp(seed, xqt, Pq, wq, uq, g_pol, vformf,vformfc)
      use random_tool
      implicit none 
      TYPE (VSL_STREAM_STATE) :: seed
      real(dp) :: xqt(3),Pq, ran
      integer :: iq 
      real(dp) :: wq(nmod)
      complex(dp) :: uq(nmod, nmod),g_pol(nmod)
      complex(dp),dimension(nbnd, nbnd, nmod, nsvd),optional :: vformf, vformfc 

   end subroutine

   !---------------------------------------------------!
   !--------------- tabulate Hamiltonian --------------!
   !---------------------------------------------------!
   subroutine Tabulate_for_DMC(el,ph,ep)
      use OMP_LIB
      use eph_svd, only : cal_vq_formf, cal_uk_formf, cal_gkq_Re0
      use Lapack_LY, only: DHEIGEN
      !use eph_holstein, only : cal_gkq_holstein
      implicit none 
      type(lattice_ifc),INTENT(IN) :: ph       !phonon dispersion 
      type(electron_wann),INTENT(IN) :: el     !el dispersion 
      type(elph_mat_wann),INTENT(IN) :: ep
      real :: max_mem = 50.d0 !GB
      real(dp) :: kp(3), qp(3), nqp(3)
      integer :: ik,iq,inq,iunit, ikpt(3)  
      CHARACTER(len=50) :: kstr,qstr  
      integer :: ig, im, ib,jb,ijb,ire,irp,icount,isvd,ika,ia  
      real(dp) :: mem_formf, abs_evq(ph%nm), Dq_re(ph%nm, ph%nm)
      integer :: kst_pool,kend_pool,qst_pool,qend_pool, i_openmp,nk_pool,nq_pool, iwork 
      complex(dp) :: mq(ph%nm,ph%nm), mnq(ph%nm,ph%nm), gkq_q(el%nb,el%nb,ph%nm),ctmp(ph%nm), uq(ph%nm,ph%nm), gpol(ph%nm)
      complex(dp) :: gtmp(el%nb,el%nb), phasefactor
      character(len=80) :: fqlist_frohlich
      logical, allocatable :: DoneList(:)
      integer :: iqpt(3), inqpt(3)

      call setup_kgrid()
      !memory check 
      max_mem = 400.d0 
      nbnd = el%nb ; nmod = ph%nm
      mem_formf = 16.d0*(nbnd**2+1)*(nmod+1)*(nsvd+1)*Nktot/1024/1024/1024  
      if(ionode) write(stdout,'(A30,f15.5,A3)')'memory for H table = ',mem_formf,'GB'
      if(mem_formf > max_mem) stop 'not enough memory; max_mem = 30 GB'



      if(.not.allocated(enk_grid)) then 
         allocate( enk_grid(el%nb, Nktot), uk_grid(el%nb, el%nb, Nktot), vnk_grid(3, el%nb, Nktot), Hk_grid(el%nb,el%nb,Nktot) )
         allocate( gq_f(nbnd,nbnd,ph%nm,Nktot), wq_grid(ph%nm,Nktot), gq_sample_f(nbnd,nbnd,ph%nm,Nktot) )
         allocate( Uk_gf(nbnd,nbnd,nsvd,Nktot), Vq_gf(nbnd,nbnd,nmod,nsvd,Nktot), grid_weight(Nktot) )
      endif 
      

      if(read_H) then 
         if(ionode) then
            write(stdout,'(A80)')' read H table from '//trim(adjustl(svd_dir))
            iunit = find_free_unit() 
            open(iunit,file=trim(adjustl(svd_dir))//'/enk.dat',action='read', form='Unformatted' )
            read(iunit) enk_grid
            close(iunit)
            open(iunit,file=trim(adjustl(svd_dir))//'/vnk.dat',action='read', form='Unformatted' )
            read(iunit) vnk_grid
            close(iunit)
            open(iunit,file=trim(adjustl(svd_dir))//'/uk.dat',action='read', form='Unformatted' )
            read(iunit) uk_grid
            close(iunit)
            open(iunit,file=trim(adjustl(svd_dir))//'/Hk.dat',action='read', form='Unformatted' )
            read(iunit) Hk_grid
            close(iunit)
            open(iunit,file=trim(adjustl(svd_dir))//'/wq.dat',action='read', form='Unformatted' )
            read(iunit) wq_grid
            close(iunit)
            open(iunit,file=trim(adjustl(svd_dir))//'/gq_f.dat',action='read', form='Unformatted' )
            read(iunit) gq_f
            close(iunit)
            open(iunit,file=trim(adjustl(svd_dir))//'/gq_sample_f.dat',action='read', form='Unformatted' )
            read(iunit) gq_sample_f
            close(iunit)
            open(iunit,file=trim(adjustl(svd_dir))//'/Uk_gf.dat',action='read', form='Unformatted' )
            read(iunit) Uk_gf
            close(iunit)
            open(iunit,file=trim(adjustl(svd_dir))//'/Vq_gf.dat',action='read', form='Unformatted' )
            read(iunit) Vq_gf
            close(iunit)
         endif
         call mp_bcast(enk_grid, meta_ionode_id, world_comm)
         call mp_bcast(uk_grid, meta_ionode_id, world_comm)
         call mp_bcast(vnk_grid, meta_ionode_id, world_comm)
         call mp_bcast(Hk_grid, meta_ionode_id, world_comm)
         call mp_bcast(wq_grid, meta_ionode_id, world_comm)
         call mp_bcast(gq_f, meta_ionode_id, world_comm)
         call mp_bcast(Uk_gf, meta_ionode_id, world_comm)
         call mp_bcast(Vq_gf, meta_ionode_id, world_comm)

         !convert to eV unit 

         complete = .true.

         write(*,'(A40,3E20.10)')'max(gq, Vq, gq_sample_f)',maxval(abs(gq_f)),maxval(abs(Vq_gf)),maxval(abs(gq_sample_f))
         return 
      endif 

      enk_grid = 0.d0; uk_grid = 0.d0; vnk_grid = 0.d0; Hk_grid= 0.d0  
      wq_grid = 0.d0;  gq_f = 0.d0; Uk_gf = 0.d0; Vq_gf = 0.d0 ; grid_weight = 0.d0  
      !mpi split of Nk point 
      call mp_split_pools(Nktot, kst_pool, kend_pool, nk_pool)
      iunit = find_free_unit()

      if(ionode)  write(*,'(A40)')'------------------------------------------'
      if(ionode)  write(*,'(A15)')'progress_k = '
      icount = 1
      !$omp parallel default(shared) private(ik,kp,i_openmp,ib)
      i_openmp = OMP_GET_THREAD_NUM()
      !$omp do schedule(guided)
      do ik = kst_pool,kend_pool
         kp = kgrid(:,ik)
         
         ! 1.electronic structure 
         call solve_band_velocity(el, kp, vnk_grid(:,:,ik),enk_grid(:,ik),uk_grid(:,:,ik))
         !! Hk = Uk Ek Uk^\dagger
         do ib = 1, el%nb 
            Hk_grid(ib,ib,ik) = enk_grid(ib,ik)
         enddo 
         Hk_grid(:,:,ik) = matmul(uk_grid(:,:,ik), matmul(Hk_grid(:,:,ik) ,transpose(conjg(uk_grid(:,:,ik))) ))
         ! 2.formfactor 
         call cal_uk_formf( kp, Uk_gf(:,:,:,ik) ) 

         !track the progress 
         !$omp critical
         icount = icount + 1
         if(ionode .and. mod(icount,nk_pool/10).eq.0 ) then 
            write(*,'(A15,i5,A3)')' ',(icount)/(nk_pool/100),'%'
         end if 
         !$omp end critical
      enddo 
      !$omp end do
      !$omp end parallel 

      call mp_sum(enk_grid, inter_pool_comm)
      call mp_sum(uk_grid, inter_pool_comm)
      call mp_sum(vnk_grid, inter_pool_comm)
      call mp_sum(Hk_grid, inter_pool_comm)
      call mp_sum(Uk_gf, inter_pool_comm)
      if(ionode) then 
         open(iunit,file=trim(adjustl(svd_dir))//'/uk.dat',action='write', form='Unformatted' )
         write(iunit) uk_grid
         close(iunit)
         open(iunit,file=trim(adjustl(svd_dir))//'/enk.dat',action='write', form='Unformatted' )
         write(iunit) enk_grid
         close(iunit)
         open(iunit,file=trim(adjustl(svd_dir))//'/vnk.dat',action='write', form='Unformatted' )
         write(iunit) vnk_grid
         close(iunit)
         open(iunit,file=trim(adjustl(svd_dir))//'/Hk.dat',action='write', form='Unformatted' )
         write(iunit) Hk_grid
         close(iunit)
         open(iunit,file=trim(adjustl(svd_dir))//'/Uk_gf.dat',action='write', form='Unformatted' )
         write(iunit) Uk_gf
         close(iunit)

         !delete 
         deallocate(uk_grid,enk_grid,vnk_grid,Hk_grid,Uk_gf )
      endif



      call mp_split_pools(Nktot, qst_pool, qend_pool, nq_pool)
      if(ionode)  write(*,'(A40)')'------------------------------------------'
      if(ionode)  write(*,'(A25)')'progress frohlich q ='
      icount = 1 
      allocate(DoneList(Nktot))
      DoneList = .False. 
      !$omp parallel default(shared) private(iq,inq,iqpt,inqpt,im,ia,ib,jb,qp,nqp,mq,mnq,i_openmp,ctmp,gtmp,phasefactor,ika,Dq_re)
      i_openmp = OMP_GET_THREAD_NUM()
      !$omp do schedule(guided)
      do iq = qst_pool,qend_pool
         if(DoneList(iq)) cycle 
         qp = kgrid(:,iq)
         iqpt = ikgrid(:,iq)
         inqpt = -iqpt 
         inq = find_ik(inqpt)    
         !nqp = kgrid(:,inq)
         nqp = -qp
         DoneList(iq) = .True. 
         DoneList(inq) = .True. 

         ! 1. phonon structure 
         call solve_phonon_modes(ph,qp,wq_grid(:,iq),mq)
         if(iq.eq.inq) then 
            write(*,'(3i4)') iqpt
            ! try to make mq as real number, multiply as factor to cancel the phase  
            Dq_re = 0.d0 
            do im = 1,nmod 
               Dq_re(im, im) = wq_grid(im,iq)
            enddo 
            do ika = 1, nmod
               ia = (ika - 1) / 3 + 1
               mq(ika, :) = mq(ika, :) * sqrt( mass(ia) )
            enddo
            !check unitary 
            !write(*,'(12f10.5)')matmul(mq, conjg(transpose(mq)) )
            mq = matmul(mq, matmul(Dq_re, conjg(transpose(mq)) ))
            Dq_re = real(mq)
            call DHEIGEN(nmod, Dq_re, wq_grid(:,iq))
            mq = Dq_re
            do ika = 1, nmod
               ia = (ika - 1) / 3 + 1
               mq(ika, :) = mq(ika, :) / sqrt( mass(ia) )
            enddo
            write(*,'(6f10.5)') wq_grid(:,iq) * ryd2mev
             
         endif    
         wq_grid(:,inq) = wq_grid(:,iq); mnq = dconjg(mq)
         do ib = 1,nbnd;do jb = 1,nbnd 
            ctmp = gq_sample_f(ib, jb, :, iq)
            gq_sample_f(ib, jb, :, iq) = matmul(ctmp, mq)
         enddo;enddo

         ! 2. e-ph, Re = 0 component  
         call cal_gkq_Re0(kp, qp, gq_sample_f(:,:,:,iq))
         do ib = 1,nbnd;do jb = 1,nbnd 
            ctmp = gq_sample_f(ib, jb, :, iq)
            gq_sample_f(ib, jb, :, iq) = matmul(ctmp, mq)
         enddo;enddo
         if(iq.ne.inq) then 
            call cal_gkq_Re0(kp, nqp, gq_sample_f(:,:,:,inq))
            do ib = 1,nbnd;do jb = 1,nbnd 
               ctmp = gq_sample_f(ib, jb, :, inq)
               gq_sample_f(ib, jb, :, inq) = matmul(ctmp, mnq)
            enddo;enddo
         else 
            !there is a bug for q = -q, that 
            do im = 1,nmod
               gtmp = gq_sample_f(:, :, im, iq)  
               if(  maxval(abs(gtmp - transpose(conjg(gtmp )))) > 1e-10) then
                  write(*,'(i5,f15.5,E20.10)') im, wq_grid(im,iq)*ryd2mev, maxval(imag(mq))
                  write(*,'(2E20.10)') maxval(abs(gtmp - transpose(conjg(gtmp )))), maxval(abs(gtmp))
                  write(*,'(3i4)') iqpt
                  write(*,'(A60)') '-q=q with degeneracy will fail the debug checking!'
                  gq_sample_f(:, :, im, iq)   = 0.5_dp*( gtmp + transpose(conjg(gtmp)) )
               endif
            enddo 
         endif 

         ! 3. e-ph formfactor 
         call cal_vq_formf( qp,  Vq_gf(:,:,:,:,iq) )
         call cal_vq_formf( nqp, Vq_gf(:,:,:,:,inq) )
         !! transform to phonon-eigen basis 
         do ib = 1,nbnd;do jb = 1,nbnd; do isvd = 1,nsvd  
            ctmp = Vq_gf(ib, jb, :, isvd, iq)
            Vq_gf(ib, jb, :, isvd, iq) = matmul(ctmp, mq)
         enddo;enddo;enddo


         if(iq.ne.inq) then 
            do ib = 1,nbnd;do jb = 1,nbnd; do isvd = 1,nsvd  
               ctmp = Vq_gf(ib, jb, :, isvd, inq)
               Vq_gf(ib, jb, :, isvd, inq) = matmul(ctmp, mnq)
            enddo;enddo;enddo
         endif 
         
         ! 3. long-ranged e-ph  
         call eph_wan_longrange(ep%pol, qp, ctmp, mq)
         !! put to matrix form 
         do ib = 1,nbnd 
            do im =  1,nmod 
               gq_f(ib,ib,im,iq) = ctmp(im) 
            enddo 
         enddo 
         call eph_wan_longrange(ep%pol, nqp, ctmp, mnq)
         !! put to matrix form 
         do ib = 1,nbnd 
            do im =  1,nmod 
               gq_f(ib,ib,im,inq) = ctmp(im)
            enddo 
         enddo 


         ! set grid_weight 
         grid_weight(iq) = (maxval(abs(ctmp))+maxval(abs(Vq_gf(:, :, :, :, iq))))**2
         grid_weight(inq) = grid_weight(iq)
         !trace the progress 
         !$omp critical
         icount = icount + 1
         if( ionode .and. mod(icount,nq_pool/10).eq.0 ) then 
            write(*,'(A15,i5,A3)')' ', icount/(nq_pool/100),'%'
         end if 

         !$omp end critical
      enddo     
      !$omp end do
      !$omp end parallel   
      call mp_sum( wq_grid, inter_pool_comm )
      call mp_sum( gq_f, inter_pool_comm )
      call mp_sum( Vq_gf, inter_pool_comm )
      call mp_sum( grid_weight, inter_pool_comm )

      if(ionode)  write(*,'(A40)')'------------------------------------------'

      
      if(ionode) then 
         open(iunit,file=trim(adjustl(svd_dir))//'/wq.dat',action='write', form='Unformatted' )
         write(iunit) wq_grid
         close(iunit)
         open(iunit,file=trim(adjustl(svd_dir))//'/gq_f.dat',action='write', form='Unformatted' )
         write(iunit) gq_f
         close(iunit)
         open(iunit,file=trim(adjustl(svd_dir))//'/gq_sample_f.dat',action='write', form='Unformatted' )
         write(iunit) gq_sample_f
         close(iunit)
         open(iunit,file=trim(adjustl(svd_dir))//'/Vq_gf.dat',action='write', form='Unformatted' )
         write(iunit) Vq_gf
         close(iunit)
         deallocate(wq_grid,gq_f,Vq_gf,gq_sample_f)
         !grid write to text for constructing alias list 
         open(iunit,file=trim(adjustl(svd_dir))//'/qgrid-dmc.dat')
         do iq = 1,Nktot
            write(iunit,'(4E20.10)') kgrid(:, iq), grid_weight(iq)
         enddo 
         close(iunit)
      endif
   end subroutine

   function find_ik(ikpt)
      integer,intent(inout) :: ikpt(3)
      integer :: find_ik 
      integer :: it 

      do it = 1,3 
         ikpt(it) = mod(ikpt(it),Nkxyz(it))
         if(ikpt(it)<0) ikpt(it) = ikpt(it) + Nkxyz(it)
      enddo      
      find_ik = IndexK( ikpt(1)+1, ikpt(2)+1, ikpt(3)+1 )
      return 
   end function

end module

subroutine tabulate_H()
   use eph_svd, only : init_eph_svd
   use HamiltonianTable, only : read_H 
   use pert_data,  only: epwan_fid, qc_dim, kc_dim
   use band_structure, only: electron_wann, init_electron_wann
   use phonon_dispersion, only: lattice_ifc, init_lattice_ifc
   use elphon_coupling_matrix, only: elph_mat_wann, init_elph_mat_wann
   use HamiltonianTable, only : tabulate_for_DMC
   implicit none 
   type(lattice_ifc)   :: phon
   type(electron_wann) :: elec
   type(elph_mat_wann) :: elph

   call init_electron_wann(epwan_fid, kc_dim, elec)
   call init_lattice_ifc(epwan_fid, qc_dim, phon)
   call init_elph_mat_wann(epwan_fid, kc_dim, qc_dim, elph)
   
   !init svd 
   call init_eph_svd( elec, phon, elph )
   
   !tabulate 
   call Tabulate_for_DMC( elec, phon, elph )

end subroutine

subroutine cal_eph_svd()
   use eph_svd, only : init_eph_svd
   use pert_data,  only: epwan_fid, qc_dim, kc_dim
   use band_structure, only: electron_wann, init_electron_wann
   use phonon_dispersion, only: lattice_ifc, init_lattice_ifc
   use elphon_coupling_matrix, only: elph_mat_wann, init_elph_mat_wann
   implicit none 
   type(lattice_ifc)   :: phon
   type(electron_wann) :: elec
   type(elph_mat_wann) :: elph

   call init_electron_wann(epwan_fid, kc_dim, elec)
   call init_lattice_ifc(epwan_fid, qc_dim, phon)
   call init_elph_mat_wann(epwan_fid, kc_dim, qc_dim, elph)
   
   !init eph_svd 
   call init_eph_svd(elec,phon,elph)
end subroutine




!check svd grerp 
!for fast evaluation of e-ph matrix, apply svd on g(r_e,r_p)
module eph_svd
   use pert_const, only : dp,twopi
   use pert_utils, only: abs2,get_exp_ikr,bose, find_free_unit, fermi
   use vector_list,only: vlist, rand_vector_list,load_vector_list
   use wigner_seitz_cell, only: ws_cell
   use pert_data,  only: mass, epwan_fid, kc_dim, qc_dim
   use pert_param, only: fqlist, fklist, prefix, band_min, band_max, phfreq_cutoff, ntemper, temper,&
            cauchy_scale, sampling, nsamples,Ryd2eV,svd_dir,read_svd,read_formf,read_svd,Nsvd,nk_svd,nq_svd,onRAM 
   use qe_mpi_mod, only: ionode, stdout, mp_sum, inter_pool_comm, npool, &
            mp_split_pools, distribute_points,ionode_id, mp_bcast, my_pool_id,mp_barrier
   use band_structure, only: electron_wann, init_electron_wann, solve_eigenvalue_vector
   use phonon_dispersion, only: lattice_ifc, init_lattice_ifc, solve_phonon_modes
   use elphon_coupling_matrix, only: eph_wannier, elph_mat_wann, init_elph_mat_wann, &
      eph_fourier_elph, eph_transform,eph_transform_fast,copy_exp_ikr, eph_wan_longrange,eph_fourier_el,eph_fourier_el_para
   implicit none 
   integer :: nmod,nbnd, nre,nrp 


   integer :: NS_tot 
   complex(dp),ALLOCATABLE :: expq(:),expk(:)
   real(dp),allocatable :: S_g(:,:) ! (Nsvd,nmod,nbnd**2) 
   real(dp), allocatable :: Re(:,:), Rp(:,:)
   complex(dp),allocatable :: Ur_g(:,:,:),Vr_g(:,:,:,:) !(nRe,Nsvd,nbnd**2), (nmod,nRp,Nsvd,nbnd**2)    
   integer :: Nsvd_tot 
   complex(dp), allocatable :: g_Rp(:,:,:,:) !(nRp, nbnd, nbnd, nmod), for treat sign problem 
   
   contains 
   
   subroutine cal_gkq_svd(xkt,xqt,gkq)
      implicit none 
      real(dp),INTENT(IN) :: xkt(3), xqt(3)
      complex(dp),INTENT(out) :: gkq(nbnd,nbnd,nmod)
      integer :: ig, im, ib,jb,ijb,ire,irp 
      complex(dp) :: Uk_g(nbnd,nbnd,Nsvd),Vq_g(nbnd,nbnd,nmod,Nsvd)
      
      if(.not.allocated(expq)) then 
         ALLOCATE(expq(nrp),expk(nre))
      end if 
      call cal_uk_formf(xkt,Uk_g)
      call cal_vq_formf(xqt,Vq_g)
      gkq = 0.d0 
      do ig = 1,nsvd    
         do im = 1,nmod 
            gkq(:,:,im) = gkq(:,:,im) + Uk_g(:,:,ig) * Vq_g(:,:,im,ig)
         enddo 
      enddo
   end subroutine

   subroutine cal_uk_formf(xk,uk)
      implicit none 
      real(dp),INTENT(IN) :: xk(3)
      complex(dp),INTENT(OUT) :: uk(nbnd,nbnd,nsvd)
      complex(dp) :: expkr(nre) 
      integer :: ib,jb,im,ire,ig,ijb
      do ire = 1,nre 
         expkr(ire) = exp( dcmplx(0.d0,1.d0)*sum(xk*Re(:,ire))* twopi )
      enddo
      do ib = 1,nbnd; do jb = 1,nbnd 
         ijb = nbnd*(ib-1) + jb
         do ig = 1,Nsvd !summation over singular mode 
            uk(ib,jb,ig) = sum(Ur_g(:,ig,ijb) * expkr)
         end do 
      enddo;enddo
   end subroutine

   subroutine cal_vq_formf(xq,vq)
      implicit none 
      real(dp),INTENT(IN) :: xq(3)
      complex(dp),INTENT(OUT) :: vq(nbnd,nbnd,nmod,nsvd)
      complex(dp) :: expqr(nrp) 
      integer :: ib,jb,im,irp,ig,ijb
      do irp = 1,nrp
         expqr(irp) = exp( dcmplx(0.d0,1.d0)*sum(xq*Rp(:,irp))* twopi )
      enddo
      do im = 1,nmod; do ib = 1,nbnd; do jb = 1,nbnd 
         ijb = nbnd*(ib-1) + jb
         do ig = 1,Nsvd !summation over singular mode 
            vq(ib,jb,im,ig) = sum(Vr_g(:,im,ig,ijb) * expqr)
         end do 
      enddo;enddo;enddo 
   end subroutine

   !calculate the Re = 0 contribution to e-ph matrix 
   subroutine cal_gkq_Re0(xkt, xqt, gkq)
      implicit none 
      real(dp),INTENT(IN) :: xkt(3), xqt(3)
      complex(dp),INTENT(out) :: gkq(nbnd,nbnd,nmod)
      complex(dp) :: gkqc(nbnd,nbnd,nmod)
      integer :: ig, im, ib,jb,ijb,ire,irp 
      complex(dp) :: expqr(nrp), cexpqr(nrp) 

      !g(q)
      do irp = 1,nrp
         expqr(irp) = exp( dcmplx(0.d0,1.d0)*sum(xqt*Rp(:,irp))* twopi )
      enddo
      gkq = 0.d0 
      do ib = 1,nbnd
         do jb = 1,nbnd 
            if(ib.ne.jb) cycle !experimental, keep only diagonal part 
            do im = 1,nmod 
               gkq(ib,jb,im) = sum(expqr*g_Rp(:,ib,jb,im)) 
            enddo 
         enddo 
      enddo 

      


   end subroutine


   !----------------------------------------------------!
   !---------------- svd on g_{ReRp} -------------------!
   !----------------------------------------------------!   
   subroutine init_eph_svd(el,ph,ep)
      use Lapack_LY, only : ZSVD, ZHEIGEN
      implicit none 
      type(lattice_ifc),INTENT(IN) :: ph       !phonon dispersion 
      type(electron_wann),INTENT(IN) :: el     !el dispersion 
      type(elph_mat_wann),INTENT(IN) :: ep     !e-ph element 
      integer :: ib, jb, ia, ix, im, ijb !choose 
      complex(dp),allocatable :: workg(:,:),g_ReRp(:,:),U(:,:),V(:,:), VT(:,:), V_reshape(:,:,:) !(max_nre,max_nrp) 
      complex(dp), allocatable :: gtmp(:,:)
      real*8,ALLOCATABLE :: S(:),AS(:,:)
      integer :: ire,irp 
      type(ws_cell), pointer :: wel,wph
      integer :: NS 
      integer :: iunit, iwork 
   
      nmod = ph%nm; nbnd = el%nb 
      
      nre = ep%nrvec_e
      nrp = ep%nrvec_p
      NS = nre 

      Nsvd_tot = nre 
      if(Nsvd_tot>100)Nsvd_tot = 100
      !if(read_formf) return !if we just readin form factor as it is tabulated before, we skip the following svd procedure 
      !allocate eph_svd 
      ALLOCATE(S_g(NS,el%nb**2), Ur_g(nre, Nsvd_tot, el%nb**2), Vr_g(nrp, ph%nm,  Nsvd_tot, el%nb**2) )
      allocate(g_Rp(nrp, nbnd, nbnd, nmod), gtmp(nbnd,nbnd))
      allocate(Re(3,nre),Rp(3,nrp))
      
   
      Ur_g = 0; Vr_g = 0; S_g = 0;
   
      if(ionode) then 
         write(*,'(A20,3i6)')'nre,nrp,Nsvd = ',nre,nrp,Nsvd 
      endif 
      iwork = 0
      if(read_svd) then 
         if(ionode) then 
            write(*,'(A20)')'read in svd g_{ReRp} '

         endif 
         call Read_eph_svd()
      
      else 
         allocate(g_ReRp(nre,nmod*nrp),workg(nre,nre), U(nre,nre), V(nre,nmod*nrp), S(nre), AS(nre,nrp))
         allocate( VT(nmod*nrp,nre), V_reshape(nrp, nmod, nre) )
         g_Rp = 0.d0
         do ib = 1,el%nb; do jb = 1,el%nb; 
            iwork = iwork + 1
            ijb = (ib-1)*el%nb + jb 
            if(mod(iwork, npool).ne.my_pool_id) cycle !mpi split of jobs 
            g_ReRp = 0  
            do ia = 1,ep%na
               wel => ep%epwan(ib,jb,ia)%ws_el
               wph => ep%epwan(ib,jb,ia)%ws_ph
               !write(*,*) wel%rvec(1)
               do ix = 1,3
                  im = ix + (ia-1)*3
                  do ire = 1,wel%nr; do irp = 1,wph%nr;  
                     g_ReRp( wel%rvec(ire), (im-1)*nrp + wph%rvec(irp) ) = ep%epwan(ib,jb,ia)%ep_hop(ix,ire,irp)
                  enddo; enddo 
                  do irp = 1,wph%nr;  
                     g_Rp( wph%rvec(irp), ib, jb, im ) = ep%epwan(ib,jb,ia)%ep_hop(ix,1,irp)
                  enddo
               enddo
            enddo
            !DO SVD 
            workg = matmul(g_ReRp, transpose(conjg(g_ReRp)))
            U = -workg 
            call ZHEIGEN(nre, U, S)
            S_g(:,ijb) = -S
            V = matmul(transpose(conjg(U)), g_ReRp)
            VT = transpose(V)
            V_reshape = RESHAPE( VT, (/nrp, nmod, nre/))
            Ur_g(:,:,ijb) = U(:,1:Nsvd_tot)
            Vr_g(:,:,:,ijb) = V_reshape(:,:,1:Nsvd_tot)
            if(ionode) write(*,'(i5,A3,i5)')iwork,'/', el%nb**2

            !hermitrize g_Rp 
            do irp = 1,nrp 
               do im = 1,nmod 
                  gtmp = g_Rp(irp,:,:,im)
                  ! version 1.0, Hermition 
                  g_Rp(irp,:,:,im) = 0.5_dp * (gtmp + transpose(conjg(gtmp)))
               enddo 
            enddo 
            

            
         enddo; enddo
         call mp_sum(S_g,inter_pool_comm)
         call mp_sum(Ur_g,inter_pool_comm)
         call mp_sum(Vr_g,inter_pool_comm)
         call mp_sum(g_Rp,inter_pool_comm)
         call Save_eph_svd()
         DEALLOCATE(g_ReRp,workg,U,V,VT,V_reshape,S,AS)
         call Save_eph_sv_dat()
      endif 
      
      do ire = 1,nre
         Re(:,ire) = ep%rvec_set_el(:,ire)
      end do  
      !write(*,*) Re(:,1)
      
      !write(*,'(A10,3f10.5)') 'Re = ',Re(:,1)
      do irp = 1,nrp
         Rp(:,irp) = ep%rvec_set_ph(:,irp)
      end do  
      if(ionode) then 
         write(*,'(A20,3E15.5)')'check max (S,U,V) : ',maxval(abs(S_g)),maxval(abs(Ur_g)),maxval(abs(Vr_g))
      endif 
     
      call mp_barrier(inter_pool_comm)
      if(ionode) write(*,'(A30)')'init_eph_svd passed' 
   end subroutine

   subroutine Save_eph_svd()
      implicit none 
      integer :: iunit,ig
      character(len=20) :: num_g  
      if(.not.ionode) return  


      iunit = find_free_unit()
      open(iunit,file=trim(adjustl(svd_dir))//'/eph_U.dat',action='WRITE',form='Unformatted')
      write(iunit)Ur_g
      close(iunit)
      open(iunit,file=trim(adjustl(svd_dir))//'/eph_V.dat',action='WRITE',form='Unformatted')
      write(iunit)Vr_g
      close(iunit)
      open(iunit,file=trim(adjustl(svd_dir))//'/eph_S.dat',action='WRITE',form='Unformatted')
      write(iunit)S_g
      close(iunit)
      open(iunit,file=trim(adjustl(svd_dir))//'/g_Rp.dat',action='WRITE',form='Unformatted')
      write(iunit)g_Rp
      close(iunit)


      !iunit = find_free_unit()
      !do ig = 1,Nsvd_tot
      !   WRITE(num_g,'(i5)') num_g
      !   open(iunit,file=trim(adjustl(svd_dir))//'/eph_U.dat-'//trim(adjustl(num_g)),action='WRITE',form='Unformatted')
      !   write(iunit)Ur_g(:,:,:,ig)
      !   close(iunit)
      !   open(iunit,file=trim(adjustl(svd_dir))//'/eph_V.dat-'//trim(adjustl(num_g)),action='WRITE',form='Unformatted')
      !   write(iunit)Vr_g
      !   close(iunit)
      !   open(iunit,file=trim(adjustl(svd_dir))//'/eph_S.dat-'//trim(adjustl(num_g)),action='WRITE',form='Unformatted')
      !   write(iunit)S_g
      !   close(iunit)
      !enddo
   end subroutine 

   subroutine Save_eph_sv_dat()
      implicit none 
      integer :: iunit,ib,jb,imod,ig, NS 
      real(dp) :: xtemp 
      if(.not.ionode) return  
      iunit = find_free_unit()
      NS = nre
      open(iunit,file='eph_S.dat',action='WRITE')
      do ib = 1,nbnd**2
         do ig = 1,NS
            xtemp = S_g(ig,ib)
            if(abs(xtemp)<1e-10) xtemp = 0
            write(iunit,'(E20.10)',advance = 'no') xtemp 
         enddo;
         write(iunit,'(A3)')' '
      enddo 
      close(iunit)
   end subroutine 

   subroutine Read_eph_svd()
      implicit none 
      integer :: iunit 
      if(ionode) then 
         iunit = find_free_unit()
         open(iunit,file=trim(adjustl(svd_dir))//'/eph_U.dat',action='read',form='Unformatted')
         read(iunit)Ur_g
         close(iunit)
         open(iunit,file=trim(adjustl(svd_dir))//'/eph_V.dat',action='read',form='Unformatted')
         read(iunit)Vr_g
         close(iunit)
         open(iunit,file=trim(adjustl(svd_dir))//'/eph_S.dat',action='read',form='Unformatted')
         read(iunit)S_g
         close(iunit)
         open(iunit,file=trim(adjustl(svd_dir))//'/g_Rp.dat',action='read',form='Unformatted')
         read(iunit)g_Rp
         close(iunit)
      end if 
      call mp_bcast(Ur_g, ionode_id, inter_pool_comm)
      call mp_bcast(Vr_g, ionode_id, inter_pool_comm)
      call mp_bcast(S_g, ionode_id, inter_pool_comm)
      call mp_bcast(g_Rp, ionode_id, inter_pool_comm)
      if(ionode) write(*,'(A10,E20.10)')'|g_Rp| = ',sum(abs(g_Rp))
   end subroutine 

end module 

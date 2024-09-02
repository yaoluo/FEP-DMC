!debug subroutine for eph svd
subroutine Test_eph_svd(el,ph,ep)
   use eph_svd
   implicit none
   type(electron_wann),INTENT(IN) :: el     !el dispersion  
   type(lattice_ifc),INTENT(IN) :: ph       !phonon dispersion 
   type(elph_mat_wann),INTENT(INOUT) :: ep     !e-ph element 
   type(vlist) :: ql,kl
   real(dp) :: xkt(3), xqt(3),dp2(ph%nm), gm2(ph%nm)
   integer :: iq,ik,im,ib,jb,i,j 
   integer :: numb 
   real(dp) :: wqt(ph%nm), ekt(el%nb), tot_mass 
   complex(dp) :: mq(ep%na*3,ep%na*3), uk(ep%nb,ep%nb), ukq(ep%nb,ep%nb)  
   complex(dp) :: gkq(el%nb,el%nb,ph%nm)
   real(dp),ALLOCATABLE :: g2(:,:,:), dpot(:,:,:), gmod(:,:,:), wq(:,:)
   character(len=80) :: fname
   
   ep%lpol = .false.
   !load in HSP 
   call load_vector_list(fqlist, ql)
   call load_vector_list(fklist, kl)

   tot_mass = sum(mass)
   numb = band_max - band_min + 1

   allocate(g2(el%nb,el%nb,ph%nm), wq(ph%nm, ql%nvec) )
   !
   allocate(dpot(ph%nm, ql%nvec, 1), gmod(ph%nm, ql%nvec, 1))

   gmod = 0.0;wq = 0.0;dpot =0.0; 
   xkt = 0; ik = 1
   call solve_eigenvalue_vector(el, xkt, ekt, uk)
   do iq = 1,ql%nvec 
      if(mod(iq,npool).ne.my_pool_id) cycle 
      xqt = ql%vec(:,iq)
      call solve_eigenvalue_vector(el, xkt+xqt, ekt, ukq)
      call solve_phonon_modes(ph, xqt, wqt, mq); wq(:,iq) = wqt(:)
      call cal_gkq_svd(xkt,xqt,gkq)
      !call eph_transform(ep, xqt, mq, uk, ukq, gkq)
      call eph_transform(ep, xqt, mq, uk, ukq, gkq)
      !
      !renormalized by wq 
      g2 = abs2(gkq)
      
      dp2 = 0.0_dp
      gm2 = 0.0_dp
      do im = 1, ph%nm
         do jb = band_min, band_max
         do ib = band_min, band_max
            dp2(im) = dp2(im) + g2(ib,jb,im)
            if(wqt(im)>phfreq_cutoff) gm2(im) = gm2(im) + g2(ib,jb,im)*0.5_dp/wqt(im)
         enddo; enddo
      enddo
      !check degenerate phonon modes
      im = 1
      do while( im < ph%nm )
         i = 0
         do j = im+1, ph%nm
            if( abs(wqt(j) - wqt(im)) > 1.0E-12_dp ) exit
            ! there are degenerate modes
            i = i + 1
         enddo
         if(i > 0) then
            dp2( im:(im+i) ) = sum(dp2( im:(im+i) )) / real(i+1, dp)
            gm2( im:(im+i) ) = sum(gm2( im:(im+i) )) / real(i+1, dp)
         endif
         !update to the next mode
         im = im + i + 1
      enddo
      
      do im = 1, ph%nm
         gmod(im, iq, ik) = sqrt( gm2(im) / numb)
         dpot(im, iq, ik) = sqrt( dp2(im) * tot_mass /  numb )
      enddo

   enddo 
   call mp_sum(dpot, inter_pool_comm)
   call mp_sum(gmod, inter_pool_comm)
   call mp_sum(wq, inter_pool_comm)

   fname = trim(prefix)//'.ephmat_svd'
   call output_ephmat(fname, kl, ql, ph%nm, wq, dpot, gmod)
   deallocate(wq, dpot, gmod)

   call mp_barrier(inter_pool_comm)
   if(ionode) write(*,'(A30)')'Test_eph_svd passed' 
end subroutine

subroutine Decompose_gr_global(el,ph,ep)
   use eph_svd
   use Lapack_LY, only : ZSVD
   implicit none 
   type(lattice_ifc),INTENT(IN) :: ph       !phonon dispersion 
   type(electron_wann),INTENT(IN) :: el     !el dispersion 
   type(elph_mat_wann),INTENT(IN) :: ep     !e-ph element 
   integer :: ib, jb, ia, ix, im, ijb !choose 
   complex*16,allocatable :: workg(:,:),g_ReRp(:,:),U(:,:),V(:,:) !(max_nre,max_nrp) 
   real*8,ALLOCATABLE :: S(:),AS(:,:)
   integer :: ire,irp 
   type(ws_cell), pointer :: wel,wph
   integer :: NS, M,N
   integer :: iunit, iwork 

   nmod = ph%nm; nbnd = el%nb 
   
   nre = ep%nrvec_e
   nrp = ep%nrvec_p
   M = nre*nbnd*nbnd; N = nrp*nmod
   NS = M; if(N < NS) NS = N
   if(read_formf) return !if we just readin form factor as it is tabulated before, we skip the following svd procedure 
   !allocate eph_svd 
   !ALLOCATE(S_g(NS,ph%nm,el%nb**2), Ur_g(M, M), Vr_g(nrp, nrp, ph%nm, el%nb**2) )
   !allocate(Re(3,nre),Rp(3,nrp))
   allocate(g_ReRp(M,N),workg(M,N), U(M,M),V(N,N),S(NS),AS(M,N))

   !Ur_g = 0; Vr_g = 0; S_g = 0;

   if(ionode) then 
      write(*,'(A20,2i6)')'nre,nrp = ',nre,nrp 
   endif 
   iwork = 0
   if(read_svd) then 
      if(ionode) then 
         write(*,'(A20)')'read in svd g_{ReRp} '
      endif 
      call Read_eph_svd()
   else 
      g_ReRp = 0.d0
      do ia = 1,ep%na; do ib = 1,el%nb; do jb = 1,el%nb; 
         !iwork = iwork + 1
         !if(mod(iwork, npool).ne.my_pool_id) cycle !mpi split of jobs 
         wel => ep%epwan(ib,jb,ia)%ws_el
         wph => ep%epwan(ib,jb,ia)%ws_ph
         ijb = (ib-1)*el%nb + jb 
         do ix = 1,3; im = ix + (ia-1)*3
            do ire = 1,wel%nr; do irp = 1,wph%nr;  
               g_ReRp( ijb+(wel%rvec(ire)-1)*nbnd*nbnd, (wph%rvec(irp)-1)*nmod+im ) = ep%epwan(ib,jb,ia)%ep_hop(ix,ire,irp)
            enddo; enddo  
            !DO SVD 
            !workg = g_ReRp
            !call ZSVD( M, N, workg, U, S, V )
            !V = transpose(conjg(V))
            !AS = 0.d0; do ire = 1,NS; AS(ire,ire) = S(ire); enddo
            !write(*,'(A15,E15.5)')'svd err = ',maxval( abs(g_ReRp - matmul(matmul(U,AS), transpose(conjg(V)))) )
            !S_g(:,im,ijb) = S
            !Ur_g(:,:,im,ijb) = U
            !Vr_g(:,:,im,ijb) = conjg(V)
         enddo
         !if(ionode) write(*,'(i5,A3,i5)')iwork,'/', ep%na*el%nb**2
      enddo; enddo; enddo
      write(*,'(A30)')'g is setup, begin svd'
      call ZSVD( M, N, g_ReRp, U, S, V )

      !call mp_sum(S_g,inter_pool_comm)
      !call mp_sum(Ur_g,inter_pool_comm)
      !call mp_sum(Vr_g,inter_pool_comm)
      !call Save_eph_svd()
      iunit = find_free_unit()
      open(iunit,file='S_tot.dat')
      write(iunit,'(E15.5)')S
      close(iunit)
   endif 

end subroutine

subroutine SVD_elph()
   use eph_svd 
   implicit none

   type(lattice_ifc) :: ph
   type(electron_wann) :: el
   type(elph_mat_wann) :: ep

   !
   if(ionode) then 
      write(stdout,'(A40)')'doing svd decomposition of g_{ReRp}'
   endif 
   !recalculate svd 
   read_svd = .false.
   read_formf = .false.

   call init_lattice_ifc(epwan_fid, qc_dim, ph)
   call init_electron_wann(epwan_fid, kc_dim, el)
   call init_elph_mat_wann(epwan_fid, kc_dim, qc_dim, ep)
   call init_eph_svd(el,ph,ep)

end subroutine 


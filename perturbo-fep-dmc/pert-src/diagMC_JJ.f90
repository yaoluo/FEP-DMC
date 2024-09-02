

subroutine diagmc_JJ()
   USE OMP_LIB
   use DiagMC 
   use diagMC_gnt_debug
   use JJ_update, only : updateJJ_drive
   
   implicit none 
   !omp shared variables 
   integer,parameter ::  Ntype = 7
   integer :: Nthread,iunit

   !omp private variables 
   integer :: i_omp, i_tmp, j_tmp, Nwork  
   integer :: i, update 
   type(diagmc_stat) :: stat(200), stat_total      
   type(fynman) :: diagrams(200) 
   type(diagmc_workspace) :: workspaces(200)

   real(dp) :: PA(Ntype), xr, process 
   integer :: finished = 0  
   !!$omp threadprivate(diagram,stat)
   call start_clock('diagmc_JJ')

   call setup_dqmc()
   !if(.not.LatticeFrohlich .and. .not.Holstein) call check_epw()

   
   !!$omp parallel default(shared)
   Nthread = OMP_GET_MAX_THREADS()
  
   !!$omp end parallel 

   write(stdout,'(A40,f6.3,2i6)')'DiagMC start : dtau, Nthread, maxOrder = ',config%dtau,Nthread, maxOrder

   call init_observable(stat_total) 
   !--------------------------------------------------------------------
   !------------------- below is openmp paralleled ---------------------
   !---------------------------------------------------------------------
   
   do i_omp = 1,Nthread
      call InitVertexList(maxN, Nph, dmc_band, nbnd,nsvd, diagrams(i_omp)%vertexList)
      call init_diagram_JJ(diagrams(i_omp))
      call init_observable(stat(i_omp))
      call init_diagmc_workspace(workspaces(i_omp))
   enddo
   call start_clock('mcmc')
   !$omp parallel default(shared) firstprivate(i_omp,i,update,PA,process,Nwork) 
      Nwork = config%Nmcmc
      i_omp = omp_get_thread_num() + 1
      call init_rand_seed_omp(i_omp + 32*my_pool_id, diagrams(i_omp)%seed)
      PA = config%PA(1:Ntype)
      do i = 1,config%Nmcmc
         ![1] select update type 
         call Discrete_omp(diagrams(i_omp)%seed,PA,Ntype,update)

         ![2] proposal new configuration according to the update type and calculate the acceptance ratio 
         call updateJJ_drive(diagrams(i_omp),stat(i_omp),update)

         ![3] print & check diagram 
         if(mod(i,10**5).eq.0.and.i_omp.eq.1 .and. ionode)  then 
            write(*,'(A40,2i10,2E15.5)')' step, order , Tau, Re(gfac) = ',&
                     i,diagrams(i_omp)%order, diagrams(i_omp)%vertexList(maxN)%tau,real(diagrams(i_omp)%gfac)
            call check_diagram(diagrams(i_omp))
         endif 

         ![4] gather data 
         if(mod(i,1000).eq.0 .and. i .ge. config%Nmcmc/4 ) then
         !if(mod(i,1000).eq.0 ) then  !debug
            call measureJJ_abinito(diagrams(i_omp), workspaces(i_omp), stat(i_omp))
         endif 
         !if(mod(i,10**5).eq.0) then 
         !   call check_phase(diagrams(i_omp))
         !end if
      end do
      !omp gather, need implement  
      !$omp critical 
         stat_total%accept = stat_total%accept + stat(i_omp)%accept 
         stat_total%update = stat_total%update + stat(i_omp)%update 
         stat_total%order = stat_total%order + stat(i_omp)%order 
         stat_total%order_g = stat_total%order_g + stat(i_omp)%order_g 
         stat_total%JJw = stat_total%JJw + stat(i_omp)%JJw
         stat_total%Njjw = stat_total%Njjw + stat(i_omp)%Njjw 
         stat_total%gtrue = stat_total%gtrue + stat(i_omp)%gtrue
      !$omp end critical 
   !$omp end parallel

   stat_total%JJw = stat_total%JJw / stat_total%Njjw
   stat_total%gtrue = stat_total%gtrue / stat_total%Njjw
   stat_total%JJw = stat_total%JJw / stat_total%gtrue 

  
   call stop_clock('mcmc')
   write(stdout, '(A20, 2E20.10)')'<sign> = ',stat_total%gtrue
   write(stdout,'(A30,i15,A7)')'total work = ', stat_total%Njjw/10000,'x1e4'
   write(stdout,'(A30,7i15)')'update try = ',stat_total%update
   write(stdout,'(A30,7E15.5)')'acceptance = ',(1.d0*stat_total%accept)/stat_total%update
   call stop_clock('diagmc_JJ')

   call output_JJ(stat_total)

   if(ionode) then 
      iunit = find_free_unit()
      open(iunit,file='parameter.dat')
      write(iunit,'(A10)') '#mu = '
      WRITE(iunit,'(E15.5)')config%mu 
      write(iunit,'(A10)') '#alpha = '
      if(Holstein) then 
         WRITE(iunit,'(E15.5)')alpha_Holstein
      else 
         WRITE(iunit,'(E15.5)')alpha_frohlich 
      endif 
      write(iunit,'(A10)') '#Ep = '
      WRITE(iunit,'(E15.5)')config%EP + config%mu
      close(iunit)
   endif 
end subroutine 

!measure current-current correlation function : fully abinitio
subroutine measureJJ_abinito(diagram, workspace,stat)
   use pert_param, only : lowest_order
   use diagmc, only : dp, diagmc_stat, fynman, diagmc_workspace,vertex, config, nbnd, dmc_band, &
             direction_mob, maxOrder,maxN,trace,band_min,band_max
   use HamiltonianTable, only : solve_electron_velocity
   use pert_const, only : ci,twopi 
   implicit none 
   type(diagmc_stat) :: stat
   type(fynman),target :: diagram 
   type(diagmc_workspace), target :: workspace
   integer :: itau,iv,order, iw, jv, itr, iel, jel, ijel, ib, ie,ix,iy       
   real(dp) :: tau1,tau2,tau3,tau4, kpt(3), vkt(dmc_band), wmata, JJw(config%Nbin), dtau
   type(vertex), pointer :: vA,vB, vi, vj
   complex(dp), pointer :: expH(:,:,:), gkq_list(:,:,:),vk_mat(:,:,:,:), right_env(:,:,:), left_env(:,:,:)
   real(dp), pointer :: vk_xyz(:,:,:)
   complex(dp) :: L(dmc_band,dmc_band), R(dmc_band,dmc_band), M(dmc_band,dmc_band)
   real(dp) :: vk_full(3,nbnd),vk(dmc_band), tau_list(maxOrder+1),Elist(dmc_band)
   complex(dp) :: vv(dmc_band,dmc_band), tt_iw(config%Nbin,2), iden(dmc_band,dmc_band), vvt(3,3),z_i, weight, weight2   
   

   expH => workspace%expH
   gkq_list => workspace%gkq
   vk_mat => workspace%Ak_mat
   left_env => workspace%left_env
   right_env => workspace%right_env
   vk_xyz => workspace%vk_xyz

   !order stat
   order = (diagram%order-1)/2 + 1
   stat%Order(order) = stat%Order(order) + 1
   stat%Njjw = stat%Njjw + 1

   !stat%order_g(order) = stat%order_g(order) + diagram%gfac
   
   iv = 1
   vi => diagram%vertexList(iv)
   if(diagram%order.eq.1) then   
      vvt = 0.d0; z_i = 0.d0  
      Elist = vi%ekout - minval(vi%ekout)

      do ib = 1,dmc_band 
         do ix = 1,3; do iy = 1,3 
            vvt(ix,iy) = vvt(ix,iy) + vi%vkout(ix, ib) * vi%vkout(iy, ib) * exp(-min( config%tauMax * Elist(ib), 100.0_dp ) )
         enddo;enddo 
         z_i = z_i + exp( -min( config%tauMax * Elist(ib), 100.0_dp ) )
      enddo 
      stat%JJw(:,:,1) = stat%JJw(:,:,1) + vvt * config%tauMax / z_i 
      stat%gtrue = stat%gtrue  + 1.d0
      return 
   endif 
   
   
   iden = 0
   do ib = 1,dmc_band
      iden(ib,ib) = 1.d0  
   enddo 
   !construct propogator and velocity matrix 
   vk_mat = 0.d0 
   left_env(:,:,1) = iden   
   do itr = 1,diagram%order 
      jv = vi%link(3)
      vj => diagram%vertexList(jv)
      !green function matrix 
      Elist = vi%ekout - minval(vi%ekout)

      tau_list(itr) = vi%tau
      dtau = vj%tau - vi%tau
      call propagator(vi%ukout, Elist, dtau, dmc_band, expH(:,:,itr) )
      do ib = 1,dmc_band
         vk_mat(ib,ib,:,itr) = vi%vkout(:,ib-1+band_min) !velocity matrix at itr-th electron line in time ordered.
      enddo 
      
      if(jv .ne. maxN) then 
         gkq_list(:,:,itr) = vj%gkq(:,:,vj%nu) / maxval(abs(vj%gkq(:,:,vj%nu))) 
      else 
         gkq_list(:,:,itr) = iden 
      endif 
      left_env(:,:,itr+1) = matmul( expH(:,:,itr), left_env(:,:,itr))
      left_env(:,:,itr+1) = matmul( gkq_list(:,:,itr), left_env(:,:,itr+1))

      iv = vi%link(3)
      vi => diagram%vertexList(iv); 
   enddo 
   tau_list(diagram%order+1) = diagram%vertexList(maxN)%tau !beta/tauMax

   !construct right env 
   right_env(:,:,diagram%order) = iden 
   do itr = diagram%order,2,-1 
      right_env(:,:,itr-1) = matmul(right_env(:,:,itr), expH(:,:,itr))
      right_env(:,:,itr-1) = matmul(right_env(:,:,itr-1), gkq_list(:,:,itr-1))
   enddo 


   !measure partition function Z : trace[G1 G2 G3 ... Gn]
   weight = trace(left_env(:,:,diagram%order+1), dmc_band) 
   weight = real(weight)
   if(real(weight) < 0.0_dp) then 
      diagram%gfac = -1.0_dp  
   else 
      diagram%gfac = 1.0_dp 
   endif 
   stat%gtrue = stat%gtrue + diagram%gfac

   !measure current-current correlator : trace[G1 G2 ... vi ... Gk ... vj ... Gn]
   do iel = 1,diagram%order     !velocity vertex 1 @ ith-electron line 
      tau1 = tau_list(iel) 
      tau2 = tau_list(iel+1)

      tt_iw(1,1) = tau2-tau1 !zero-freq use a different formula 
      tt_iw(2:,1) = (exp(ci*stat%wmata(2:)*tau2)-exp(ci*stat%wmata(2:)*tau1) ) / (ci*stat%wmata(2:))
      if(iel.eq.1) then 
         L = iden 
      else 
         L = left_env(:,:,iel)
      endif    

      M = expH(:,:,iel)  
      do jel = iel,diagram%order !velocity vertex 2 @ jth-electron line 
         tau3 = tau_list(jel) 
         tau4 = tau_list(jel+1)
         tt_iw(1,2) = tau4-tau3 
         tt_iw(2:,2) = (exp(-ci*stat%wmata(2:)*tau4)-exp(-ci*stat%wmata(2:)*tau3) ) / (-ci*stat%wmata(2:))
         R = right_env(:,:,jel)

         ! ! ! debug for weight ! 1 ! 
         vv = matmul(R, M)
         vv = matmul(vv, L)
         weight2 = trace(vv, dmc_band) 
         if(abs(real(weight2) - weight)/abs(weight) > 1e-10) then  
            write(*,'(A40,2E30.20)')'weight inconsistent@JJ measure',real(weight2) , real(weight)
         endif 

         !vv = Tr[R*Vj*M*Vi*L]
         do ix = 1,3
            iy = ix 
            vv = matmul(R,  vk_mat(:, :, ix, jel))
            vv = matmul(vv, M)
            vv = matmul(vv, vk_mat(:, :, iy, iel))
            vv = matmul(vv, L)
            vvt(ix,iy) = trace(vv,dmc_band) / abs(real(weight))
            if(iel.eq.jel) then 
               stat%JJw(ix,iy,:) = stat%JJw(ix,iy,:) + real(vvt(ix,iy)) * real(tt_iw(:,1)*tt_iw(:,2)) / config%tauMax
            else 
               stat%JJw(ix,iy,:) = stat%JJw(ix,iy,:) + 2.0_dp * real(vvt(ix,iy)) * real(tt_iw(:,1)*tt_iw(:,2)) / config%tauMax
            endif 
         enddo


         if(maxval(abs(vvt))>1e10) then 
            write(*,'(A30)') "measure JJ error : vv > 1e10" 
            stop 
         endif  
         !update M 
         if(jel.eq.diagram%order) cycle  !no need to update M for the last iteration 
         M = matmul(gkq_list(:,:,jel),M) !jel == diagram%order, gkq = identity matrix 
         M = matmul(expH(:,:,jel+1),M)
      enddo 

     !if(iel<=1) cycle  
     !M = matmul(expH(:,:,iel),gkq_list(:,:,iel-1) ) 
     !M = matmul(M, expH(:,:,iel-1)) 
     !R = right_env(:,:,iel)
     !do jel = iel-1, 1, -1 
     !   tau3 = tau_list(jel) 
     !   tau4 = tau_list(jel+1)
     !   tt_iw(1,2) = tau4-tau3 
     !   tt_iw(2:,2) = (exp(-ci*stat%wmata(2:)*tau4)-exp(-ci*stat%wmata(2:)*tau3) ) / (-ci*stat%wmata(2:))
     !   if(jel.eq.1) then 
     !      L = iden 
     !   else 
     !      L = left_env(:,:,jel)
     !   endif    
     !   ! ! ! debug for weight ! 1 ! 
     !   !vv = matmul(R, M)
     !   !vv = matmul(vv, L)
     !   !weight2 = trace(vv, dmc_band) 
     !   !write(*,'(3E30.20)')abs(real(weight2) - weight)/abs(weight), real(weight2) , real(weight)
     !   !if(abs(real(weight2) - weight)/abs(weight) > 1e-4) then  
     !   !   write(*,'(A40,2E30.20)')'weight inconsistent@JJ measure',real(weight2) , real(weight)
     !   !endif 

     !   !vv = Tr[R*Vj*M*Vi*L]
     !   do ix = 1,3
     !       iy = ix 
     !       vv = matmul(R,  vk_mat(:, :, iy, iel))
     !       vv = matmul(vv, M)
     !       vv = matmul(vv, vk_mat(:, :, ix, jel))
     !       vv = matmul(vv, L)
     !       vvt(ix,iy) = trace(vv,dmc_band) / abs(real(weight))
     !       stat%JJw(ix,iy,:) = stat%JJw(ix,iy,:) + real(vvt(ix,iy)) * real(tt_iw(:,1)*tt_iw(:,2)) / config%tauMax
     !   enddo

     !   if(maxval(abs(vvt))>1e10) then 
     !      write(*,'(A30)') "measure JJ error : vv > 1e10" 
     !      stop 
     !   endif  
     !   !update M 
     !   if(jel.eq.1) cycle !no need for the last iteration 
     !   M = matmul(M, gkq_list(:,:,jel-1))
     !   M = matmul(M, expH(:,:,jel-1)) !jel == diagram%order, gkq = identity matrix 
     !enddo 
   enddo 
end subroutine 

!output 
subroutine output_JJ(stat,suffix_name)
   use diagMC
   use pert_param, only : Nsample_MC
   implicit none
   type(diagmc_stat) :: stat
   character(len=100), optional :: suffix_name
   character(len=100) :: str  
   integer :: iunit, it, iw, order,ib,jb,job    
   real(dp) :: dtau, tau_ev
   INTEGER*4  :: access,status

   do job = 1,Nsample_MC
      write(str,'(i10)') job 
      status = access( 'JJ.dat-'//trim(adjustL(str)),  ' ' )    ! blank mode
      if ( status .ne. 0 ) then 
         exit 
      endif 
   enddo 

   iunit = find_free_unit()

   open(iunit,file='JJ.dat-'//trim(adjustL(str)))
   write(iunit,'(3A30)')'# mata freq' , '<J(iw_n)J(-iw_n)>'
   do iw = 1,config%Nbin
      write(iunit,'(E30.20,9E30.20)') stat%wmata(iw), stat%JJw(:,:,iw)
   enddo 
   close(iunit)    

   open(iunit,file='orderStat.dat-'//trim(adjustL(str)))
   do it = 1,maxOrder/2+1
      write(iunit,'(3E30.20)') real(stat%order(it),dp), stat%order_g(it)
   enddo 
   close(iunit)  

   open(iunit,file='sign.dat-'//trim(adjustL(str)))
   write(iunit,'(2E30.20)') stat%gtrue
   close(iunit)   
end subroutine 


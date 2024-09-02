

subroutine diagmc_gt_matrix()
   USE OMP_LIB
   use DiagMC 
   use diagMC_gnt_debug
   use zerophonon_update_matrix, only : update_drive
   implicit none 
   !omp shared variables 
   integer,parameter ::  Ntype = 7
   integer :: Nthread,iunit

   !omp private variables 
   integer :: i_omp, i_tmp, j_tmp, Nwork  
   integer :: i, update 
   type(diagmc_stat) :: stat(200), stat_total      
   type(fynman) :: diagrams(200) 
   type(diagmc_workspace) :: workspace(200) 

   real(dp) :: PA(Ntype), xr, process 
   integer :: finished = 0, job_process, total_work    

   multiband = .true.
   !!$omp threadprivate(diagram,stat)
   call start_clock('diagmc_gt')

   call setup_dqmc()
   !if(.not.LatticeFrohlich .and. .not.Holstein) call check_epw()
   !call ReSelfenergy_degenerate_frolich(el,ph,ep)
   !stop 
   
   !!$omp parallel default(shared)
   Nthread = OMP_GET_MAX_THREADS()
   !!$omp end parallel 

   write(stdout,'(A40,f6.3,2i6)')'DiagMC start : dtau, Nthread, maxOrder = ',config%dtau,Nthread, maxOrder

   call init_observable(stat_total) 
   !--------------------------------------------------------------------
   !------------------- below is openmp paralleled ---------------------
   !---------------------------------------------------------------------
   
   do i_omp = 1,Nthread
      call InitVertexList(maxN, Nph, dmc_band, nbnd, nsvd, diagrams(i_omp)%vertexList)
      call init_diagram(diagrams(i_omp))
      call init_observable(stat(i_omp))
      call init_diagmc_workspace(workspace(i_omp))
   enddo
   call start_clock('mcmc')

   job_process = 0
   total_work =  config%Nmcmc / 10000 * Nthread 
   write(*,'(A20,i10,A10)')'total_work = ',total_work, '( x 1e4)'
   !$omp parallel default(shared) firstprivate(i_omp,i,update,PA,process,Nwork) 
      Nwork = config%Nmcmc
      i_omp = omp_get_thread_num() + 1
      
      call init_rand_seed_omp(i_omp, diagrams(i_omp)%seed)
      PA = config%PA(1:Ntype)
      diagrams(i_omp)%vertexList(maxN)%tau = config%tauMax
      i = 0 
      do while (job_process<total_work)
         i = i + 1
         ![1] select update type 
         call Discrete_omp(diagrams(i_omp)%seed,PA,Ntype,update)
         ![2] proposal new configuration according to the update type and calculate the acceptance ratio 
         call update_drive(diagrams(i_omp),stat(i_omp),update)
         ![3] print & check diagram 
         if(mod(i,10**5).eq.0) then 
            call check_diagram(diagrams(i_omp))
         endif

         ![4] gather data 
         !WRITE(*,*) '1-measure'
         if(mod(i,20).eq.0 .and. i.ge.Nwork/4) then  !drop the first one fourth 
            call measureG_matrix(diagrams(i_omp), workspace(i_omp),stat(i_omp)) !measure first, recaluclate gfac!
         endif 

         if(mod(i,10**5).eq.0.and.i_omp.eq.1 )  then 
            write(*,'(A40,2i10,2E15.5)')' step, order , Tau, Re(gfac) = ',&
                     i,diagrams(i_omp)%order, diagrams(i_omp)%vertexList(maxN)%tau,real(diagrams(i_omp)%gfac)
         endif 
         if(mod(i,10000).eq.9999) then 
            !$omp atomic 
            job_process = job_process + 1 
         endif 
      end do
      !omp gather data
      !stat(i_omp)%GreenFunc_matrix = stat(i_omp)%GreenFunc_matrix 
      !stat(i_omp)%trG_abs = stat(i_omp)%trG_abs 
      !$omp critical 
         stat_total%Njjw = stat_total%Njjw + stat(i_omp)%Njjw 
         stat_total%accept = stat_total%accept + stat(i_omp)%accept 
         stat_total%update = stat_total%update + stat(i_omp)%update 
         stat_total%Tstat = stat_total%Tstat + stat(i_omp)%Tstat  
         stat_total%order = stat_total%order + stat(i_omp)%order  
         stat_total%order_g = stat_total%order_g + stat(i_omp)%order_g  
         stat_total%GreenFunc_matrix = stat_total%GreenFunc_matrix + stat(i_omp)%GreenFunc_matrix  
         stat_total%trG_abs = stat_total%trG_abs + stat(i_omp)%trG_abs  
      !$omp end critical 
      finished = 1
   !$omp end parallel


   call stop_clock('mcmc')

   write(stdout,'(A30,i15,A7)')'total work = ',sum(stat_total%Tstat)/10000,'x1e4'
   write(stdout,'(A30,7i15)')'update try = ',stat_total.update
   write(stdout,'(A30,7E15.5)')'acceptance = ',(1.d0*stat_total.accept)/stat_total.update
   call output_gt(stat_total)
   
   call stop_clock('diagmc_gt')

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
   endif 

end subroutine 


subroutine diagmc_EZ_matrix()
   USE OMP_LIB
   use DiagMC 
   use diagMC_gnt_debug
   use multiphonon_update_matrix, only : update_drive
   implicit none 
   !omp shared variables 
   integer,parameter ::  Ntype = 7
   integer :: Nthread,iunit

   !omp private variables 
   integer :: i_omp, i_tmp, j_tmp, Nwork  
   integer :: i, ib, update 
   type(diagmc_stat) :: stat(200), stat_total      
   type(fynman) :: diagrams(200) 
   type(diagmc_workspace) :: workspaces(200)

   real(dp) :: PA(Ntype), xr, process,njjw 
   integer :: finished = 0, job_process, total_work    

   multiband = .true.
   !!$omp threadprivate(diagram,stat)
   call start_clock('diagmc_gt')

   call setup_dqmc()
   !if(.not.LatticeFrohlich .and. .not.Holstein) call check_epw()
   !call ReSelfenergy_degenerate_frolich(el,ph,ep)
   !stop 
   
   !!$omp parallel default(shared)
   Nthread = OMP_GET_MAX_THREADS()
   !!$omp end parallel 

   write(stdout,'(A40,f6.3,2i6)')'DiagMC start : dtau, Nthread, maxOrder = ',config%dtau,Nthread, maxOrder

   call init_observable(stat_total) 
   !--------------------------------------------------------------------
   !------------------- below is openmp paralleled ---------------------
   !---------------------------------------------------------------------
   
   do i_omp = 1,Nthread
      call InitVertexList(maxN, Nph, dmc_band, nbnd, nsvd, diagrams(i_omp)%vertexList)
      !call init_diagram(diagrams(i_omp))
      call init_diagram(diagrams(i_omp), config%tauMax)
      call init_observable(stat(i_omp))
      call init_diagmc_workspace(workspaces(i_omp))
   enddo
   call start_clock('mcmc')


   job_process = 0
   total_work =  config%Nmcmc / 10000 * Nthread 
   write(*,'(A20,i10,A10)')'total_work = ',total_work, '( x 1e4)'
   !$omp parallel default(shared) firstprivate(i_omp,i,update,PA,process,Nwork) 
      Nwork = config%Nmcmc
      i_omp = omp_get_thread_num() + 1
      
      call init_rand_seed_omp(i_omp, diagrams(i_omp)%seed)
      PA = config%PA(1:Ntype)
      PA = PA/sum(PA)
      diagrams(i_omp)%vertexList(maxN)%tau = config%tauMax
      i = 0 
      do while (job_process<total_work)
         i = i + 1
         ![1] select update type 
         call Discrete_omp(diagrams(i_omp)%seed,PA,Ntype,update)
         !WRITE(*,*) '1-',UPDATE 
         ![2] proposal new configuration according to the update type and calculate the acceptance ratio 
         call update_drive(diagrams(i_omp),stat(i_omp),update)
         ![3] print & check diagram 
         if(mod(i,10**5).eq.0) call check_diagram(diagrams(i_omp))

         ![4] gather data 
         if(mod(i,100).eq.0 .and. i.ge.config%Nmcmc/4) then 
            call measure_EZ_wfn(diagrams(i_omp),workspaces(i_omp),stat(i_omp))
            !call measure_EZ_toy(diagrams(i_omp),stat(i_omp))
         endif 
         if(mod(i,10**5).eq.0.and.i_omp.eq.1 .and. ionode)  then 
            write(*,'(A50,i5,2i10,3E15.5)')' update, step, order , Tau, s/gfac = ',&
                     update, i, diagrams(i_omp)%order, diagrams(i_omp)%vertexList(maxN)%tau, &
                     real(diagrams(i_omp)%sfac), real(diagrams(i_omp)%gfac)
            call check_diagram(diagrams(i_omp))
         endif  
         if(mod(i,10000).eq.9999) then 
            !$omp atomic 
            job_process = job_process + 1 
         endif 
      end do
      !omp gather data
     
      !$omp critical 
         stat_total%Njjw = stat_total%Njjw + stat(i_omp)%Njjw 
         stat_total%Znph = stat_total%Znph + stat(i_omp)%Znph 
         stat_total%PolaronWF%psi = stat_total%PolaronWF%psi + stat(i_omp)%PolaronWF%psi
         stat_total%PolaronLat%psi = stat_total%PolaronLat%psi + stat(i_omp)%PolaronLat%psi
         stat_total%Etrue = stat_total%Etrue + stat(i_omp)%Etrue 
         stat_total%gtrue = stat_total%gtrue + stat(i_omp)%gtrue
         stat_total%accept = stat_total%accept + stat(i_omp)%accept 
         stat_total%update = stat_total%update + stat(i_omp)%update 
         stat_total%Tstat = stat_total%Tstat + stat(i_omp)%Tstat  
         stat_total%order = stat_total%order + stat(i_omp)%order  
         stat_total%order_g = stat_total%order_g + stat(i_omp)%order_g  
         stat_total%order_sample = stat_total%order_sample + stat(i_omp)%order_sample  
      !$omp end critical 
   !$omp end parallel

   
   stat_total%Etrue = stat_total%Etrue / config%tauMax
   stat_total%Etrue = stat_total%Etrue / stat_total%gtrue
   stat_total%gtrue = stat_total%gtrue / stat_total%Njjw
   call output_EZ(stat_total)
   stat_total%Znph = stat_total%Znph / sum( stat_total%Znph )
   

   call stop_clock('mcmc')
   if(ionode) then 
      write(stdout, '(A20, 2E20.10)')'Etrue = ',stat_total%Etrue
      write(stdout, '(A20, 2E20.10)')'<g/Z> = ',stat_total%gtrue
      write(stdout, '(A20, 3E20.10)')'<Z_0> = ',real(stat_total%Znph(:,1))
      write(stdout, '(A20, 3E20.10)')'<Z_1> = ',real(stat_total%Znph(:,2))
      write(stdout,'(A30,i15,A7)')'total work = ',sum(stat_total%Tstat)/10000,'x1e4'
      write(stdout,'(A30,7i15)')'update try = ',stat_total.update
      write(stdout,'(A30,7E15.5)')'acceptance = ',(1.d0*stat_total.accept)/stat_total.update
   endif 
   


   call mp_sum(stat_total%Etrue, inter_pool_comm)
   call mp_sum(stat_total%gtrue, inter_pool_comm)
   call mp_sum(stat_total%Znph, inter_pool_comm)

   stat_total%Znph = stat_total%Znph / npool
   stat_total%Etrue = stat_total%Etrue / npool
   stat_total%gtrue = stat_total%gtrue / npool

   call stop_clock('diagmc_gt')

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
      write(iunit,'(A10)') '#Etrue = '
      WRITE(iunit,'(E15.5)')real(stat_total%Etrue)
      write(iunit,'(A10)') '#<g/Z> = '
      WRITE(iunit,'(E15.5)')real(stat_total%gtrue)
      write(iunit,'(A10)') '#<Z0> = '
      WRITE(iunit,'(E15.5)')stat_total%Znph(:,1)
      close(iunit)
   endif 

end subroutine 

!measure GreenFunc
subroutine measureG(diagram,stat)
   use diagmc, only : diagmc_stat, fynman, maxN,config
   implicit none 
   type(diagmc_stat) :: stat
   type(fynman) :: diagram 
   integer :: itau, order 
   order = diagram%order/2  

   
   !write(*,*) diagram%vertexList(maxN)%tau
   itau = floor( (diagram%vertexList(maxN)%tau - config%tauMin) /config%dtau)+1
   if(itau>config%Nbin) itau = config%Nbin
   !write(*,*) itau,order+1

   stat%Tstat(itau) = stat%Tstat(itau) + 1
   stat%GreenFunc(itau,order+1) = stat%GreenFunc(itau,order+1) + diagram%gfac
   stat%Order(order+1) = stat%Order(order+1) + 1
   stat%order_g(order+1) = stat%order_g(order+1) + diagram%gfac
end subroutine 


!measure GreenFunc
subroutine measureG_matrix(diagram,workspace,stat)
   use diagmc, only : Gel,diagmc_stat, vertex,fynman,diagmc_workspace, maxN,config, dp,cone, dmc_band,trace
   use pert_param, only : print_sign, sample_gt, ib_sample
   implicit none 
   type(diagmc_stat) :: stat
   type(fynman),target :: diagram 
   type(diagmc_workspace), target :: workspace
   integer :: itau, order, ib, jb  
   real(dp) :: tau1,tau2, Elist(dmc_band),dtau
   complex(dp) :: g_matrix(dmc_band,dmc_band), G_amp(dmc_band,dmc_band), expH(dmc_band,dmc_band), weight 
   type(vertex), pointer :: v_head,v_tail

   order = diagram%order/2  
   v_head => diagram%vertexList(1)
   v_tail => diagram%vertexList(maxN)
   
   !write(*,*) diagram%vertexList(maxN)%tau
   itau = floor( (v_tail%tau - config%tauMin) /config%dtau)+1
   if(itau>config%Nbin) itau = config%Nbin

   stat%Tstat(itau) = stat%Tstat(itau) + 1
   stat%Order(order+1) = stat%Order(order+1) + 1
   stat%trG_abs(itau) = stat%trG_abs(itau) + 1.d0 

   Elist = v_tail%ekin - minval(v_tail%ekin)
   if(diagram%order.eq.1) then 
      g_matrix = 0.d0 
      weight = 0.d0 
      do ib = 1,dmc_band 
         g_matrix(ib,ib) = Gel(Elist(ib), v_tail%tau)
      enddo 
      weight = trace(g_matrix,dmc_band)
      if(sample_gt) weight = g_matrix(ib_sample, ib_sample)
      g_matrix = g_matrix / weight
      stat%order_g(order+1) = stat%order_g(order+1) + 1.d0 
      stat%GreenFunc_matrix(:,:,itau) = stat%GreenFunc_matrix(:,:,itau) + real(g_matrix)
      return 
   endif 
   tau1 = diagram%vertexList(v_head%link(3))%tau
   tau2 = diagram%vertexList(v_tail%link(1))%tau

   !get left env of the last 
   call left_environment_matrix(diagram, diagram%order, g_matrix)
   dtau = v_tail%tau - diagram%vertexList(v_tail%link(1))%tau
   call propagator(v_tail%ukout, Elist, dtau, dmc_band, expH)
   g_matrix = matmul(expH, g_matrix)

   weight = trace(g_matrix, dmc_band)
   if(sample_gt) weight = g_matrix(ib_sample, ib_sample)
   if(real(weight)<0) then 
      diagram%gfac = -1.d0 
   else 
      diagram%gfac = 1.d0
   endif
   g_matrix =  g_matrix / abs(real(weight))
   
   stat%order_g(order+1) = stat%order_g(order+1) + diagram%gfac
   stat%GreenFunc_matrix(:,:,itau) = stat%GreenFunc_matrix(:,:,itau) + real(g_matrix)

end subroutine 


!measure Epol 
subroutine measure_EZ(diagram, stat)
   use diagmc, only : Gel,diagmc_stat, vertex,fynman, maxN,config, dp,cone, dmc_band,trace
   use pert_param, only : print_sign, nk_svd
   implicit none 
   type(fynman),target, intent(inout) :: diagram 
   type(diagmc_stat),intent(inout) :: stat
   !local 
   integer :: itau, order, ib, jb , itr, iv,jv   
   real(dp) :: tau1,tau2, Elist(dmc_band), iden(dmc_band,dmc_band), dtau, mu, Wdt, tau_barek,  tauList(maxN)
   real(dp) :: dkt(3, maxN), tau_ratio
   integer :: i_dkt(3,maxN)
   integer :: nph_line(maxN), iph, n1,n2
   complex(dp) ::  weight, weight_2, gtrue  
   type(vertex), pointer :: v_head,v_tail, vi,vj
   complex(dp), dimension(dmc_band,dmc_band,maxN) :: expH, gkq,left_env, HexpH
   complex(dp), dimension(dmc_band,dmc_band) :: right_env, Hdt_matrix, G_amp,g_matrix

   order = diagram%order/2  
   v_head => diagram%vertexList(1)
   v_tail => diagram%vertexList(maxN)
   tau_barek = 0.d0 

   !write(*,*) diagram%vertexList(maxN)%tau
   itau = floor( (v_tail%tau - config%tauMin) /config%dtau)+1
   if(itau>config%Nbin) itau = config%Nbin

   stat%Tstat(itau) = stat%Tstat(itau) + 1
   stat%Order(order+1) = stat%Order(order+1) + 1
   stat%Njjw = stat%Njjw + 1
   Elist = v_tail%ekin - minval(v_tail%ekin) 
   if(diagram%order.eq.1) then 
      g_matrix = 0.d0 
      Hdt_matrix = 0.d0 
      weight = 0.d0 
      do ib = 1,dmc_band 
         g_matrix(ib,ib) = Gel(Elist(ib), v_tail%tau)
         Hdt_matrix(ib,ib) = Gel(Elist(ib), v_tail%tau) * (v_tail%ekin(ib)) * config%tauMax
      enddo 
      weight = trace(g_matrix,dmc_band)
      stat%gtrue = stat%gtrue + 1.d0 
      stat%Etrue = stat%Etrue + trace(Hdt_matrix,dmc_band) /abs(weight)

      diagram%gfac = 1.d0 
      do ib = 1,dmc_band 
         stat%Znph(ib,1) = stat%Znph(ib,1) + g_matrix(ib,ib) / weight 
      enddo 
      return 
   endif 

  
   
   !Hdtau contribution 
   !1. prepare propogator along the imaginary time 
   iden = 0.d0
   do ib = 1,dmc_band; iden(ib,ib) = 1.d0; enddo 

   vi => diagram%vertexList(1)
   jv = vi%link(3)
   vj => diagram%vertexList(jv)
   left_env(:,:,1) = iden   

   nph_line = diagram%nph_ext
   n1 = 0; n2 = 0 
   do itr = 1,diagram%order    
      Elist = vi%ekout - minval(vi%ekout)
      mu = minval(vi%ekout)  
      dkt(:,itr) = vi%kout - config%P
      i_dkt(:,itr) = mod(vi%i_kout - config%i_P, nk_svd)
      
      !call To1BZ(dkt(:,itr))
      dtau = vj%tau-vi%tau
      tauList(itr) = vi%tau
      tauList(itr+1) = vj%tau
      call propagator(vi%ukout, Elist, dtau, dmc_band, expH(:,:,itr))
      call Hpropagator(vi%ukout, Elist, mu, dtau, dmc_band, HexpH(:,:,itr))
      if(itr.eq.diagram%order)then 
         gkq(:,:,itr) = iden
      else 
         gkq(:,:,itr) = vj%gkq(:,:,vj%nu) / maxval(abs(vj%gkq(:,:,vj%nu)))
      endif
 
      !gkq(:,:,itr) = iden

      left_env(:,:,itr) = matmul(expH(:,:,itr), left_env(:,:,itr))
      left_env(:,:,itr) = matmul(gkq(:,:,itr), left_env(:,:,itr)) 
      left_env(:,:,itr+1) = left_env(:,:,itr)
      if(itr>1) then 
         if(diagram%vertexList(vi%link(2))%tau > vi%tau) then 
            nph_line(itr) = nph_line(itr-1) + 1
         else 
            nph_line(itr) = nph_line(itr-1) - 1
         endif 
      endif 
      if(jv.eq.maxN) exit 
      if(vj%link(2).eq.1)n1 = n1 + 1
      if(vj%link(2).eq.maxN)n2 = n2 + 1

      iv = vi%link(3)
      vi => diagram%vertexList(iv)
      jv = vi%link(3)
      vj => diagram%vertexList(jv)
   enddo  
   if(diagram%nph_ext.ne.n1) stop '@measure Z : nph_ext in inconsistent'
   if(diagram%nph_ext.ne.n2) stop '@measure Z : nph_ext out inconsistent'
   g_matrix = left_env(:,:,diagram%order+1)

   weight = trace(g_matrix,dmc_band)
   weight = real(weight)
   gtrue = weight / abs(real(weight))
   stat%gtrue = stat%gtrue + gtrue
   
   if(real(weight)>0.d0) then 
      diagram%gfac = 1.d0 
   else 
      diagram%gfac = -1.d0 
   endif 
   
   !Energy : order contribution 
   stat%Etrue = stat%Etrue - (diagram%order-1) * gtrue
   !Energy : phonon energy part 
   Wdt = 0.d0 
   do iv = 2,diagram%order
      vi => diagram%vertexList(iv)
      jv = vi%link(2)
      vj => diagram%vertexList(jv)
      if(jv.ne.1 .and. jv .ne. maxN) then 
         Wdt = Wdt + 0.5 * vi%wq(vi%nu) * abs(vi%tau-vj%tau) !0.5 for double counting 
      else 
         Wdt = Wdt + vi%wq(vi%nu) * abs(vi%tau-vj%tau) !no double counting for external phonon 
      endif 
   enddo 
   stat%Etrue = stat%Etrue + Wdt * gtrue
   !2. do measurement 
   right_env = iden
  
   do itr = diagram%order,1,-1 
      Hdt_matrix = matmul(right_env, HexpH(:,:,itr))
      g_matrix =  matmul(right_env, expH(:,:,itr))
      if(itr>1) then 
         Hdt_matrix = matmul(Hdt_matrix, left_env(:,:,itr-1))
         g_matrix = matmul(g_matrix, left_env(:,:,itr-1))
      endif 

      weight_2 = real(trace(g_matrix, dmc_band))
      if(abs((weight_2-weight)/weight) > 0.01d0) then 
         write(*,'(3E20.10)')abs((weight_2-weight)/weight), real(weight), real(weight_2)
         write(*,'(A20,2i6)')'itr, order = ',itr,diagram%order
         !stop 'measurement error'
      endif 
      stat%Etrue = stat%Etrue + trace(Hdt_matrix,dmc_band) /abs(weight)
     
      
      dtau = tauList(itr+1) - tauList(itr)
      tau_ratio = dtau / config%tauMax
      if(nph_line(itr).eq.0) then 
         if(maxval(abs(i_dkt(:,itr)))>1e-10) then 
            write(*,'(3i10)')i_dkt(:,itr)
            stop "momentum error @ Z-zero"
         endif 
         !back to bloch 
         !g_matrix = matmul(conjg(transpose(v_head%ukout)) ,matmul(g_matrix, v_head%ukout))
         do ib = 1,dmc_band 
            stat%Znph(ib,1) = stat%Znph(ib,1) + real(g_matrix(ib,ib)) / abs(real(weight)) * tau_ratio
         enddo 
         !stat%Znph(1,1) = stat%Znph(1,1) +  tau_ratio * real(diagram%gfac)
      else 
         !if(maxval(abs(dkt(:,itr)))<1e-15) stop "momentum error @ Z-multi"
         iph = nph_line(itr) + 1
         if(iph>100) iph = 100 
         do ib = 1,dmc_band 
            stat%Znph(ib,iph) = stat%Znph(ib,iph) + real(g_matrix(ib,ib)) / abs(real(weight)) * tau_ratio
         enddo 
         !stat%Znph(1,iph) = stat%Znph(1,iph) + tau_ratio * real(diagram%gfac)
      endif 

      right_env = matmul(right_env, expH(:,:,itr))
      if(itr.ne.1) right_env = matmul(right_env, gkq(:,:,itr-1))
   enddo 

   contains 
   function overlap(tau1,tau2,tau3,tau4)
      real(dp) :: tau1,tau2,tau3,tau4, overlap
      real(dp) :: left,right 
      left = max( tau1,tau3 )
      right = min(tau2, tau4) 
      overlap = 0.d0 
      if(left<right) overlap = right - left
   end function
   
end subroutine 

!measure Epol + polaron wfn 
! this subroutine may cause momory leaking, solution = memory allocation by hand 
subroutine measure_EZ_wfn(diagram, workspace, stat)
   
   use diagmc, only : Gel,diagmc_stat, vertex,fynman, diagmc_workspace,maxN,config, dp,cone, dmc_band,trace,get_grid_index
   use pert_param, only : print_sign, nk_svd
   implicit none 
   type(fynman),target, intent(inout) :: diagram 
   type(diagmc_stat),intent(inout) :: stat
   type(diagmc_workspace), target :: workspace
   !local 
   integer :: itau, order, ib, jb , itr, iv,jv, kv, lv, icount    
   real(dp) :: tau1,tau2, Elist(dmc_band), iden(dmc_band,dmc_band), dtau, mu, Wdt, tau_barek,  tauList(maxN)
   real(dp) :: dkt(3, maxN), tau_ratio, pt(3), qp(3)
   integer :: i_dkt(3, maxN)
   integer :: nph_line(maxN), iph, n1, n2, grid_index(3), nu_list(100,maxN), qpt_list(3,100,maxN), nph2_list(maxN)
   integer :: grid_index_el(3),  grid_index_ph(3)
   complex(dp) ::  weight, weight_2, gtrue
   type(vertex), pointer :: v_head,v_tail, vi,vj, vk, vl 
   !complex(dp), dimension(dmc_band,dmc_band,maxN) :: expH, gkq,left_env, HexpH, unk
   complex(dp), pointer :: expH(:,:,:), gkq(:,:,:),left_env(:,:,:), HexpH(:,:,:)
   complex(dp), dimension(dmc_band,dmc_band) :: right_env, Hdt_matrix, G_amp,g_matrix
   real(dp) :: denminator

   expH => workspace%expH
   gkq => workspace%gkq
   left_env => workspace%left_env
   HexpH => workspace%HexpH

   order = diagram%order/2  
   v_head => diagram%vertexList(1)
   v_tail => diagram%vertexList(maxN)
   tau_barek = 0.d0 

   !write(*,*) diagram%vertexList(maxN)%tau
   itau = floor( (v_tail%tau - config%tauMin) /config%dtau)+1
   if(itau>config%Nbin) itau = config%Nbin

   stat%Tstat(itau) = stat%Tstat(itau) + 1
   stat%Order(order+1) = stat%Order(order+1) + 1
   stat%Njjw = stat%Njjw + 1
   Elist = v_tail%ekin - minval(v_tail%ekin) 
   if(diagram%order.eq.1) then 
      g_matrix = 0.d0 
      Hdt_matrix = 0.d0 
      weight = 0.d0 
      do ib = 1,dmc_band 
         g_matrix(ib,ib) = Gel(Elist(ib), v_tail%tau)
         Hdt_matrix(ib,ib) = Gel(Elist(ib), v_tail%tau) * (v_tail%ekin(ib)) * config%tauMax
      enddo 
      weight = trace(g_matrix,dmc_band)
      stat%gtrue = stat%gtrue + 1.d0 
      stat%Etrue = stat%Etrue + trace(Hdt_matrix,dmc_band) /abs(weight)

      diagram%gfac = 1.d0 
      pt = config%i_P / real(nk_svd) 
      call get_grid_index( pt, stat%PolaronWF%ngrid, grid_index )
      do ib = 1,dmc_band 
         stat%Znph(ib,1) = stat%Znph(ib,1) + g_matrix(ib,ib) / weight 
         stat%PolaronWF%psi(ib,grid_index(1),grid_index(2),grid_index(3)) = &
            stat%PolaronWF%psi(ib,grid_index(1),grid_index(2),grid_index(3)) + real(g_matrix(ib,ib)) / abs(real(weight)) 
      enddo 
      stat%order_g(order+1) = stat%order_g(order+1) + diagram%gfac
      stat%order_sample(order+1) = stat%order_sample(order+1) + diagram%gfac
      return 
   endif 

   ! calculate deminator 
   call measure_denminator(diagram, workspace, denminator)
   weight = denminator
   diagram%sfac = denminator / abs(denminator) 
   stat%order_sample(order+1) = stat%order_sample(order+1) + denminator / abs(denminator)
   
   !Hdtau contribution 
   !1. prepare propogator along the imaginary time 
   iden = 0.d0
   do ib = 1,dmc_band; iden(ib,ib) = 1.d0; enddo 

   vi => diagram%vertexList(1)
   jv = vi%link(3)
   vj => diagram%vertexList(jv)
   left_env(:,:,1) = iden   

   nph_line = diagram%nph_ext
   nph2_list = 0
   n1 = 0; n2 = 0 
   do itr = 1,diagram%order    
      !collect phonon line 
      tau1 = 0.5*(vi%tau + vj%tau)
      icount = 0
      do kv = 2,diagram%order
         vk => diagram%vertexList(kv)
         if(vk%link(2).eq.maxN) then 
            if(vk%tau < tau1) then 
               icount = icount + 1
               if(icount>99) cycle 
               qp = -vk%i_q / real(nk_svd)
               call get_grid_index(qp, stat%PolaronLat%ngrid, grid_index)
               qpt_list(:, icount, itr) = grid_index
               nu_list(icount, itr) = vk%nu 
            endif 
         else 
            if(vk%tau < tau1) cycle !kill double couting 
            vl => diagram%vertexList(vk%link(2))
            if( vl%tau < tau1 ) then 
               icount = icount + 1
               if(icount>99) cycle 
               qp = vk%i_q / real(nk_svd) 
               call get_grid_index(qp, stat%PolaronLat%ngrid, grid_index)
               qpt_list(:,icount, itr) = grid_index
               nu_list(icount, itr) = vk%nu 
            endif 
         endif 
      enddo 
      nph2_list(itr) = icount 
      !collect el info 
      Elist = vi%ekout - minval(vi%ekout)
      mu = minval(vi%ekout)  

      i_dkt(:,itr) = mod(vi%i_kout - config%i_P, nk_svd)
      !call To1BZ(dkt(:,itr))
      dtau = vj%tau-vi%tau
      tauList(itr) = vi%tau

      call propagator(vi%ukout, Elist, dtau, dmc_band, expH(:,:,itr))
      call Hpropagator(vi%ukout, Elist, mu, dtau, dmc_band, HexpH(:,:,itr))
      if(itr.eq.diagram%order)then 
         gkq(:,:,itr) = iden
      else 
         gkq(:,:,itr) = vj%gkq_full(:,:,vj%nu) / maxval(abs(vj%gkq(:,:,vj%nu)))
      endif

      left_env(:, :, itr+1) = matmul( expH(:, :,itr), left_env(:, :, itr) )
      left_env(:, :, itr+1) = matmul( gkq(:, :, itr), left_env(:, :, itr+1) ) 
      if(itr>1) then 
         if(diagram%vertexList(vi%link(2))%tau > vi%tau) then 
            nph_line(itr) = nph_line(itr-1) + 1
         else 
            nph_line(itr) = nph_line(itr-1) - 1
         endif 
      endif 
      if(jv.eq.maxN) exit 
      if(vj%link(2).eq.1)n1 = n1 + 1
      if(vj%link(2).eq.maxN)n2 = n2 + 1

      iv = vi%link(3)
      vi => diagram%vertexList(iv)
      jv = vi%link(3)
      vj => diagram%vertexList(jv)
      if(nph2_list(itr).ne.nph_line(itr)) stop 'ph line inconsistent'
   enddo  
   tauList(diagram%order+1) = diagram%vertexList(maxN)%tau

   if(diagram%nph_ext.ne.n1) stop '@measure Z : nph_ext in inconsistent'
   if(diagram%nph_ext.ne.n2) stop '@measure Z : nph_ext out inconsistent'
   g_matrix = left_env(:,:,diagram%order+1)

   weight_2 = trace(g_matrix,dmc_band)
   weight_2 = real(weight_2)
   gtrue = weight_2 / abs(real(weight))
   stat%gtrue = stat%gtrue + gtrue
   stat%order_g(order+1) = stat%order_g(order+1) + gtrue

   diagram%gfac = gtrue

   !Energy : order contribution 
   stat%Etrue = stat%Etrue - (diagram%order-1) * gtrue
   !Energy : phonon energy part 
   Wdt = 0.d0 
   do iv = 2,diagram%order
      vi => diagram%vertexList(iv)
      jv = vi%link(2)
      vj => diagram%vertexList(jv)
      if(jv.ne.1 .and. jv .ne. maxN) then 
         Wdt = Wdt + 0.5 * vi%wq(vi%nu) * abs(vi%tau-vj%tau) !0.5 for double counting 
      else 
         Wdt = Wdt + vi%wq(vi%nu) * abs(vi%tau-vj%tau) !no double counting for external phonon 
      endif 
   enddo 
   stat%Etrue = stat%Etrue + Wdt * gtrue

   !2. do measurement 
   right_env = iden
   do itr = diagram%order,1,-1 
      Hdt_matrix = matmul(right_env, HexpH(:,:,itr))
      g_matrix =  matmul(right_env, expH(:,:,itr))
      if(itr>1) then 
         Hdt_matrix = matmul(left_env(:,:,itr), Hdt_matrix)
         g_matrix = matmul(left_env(:,:,itr), g_matrix)
      endif 
      weight_2 = trace(g_matrix,dmc_band)
      if( abs(real(weight_2)-real(weight)) > 1e-4*abs(real(weight)) )  then 
         stop 'weight inconsistent'
      endif 

      stat%Etrue = stat%Etrue + trace(Hdt_matrix,dmc_band) / abs(real(weight))

      dtau = tauList(itr+1) - tauList(itr)
      tau_ratio = dtau / config%tauMax
      if(tau_ratio<0) stop 'dtau<0'
      grid_index_el = mod( mod(i_dkt(:,itr) + config%i_P, nk_svd) + nk_svd, nk_svd ) + 1 !do it two times to ensure positivity
      if(nph_line(itr).eq.0) then 
         if(maxval(abs(i_dkt(:,itr)))>1e-10) then 
            write(*,'(3E20.10)')i_dkt(:,itr)
            stop "momentum error @ Z-zero"
         endif 
         !back to bloch
         do ib = 1,dmc_band 
            stat%Znph(ib,1) = stat%Znph(ib,1) + real(g_matrix(ib,ib)) / abs(real(weight)) * tau_ratio
            stat%PolaronWF%psi(ib,grid_index_el(1),grid_index_el(2),grid_index_el(3)) = &
               stat%PolaronWF%psi(ib,grid_index_el(1),grid_index_el(2),grid_index_el(3)) + real(g_matrix(ib,ib)) / abs(real(weight)) * tau_ratio
         enddo 

      else 
         !collect ph statistics 
         do iph = 1,nph_line(itr)
            if(iph>99) cycle 
            grid_index = qpt_list(:, iph, itr)
            stat%PolaronLat%psi(nu_list(iph,itr),grid_index(1),grid_index(2),grid_index(3) ) = & 
               stat%PolaronLat%psi(nu_list(iph,itr),grid_index(1),grid_index(2),grid_index(3) ) + diagram%gfac * tau_ratio
         enddo 
         !if(maxval(abs(dkt(:,itr)))<1e-15) stop "momentum error @ Z-multi"
         iph = nph_line(itr) + 1
         if(iph>100) iph = 100 
         do ib = 1,dmc_band 
            stat%Znph(ib,iph) = stat%Znph(ib,iph) + real(g_matrix(ib,ib)) / abs(real(weight)) * tau_ratio
            stat%PolaronWF%psi(ib,grid_index_el(1),grid_index_el(2),grid_index_el(3)) = &
               stat%PolaronWF%psi(ib,grid_index_el(1),grid_index_el(2),grid_index_el(3)) + real(g_matrix(ib,ib)) / abs(real(weight)) * tau_ratio
         enddo 
      endif 

      right_env = matmul(right_env, expH(:,:,itr))
      if(itr.ne.1) right_env = matmul(right_env, gkq(:,:,itr-1))
   enddo 

   contains 
   function overlap(tau1,tau2,tau3,tau4)
      real(dp) :: tau1,tau2,tau3,tau4, overlap
      real(dp) :: left,right 
      left = max( tau1,tau3 )
      right = min(tau2, tau4) 
      overlap = 0.d0 
      if(left<right) overlap = right - left
   end function
   
end subroutine 

subroutine measure_denminator(diagram, workspace, denminator )
   use diagmc, only : Gel,diagmc_stat, vertex,fynman, diagmc_workspace,maxN,config, dp,cone, dmc_band,trace,get_grid_index
   use pert_param, only : print_sign, nk_svd
   implicit none 
   type(fynman),target, intent(inout) :: diagram 
   type(diagmc_workspace), target :: workspace
   real(dp) :: denminator
   !local 
   integer :: itau, order, ib, jb , itr, iv,jv, kv, lv, icount    
   real(dp) :: tau1,tau2, Elist(dmc_band), iden(dmc_band,dmc_band), dtau, mu, Wdt, tau_barek,  tauList(maxN)
   real(dp) :: dkt(3, maxN), tau_ratio, pt(3), qp(3)
   integer :: i_dkt(3, maxN)
   integer :: nph_line(maxN), iph, n1, n2, grid_index(3), nu_list(100,maxN), qpt_list(3,100,maxN), nph2_list(maxN)
   integer :: grid_index_el(3),  grid_index_ph(3)
   complex(dp) ::  weight, weight_2, gtrue
   type(vertex), pointer :: v_head,v_tail, vi,vj, vk, vl 
   !complex(dp), dimension(dmc_band,dmc_band,maxN) :: expH, gkq,left_env
   complex(dp), pointer :: expH(:,:,:), gkq(:,:,:),left_env(:,:,:)
   complex(dp), dimension(dmc_band,dmc_band) :: right_env, Hdt_matrix, G_amp,g_matrix

   expH => workspace%expH
   gkq => workspace%gkq
   left_env => workspace%left_env

   if(diagram%order.eq.1) return 

   !Hdtau contribution 
   !1. prepare propogator along the imaginary time 
   iden = 0.d0
   do ib = 1,dmc_band; iden(ib,ib) = 1.d0; enddo 

   vi => diagram%vertexList(1)
   jv = vi%link(3)
   vj => diagram%vertexList(jv)
   left_env(:,:,1) = iden   

   nph_line = diagram%nph_ext
   nph2_list = 0
   n1 = 0; n2 = 0 
   do itr = 1,diagram%order    
      !collect el info 
      Elist = vi%ekout - minval(vi%ekout)
      dtau = vj%tau-vi%tau
      call propagator(vi%ukout, Elist, dtau, dmc_band, expH(:,:,itr))
      if(itr.eq.diagram%order)then 
         gkq(:,:,itr) = iden
      else 
         gkq(:,:,itr) = vj%gkq(:,:,vj%nu) / maxval(abs(vj%gkq(:,:,vj%nu)))
      endif

      left_env(:,:,itr+1) = matmul(expH(:,:,itr), left_env(:,:,itr))
      left_env(:,:,itr+1) = matmul(gkq(:,:,itr), left_env(:,:,itr+1)) 

      if(jv.eq.maxN) exit 
      iv = vi%link(3)
      vi => diagram%vertexList(iv)
      jv = vi%link(3)
      vj => diagram%vertexList(jv)
   enddo  

   g_matrix = left_env(:,:,diagram%order+1)
   weight = trace(g_matrix,dmc_band)
   denminator = real(weight)
end subroutine

!measure Epol + polaron wfn 
subroutine measure_EZ_wfn_backup(diagram, stat)
   use diagmc, only : Gel,diagmc_stat, vertex,fynman, maxN,config, dp,cone, dmc_band,trace,get_grid_index
   use pert_param, only : print_sign, nk_svd
   implicit none 
   type(fynman),target, intent(inout) :: diagram 
   type(diagmc_stat),intent(inout) :: stat
   !local 
   integer :: itau, order, ib, jb , itr, iv,jv, kv, lv, icount    
   real(dp) :: tau1,tau2, Elist(dmc_band), iden(dmc_band,dmc_band), dtau, mu, Wdt, tau_barek,  tauList(maxN)
   real(dp) :: dkt(3, maxN), tau_ratio, pt(3), qp(3)
   integer :: i_dkt(3, maxN)
   integer :: nph_line(maxN), iph, n1, n2, grid_index(3), nu_list(100,maxN), qpt_list(3,100,maxN), nph2_list(maxN)
   integer :: grid_index_el(3),  grid_index_ph(3)
   complex(dp) ::  weight, weight_2, gtrue
   type(vertex), pointer :: v_head,v_tail, vi,vj, vk, vl 
   complex(dp), dimension(dmc_band,dmc_band,maxN) :: expH, gkq,left_env, HexpH, unk
   complex(dp), dimension(dmc_band,dmc_band) :: right_env, Hdt_matrix, G_amp,g_matrix

   order = diagram%order/2  
   v_head => diagram%vertexList(1)
   v_tail => diagram%vertexList(maxN)
   tau_barek = 0.d0 

   !write(*,*) diagram%vertexList(maxN)%tau
   itau = floor( (v_tail%tau - config%tauMin) /config%dtau)+1
   if(itau>config%Nbin) itau = config%Nbin

   stat%Tstat(itau) = stat%Tstat(itau) + 1
   stat%Order(order+1) = stat%Order(order+1) + 1
   stat%Njjw = stat%Njjw + 1
   Elist = v_tail%ekin - minval(v_tail%ekin) 
   if(diagram%order.eq.1) then 
      g_matrix = 0.d0 
      Hdt_matrix = 0.d0 
      weight = 0.d0 
      do ib = 1,dmc_band 
         g_matrix(ib,ib) = Gel(Elist(ib), v_tail%tau)
         Hdt_matrix(ib,ib) = Gel(Elist(ib), v_tail%tau) * (v_tail%ekin(ib)) * config%tauMax
      enddo 
      weight = trace(g_matrix,dmc_band)
      stat%gtrue = stat%gtrue + 1.d0 
      stat%Etrue = stat%Etrue + trace(Hdt_matrix,dmc_band) /abs(weight)

      diagram%gfac = 1.d0 
      pt = config%i_P / real(nk_svd) 
      call get_grid_index( pt, stat%PolaronWF%ngrid, grid_index )
      do ib = 1,dmc_band 
         stat%Znph(ib,1) = stat%Znph(ib,1) + g_matrix(ib,ib) / weight 
         stat%PolaronWF%psi(ib,grid_index(1),grid_index(2),grid_index(3)) = &
            stat%PolaronWF%psi(ib,grid_index(1),grid_index(2),grid_index(3)) + real(g_matrix(ib,ib)) / abs(real(weight)) 
      enddo 
      stat%order_g(order+1) = stat%order_g(order+1) + diagram%gfac
      return 
   endif 

   !Hdtau contribution 
   !1. prepare propogator along the imaginary time 
   iden = 0.d0
   do ib = 1,dmc_band; iden(ib,ib) = 1.d0; enddo 

   vi => diagram%vertexList(1)
   jv = vi%link(3)
   vj => diagram%vertexList(jv)
   left_env(:,:,1) = iden   

   nph_line = diagram%nph_ext
   nph2_list = 0
   n1 = 0; n2 = 0 
   do itr = 1,diagram%order    
      !collect phonon line 
      tau1 = 0.5*(vi%tau + vj%tau)
      icount = 0
      do kv = 2,diagram%order
         vk => diagram%vertexList(kv)
         if(vk%link(2).eq.maxN) then 
            if(vk%tau < tau1) then 
               icount = icount + 1
               if(icount>99) cycle 
               qp = -vk%i_q / real(nk_svd)
               call get_grid_index(qp, stat%PolaronLat%ngrid, grid_index)
               qpt_list(:, icount, itr) = grid_index
               nu_list(icount, itr) = vk%nu 
            endif 
         else 
            if(vk%tau < tau1) cycle !kill double couting 
            vl => diagram%vertexList(vk%link(2))
            if( vl%tau < tau1 ) then 
               icount = icount + 1
               if(icount>99) cycle 
               qp = vk%i_q / real(nk_svd) 
               call get_grid_index(qp, stat%PolaronLat%ngrid, grid_index)
               qpt_list(:,icount, itr) = grid_index
               nu_list(icount, itr) = vk%nu 
            endif 
         endif 
      enddo 
      nph2_list(itr) = icount 
      !collect el info 
      Elist = vi%ekout - minval(vi%ekout)
      mu = minval(vi%ekout)  
      dkt(:,itr) = vi%kout - config%P
      i_dkt(:,itr) = mod(vi%i_kout - config%i_P, nk_svd)
      !call To1BZ(dkt(:,itr))
      dtau = vj%tau-vi%tau
      tauList(itr) = vi%tau
      tauList(itr+1) = vj%tau
      unk(:,:,itr) = vi%ukout
      call propagator(vi%ukout, Elist, dtau, dmc_band, expH(:,:,itr))
      call Hpropagator(vi%ukout, Elist, mu, dtau, dmc_band, HexpH(:,:,itr))
      if(itr.eq.diagram%order)then 
         gkq(:,:,itr) = iden
      else 
         gkq(:,:,itr) = vj%gkq(:,:,vj%nu) / maxval(abs(vj%gkq(:,:,vj%nu)))
      endif
 
      !gkq(:,:,itr) = iden

      left_env(:,:,itr) = matmul(expH(:,:,itr), left_env(:,:,itr))
      left_env(:,:,itr) = matmul(gkq(:,:,itr), left_env(:,:,itr)) 
      left_env(:,:,itr+1) = left_env(:,:,itr)
      if(itr>1) then 
         if(diagram%vertexList(vi%link(2))%tau > vi%tau) then 
            nph_line(itr) = nph_line(itr-1) + 1
         else 
            nph_line(itr) = nph_line(itr-1) - 1
         endif 
      endif 
      if(jv.eq.maxN) exit 
      if(vj%link(2).eq.1)n1 = n1 + 1
      if(vj%link(2).eq.maxN)n2 = n2 + 1

      iv = vi%link(3)
      vi => diagram%vertexList(iv)
      jv = vi%link(3)
      vj => diagram%vertexList(jv)

      if(nph2_list(itr).ne.nph_line(itr)) then 
         write(*,'(A20,2i4)') 'so many ph line : ',nph2_list(itr), nph_line(itr)
         if(nph_line(itr)<99) stop 'ph line inconsistent'
      endif 

   enddo  

   if(diagram%nph_ext.ne.n1) stop '@measure Z : nph_ext in inconsistent'
   if(diagram%nph_ext.ne.n2) stop '@measure Z : nph_ext out inconsistent'
   g_matrix = left_env(:,:,diagram%order+1)

   weight = trace(g_matrix,dmc_band)
   weight = real(weight)
   gtrue = weight / abs(real(weight))
   stat%gtrue = stat%gtrue + gtrue
   

   if(real(weight)>0.d0) then 
      diagram%gfac = 1.d0 
   else 
      diagram%gfac = -1.d0 
   endif 
   stat%order_g(order+1) = stat%order_g(order+1) + diagram%gfac
   
   !Energy : order contribution 
   stat%Etrue = stat%Etrue - (diagram%order-1) * gtrue
   !Energy : phonon energy part 
   Wdt = 0.d0 
   do iv = 2,diagram%order
      vi => diagram%vertexList(iv)
      jv = vi%link(2)
      vj => diagram%vertexList(jv)
      if(jv.ne.1 .and. jv .ne. maxN) then 
         Wdt = Wdt + 0.5 * vi%wq(vi%nu) * abs(vi%tau-vj%tau) !0.5 for double counting 
      else 
         Wdt = Wdt + vi%wq(vi%nu) * abs(vi%tau-vj%tau) !no double counting for external phonon 
      endif 
   enddo 
   stat%Etrue = stat%Etrue + Wdt * gtrue
   !2. do measurement 
   right_env = iden
  
   do itr = diagram%order,1,-1 
      Hdt_matrix = matmul(right_env, HexpH(:,:,itr))
      g_matrix =  matmul(right_env, expH(:,:,itr))
      if(itr>1) then 
         Hdt_matrix = matmul(Hdt_matrix, left_env(:,:,itr-1))
         g_matrix = matmul(g_matrix, left_env(:,:,itr-1))
      endif 

      weight_2 = real(trace(g_matrix, dmc_band))
      if(abs((weight_2-weight)/weight) > 0.01d0) then 
         write(*,'(3E20.10)')abs((weight_2-weight)/weight), real(weight), real(weight_2)
         write(*,'(A20,2i6)')'itr, order = ',itr,diagram%order
         !stop 'measurement error'
      endif 
      stat%Etrue = stat%Etrue + trace(Hdt_matrix,dmc_band) /abs(weight)
     
      
      dtau = tauList(itr+1) - tauList(itr)
      tau_ratio = dtau / config%tauMax
      
      pt = (i_dkt(:,itr) + config%i_P) / real(nk_svd)
      call get_grid_index( pt, stat%PolaronWF%ngrid, grid_index_el )

      if(nph_line(itr).eq.0) then 
         if(maxval(abs(dkt(:,itr)))>1e-10) then 
            write(*,'(3E20.10)')dkt(:,itr)
            stop "momentum error @ Z-zero"
         endif 
         !back to bloch 
         g_matrix = matmul(conjg(transpose(v_head%ukout)) ,matmul(g_matrix, v_head%ukout))
         do ib = 1,dmc_band 
            stat%Znph(ib,1) = stat%Znph(ib,1) + real(g_matrix(ib,ib)) / abs(real(weight)) * tau_ratio
            stat%PolaronWF%psi(ib,grid_index_el(1),grid_index_el(2),grid_index_el(3)) = &
               stat%PolaronWF%psi(ib,grid_index_el(1),grid_index_el(2),grid_index_el(3)) + real(g_matrix(ib,ib)) / abs(real(weight)) * tau_ratio
         enddo 
      else 
         !collect ph statistics 
         do iph = 1,nph_line(itr)
            if(iph>99) cycle 
            grid_index = qpt_list(:, iph, itr)
            stat%PolaronLat%psi(nu_list(iph,itr),grid_index(1),grid_index(2),grid_index(3) ) = & 
               stat%PolaronLat%psi(nu_list(iph,itr),grid_index(1),grid_index(2),grid_index(3) ) + diagram%gfac * tau_ratio
         enddo 
         !if(maxval(abs(dkt(:,itr)))<1e-15) stop "momentum error @ Z-multi"
         iph = nph_line(itr) + 1
         if(iph>100) iph = 100 
         !back to bloch 
         g_matrix = matmul( conjg(transpose(unk(:,:,itr))), matmul(g_matrix, unk(:,:,itr)) )
         do ib = 1,dmc_band 
            stat%Znph(ib,iph) = stat%Znph(ib,iph) + real(g_matrix(ib,ib)) / abs(real(weight)) * tau_ratio
            !if(maxval())
            stat%PolaronWF%psi(ib,grid_index_el(1),grid_index_el(2),grid_index_el(3)) = & 
               stat%PolaronWF%psi(ib,grid_index_el(1),grid_index_el(2),grid_index_el(3)) + real(g_matrix(ib,ib)) / abs(real(weight)) * tau_ratio
         enddo 
      endif 

      right_env = matmul(right_env, expH(:,:,itr))
      if(itr.ne.1) right_env = matmul(right_env, gkq(:,:,itr-1))
   enddo 

   contains 
   function overlap(tau1,tau2,tau3,tau4)
      real(dp) :: tau1,tau2,tau3,tau4, overlap
      real(dp) :: left,right 
      left = max( tau1,tau3 )
      right = min(tau2, tau4) 
      overlap = 0.d0 
      if(left<right) overlap = right - left
   end function
   
end subroutine 

subroutine measure_EZ_toy(diagram, stat)
   use diagmc, only : Gel,diagmc_stat, vertex,fynman, maxN,config, dp,cone, dmc_band,trace
   use pert_param, only : print_sign
   implicit none 
   type(fynman),target, intent(inout) :: diagram 
   type(diagmc_stat),intent(inout) :: stat
   !local 
   integer :: itau, order, ib, jb , itr, iv,jv   
   real(dp) :: tau1,tau2, Elist(dmc_band), iden(dmc_band,dmc_band), dtau, mu, Wdt, tau_barek,  tauList(maxN)
   real(dp) :: dkt(3, maxN), tau_ratio
   integer :: nph_line(maxN), iph, n1,n2
   complex(dp) ::  weight, weight_2, gtrue  
   type(vertex), pointer :: v_head,v_tail, vi,vj
   complex(dp), dimension(dmc_band,dmc_band,maxN) :: expH, gkq,left_env, HexpH
   complex(dp), dimension(dmc_band,dmc_band) :: right_env, Hdt_matrix, G_amp,g_matrix

   order = diagram%order/2  
   v_head => diagram%vertexList(1)
   v_tail => diagram%vertexList(maxN)
   tau_barek = 0.d0 

   !write(*,*) diagram%vertexList(maxN)%tau
   itau = floor( (v_tail%tau - config%tauMin) /config%dtau)+1
   if(itau>config%Nbin) itau = config%Nbin

   stat%Tstat(itau) = stat%Tstat(itau) + 1
   stat%Order(order+1) = stat%Order(order+1) + 1
   stat%Njjw = stat%Njjw + 1
   Elist = v_tail%ekin - minval(v_tail%ekin) 
   if(diagram%order.eq.1) then 
      stat%gtrue = stat%gtrue + 1.d0 
      stat%Etrue = stat%Etrue + v_tail%ekin(1) * config%tauMax
      diagram%gfac = 1.d0 
      do ib = 1,dmc_band 
         stat%Znph(ib,1) = stat%Znph(ib,1) + 1.d0 
      enddo 
      return 
   endif 

  
   
   !Hdtau contribution 
   !1. prepare propogator along the imaginary time 
   vi => diagram%vertexList(1)
   jv = vi%link(3)
   vj => diagram%vertexList(jv)

   nph_line = diagram%nph_ext
   n1 = 0; n2 = 0 
   do itr = 1,diagram%order    
      dtau = vj%tau-vi%tau
      tauList(itr) = vi%tau
      tauList(itr+1) = vj%tau
     
      stat%Etrue = stat%Etrue + vi%ekout(1) * dtau
      if(itr>1) then 
         if(diagram%vertexList(vi%link(2))%tau > vi%tau) then 
            nph_line(itr) = nph_line(itr-1) + 1
         else 
            nph_line(itr) = nph_line(itr-1) - 1
         endif 
      endif 
      if(jv.eq.maxN) exit 
      if(vj%link(2).eq.1)n1 = n1 + 1
      if(vj%link(2).eq.maxN)n2 = n2 + 1

      iv = vi%link(3)
      vi => diagram%vertexList(iv)
      jv = vi%link(3)
      vj => diagram%vertexList(jv)
   enddo  
   if(diagram%nph_ext.ne.n1) stop '@measure Z : nph_ext in inconsistent'
   if(diagram%nph_ext.ne.n2) stop '@measure Z : nph_ext out inconsistent'

   stat%gtrue = stat%gtrue + 1.d0 
   diagram%gfac = 1.d0 

   !Energy : order contribution 
   stat%Etrue = stat%Etrue - (diagram%order-1) 
   !Energy : phonon energy part 
   Wdt = 0.d0 
   do iv = 2,diagram%order
      vi => diagram%vertexList(iv)
      jv = vi%link(2)
      vj => diagram%vertexList(jv)
      if(jv.ne.1 .and. jv .ne. maxN) then 
         Wdt = Wdt + 0.5 * vi%wq(vi%nu) * abs(vi%tau-vj%tau) !0.5 for double counting 
      else 
         Wdt = Wdt + vi%wq(vi%nu) * abs(vi%tau-vj%tau) !no double counting for external phonon 
      endif 
   enddo 
   stat%Etrue = stat%Etrue + Wdt 

   !2. do measurement 
   do itr = diagram%order,1,-1 
      dtau = tauList(itr+1) - tauList(itr)
      tau_ratio = dtau / config%tauMax
      if(nph_line(itr).eq.0) then 
         if(maxval(abs(dkt(:,itr)))>1e-10) then 
            write(*,'(3E20.10)')dkt(:,itr)
            stop "momentum error @ Z-zero"
         endif 
         !back to bloch 
         do ib = 1,dmc_band 
            stat%Znph(ib,1) = stat%Znph(ib,1) +  tau_ratio
         enddo 
      else 
         iph = nph_line(itr) + 1
         if(iph>100) iph = 100 
         do ib = 1,dmc_band 
            stat%Znph(ib,iph) = stat%Znph(ib,iph) +  tau_ratio
         enddo 
      endif 
   enddo 

   contains 
   function overlap(tau1,tau2,tau3,tau4)
      real(dp) :: tau1,tau2,tau3,tau4, overlap
      real(dp) :: left,right 
      left = max( tau1,tau3 )
      right = min(tau2, tau4) 
      overlap = 0.d0 
      if(left<right) overlap = right - left
   end function
   
end subroutine 

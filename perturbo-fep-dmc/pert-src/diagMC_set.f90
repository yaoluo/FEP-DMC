subroutine diagmc_SEt()
   USE OMP_LIB
   use DiagMC 
   use diagMC_set_debug
   use selfE_update, only : update_drive, init_selfe_diagram
   implicit none 
   !omp shared variables 
   integer,parameter ::  Ntype = 8
   integer :: Nthread, iunit

   !omp private variables 
   integer :: i_omp, i_tmp, j_tmp, Nwork  
   integer :: i, update 
   type(diagmc_stat) :: stat(200), stat_total      
   type(fynman) :: diagrams(200) 

   real(dp) :: PA(Ntype), process 
   integer :: finished = 0  
   !!$omp threadprivate(diagram,stat)
   multiband = .false. 
   call start_clock('diagmc_SEt')

   call setup_dqmc()
   !do sanity checking 
   if(se_check.eq.1) then 
      call ReSelfenergy()
      !call selfenergy_fan_time()
      !call ReSelfenergy_degenerate()
      return 
   end if 
   call check_epw()

   
   !!$omp parallel default(shared)
   Nthread = OMP_GET_MAX_THREADS()
   !!$omp end parallel 

   write(stdout,'(A40,f6.3,2i6)')'DiagMC-SelfE start : dtau, Nthread, maxOrder = ',config%dtau,Nthread, (maxOrder-1)/2
   call Choose_FakeRatio()
   !FakeRatio = 1e-5
   call init_observable(stat_total) 
   !--------------------------------------------------------------------
   !------------------- below is openmp paralleled ---------------------
   !---------------------------------------------------------------------
   
   do i_omp = 1,Nthread
      call InitVertexList(maxN, Nph, dmc_band,,nbnd, diagrams(i_omp)%vertexList)
      call init_selfe_diagram(diagrams(i_omp))
      call init_observable(stat(i_omp))
   enddo
   config%PA = config%PA/sum(config%PA)

   call start_clock('mcmc')
   !$omp parallel default(shared) firstprivate(i_omp,i,update,PA,process,Nwork) 
      Nwork = config%Nmcmc
      i_omp = omp_get_thread_num() + 1
      call init_rand_seed_omp(i_omp, diagrams(i_omp)%seed)
      PA = config%PA(1:Ntype)
      do i = 1,config%Nmcmc
         ![1] select update type 
         call Discrete_omp(diagrams(i_omp)%seed,PA,Ntype,update)

         ![2] proposal new configuration according to the update type and calculate the acceptance ratio 
         call update_drive(diagrams(i_omp),stat(i_omp),update)
      
         ![3] print & check diagram 
         if(mod(i,10**5).eq.0.and.i_omp.eq.1 )  then 
            write(*,'(A50,3i10,2E15.5)')' step, order, update , Tau, Re(gfac) = ',&
            i,diagrams(i_omp)%order, update, diagrams(i_omp)%vertexList(diagrams(i_omp)%tail)%tau,real(diagrams(i_omp)%gfac)
         endif 
         if(mod(i,10**5).eq.0) then 
            call check_diagram(diagrams(i_omp))
         endif
         if(mod(i,10**5).eq.0) then 
            call check_phase(diagrams(i_omp))
         end if 
         ![4] gather data 
         call measureSE(diagrams(i_omp),stat(i_omp))
         
         if(mod(i,10**5).and.finished.eq.1) then  !if one core is finished, this core stops but remember how much I did
            Nwork = i; exit 
         endif 
      end do
      !omp gather, need implement  
      stat(i_omp)%GreenFunc = stat(i_omp)%GreenFunc / sum(abs(stat(i_omp)%GreenFunc)/config%Nbin)
      process = real(Nwork,dp)/config%Nmcmc
      !$omp critical 
         stat_total%accept = stat_total%accept + stat(i_omp)%accept * process
         stat_total%update = stat_total%update + stat(i_omp)%update * process
         stat_total%Tstat = stat_total%Tstat + stat(i_omp)%Tstat * process 
         stat_total%order = stat_total%order + stat(i_omp)%order * process 
         stat_total%order_g = stat_total%order_g + stat(i_omp)%order_g * process 
         stat_total%GreenFunc = stat_total%GreenFunc + stat(i_omp)%GreenFunc * process 
      !$omp end critical 
      finished = 1
   !$omp end parallel

   call stop_clock('mcmc')
   write(stdout,'(A30,i15,A7)')'total work = ',sum(stat_total%Tstat)/10000,'x1e4'
   write(stdout,'(A30,7i15)')'update try = ',stat_total.update
   write(stdout,'(A30,7E15.5)')'acceptance = ',(1.d0*stat_total.accept)/stat_total.update
   call output_gt(stat_total)
   
   call stop_clock('diagmc_SEt')
   if(ionode) then 
      iunit = find_free_unit()
      open(iunit,file='parameter.dat')
      write(iunit,'(A10)') '#FakeRatio = '
      WRITE(iunit,'(E15.5)')FakeRatio
      write(iunit,'(A10)') '#alpha = '
      WRITE(iunit,'(E15.5)')alpha_frolich 
      write(iunit,'(A10)') '#Eshift = '
      WRITE(iunit,'(E15.5)')config%Eshift 
      write(iunit,'(A10)') '#Ep = '
      WRITE(iunit,'(E15.5)')config%EP + config%mu
      close(iunit)
   endif 
end subroutine 

subroutine diagmc_SEt_multiband()
   USE OMP_LIB
   use DiagMC 
   use diagMC_set_debug
   use selfE_update_multiband, only : update_drive, init_selfe_diagram
   implicit none 
   !omp shared variables 
   integer,parameter ::  Ntype = 8
   integer :: Nthread, iunit

   !omp private variables 
   integer :: i_omp, i_tmp, j_tmp, Nwork  
   integer :: i, update 
   type(diagmc_stat) :: stat(200), stat_total      
   type(fynman) :: diagrams(200) 

   real(dp) :: PA(Ntype), process, target_ratio  
   integer :: finished = 0  
   !!$omp threadprivate(diagram,stat)

   multiband = .true. 
   call start_clock('diagmc_SEt_multi')

   call setup_dqmc()
   !do sanity checking 
   

   select case(se_check)
      case(1)
         call ReSelfenergy()
      case(2)
         call ReSelfenergy_degenerate()
      case(3)
         call selfenergy_fan_time_matrix()
      case(4)
         call MataSelfEnergy_Matrix()
   end select 
   if(se_check.ne.0) return  
   
   call check_epw()

   !!$omp parallel default(shared)
   Nthread = OMP_GET_MAX_THREADS()
   !!$omp end parallel 

   write(stdout,'(A40,f6.3,2i6)')'DiagMC-SelfE start : dtau, Nthread, maxOrder = ',config%dtau,Nthread, (maxOrder-1)/2
   target_ratio = 0.5 
   call Choose_FakeRatio(target_ratio)
   call init_observable(stat_total) 
   !--------------------------------------------------------------------
   !------------------- below is openmp paralleled ---------------------
   !---------------------------------------------------------------------
   
   do i_omp = 1,Nthread
      call InitVertexList(maxN, Nph, dmc_band,nbnd, diagrams(i_omp)%vertexList)
      call init_selfe_diagram(diagrams(i_omp))
      call init_observable(stat(i_omp))
   enddo
   config%PA = config%PA/sum(config%PA)

   call start_clock('mcmc')
   !$omp parallel default(shared) firstprivate(i_omp,i,update,PA,process,Nwork) 
      Nwork = config%Nmcmc
      i_omp = omp_get_thread_num() + 1
      call init_rand_seed_omp(i_omp, diagrams(i_omp)%seed)
      PA = config%PA(1:Ntype)
      do i = 1,config%Nmcmc
         ![1] select update type 
         call Discrete_omp(diagrams(i_omp)%seed,PA,Ntype,update)

         ![2] proposal new configuration according to the update type and calculate the acceptance ratio 
         call update_drive(diagrams(i_omp),stat(i_omp),update)
      
         ![3] print & check diagram 
         if(mod(i,10**5).eq.0.and.i_omp.eq.1 )  then 
            write(*,'(A50,3i10,2E15.5)')' step, order, update , Tau, Re(gfac) = ',&
            i,diagrams(i_omp)%order, update, diagrams(i_omp)%vertexList(diagrams(i_omp)%tail)%tau,real(diagrams(i_omp)%gfac)
         endif 
         if(mod(i,10**5).eq.0) then 
            call check_diagram(diagrams(i_omp))
         endif
         if(mod(i,10**5).eq.0) then 
            call check_phase(diagrams(i_omp))
         end if 
         ![4] gather data 
         call measureSE_multiband(diagrams(i_omp),stat(i_omp))
         
         if(mod(i,10**5).and.finished.eq.1) then  !if one core is finished, this core stops but remember how much I did
            Nwork = i; exit 
         endif 
      end do
      !omp gather, need implement  
      stat(i_omp)%GreenFunc_matrix = stat(i_omp)%GreenFunc_matrix / sum(abs(stat(i_omp)%GreenFunc_matrix)/config%Nbin)
      process = real(Nwork,dp)/config%Nmcmc
      !$omp critical 
         stat_total%accept = stat_total%accept + stat(i_omp)%accept * process
         stat_total%update = stat_total%update + stat(i_omp)%update * process
         stat_total%Tstat = stat_total%Tstat + stat(i_omp)%Tstat * process 
         stat_total%order = stat_total%order + stat(i_omp)%order * process 
         stat_total%order_g = stat_total%order_g + stat(i_omp)%order_g * process 
         stat_total%GreenFunc_matrix = stat_total%GreenFunc_matrix + stat(i_omp)%GreenFunc_matrix * process 
      !$omp end critical 
      finished = 1
   !$omp end parallel

   call stop_clock('mcmc')
   write(stdout,'(A30,i15,A7)')'total work = ',sum(stat_total%Tstat)/10000,'x1e4'
   write(stdout,'(A30,7i15)')'update try = ',stat_total.update
   write(stdout,'(A30,7E15.5)')'acceptance = ',(1.d0*stat_total.accept)/stat_total.update
   call NormalizeSigma(stat_total)
   call output_gt(stat_total)
   
   call stop_clock('diagmc_SEt_multi')
   if(ionode) then 
      iunit = find_free_unit()
      open(iunit,file='parameter.dat')
      write(iunit,'(A10)') '#FakeRatio = '
      WRITE(iunit,'(E15.5)')FakeRatio 
      write(iunit,'(A10)') '#alpha = '
      if(Holstein) then 
         WRITE(iunit,'(E15.5)')alpha_Holstein
      else 
         WRITE(iunit,'(E15.5)')alpha_frolich 
      endif
      write(iunit,'(A10)') '#mu = '
      WRITE(iunit,'(E15.5)')config%mu 
      write(iunit,'(A10)') '#Ep = '
      WRITE(iunit,'(E15.5)')config%EP + config%mu
      close(iunit)
   endif 
   !WRITE(*,*) config%P
   
end subroutine 

subroutine diagmc_SEt_DMFT()
   USE OMP_LIB
   use DiagMC 
   use diagMC_set_debug
   use selfE_update_DMFT, only : update_drive, init_selfe_diagram
   use HamiltonianTable, only : FineK
   implicit none 
   !omp shared variables 
   integer,parameter ::  Ntype = 8
   integer :: Nthread, iunit

   !omp private variables 
   integer :: i_omp, i_tmp, j_tmp, Nwork  
   integer :: i, update 
   type(diagmc_stat) :: stat(200), stat_total      
   type(fynman) :: diagrams(200) 

   real(dp) :: PA(Ntype), process, target_ratio, tmp   
   integer :: finished = 0, maxik(1)  
   !!$omp threadprivate(diagram,stat)

   multiband = .true. 
   call start_clock('diagmc_SEt_DMFT')

   call setup_dqmc()
   !do sanity checking 
   if(se_check.eq.1) then 
      !call ReSelfenergy()
      call selfenergy_fan_time()
      !call ReSelfenergy_degenerate()
      return 
   end if 
   call check_epw()

   !!$omp parallel default(shared)
   Nthread = OMP_GET_MAX_THREADS()
   !!$omp end parallel 

   write(stdout,'(A40,f6.3,2i6)')'DiagMC-SelfE start : dtau, Nthread, maxOrder = ',config%dtau,Nthread, (maxOrder-1)/2
   target_ratio = 1.0 
   call Choose_FakeRatio(target_ratio)
   call init_observable(stat_total) 
   !--------------------------------------------------------------------
   !------------------- below is openmp paralleled ---------------------
   !---------------------------------------------------------------------
   
   do i_omp = 1,Nthread
      call InitVertexList(maxN, Nph, dmc_band,nbnd, diagrams(i_omp)%vertexList)
      call init_selfe_diagram(diagrams(i_omp))
      call init_observable(stat(i_omp))
   enddo
   config%PA = config%PA/sum(config%PA)

   call start_clock('mcmc')
   !$omp parallel default(shared) firstprivate(i_omp,i,update,PA,process,Nwork,tmp) 
      Nwork = config%Nmcmc
      i_omp = omp_get_thread_num() + 1
      call init_rand_seed_omp(i_omp, diagrams(i_omp)%seed)
      PA = config%PA(1:Ntype)
      do i = 1,config%Nmcmc
         ![1] select update type 
         call Discrete_omp(diagrams(i_omp)%seed,PA,Ntype,update)

         ![2] proposal new configuration according to the update type and calculate the acceptance ratio 
         call update_drive(diagrams(i_omp),stat(i_omp),update)
      
         ![3] print & check diagram 
         if(mod(i,10**5).eq.0.and.i_omp.eq.1 )  then 
            write(*,'(A50,3i10,2E15.5)')' step, order, update , Tau, Re(gfac) = ',&
            i,diagrams(i_omp)%order, update, diagrams(i_omp)%vertexList(diagrams(i_omp)%tail)%tau,real(diagrams(i_omp)%gfac)
         endif 
         if(mod(i,10**5).eq.0) then 
            call check_diagram(diagrams(i_omp))
         endif
         if(mod(i,10**5).eq.0) then 
            call check_phase(diagrams(i_omp))
         end if 
         ![4] gather data 
         call measureSE_DMFT(diagrams(i_omp),stat(i_omp))
         
         if(mod(i,10**5).and.finished.eq.1) then  !if one core is finished, this core stops but remember how much I did
            Nwork = i; exit 
         endif 
      end do
      !omp gather, need implement  
      tmp = 1.d0 / sum(abs(stat(i_omp)%GreenFunc_matrix)/config%Nbin)
      stat(i_omp)%GreenFunc_matrix = stat(i_omp)%GreenFunc_matrix * tmp
      stat(i_omp)%SigmaR = stat(i_omp)%SigmaR * tmp
      process = real(Nwork,dp)/config%Nmcmc
      !$omp critical 
         stat_total%accept = stat_total%accept + stat(i_omp)%accept * process
         stat_total%update = stat_total%update + stat(i_omp)%update * process
         stat_total%Tstat = stat_total%Tstat + stat(i_omp)%Tstat * process 
         stat_total%order = stat_total%order + stat(i_omp)%order * process 
         stat_total%order_g = stat_total%order_g + stat(i_omp)%order_g * process 
         stat_total%GreenFunc_matrix = stat_total%GreenFunc_matrix + stat(i_omp)%GreenFunc_matrix * process 
         stat_total%SigmaR = stat_total%SigmaR + stat(i_omp)%SigmaR * process 
         stat_total%Pstat = stat_total%Pstat + stat(i_omp)%Pstat * process 
      !$omp end critical 
      finished = 1
   !$omp end parallel
      
   call print_vertex(diagrams(1))
   call stop_clock('mcmc')
   write(stdout,'(A30,i15,A7)')'total work = ',sum(stat_total%Tstat)/10000,'x1e4'
   write(stdout,'(A30,7i15)')'update try = ',stat_total.update
   write(stdout,'(A30,7E15.5)')'acceptance = ',(1.d0*stat_total.accept)/stat_total.update
   call NormalizeSigma(stat_total)
   
   call output_gt(stat_total)
   call output_SigmaR(stat_total)

   maxik = maxloc(stat_total%Pstat)
   write(*,'(A20,3f10.5)') 'Hot K = ',FineK(:,maxik(1))

   call stop_clock('diagmc_SEt_DMFT')
   if(ionode) then 
      iunit = find_free_unit()
      open(iunit,file='parameter.dat')
      write(iunit,'(A10)') '#FakeRatio = '
      WRITE(iunit,'(E15.5)')FakeRatio 
      write(iunit,'(A10)') '#alpha = '
      WRITE(iunit,'(E15.5)')alpha_frolich 
      write(iunit,'(A10)') '#mu = '
      WRITE(iunit,'(E15.5)')config%mu 
      write(iunit,'(A10)') '#Ep = '
      WRITE(iunit,'(E15.5)')config%EP + config%mu
      close(iunit)

      open(iunit,file='kstat.dat')
      write(iunit,*) stat_total%Pstat
      close(iunit)
   endif 
end subroutine 


subroutine Choose_FakeRatio(target_ratio_in)
   USE OMP_LIB
   use DiagMC 
   use diagMC_set_debug
   use selfE_update_multiband, only : update_drive, init_selfe_diagram

   implicit none 
   !omp shared variables
   real(dp),optional :: target_ratio_in
   integer,parameter ::  Ntype = 8
   integer :: Nthread
   real(dp), parameter :: target_ratio_default = 0.3
   !omp private variables 
   integer :: i_omp, i_tmp, j_tmp, Nwork  
   integer :: i, update 
   type(diagmc_stat) :: stat(200), stat_total      
   type(fynman) :: diagrams(200) 

   real(dp) :: PA(Ntype), process, Ratio ,target_ratio
   integer :: finished = 0, iunit  
   !!$omp threadprivate(diagram,stat)

   target_ratio = target_ratio_default
   if(present(target_ratio_in))target_ratio = target_ratio_in
   call start_clock('FakeRatio')

   
   !!$omp parallel default(shared)
   Nthread = OMP_GET_MAX_THREADS()
   !!$omp end parallel 


   write(stdout,'(A40)')'Measure FakeRatio ----- '

   call init_observable(stat_total) 
   !--------------------------------------------------------------------
   !------------------- below is openmp paralleled ---------------------
   !---------------------------------------------------------------------
   
   do i_omp = 1,Nthread
      call InitVertexList(maxN, Nph, dmc_band,nbnd, diagrams(i_omp)%vertexList)
      call init_selfe_diagram(diagrams(i_omp))
      call init_observable(stat(i_omp))
   enddo
   config%PA = config%PA/sum(config%PA)

   call start_clock('mcmc')
   !$omp parallel default(shared) firstprivate(i_omp,i,update,PA,process,Nwork) 
      Nwork = config%Nmcmc
      i_omp = omp_get_thread_num() + 1
      call init_rand_seed_omp(i_omp, diagrams(i_omp)%seed)
      PA = config%PA(1:Ntype)
      PA(1) = 0; PA(2) = 0; PA(5) = 0; !turn off add/remove/move internal vertex: so that only sample between n = 0/1
      PA = PA/sum(PA)
      do i = 1,config%Nmcmc/100
         ![1] select update type 
         call Discrete_omp(diagrams(i_omp)%seed,PA,Ntype,update)

         ![2] proposal new configuration according to the update type and calculate the acceptance ratio 
         call update_drive(diagrams(i_omp),stat(i_omp),update)
      
         ![4] gather data 
         call measureSE(diagrams(i_omp),stat(i_omp))
         
         if(mod(i,10**5).and.finished.eq.1) then  !if one core is finished, this core stops but remember how much I did
            Nwork = i; exit 
         endif 
      end do
      !omp gather, need implement  
      process = real(Nwork,dp)/config%Nmcmc
      !$omp critical 
         stat_total%order = stat_total%order + stat(i_omp)%order * process 
      !$omp end critical 
      finished = 1
   !$omp end parallel
   stat_total%order(1) = int(stat_total%order(1) / FakeRatio)
   Ratio = real(stat_total%order(1)) / stat_total%order(2)
   !target_ratio = FakeRatio*order(1)/order(2) 
   FakeRatio = target_ratio/ratio 
   if(multiband) FakeRatio = FakeRatio 
   call stop_clock('FakeRatio')


end subroutine 

!measure 
subroutine measureSE(diagram,stat)
   use diagmc, only : dp, diagmc_stat, fynman, maxN, config
   implicit none 
   type(diagmc_stat) :: stat
   type(fynman) :: diagram 
   integer :: itau,iv1,iv2,order  
   real(dp) :: tau

   order = diagram%order/2
   stat%Order(order+1) = stat%Order(order + 1) + 1
   stat%order_g(order+1) = stat%order_g(order + 1) + diagram%gfac
   if(order.eq.0) return 
   tau = diagram%vertexList(diagram%tail)%tau 
   itau = floor( tau/config%dtau )+1
   if(itau>config%Nbin) itau = config%Nbin
   stat%Tstat(itau) = stat%Tstat(itau) + 1
   stat%GreenFunc(itau,order) = stat%GreenFunc(itau,order) + diagram%gfac
end subroutine 

!measure 
subroutine measureSE_multiband(diagram,stat)
   use pert_param, only : lowest_order
   use diagmc, only : dp, diagmc_stat, fynman, maxN, config
   use HamiltonianTable, only : IndexK,Nkxyz
   implicit none 
   type(diagmc_stat) :: stat
   type(fynman) :: diagram 
   integer :: itau,iv1,iv2,order, m,n, ik  
   real(dp) :: tau, K_ex(3)
   
   !sigma_{mn}
   m = diagram%vertexList(diagram%tail)%n_out 
   n = diagram%vertexList(diagram%head)%n_in 

   order = diagram%order/2
   stat%Order(order+1) = stat%Order(order + 1) + 1
   stat%order_g(order+1) = stat%order_g(order + 1) + diagram%gfac
   
   if(order.eq.0) return 

   tau = diagram%vertexList(diagram%tail)%tau 
   itau = floor( tau/config%dtau )+1
   if(itau>config%Nbin) itau = config%Nbin
   stat%Tstat(itau) = stat%Tstat(itau) + 1
   stat%GreenFunc_matrix(itau,order,m,n) = stat%GreenFunc_matrix(itau,order,m,n) + diagram%gfac

   !DMFT measure
   K_ex = diagram%vertexList(diagram%tail)%kout  
   ik = IndexK(K_ex,Nkxyz)
   stat%Pstat(ik) = stat%Pstat(ik) + 1

end subroutine 

!
subroutine measureSE_DMFT(diagram,stat)
   use pert_param, only : lowest_order
   use diagmc, only : dp, diagmc_stat, fynman, maxN, config
   use HamiltonianTable, only : IndexK,Nkxyz
   use pert_const, only : ci,twopi 
   implicit none 
   type(diagmc_stat) :: stat
   type(fynman) :: diagram 
   integer :: itau,iv1,iv2,order, m,n, ik, ir, isym   
   real(dp) :: tau, K_ex(3)
   real(dp) :: Rot90(3,3)
   !sigma_{mn}
   m = diagram%vertexList(diagram%tail)%n_out 
   n = diagram%vertexList(diagram%head)%n_in 

   order = diagram%order/2
   stat%Order(order+1) = stat%Order(order + 1) + 1
   stat%order_g(order+1) = stat%order_g(order + 1) + diagram%gfac
   
   if(order<lowest_order) return 

   tau = diagram%vertexList(diagram%tail)%tau 
   itau = floor( tau/config%dtau )+1
   if(itau>config%Nbin) itau = config%Nbin
   stat%Tstat(itau) = stat%Tstat(itau) + 1
   stat%GreenFunc_matrix(itau,order,m,n) = stat%GreenFunc_matrix(itau,order,m,n) + diagram%gfac

   !DMFT measure
   K_ex = diagram%vertexList(diagram%tail)%kout  
   ik = IndexK(K_ex,Nkxyz)
   stat%Pstat(ik) = stat%Pstat(ik) + 1
   Rot90 = 0;Rot90(1,2) = -1; Rot90(2,1) = 1; Rot90(3,3) = 1
   !measure sigma_R 
   do ir = 1,stat%NR 
      do isym = 1,4 
         K_ex = matmul(Rot90, K_ex)
         stat%SigmaR(itau,m,n,ir) = stat%SigmaR(itau,m,n,ir) & 
                           + diagram%gfac * exp(ci*twopi* sum( stat%Rmesh(:,ir)*K_ex) )
      enddo 
   enddo 

   

end subroutine 

!
subroutine NormalizeSigma(stat)
   use pert_param, only : lowest_order
   use diagmc, only : dp, diagmc_stat, config,FakeRatio,multiband
   implicit none 
   type(diagmc_stat) :: stat
   real(dp) :: I0,P0, sumSigma,dtau, temp  

   I0 = FakeRatio * config%tauMax
   if(lowest_order<=1)  P0 = real(stat%order(1)) / (sum(stat%order)-stat%order(1))
   if(lowest_order.eq.2) P0 = real(stat%order(1)) / ( sum(stat%order)-sum(stat%order(1:2)) )

   sumSigma = I0/P0 
   dtau = config%dtau

   if(multiband) then 
      temp = 1.d0 / sum(stat%GreenFunc_matrix*dtau) * sumSigma
      stat%GreenFunc_matrix = stat%GreenFunc_matrix * temp
      stat%SigmaR = stat%SigmaR * temp
   else 
      stat%GreenFunc = stat%GreenFunc / sum(stat%GreenFunc*dtau) * sumSigma
   endif 

end subroutine


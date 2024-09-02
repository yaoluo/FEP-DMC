subroutine diagmc_JJt()
   
   USE OMP_LIB
   use DiagMC 
   use diagMC_gt_debug
   use JJt_update, only : update_drive
   implicit none 
   !omp shared variables 
   integer,parameter ::  Ntype = 6
   integer :: Nthread,iunit

   !omp private variables 
   integer :: i_omp, i_tmp, j_tmp, Nwork  
   integer :: i, update 
   type(diagmc_stat) :: stat(200), stat_total      
   type(fynman) :: diagrams(200) 

   real(dp) :: PA(Ntype), xr, process 
   integer :: finished = 0  
   !!$omp threadprivate(diagram,stat)
   call start_clock('diagmc_gt')

   call setup_dqmc()
   if(.not.LatticeFrohlich .and. .not.Holstein) call check_epw()

   
   !!$omp parallel default(shared)
   Nthread = OMP_GET_MAX_THREADS()
   !!$omp end parallel 

   write(stdout,'(A40,f6.3,2i6)')'DiagMC start : dtau, Nthread, maxOrder = ',config%dtau,Nthread, maxOrder

   call init_observable(stat_total) 
   !--------------------------------------------------------------------
   !------------------- below is openmp paralleled ---------------------
   !---------------------------------------------------------------------
   
   do i_omp = 1,Nthread
      call InitVertexList(maxN, Nph,dmc_band, nbnd,nsvd, diagrams(i_omp)%vertexList)
      call init_diagram_JJ(diagrams(i_omp))
      call init_observable(stat(i_omp))
   enddo
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
            write(*,'(A40,2i10,2E15.5)')' step, order , Tau, Re(gfac) = ',&
                     i,diagrams(i_omp)%order, diagrams(i_omp)%vertexList(maxN)%tau,real(diagrams(i_omp)%gfac)
         endif 
         if(mod(i,10**5).eq.0) then 
            call check_diagram(diagrams(i_omp))
         endif
         if(mod(i,10**5).eq.0) then 
            call check_phase(diagrams(i_omp))
         end if 
         ![4] gather data 
         call measureG(diagrams(i_omp),stat(i_omp))
         
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
         WRITE(iunit,'(E15.5)')alpha_frolich 
      endif 
      write(iunit,'(A10)') '#Ep = '
      WRITE(iunit,'(E15.5)')config%EP + config%mu
      close(iunit)
   endif 

end subroutine 

!init self-energy diagram, start form the fake one. 
subroutine init_diagram_JJ(diagram)
   use diagMC 
      implicit none 
      type(fynman), target :: diagram
      type(vertex), pointer :: vi

  
   
      !init diagram 
      diagram%order = 1; 
      diagram%gfac = 1.d0; 
      vi => diagram%vertexList(1)
      vi%kout = config%P;  vi%ekout = config%EP; vi%tau = 0.d0
      vi%link(3) = maxN
      call cal_ek(vi%kout,vi%ekout,vi%ukout)
      vi%n_out = config%ib 

      vi => diagram%vertexList(maxN)
      vi%kin = config%P; vi%ekin = config%EP; vi%tau = config%tauMax
      vi%link(1) = 1
      call cal_ek(vi%kin,vi%ekin,vi%ukin) 
      vi%n_in = config%ib 
end subroutine


!
subroutine measureJJ(diagram,stat)
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

   order = (diagram%order-3)/2
   stat%Order(order+1) = stat%Order(order + 1) + 1
   stat%order_g(order+1) = stat%order_g(order + 1) + diagram%gfac
   
   if(order<lowest_order) return 

   tau = diagram%vertexList(diagram%tail)%tau 
   itau = floor( tau/config%dtau )+1
   if(itau>config%Nbin) itau = config%Nbin
   stat%Tstat(itau) = stat%Tstat(itau) + 1
   stat%GreenFunc_matrix(itau,order,m,m) = stat%GreenFunc_matrix(itau,order,m,m) + diagram%gfac

end subroutine 

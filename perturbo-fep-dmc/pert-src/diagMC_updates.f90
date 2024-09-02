!------------------------------------------------------------------------
!diagrammatic monte carlo for electron phonon Hamiltonian
!key subroutines are the 7 updates: 
!1. add_ph : add phonon line on one electron line 
!2. remove_ph : the complementary update of add_ph
!3. change_ph_mode: change the phonon mode 
!4. change_band : change internal band index 
!5. update_vertex_time : move the time of e-ph vertex 
!6. update_tau: change the external time 
!7. swap: swap nearby two vertex in time 
!author : Yao Luo, Email: yluo7@caltech.edu
!
!add_ph() & swap() are the most expensive (other updates cost no time in comparsion)  
! because it needs to evaluate new e-ph matrix. 
!For efficiency, the prob for the two updates should be comparablely lower than the others. 
!Version 2.0 @ April 2022 : 
!     changed every updates with fynman diagrams and statistics as input for openmp parallelization. 
!Version 3.0 
!     matrix version for multiband system. the band index of each electron line is summed using matrix product.
!------------------------------------------------------------------------
module ph_mcsample 
   use pert_const, only: dp
   implicit none 
    
   real(dp) :: qa(3) = -(/ 0.5,0.5,0.5 /)
   real(dp) :: qb(3) = (/ 0.5,0.5,0.5 /)
   real(dp) :: alpha_ = 200.0_dp
end module 


module zerophonon_update_matrix_backup
   use DiagMC, only : config,fynman,  diagmc_stat, vertex, nsvd,copy_vtex,maxN,abinitio,kbT ,trace
   use pert_param, only : LatticeFrohlich,gauss_frolich
   use HamiltonianTable, only : random_qfromlist_omp
   use diagMC_gt_debug, only :check_phase
   contains 

   subroutine update_drive(diagram,stat,update)
      
      implicit none 
      type(fynman) :: diagram
      type(diagmc_stat) :: stat
      integer :: update , j_tmp 

      do j_tmp = 1,maxN
         if(.not.allocated(diagram%vertexList(j_tmp)%wq)) then 
            write(*,'(A20,i5)')'j_tmp = ',j_tmp
            stop '@update drive'
         endif 
      enddo

      select case(update)
      case(1)  
         call add_ph(diagram,stat)
      case(2)  
         call remove_ph(diagram,stat)
      case(3)  
         call change_ph_mode(diagram,stat)
      case(5)
         call update_vertex_time(diagram,stat)
      case(6)
         call update_tau(diagram,stat)
      case(7) 
         call update_swap(diagram,stat)
      end select 
   end subroutine 
    
   subroutine change_internal(diagram,stat,N)
      !change internal variable 
      use diagmc 
      implicit none 
      type(fynman), target :: diagram
      type(diagmc_stat) :: stat
      integer :: i,N, update 
      do i = 1,N
         call Discrete_omp(diagram%seed,config%PA,7,update)
         if(update.eq.6) cycle
         call update_drive(diagram,stat,update)
      end do 
   end subroutine 
   
   !1 
   subroutine add_ph(diagram, stat)
      use ph_mcsample
      use DiagMC
      implicit none 
      type(fynman), target :: diagram
      type(diagmc_stat) :: stat

      integer :: ivL, ivR, ivn1, ivn2, nu, ib,jb 
      real(dp) :: tauL, tauR, tauLR,tau1, tau2, tau12
      real(dp) :: Pq, Pnu(Nph),Ptau12, Pib, Pel_line, temp  !P(B)
      real(dp) :: P_A,P_B,D_AB, P_accept
      real(dp) :: decay, decay_exp
      real(dp) :: ran, err_add,sigma_qm, factor, Elist(nbnd), dtau 
      complex(dp) :: left_env(nbnd,nbnd), right_env(nbnd,nbnd), A(nbnd, nbnd), left_new(nbnd,nbnd), mat_old, mat_new 
      complex(dp) :: gkq(nbnd,nbnd,Nph), expH(nbnd,nbnd)
      type(vertex),pointer :: v1,v2
      type(vertex),pointer :: vL,vR, vn1, vn2 

      if(diagram%order+2.gt.maxOrder) return 
      !count update   
      stat%update(1) = stat%update(1) + 1 
      ivn1 = diagram%order + 1; ivn2 = diagram%order + 2
      vn1 => diagram%vertexList(ivn1); vn2 => diagram%vertexList(ivn2)
      

      !randomly choose an electron line 
      ivL = 1
      call uniform_int_omp(diagram%seed, 1, diagram%order, ivL )  
      vL => diagram%vertexList(ivL)
      tauL = vL%tau; 
      ivR = vL%link(3)
      vR => diagram%vertexList(ivR)
      !propose ph-momentum 

      !call sample_q_omp(diagram%seed,vn1%i_q,vn1%Pq,vn1); Pq = vn1%Pq
      call sample_q_omp_int(diagram%seed, vn1%i_q, vn1%Pq, vn1)
      
      vn1%i_kin = vL%i_kout; vn1%ekin = vL%ekout; vn1%ukin = vL%ukout
      vn1%i_kout = vn1%i_kin + vn1%i_q 
      call cal_ek_int( vn1%i_kout,  vn1%ekout,  vn1%ukout )
      call cal_wq_int( vn1%i_q,     vn1%wq    )
      call cal_gkq_vtex_int(vn1, vn1%gkq)

      do nu = 1,Nph 
         vn2%gkq(:,:,nu) = conjg(transpose(vn1%gkq(:,:,nu)))
      enddo 
      do nu = 1,Nph 
         Pnu(nu) = maxval(abs(vn1%gkq(:,:,nu)))**2 !* sqrt(1.0_dp + bose(kbT, vn1%wq(nu)))
      enddo 
      Pnu = Pnu / sum(Pnu) 
      call Discrete_omp(diagram%seed,Pnu, Nph, nu)

      if(nu>Nph) return    
      vn1%nu = nu 

      !choose tau first P(t1,t2) ~ exp(-decay_exp*(t2-t1)), t2>t1
      decay_exp = vn1%wq(nu) + minval(vn1%ekout) - minval(vn1%ekin)
      tauR = diagram%vertexList(ivR)%tau
      call uniform_real_omp(diagram%seed,tauL,tauR,tau1,temp,1)
      call exp_sample_omp(diagram%seed,tau1,tauR,decay_exp,tau2,Ptau12,1)
      Ptau12 = Ptau12 / (tauR-tauL)
      tau12 = tau2 - tau1
      if(tau12<0) then 
         write(*,'(i5,2f15.5)') diagram%order,tauL,tauR
         write(*,'(3f15.5)')tau1,tau2
         write(*,'(4E15.5)')vn1%wq(nu) , vn1%ekout(vn1%n_out) , vn1%ekin(vn1%n_in) ,decay_exp
         stop 'tau12 < 0'
      endif 


      !prob for choosing to remove this phonon line 
      P_A = 1.d0/( (diagram%order-1)/2 + 1 ) * config%PA(2)             
      !prob for choosing to add this phonon line 
      Pel_line = 1.d0 / diagram%order
      P_B = Pq * Pel_line  * Pnu(vn1%nu) * Ptau12 * config%PA(1) 
      !P_B = Pq * Pel_line * Ptau12 * config%PA(1)         
      !prob ratio W(B)/W(A) 
      D_AB = Gel(minval(vn1%ekout),tau12) / Gel(minval(vn1%ekin),tau12) !electron gf change 
      D_AB = D_AB * Dph(vn1%wq(nu),tau12)                                           !phonon gf change 


      !cal env 
      call left_environment_matrix( diagram,  ivL, left_env ); left_new = left_env 
      call right_environment_matrix( diagram,  ivR, right_env )
      !write(*,*)ivR, maxval(abs(right_env))
      !contract between tauL - tauR 
      Elist = vL%ekout - minval(vL%ekout)
      dtau = tauR-tauL
      call propagator(vL%ukout,Elist, dtau, nbnd, expH)
      left_env = matmul(expH,left_env)
      A = matmul(right_env,left_env)
      mat_old = trace(A,nbnd)

      !new : vL -> v1 
      Elist = vL%ekout - minval(vL%ekout)
      dtau = tau1-tauL
      call propagator(vL%ukout,Elist, dtau, nbnd, expH)
      left_new = matmul(expH,left_new)
      !v1 e-ph 
      left_new = matmul(vn1%gkq(:,:,vn1%nu), left_new)
      !v1 -> v2 
      Elist = vn1%ekout - minval(vn1%ekout)
      dtau = tau2-tau1
      call propagator(vn1%ukout, Elist, dtau, nbnd, expH)
      left_new = matmul(expH,left_new) 
      !v2 e-ph 
      left_new = matmul(vn2%gkq(:,:,vn1%nu), left_new)
      !v2 -> vR
      Elist = vR%ekin - minval(vR%ekin)
      dtau = tauR-tau2
      call propagator(vR%ukin, Elist, dtau, nbnd, expH)
      left_new = matmul(expH,left_new) 
      mat_new = trace(matmul(right_env, left_new),nbnd) 

      factor = real(mat_new) / real(mat_old)
      D_AB = D_AB * abs(factor) 


      P_accept = D_AB * P_A / P_B 
      call random_number_omp(diagram%seed,ran)
      if(ran<P_accept) then 
         v1 => diagram%vertexList(ivn1) 
         v2 => diagram%vertexList(ivn2) 
         v1%tau   = tau1; v2%tau = tau2
         v2%i_kin   = v1%i_kout   
         v2%i_kout  = v1%i_kin 
         v2%i_q     = -v1%i_q
         v2%Pq    = v1%Pq !for ab-initio
         v2%nu    = v1%nu 
         v2%wq    = v1%wq; 
         v2%ekin  = v1%ekout; v2%Einmin = minval(v2%ekin)
         
         v2%ekout = v1%ekin; v2%Eoutmin = minval(v2%ekout)  
         v2%ukin  = v1%ukout  
         v2%ukout = v1%ukin 
         do nu = 1,Nph 
            v2%gkq(:,:,nu) = conjg(transpose(v1%gkq(:,:,nu)))
         end do 

         !copy of formfactor, maybe not needed 
         !v2%uk_f = v1%uk_f
         !v2%ukq_f = v1%ukq_f
         !v2%vq_f = v1%vq_f
         !v2%vnq_f = v1%vnq_f
        
         
         !link info 
         v1%link = (/ ivL, ivn2, ivn2/)
         v2%link = (/ ivn1, ivn1, ivR /)
         diagram%vertexList(ivL)%link(3) = ivn1
         diagram%vertexList(ivR)%link(1) = ivn2
         stat%accept(1) = stat%accept(1) + 1
         diagram%order = diagram%order + 2
      end if 

      return 
   end subroutine 
    
   !2
   subroutine remove_ph(diagram, stat)
      use ph_mcsample
      use DiagMC 
      implicit none 
      type(fynman), target :: diagram
      type(diagmc_stat) :: stat
      integer :: iv1, iv2,ivL,ivR, nu,nup, ib
      type(vertex), pointer :: vL,vR, vLL, vRR, v1,v2

      real(dp) :: tauL, tauR,tauLR, tau1, tau2, tau12, temp, Ptau12
      real(dp) :: xqt(3)
      real(dp) :: Pq,Pnu(Nph), Pib, Pel_line          !P(B)
      real(dp) :: Pph_line                            !Prob for selecting one ph line 
      real(dp) :: ek, ekq, Elist(nbnd), dtau
      real(dp) :: P_accept, P_A, P_B,D_AB 
      real(dp) :: wq, decay, decay_exp,errG
      real(dp) :: ran, sigma_q, factor 
      complex(dp) :: left_env(nbnd,nbnd), right_env(nbnd,nbnd), A(nbnd, nbnd), left_new(nbnd,nbnd), mat_old, mat_new 
      complex(dp) :: expH(nbnd,nbnd)
      if(diagram%order.eq.1)return 

      !randomly choose a ph-line 
      call uniform_int_omp(diagram%seed,2,diagram%order,iv1)
      v1 => diagram%vertexList(iv1);  iv2 = v1%link(2)
      v2 => diagram%vertexList(iv2)
      if(iv2.eq.1 .or. iv2.eq.maxN) return 
      if(v2%tau<v1%tau) then  !make sure vR vertex is later in time 
         iv1 = iv2
         v1 => diagram%vertexList(iv1); iv2 = v1%link(2)
         v2 => diagram%vertexList(iv2)
      end if 
      if( v1%link(3).ne.v1%link(2) ) return      !vL vR should be connected by a single electron line 

      
      vL => diagram%vertexList(v1%link(1)) 
      vR => diagram%vertexList(v2%link(3)) 
      !count update 
      stat%update(2) = stat%update(2) + 1


      tauL = vL%tau 
      tauR = vR%tau 
      tau1 = v1%tau; 
      tau2 = v2%tau 
      nu = v1%nu; wq = v1%wq(nu)

      !tau1 tau2 distribution
      tau12 = tau2 - tau1  
      if(tau12<0) then 
         write(*,*)tau12
         stop 'time disorder'
         call mp_global_end() 
      end if 
      decay_exp = wq + minval(v1%ekout) - minval(v1%ekin)
      call exp_sample_omp(diagram%seed,tau1,tauR,decay_exp,tau2,Ptau12)
      Ptau12 = Ptau12 / (tauR-tauL)

      !nu probility
      do nu = 1,Nph 
         Pnu(nu) = maxval(abs(v1%gkq(:,:,nu)))**2 !* sqrt(1.0_dp + bose(kbT, v1%wq(nu)))
      enddo 
      Pnu = Pnu / sum(Pnu) 

      ![2.1] prob for adding ph line  
      if(abinitio) then 
         Pq = v1%Pq
      else 
         Pq = 1.d0 
      endif 
      Pel_line = 1.d0 / (diagram%order - 2); 
      P_A = Pq * Pel_line  * Pnu(v1%nu) * Ptau12 * config%PA(1)
      ![2.2] prob for removing ph line 
      Pph_line = 1.d0 / ( (diagram%order-1)/2 )
      P_B = Pph_line * config%PA(2)      !there are (order-1)/2 phonon line 
      D_AB = 1.d0/( Gel(minval(v1%ekout),tau12)/Gel(minval(v1%ekin),tau12) * Dph(wq, tau12 )  )

      !matrix term 
      call left_environment_matrix( diagram, v1%link(1), left_env ); left_new = left_env 
      call right_environment_matrix(diagram, v2%link(3), right_env )
      Elist = vL%ekout - minval(vL%ekout) 
      dtau = vR%tau - vL%tau 
      call propagator(vL%ukout, Elist, dtau, nbnd, expH)
      left_env = matmul(expH,left_env)
      mat_old = trace(matmul(right_env, left_env),nbnd)
      !new vL -> v1
      Elist = vL%ekout - minval(vL%ekout)
      dtau = tau1-tauL
      call propagator(vL%ukout,Elist, dtau,nbnd, expH)
      left_new = matmul(expH,left_new)
      !v1 e-ph 
      left_new = matmul(v1%gkq(:,:,v1%nu), left_new)
      !v1 -> v2 
      Elist = v1%ekout - minval(v1%ekout)
      dtau = tau2-tau1
      call propagator(v1%ukout, Elist, dtau, nbnd, expH)
      left_new = matmul(expH,left_new) 
      ! v2 e-ph 
      left_new = matmul(v2%gkq(:,:,v2%nu),left_new)

      Elist = vR%ekin - minval(vR%ekin)
      dtau = tauR-tau2
      call propagator(vR%ukin, Elist, dtau, nbnd, expH)
      left_new = matmul(expH,left_new) 
      mat_new = trace(matmul(right_env, left_new), nbnd)
     
      factor = real(mat_old) / real(mat_new)  
      D_AB = D_AB * abs(factor) 

      P_accept = D_AB * P_A / P_B 
      !write(*,'(4E15.5)') P_accept, D_AB,  P_A / P_B
      call random_number_omp(diagram%seed, ran )
      if(ran<P_accept) then 
   
         !removal also create a new link, make sure the new electron line is consistent 
         errG =  maxval(abs(vR%ukin - vL%ukout)) / maxval(abs(vR%ukin))
         if(errG>1e-10) then 
            if(ivL.ne.1) then 
               vL%i_kout = vR%i_kin 
               vL%ekout = vR%ekin 
               vL%vkout = vR%vkin
               vL%ukout = vR%ukin
               !call cal_gkq_vtex(vL,vL%gkq)
            else 
               if(ivR.ne.maxN) then 
                  vR%i_kin = vL%i_kout  
                  vR%ekin  = vL%ekout 
                  vR%vkin = vL%vkout
                  vR%ukin = vL%ukout 
                  !call cal_gkq_vtex(vR,vR%gkq)
               endif 
            endif 
         endif 

         call delete_vertex(diagram,iv1,iv2)
         stat%accept(2) = stat%accept(2) + 1
      end if 
   end subroutine 
    
   !3 change the phonon modes, this may change the phase  
   subroutine change_ph_mode(diagram, stat)
      use DiagMC 
      implicit none 
      type(fynman), target :: diagram
      type(diagmc_stat) :: stat
      integer :: iv1,iv2
      integer :: nu,m1,n1,m2,n2, nup, ib
      real(dp) :: Pnu(Nph),absPnu(Nph), tau, factor, dtau, Elist(nbnd)
      type(vertex), pointer :: v1,v2,vL,vR
      complex(dp) :: Dratio
      complex(dp) :: left_env(nbnd,nbnd), right_env(nbnd,nbnd),LR(nbnd,nbnd), mid_env(nbnd, nbnd), A(nbnd, nbnd), left_new(nbnd,nbnd)
      complex(dp) :: expH(nbnd,nbnd), mat_old, mat_new 
      
      if(diagram%order.eq.1) return 
      if(Nph.eq.1) return 
      
      !randomly choose ph line 
      call uniform_int_omp(diagram%seed,2,diagram%order,iv1)
      v1 => diagram%vertexList(iv1); iv2 = v1%link(2)
      if(iv2.eq.1 .or. iv2.eq.maxN) then 
         !call change_extph_mode(iv1,diagram, stat)
         return 
      endif 
      v2 => diagram%vertexList(iv2)
      if(v1%tau>v2%tau) then 
         iv1 = iv2 
         v1 => diagram%vertexList(iv1); iv2 = v1%link(2)
         v2 => diagram%vertexList(iv2)
      endif 
      tau =  v2%tau - v1%tau 
      if(tau<0) stop 'tau disorder in change_ph_mode '
      vL => diagram%vertexList(v1%link(1))
      vR => diagram%vertexList(v2%link(3))

      call left_environment_matrix(diagram, v1%link(1), left_env)
      !add time propagate : v1$link(1) -> v1 
      Elist = v1%ekin - minval(v1%ekin)
      dtau = v1%tau - vL%tau
      call propagator(v1%ukin, Elist, dtau, nbnd, expH)
      left_env = matmul(expH,left_env)
      call right_environment_matrix(diagram, v2%link(3), right_env)
      !add time propagate : v2 -> v2$link(3)  
      Elist = v2%ekout - minval(v2%ekout)
      dtau = vR%tau - v2%tau
      call propagator(v2%ukout, Elist, dtau, nbnd, expH)
      right_env = matmul(right_env,expH)

      call mid_environment_matrix(diagram, iv1, iv2, mid_env)
      !weight = trace(LR * g2 *M * g1)
      do nu = 1,Nph 
         A = matmul(v1%gkq(:,:,nu), left_env)
         A = matmul(mid_env, A)
         A = matmul(v2%gkq(:,:,nu),A)
         A = matmul(right_env,A)
         Pnu(nu) =  real(trace(A, nbnd)) * Dph(v1%wq(nu),tau)
      end do 
      absPnu = abs(Pnu) 
      absPnu = absPnu / sum(absPnu)
      call Discrete_omp(diagram%seed,absPnu,Nph,nu)           !sampling wrt |D|, acceptance is 1 
      if(nu>Nph) return
      if(nu.eq.v1%nu) return 
      v1%nu = nu; v2%nu = nu
      stat%update(3) = stat%update(3) + 1
      stat%accept(3) = stat%accept(3) + 1

   end subroutine 

    
   !5
   subroutine update_vertex_time(diagram, stat)
      use DiagMC 
      implicit none 
      type(fynman), target :: diagram
      type(diagmc_stat) :: stat
      integer :: iv1, ib 
      type(vertex), pointer :: vA
      real(dp) :: tauL, tauR, tauQ,tau1, decay, Pt, sign_wq, Elist(nbnd), dtau
      real(dp) :: wq,ran, Paccept, factor 
      complex(dp) :: expH(nbnd,nbnd), left_env(nbnd,nbnd), right_env(nbnd,nbnd), A(nbnd, nbnd), left_new(nbnd,nbnd)
      complex(dp) :: mat_old, mat_new 

      if(diagram%order.eq.1) return 

      stat%update(5) = stat%update(5) + 1
      call uniform_int_omp(diagram%seed,2,diagram%order,iv1) !choose a vertex 
      vA => diagram%vertexList(iv1)
      tauL = diagram%vertexList(vA%link(1))%tau
      tauR = diagram%vertexList(vA%link(3))%tau
      tauQ = diagram%vertexList(vA%link(2))%tau

      wq = vA%wq(vA%nu)
      sign_wq = -1 
      if(vA%tau>tauQ) sign_wq = 1 
      decay = minval(vA%ekin) - minval(vA%ekout) + sign_wq * wq

      if(tauL .eq. tauR) stop "tauLR equal @ update_vertex_time"
      if(decay .eq. 0.d0) stop "decay .eq. 0 @ update_vertex_time"
      call exp_sample_omp(diagram%seed,tauL,tauR,decay,tau1,Pt,1)

      Paccept = Dph( wq, abs(tau1-tauQ) ) / Dph( wq, abs(vA%tau-tauQ) ) !ph gf change 
      Paccept = Paccept * exp( sign_wq*wq*(tau1-vA%tau) ) 
      !matrix term 
      call left_environment_matrix(diagram,  vA%link(1),left_env)
      call right_environment_matrix(diagram,  vA%link(3),right_env)
      left_new = left_env 
   
      !old matrix 
      Elist = vA%ekin-minval(vA%ekin)
      dtau =  vA%tau - tauL
      call propagator(vA%ukin, Elist, dtau, nbnd, expH)
      left_env = matmul(expH, left_env)
      !vA e-ph 
      left_env = matmul(vA%gkq(:,:,vA%nu), left_env)

      Elist = vA%ekout-minval(vA%ekout)
      dtau =  tauR-vA%tau
      call propagator(vA%ukout, Elist, dtau, nbnd, expH)
      left_env = matmul(expH, left_env)
      mat_old = trace( matmul(right_env, left_env), nbnd)

      !new matrix 
      Elist = vA%ekin-minval(vA%ekin)
      dtau =  tau1 - tauL
      call propagator(vA%ukin, Elist, dtau, nbnd, expH)
      left_new = matmul(expH, left_new)
      !vA e-ph 
      left_new = matmul(vA%gkq(:,:,vA%nu), left_new)

      Elist = vA%ekout-minval(vA%ekout)
      dtau =  tauR-tau1
      call propagator(vA%ukout, Elist, dtau, nbnd, expH)
      left_new = matmul(expH, left_new)
      mat_new = trace( matmul(right_env, left_new), nbnd)

      factor = real(mat_new) / real(mat_old)
      Paccept = Paccept * abs(factor) 

      !write(*,*) Paccept
      call random_number_omp(diagram%seed,ran)
      if(ran < Paccept) then 
         vA%tau = tau1
         stat%accept(5) = stat%accept(5) + 1
      end if 
   end subroutine 
   !6 
   subroutine update_tau(diagram, stat)
      use DiagMC 
      implicit none 
      type(fynman), target :: diagram
      type(diagmc_stat) :: stat
      real(dp) :: tauL, decay, newT, Pt 
      type(vertex), pointer :: vR 
      real(dp) :: Paccept, ran, factor, Elist(nbnd), dtau 
      integer :: ib 
      complex(dp) :: left_env(nbnd,nbnd), expH(nbnd,nbnd), mat_old, mat_new 

      vR =>  diagram%vertexList(maxN)
      decay = minval(vR%ekin) 

      tauL = diagram%vertexList(vR%link(1))%tau 
      if(tauL<config%tauMin) tauL = config%tauMin
      call exp_sample_omp(diagram%seed,tauL,config%tauMax, decay, newT, Pt,1)
      stat%update(6) = stat%update(6) + 1

      !matrix term 
      call left_environment_matrix(diagram, vR%link(1), left_env )
      Elist = vR%ekin - decay
      dtau = vR%tau - tauL
      call propagator(vR%ukin, Elist, dtau, nbnd, expH )
      mat_old = trace( matmul(expH,left_env), nbnd )
      dtau = newT - tauL
      call propagator(vR%ukin, Elist, dtau, nbnd, expH )
      mat_new = trace( matmul(expH,left_env), nbnd )
   
      factor = real(mat_new) / real(mat_old) 
      Paccept = abs(factor) 
      call random_number_omp(diagram%seed, ran) 
      if(ran<Paccept) then 
         stat%accept(6) = stat%accept(6) + 1      
         vR%tau = newT
      endif 
   end subroutine 
    
   !7 
   subroutine update_swap(diagram, stat)
      use DiagMC 

      implicit none 
      type(fynman), target :: diagram
      type(diagmc_stat) :: stat
      integer :: iv1, iv2, ib
      type(vertex), pointer :: vL,vR, vLL, vRR

      real(dp) :: tau1, tau2, tau12,tauL,tauR,tauLL,tauRR, decay, wq1, wq2,Elist(nbnd), dtau
      real(dp) :: P_accept, P_kchange
   
      real(dp) :: ran
      real(dp) ::  factor 
      type(vertex),pointer :: vLn,vRn
      complex(dp) :: left_env(nbnd,nbnd), right_env(nbnd,nbnd), A(nbnd, nbnd), left_new(nbnd,nbnd), mat_old, mat_new 
      complex(dp) :: expH(nbnd,nbnd)
      if(diagram%order<=3)return 

      !randomly choose a pair of vertex 
      call uniform_int_omp(diagram%seed,2,diagram%order,iv1)
      vL => diagram%vertexList(iv1);  iv2 = vL%link(3)
      if(vL%link(2).eq.iv2) return            !the two vertex are not phonon linked 
      if(iv2.eq.maxN)       return            !not the last vertex   
      !count update 
      stat%update(7) = stat%update(7) + 1

      vR => diagram%vertexList(iv2)
      vLL => diagram%vertexList(vL%link(1))
      vRR => diagram%vertexList(vR%link(3))
      tauL = vL%tau; tauLL = vLL%tau; 
      tauR = vR%tau; tauRR = vRR%tau

      if(tauL>tauR) stop ' time disordered @ update_swap'

      vLn => diagram%vertexList(diagram%order+1) 
      vRn => diagram%vertexList(diagram%order+2) 

      vRn%tau = vL%tau
      !q depdent quantity 
      vRn%nu = vR%nu; vRn%i_q = vR%i_q; vRn%Pq = vR%Pq
      vRn%wq = vR%wq; 
      !vRn%uk_f = vR%uk_f;  vRn%ukq_f = vR%ukq_f
      !vRn%vq_f = vR%vq_f;  vRn%vnq_f = vR%vnq_f
      !kin depdent quantity 
      vRn%i_kin = vL%i_kin; vRn%ekin = vL%ekin; vRn%ukin = vL%ukin;  
      !kout depdent quantity 
      vRn%i_kout = vRn%i_kin + vRn%i_q
      call cal_ek_int(vRn%i_kout, vRn%ekout,vRn%ukout)
      !new e-ph vertex@vRn 
      call cal_gkq_vtex( vRn, vRn%gkq ) 

      vLn%tau = vR%tau
      !q depdent quantity 
      vLn%nu = vL%nu; vLn%i_q = vL%i_q; vLn%Pq = vL%Pq
      vLn%wq = vL%wq; 
      !vLn%uk_f = vL%uk_f; vLn%ukq_f = vL%ukq_f
      !vLn%vq_f = vL%vq_f; vLn%vnq_f = vL%vnq_f
      !kin depdent quantity 
      vLn%i_kin = vRn%i_kout; vLn%ekin = vRn%ekout; vLn%ukin = vRn%ukout 
      !kout depdent quantity 
      vLn%i_kout = vR%i_kout; vLn%ekout = vR%ekout; vLn%ukout = vR%ukout 
      !new e-ph vertex@vLn
      call cal_gkq_vtex( vLn, vLn%gkq ) 

      !electron propagator ratio 
      tau12 =  vR%tau - vL%tau
      P_kchange = Gel( minval(vRn%ekout),tau12 ) / Gel( minval(vL%ekout),tau12 )

      !ph propagator ratio 
      wq1 = vL%wq(vL%nu); wq2 = vR%wq(vR%nu)
      P_kchange = P_kchange * Dph( wq1, abs( vR%tau - diagram%vertexList(vL%link(2))%tau ) ) /&
                  Dph( wq1, abs( vL%tau - diagram%vertexList(vL%link(2))%tau ) )
      P_kchange = P_kchange * Dph( wq2, abs( vL%tau - diagram%vertexList(vR%link(2))%tau ) ) /&
                  Dph( wq2, abs( vR%tau - diagram%vertexList(vR%link(2))%tau ) )
       
      !matrix term for new e-ph and propagator 
      call left_environment_matrix(diagram, vL%link(1), left_env);
      Elist = vL%ekin - minval(vL%ekin)
      dtau =  vL%tau - vLL%tau
      call propagator(vL%ukin, Elist, dtau, nbnd, expH)
      left_env = matmul(expH, left_env)
      left_new = left_env

      call right_environment_matrix(diagram,  vR%link(3), right_env)
      Elist = vR%ekout - minval(vR%ekout)
      dtau =  vRR%tau - vR%tau
      call propagator(vR%ukout, Elist, dtau, nbnd, expH)
      right_env = matmul(right_env, expH) 

      !old mat 
      !e-ph  @vL
      left_env = matmul(vL%gkq(:,:,vL%nu), left_env)
      !propagator : vL -> vR
      Elist = vL%ekout - minval(vL%ekout)
      dtau =  vR%tau - vL%tau
      call propagator(vL%ukout, Elist, dtau, nbnd, expH)
      left_env = matmul(expH, left_env)
      !e-ph @ vR
      left_env = matmul(vR%gkq(:,:,vR%nu), left_env)

      A = matmul(right_env,left_env)
      mat_old = trace(A, nbnd)

      !new mat 
      !e-ph @vRn
      left_new = matmul(vRn%gkq(:,:,vRn%nu), left_new)
      !propagator : vRn -> vLn 
      Elist = vRn%ekout - minval(vRn%ekout)
      call propagator(vRn%ukout, Elist, dtau, nbnd, expH)
      left_new = matmul(expH, left_new)
      !e-ph @vLn 
      left_new = matmul(vLn%gkq(:,:,vLn%nu), left_new)

      A = matmul(right_env, left_new)
      mat_new = trace(A, nbnd)

      
      factor = real(mat_new) / real(mat_old)
      !write(*,'(3E20.10)')  maxval(abs(left_env)),maxval(abs(right_env)),P_kchange
      P_accept = abs(factor) * P_kchange  
      
      call random_number_omp(diagram%seed,ran)
      
      if(ran<P_accept) then 
         vLn%link = vL%link;       vRn%link = vR%link
         vLn%plink = vL%plink;     vRn%plink = vR%plink 
         vLn%link(1) = iv2;        vLn%link(3) = vR%link(3) 
         vRn%link(1) = vL%link(1); vRn%link(3) = iv1  
         diagram%vertexList(vL%link(1))%link(3) = iv2
         diagram%vertexList(vR%link(3))%link(1) = iv1
         
         call copy_vtex(Nph, dmc_band, nbnd,nsvd,vLn,diagram%vertexList(iv1))
         call copy_vtex(Nph, dmc_band, nbnd,nsvd,vRn,diagram%vertexList(iv2))
         stat%accept(7) = stat%accept(7) + 1
      end if 

   end subroutine 
    
    
   !----------------------------------------------------------------------------
   !-----------  driver of vertex storage management ---------------------------
   !----------------------------------------------------------------------------
   subroutine insert_vertex(diagram,v1, v2) !push_back 
      implicit none 
      type(fynman), target :: diagram
      type(vertex) :: v1,v2 
      integer :: v1_in, v2_out
      integer :: iv1,iv2  
      if(v1%tau>v2%tau) stop 'time disorder @ insert_vertex : v1%tau should < v2%tau'
      v1_in = v1%link(1); v2_out = v2%link(3)
      iv1 = diagram%order + 1; iv2 = diagram%order + 2
      v1%link(2:3) = iv2; v2%link(1:2) = iv1 
      diagram%vertexList(v1_in)%link(3) = iv1
      diagram%vertexList(v2_out)%link(1) = iv2 
      diagram%vertexList(iv1) = v1
      diagram%vertexList(iv2) = v2
      diagram%order = diagram%order + 2
   end subroutine 

   subroutine swap_vertex(diagram,iv1,iv2)
      implicit none 
      type(fynman) :: diagram
      integer :: iv1,iv2 
      if(iv1.eq.iv2) return 
      call move_vertex(diagram,iv1,diagram%order+1) ! 1 -> order+1
      call move_vertex(diagram,iv2,iv1)             ! 2 -> 1
      call move_vertex(diagram,diagram%order+1,iv2) !order + 1 -> 2
   end subroutine 

   subroutine move_vertex(diagram,iv1,iv2)
      !move vertex iv1 to iv2 -th place 
      use DiagMC, only : dmc_band,nbnd,Nph
      implicit none 
      type(fynman) :: diagram
      integer :: iv1,iv2,j_tmp
      integer :: link(3)

      call copy_vtex(Nph, dmc_band, nbnd,nsvd,diagram%vertexList(iv1),diagram%vertexList(iv2))

      link = diagram%vertexList(iv1)%link
      diagram%vertexList(link(1))%link(3) = iv2
      diagram%vertexList(link(2))%link(2) = iv2
      diagram%vertexList(link(3))%link(1) = iv2

   end subroutine 

   subroutine delete_vertex(diagram, iv1,iv2)
      implicit none 
      integer :: iv1,iv2 
      type(fynman), target :: diagram
      type(vertex),pointer :: v1,v2 
     
      integer :: v1_in, v2_out, j_tmp 
   
      v1 => diagram%vertexList(iv1)
      v2 => diagram%vertexList(iv2)
      if( v1%link(3).ne.v1%link(2) ) stop 'link err@delete_vertex 1'

      if(iv1.eq.diagram%order .and. iv2.eq.diagram%order-1) then 
         call swap_vertex(diagram,iv2,iv1)
      else 
         call swap_vertex(diagram,iv1,diagram%order-1)
         call swap_vertex(diagram,iv2,diagram%order)
      endif 


      v1 => diagram%vertexList(diagram%order-1)
      v2 => diagram%vertexList(diagram%order)
      if( v1%link(3).ne.v1%link(2) ) then 
         write(*,*) iv1,iv2,diagram%order 
         write(*,'(3i5)') v1%link
         write(*,'(3i5)') v2%link

         write(*,'(3i5)')  diagram%vertexList(iv1)%link
         write(*,'(3i5)')  diagram%vertexList(iv2)%link
         stop 'link err@delete_vertex 2'
      endif 
   
      v1_in = v1%link(1); v2_out = v2%link(3)
      diagram%vertexList(v1_in)%link(3) = v2_out
      diagram%vertexList(v2_out)%link(1) = v1_in

      diagram%order = diagram%order - 2           !then the last two elements are forgetten 
      diagram%vertexList(diagram%order+1)%tau = 0.0
      diagram%vertexList(diagram%order+2)%tau = 0.0
        
   end subroutine 

end module 






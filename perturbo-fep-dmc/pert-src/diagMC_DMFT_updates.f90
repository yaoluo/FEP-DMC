module selfE_update_DMFT  
   use DiagMC, only : FakeRatio,dp, config,fynman, diagmc_stat, vertex, copy_vtex,maxN,&
                     Dph,Gel,nbnd,Nph,sample_q_omp,cal_wq,cal_ek,cal_gkq_symmetrize
   use random_tool, only : uniform_int_omp,Discrete_omp,random_number_omp
   

   contains 

   function Equal(k1,k2)
      use mod_HTable, only : IndexK,Nkxyz
      use DiagMC, only : Holstein,dim_Holstein
      implicit none 
      Logical :: Equal 
      real(dp) :: k1(3),k2(3), k12(3)
      integer :: ik1,ik2 
      
      Equal = .false.
      if(Holstein) then 
         k12 = k1-k2 
         call pbc(k12)
         if(sum(abs(k12(1:dim_Holstein)))<1e-8) Equal = .true. 
      else 
         ik1 = IndexK(k1,Nkxyz)
         ik2 = IndexK(k2,Nkxyz)
         if(ik1.eq.ik2) Equal = .true. 
      endif 
      contains 
      subroutine pbc(kpt)
         implicit none 
         real(dp) :: kpt(3)
         integer :: i 
         do i = 1,3 
            if(kpt(i)<0) kpt(i) = kpt(i) + int(abs(kpt(i))+1)
            if(kpt(i)>1) kpt(i) = kpt(i) - int(kpt(i))
         enddo 
      end subroutine
   end function
   subroutine update_drive(diagram,stat,update)
      
      implicit none 
      type(fynman) :: diagram
      type(diagmc_stat) :: stat
      integer :: update , j_tmp 

      select case(update)
      case(1)  
         call add_ph(diagram,stat)
      case(2)  
         call remove_ph(diagram,stat)
      case(3)  
         call change_ph_mode(diagram,stat)
      case(4)  
         call change_band(diagram,stat)
      case(5)
         call update_vertex_time(diagram,stat)
      case(6)
         call update_external(diagram,stat)
      case(7) 
         call update_swap(diagram,stat)
      case(8)
         call ForNormalization(diagram,stat)
      end select 
      if(diagram%order.eq.0) call Change_K(diagram, stat)
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

      integer :: ivL, ivR, ivn1, ivn2, nu 
      real(dp) :: tauL, tauR, tauLR,tau1, tau2, tau12
      real(dp) :: Pq, Pnu(Nph),Ptau12, Pib, Pel_line, temp  !P(B)
      real(dp) :: P_A,P_B,D_AB, P_accept
      real(dp) :: decay, decay_exp
      real(dp) :: ran, err_add
      complex(dp) :: gkq(nbnd,nbnd,Nph)
      type(vertex),pointer :: v1,v2
      type(vertex),pointer :: vL, vn1, vn2 

      if(diagram%order+2.gt.maxOrder) return 
      if(diagram%order .eq. 0) return 
      !count update   
      stat%update(1) = stat%update(1) + 1 

      ivn1 = diagram%order + 1; ivn2 = diagram%order + 2
      vn1 => diagram%vertexList(ivn1); vn2 => diagram%vertexList(ivn2)
      ![1.1] randomly choose an electron line 
      ivL = 1
      call uniform_int_omp(diagram%seed, 1, diagram%order, ivL )  !randomly choose a vertex 
      if( ivL.eq.diagram%tail ) return !give up if select the last one 
      vL => diagram%vertexList(ivL)
      tauL = vL%tau; 
      ivR = vL%link(3)
      
      ![1.2]propose ph-momentum 
      call sample_q_omp(diagram%seed,vn1%q); Pq = 1.d0 !choose from pre-tabulated mesh 
      
      vn1%kin = vL%kout; 
      vn1%kout = vn1%kin + vn1%q
      call cal_ek( vn1%kin,   vn1%ekin ,  vn1%ukin  )
      !vn1%ekin = vL%ekout; vn1%ukin = vL%ukout
      call cal_ek( vn1%kout,  vn1%ekout,  vn1%ukout )
      call cal_wq( vn1%q,     vn1%wq,     vn1%uq    )
      call cal_gkq_symmetrize(vn1%kin, vn1%q, vn1%wq, vn1%uq, vn1%ukin, vn1%ukout, vn1%gkq)

      ![1.3] propose el/ph-index ~ abs(gkq)
      vn1%n_in = vL%n_out
      call uniform_int_omp(diagram%seed,1,nbnd,vn1%n_out) !randomly choose a band index 
      do nu = 1,Nph 
         Pnu(nu) = abs(vn1%gkq(vn1%n_out,vn1%n_in,nu))**2 
      enddo 
      Pnu = Pnu / sum(Pnu) 
      call Discrete_omp(diagram%seed,Pnu, Nph, nu);
      if(nu>Nph) return    
      vn1%nu = nu 
         

      ![1.4] choose tau first P(t1,t2) ~ exp(-decay_exp*(t2-t1)), t2>t1
      decay_exp = vn1%wq(nu) + vn1%ekout(vn1%n_out) - vn1%ekin(vn1%n_in) 
      tauR = diagram%vertexList(ivR)%tau
      tauLR = tauR - tauL; 
      !call exp_sample_omp( diagram%seed, 0.d0, tauLR, decay_exp, tau12, Ptau12, 1 )
      !if(tau12 > tauR-tauL) stop 'tau sample is out of range @ add_ph()' 
      !call uniform_real_omp(diagram%seed,tauL,tauR-tau12,tau1,temp,1); 
      !tau2 = tau1 + tau12 
      !Ptau12 = Ptau12*temp
      call uniform_real_omp(diagram%seed,tauL,tauR,tau1,temp,1)
      call exp_sample_omp(diagram%seed,tau1,tauR,decay_exp,tau2,Ptau12,1)
      Ptau12 = Ptau12 * temp 
      tau12 = tau2 - tau1

      ![2.1]prob for choosing to remove this phono line 
      P_A = 1.d0/( diagram%order/2 + 1 ) * config%PA(2)                       
      ![2.2]prob for choosing to add this phono line 
      Pib = 1.d0 / nbnd
      Pel_line = 1.d0 / diagram%order !electron line selected as the outcome leg of a vertex
      P_B = Pq * Pib * Pnu(nu) * Ptau12  * Pel_line * config%PA(1)         
      ![2.3]prob ratio W(B)/W(A) 
      D_AB = abs(vn1%gkq(vn1%n_out,vn1%n_in,nu))**2
      D_AB = D_AB * Gel(vn1%ekout(vn1%n_out),tau12)/Gel(vn1%ekin(vn1%n_in),tau12) !electron gf change 

      D_AB = D_AB * Dph(vn1%wq(nu),tau12)                                           !phonon gf change 
      P_accept = D_AB * P_A / P_B 


      call random_number_omp(diagram%seed,ran)
      if(ran<P_accept) then 

         v1 => diagram%vertexList(ivn1) 
         v2 => diagram%vertexList(ivn2) 
         v1%tau   = tau1; v2%tau = tau2
         v2%kin   = v1%kout   
         v2%kout  = v1%kin 
         v2%q     = -v1%q
         v2%n_in  = v1%n_out 
         v2%n_out = v1%n_in 
         v2%nu    = v1%nu 
         v2%wq    = v1%wq; 
         v2%ekin  = v1%ekout
         
         v2%ekout = v1%ekin  
         v2%uq    = conjg(v1%uq)
         v2%ukin  = v1%ukout  
         v2%ukout = v1%ukin 
         do nu = 1,Nph 
            v2%gkq(:,:,nu) = conjg(transpose(v1%gkq(:,:,nu)))
         end do 
         
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
      integer :: iv1, iv2,ivL,ivR, nu,nup
      type(vertex), pointer :: vL,vR

      real(dp) :: tauL, tauR,tauLR, tau1, tau2, tau12, temp, Ptau12
      real(dp) :: xqt(3)
      real(dp) :: Pq,Pnu(Nph), Pib, Pel_line          !P(B)
      real(dp) :: Pph_line                            !Prob for selecting one ph line 
      real(dp) :: ek, ekq
      real(dp) :: P_accept, P_A, P_B,D_AB 
      real(dp) :: wq, decay, decay_exp
      real(dp) :: ran

      if(diagram%order<=2) return 
      
      ![1.1]randomly choose a ph-line 
      call uniform_int_omp(diagram%seed,1,diagram%order,iv1)
      !if the one of the two vertex are head or tail, drop the update.  
      if(iv1.eq.diagram%head) return 
      if(iv1.eq.diagram%tail) return 
      vL => diagram%vertexList(iv1);  iv2 = vL%link(2)
      if(iv2.eq.diagram%head) return 
      if(iv2.eq.diagram%tail) return 

      vR => diagram%vertexList(iv2)
      if(vR%tau<vL%tau) then  !make sure vR vertex is later in time 
          vL => diagram%vertexList(iv2)
          vR => diagram%vertexList(iv1)
      end if 
      ivL = vR%link(2)
      ivR = vL%link(2)

      if( vL%link(3).ne.ivR ) return      !vL vR should be connected by a single electron line 
      if( vL%n_in .ne. vR%n_out ) return  !for multi-band : the band index of the incoming & outgoing band index is the same 
      
      !count update 
      stat%update(2) = stat%update(2) + 1
      
      !local configuration info 
      xqt = vL%q 
      tauL = diagram%vertexList(vL%link(1))%tau 
      tauR = diagram%vertexList(vR%link(3))%tau 
      tau1 = vL%tau; tau2 = vR%tau 
      ek = vL%ekin(vL%n_in); ekq = vL%ekout(vL%n_out); 
      nu = vL%nu; wq = vL%wq(nu)

      ![2.1] Priori prob of tau1 tau2 : Ptau12
      tau12 = tau2 - tau1  
      if(tau12<0) then 
         write(*,*)tau12
         stop 'time disorder'
         call mp_global_end() 
      end if 
      decay_exp = wq + ekq - ek
      tauLR = tauR - tauL
      !call exp_sample(0.d0,tauLR,decay_exp,tau12,Ptau12)
      !Ptau12 = Ptau12 / (tauLR - tau12)
      call uniform_real_omp(diagram%seed,tauL,tauR,tau1,temp)
      call exp_sample_omp(diagram%seed,tau1,tauR,decay_exp,tau2,Ptau12)
      Ptau12 = Ptau12 * temp 

      ![2.2] Priori prob of nu : Pnu
      !decay =  ekq - ek
      do nup = 1,Nph 
         Pnu(nup) = abs(vL%gkq(vL%n_out,vL%n_in,nup))**2
      end do 
      Pnu = Pnu / sum(Pnu)

      ![2.3] prob of adding this ph line : P_A
      Pq = 1.d0; 
      Pib = 1.d0/nbnd; 
      Pel_line = 1.d0 / (diagram%order-2) !there are order-2 vertex of the lower order diagram 
      P_A = Pq * Pib * Pnu(nu) * Ptau12  * Pel_line * config%PA(1)
      P_A =  Ptau12 * Pel_line * config%PA(1) 

      ![2.4] prob of removing this ph line : P_B 
      Pph_line = 2.d0/ diagram%order      !there are order/2 phonon line, randomly choose one  
      P_B = Pph_line * config%PA(2)      
      ![2.5] ratio of weight 
      D_AB = 1.d0/( abs(vL%gkq(vL%n_out,vL%n_in,nu))**2 * Gel(ekq,tau12)/Gel(ek,tau12) * Dph(wq, tau12 ) )
      ![2.6] acceptance according to importance sampling 
      P_accept = D_AB * P_A/P_B 
      
      call random_number_omp(diagram%seed, ran )
      if(ran<P_accept) then 
         call delete_vertex(diagram,ivL,ivR)
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
      integer :: nu,m1,n1,m2,n2, nup
      real(dp) :: Pnu(Nph), tau
      type(vertex), pointer :: v1,v2
      complex(dp) :: Dratio

      if(diagram%order.eq.0) return 
      !randomly choose ph line 
      call uniform_int_omp(diagram%seed,1,diagram%order,iv1)
      v1 => diagram%vertexList(iv1); iv2 = v1%link(2)
      v2 => diagram%vertexList(iv2)
      tau = abs( v1%tau - v2%tau )

      m1 = v1%n_out; n1 = v1%n_in 
      m2 = v2%n_out; n2 = v2%n_in 
      do nu = 1,Nph 
         Pnu(nu) = abs( v1%gkq(m1,n1,nu)*v2%gkq(m2,n2,nu)*Dph(v1%wq(nu),tau) )
      end do 
      Pnu = Pnu / sum(Pnu)
      call Discrete_omp(diagram%seed,Pnu,Nph,nu)           !sampling wrt |D|, acceptance is 1 
      if(nu>Nph) return
      if(nu.eq.v1%nu) return 
      nup = v1%nu
      !update phase change 
      Dratio = v1%gkq(m1,n1,nu)*v2%gkq(m2,n2,nu) / (v1%gkq(m1,n1,nup)*v2%gkq(m2,n2,nup))

      diagram%gfac = diagram%gfac * Dratio / abs(Dratio)
      v1%nu = nu; v2%nu = nu

      stat%update(3) = stat%update(3) + 1
      stat%accept(3) = stat%accept(3) + 1

   end subroutine 

   !4 this may change the phase  
   subroutine change_band(diagram, stat)
      use DiagMC 
      implicit none 
      type(fynman), target :: diagram
      type(diagmc_stat) :: stat
      integer :: iv1,iv2
      integer :: nu1,nu2,m1,n1,m2,n2, ib
      real(dp) :: Pb(nbnd), tau
      type(vertex), pointer :: v1,v2
      complex(dp) :: Dratio

      !randomly choose an internal electron line 
      if(diagram%order.eq.0) return

      call uniform_int_omp(diagram%seed,1,diagram%order,iv1) 
      if(iv1.eq.diagram%tail) return  !external line is unchanged 
      v1 => diagram%vertexList(iv1); iv2 = v1%link(3)

      v2 => diagram%vertexList(iv2)
      tau = v2%tau - v1%tau 
      m1 = v1%n_out; n1 = v1%n_in 
      m2 = v2%n_out; n2 = v2%n_in 
      nu1 = v1%nu; nu2 = v2%nu 
      if(m1.ne.n2) stop 'band inconsistent @ change_band()'
      !write(*,*)m1,n1,nu1,m2,n2,nu2
      do ib = 1,nbnd 
         Pb(ib) = Gel(v1%ekout(ib),tau) * abs( v1%gkq(ib,n1,nu1)*v2%gkq(m2,ib,nu2) )
      end do 
      Pb = Pb / sum(Pb); ib = 1;
      call Discrete_omp(diagram%seed, Pb, nbnd, ib)
      v1%n_out = ib; v2%n_in = ib 
      !update phase change 
      Dratio = v1%gkq(ib,n1,nu1)*v2%gkq(m2,ib,nu2) / ( v1%gkq(m1,n1,nu1)*v2%gkq(m2,n2,nu2) )
      diagram%gfac = diagram%gfac * (Dratio / abs(Dratio))

      stat%update(4) = stat%update(4) + 1
      stat%accept(4) = stat%accept(4) + 1
   end subroutine
    
   !5
   subroutine update_vertex_time(diagram, stat)
      use DiagMC 
      implicit none 
      type(fynman), target :: diagram
      type(diagmc_stat) :: stat
      integer :: iv1
      type(vertex), pointer :: vA
      real(dp) :: tauL, tauR, tauQ,tau1, decay, Pt
      real(dp) :: wq,ran, Paccept, sign_wq 

      if(diagram%order<=2) return 

      stat%update(5) = stat%update(5) + 1
      call uniform_int_omp(diagram%seed,1,diagram%order,iv1) !choose a vertex    
      if(iv1.eq.diagram%tail) return 
      if(iv1.eq.diagram%head) return 

      vA => diagram%vertexList(iv1)
      tauL = diagram%vertexList(vA%link(1))%tau
      tauR = diagram%vertexList(vA%link(3))%tau
      tauQ = diagram%vertexList(vA%link(2))%tau

      wq = vA%wq(vA%nu)
      sign_wq = -1 
      if(vA%tau>tauQ) sign_wq = 1 
      decay = vA%ekin(vA%n_in) - vA%ekout(vA%n_out) + sign_wq * wq 
      call exp_sample_omp(diagram%seed,tauL,tauR,decay,tau1,Pt,1)
      Paccept = Dph( wq, abs(tau1-tauQ) ) / Dph( wq, abs(vA%tau-tauQ) ) !ph gf change 
      Paccept = Paccept * exp( sign_wq*wq*(tau1-vA%tau) ) 
      call random_number_omp(diagram%seed,ran)
      if(ran < Paccept) then 
         vA%tau = tau1
         stat%accept(5) = stat%accept(5) + 1
      end if  


   end subroutine 
    
   !6 
   subroutine update_external(diagram, stat)
      use DiagMC 
      implicit none 
      type(fynman), target :: diagram
      type(diagmc_stat) :: stat
      real(dp) :: tauL, decay, newT, oldT, Pt,ek,wq,ran,P_accept 
      type(vertex), pointer :: vR,vL 
      integer :: ib 
      real(dp) :: P_band(nbnd),tau_ph
      complex(dp) :: phase 

      ![1] update external time 
      vR => diagram%vertexList(diagram%tail)
      if(diagram%order.eq.0) then 
         call random_number_omp(diagram%seed, vR%tau)
         vR%tau = vR%tau * config%tauMax
         return 
      endif 
      ek = vR%ekin(vR%n_in)
      wq = vR%wq(vR%nu) 
      decay = ek + wq 
      tau_ph = diagram%vertexList(vR%link(2))%tau  
      tauL = diagram%vertexList(vR%link(1))%tau 
      if(tauL<config%tauMin) tauL = config%tauMin
      call exp_sample_omp(diagram%seed,tauL,config%tauMax, decay, newT, Pt,1)

      oldT = diagram%vertexList(diagram%tail)%tau
      P_accept = Dph(wq,newT-tau_ph) / Dph(wq,oldT-tau_ph) * exp( wq*(newT-oldT) )
      call random_number_omp(diagram%seed, ran)
      stat%update(6) = stat%update(6) + 1
      if(ran < P_accept ) then 
         stat%accept(6) = stat%accept(6) + 1
         diagram%vertexList(diagram%tail)%tau = newT
      endif
           

      ![2] change tail band index 
      if(nbnd.eq.1) return 
      do ib = 1,nbnd
         P_band(ib) = abs(vR%gkq(ib,vR%n_in,vR%nu))
      enddo 
      P_band = P_band/sum(P_band)
      call Discrete_omp(diagram%seed,P_band,nbnd,ib)
      if(ib>nbnd) return 
      !update phase 
      phase = vR%gkq(ib,vR%n_in,vR%nu) / vR%gkq(vR%n_out,vR%n_in,vR%nu)
      diagram%gfac = diagram%gfac * phase / abs(phase)
      vR%n_out = ib

      ![3] change head band index  
      vL => diagram%vertexList(diagram%head)
      do ib = 1,nbnd
         P_band(ib) = abs(vL%gkq(vL%n_out,ib,vL%nu))
      enddo 
      P_band = P_band/sum(P_band)
      call Discrete_omp(diagram%seed,P_band,nbnd,ib)
      if(ib>nbnd) return 
      !update phase 
      phase = vL%gkq(vL%n_out,ib,vL%nu) / vL%gkq(vL%n_out,vL%n_in,vL%nu)
      diagram%gfac = diagram%gfac * phase / abs(phase)
      vL%n_in = ib




   end subroutine 

   !7 this update may alternate the sign
   !this update can change gfac to complex number 
   !this may change the head & tail vertex 
   subroutine update_swap(diagram, stat)
      use DiagMC 
      use diagMC_set_debug, only : check_diagram
      implicit none 
      type(fynman), target :: diagram
      type(diagmc_stat) :: stat
      integer :: iv1, iv2, newtail, newhead
      type(vertex), pointer :: vL,vR

      real(dp) :: tau1, tau2, tau12, decay, wq1, wq2
      real(dp) :: P_accept, P_kchange
   
      real(dp) :: ran, k_ex(3)
      complex(dp) ::  factor 
      type(vertex),pointer :: vLn,vRn

      if(diagram%order<=2)return 

      !randomly choose a pair of vertex 
      call uniform_int_omp(diagram%seed,1,diagram%order,iv1)
      if(iv1.eq.diagram%tail) return          !no electron line after the tail vertex 
      if(iv1.eq.diagram%head) return          !no electron line after the tail vertex 
      vL => diagram%vertexList(iv1);  iv2 = vL%link(3)
      if(vL%link(2).eq.iv2) return            !the two vertex are not phonon linked 
      

      !count update 
      stat%update(7) = stat%update(7) + 1

      vR => diagram%vertexList(iv2)
      !if(vL%n_out.ne.vR%n_out) return            !the two vertex are not phonon linked 

      tau1 = vL%tau 
      tau2 = vR%tau 
      if(tau1>tau2) stop ' time disordered @ update_swap'
      tau12 = tau2 -tau1
      vLn => diagram%vertexList(diagram%order+1) 
      vRn => diagram%vertexList(diagram%order+2) 

      

      vRn%tau = vL%tau
      vRn%n_in = vL%n_in; vRn%n_out = vL%n_out; vRn%nu = vR%nu
      vRn%kin = vL%kin; vRn%q = vR%q
      vRn%kout = vRn%kin + vRn%q
      vRn%ekin = vL%ekin; vRn%ukin = vL%ukin;   
      
      k_ex = diagram%vertexList(diagram%head)%kin
      if(Equal(vRn%kout,k_ex)) return !check whether there are reducible diagram 

      call cal_ek(vRn%kout, vRn%ekout,vRn%ukout)
      vRn%wq = vR%wq; vRn%uq = vR%uq
      call cal_gkq_vtex( vRn, vRn%gkq ) !new e-ph vertex 

      vLn%tau = vR%tau
      vLn%n_in = vR%n_in; vLn%n_out = vR%n_out; vLn%nu = vL%nu
      vLn%kin = vRn%kout; vLn%q = vL%q
      vLn%kout = vLn%kin + vLn%q
      vLn%ekin = vRn%ekout; vLn%ukin = vRn%ukout 
      vLn%ekout = vR%ekout; vLn%ukout = vR%ukout 
      vLn%wq = vL%wq; vLn%uq = vL%uq
      call cal_gkq_vtex( vLn, vLn%gkq ) !new e-ph vertex

      !electron propagator ratio 
      P_kchange = Gel(vRn%ekout(vRn%n_out),tau12 ) / Gel(vL%ekout(vL%n_out),tau12 )

      !ph propagator ratio 
      wq1 = vL%wq(vL%nu); wq2 = vR%wq(vR%nu)
      
      P_kchange = P_kchange * Dph( wq1, abs( vR%tau - diagram%vertexList(vL%link(2))%tau ) ) /&
                  Dph( wq1, abs( vL%tau - diagram%vertexList(vL%link(2))%tau ) )
      P_kchange = P_kchange * Dph( wq2, abs( vL%tau - diagram%vertexList(vR%link(2))%tau ) ) /&
                  Dph( wq2, abs( vR%tau - diagram%vertexList(vR%link(2))%tau ) )
      factor = ( vRn%gkq(vRn%n_out,vRn%n_in,vRn%nu) * vLn%gkq(vLn%n_out,vLn%n_in,vLn%nu) ) & 
               / ( vL%gkq(vL%n_out,vL%n_in,vL%nu) * vR%gkq(vR%n_out,vR%n_in,vR%nu) )    

      P_accept = abs(factor) * P_kchange 
      !write(*,'(3f10.5)') P_accept, abs(factor),P_kchange
      call random_number_omp(diagram%seed,ran)
      
      if(ran<P_accept) then 
         !check this updates, it may be wrong because of contratry trend 
         factor = factor / abs(factor)
         diagram%gfac = diagram%gfac * factor 

         vLn%link = vL%link;       vRn%link = vR%link 
         vLn%link(1) = iv2;        vLn%link(3) = vR%link(3) 
         vRn%link(1) = vL%link(1); vRn%link(3) = iv1  
         diagram%vertexList(vL%link(1))%link(3) = iv2
         diagram%vertexList(vR%link(3))%link(1) = iv1
         
         call copy_vtex(Nph, dmc_band, nbnd,vLn,diagram%vertexList(iv1))
         call copy_vtex(Nph, dmc_band, nbnd,vRn,diagram%vertexList(iv2))
         !update tail/head vertex
         newtail = diagram%tail; newhead = diagram%head

         if(iv1.eq.diagram%tail) newtail = iv2
         if(iv2.eq.diagram%tail) newtail = iv1

         if(iv1.eq.diagram%head) newhead = iv2
         if(iv2.eq.diagram%head) newhead = iv1

         diagram%tail = newtail ; diagram%head = newhead
         stat%accept(7) = stat%accept(7) + 1
         !write(*,*) 'check phase @ swap'
         !call check_diagram(diagram)
         !write(*,*) 'check phase @ swap pass'
      end if 
      !call CleanVertex(vLn)
      !call CleanVertex(vRn)
   end subroutine 
    
   !ForNormalization 
   subroutine ForNormalization(diagram, stat)
      use random_tool, only : random_number_omp
      implicit none 
      type(fynman), target :: diagram
      type(diagmc_stat) :: stat
      real(dp) :: ran

      call random_number_omp(diagram%seed,ran)
      if(ran>0.5d0) then 
         call ToFake(diagram,stat)
      else 
         call ToReal(diagram,stat)
      end if
   end 

   !from 1st order self-energy to bare green 
   subroutine ToFake(diagram, stat)
      implicit none 
      type(fynman), target :: diagram
      type(diagmc_stat) :: stat

      integer :: ivL, ivR, ivn1, ivn2, nu 
      real(dp) :: tauLR,tau1, tau2, tau12
      real(dp) :: Pq, Pnu(Nph),Ptau12, Pib, Pel_line, temp  !P(B)
      real(dp) :: P_A,P_B,D_AB, P_accept
      real(dp) :: decay, decay_exp
      real(dp) :: ran, err_add
      complex(dp) :: gkq(nbnd,nbnd,Nph)
      type(vertex),pointer :: v1,v2
      type(vertex),pointer :: vL, vR

      if(diagram%order.ne.2) return  

      vL => diagram%vertexList(diagram%head)
      vR => diagram%vertexList(diagram%tail)
      if(vL%n_in.ne.vR%n_out) return 
      if(vL%n_in.ne.1) return  !only the lower band 
      tauLR = vR%tau 
      if(vL%tau>1e-10) stop 'tail time not equal 0'
      ![1.1] P_B 
      Pq = 1.d0 
      do nu = 1,Nph 
         Pnu(nu) = abs(vL%gkq(vL%n_out,vL%n_in,nu))**2
      enddo 
      Pnu = Pnu / sum(Pnu)
      P_B = Pq * Pnu(vL%nu) / nbnd 
      ![1.2] D_AB
      D_AB = 1.d0 / ( abs(vL%gkq(vL%n_out,vL%n_in,vL%nu))**2 * Dph( vL%wq(vL%nu),tauLR )  )
      D_AB =  D_AB /( Gel(vL%ekout(vL%n_out),tauLR) / FakeRatio )
      P_accept = D_AB * P_B 
      !write(*,*) P_accept,D_AB,P_B
      call random_number_omp(diagram%seed,ran)
      if(ran < P_accept) then 
         vL%n_out = vL%n_in; vR%n_in = vL%n_out
         vL%q = 0.d0;  vR%q = 0.d0; vL%nu = -1; vR%nu = -1; 
         vL%kout = vL%kin; vL%ekout = vL%ekin; vL%ukout = vL%ukin       !el eigen mode
         vR%kin  = vR%kout;  vR%ekin = vR%ekout;  vR%ukin = vR%ukout 
         diagram%order = 0
      endif 
       

   end subroutine

   !from bare green to 1st order self-energy 
   subroutine ToReal(diagram, stat)
      implicit none 
      type(fynman), target :: diagram
      type(diagmc_stat) :: stat

      integer :: ivL, ivR, ivn1, ivn2, nu 
      real(dp) :: tauLR,tau1, tau2, tau12
      real(dp) :: Pq, Pnu(Nph),Ptau12, Pib, Pel_line, temp  !P(B)
      real(dp) :: P_A,P_B,D_AB, P_accept
      real(dp) :: decay, decay_exp
      real(dp) :: ran, err_add
      complex(dp) :: gkq(nbnd,nbnd,Nph)
      type(vertex),pointer :: v1,v2
      type(vertex),pointer :: vL, vR, vn1, vn2 


      if(diagram%order.ne.0) return 

      vn1 => diagram%vertexList(diagram%order+1)
      vL => diagram%vertexList(diagram%head)
      vR => diagram%vertexList(diagram%tail)
      tauLR = vR%tau - vL%tau

      ![1.1] choose ph momentum 
      call sample_q_omp(diagram%seed,vn1%q); Pq = 1.d0 !choose from pre-tabulated mesh 
      vn1%kin = vL%kin 
      vn1%kout = vn1%kin + vn1%q
      call cal_wq(vn1%q,vn1%wq,vn1%uq)
      vn1%ukin = vL%ukin; vn1%ekin = vL%ekin 
      call cal_ek(vn1%kout,vn1%ekout,vn1%ukout)
      call cal_gkq_symmetrize(vn1%kin, vn1%q, vn1%wq, vn1%uq, vn1%ukin, vn1%ukout, vn1%gkq)

      ![1.2] propose el/ph-index ~ abs(gkq)
      vn1%n_in = vL%n_in
      call uniform_int_omp(diagram%seed,1,nbnd,vn1%n_out)
      do nu = 1,Nph 
         Pnu(nu) = abs(vn1%gkq(vn1%n_out,vn1%n_in,nu))**2
      enddo 
      Pnu = Pnu / sum(Pnu) 
      call Discrete_omp(diagram%seed,Pnu, Nph, nu);
      if(nu>Nph) return    
      vn1%nu = nu    
      
      ![2.1] propose prob for q,m,nu
      P_A = Pq * Pnu(nu) / nbnd  
      ![2.2] 
      D_AB = abs(vn1%gkq(vn1%n_out,vn1%n_in,nu))**2 * Dph( vn1%wq(nu),tauLR ) 
      D_AB = D_AB * Gel(vn1%ekout(vn1%n_out),tauLR) / FakeRatio
      P_accept = D_AB / P_A
      !write(*,*) D_AB,P_A
      call random_number_omp(diagram%seed,ran)
      if(ran < P_accept) then 
         !write(*,*) 'nu = ',nu
         vL%link(2) = diagram%tail; vR%link(2) = diagram%head 
         vL%n_out = vn1%n_out; vL%nu = vn1%nu                           !el/ph index : m,n,\nu  
         vR%n_in = vn1%n_out;  vR%nu = vn1%nu 
         vL%q = vn1%q;  vL%wq = vn1%wq; vL%uq = vn1%uq                  !ph eigen mode
         vR%q = -vn1%q; vR%wq = vn1%wq; vR%uq = conjg(vn1%uq)
         vL%kout = vn1%kout; vL%ekout = vn1%ekout; vL%ukout = vn1%ukout !el eigen mode
         vR%kin  = vn1%kout;  vR%ekin = vn1%ekout;  vR%ukin = vn1%ukout 
         vL%gkq = vn1%gkq                                               !e-ph matrix     
         do nu = 1,Nph 
            vR%gkq(:,:,nu) = conjg(transpose(vL%gkq(:,:,nu)))
         end do 
         diagram%order = 2
      endif 
      

   end subroutine


   !9 change momentum  
   subroutine Change_K(diagram, stat)
      use mod_HTable, only : random_Q_omp
      use DiagMC, only : Holstein
      implicit none 
      type(fynman), target :: diagram
      type(diagmc_stat) :: stat
      type(vertex),pointer :: vL, vR

      if(diagram%order.ne.0) return  

      vL => diagram%vertexList(diagram%head)
      vR => diagram%vertexList(diagram%tail)      

      call sample_q_omp(diagram%seed,vL%kin)

      vL%kout = vL%kin
      vR%kin = vL%kin
      vR%kout = vL%kin
      call random_number_omp(diagram%seed, vR%tau)
      vR%tau = vR%tau * config%tauMax
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
       
      call move_vertex(diagram,iv1,diagram%order+1) ! 1 -> order+1
      call move_vertex(diagram,iv2,iv1)             ! 2 -> 1
      call move_vertex(diagram,diagram%order+1,iv2) !order + 1 -> 2
      
      !swap may change the position of the tail/head, keep track of that 
      if(iv1.eq.diagram%tail) diagram%tail = iv2
      if(iv1.eq.diagram%head) diagram%head = iv2
      if(iv2.eq.diagram%tail) diagram%tail = iv1
      if(iv2.eq.diagram%head) diagram%head = iv1

   end subroutine 

   subroutine move_vertex(diagram,iv1,iv2)
      !move vertex iv1 to iv2 -th place 
      use DiagMC, only : nbnd,Nph
      implicit none 
      type(fynman) :: diagram
      integer :: iv1,iv2,j_tmp
      integer :: link(3)

      call copy_vtex(Nph, dmc_band, nbnd,diagram%vertexList(iv1),diagram%vertexList(iv2))

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
   
      call swap_vertex(diagram,iv1,diagram%order-1)
      call swap_vertex(diagram,iv2,diagram%order)

      v1 => diagram%vertexList(diagram%order-1)
      v2 => diagram%vertexList(diagram%order)
      if( v1%link(3).ne.v1%link(2) ) stop 'link err@delete_vertex'
   
      v1_in = v1%link(1); v2_out = v2%link(3)
      diagram%vertexList(v1_in)%link(3) = v2_out
      diagram%vertexList(v2_out)%link(1) = v1_in

      diagram%order = diagram%order - 2           !then the last two elements are forgetten 
      diagram%vertexList(diagram%order+1)%tau = 0.0
      diagram%vertexList(diagram%order+2)%tau = 0.0
        
   end subroutine 

   !init self-energy diagram, start form the fake one. 
   subroutine init_selfe_diagram(diagram)
      implicit none 
      type(fynman), target :: diagram
      type(vertex), pointer :: vi 
      real(dp) :: tau 

      tau = 1.0    
      if(tau<config%tauMin .or. tau>config%tauMax) tau = 0.5*(config%tauMin + config%tauMax)

      !init diagram 
      diagram%head = 1; diagram%tail = 2;
      diagram%order = 0; 
      diagram%gfac = 1.d0; 

      vi => diagram%vertexList(1); vi%tau = 0.d0
      vi%kin  = config%P;   vi%ekin = config%EP
      vi%kout = config%P;  vi%ekout = config%EP
      vi%link = (/ maxN, -1, 2 /)
      call cal_ek(vi%kout,vi%ekout,vi%ukout); vi%ukin = vi%ukout 
      vi%n_in = 1; vi%n_out = 1

      vi => diagram%vertexList(2); vi%tau = tau
      vi%kin  = config%P;  vi%ekin = config%EP 
      vi%kout = config%P; vi%ekout = config%EP
      vi%link = (/ 1,-1,maxN /)
      call cal_ek(vi%kin,vi%ekin,vi%ukin); vi%ukout = vi%ukin 
      vi%n_in = 1 ; vi%n_out = 1

   end subroutine

end module 

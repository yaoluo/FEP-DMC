!------------------------------------------------------------------------
!diagrammatic monte carlo for electron phonon Hamiltonian
!key subroutines are the 7 updates: 
!1. add_ph : add phonon line on one electron line 
!2. remove_ph : the complementary update of add_ph
!3. change_ph_mode: change the phonon mode 
!4. change_band : change internal band index 
!5. update_vertex_time : move the time of e-ph vertex 
!6. swap: swap nearby two vertex in time 
!author : Yao Luo, Email: yluo7@caltech.edu
!
!add_ph() & swap() are the most expensive (other updates cost no time in comparsion)  
! because it needs to evaluate new e-ph matrix. 
!For efficiency, the prob for the two updates should be comparablely lower than the others. 
!Version 2.0 @  April 2022 : 
!     changed every updates with fynman diagrams and statistics as input for openmp parallelization. 
!Version 3.0 @ May 2022 : 
!     fixed external tau = \beta for current-current correlator 
!------------------------------------------------------------------------


module JJt_update  
   use DiagMC, only : config,fynman,  diagmc_stat, vertex, copy_vtex,maxN 

   contains 

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
      case(5)
         call update_vertex_time(diagram,stat)
      case(6) 
         call update_swap(diagram,stat)
      end select 
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
      !propose ph-momentum 
      call sample_q_omp(diagram%seed,vn1%q); Pq = 1.d0 !choose from pre-tabulated mesh 
      
     
      vn1%kin = vL%kout; 
      vn1%kout = vn1%kin + vn1%q
      call cal_ek( vn1%kin,   vn1%ekin ,  vn1%ukin  )
      
      !vn1%ekin = vL%ekout; vn1%ukin = vL%ukout
      call cal_ek( vn1%kout,  vn1%ekout,  vn1%ukout )
      call cal_wq( vn1%q,     vn1%wq,     vn1%uq    )
      call cal_gkq_symmetrize(vn1%kin, vn1%q, vn1%wq, vn1%uq, vn1%ukin, vn1%ukout, vn1%gkq)
      !propose el/ph-index ~ abs(gkq)
      vn1%n_in = vL%n_out
      call uniform_int_omp(diagram%seed,1,nbnd,vn1%n_out)
      do nu = 1,Nph 
         Pnu(nu) = abs(vn1%gkq(vn1%n_out,vn1%n_in,nu))**2
      enddo 
      Pnu = Pnu / sum(Pnu) 
      call Discrete_omp(diagram%seed,Pnu, Nph, nu);
      if(nu>Nph) return    
      vn1%nu = nu 

      !choose tau first P(t1,t2) ~ exp(-decay_exp*(t2-t1)), t2>t1
      decay_exp = vn1%wq(nu) + vn1%ekout(vn1%n_out) - vn1%ekin(vn1%n_in) 
      tauR = diagram%vertexList(ivR)%tau
      call uniform_real_omp(diagram%seed,tauL,tauR,tau1,temp,1)
      call exp_sample_omp(diagram%seed,tau1,tauR,decay_exp,tau2,Ptau12,1)
      Ptau12 = Ptau12 / (tauR-tauL)
      tau12 = tau2 - tau1

      !prob for choosing to remove this phonon line 
      P_A = 1.d0/( (diagram%order-1)/2 + 1 ) * config%PA(2)             
      !prob for choosing to add this phonon line 
      Pel_line = 1.d0 / diagram%order; 
      Pib = 1.d0 / nbnd
      P_B = Pq * Pel_line * Pib * Pnu(nu) * Ptau12 * config%PA(1)         
      !prob ratio W(B)/W(A) 
      D_AB = abs(vn1%gkq(vn1%n_out,vn1%n_in,nu))**2
      D_AB = D_AB * Gel(vn1%ekout(vn1%n_out),tau12) / Gel(vn1%ekin(vn1%n_in),tau12) !electron gf change 
      D_AB = D_AB * Dph(vn1%wq(nu),tau12)                                           !phonon gf change 
      P_accept = D_AB * P_A / P_B 

      !write(*,*) decay_exp,P_A,P_B
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

      if(diagram%order.eq.1)return 
      
      !randomly choose a ph-line 
      call uniform_int_omp(diagram%seed,2,diagram%order,iv1)
      vL => diagram%vertexList(iv1);  iv2 = vL%link(2)
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
      
      xqt = vL%q 
      tauL = diagram%vertexList(vL%link(1))%tau 
      tauR = diagram%vertexList(vR%link(3))%tau 
      tau1 = vL%tau; tau2 = vR%tau 
      ek = vL%ekin(vL%n_in); ekq = vL%ekout(vL%n_out); 
      nu = vL%nu; wq = vL%wq(nu)

      !tau1 tau2 distribution
      tau12 = tau2 - tau1  
      if(tau12<0) then 
         write(*,*)tau12
         stop 'time disorder'
         call mp_global_end() 
      end if 
      decay_exp = wq + ekq - ek
      call exp_sample_omp(diagram%seed,tau1,tauR,decay_exp,tau2,Ptau12)
      Ptau12 = Ptau12 / (tauR-tauL)

      !nu probility
      do nup = 1,Nph 
         Pnu(nup) = abs(vL%gkq(vL%n_out,vL%n_in,nup))**2
      end do 
      Pnu = Pnu / sum(Pnu)

      ![2.1] prob for adding ph line  
      Pq = 1.d0; 
      Pel_line = 1.d0 / (diagram%order - 2); 
      !Pel_line = 1.d0 / 3; 
      Pib = 1.d0/nbnd; 
      P_A = Pq * Pel_line * Pib * Pnu(nu) * Ptau12 * config%PA(1)
      ![2.2] prob for removing ph line 
      Pph_line = 1.d0 / ( (diagram%order-1)/2 )
      !Pph_line = 1.d0 / 2
      P_B = Pph_line * config%PA(2)      !there are (order-1)/2 phonon line 
      D_AB = 1.d0/( abs(vL%gkq(vL%n_out,vL%n_in,nu))**2 * Gel(ekq,tau12)/Gel(ek,tau12) * Dph(wq, tau12 ) )

      P_accept = D_AB * P_A / P_B 
      
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

      if(diagram%order.eq.1) return 
      !randomly choose ph line 
      call uniform_int_omp(diagram%seed,2,diagram%order,iv1)
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
      if(diagram%order.eq.1) return

      call uniform_int_omp(diagram%seed,2,diagram%order,iv1) 
      v1 => diagram%vertexList(iv1); iv2 = v1%link(3)
      if(iv2.eq.maxN) return !external line is unchanged 

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
      Pb = Pb / sum(Pb)
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
      !because of cyclic symmetry, vertex connected to tail or head needs special care
      !connectivity may change if vertex move across the time boundary 
      use DiagMC 
      implicit none 
      type(fynman), target :: diagram
      type(diagmc_stat) :: stat
      integer :: ivA,ivL,ivR
      type(vertex), pointer :: vA, vL,vR, vH,vT
      real(dp) :: tauL, tauR, tauLR, tauQ,tau1, decay, Pt, sign_wq
      real(dp) :: wq,ran, Paccept 
      logical :: cross

      if(diagram%order.eq.1) return 

      stat%update(5) = stat%update(5) + 1
      call uniform_int_omp(diagram%seed,2,diagram%order,ivA) !choose a vertex 
      vA => diagram%vertexList(ivA)
      wq = vA%wq(vA%nu)
      decay = vA%ekin(vA%n_in) - vA%ekout(vA%n_out)
      tauQ = diagram%vertexList(vA%link(2))%tau
      cross = .False.
      if(vA%link(1).eq.diagram%head) then 
         vH => diagram%vertexList(diagram%head)
         vT => diagram%vertexList(diagram%tail)
         ivR = vA%link(3)
         vR => diagram%vertexList(ivR)
         ivL = vT%link(1)
         vL => diagram%vertexList(ivL)

         tauL = vL%tau 
         tauR = vR%tau 
         tauLR = tauR - tauL + config%tauMax 
         call exp_sample_omp(diagram%seed,0,tauLR,decay,tau1,Pt,1)
         tau1 = tauL + tau1 
         if(tau1<config%tauMax) then 
            cross = .True.
         else 
            tau1 = tau1-config%tauMax 
         endif 
         Paccept = Dph( wq, abs(tau1-tauQ) ) / Dph( wq, abs(vA%tau-tauQ) )
         call random_number_omp(diagram%seed,ran)
         if(ran < Paccept) then 
            vA%tau = tau1
            stat%accept(5) = stat%accept(5) + 1
            !change connetivity and update head and tail 
            if(corss) then 
               vA%link(1) = ivL;  vL%link(3) = ivA
               vR%link(1) = diagram%head ; vH%link(3) = ivR
               vT%link(1) = ivA; vA%link(3) = diagram%tail
               vH%ekout = vA%ekout; vH%ukout = vA%ukout; vH%n_in =  vA%n_out
               vH%ekin =  vA%ekout; vH%ukin =  vA%ukout; vH%n_out = vA%n_out
               vT%ekout = vA%ekout; vT%ukout = vA%ukout; vT%n_in =  vA%n_out
               vT%ekin =  vA%ekout; vT%ukin =  vA%ukout; vT%n_out = vA%n_out
            endif
         end if 
         return 
      endif 

      if(vA%link(3).eq.diagram%tail) then 
         vH => diagram%vertexList(diagram%head)
         vT => diagram%vertexList(diagram%tail)
         ivR = vH%link(3)
         vR => diagram%vertexList(ivR)
         ivL = vA%link(1)
         vL => diagram%vertexList(ivL)

         tauL = vL%tau
         tauR = vR%tau 
         tauLR = tauR - tauL + config%tauMax
         call exp_sample_omp(diagram%seed,0,tauLR,decay,tau1,Pt,1)
         tau1 = tauL + tau1 
         if(tau1>config%tauMax) then 
            cross = .True.
            tau1 = tau1-config%tauMax 
         endif 
         Paccept = Dph( wq, abs(tau1-tauQ) ) / Dph( wq, abs(vA%tau-tauQ) )
         call random_number_omp(diagram%seed,ran)
         if(ran < Paccept) then 
            vA%tau = tau1
            stat%accept(5) = stat%accept(5) + 1
            !change connetivity and update head and tail 
            if(corss) then 
               vA%link(1) = diagram%head;  vH%link(3) = ivA
               vR%link(1) = ivA ; vA%link(3) = ivR
               vT%link(1) = ivL; vL%link(3) = diagram%tail
               vH%ekout = vA%ekin; vH%ukout = vA%ukin; vH%n_in = vA%n_in
               vH%ekin = vA%ekin; vH%ukin = vA%ukin; vH%n_out = vA%n_in
               vT%ekout = vA%ekin; vT%ukout = vA%ukin; vT%n_in = vA%n_in
               vT%ekin = vA%ekin; vT%ukin = vA%ukin; vT%n_out = vA%n_in
            endif
         end if 
         return  
      endif 

      tauL = diagram%vertexList(vA%link(1))%tau
      tauR = diagram%vertexList(vA%link(3))%tau
   
      if(tauL .eq. tauR) stop "tauLR equal @ update_vertex_time"
      if(decay .eq. 0.d0) stop "decay .eq. 0 @ update_vertex_time"
      call exp_sample_omp(diagram%seed,tauL,tauR,decay,tau1,Pt,1)


      Paccept = Dph( wq, abs(tau1-tauQ) ) / Dph( wq, abs(vA%tau-tauQ) ) !ph gf change 
      Paccept = Paccept !* exp( sign_wq*wq*(tau1-vA%tau) ) 

      !write(*,*) Paccept
      call random_number_omp(diagram%seed,ran)
      if(ran < Paccept) then 
         vA%tau = tau1
         stat%accept(5) = stat%accept(5) + 1
      end if 
   end subroutine 
    
    
   !6 this update may alternate the sign
   !this update can change gfac to complex number 
   subroutine update_swap(diagram, stat)
      use DiagMC 
      implicit none 
      type(fynman), target :: diagram
      type(diagmc_stat) :: stat
      integer :: iv1, iv2
      type(vertex), pointer :: vL,vR

      real(dp) :: tau1, tau2, tau12, decay, wq1, wq2
      real(dp) :: P_accept, P_kchange
   
      real(dp) :: ran
      complex(dp) ::  factor 
      type(vertex),pointer :: vLn,vRn

      if(diagram%order<=3)return 

      !randomly choose a pair of vertex 
      call uniform_int_omp(diagram%seed,2,diagram%order,iv1)
      vL => diagram%vertexList(iv1);  iv2 = vL%link(3)
      if(vL%link(2).eq.iv2) return            !the two vertex are not phonon linked 
      if(iv2.eq.maxN)       return            !not the last vertex   
      !count update 
      stat%update(7) = stat%update(7) + 1

      vR => diagram%vertexList(iv2)
      tau1 = vL%tau 
      tau2 = vR%tau 
      if(tau1>tau2) stop ' time disordered @ update_swap'
      tau12 = tau2 -tau1
      vLn => diagram%vertexList(diagram%order+1) 
      vRn => diagram%vertexList(diagram%order+2) 
      !call CreateVertex(Nph,nbnd,vLn)
      !call CreateVertex(Nph,nbnd,vRn)

      vRn%tau = vL%tau
      vRn%n_in = vL%n_in; vRn%n_out = vL%n_out; vRn%nu = vR%nu
      vRn%kin = vL%kin; vRn%q = vR%q
      vRn%kout = vRn%kin + vRn%q
      vRn%ekin = vL%ekin; vRn%ukin = vL%ukin;     
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
      factor = (vRn%gkq(vRn%n_out,vRn%n_in,vRn%nu) * vLn%gkq(vLn%n_out,vLn%n_in,vLn%nu) ) & 
               / ( vL%gkq(vL%n_out,vL%n_in,vL%nu) * vR%gkq(vR%n_out,vR%n_in,vR%nu) )    

      P_accept = abs(factor) * P_kchange  

      call random_number_omp(diagram%seed,ran)
      
      if(ran<P_accept) then 
         factor = factor / abs(factor)
         diagram%gfac = diagram%gfac * factor 

         vLn%link = vL%link;       vRn%link = vR%link 
         vLn%link(1) = iv2;        vLn%link(3) = vR%link(3) 
         vRn%link(1) = vL%link(1); vRn%link(3) = iv1  
         diagram%vertexList(vL%link(1))%link(3) = iv2
         diagram%vertexList(vR%link(3))%link(1) = iv1
         
         !diagram%vertexList(iv1) = vLn 
         !diagram%vertexList(iv2) = vRn 
         call copy_vtex(Nph, dmc_band, nbnd,vLn,diagram%vertexList(iv1))
         call copy_vtex(Nph, dmc_band, nbnd,vRn,diagram%vertexList(iv2))
         stat%accept(7) = stat%accept(7) + 1
         !write(*,*) 'check phase @ swap'
         !call check_phase(diagram)
         !write(*,*) 'check phase @ swap pass'
      end if 
      !call CleanVertex(vLn)
      !call CleanVertex(vRn)
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
   end subroutine 

   subroutine move_vertex(diagram,iv1,iv2)
      !move vertex iv1 to iv2 -th place 
      use DiagMC, only : nbnd,Nph,dmc_band
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

end module 


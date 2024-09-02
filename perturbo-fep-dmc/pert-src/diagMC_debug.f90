
module diagMC_gt_debug
   use diagmc
   contains 
   
   subroutine check_diagram(diagram)
       use diagMC 
       implicit none
       type(fynman), target :: diagram  
       integer :: iv, ib, nu 
       type(vertex), pointer :: v1 , v2 
       character(len=50) :: line 
       real(dp) :: ek(dmc_band), wq(Nph), ks(3), errU, errG, vnk(3,dmc_band)
       complex(dp) :: uq(Nph,Nph), gkq(dmc_band,dmc_band,Nph),iden(dmc_band,dmc_band)

      iden = 0.d0 
      do ib = 1,dmc_band;
         iden(ib,ib) = 1.d0 
      enddo 
   
       
       !check link 
       do iv = 2,diagram%order 
           v1 => diagram%vertexList(iv)
           v2 => diagram%vertexList(v1%link(1))
           if( v2%link(3).ne.iv ) stop 'diagram err : link err 1'
           v2 => diagram%vertexList(v1%link(2))
           if( v2%link(2).ne.iv ) stop 'diagram err : link err 2'
           v2 => diagram%vertexList(v1%link(3))
           if( v2%link(1).ne.iv ) stop 'diagram err : link err 3'
           if( maxval(v1%link)>diagram%order .and. maxval(v1%link).ne.maxN ) stop 'link num > diagram%order'
       end do 
       
       !check tau diagram%order 
       do iv = 1,diagram%order 
           v1 => diagram%vertexList(iv)
           v2 => diagram%vertexList(v1%link(3)) 
           if(v2%tau < v1%tau ) stop 'diagram err : vertex not time diagram%ordered'
           if(iv.eq.1) cycle
           v2 => diagram%vertexList(v1%link(1)) 
           if(v2%tau > v1%tau ) stop 'diagram err : vertex not time diagram%ordered'
       end do 
   
       !check mometum conservation 
       do iv = 2,diagram%order 
           v1 => diagram%vertexList(iv)
           ks = v1%i_kin+v1%i_q-v1%i_kout
           !call To1BZ(ks)
           if( norm2(ks)>1e-8 ) stop 'diagram err : momentum conservation breaks at vertex'
           v2 => diagram%vertexList(v1%link(1)) 
           ks = v2%i_kout-v1%i_kin
           !call To1BZ(ks)
           if( norm2(ks)>1e-8 ) stop 'diagram err : momentum conservation breaks at electron bond'
           v2 => diagram%vertexList(v1%link(2)) 
           ks = v2%i_q+v1%i_q
           !call To1BZ(ks)
           if( norm2(ks)>1e-8 ) stop 'diagram err : momentum conservation breaks at phonon bond'
       end do 
   
   
       !check dispersion , el-ph vertex
       do iv = 1,diagram%order
            v1 => diagram%vertexList(iv)
            call cal_ek_int(v1%i_kout, v1%ekout, v1%ukout, vnk)
            !if(maxval(abs(vnk-v1%vkout))>1e-10) then 
            !   stop 'diagram err : velocity inconsistent'
            !endif  
            v2 => diagram%vertexList(v1%link(3)) 
            errU = maxval(abs(v1%ukout-v2%ukin))
            if(errU>1E-10) then
               v2%ukin = v1%ukout
               v2%i_kin = v1%i_kout
               if(errU>1e-3) then    
                  write(*,'(3i5,E20.10)') iv, v1%link(3),diagram%order,errU
                  write(*,'(3E15.5)')  v1%ekout
                  write(*,'(3E15.5)')  v2%ekin
               stop 'diagram err : uk inconsistent'
               endif 
            endif 
            errU = maxval(abs(matmul(transpose(conjg(v1%ukout)),v1%ukout) - iden))
            if(errU>1E-7) stop 'diagram err : uk not unitary!'
       end do 

      !check gkq 
      do iv = 2, diagram%order
         v1 => diagram%vertexList(iv)
         call cal_gkq_vtex_int(v1, gkq)
         errU = maxval(abs(gkq - v1%gkq))
         if(errU>1e-6) then
            write(*,*) maxval(abs(gkq)), maxval(abs(v1%gkq))
            write(*,*) errU 
            stop 'gkq err'
         endif 
      enddo 
      
      !check gkq and gkq conj
      do iv = 2,diagram%order
         v1 => diagram%vertexList(iv)
         v2 => diagram%vertexList(v1%link(3)) 
         if(v1%link(3).eq.v1%link(2)) then 
            do nu = 1,Nph 
               errU = maxval(abs(v1%gkq(:,:,nu)-transpose(conjg(v2%gkq(:,:,nu)))))
               if(errU>1e-5) then 
                  if(errU>1E-7) stop 'diagram err : gkq, gk+q,-q not conjguate to each other'
               endif 
            enddo 
         endif 
      end do 


   
   end subroutine 
   
   subroutine check_epw()
      use diagMC 
      implicit none 
      real(dp) :: xkt(3), xqt(3), errg
      complex(dp) ::  gkq(dmc_band,dmc_band,Nph), gkq_2(dmc_band,dmc_band,Nph)
      integer :: i,j,itr  
      type(vertex) :: v1,v2 
   
      call start_clock('check_epw')
      call CreateVertex(Nph,dmc_band,nbnd,nsvd,v1); call CreateVertex(Nph,dmc_band,nbnd,nsvd,v2)
      do itr = 1,100 
   
         v1%i_kin = xkt!(/-0.40796E-01,   -0.12336E+00   ,-0.17060E-02/)
         v1%i_q = xqt!(/0.20337E-03  , -0.94166E-03 ,  -0.88046E-03/)
         v1%i_kout = v1%i_kin + v1%i_q
         
         v2%i_kin = v1%i_kout; v2%i_q = -v1%i_q; v2%i_kout = v1%i_kin 
   
         call cal_wq_int(v1%i_q,v1%wq)
         call cal_ek_int(v1%i_kin,v1%ekin,v1%ukin)
         call cal_ek_int(v1%i_kout,v1%ekout,v1%ukout)
   
         call cal_gkq_vtex_int(v1,gkq)
   
         v2%wq = v1%wq; 
         v2%ukin = v1%ukout; v2%ukout = v1%ukin; v2%uq = conjg(v1%uq)
         call cal_gkq_vtex_int(v2,gkq_2)
         do i = 1,Nph 
            gkq_2(:,:,i) = conjg(transpose(gkq_2(:,:,i)))
            !the mesh for uq has inversion-symemtry related ones. 
            !This can introduce discrependcies 
            !need to check 
         end do 
         errg = maxval(abs(gkq - gkq_2 ))/maxval(abs(gkq))
         !write(*,'(A15,2E15.5)')'maxg = ',maxval(abs(gkq))**2,maxval(abs(gkq_2))**2
         if(errg > 1e-5) then 
            write(*,'(A15,2E15.5)')'maxg = ',maxval(abs(gkq))**2,maxval(abs(gkq_2))**2
            write(*,'(9f10.5)')v1%i_kin,v1%i_kout,v1%i_q
            write(*,'(9f10.5)')v2%i_kin,v2%i_kout,v2%i_q
            write(*,'(A15,2E15.5)')'maxg = ',maxval(abs(gkq)),maxval(abs(gkq_2))
            write(*,'(A15,E15.5)')'errg = ',maxval(abs(gkq - gkq_2 ))/maxval(abs(gkq))
            write(*,'(A15,E15.5)')'errabsg = ',maxval(abs(abs(gkq) - abs(gkq_2)) )/maxval(abs(gkq))
            stop 'err@check_epw()'
         endif 
      enddo
     
      
      !write(*,'(A15,2E15.5)')'maxg = ',maxval(abs(gkq)),maxval(abs(gkq_2))
      call stop_clock('check_epw')
      if(maxval(abs(gkq - gkq_2 ))/maxval(abs(gkq)) < 1e-10) then 
         call CleanVertex(v1); call CleanVertex(v2)
         return 
      end if 
      if(Nph<=6) then 
         write(*,'(A15,6E15.5)')'w = ',v1%wq
         write(*,'(A15,6E15.5)')'imag g1g2 = ',imag(gkq*gkq_2)
         write(*,'(A15,6E15.5)')'imag g1 = ',imag(gkq)
         write(*,'(A15,6E15.5)')'imag g2 = ',imag(gkq_2)
         write(*,'(A15,6E15.5)')'real g1 = ',real(gkq)
         write(*,'(A15,6E15.5)')'real g2 = ',real(gkq_2)
         write(*,'(A15,6E15.5)')'abs g1 = ',abs(gkq)
         write(*,'(A15,6E15.5)')'abs g2 = ',abs(gkq_2)
      end if
      stop 'g_{k,q} != conjg( g_{k+q,-q} )' 
   
   end subroutine 
   
   
   subroutine check_phase(diagram)
      implicit none 
      type(fynman),target :: diagram
      integer :: iv, ie, jv, itr ,final_vertex 
      real(dp) :: dtau, factor 
      complex(dp) :: gfac_exact, vn(dmc_band), right_env(dmc_band, dmc_band), left_env(dmc_band,dmc_band),A(dmc_band,dmc_band), x
      type(vertex), pointer :: vi,vj

      
      vi => diagram%vertexList(1)
      vj => diagram%vertexList(vi%link(3))
      call right_environment_matrix(diagram, 1, right_env )
      

      x = 0.d0 
      do ie = 1,dmc_band 
         x = x + right_env(ie,ie)
      enddo 
      factor = real(x)
      factor = factor / abs(factor)
      if( abs(factor-diagram%gfac)>1e-5 ) then 
         write(*,*)'sign miss match!'
         write(*,'(2f10.5)')factor, real(diagram%gfac)
         stop 
      endif 

    

   end subroutine 
   
   !for sanity checking, direct calcuate the green function at leading diagram%order 
   !do this only for a single k-point
   !compute the lowest diagram%order self-energy in matasubrara formulism at a single kpoints 
   subroutine MataSelfEnergy()
      use DiagMC, only : config,ph,el,ep,kbT,nq_se_check
      use selfenergy
      use pert_utils, only : random_uniform, find_free_unit 
      use pert_param, only : ryd2ev
      implicit none 
      integer,parameter :: Nmata = 10000  
      integer :: nq, bmin  
      ! selfe_i(Nmata, ntmpr, nband, nkpt)
      complex(dp),allocatable :: selfE(:,:,:,:)
      
      
      real(dp) :: kpts(3,1), tmpr(1), ef(1), egrid(Nmata)
      real(dp),allocatable :: qset(:,:), qwt(:)
      real(DP) :: enk(el%nb,1)
      integer :: iw, iunit  
   
   
      kpts(:,1) = config%P 
      do iw = -Nmata/2,Nmata/2-1 
         egrid(iw+1+Nmata/2) = (2*iw+1) * pi*kbT
      enddo 
      nq = nq_se_check
      ALLOCATE(qset(3,nq),qwt(nq),selfE(Nmata,1,el%nb,1))
      CALL random_uniform(qset,qwt)
      
      tmpr(1) = kbT 
      ef(1) = config%mu / ryd2ev 
      bmin = 1
      call selfenergy_mata(ep, el, ph, kpts, egrid, qset, qwt, tmpr, ef, bmin, enk, selfE )  
      !output 
      if(ionode) then 
         iunit = find_free_unit()
         open(iunit,file = 'SelfE_sanity-check.dat')
         write(iunit,'(3A15)') '#Mata freq','Re SelfE','Im SelfE'
         do iw = 1,Nmata 
            write(iunit,'(3E15.5)')egrid(iw),selfE(iw,1,config%ib,1)
         enddo  
         close(iunit)
      endif 
   
      DEALLOCATE(qset,qwt,selfE)
   
   end subroutine


   
end 


module diagMC_gnt_debug
   use diagmc
   contains 
   
   subroutine check_diagram(diagram)
      use diagMC 
      use pert_param, only : nk_svd
      use HamiltonianTable, only: find_ik
      implicit none
      type(fynman), target :: diagram  
      integer :: iv, ib, nu 
      type(vertex), pointer :: v1 , v2 
      character(len=50) :: line 
      real(dp) :: ek(dmc_band), wq(Nph), ks(3), errU, errG, vnk(3, dmc_band)
      integer :: i_ks(3)
      complex(dp) :: uq(Nph,Nph), gkq(dmc_band,dmc_band,Nph),iden(nbnd,nbnd)
      integer :: n1,n2, inqt(3), iqt(3)
      
      iden = 0.d0 
      do ib = 1,nbnd;
         iden(ib,ib) = 1.d0 
      enddo 
   
       
       !check link 
       do iv = 2,diagram%order 
           v1 => diagram%vertexList(iv)
           v2 => diagram%vertexList(v1%link(1))
           if( v2%link(3).ne.iv ) stop 'diagram err : link err 1'
           if(v1%link(2).ne.1 .and. v1%link(2).ne.maxN) then 
               v2 => diagram%vertexList(v1%link(2))
               if( v2%link(2).ne.iv ) stop 'diagram err : link err 2'
            else 
               v2 => diagram%vertexList(v1%plink)
               if( v2%plink.ne.iv ) stop 'diagram err : plink err'
           endif 
           v2 => diagram%vertexList(v1%link(3))
           if( v2%link(1).ne.iv ) stop 'diagram err : link err 3'
           if( maxval(v1%link)>diagram%order .and. maxval(v1%link).ne.maxN ) stop 'link num > diagram%order'
       end do 
       
       !check tau diagram%order 
       do iv = 1,diagram%order 
           v1 => diagram%vertexList(iv)
           v2 => diagram%vertexList(v1%link(3)) 
           if(v2%tau < v1%tau ) stop 'diagram err : vertex not time diagram%ordered'
           if(iv.eq.1) cycle
           v2 => diagram%vertexList(v1%link(1)) 
           if(v2%tau > v1%tau ) stop 'diagram err : vertex not time diagram%ordered'
       end do 
   
       !check mometum conservation 
       do iv = 2,diagram%order 
           v1 => diagram%vertexList(iv)
           i_ks = v1%i_kin+v1%i_q-v1%i_kout
           i_ks = mod(i_ks, nk_svd)
           if( maxval(abs(i_ks))>1e-8 ) then 
               write(*,*) i_ks 
               stop 'diagram err : momentum conservation breaks at vertex'
            endif 
           v2 => diagram%vertexList(v1%link(1)) 
           i_ks = v2%i_kout-v1%i_kin
           i_ks = mod(i_ks, nk_svd)
           if( maxval(abs(i_ks))>1e-8 ) stop 'diagram err : momentum conservation breaks at electron bond'
           if(v1%link(2).ne.1 .and. v1%link(2).ne.maxN) then 
               v2 => diagram%vertexList(v1%link(2)) 
               i_ks = v2%i_q+v1%i_q
               i_ks = mod(i_ks, nk_svd)
               if( maxval(abs(i_ks))>1e-8 ) stop 'diagram err : momentum conservation breaks at phonon bond'
           else 
               v2 => diagram%vertexList(v1%plink) 
               i_ks = v2%i_q+v1%i_q
               i_ks = mod(i_ks, nk_svd)
               if( maxval(abs(i_ks))>1e-8 ) stop 'diagram err : momentum conservation breaks at ext phonon bond'
           endif 
       end do 
   
   
       !check dispersion , el-ph vertex
       do iv = 1,diagram%order
            line = 'diagram err : uk inconsistent'
            v1 => diagram%vertexList(iv)
            v2 => diagram%vertexList(v1%link(3)) 
            errU = maxval(abs(v1%ukout-v2%ukin))
            if(errU>1E-10) then
               v2%ukin = v1%ukout
               v2%i_kin = v1%i_kout
               if(errU>1e-3) then    
                  write(*,'(3i5,E20.10)') iv, v1%link(3),diagram%order,errU
                  write(*,'(3E15.5)')  v1%ekout
                  write(*,'(3E15.5)')  v2%ekin
               stop 'diagram err : uk inconsistent'
               endif 
            endif 
            errU = maxval(abs(matmul(transpose(conjg(v1%ukout)),v1%ukout) - iden))
            if(errU>1E-7) then 
               
               stop 'diagram err : uk not unitary!'
            endif 
       end do 

      !check gkq and gkq conj
      do iv = 2,diagram%order
         v1 => diagram%vertexList(iv)
         v2 => diagram%vertexList(v1%link(3)) 
         if(v1%link(3).eq.v1%link(2) .and. v1%link(2).ne.maxN) then 
            do nu = 1,Nph 
               errU = maxval(abs(v1%gkq(:,:,nu)-transpose(conjg(v2%gkq(:,:,nu)))))
               if(errU>1e-5) then 
                  write(*,'(i4)') diagram%order
                  write(*,'(3i4)') v1%i_kin
                  write(*,'(3i4)') v1%i_kout
                  write(*,'(3i4)') v1%i_q
                  write(*,*) errU
                  if(errU>1E-7) stop 'diagram err : gkq, gk+q,-q not conjguate to each other'
               endif 
            enddo 
         endif 
      end do 

      !check gkq 
      do iv = 2, diagram%order
         v1 => diagram%vertexList(iv)
         call cal_gkq_vtex_int(v1, gkq)
         errU = maxval(abs(gkq - v1%gkq)) / maxval(abs(gkq))
         if(errU>1e-4) then
            write(*,'(i4)') diagram%order
            write(*,'(3i4)') v1%i_kin
            write(*,'(3i4)') v1%i_kout
            write(*,'(3i4)') v1%i_q
            write(*,*) maxval(abs(gkq)), maxval(abs(v1%gkq))
            write(*,*) errU 
            stop 'gkq err'
         endif 
      enddo 
      
      !for JJ, velocity
      do iv = 1,diagram%order
         v1 => diagram%vertexList(iv)
         v2 => diagram%vertexList(v1%link(3))
         !check velocity  
         call cal_ek_int(v1%i_kout, v1%ekout, v1%ukout, vnk)
         errU = maxval(abs(vnk - v1%vkout)) 
         if(errU>1e-10) then
            write(*,'(i4,A20,3E20.10)')iv,'err Vnk 1 = ',maxval(abs(vnk)), maxval(abs(v1%vkout)) ,errU 
            v1%vkout = vnk
         endif 

      end do 

   end subroutine 
   
   subroutine check_epw(seed)
      use diagMC 
      implicit none 
      real(dp) :: xkt(3), xqt(3), errg,Pq
      complex(dp) ::  gkq(dmc_band,dmc_band,Nph), gkq_2(dmc_band,dmc_band,Nph)
      integer :: i,j,itr  
      type(vertex) :: v1,v2 
      TYPE(VSL_STREAM_STATE) :: seed
   
      call start_clock('check_epw')
   
     !call CreateVertex(Nph,dmc_band,nsvd,v1); call CreateVertex(Nph,dmc_band,nsvd,v2)
     !do itr = 1,100  
     !   call sample_q_omp(seed, xqt, Pq, v1)
     !   call random_number(v1%i_kin)
     !   v1%i_q = xqt
     !   v1%i_kout = v1%i_kin + v1%i_q
     !   
     !   v2%i_kin = v1%i_kout; v2%i_q = -v1%i_q; v2%i_kout = v1%i_kin 
     !   v2%gReQ = v1%gRenQ; v2%gRenQ = v1%gReQ
     ! 

     !   call cal_gkq_vtex(v1,gkq)
     !   call cal_gkq_vtex(v2,gkq_2)

     !   do i = 1,Nph 
     !      gkq_2(:,:,i) = conjg(transpose(gkq_2(:,:,i)))
     !   end do 
     !   errg = maxval(abs(gkq - gkq_2 ))/maxval(abs(gkq))
     !   write(*,'(A15, 2E15.5)')'errg = ',errg, maxval(abs(gkq))
     !   if(errg > 1e-5) then 
     !      write(*,'(A15,2E15.5)')'maxg = ',maxval(abs(gkq))**2,maxval(abs(gkq_2))**2
     !      write(*,'(9f10.5)')v1%i_kin,v1%i_kout,v1%i_q
     !      write(*,'(9f10.5)')v2%i_kin,v2%i_kout,v2%i_q
     !      write(*,'(A15,2E15.5)')'maxg = ',maxval(abs(gkq)),maxval(abs(gkq_2))
     !      write(*,'(A15,E15.5)')'errg = ',maxval(abs(gkq - gkq_2 ))/maxval(abs(gkq))
     !      write(*,'(A15,E15.5)')'errabsg = ',maxval(abs(abs(gkq) - abs(gkq_2)) )/maxval(abs(gkq))
     !      stop 'err@check_epw()'
     !   endif 
     !enddo
     

      call stop_clock('check_epw')

   end subroutine 
   
   subroutine check_phase(diagram)
      implicit none 
      type(fynman),target :: diagram
      integer :: iv, ie, jv, itr ,final_vertex 
      real(dp) :: dtau, factor 
      complex(dp) :: gfac_exact, vn(dmc_band), right_env(dmc_band, dmc_band), left_env(dmc_band,dmc_band),A(dmc_band,dmc_band), x
      type(vertex), pointer :: vi,vj

      
      vi => diagram%vertexList(1)
      vj => diagram%vertexList(vi%link(3))
      call right_environment_matrix(diagram, 1, right_env )
      

      x = 0.d0 
      do ie = 1,dmc_band 
         x = x + right_env(ie,ie)
      enddo 
      factor = real(x)
      factor = factor / abs(factor)
      if( abs(factor-diagram%gfac)>1e-5 ) then 
         write(*,*)'sign miss match!'
         write(*,'(2f10.5)')factor, real(diagram%gfac)
         stop 
      endif 

    

   end subroutine 
   

end 


!for matrix-diagmc 
subroutine left_environment_vector(diagram, init_band, final_vertex,alpha)
   use diagMC, only : fynman, dmc_band, vertex,maxN
   use pert_const, only : dp
   implicit none 
   type(fynman), target :: diagram 
   complex(dp),intent(out) :: alpha(dmc_band)
   integer :: init_band, final_vertex
   ! 
   type(vertex), pointer :: vi,vj 
   real(dp) :: Emin, dtau 
   integer :: iv,jv,ie,itr 
   complex(dp) :: A(dmc_band,dmc_band) 
   
   !begin contract environment 
   !check matrix phase 
   alpha = 0.d0; alpha(init_band) = 1.d0 
   if(final_vertex.eq.1) return 

   vi => diagram%vertexList(1)
   vj => diagram%vertexList(vi%link(3))
   do itr = 1,diagram%order-1 
      iv = vi%link(3)
      vi => diagram%vertexList(iv)
      jv = vi%link(3)
      vj => diagram%vertexList(jv)
      
      !vertex action 
      A = vi%gkq(:,:,vi%nu) 
      A = A / maxval(abs(A))
      alpha = matmul( A, alpha)
      !alpha = alpha / vi%norm_gkq
      !alpha = alpha / maxval(abs(alpha))
      if(iv.eq.final_vertex)    return 
      

      !propagator action 
      Emin = minval(vi%ekout)
      dtau = vj%tau-vi%tau
      do ie = 1,dmc_band 
         alpha(ie) = alpha(ie) * exp(-dtau*(vi%ekout(ie)-Emin))
      enddo 
      alpha = alpha / maxval(abs(alpha)) !rescale 
   enddo  
   write(*,*) 'error left_environment_vector'
   stop 
end subroutine 

subroutine left_environment_vector_hermi(diagram, init_band, final_vertex,alpha)
   use diagMC, only : fynman, dmc_band, vertex,maxN
   use pert_const, only : dp
   implicit none 
   type(fynman), target :: diagram 
   complex(dp),intent(out) :: alpha(dmc_band)
   integer :: init_band, final_vertex
   ! 
   type(vertex), pointer :: vi,vj 
   real(dp) :: Emin, dtau,Elist(dmc_band), errH 
   integer :: iv,jv,ie,itr 
   complex(dp) :: A(dmc_band,dmc_band), expH(dmc_band,dmc_band), Gt1t2(dmc_band,dmc_band) 
   
   !begin contract environment 
   vi => diagram%vertexList(1)
   alpha = 0.d0; alpha =  vi%ukout(:,init_band)
   if(final_vertex.eq.1) return 
   expH = 0.d0 
   do ie = 1,dmc_band 
      expH(ie,ie) = 1.d0 
   enddo 
   vj => diagram%vertexList(vi%link(3))
   do itr = 1,diagram%order-1 
      iv = vi%link(3)
      vi => diagram%vertexList(iv)
      jv = vi%link(3)
      vj => diagram%vertexList(jv)
      
      !vertex action 
      dtau = vj%tau -vi%tau
      Elist = vi%ekout - minval(vi%ekout)
      call propagator(vi%ukout, Elist, dtau, dmc_band, Gt1t2 )
      errH = maxval(abs(Gt1t2 - transpose(conjg(Gt1t2))))/ maxval(abs(Gt1t2))
      if(errH>1e-8) then 
         write(*,*) 'G hermicity wrong!'
         stop 
      endif 
      
      !re-hermi 
      expH = matmul(Gt1t2,expH) 
      errH = maxval(abs(expH - transpose(conjg(expH))))/ maxval(abs(expH))
      if(errH>1e-8) then 
         write(*,*) 'expH hermicity wrong!',errH
         write(*,'(A20)')'real, imag = '
         write(*,'(3f10.5)') real(expH)
         write(*,'(3f10.5)') imag(expH)
         stop 
      endif
      expH = expH / maxval(abs(expH))
      if(jv.eq.final_vertex) then 
         !check hermicity 
         write(*,'(A20)')'real, imag = '
         write(*,'(3f10.5)') real(expH)
         write(*,'(3f10.5)') imag(expH)
         return 
      endif 
   enddo  
   write(*,*) 'error left_environment_vector'
   stop 
end subroutine 


subroutine right_environment_vector(diagram, init_band, final_vertex,alpha)
   use diagMC, only : fynman, dmc_band, vertex,maxN
   use pert_const, only : dp
   implicit none 
   type(fynman), target :: diagram 
   complex(dp),intent(out) :: alpha(dmc_band)
   integer :: init_band, final_vertex
   ! 
   type(vertex), pointer :: vi,vj 
   real(dp) :: Emin, dtau 
   integer :: iv,jv,ie,itr
   complex(dp) :: A(dmc_band,dmc_band) 
   
   
   !begin contract environment 
   !check matrix phase 
   alpha = 0.d0; alpha(init_band) = 1.d0 
   if(final_vertex.eq.maxN) return 
   if(final_vertex.eq.1) stop 'head can not be on the right' 
   iv = maxN
   vi => diagram%vertexList(iv)
   jv = vi%link(1)
   vj => diagram%vertexList(jv)
   do itr = 1,diagram%order-1 
      iv = vi%link(1)
      vi => diagram%vertexList(iv)
      jv = vi%link(1)
      vj => diagram%vertexList(jv)
      
      !vertex action 
      A = transpose(vi%gkq(:,:,vi%nu)) 
      A = A / maxval(abs(A))
      alpha = matmul( A,alpha)
      alpha = alpha / maxval(abs(alpha))
      !alpha = alpha / vi%norm_gkq
      if(iv.eq.final_vertex) then 
         !alpha = alpha / maxval(abs(alpha)) !rescale 
         return 
      endif 

      !propagator action 
      Emin = minval(vi%ekin)
      dtau = vi%tau-vj%tau
      do ie = 1,dmc_band 
         alpha(ie) = alpha(ie) * exp(-dtau*(vi%ekin(ie)-Emin))
      enddo 
      !alpha = alpha / maxval(abs(alpha)) !rescale 
   enddo  
   write(*,*) 'error right_environment_vector'
   stop 
end subroutine 

subroutine right_environment_vector_hermi(diagram, init_band, final_vertex,alpha)
   use diagMC, only : fynman, dmc_band, vertex,maxN
   use pert_const, only : dp
   implicit none 
   type(fynman), target :: diagram 
   complex(dp),intent(out) :: alpha(dmc_band)
   integer :: init_band, final_vertex
   ! 
   type(vertex), pointer :: vi,vj 
   real(dp) :: Emin, dtau 
   integer :: iv,jv,ie,itr
   complex(dp) :: A(dmc_band,dmc_band),  expH(dmc_band,dmc_band)  
   
   
   !begin contract environment 
   !check matrix phase 
   alpha = 0.d0; alpha(init_band) = 1.d0 
   if(final_vertex.eq.maxN) return 
   if(final_vertex.eq.1) stop 'head can not be on the right' 
   iv = maxN
   vi => diagram%vertexList(iv)
   jv = vi%link(1)
   vj => diagram%vertexList(jv)
   do itr = 1,diagram%order-1 
      iv = vi%link(1)
      vi => diagram%vertexList(iv)
      jv = vi%link(1)
      vj => diagram%vertexList(jv)
      
      !vertex action 
      call propagator(vi%ukout, vi%ekout, vj%tau-vi%tau,dmc_band,expH )
      alpha = matmul(transpose(expH), alpha)
      alpha = alpha / maxval(abs(alpha))
      if(jv.eq.final_vertex) then 
         return 
      endif 
   enddo  
   write(*,*) 'error right_environment_vector'
   stop 
end subroutine 

subroutine left_environment_matrix(diagram, final_vertex,alpha)
   use diagMC, only : fynman, dmc_band, vertex,maxN,Gel
   use pert_const, only : dp
   implicit none 
   type(fynman), target :: diagram 
   complex(dp),intent(out) :: alpha(dmc_band,dmc_band)
   integer :: init_band, final_vertex
   ! 
   type(vertex), pointer :: vi,vj 
   real(dp) :: Emin, dtau,Elist(dmc_band) 
   integer :: iv,jv,ie,itr 
   complex(dp) :: A(dmc_band,dmc_band), expH(dmc_band,dmc_band) 
   
   !begin contract environment 
   !check matrix phase 
   alpha = 0.d0; 
   do ie = 1,dmc_band 
      alpha(ie,ie) = 1.d0 
   enddo 
   if(final_vertex.eq.1) return 

   vi => diagram%vertexList(1)
   jv = vi%link(3)
   vj => diagram%vertexList(jv)
   do itr = 1,diagram%order
      !propagator action 
      Elist = vi%ekout - minval(vi%ekout)
      dtau = vj%tau-vi%tau
      call propagator(vi%ukout,Elist, dtau, dmc_band, expH )
      alpha = matmul( expH ,alpha)
      alpha = matmul(vj%gkq(:,:,vj%nu), alpha) 
      alpha = alpha / maxval(abs(alpha))
      if(jv.eq.final_vertex) then 
         return 
      endif 

      iv = vi%link(3)
      vi => diagram%vertexList(iv)
      jv = vi%link(3)
      vj => diagram%vertexList(jv)
      !alpha = alpha / maxval(abs(alpha)) !rescale 
   enddo  
end subroutine 

subroutine right_environment_matrix(diagram, final_vertex,alpha)
   use diagMC, only : fynman, dmc_band, vertex,maxN
   use pert_const, only : dp
   implicit none 
   type(fynman), target :: diagram 
   complex(dp),intent(out) :: alpha(dmc_band,dmc_band)
   integer :: final_vertex
   ! 
   type(vertex), pointer :: vi,vj 
   real(dp) :: Emin, dtau,Elist(dmc_band) 
   integer :: iv,jv,ie,itr
   complex(dp) :: A(dmc_band,dmc_band), expH(dmc_band,dmc_band) 
   
   
   !begin contract environment 
   !check matrix phase 
   alpha = 0.d0 
   do ie = 1,dmc_band 
      alpha(ie,ie) = 1.d0 
   enddo 

   if(final_vertex.eq.maxN) return 
   !if(final_vertex.eq.1) stop 'head can not be on the right' 
   iv = maxN
   vi => diagram%vertexList(iv)
   jv = vi%link(1)
   vj => diagram%vertexList(jv)
   do itr = 1,diagram%order
      Elist = vi%ekin - minval(vi%ekin)
      dtau = vi%tau-vj%tau
      call propagator(vi%ukin, Elist, dtau, dmc_band, expH)
      alpha = matmul( alpha, expH)
      alpha = matmul(alpha, vj%gkq(:,:,vj%nu))
      alpha = alpha / maxval(abs(alpha))
      if(jv.eq.final_vertex) then 
         return 
      endif 
      iv = vi%link(1)
      vi => diagram%vertexList(iv)
      jv = vi%link(1)
      vj => diagram%vertexList(jv) 
   enddo  
end subroutine 


subroutine mid_environment_matrix(diagram, init_vertex, final_vertex, alpha)
   use diagMC, only : Gel,fynman, dmc_band, vertex,maxN
   use pert_const, only : dp
   implicit none 
   type(fynman), target :: diagram 
   complex(dp),intent(out) :: alpha(dmc_band,dmc_band)
   integer :: init_vertex,final_vertex
   ! 
   type(vertex), pointer :: vi,vj 
   real(dp) :: Emin, dtau,Elist(dmc_band) 
   integer :: iv,jv,ie,itr
   complex(dp) :: A(dmc_band,dmc_band), expH(dmc_band,dmc_band) 
   
   
   !begin contract environment 
   alpha = 0.d0 
   do ie = 1,dmc_band 
      alpha(ie,ie) = 1.d0 
   enddo 
   if(init_vertex.eq.final_vertex) return 
   if(final_vertex.eq.maxN) stop  'final = maxN @ mid_environment_matrix' 
   if(init_vertex.eq.1) stop 'init = 1 @ mid_environment_matrix' 
   
   vi => diagram%vertexList(init_vertex)
   jv = vi%link(3)
   vj => diagram%vertexList(jv)
   do itr = 1,diagram%order-1 
      !propagator action 

      Elist = vi%ekout-minval(vi%ekout)
      dtau = vj%tau-vi%tau
      call propagator(vi%ukout, Elist, dtau, dmc_band,expH)
      alpha = matmul(expH,alpha)
      if(jv.eq.final_vertex) return 
      
      alpha = matmul( vj%gkq(:,:,vj%nu),alpha )
      alpha = alpha / maxval(abs(alpha))
      iv = vi%link(3)
      vi => diagram%vertexList(iv)
      jv = vi%link(3)
      vj => diagram%vertexList(jv)
      !alpha = alpha / maxval(abs(alpha)) !rescale 
   enddo   
   write(*,*) 'error mid_environment_matrix'
   write(*,'(4i5)') iv,jv,final_vertex,diagram%order
   stop 
end subroutine 


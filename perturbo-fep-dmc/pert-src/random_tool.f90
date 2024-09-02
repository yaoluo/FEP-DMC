include 'mkl_vsl.f90'

module random_tool
   use mkl_vsl_type
   use mkl_vsl
   !use special_func
   use pert_const, only: dp,twopi,czero,pi
   implicit none 
   integer, private, parameter :: Nran = 1000
   integer, private :: iran = 1
   real*8, private :: rans(Nran) 

   interface random_number_omp
      module procedure random_number_0d_omp
      module procedure random_number_nd_omp
   end interface

   contains 

   subroutine reload_random()
       implicit none 
       iran = 1
       call random_number(rans)
   end subroutine 

   subroutine init_rand_seed(id)
       implicit none
       integer :: id 
       integer :: n,mytime(8)
       integer, dimension(:),allocatable :: seed
       call random_seed(size=n)
       allocate(seed(n))
       call date_and_time(VALUES=mytime)
       seed = 1000*(mytime(7)*1000+mytime(8)) + id*333333!+ 37 *(/(i-1,i=1,n)/)
       call random_seed(put=seed)
       !write(*,*)seed
       deallocate(seed)
   end subroutine

   subroutine uniform_int(a,b,n,pn)
       implicit none 
       integer, intent(in) ::a,b
       integer, intent(out) :: n  
       real(dp),optional, intent(out) :: pn 
       real(dp) ::  r
       call random_number(r) 
       n = int( r*(b-a+1) ) + a
       if(n>b) n = b
       if(.not.present(pn)) return 
       pn = 1.d0/(b-a+1)
   end subroutine uniform_int

   subroutine uniform_real(a,b,x,px,sample)
       implicit none 
       real(dp), intent(in) ::a,b
       real(dp), intent(out) :: px 
       real(dp) :: x
       real(dp) ::  r
       integer, optional :: sample 

       if( present(sample)) then   
           call random_number(r) 
           x = a + r*(b-a)
       end if 

       px = 1.d0/abs(b-a)
   end subroutine uniform_real

   subroutine uniform_real_3d(a,b,x)
     implicit none 
     real(dp), intent(in) ::a,b 
     real(dp), intent(out) :: x(3)
     call random_number(x) 
     x = a + x*(b-a)
 end subroutine uniform_real_3d

   subroutine Discrete(P,Ntype,i)
       implicit none 
       integer :: Ntype,i
       real(dp) :: P(Ntype), r, cumulateP
       if(Ntype.eq.1) then 
           i = 1
           return 
       end if
       call random_number(r)
       cumulateP = 0.d0
       do i = 1,Ntype
           cumulateP = cumulateP + P(i)
           if (cumulateP>r ) return;
       end do
   end subroutine Discrete


   subroutine exp_sample(a,b,beta,x,px,sample)
      implicit none 
      real(dp), intent(in) :: a,b,beta
      real(dp) :: x, px 
      integer, optional :: sample
      real(dp) :: L, xdis, absBeta 
      real(dp) :: xt, ya, yb

      if(a.eq.b) stop 'error a=b @ exp_sample()'
   
      !if(beta.eq.0) stop 'error beta=0 @ exp_sample()'
      absBeta = abs(beta)
      !if(absBeta>10) absBeta = 10.0_dp 
      L = b - a
      ya = 1.d0; yb = exp(-absBeta*L)
      if( present(sample) ) then 
         call random_number(xt)
            xdis = abs(log( xt + yb*(1.d0-xt) )/absBeta)
         if(beta>0) then 
            x = xdis + a
         else 
           x = b - xdis
         end if 
      else 
         if(beta>0) then 
           xdis =  x - a
         else 
           xdis = b - x
         end if 
      end if 

      px = absBeta*exp(-absBeta*xdis)/( 1.d0 - yb )
      !write(*,'(4E15.5)') absBeta, xdis, a,b
   end subroutine 


   !--------------------------------------------------!
   !-------------- openmp random tool ----------------!
   !--------------------------------------------------!
   subroutine init_rand_seed_omp(id,seed)
      implicit none
      integer, parameter ::  M = 2**31
      integer :: id 
      integer :: n,mytime(8), seed_num, stat
      TYPE (VSL_STREAM_STATE) :: seed
      call date_and_time(VALUES=mytime)
      seed_num = 71*(mytime(7)*13+mytime(8)) + id*3499
      stat = vslnewstream(seed, VSL_BRNG_MT19937, seed_num)
  end subroutine

   subroutine random_number_nd_omp(seed,x)
      implicit none 
      integer, parameter :: A = 1103515245 , B = 12345, M = 2**31
      real(dp) :: x(:)
      TYPE (VSL_STREAM_STATE) :: seed
      integer :: N,stat
      N = size(x)
      stat = vdrnguniform( 0, seed, N, x, 0.d0, 1.d0 )
   end subroutine

   subroutine random_number_0d_omp(seed,x)
      implicit none 
      integer, parameter :: A = 1103515245 , B = 12345, M = 2**31
      real(dp) :: x,xt(1)
      TYPE (VSL_STREAM_STATE) :: seed
      integer :: i,N,stat
      stat = vdrnguniform( 0, seed, 1, xt, 0.d0, 1.d0 )
      x = xt(1)
   end subroutine

   subroutine uniform_int_omp(seed,a,b,n,pn)
      implicit none 
      TYPE (VSL_STREAM_STATE) :: seed
      integer, intent(in) ::a,b
      integer, intent(out) :: n  
      real(dp),optional, intent(out) :: pn 
      real(dp) ::  r
      call random_number_omp(seed,r) 
      n = int( r*(b-a+1) ) + a
      if(n>b) n = b
      if(.not.present(pn)) return 
      pn = 1.d0/(b-a+1)
   end subroutine 

   subroutine uniform_real_omp(seed,a,b,x,px,sample)
      implicit none 
      TYPE (VSL_STREAM_STATE) :: seed
      real(dp), intent(in) ::a,b
      real(dp), intent(out) :: px 
      real(dp) :: x
      real(dp) ::  r
      integer, optional :: sample 

      if( present(sample)) then   
          call random_number_omp(seed,r) 
          x = a + r*(b-a)
      end if 

      px = 1.d0/abs(b-a)
   end subroutine 

   subroutine Discrete_omp(seed,P,Ntype,i)
      implicit none 
      TYPE (VSL_STREAM_STATE) :: seed
      integer :: Ntype,i
      real(dp) :: P(Ntype), r, cumulateP
      if(Ntype.eq.1) then 
          i = 1
          return 
      end if
      call random_number_omp(seed,r)
      !call random_number(r)
      cumulateP = 0.d0
      do i = 1,Ntype
          cumulateP = cumulateP + P(i)
          if (cumulateP>r ) return;
      end do
  end subroutine 

   subroutine exp_sample_omp(seed,a,b,beta,x,px,sample)
      implicit none 
      TYPE (VSL_STREAM_STATE) :: seed
      real(dp), intent(in) :: a,b,beta
      real(dp) :: x, px 
      integer, optional :: sample
      real(dp) :: L, xdis, absBeta 
      real(dp) :: xt, ya, yb

      if(a.eq.b) then 
         !stop 'error a=b @ exp_sample_omp()'
         px = 100.d0 
         x = a
         return 
      endif 
      !if(beta.eq.0) stop 'error beta=0 @ exp_sample()'
      absBeta = abs(beta)
      !if(absBeta>10) absBeta = 10.0_dp 
      L = b - a
      ya = 1.d0; yb = exp(-absBeta*L)
      if( present(sample) ) then 
         call random_number_omp(seed,xt)
            xdis = abs(log( xt + yb*(1.d0-xt) )/absBeta)
         if(beta>0) then 
            x = xdis + a
         else 
           x = b - xdis
         end if 
      else 
         if(beta>0) then 
           xdis =  x - a
         else 
           xdis = b - x
         end if 
      end if 

      px = absBeta*exp(-absBeta*xdis)/( 1.d0 - yb )
      !write(*,'(4E15.5)') absBeta, xdis, a,b
   end subroutine 

   !statisticals 


end module 


module AliasMethod
   use pert_const, only : dp
   use pert_param, only : svd_dir
   use qe_mpi_mod, only : ionode,mp_bcast,ionode_id,inter_pool_comm
   implicit none 
   integer :: Npts
   real(dp),allocatable :: grid_weight(:), onsite_weight(:)
   integer,allocatable :: alias(:)
   
   contains 

      subroutine init_alias(N,weight)
         implicit none 
         integer :: N
         real(dp) :: weight(N)
         Npts = N 
         if(allocated(grid_weight)) deallocate(grid_weight,onsite_weight,alias)
         allocate(grid_weight(Npts),onsite_weight(Npts),alias(Npts))
         grid_weight = weight
         grid_weight = grid_weight / sum(grid_weight) * Npts
      end subroutine

      !alias table is constructed using python, read in alias and onsite_weight
      subroutine read_alias()
         implicit none 
         integer :: ip 
         if (ionode) then 
            open(202, file=trim(adjustl(svd_dir))//'/alias.dat',status='old')
            do ip = 1,Npts
               read(202,*) alias(ip), onsite_weight(ip), grid_weight(ip)
               !write(*,*) alias(ip), onsite_weight(ip)
            enddo 
            close(202)
         endif 
         call mp_bcast(alias , ionode_id, inter_pool_comm)
         call mp_bcast(onsite_weight , ionode_id, inter_pool_comm)
         call mp_bcast(grid_weight , ionode_id, inter_pool_comm)
         grid_weight = grid_weight / sum(grid_weight) * Npts
      end subroutine

      subroutine sample_AlisaMethod(rnd_in, ip, Pq)
         implicit none 
         real(dp), intent(in) :: rnd_in 
         integer, intent(out)  :: ip 
         real(dp),intent(out) :: Pq
         !
         real(dp) :: rnd_res 

         ip = int(rnd_in*Npts) + 1
         rnd_res = rnd_in*Npts - (ip - 1)
         if(rnd_res>onsite_weight(ip)) ip = alias(ip)
         Pq = grid_weight(ip)

         return 
      end subroutine
end module 
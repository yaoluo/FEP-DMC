

module FynDiagram
   use pert_const, only: dp,twopi,czero,cone
   use pert_param, only : nsvd 
   !storage of diagrams 
   implicit none 

   type vertex
      integer :: link(3) = (/-1,-1,-1/)    !k(e in),q(ph),k+q(e out)
      integer :: plink 
      integer :: ikin,ikout,iq    !index 
      real(dp) :: tau   = 0.d0    !time 
      real(dp) :: kin(3)  = 0.d0  !electron in momentum / annihilation 
      real(dp) :: kout(3) = 0.d0  !electron out momentum / creation 
      real(dp) :: q(3)  = 0.d0    !phonon wave vector kout = kin + q
      integer :: i_kin(3), i_kout(3), i_q(3) !el/ph momentum index in the cubic grid  

      real(dp) :: Pq = 1.d0       !the sampling probility of q 
      real(dp) :: Einmin = 1.d0, Eoutmin = 1.d0 
      integer :: n_in, n_out, nu           !the phonon/band index m for kout, n for kin 
      real(dp),allocatable :: wq(:),ekin(:),ekout(:)                   !wq(nat*3),ek(nb)
      real(dp), allocatable :: vkin(:,:), vkout(:,:)
      complex(dp),allocatable :: uq(:,:), ukin(:,:), ukout(:,:)        !uq(nat*3, nat*3), uk(nb,nb)
      complex(dp),allocatable ::  gkq(:,:,:)               !gkg(nb,nb,nat*3), g_pol(3*nat)
      complex(dp),allocatable :: gq(:,:,:)                  !(nb, nb, nmod)
      complex(dp),allocatable :: gkq_full(:,:,:)
      !complex(dp), allocatable :: uk_f(:,:,:), ukq_f(:,:,:) !(nb,nb,nsvd)
      !complex(dp), allocatable :: vq_f(:,:,:,:), vnq_f(:,:,:,:) !(nb,nb,nmod,nsvd)
   end type vertex 

   type qmcconfig
      integer :: Nmcmc, Nbin, Ntype  
      integer :: ib               !band index
      real(dp) :: P(3)              !in crystal coord
      integer :: i_P(3)
      real(dp),allocatable :: EP(:)
      real(dp) :: mu, mu_fictious, Eshift, E0, W0  !mu is for control num, mu_fictious is for numerical accuracy 
      real(dp) :: sqmass = 1.d0 !effective mass 
      real(dp) :: tauMin = 10, tauMax = 50, dtau  
      real(dp) :: PA(10)  !probalilty for the different update 
      real(dp) :: alpha 
   end type qmcconfig 

   contains 
   subroutine InitVertexList(N, nmode, nb, nb_wan, nsvd, vertexL)
      integer,intent(in) :: N, nmode,nb,nsvd, nb_wan    
      type(vertex),intent(inout) :: vertexL(N)
      !local variable 
      integer :: iv 
      do iv = 1,N 
         !write(*,*)iv,shape(vertexL(iv)%ekin)
         call CreateVertex(nmode,nb, nb_wan,nsvd,vertexL(iv))
      end do 
   end subroutine 

   subroutine CreateVertex(nmode, nb, nb_wan, nsvd, vtex)
      integer,intent(in) :: nmode,nb,nb_wan,nsvd    
      type(vertex),intent(inout) :: vtex

      
      if(allocated(vtex%ekin)) then 
         write(*,*)shape(vtex%ekin),shape(vtex%ekout),shape(vtex%wq)
         write(*,*)shape(vtex%ukin),shape(vtex%ukout),shape(vtex%uq)
         write(*,*)shape(vtex%gkq)
         stop 'vertexlist error'
      endif 
      allocate( vtex%ekin(nb),vtex%ekout(nb),vtex%vkin(3,nb),vtex%vkout(3,nb) )
      allocate( vtex%ukin(nb_wan,nb_wan),vtex%ukout(nb_wan,nb_wan) )
      allocate( vtex%wq(nmode),vtex%uq(nmode,nmode) )
      allocate( vtex%gkq(nb,nb,nmode) )
      allocate( vtex%gkq_full(nb,nb,nmode) )
      !allocate( vtex%uk_f(nb,nb,nsvd), vtex%ukq_f(nb,nb,nsvd) )
      !allocate( vtex%vq_f(nb,nb,nmode,nsvd), vtex%vnq_f(nb,nb,nmode,nsvd) )
     
   end subroutine

   subroutine CleanVertex(vtex) 
      type(vertex),intent(inout) :: vtex
      deallocate(vtex%ekin,vtex%ekout,vtex%vkin,vtex%vkout,vtex%ukin,vtex%ukout,vtex%wq,vtex%uq)
      deallocate(vtex%gq,vtex%gkq,vtex%gkq_full)
      !deallocate(vtex%uk_f,vtex%ukq_f,vtex%vq_f,vtex%vnq_f)
   end subroutine

   subroutine copy_vtex(nmod, nbnd, nb_wan,nsvd,v1,v2)
      ! v2 = v1 
      type(vertex) :: v1
      type(vertex) ::v2 
      integer :: nbnd, nb_wan, nmod, nsvd
      !topo
      v2%link = v1%link  
      v2%plink = v1%plink    
      v2%tau = v1%tau 
      !q-dependent 
      v2%q = v1%q; v2%nu = v1%nu; v2%i_q = v1%i_q; 
      v2%Pq = v1%Pq
      v2%wq(1:nmod) = v1%wq(1:nmod)       
      !kin-dependent 
      v2%kin = v1%kin; v2%i_kin = v1%i_kin
      v2%ekin(1:nbnd)        = v1%ekin(1:nbnd)
      v2%vkin(1:3,1:nbnd)    = v1%vkin(1:3,1:nbnd)
      v2%ukin(1:nb_wan,1:nb_wan) = v1%ukin(1:nb_wan,1:nb_wan)
      !kout-dependent 
      v2%kout = v1%kout; v2%i_kout = v1%i_kout 
      v2%ekout(1:nbnd)        = v1%ekout(1:nbnd)
      v2%vkout(1:3,1:nbnd)    = v1%vkout(1:3,1:nbnd)
      v2%ukout(1:nb_wan,1:nb_wan) = v1%ukout(1:nb_wan,1:nb_wan)
      !e-ph 
      v2%gkq(1:nbnd,1:nbnd,1:nmod)          = v1%gkq(1:nbnd,1:nbnd,1:nmod)  ! maybe, I only need this one. 
      v2%gkq_full(1:nbnd,1:nbnd,1:nmod)     = v1%gkq_full(1:nbnd,1:nbnd,1:nmod)

   end subroutine
end module FynDiagram 

module PolaronWF
   use pert_const, only : dp,twopi,czero,cone
   use pert_param, only : nk_svd  
   type Psi_nk
      integer :: ngrid(3)
      real(dp),allocatable :: psi(:,:,:,:) !(n,ik1,ik2,ik3)
   end type Psi_nk

   contains 
   subroutine get_grid_index(pt, ngrid, grid_index)
      real(dp) :: pt(3)
      integer :: ngrid(3)
      integer :: grid_index(3)

      pt = pt + int(1+maxval(abs(pt)))
      !shifted center 
      grid_index = mod(int(pt*nk_svd+1e-8), nk_svd) + 1

   end subroutine 

end module 

module DiagMC
   use FynDiagram
   use PolaronWF
   use random_tool
   use pert_const, only : ryd2ev, ryd2mev
   use pert_utils, only: abs2,get_exp_ikr,bose, find_free_unit, fermi
   use vector_list,only: vlist, rand_vector_list, load_vector_list
   use wigner_seitz_cell, only: ws_cell
   use pert_data,  only: mass, epwan_fid, kc_dim, qc_dim, alat,at
   use pert_param, only: hole, fqlist, fklist, prefix, band_min, band_max, phfreq_cutoff, ntemper, temper,calc_mode,&
            cauchy_scale, sampling, nsamples,FakeRatio,Holstein,w0_Holstein,alpha_Holstein,LatticeFrohlich,alpha_frohlich,& 
            tauMin,tauMax,Nbin,dim_Holstein, gauss_frolich, spinHolstein, flip_strength, zeroTMC, read_H,Nsample_MC,frohlich_only, &
            DMC_Method, sample_gt, ib_sample   
   use qe_mpi_mod, only: ionode, stdout, mp_sum, inter_pool_comm, npool, &
            mp_split_pools, distribute_points,ionode_id, mp_bcast, my_pool_id,mp_barrier,mp_global_end
   use band_structure, only: electron_wann, init_electron_wann, solve_eigenvalue_vector,solve_band_velocity
   use phonon_dispersion, only: lattice_ifc, init_lattice_ifc, solve_phonon_modes
   use elphon_coupling_matrix, only: eph_wannier, elph_mat_wann, init_elph_mat_wann, &
      eph_fourier_elph, eph_transform,eph_transform_fast,copy_exp_ikr, eph_wan_longrange,eph_fourier_el,eph_fourier_el_para
   use AliasMethod, only : init_alias, read_alias, sample_AlisaMethod

  implicit none 

   !Notes : 
   !      variable defined in module files are always shared by different omp threads
   !logical,parameter :: LatticeFrohlich = .true. 
   !real(dp) :: alpha_frohlich = 0.1d0 

   logical ::  multiband = .false.
   logical :: abinitio = .true.   

   type(qmcconfig) :: config        !configuration
   type(lattice_ifc) :: ph          !phonon dispersion 
   type(electron_wann) :: el        !el dispersion 
   type(elph_mat_wann) :: ep        !e-ph element 


   integer :: ndim = 3              !dimension of the lattice       
   integer :: Nph = 1               !# of ph modes 
   integer :: nbnd = 1              !# of bands
   integer :: dmc_band = 1          !# of bands for DMC
   integer :: nbnd_dmc = 1          !# of bands included in DMC 
   real(dp) :: kbT                  !temperature in unit of ryd 
   real(dp), dimension(3) :: direction_mob = (/1.0,0.0,0.0/)
   integer,parameter :: maxN = 1000  !length of vertex list
   integer ::  maxOrder             !max # of vertices, maxOrder = maxN - 3, two rest are for temporary storage   

   real(dp) :: n_holstein 
   !---- fynman diagram 
   type fynman 
      TYPE (VSL_STREAM_STATE) :: seed        !random num seed 
      integer :: head, tail                  !head & tail for self-energy update 
      integer :: order = 0                   !order of the present diagrams 
      complex(dp) :: gfac = 1.d0 ,  sfac = 1.d0            !D/abs(D)  
      type(vertex) :: vertexList(maxN)       !fynman diagram is stored as a list of vertex ----
      integer :: nph_ext = 0 
   end type 
   !---- observables and statistics of updates 
   type diagmc_stat
      integer*8, dimension(:),allocatable :: update, accept, order, Tstat   !Tstat = GF if gfac always = 1                    
      complex(dp), dimension(:,:),allocatable :: GreenFunc                  !(config%Nbin,order) : sum of D/abs(D)
      complex(dp), dimension(:,:,:),allocatable :: GreenFunc_matrix       !(nb,nb,config%Nbin) : sum of D/abs(D)
      real(dp), dimension(:),allocatable :: trG_abs       !(nb,nb,config%Nbin) : sum of D/abs(D)
      real(dp), dimension(:),allocatable :: wmata                      !current-current correlator and its Mata Freq Mesh. 
      real(dp), dimension(:,:,:),allocatable :: JJw, std_JJw           !(3,3,Nw)

      complex(dp), dimension(:), allocatable :: order_g, order_sample        
      real(dp), dimension(:),allocatable :: stdGF                           !standard deviation, std can be evaluated by different mc-run  
      complex(dp) :: Etrue, Ztrue, gtrue, gfac, sfac                                            
      real(dp) :: Z, meanOrder                                           !quasi-particle weight & average order  
      real(dp),allocatable :: Znph(:,:)                       !Z for multi-phonon 
      integer*8 :: Ne
      integer,ALLOCATABLE :: Pstat(:) 
      integer :: NR = 4 
      complex(dp), dimension(:,:,:,:), allocatable :: SigmaR
      real(dp), dimension(:,:),allocatable :: Rmesh
      complex(dp), dimension(:),allocatable :: cwork
      integer*8 :: Njjw = 0 
      !polaron wfn 
      type(Psi_nk) ::  PolaronWF       !polaron wfn-electron part  
      type(Psi_nk) ::  PolaronLat      !polaron wfn-phonon part 
   end type 

   !do the memory allocation by hand to avoid memory leaking 
   type diagmc_workspace 
      complex(dp), allocatable :: expH(:,:,:), gkq(:,:,:),left_env(:,:,:),right_env(:,:,:), HexpH(:,:,:), Ak_mat(:,:,:,:),unk(:,:,:) 
      real(dp), allocatable :: vk_xyz(:,:,:) !(3,nbnd,Nmax)
   end type 

   !---- debug setup 
   integer :: debug = 1
   integer :: se_check = 1, nq_se_check = 10000
   contains 

   !set up global variable for diagmc, shared by omp thread  
   subroutine setup_dqmc()
      use pert_utils, only: abs2,get_exp_ikr,bose, find_free_unit
      use pert_param, only : Nsvd,read_svd,read_formf,svd_dir
      use eph_svd, only : init_eph_svd
      use HamiltonianTable, only : tabulate_for_DMC, grid_weight, Nktot
      implicit none
       
      real(dp) :: tau
      integer :: iunit 
      character(len=100) :: line
      
      abinitio = .false.
      if(.not. Holstein .and. .not. LatticeFrohlich .and. .not. spinHolstein) then 
         abinitio = .true.
      endif
      !write(*,*) 111 
      kbT = temper(1)*ryd2ev
      config%ib = 1
      if(ionode) then 
         iunit = find_free_unit()
         open(iunit,file='diagMC.in')
         read(iunit,*) line  
         read(iunit,*) config%Nmcmc, config%i_P, config%mu, se_check, nq_se_check, maxOrder
         read(iunit,*) line 
         read(iunit,*) config%Ntype
         read(iunit,*) line 
         read(iunit,*) config%PA(1:config%Ntype)
         close(iunit)
         config%P = config%P 
         config%PA = config%PA / sum(config%PA) !normlize 
         config%Nmcmc = config%Nmcmc * 10000
         
      end if 
     
      config%tauMin = tauMin; config%tauMax = tauMax; config%Nbin = Nbin
      if( abinitio ) config%tauMax = 1.d0/kbT
      tau = 0.1*(config%tauMax - config%tauMin) + config%tauMin
      call mp_bcast(config%Nmcmc, ionode_id, inter_pool_comm)
      call mp_bcast(config%i_P, ionode_id, inter_pool_comm)
      call mp_bcast(config%mu , ionode_id, inter_pool_comm)
      call mp_bcast(se_check , ionode_id, inter_pool_comm)
      call mp_bcast(nq_se_check , ionode_id, inter_pool_comm)
      call mp_bcast(maxOrder , ionode_id, inter_pool_comm)
      call mp_bcast(config%Nbin , ionode_id, inter_pool_comm)
      call mp_bcast(config%PA, ionode_id, inter_pool_comm)
      call mp_bcast(config%Ntype, ionode_id, inter_pool_comm)
      call mp_bcast(Tau, ionode_id, inter_pool_comm)
      
      if(maxOrder.gt.maxN) then 
         stop 'Too large maxOrder( it should <= 300)'
      endif
      !finite temeprature time 

      !statisticals of Green function and mcmc 
      config%dtau = (config%TauMax-config%TauMin)/(config%Nbin-1)

      !init DFT/DFPT input 
      call init_lattice_ifc(epwan_fid, qc_dim, ph)
      call init_electron_wann(epwan_fid, kc_dim, el)
      call init_elph_mat_wann(epwan_fid, kc_dim, qc_dim, ep)
      Nph = ep%na*3; nbnd = ep%nb 
      dmc_band = band_max - band_min + 1

      !setting for Holstein model 
      if(Holstein .or. LatticeFrohlich .or. spinHolstein) then 
         if(Holstein .or. LatticeFrohlich) then 
            Nph = 1; dmc_band = 1
         endif 
         if(spinHolstein) then 
            Nph = 1; dmc_band = 2
         endif 
         n_holstein = config%tauMax*w0_Holstein
         if(n_holstein<40) then 
            n_holstein = 1.d0/( exp(n_holstein) - 1.d0 )
         else 
            n_holstein = 0.d0 
         endif 
      end if 
      if(zeroTMC)n_holstein = 0.d0 
      
      write(*,'(A20,f10.5)') 'N_holstein = ',n_holstein

      if(Holstein.and.LatticeFrohlich) stop 'model ill defined, whether it is holstein or LatticeFrohlich'
      !---- eph-svd + Htable 
      
      if(abinitio) then 
         read_H = .true. 
         call tabulate_for_DMC(el,ph,ep)
         call init_alias(Nktot, grid_weight)
         call read_alias()
      endif

      allocate(config%Ep(dmc_band))

      call cal_ek_int(config%i_P,config%EP)
      config%EP = config%EP 
      config%PA = config%PA / sum(config%PA)
      if(ionode) then 
         write(stdout, '(A60)')'-----------------------------------------------------------------'
         write(stdout, '(A60)')'Please check the unit below! diagMC works in eV unit!'

         write(stdout,'(A20,f12.3,A5)')'Temperature  = ',kbT*1000.d0, '(meV)'
         write(stdout,'(A20,f12.5)')'alat = ',alat
         write(stdout,'(A20,3f12.5)')'A = ',at(:,1)
         write(stdout,'(A20,3f12.5)')'B = ',at(:,2)
         write(stdout,'(A20,3f12.5)')'C = ',at(:,3)
         write(stdout,'(A20,f12.3,A5)')'Ek-mu  = ',minval(config%EP)*1000.d0, '(meV)'
         write(stdout,'(A20,f12.3,A5)')'phfreq_cut = ',phfreq_cutoff*ryd2ev*1000.d0, '(meV)'
         write(stdout,'(A20,i10,A7, i10)')'Nmc, Nbin=',int(config%Nmcmc/1E4), ' x10^4, ', config%Nbin
         write(stdout,'(A20,f15.5,A8)')'beta =', config%tauMax, '(eV^-1)'
         write(stdout,'(A20,i5,4f8.3)')'ib,P,mu = ',config%ib,config%P,config%mu
         write(stdout,'(A30,3i5)')'Nph,nbnd,dmc_band = ',Nph,nbnd,dmc_band
         write(stdout, '(A15, L3,A15,i5)')'sample_gt = ',sample_gt,';  ib_sample', ib_sample
         write(stdout,'(A40)')'first-principle diagMC init finished'
         write(stdout, '(A60)')'-----------------------------------------------------------------'
         iunit = find_free_unit()
         open(iunit,file="chemical_pot.dat")
         WRITE(iunit,'(E15.5,A20)') config%mu  ,'# (eV)'
         close(iunit)
      end if 

      !stop 
   end subroutine setup_dqmc


   !init workspace 
   subroutine init_diagmc_workspace(diagram_wspace)
      implicit none 
      type(diagmc_workspace) :: diagram_wspace
      if(allocated(diagram_wspace%expH)) stop 'err init_diagmc_workspace @ diagmc_workspace allocated!'
      allocate(diagram_wspace%expH(dmc_band, dmc_band, maxN))
      allocate(diagram_wspace%gkq(dmc_band, dmc_band, maxN))
      allocate(diagram_wspace%left_env(dmc_band, dmc_band, maxN+1))
      allocate(diagram_wspace%right_env(dmc_band, dmc_band, maxN+1))
      allocate(diagram_wspace%Ak_mat(dmc_band, dmc_band, 3, maxN))
      allocate(diagram_wspace%HexpH(dmc_band, dmc_band, maxN))
      allocate(diagram_wspace%unk(dmc_band, dmc_band, maxN))
      allocate(diagram_wspace%vk_xyz(3, dmc_band, maxN))
   end subroutine

   !init observable 
   subroutine init_observable(stat)
      use HamiltonianTable , only : Nktot
      implicit none 
      integer :: iw, Npol_grid 
      type(diagmc_stat) :: stat
      if(.not.allocated(stat%update)) then 
         allocate(stat%update(7),stat.accept(7))
         allocate(stat%order(maxOrder/2+1), stat%order_g(maxOrder/2+1), stat%order_sample(maxOrder/2+1))
         allocate(stat%Tstat(config%Nbin))
         allocate(stat%stdGF(config%Nbin))
         allocate(stat%GreenFunc_matrix(dmc_band,dmc_band,config%Nbin))
         allocate(stat%trG_abs(config%Nbin))
         allocate(stat%GreenFunc(config%Nbin,maxOrder/2+1))
         allocate(stat%wmata(config%Nbin))
         allocate(stat%JJw(3,3,config%Nbin),stat%std_JJw(3,3,config%Nbin))
         allocate(stat%cwork(config%Nbin))
         allocate(stat%SigmaR(config%Nbin,dmc_band,dmc_band,stat%NR))
         allocate(stat%Rmesh(3,stat%NR))
         allocate(stat%Pstat(Nktot))
         allocate(stat%Znph(dmc_band,100))
         stat%Rmesh(:,1) = (/0,0,0/)
         stat%Rmesh(:,2) = (/1.0,0.0,0.0/)
         stat%Rmesh(:,3) = (/1.0,1.0,0.0/)
         stat%Rmesh(:,4) = (/2.0,0.0,0.0/)
      endif 
      

      stat%update = 0;    stat%accept = 0
      stat%Order = 0;     stat%order_g = 0.d0 ; stat%order_sample = 0.d0 
      stat%Tstat = 0
      stat%stdGF = 0 
      stat%Etrue = 0;stat%Ztrue = 0;stat%gtrue = 0.d0
      stat%Ne = 0; stat%gfac = 0.d0 ; stat%sfac = 0.d0 
      stat%JJw = 0; stat%std_JJw = 0 
      stat%Znph = 0.d0 

      if(multiband) then 
         stat%GreenFunc_matrix = 0;
         stat%trG_abs = 0;
      else 
         stat%GreenFunc = 0;
      endif
      stat%Pstat = 0
      stat%SigmaR = 0
      do iw = 1,config%Nbin 
         stat%wmata(iw) = ( iw - 1 ) / config%tauMax * twopi 
      enddo 
      stat%Njjw = 0 

      !initialize polaron wfn 
      Npol_grid = nk_svd 
      stat%PolaronWF%ngrid = Npol_grid
      allocate( stat%PolaronWF%psi(dmc_band,Npol_grid,Npol_grid,Npol_grid) )
      stat%PolaronWF%psi = 0.d0
      stat%PolaronLat%ngrid =  Npol_grid
      allocate( stat%PolaronLat%psi(nph,Npol_grid,Npol_grid,Npol_grid) )
      stat%PolaronLat%psi = 0.d0 
   end subroutine

   !init diagram for green function 
   subroutine init_diagram(diagram, tau)
      implicit none 
      type(fynman), target :: diagram
      type(vertex), pointer :: vi
      real(dp), optional :: tau 
      real(dp) :: tau_, iden(dmc_band,dmc_band)    
      integer :: ib, nu 

  
      tau_ = 1.0    
      if(present(tau)) tau_ = tau 
      if(tau_<config%tauMin .or. tau_>config%tauMax) tau_ =0.5*(config%tauMin + config%tauMax)

      iden = 0.d0 
      do ib = 1,dmc_band 
         iden(ib,ib) = 1.d0 
      enddo 
      !init diagram 
      diagram%order = 1; 
      diagram%gfac = 1.d0; diagram%sfac = 1.d0; 
      vi => diagram%vertexList(1)
      vi%i_kout = config%i_P
      vi%ekout = config%EP; vi%tau = 0.d0
      vi%link(3) = maxN
      call cal_ek_int(vi%i_kout,vi%ekout,vi%ukout,vi%vkout); 
      vi%n_out = config%ib 
      vi%nu = 1
      do nu = 1,Nph 
         vi%gkq(:,:,nu) = iden 
         vi%gkq_full(:,:,nu) = iden 
      enddo 

      vi => diagram%vertexList(maxN)
      vi%i_kin = config%i_P;  vi%tau = tau_
      vi%link(1) = 1
      call cal_ek_int(vi%i_kin,vi%ekin,vi%ukin,vi%vkin) 
      vi%n_in = config%ib 
      vi%nu = 1
      do nu = 1,Nph 
         vi%gkq(:,:,nu) = iden
         vi%gkq_full(:,:,nu) = iden 
      enddo 
   end subroutine 

   !init diagram for current correlator
   subroutine init_diagram_JJ(diagram, tau)
      implicit none 
      type(fynman), target :: diagram
      type(vertex), pointer :: vi
      real(dp), optional :: tau 
      real(dp) :: tau_ 

      tau_ = 1.0    
   
      !init diagram 
      diagram%order = 1; 
      diagram%gfac = 1.d0; diagram%sfac = 1.d0; 
      vi => diagram%vertexList(1)
      vi%i_kout = 0; vi%tau = 0.d0
      vi%link(3) = maxN
      call cal_ek_int(vi%i_kout,vi%ekout,vi%ukout,vi%vkout)

      vi => diagram%vertexList(maxN)
      vi%i_kin = 0; vi%tau = config%tauMax 
      vi%link(1) = 1
      call cal_ek_int(vi%i_kin,vi%ekin,vi%ukin,vi%vkin) 
   end subroutine 

   !phonon propagator 
   function Dph(w,tau) 
      implicit none 
      real(dp) :: Dph, D1, D2  
      real(dp),INTENT(IN) :: w,tau 
      real(dp) :: nq,tau_abs
      
      tau_abs = abs(tau)
      if(w<phfreq_cutoff*ryd2ev) then 
         Dph = 1e-15_dp
         return  
      end if 
      if(zeroTMC) then 
         Dph =  exp(-tau_abs*w)
         return 
      endif 
     
      !numerical stable form 
      D1 = exp(-min(tau_abs*w, 100.0_dp)) / (1.0_dp - exp( -min(config%tauMax*w, 100.0_dp) ) )
      D2 = exp(-min((config%tauMax - tau_abs)*w, 100.0_dp)) / (1.0_dp - exp( -min(config%tauMax*w, 100.0_dp) ) )
      Dph = D1 + D2
   end function

   !electron propagtor 
   function Gel(ek,tau)  
      implicit none
      real(dp) :: Gel  
      real(dp),INTENT(IN) :: ek,tau
      real(dp) :: fk 
      !fk = fermi( kbT, ek )
      Gel = exp(-ek*tau) 
   end function

   subroutine sample_q_omp(seed, xqt, Pq, vtex)
      use HamiltonianTable, only : random_qfromlist_omp
      implicit none 
      TYPE(VSL_STREAM_STATE) :: seed
      real(dp) :: xqt(3), rnd  
      real(dp), optional :: Pq
      type(vertex), optional :: vtex
      !
      integer :: iq 
      
      if(Holstein) then 
         call random_number_omp(seed,xqt)
         xqt = xqt !- 0.5  
      endif 
      
   end subroutine

   subroutine sample_q_omp_int(seed, ixqt, Pq, vtex)
      use HamiltonianTable, only : random_qfromlist_omp, ikgrid, grid_weight
      implicit none 
      TYPE(VSL_STREAM_STATE) :: seed
      real(dp) ::  rnd  
      integer :: ixqt(3)

      real(dp), optional :: Pq
      type(vertex), optional :: vtex
      !
      integer :: iq 
      
      if(present(vtex)) then  
         call random_number_omp(seed,rnd)
         call sample_AlisaMethod(rnd, iq, Pq)
         ixqt = ikgrid(:,iq)
         return 
      endif 

   end subroutine

   !----------------------------------------------------------------------------
   !----------    driver for abinit e-ph Hamiltonian ---------------------------
   !----------------------------------------------------------------------------
   subroutine cal_ek(xpt,ek,uk,vk)
      use HamiltonianTable, only : solve_electron_fast
      use pert_param, only : exact_electron 
      implicit none 
      real(dp), intent(in):: xpt(3)
      real(dp), intent(out) :: ek(dmc_band)
      real(dp) :: xcop(3)
      complex(dp),optional,intent(out) :: uk(nbnd,nbnd)
      real(dp), optional, intent(out) :: vk(3,dmc_band)
      real(dp) :: vkt(3,nbnd)
      complex(dp) :: ukt(nbnd,nbnd)
      !write(*,'(A20,i10)') 'size ek = ',size(ek)

      
      if(spinHolstein) then 
         ek= -2.d0*sum(cos(twopi*xpt(1:dim_Holstein))) - config%mu !1D Holstein 
         if( present(uk) ) then 
            uk(1,1)=1.d0 ; uk(2,2)=1.d0 
         endif 
         return 
      endif

      if(LatticeFrohlich) then 
         ek(1) = 0.5d0 * sum(xpt**2) - config%mu !3D Frolich 
         if( present(uk) ) uk = 1.d0 
         return 
      endif

      if(Holstein) then 
         ek(1) = -2.d0*sum(cos(twopi*xpt(1:dim_Holstein))) - config%mu !1D Holstein 
         
         if( present(uk) ) uk=1.d0 
         if( present(vk)) then 
            vk = 0 
            vk(1,1) =  2.d0*sin(twopi*xpt(1))
         endif 
         return 
      endif

      !first principle electron

      xcop = xpt 
      if(hole) xcop = -xcop 
      !call To1BZ(xcop)
      if( present(uk) )  then 
         if(present(vk)) then 
            call solve_band_velocity(el, xcop, vkt, ek, uk)
            ek = ek * ryd2ev  
            vk = vk * ryd2ev * alat
            if(hole) vk = -vk 
         else 
            call solve_eigenvalue_vector(el, xcop, ek, uk); 
            ek = ek * ryd2ev 
         endif 
         
      else 
         call solve_eigenvalue_vector(el, xcop, ek); ek = ek * ryd2ev 
      endif 
      ek = ek - config%mu

      if(hole) ek = -ek  
   end subroutine 
   
   subroutine cal_wq(xqt,wq,uq)

      implicit none 
      real(dp), intent(in):: xqt(3)
      real(dp), intent(out) :: wq(Nph) 
      complex(dp), intent(out) :: uq(Nph,Nph)
      !call solve_phonon_modes(ph, xqt, wq, uq)
      if(Holstein .or. spinHolstein) then 
         uq = 1.d0; wq = w0_Holstein; 
         return 
      endif
      if(LatticeFrohlich) then 
         uq = 1.d0; wq = 1.0d0; 
         return 
      endif
      call solve_phonon_modes(ph, xqt, wq,  uq)
   end subroutine 
   
   !bare vertex at k  
   subroutine cal_vk(xpt,vk)
      !v_kx = d E/ dk_x 
      !use HamiltonianTable, only 
      implicit none 
      real(dp), intent(in):: xpt(3)
      real(dp), intent(out) :: vk(3,nbnd)
      if(Holstein) then 
         vk(1,1) = 2.d0*sin( twopi*xpt(1) ) !bare (velocity) vertex at k  
         return 
      endif
      if(spinHolstein) then 
         vk(1,1) = 1; vk(1,2) = -1 !bare (spin) vertex at k  
         return 
      endif

   end subroutine 

   !g_kq_sym = 0.5*(g_kq + conjg(g_{k+q,-q})
   !ensure g_kq_sym =conjg( g_{k+q -q}_sym )
   subroutine cal_gkq_symmetrize(kpt,xqt,wq,uq,uk,ukq,g_kq, g_pol, vformf, vformfc)
      use HamiltonianTable, only : solve_gkq_fast
      use pert_param, only :exact_gkq, drop_frohlich
      implicit none 
      real(dp),intent(in) :: kpt(3), xqt(3),wq(Nph)                           
      complex(dp),intent(in) :: uq(Nph,Nph), uk(nbnd,nbnd), ukq(nbnd,nbnd) !eigenvector for el and ph 
      complex(dp),intent(out) ::  g_kq(nbnd, nbnd, Nph)
      complex(dp),intent(in), dimension(nbnd,nbnd,Nph,nsvd),optional :: vformf, vformfc

      complex(dp),optional :: g_pol(Nph) 
      complex(dp) :: ephlr(Nph), iden(nbnd,nbnd),ctmp(Nph)
      integer :: im, ib, jb 
      real(dp) :: q_1bz(3), g1,g2 

      g_kq = 0.d0
      if(spinHolstein) then 
         g1 = sqrt(dim_Holstein * w0_Holstein * alpha_Holstein ) 
         g2 = sqrt(dim_Holstein * w0_Holstein * flip_strength * alpha_Holstein ) 
         g_kq(1,1,1) = g1; g_kq(2,2,1) = g1
         g_kq(1,2,1) = g2; g_kq(2,1,1) = g2
         return 
      endif 

      if(Holstein) then
         g_kq = sqrt(dim_Holstein * w0_Holstein * alpha_Holstein ) 
         return 
      endif
      
      if(LatticeFrohlich) then 
         !q_1bz = xqt; 
         !call To1BZ(q_1bz)
         if(sum(xqt**2)>1e-9) then 
            !g_kq = sqrt(alpha_frohlich * sqrt(2.d0)*twopi/ sum(xqt**2)/ (twopi**3) )
            g_kq = sqrt(alpha_frohlich * sqrt(2.d0)/ sum(xqt**2)/ (twopi**2) ) ! alpha/q^2, rescale for convenience 
         endif 
         return
      endif 
   end subroutine


   subroutine cal_gkq_vtex(vtex,g_kq)
      !use eph_ReQ, only : cal_symmetric_gkq_from_gReQ_gRenQ, cal_symmetric_gkq_from_gReQ_gRenQ_hole_version
      implicit none 
      type(vertex) :: vtex
      complex(dp),intent(out) ::  g_kq(nbnd, nbnd, Nph)
      !work 
      integer :: im 
      if(abinitio) then 
         if(hole) then 
            !call cal_symmetric_gkq_from_gReQ_gRenQ_hole_version(vtex%kin,vtex%q,vtex%gReQ, vtex%gRenQ, g_kq )
         else 
            !call cal_symmetric_gkq_from_gReQ_gRenQ(vtex%kin,vtex%q,vtex%gReQ, vtex%gRenQ, g_kq )     
         endif      
      else 
         !call cal_gkq_symmetrize(vtex%kin,vtex%q,vtex%wq,vtex%uq,vtex%ukin,vtex%ukout, g_kq) 
      endif 

      g_kq = g_kq * sqrt(alpha_frohlich)

   end subroutine

   !ixpt 
   subroutine cal_ek_int(ixpt,ek,uk,vk)
      use HamiltonianTable, only : solve_electron_fast
      use pert_param, only : exact_electron 
      implicit none 
      integer, intent(inout):: ixpt(3)
      real(dp), intent(out) :: ek(dmc_band)
      integer :: ixcop(3), ib 
      complex(dp),optional,intent(out) :: uk(nbnd,nbnd)
      real(dp), optional, intent(out) :: vk(3,dmc_band)
      real(dp) :: vkt(3,nbnd)
      real(dp) :: ekt(nbnd)
      !write(*,'(A20,i10)') 'size ek = ',size(ek)

      ixcop = ixpt 
      if(hole) ixcop = -ixpt 
      if( present(uk) )  then 
         if(present(vk)) then 
            call solve_electron_fast(ixcop, ekt, uk=uk, vk=vkt) !uk=, vk= for optional variable 
            vkt = vkt * ryd2ev * alat
            if(hole) vkt = -vkt 
            do ib = 1,dmc_band
               vk(:, ib) = vkt(:, ib-1+band_min)
            enddo 
         else 
            call solve_electron_fast(ixcop, ekt, uk=uk)  
         endif 
      else 
         call solve_electron_fast(ixcop, ekt)
      endif 

      ekt = ekt * ryd2ev 
      ekt = ekt - config%mu
      if(hole) ekt = -ekt  
      do ib = 1,dmc_band
         ek(ib) = ekt(ib-1+band_min)
      enddo 
   end subroutine 

   subroutine cal_wq_int(ixqt,wq)
      use HamiltonianTable, only : solve_phonon_fast
      implicit none 
      integer, intent(inout):: ixqt(3)
      real(dp), intent(out) :: wq(Nph) 

      if(Holstein .or. spinHolstein) then 
        wq = w0_Holstein; 
         return 
      endif
      if(LatticeFrohlich) then 
         wq = 1.0d0; 
         return 
      endif
      call solve_phonon_fast(ixqt, wq)

      wq = wq * ryd2ev 
   end subroutine 
   
   subroutine cal_gkq_vtex_int(vtex,g_kq)
      use HamiltonianTable, only : solve_gkq_fast, solve_gkq_full_fast
      implicit none 
      type(vertex) :: vtex
      complex(dp),intent(out) ::  g_kq(dmc_band, dmc_band, Nph)
      !work 
      integer :: im, ikt(3), iqt(3) 
      complex(dp) ::  gtmp(nbnd, nbnd, Nph), ukqc(dmc_band, nbnd), uk(nbnd, dmc_band)

      !put i_kin, i_kout, i_q
      vtex%i_kin = mod(vtex%i_kin, nk_svd)
      vtex%i_kout = mod(vtex%i_kout, nk_svd)
      if(abinitio) then    
         if(hole) then 
            !g_{nm}(-k-q,q) = g^*_{mn}(-k,-q)
            ikt = -vtex%i_kin 
            iqt = -vtex%i_q
            if(DMC_Method.eq.0) then 
               call solve_gkq_full_fast( ikt, iqt, gtmp ) 
            else 
               call solve_gkq_fast( ikt, iqt, gtmp ) 
            endif 
            gtmp = dconjg(gtmp)
         else 
            ikt = vtex%i_kin 
            iqt = vtex%i_q
            if(DMC_Method.eq.0) then 
               call solve_gkq_full_fast( ikt, iqt, gtmp ) 
            else 
               call solve_gkq_fast( ikt, iqt, gtmp ) 
            endif 
         endif    
         
      else 
         call cal_gkq_symmetrize(vtex%kin,vtex%q,vtex%wq,vtex%uq,vtex%ukin,vtex%ukout, gtmp) 
      endif 
      gtmp = gtmp * ryd2ev
      gtmp = gtmp * sqrt(alpha_frohlich)

      ! put to bloch 
      ukqc = transpose(conjg(vtex%ukout(:,band_min:band_max)))
      uk = vtex%ukin(:,band_min:band_max)
      do im = 1,Nph 
         g_kq(:,:,im) = matmul( ukqc, matmul(gtmp(:,:,im), uk) )
      enddo 
   end subroutine

   subroutine cal_gkq_full_vtex_int(vtex,g_kq)
      use HamiltonianTable, only : solve_gkq_full_fast
      implicit none 
      type(vertex) :: vtex
      complex(dp),intent(out) ::  g_kq(dmc_band, dmc_band, Nph)
      !work 
      integer :: im, ikt(3), iqt(3) 
      complex(dp) ::  gtmp(nbnd, nbnd, Nph), ukqc(dmc_band, nbnd), uk(nbnd, dmc_band)

      !put i_kin, i_kout, i_q
      vtex%i_kin = mod(vtex%i_kin, nk_svd)
      vtex%i_kout = mod(vtex%i_kout, nk_svd)
      if(abinitio) then    
         if(hole) then 
            !g_{nm}(-k-q,q) = g^*_{mn}(-k,-q)
            ikt = -vtex%i_kin 
            iqt = -vtex%i_q
            call solve_gkq_full_fast( ikt, iqt, gtmp )    
            gtmp = dconjg(gtmp)
         else 
            ikt = vtex%i_kin 
            iqt = vtex%i_q
            call solve_gkq_full_fast( ikt, iqt, gtmp )    
         endif    
         
      else 
         call cal_gkq_symmetrize(vtex%kin,vtex%q,vtex%wq,vtex%uq,vtex%ukin,vtex%ukout, gtmp) 
      endif 
      gtmp = gtmp * ryd2ev
      gtmp = gtmp * sqrt(alpha_frohlich)


      ! put to bloch 
      ukqc = transpose(conjg(vtex%ukout(:,band_min:band_max)))
      uk = vtex%ukin(:,band_min:band_max)
      do im = 1,Nph 
         g_kq(:,:,im) = matmul( ukqc, matmul(gtmp(:,:,im), uk) )
      enddo 
   end subroutine

   !----------------------------------------------------------------------------
   !---------------------- debug -----------------------------------------------
   !---------------------------------------------------------------------------- 
   subroutine print_vertex(diagram)
      implicit none 
      integer :: iv 
      type(fynman), target :: diagram
      type(vertex), pointer :: vi
      do iv = 1,diagram%order 
         vi => diagram%vertexList(iv)
         write(*,'(i5,A4,3i5,4E15.5, 3E15.5)')iv,'|', vi%link,vi%tau, vi%ekin(vi%n_in), &
                     vi%wq(vi%nu)  ,vi%ekout(vi%n_out), vi%q
         !write(*,'(i5,A4,3i5,4E15.5, 3E15.5)')iv,'|', vi%link,vi%tau, vi%ekin(vi%n_in), &
         !            vi%wq(vi%nu)  ,vi%ekout(vi%n_out), vi%q
      end do 
      vi => diagram%vertexList(maxN)
      write(*,'(i5,A4,3i5,f10.5)') maxN,'|',vi%link,vi%tau 
   end subroutine 

   !output 
   subroutine output_gt(stat,suffix_name)

      implicit none
      type(diagmc_stat) :: stat
      character(len=100), optional :: suffix_name
      integer :: iunit, it, order,ib,jb, job   
      real(dp) :: dtau, tau_ev
      real(dp) :: T_real(config%Nbin)
      INTEGER*4  :: access,status
      logical :: res
      character(len=100) :: str  

      do job = 1,Nsample_MC
         write(str,'(i10)') job 
         !status = access( 'green_Matrix.dat-'//trim(adjustL(str)), ' ')    ! blank mode
         inquire( file='green_Matrix.dat-'//trim(adjustL(str)), exist=res )
         if ( .not. res ) then 
            exit 
         endif 
      enddo 
      T_real = stat%Tstat
      T_real = T_real / sum( T_real/config%Nbin )


      dtau = (config%TauMax-config%TauMin)/config%Nbin
     
      iunit = find_free_unit() 
      open(iunit,file='green_Matrix.dat-'//trim(adjustL(str)))
      write(iunit,'(A20,2A30)')'# tau', 'tau_stats','Real GF'
      do it = 1,config%Nbin 
         tau_ev = ((it-1)*dtau+config%tauMin)
         write(iunit,'(E20.10,E30.20)',advance='no') tau_ev, T_real(it)
         do ib = 1,dmc_band
            do jb = 1,dmc_band 
                  write(iunit,'(E30.20)',advance='no') real(stat%GreenFunc_matrix(ib,jb,it))
            enddo
         enddo
         write(iunit,'(A2)') ' '
      end do  
      close(iunit)

      open(iunit,file='trG_abs.dat-'//trim(adjustL(str)))
      write(iunit,'(A20,2A30)')'# tau', 'tau_stats','Real GF'
      do it = 1,config%Nbin 
         tau_ev = ((it-1)*dtau+config%tauMin)
         write(iunit,'(E20.10, E30.20)') tau_ev, stat%trG_abs(it)
      end do  
      close(iunit)

      

      
    
      iunit = find_free_unit()
      open(iunit,file='orderStat.dat-'//trim(adjustL(str)))
      do it = 1,maxOrder/2+1
         write(iunit,'(5E30.20)') real(stat%order(it),dp), stat%order_g(it) , stat%order_sample(it)
      enddo 
      close(iunit)    


    
   end subroutine 

   !output 
   subroutine output_EZ(stat,suffix_name)
      use HamiltonianTable, only : find_ik, enk_grid, wq_grid
      implicit none
      type(diagmc_stat) :: stat
      character(len=100), optional :: suffix_name
      integer :: iunit, it, order,ib,jb, job, i1,i2,i3,iph, ikpt(3), ik    
      real(dp) :: dtau, tau_ev, ekt(nbnd)
      INTEGER*4  :: access,status
      logical :: res
      character(len=100) :: str  

      do job = 1,Nsample_MC
         write(str,'(i10)') job 
         !status = access( 'green_Matrix.dat-'//trim(adjustL(str)), ' ')    ! blank mode
         inquire( file='Znph.dat-'//trim(adjustL(str)), exist=res )
         if ( .not. res ) then 
            exit 
         endif 
      enddo 

      iunit = find_free_unit() 
      open(iunit,file='Znph.dat-'//trim(adjustL(str)))
      write(iunit,'(A10)') '#<Z0> = '
      do it = 1,100
         do ib = 1,dmc_band
            WRITE(iunit,'(E30.10)',advance='no') stat%Znph(ib,it)
         enddo 
         WRITE(iunit,'(A2)') ' '
      enddo
      write(iunit,'(A10)') '#<E> = '
      do ib = 1,dmc_band
         WRITE(iunit,'(E30.10)',advance='no') real(stat%Etrue)
      enddo 
      WRITE(iunit,'(A2)') ' '
      do ib = 1,dmc_band
         WRITE(iunit,'(E30.10)',advance='no') imag(stat%Etrue)
      enddo 
      close(iunit)

      open(iunit,file='orderStat.dat-'//trim(adjustL(str)))
      do it = 1,maxOrder/2+1
         write(iunit,'(5E30.20)') real(stat%order(it),dp), stat%order_g(it), stat%order_sample(it)
      enddo 
      close(iunit)

      open(iunit,file='sign.dat-'//trim(adjustL(str)))
      write(iunit,'(2E30.20, i30)') stat%gtrue, stat%Njjw
      close(iunit)    

      !output wfn 
      open(iunit, file='PolaronWF.dat-'//trim(adjustL(str)))
      do i1 = 1,stat%PolaronWF%ngrid(1)
         do i2 = 1,stat%PolaronWF%ngrid(2)
            do i3 = 1,stat%PolaronWF%ngrid(3)
               do ib = 1,dmc_band 
                  write(iunit,'(E20.10)',advance='no') stat%PolaronWF%psi(ib,i1,i2,i3) 
               enddo 
               !energy 
               ikpt = (/i1-1, i2-1, i3-1/)
               ik = find_ik(ikpt)
               do ib = 1,dmc_band 
                  write(iunit,'(E20.10)',advance='no') enk_grid(ib-1+band_min,ik) * ryd2ev
               enddo 
               write(iunit,'(A2)') ' '
      enddo;enddo;enddo
      close(iunit)    

      open(iunit, file='PolaronLat.dat-'//trim(adjustL(str)))
      do i1 = 1,stat%PolaronLat%ngrid(1)
         do i2 = 1,stat%PolaronLat%ngrid(2)
            do i3 = 1,stat%PolaronLat%ngrid(3)
               do iph = 1,nph 
                  write(iunit,'(E20.10)',advance='no') stat%PolaronLat%psi(iph,i1,i2,i3) 
               enddo 
               !energy 
               ikpt = (/i1-1, i2-1, i3-1/)
               ik = find_ik(ikpt)
               do iph = 1,nph 
                  write(iunit,'(E20.10)',advance='no') wq_grid(iph,ik) * ryd2mev 
               enddo 
               write(iunit,'(A2)') ' '
      enddo;enddo;enddo
      close(iunit)    



      open(iunit,file='parameter.dat-'//trim(adjustL(str)))
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
      WRITE(iunit,'(E20.10)')real(stat%Etrue)
      write(iunit,'(A10)') '#<g/Z> = '
      WRITE(iunit,'(E20.10)')real(stat%gtrue)
      write(iunit,'(A10)') '#<Z0> = '
      WRITE(iunit,'(E20.10)') stat%Znph(:,1) / sum( stat%Znph )
      close(iunit)



   end subroutine 



   !output 
   subroutine output_SigmaR(stat,suffix_name)
      implicit none
      type(diagmc_stat) :: stat
      character(len=100), optional :: suffix_name
      integer :: iunit, it, order,ib,jb,ir   
      real(dp) :: dtau, tau_ev
      real(dp) :: T_real(config%Nbin)

      dtau = (config%TauMax-config%TauMin)/config%Nbin
      if(.not.ionode) return 

      iunit = find_free_unit()
      !sigmaR

      if( present(suffix_name) ) then 
         open(iunit,file='sigmaR.dat-'//trim(adjustL(suffix_name)))
      else 
         open(iunit,file='sigmaR.dat')
      end if 
      write(iunit,'(2A20,A50)')'# tau', 'Real sigma_R'
      do it = 1,config%Nbin 
         tau_ev = ((it-1)*dtau+config%tauMin)
         write(iunit,'(2E20.10)',advance='no') tau_ev
         do ir = 1,stat%NR
            do ib = 1,dmc_band; do jb = 1,dmc_band 
               write(iunit,'(2E30.20)',advance='no') real(stat%SigmaR(it,ib,jb,ir)),imag(stat%SigmaR(it,ib,jb,ir))
            enddo;enddo
         enddo
         write(iunit,'(A2)') ' '
      end do 
      close(iunit)  
   end subroutine 
   
   function trace(A,n)
      use pert_const, only : dp
      integer, intent(in) :: n 
      complex(dp), intent(in) :: A(n,n)
      complex(dp) :: trace 
      integer :: i 
      trace = 0.d0 
      do i = 1,n 
         trace = trace + A(i,i)
      enddo 
      return 
   end function

end module DiagMC



!measure 
subroutine measureE(diagram,stat)
   use diagmc, only : dp,diagmc_stat, fynman, maxN,config,vertex
   implicit none 
   type(diagmc_stat) :: stat
   type(fynman),target :: diagram 
   integer :: iv,m,n,nu 
   type(vertex), pointer :: vi 
   real(dp) :: Edt, tau 
   Edt = 0.d0 
   do iv = 2,diagram%order 
      vi => diagram%vertexList(iv)
      m = vi%n_out; n = vi%n_in; nu = vi%nu 
      tau = vi%tau
      Edt = Edt + (vi%ekin(n)-vi%ekout(m))*tau
      Edt = Edt + 0.5d0 * abs(tau - diagram%vertexList(vi%link(2))%tau) * vi%wq(nu) 
   enddo 
   vi => diagram%vertexList(maxN)
   Edt = Edt + vi%ekin(vi%n_in) * vi%tau
   stat%gfac = stat%gfac + diagram%gfac
   stat%Etrue = stat%Etrue + diagram%gfac * (Edt - (diagram%order-1) ) / vi%tau
   stat%Ne = stat%Ne + 1
end subroutine 


!do some stat 
subroutine post_process(stat)
   use diagmc, only : maxOrder,dp,diagmc_stat, maxN,config,vertex
   implicit none 
   type(diagmc_stat) :: stat
   integer :: i 
   real(dp) :: meanOrder, orderState(maxN)

   orderState(1:maxOrder/2) = stat%order(1:maxOrder/2)
   orderState = orderState / sum(orderState)
   meanOrder = 0
   do i = 1,maxOrder/2
      meanOrder = meanOrder + orderState(i) * i 
   enddo
   stat%meanOrder = meanOrder
end subroutine



subroutine propagator(uk,ek,tau,nb,expH)
   use pert_const, only : dp 
   use diagMC, only : nbnd 
   implicit none 
   complex(dp) :: uk(nbnd,nbnd), expH(nb,nb)
   real(dp) :: ek(nb), tau 
   integer :: nb
   !
   integer :: ib

   expH = 0.d0 
   do ib =  1,nb 
      expH(ib,ib) =  exp(-min(tau*ek(ib),100.0_dp))
   enddo 
   return 

   !A = matmul(transpose(conjg(uk)),uk)
   !err = 0 
   !do ib = 1,nb 
   !   err = err + abs(A(ib,ib)-1.d0)
   !enddo 
   !if(err>1e-10) stop 'not unitary'
   !expH = 0.5d0 * ( expH + transpose(conjg(expH)) )
end subroutine


subroutine Hpropagator(uk,ek, mu,tau,nb,HexpH)
   use pert_const, only : dp 
   use diagMC, only : nbnd 
   implicit none 
   complex(dp) :: uk(nbnd,nbnd), HexpH(nb,nb)
   real(dp) :: ek(nb), mu, tau 
   integer :: nb
   !
   integer :: ib
   HexpH = 0.d0 
   do ib =  1,nb 
      HexpH(ib,ib) = exp(-tau*ek(ib)) * (ek(ib)+mu) * tau
    enddo 
   return
end subroutine

subroutine Apropagator(uk,ek, vk, mu,tau,nb,HexpH)
   use pert_const, only : dp 
   implicit none 
   complex(dp) :: uk(nb,nb), HexpH(nb,nb)
   real(dp) :: ek(nb), vk(nb), mu, tau 
   integer :: nb
   !
   integer :: ib
   do ib =  1,nb 
      HexpH(:,ib) = uk(:,ib)  * exp(-tau*ek(ib)) * vk(ib) * tau
    enddo 
   HexpH = matmul(HexpH,transpose(conjg(uk)))
   !expH = 0.5d0 * ( expH + transpose(conjg(expH)) )
end subroutine



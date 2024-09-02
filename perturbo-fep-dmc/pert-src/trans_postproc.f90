!===============================================================================
! Copyright (C) 2016-2020 Jin-Jian Zhou, Jinsoo Park, I-Te Lu, Marco Bernardi
!
! This program is distributed under the terms of the GNU General Public License.
! See the file `LICENSE' in the root directory of this distribution, or obtain 
! a copy of the License at <https://www.gnu.org/licenses/gpl-3.0.txt>.
!
! Author: jjzhou <jjchou.comphy@gmail.com>
! Comment:
!   compute more transport coefficients starting from the pre-calculated TDF
!
! Maintenance:
!===============================================================================

subroutine trans_postproc()
   use pert_const, only: dp, ryd2ev
   use pert_data,  only: volume
   use pert_utils, only: mfermi_deriv, fermi
   use qe_mpi_mod, only: stdout
   use pert_param, only: ntemper, temper, efermi, ftemper, prefix
   use boltz_trans_mod, only: extract_trans_coeff
   use boltz_trans_output, only: output_trans_coef
   use hdf5_utils
   implicit none
   logical :: has_file, is_hole
   integer(HID_T) :: file_id
   integer :: ntmp, nestep, it, i, itmp, ncomp, tdf_size(2)
   real(dp) :: de, tp, ef, cond_seebeck(6), kk_coef(6), alpha(6)
   character(len=120) :: dset_name, fname, msg
   integer, allocatable :: niter(:)
   real(dp), allocatable :: tmpr(:), efem(:), density(:), occup(:)
   real(dp), allocatable :: ene(:), tdf(:,:), dos(:), tdf_tmp(:,:)
   real(dp), allocatable :: cond(:,:), seebeck(:,:), therm_cond(:,:)

   if(ntemper .eq. 0) call errore('tranpsort', &
      'ftemper is not specified, missing temperatures and chemical potential',1)

   !load data from hdf5 file
   fname = trim(prefix)//"_tdf.h5"
   inquire(file=trim(fname), exist=has_file)
   if(.not. has_file) call errore('trans_postproc','missing '// trim(fname), 1)
   call hdf_open_file(file_id, trim(fname), status='OLD', action='READ')
   
   !read temperature
   call hdf_read_attribute(file_id, 'energy_grid', 'ne', nestep)
   call hdf_read_attribute(file_id, 'temperature', 'nt', ntmp)
   call hdf_read_attribute(file_id, 'dos', 'hole', itmp)
   is_hole = merge(.true., .false., itmp .eq. 1)
   allocate( tmpr(ntemper), efem(ntemper), niter(ntemper), ene(nestep), dos(nestep) )
   !
   write(msg,'(a)') "provided temperatures are different from those in " // trim(fname)
   if(ntmp .ne. ntemper) call errore('trans_postproc', trim(msg), 1)
   
   call hdf_read_dataset(file_id, 'temperature', tmpr)
   if( any(abs(tmpr-temper) > 1.0E-9_dp) ) call errore('trans_postproc', trim(msg), 1)
    
   call hdf_read_dataset(file_id, 'efermi', efem)
   write(msg,'(a)') "input chemical potentials are different from those in " // trim(fname)
   if( any(abs(efem - efermi) > 1.0E-6_dp) ) then
      write(stdout, '(5x, a)') "Warn (trans_postproc): "
      write(stdout, '(9x, a)') trim(msg)
   endif
   !
   deallocate( tmpr, efem )
   
   call hdf_read_dataset(file_id, 'dos', dos)
   call hdf_read_dataset(file_id, 'num_iter', niter)
   call hdf_read_dataset(file_id, 'energy_grid', ene)
   call hdf_get_dims(file_id, 'tdf_t1_i1', tdf_size)
   de = (ene(nestep)-ene(1)) / (nestep-1.0_dp)
   ncomp = tdf_size(1)
   if(ncomp.ne.6 .and. ncomp.ne.12) call errore('trans_postproc','wrong dimension in the tdf data',1)
   !
   allocate( tdf_tmp(6, nestep), tdf(ncomp, nestep), occup(nestep), density(ntemper) )
   allocate( cond(6, ntemper), seebeck(6, ntemper), therm_cond(6, ntemper) )

   do it = 1, ntemper
      tp = temper(it) !temperature
      ef = efermi(it) !chemical potential
      !load tdf in Rydberg atomic unit
      ! TDF_ij(E) = spinor/(N*V) sum_{nk} vnk*vnk*tau*delta(E-e_nk)
      call read_tdf_hdf5(file_id, it, niter(it), nestep, tdf)
      
      !compute carrier concentration: #. carrier / bohr^3
      occup(:) = fermi(tp, merge(ef-ene(:), ene(:)-ef, is_hole) )
      density(it) = dot_product(dos, occup) * de / volume
      
      !compute conductivity
      do i = 1, 6
         tdf_tmp(i,:)  = tdf(i,:) * mfermi_deriv(tp, (ene(:) - ef))
      enddo
      !\sigma = q^2 * Int tdf(E) * (-df/dE) * dE
      cond(1:6, it) = sum(tdf_tmp, dim=2) * de

      ! alpha corresponding to T * \sigma * seebeck, using tdf of E-field
      !  = - q Int tdf(E) * (-df/dE) * (E-\mu) * dE   (Eq.7c in PRB 94, 085204, 2016)
      do i = 1, 6
         tdf_tmp(i,:) = tdf_tmp(i,:) * (ene(:) - ef)
      enddo
      alpha(1:6) = - sum(tdf_tmp, dim=2) * de
   
      !
      ! \sigma * seebeck = - q/T Int tdf(E) * (-df/dE) * (E-\mu) * dE
      if(ncomp > 6) then
         !compute \sigma*seebeck, using tdf of T-field (note the minus sign below)
         do i = 1, 6
            tdf_tmp(i,:)  = tdf(i+6,:) * mfermi_deriv(tp, (ene(:) - ef))
         enddo
         cond_seebeck(1:6) = - sum(tdf_tmp, dim=2) * ( de / tp )
      else
         !!\sigma * seebeck = - q/T Int tdf(E) * (-df/dE) * (E-\mu) * dE
         !!cond_seebeck(1:6) = - sum(tdf_tmp, dim=2) * ( de / tp )
         cond_seebeck = alpha / tp
      endif

      !compute K for thermal conductivity, K = 1/T Int tdf(E) * (-df/dE) * (E-\mu)^2 * dE
      ! (using tdf of T-field (beta in Eq.7d PRB 94,085204,2016) if available.)
      do i = 1, 6
         tdf_tmp(i,:) = tdf_tmp(i,:) * (ene(:) - ef)
      enddo
      kk_coef(1:6) = sum(tdf_tmp, dim=2) * ( de / tp )

      ! seebeck = sigma^-1 * [sigma*seebeck]
      ! 1, K_B/q is omitted, q is electron charge and K_b is boltzmann constant 
      ! 2, the vlue of seebeck here is dimensionless.
      !
      ! thermal conductivity \kappa = K - alpha * seebeck (or K - T*[sigma*seebeck]*seebeck)
      call extract_trans_coeff &
         (tp, cond(:,it), cond_seebeck, kk_coef, seebeck(:,it), therm_cond(:,it), alpha)
   enddo
   call output_trans_coef(temper, efermi, density, cond, seebeck, therm_cond)

   deallocate( tdf_tmp, tdf, occup, density, cond, seebeck, therm_cond )
   return

   contains

      subroutine read_tdf_hdf5(fid, tid, iter, nstep, tdf)
         use pert_const, only: dp
         use HDF5_utils
         implicit none
         integer(HID_T), intent(in) :: fid
         integer, intent(in) :: tid, iter, nstep
         real(dp), intent(out) :: tdf(:,:)
         !
         character(len=120) :: dset_name
         character(len=6), external :: int_to_char

         dset_name = 'tdf_t' // trim(int_to_char(it)) // '_i' // trim(int_to_char(iter))
         call hdf_read_dataset(fid, trim(dset_name), tdf(:, 1:nstep))
      end subroutine read_tdf_hdf5
      !
end subroutine trans_postproc

!
! Copyright (C) 2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! adapted from QE/PW/src/pw_restart_new.f90/read_collected_to_evc 
!----------------------------------------------------------------------------
! adapted from QE/PW/src/pw_restart_new.f90/read_collected_to_evc 
!  modified to skip save_buffer, and only read one kpoint a time.
!------------------------------------------------------------------------
SUBROUTINE pw_read_collected_to_evc(dirname, iunit, ik, ik_g, isk, ngk_g, evc)
   !------------------------------------------------------------------------
   !
   ! ... This routines reads wavefunctions from the new file format and
   ! ... writes them into the old format
   !
   ! jjzhou: modified to skip save_buffer, and only read one kpoint a time.
   !
   use kinds,  only: dp
   !USE control_flags,        ONLY : gamma_only
   USE lsda_mod,             ONLY : nspin !, isk
   USE klist,                ONLY : nkstot, wk, nks, xk, ngk, igk_k
   USE wvfct,                ONLY : npwx, g2kin, et, wg, nbnd
   !USE wavefunctions, ONLY : evc
   !USE io_files,             ONLY : nwordwfc, iunwfc
   !USE buffers,              ONLY : save_buffer
   USE gvect,                ONLY : ig_l2g
   USE noncollin_module,     ONLY : noncolin, npol
   USE mp_bands,             ONLY : nbgrp, root_bgrp, intra_bgrp_comm
   USE mp_pools,             ONLY : me_pool, root_pool, &
                                    intra_pool_comm, inter_pool_comm
   USE mp,                   ONLY : mp_sum, mp_max
   USE io_base,              ONLY : read_wfc
   !
   IMPLICIT NONE
   !
   CHARACTER(LEN=*), INTENT(IN)  :: dirname
   INTEGER, INTENT(in) :: iunit, ik, ik_g, ngk_g, isk
   complex(dp), INTENT(out) :: evc(npwx*npol, nbnd)
   ! just work arry
   !INTEGER, intent(out) :: igk_l2g(npwx), igk_l2g_kdip(npwx), mill_k(3, npwx)
   !
   !
   CHARACTER(LEN=2), DIMENSION(2) :: updw = (/ 'up', 'dw' /)
   CHARACTER(LEN=320)   :: filename, msg
   !INTEGER              :: i, ik, ik_g, ig, ipol, ik_s
   INTEGER              :: i, ig, ipol, ik_s
   INTEGER              :: npol_, npwx_g, nbnd_
   INTEGER              :: nupdwn(2), ike, iks, npw_g, ispin
   INTEGER, EXTERNAL    :: global_kpoint_index
   !INTEGER, ALLOCATABLE :: ngk_g(:) ! the input argument ngk_g corresponds to ngk_g(ik_g)
   INTEGER, ALLOCATABLE :: igk_l2g(:), igk_l2g_kdip(:), mill_k(:,:)
   LOGICAL              :: opnd, ionode_k, gamma_only
   REAL(DP)             :: scalef, xk_(3), b1(3), b2(3), b3(3)
   character(len=6), external :: int_to_char
   !
   !iks = global_kpoint_index (nkstot, 1)
   !ike = iks + nks - 1
   !
   ! ... ngk_g: global number of k+G vectors for all k points
   !
   !ALLOCATE( ngk_g( nkstot ) )
   !ngk_g = 0
   !ngk_g(iks:ike) = ngk(1:nks)
   !CALL mp_sum( ngk_g, inter_pool_comm)
   !CALL mp_sum( ngk_g, intra_pool_comm)
   !ngk_g = ngk_g / nbgrp
   !
   ! ... npwx_g: maximum number of G vector among all k points
   !
   !npwx_g = MAXVAL( ngk_g(1:nkstot) )
   !
   ! ... the root processor of each pool reads
   !
   !ionode_k = (me_pool == root_pool)
   !
   ! ... The igk_l2g array yields the correspondence between the
   ! ... local k+G index and the global G index
   !
   ALLOCATE ( igk_l2g( npwx ) )
   !
   ! ... the igk_l2g_kdip local-to-global map is needed to read wfcs
   !
   ALLOCATE ( igk_l2g_kdip( npwx ) )
   !
   ALLOCATE( mill_k ( 3,npwx ) )
   !
   !k_points_loop: DO ik = 1, nks
      !
      ! index of k-point ik in the global list
      !
      !ik_g = ik + iks - 1
      !
      ! ... Compute the igk_l2g array from previously computed arrays
      ! ... igk_k (k+G indices) and ig_l2g (local to global G index map)
      !
      igk_l2g = 0
      DO ig = 1, ngk(ik)
         igk_l2g(ig) = ig_l2g(igk_k(ig,ik))
      END DO
      !
      ! ... npw_g: the maximum G vector index among all processors
      !
      npw_g = MAXVAL( igk_l2g(1:ngk(ik)) )
      CALL mp_max( npw_g, intra_pool_comm )
      !
      igk_l2g_kdip = 0
      CALL pw_gk_l2gmap_kdip( npw_g, ngk_g, ngk(ik), igk_l2g, &
                           igk_l2g_kdip )
      !
      evc=(0.0_DP, 0.0_DP)
      
      IF ( nspin == 2 ) THEN
         !
         ! ... LSDA: spin mapped to k-points, isk(ik) tracks up and down spin
         !
         !ik_g = MOD ( ik_g-1, nkstot/2 ) + 1 
         ispin = isk !(ik)
         filename = TRIM(dirname) // 'wfc' // updw(ispin) // &
              & TRIM(int_to_char(ik_g))
         !
      ELSE
         !
         filename = TRIM(dirname) // 'wfc' // TRIM(int_to_char(ik_g))
         !
      ENDIF
      !
      CALL read_wfc( iunit, filename, root_bgrp, intra_bgrp_comm, &
           ik_g, xk_, ispin, npol_, evc, npw_g, gamma_only, nbnd_, &
           igk_l2g_kdip(:), ngk(ik), b1, b2, b3, mill_k, scalef )
      !
      ! ... here one should check for consistency between what is read
      ! ... and what is expected
      !
      IF ( nbnd_ < nbnd ) THEN
         WRITE (msg,'("The number of bands for this run is",I6,", but only",&
              & I6," bands were read from file")')  nbnd, nbnd_  
         CALL errore ('pw_restart - read_collected_to_evc', msg, 1 )
      END IF
      !CALL save_buffer ( evc, nwordwfc, iunwfc, ik )
      ! 
   !END DO k_points_loop
   !
   DEALLOCATE ( mill_k )
   DEALLOCATE ( igk_l2g )
   DEALLOCATE ( igk_l2g_kdip )
   !!
   RETURN
   !
END SUBROUTINE pw_read_collected_to_evc


SUBROUTINE pw_gk_l2gmap_kdip( npw_g, ngk_g, ngk, igk_l2g, igk_l2g_kdip)
   !-----------------------------------------------------------------------
   !
   ! ... This subroutine maps local G+k index to the global G vector index
   ! ... the mapping is used to collect wavefunctions subsets distributed
   ! ... across processors.
   ! ... This map is used to obtained the G+k grids related to each kpt
   !
   USE mp_bands,             ONLY : intra_bgrp_comm
   USE mp,                   ONLY : mp_sum
   !
   IMPLICIT NONE
   !
   ! ... Here the dummy variables
   !
   INTEGER, INTENT(IN)  :: npw_g, ngk_g, ngk
   INTEGER, INTENT(IN)  :: igk_l2g(ngk)
   INTEGER, INTENT(OUT) :: igk_l2g_kdip(ngk)
   !
   INTEGER, ALLOCATABLE :: igwk_(:), itmp(:), igwk_lup(:)
   INTEGER              :: ig, ig_, ngg
   !
   !
   ALLOCATE( itmp( npw_g ) )
   ALLOCATE( igwk_( ngk_g ) )
   !
   itmp(:)  = 0
   igwk_(:) = 0
   !
   DO ig = 1, ngk
      itmp(igk_l2g(ig)) = igk_l2g(ig)
   END DO
   !
   CALL mp_sum( itmp, intra_bgrp_comm )
   !
   ngg = 0
   DO ig = 1, npw_g
      !
      IF ( itmp(ig) == ig ) THEN
         !
         ngg = ngg + 1
         igwk_(ngg) = ig
         !
      END IF
      !
   END DO
   !
   IF ( ngg /= ngk_g ) &
      CALL errore( 'gk_l2gmap_kdip', 'unexpected dimension in ngg', 1 )
   !
   ALLOCATE( igwk_lup( npw_g ) )
   !
!$omp parallel private(ig_, ig)
!$omp workshare
   igwk_lup = 0
!$omp end workshare
!$omp do
   DO ig_ = 1, ngk_g
      igwk_lup(igwk_(ig_)) = ig_
   END DO
!$omp end do
!$omp do
   DO ig = 1, ngk
      igk_l2g_kdip(ig) = igwk_lup(igk_l2g(ig))
   END DO
!$omp end do
!$omp end parallel
   !
   DEALLOCATE( igwk_lup )
   !
   DEALLOCATE( itmp, igwk_ )
   !
   RETURN
   !
END SUBROUTINE pw_gk_l2gmap_kdip

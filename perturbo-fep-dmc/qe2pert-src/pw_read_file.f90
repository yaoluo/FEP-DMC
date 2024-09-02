!
! Copyright (C) 2016-2019 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! adapted from PW/src/read_file_new.f90
! Wrapper routine to read xml data file and initialze, 
!  without allocate space for wavefunction
!----------------------------------------------------------------------------

#if defined(__QE64)

!work with QE version 6.4
SUBROUTINE pw_read_file()
  !----------------------------------------------------------------------------
  !
  ! Read data produced by pw.x or cp.x - new xml file and binary files
  ! Wrapper routine for backwards compatibility
  ! 
  ! jjzhou: modified to skip generating prefix.wfc_xx files
  !
  USE io_files,             ONLY : nwordwfc, iunwfc, prefix, tmp_dir, wfc_dir
  USE io_global,            ONLY : stdout, ionode
  USE buffers,              ONLY : open_buffer, close_buffer
  USE wvfct,                ONLY : nbnd, npwx
  USE noncollin_module,     ONLY : npol
  USE paw_variables,        ONLY : okpaw, ddd_PAW
  USE paw_onecenter,        ONLY : paw_potential
  USE uspp,                 ONLY : becsum
  USE scf,                  ONLY : rho
  USE realus,               ONLY : betapointlist, &
                                   init_realspace_vars,real_space
  USE dfunct,               ONLY : newd
  USE ldaU,                 ONLY : lda_plus_u, U_projection
  USE pw_restart_new !,     ONLY : read_xml_file
  USE io_files,             ONLY : tmp_dir, prefix, postfix
  USE control_flags,        ONLY : io_level
  USE klist,                ONLY : init_igk
  USE gvect,                ONLY : ngm, g
  USE gvecw,                ONLY : gcutw
  USE qes_types_module,     ONLY : output_type
  !
  IMPLICIT NONE 
  TYPE ( output_type) :: output_obj 
  INTEGER :: ierr
  LOGICAL :: exst, wfc_is_collected
  CHARACTER( LEN=256 )  :: dirname
  !
  !
  ierr = 0 
  !
  ! ... Read the contents of the xml data file
  !
  dirname = TRIM( tmp_dir ) // TRIM( prefix ) // postfix
  IF ( ionode ) WRITE( stdout, '(/,5x,A,/,5x,A)') &
     'Reading data from directory:', TRIM( dirname )
  !
  CALL read_xml_file ( wfc_is_collected )
  !
  ! ... Open unit iunwfc, for Kohn-Sham orbitals - we assume that wfcs
  ! ... have been written to tmp_dir, not to a different directory!
  ! ... io_level = 1 so that a real file is opened
  !
  wfc_dir = tmp_dir
  nwordwfc = nbnd*npwx*npol
  io_level = 1
  !CALL open_buffer ( iunwfc, 'wfc', nwordwfc, io_level, exst )
  !
  ! ... Allocate and compute k+G indices and number of plane waves
  ! ... FIXME: should be read from file, not re-computed
  !
  CALL init_igk ( npwx, ngm, g, gcutw ) 
  !
  ! ... FIXME: this should be taken out from here
  !
  !IF ( wfc_is_collected ) CALL read_collected_to_evc(dirname) 
  !
  ! ... Assorted initialization: pseudopotentials, PAW
  ! ... Not sure which ones (if any) should be done here
  !
  CALL init_us_1()
  !
  IF (lda_plus_u .AND. (U_projection == 'pseudo')) CALL init_q_aeps()
  !
  IF (okpaw) THEN
     becsum = rho%bec
     CALL PAW_potential(rho%bec, ddd_PAW)
  ENDIF 
  !
  IF ( real_space ) THEN
    CALL betapointlist()
    CALL init_realspace_vars()
    IF( ionode ) WRITE(stdout,'(5x,"Real space initialisation completed")')
  ENDIF
  CALL newd()
  !
  !CALL close_buffer  ( iunwfc, 'KEEP' )
  !
END SUBROUTINE pw_read_file

#else

! work with QE version 6.5
SUBROUTINE pw_read_file()
   USE io_files,  ONLY : wfc_dir, tmp_dir
   USE klist,     ONLY : nkstot, nks, xk, wk, init_igk
   USE lsda_mod,  ONLY : isk
   USE wvfct,     ONLY : nbnd, et, wg, npwx
   USE gvecw,     ONLY : gcutw
   USE gvect,     ONLY : ngm, g
   !
   IMPLICIT NONE
   ! 
   LOGICAL :: wfc_is_collected
   INTEGER, EXTERNAL :: n_plane_waves

   wfc_is_collected = .false.
   CALL read_file_new( wfc_is_collected )

   if( .not. wfc_is_collected ) call errore('pw_read_file', &
     'Wavefunctions in collected format not available', 1)

   !performs wavefunction-related initialization without allocate space
   !
   ! ... initialization of KS orbitals
   !
   wfc_dir = tmp_dir ! this is likely obsolete and no longer used
   !
   ! ... distribute across pools k-points and related variables.
   ! ... nks is defined by the following routine as the number 
   ! ... of k-points in the current pool
   !
   CALL divide_et_impera( nkstot, xk, wk, isk, nks )
   CALL poolscatter( nbnd, nkstot, et, nks, et )
   CALL poolscatter( nbnd, nkstot, wg, nks, wg )

   !
   !   calculate number of PWs for all kpoints
   !
   npwx = n_plane_waves( gcutw, nks, xk, g, ngm )
   !
   !   compute indices j=igk(i) such that (k+G)_i = k+G_j, for all k
   !   compute number of plane waves ngk(ik) as well
   !
   CALL init_igk( npwx, ngm, g, gcutw )
   
   RETURN
END SUBROUTINE pw_read_file

#endif

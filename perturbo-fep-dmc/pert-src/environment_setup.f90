!
!!set up enviroment, adapted from QE/Modules/environment_start
!-------------------------------------------------------------
! Copyright (C) 2002-2011 Quantum ESPRESSO groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .

subroutine environment_setup(code, lclock)
   use io_files,  only: crash_file, nd_nmbr
   use mp_images, only: me_image, root_image, my_image_id
   use mp_pools,  only: npool
   use mp_world,  only: nproc
   USE io_global, ONLY: stdout, meta_ionode
   use fox_init_module, ONLY: fox_init
   !
#if defined(__HDF5)
   use qeh5_base_module, only: initialize_hdf5
#else
   use hdf5_utils, only: hdf_init
#endif

   implicit none
   CHARACTER(LEN=*), INTENT(IN) :: code
   LOGICAL, intent(in) :: lclock
  
   LOGICAL :: exst, debug = .false.
   INTEGER :: ios, crashunit
   CHARACTER(LEN=80) :: uname
   CHARACTER(LEN=9)  :: cdate, ctime
   
   INTEGER, EXTERNAL :: find_free_unit
   CHARACTER(LEN=6), EXTERNAL :: int_to_char
   
   ! ... The Intel compiler allocates a lot of stack space
   ! ... Stack limit is often small, thus causing SIGSEGV and crash
   ! ... One may use "ulimit -s unlimited" but it doesn't always work
   ! ... The following call does the same and always works
   ! 
#if defined(__INTEL_COMPILER)
   CALL remove_stack_limit()
#endif

   ! ... use ".FALSE." to disable all clocks except the total cpu time clock
   ! ... use ".TRUE."  to enable clocks
   CALL init_clocks( lclock )
   CALL start_clock( trim(code) )

   IF( meta_ionode ) THEN
      ! ...  search for file CRASH and delete it
      INQUIRE( FILE=TRIM(crash_file), EXIST=exst )
      IF( exst ) THEN
         crashunit = find_free_unit()
         OPEN( UNIT=crashunit, FILE=TRIM(crash_file), STATUS='OLD',IOSTAT=ios )
         IF (ios==0) THEN
            CLOSE( UNIT=crashunit, STATUS='DELETE', IOSTAT=ios )
         ELSE
            WRITE(stdout,'(5x,"Remark: CRASH file could not be deleted")')
         END IF
      END IF

   ELSE
      ! ... one processor per image (other than meta_ionode)
      ! ... or, for debugging purposes, all processors,
      ! ... open their own standard output file
!#define DEBUG
#if defined(DEBUG)
      debug = .true.
#endif
      IF (me_image == root_image .OR. debug ) THEN
         uname = 'out.' // trim(int_to_char( my_image_id )) // '_' // &
              trim(int_to_char( me_image))
         OPEN ( unit = stdout, file = TRIM(uname),status='unknown')
      ELSE
#if defined(_WIN32)
         OPEN ( unit = stdout, file='NUL:', status='unknown' )
#else
         OPEN ( unit = stdout, file='/dev/null', status='unknown' )
#endif
      END IF

   END IF

   !print some message
   call date_and_tim(cdate, ctime)
   WRITE( stdout, '(/3X,"Program ",A," starts on ",A9," at ",A9)' ) &
      TRIM(code), cdate, ctime

   ! ... for compatibility with PWSCF
#if defined(__MPI)
   nd_nmbr = TRIM ( int_to_char( me_image+1 ))
   CALL parallel_info(code, stdout)
#else
   nd_nmbr = ' '
   CALL serial_info(stdout)
#endif
   
   !open HDF5 interface
   !if __HDF5 is used in QE/make.inc, 
   !  then we initialized the HDF5 interface with qeh5_module subroutine
   !if __HDF5 is not used in QE, then we use our own HDF5 wrapper.
#if defined(__HDF5)
   call initialize_hdf5()
#else
   call hdf_init()
#endif

   CALL fox_init()
   !check if npool == nproc
   if (nproc /= npool) &
      CALL errore('qe2pert','npools must match the number of MPI processes',1)
end subroutine environment_setup


subroutine environment_stop( code )
   USE io_global, ONLY: stdout, meta_ionode
   !
#if defined(__HDF5)
   use qeh5_base_module, only: finalize_hdf5
#else
   use hdf5_utils, only: hdf_finalize
#endif

   implicit none
   CHARACTER(LEN=*), INTENT(IN) :: code
   !
   CHARACTER(LEN=9)  :: cdate, ctime
   CHARACTER(LEN=80) :: time_str

   CALL stop_clock(  TRIM(code) )
   !CALL print_clock( TRIM(code) )
   !print out all the clock
   CALL print_clock(' ')
   
   !close HDF5 interface
#if defined(__HDF5)
   call finalize_hdf5()
#else
   call hdf_finalize()
#endif

   call date_and_tim( cdate, ctime )
   time_str = 'Program was terminated on:  ' // ctime // ' ' // cdate

   if(meta_ionode) write(stdout,'(/3X, A60)') time_str
end subroutine environment_stop


  !==-----------------------------------------------------------------------==!
  SUBROUTINE parallel_info(code, stdout)
    use mp_world,  only: nproc, nnode
    implicit none
    CHARACTER(LEN=*), INTENT(IN) :: code
    integer, intent(in) :: stdout

#if defined(_OPENMP)
    INTEGER, EXTERNAL :: omp_get_max_threads
    !
    WRITE( stdout, '(/5X,"Parallel version (MPI & OpenMP), running on ", &
         &I7," processor cores")' ) nproc * omp_get_max_threads()
    !
    WRITE( stdout, '(5X,"Number of MPI processes:           ",I7)' ) nproc
    !
    WRITE( stdout, '(5X,"Threads/MPI process:               ",I7/)' ) &
         omp_get_max_threads()
#else
    WRITE( stdout, '(/5X,"Parallel version (MPI), running on ",&
         &I5," processors")' ) nproc 
#endif
    !
#if !defined(__GFORTRAN__) ||  ((__GNUC__>4) || ((__GNUC__==4) && (__GNUC_MINOR__>=8)))
    WRITE( stdout, '(5X,"MPI processes distributed on ",&
         &I5," nodes",/)' ) nnode
#endif
    !
  END SUBROUTINE parallel_info

  !==-----------------------------------------------------------------------==!
  SUBROUTINE serial_info ( stdout )
    implicit none
    integer, intent(in) :: stdout
    !
#if defined(_OPENMP)
    INTEGER, EXTERNAL :: omp_get_max_threads
#endif
    !
#if defined(_OPENMP)
    WRITE( stdout, '(/5X,"Serial multi-threaded version, running on ",&
         &I4," processor cores")' ) omp_get_max_threads()
    !
#else
    WRITE( stdout, '(/5X,"Serial version")' )
#endif
    !
  END SUBROUTINE serial_info

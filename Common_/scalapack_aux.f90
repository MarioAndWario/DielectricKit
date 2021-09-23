!================================================================================
!
! Modules:
!
! (1) scalapack_aux_m      Originally by FHJ 07/24/2015
!
!     Defines functions associated to matrix distributions. These functions
!     were copied from netlib ScaLAPACK (BSD-licensed), and are included here
!     as they can be used on builds with and without ScaLAPACK. None of these
!     functions require ScaLAPACK.
!
!================================================================================

#include "f_defs.h"

module scalapack_aux_m

  implicit none

  private

  public ::      &
    numroc,           &
    indxl2g,          &
    indxg2l,          &
    indxg2p


contains

! -- ScaLAPACK tools routine (version 1.7) --

  INTEGER FUNCTION NUMROC( N, NB, IPROC, ISRCPROC, NPROCS )
    integer, intent(in) :: N, NB, IPROC, ISRCPROC, NPROCS

    integer :: EXTRABLKS, MYDIST, NBLOCKS
    MYDIST = MOD( NPROCS+IPROC-ISRCPROC, NPROCS )
    NBLOCKS = N / NB
    NUMROC = (NBLOCKS/NPROCS) * NB
    EXTRABLKS = MOD( NBLOCKS, NPROCS )
    IF( MYDIST.LT.EXTRABLKS ) THEN
      NUMROC = NUMROC + NB
    ELSE IF( MYDIST.EQ.EXTRABLKS ) THEN
      NUMROC = NUMROC + MOD( N, NB )
    END IF

  END FUNCTION NUMROC

  INTEGER FUNCTION INDXL2G( INDXLOC, NB, IPROC, ISRCPROC, NPROCS )
    integer, intent(in) :: INDXLOC, IPROC, ISRCPROC, NB, NPROCS

    INDXL2G = NPROCS*NB*((INDXLOC-1)/NB) + MOD(INDXLOC-1,NB) + &
      MOD(NPROCS+IPROC-ISRCPROC, NPROCS)*NB + 1

  END FUNCTION INDXL2G

  INTEGER FUNCTION INDXG2L( INDXGLOB, NB, IPROC, ISRCPROC, NPROCS )
    integer, intent(in) :: INDXGLOB, IPROC, ISRCPROC, NB, NPROCS

    INDXG2L = NB*((INDXGLOB-1)/(NB*NPROCS))+MOD(INDXGLOB-1,NB)+1

  END FUNCTION INDXG2L

  INTEGER FUNCTION INDXG2P( INDXGLOB, NB, IPROC, ISRCPROC, NPROCS )
    integer, intent(in) :: INDXGLOB, IPROC, ISRCPROC, NB, NPROCS

    INDXG2P = MOD( ISRCPROC + (INDXGLOB - 1) / NB, NPROCS )

  END FUNCTION INDXG2P
  
end module scalapack_aux_m

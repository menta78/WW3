!> @file
!> @brief Calculate ice source term S_{ice} using UOST method.
!>
!> @author Lorenzo Mentaschi
!> @date   14-Feb-2026
!>

#include "w3macros.h"
!/ ------------------------------------------------------------------- /
!>
!> @brief Calculate ice source term S_{ice} using UOST method.
!>
!> @author Lorenzo Mentaschi
!> @date   14-Feb-2026
!>
!> @copyright Copyright 2009-2022 National Weather Service (NWS),
!>       National Oceanic and Atmospheric Administration.  All rights
!>       reserved.  WAVEWATCH III is a trademark of the NWS.
!>       No unauthorized use without permission.
!>
MODULE W3SIC6MD
  !/
  !/                  +-----------------------------------+
  !/                  | WAVEWATCH III           NOAA/NCEP |
  !/                  |      Lorenzo Mentaschi            |
  !/                  |                        FORTRAN 90 |
  !/                  | Last update :         14-Feb-2026 |
  !/                  +-----------------------------------+
  !/
  !/    For updates see W3SIC6 documentation.
  !/
  !  1. Purpose :
  !
  !     Calculate ice source term S_{ice} using UOST (Unresolved 
  !     Obstacles Source Term) method.
  !
  !     First-order approximation of the wave attenuation due to ice 
  !     concentration. ICEUODIS stands for Ice Unresolved Obstacles 
  !     Dissipation, based on UOST (Mentaschi et al. 2015, 2018, 2020).
  !     
  !     From the ice concentration an isotropic transparency coefficient 
  !     is estimated. To simplify the problem, for now it is assumed 
  !     that beta==alpha which means that all the energy is dissipated 
  !     in the current cell. This gets rid of the shadow problems.
  !     This is inaccurate for the current cell, but should dissipate 
  !     all the energy needed. Furthermore, local wave growth (which 
  !     reduces the effect of unresolved obstacles) is neglected.
  !
  !  2. Variables and types :
  !
  !  3. Subroutines and functions :
  !
  !      Name      Type  Scope    Description
  !     ----------------------------------------------------------------
  !      W3SIC6    Subr. Public   Ice source term using UOST.
  !     ----------------------------------------------------------------
  !
  !  4. Subroutines and functions used :
  !
  !     See subroutine documentation.
  !
  !  5. Remarks :
  !
  !     In the future, if this approximation will not be enough, 
  !     some improvement could be introduced:
  !     - beta different from alpha could be estimated by assuming 
  !       a uniform distribution of the ice in the cell, or by loading
  !       the distribution of the ice, and the shadow could be estimated 
  !       for the neighboring cells.
  !     - local wave growth could be taken into account (see the psi 
  !       function in UOST)
  !     - if there is information on the size of the ice floes, the 
  !       transparency coeff. could be made frequency-dependent.
  !
  !     Reference: Mentaschi et al. (2015, 2018, 2020)
  !
  !  6. Switches :
  !
  !     See subroutine documentation.
  !
  !  7. Source code :
  !/
  !/ ------------------------------------------------------------------- /
  !/
  PUBLIC :: W3SIC6
  !/
CONTAINS
  !/ ------------------------------------------------------------------- /
  !>
  !> @brief S_{ice} source term using UOST method for unresolved obstacles.
  !>
  !> @details First-order approximation of wave attenuation due to ice
  !>  concentration using the transparency coefficient approach.
  !>
  !> @param[in]  A      Action density spectrum (1-D).
  !> @param[in]  CG     Group velocities.
  !> @param[in]  IX     Grid index.
  !> @param[in]  IY     Grid index.
  !> @param[out] S      Source term (1-D version).
  !> @param[out] D      Diagonal term of derivative (1-D version).
  !>
  !> @author Lorenzo Mentaschi
  !> @date   14-Feb-2026
  !>
  SUBROUTINE W3SIC6 (A, CG, IX, IY, S, D)
    !/
    !/                  +-----------------------------------+
    !/                  | WAVEWATCH III           NOAA/NCEP |
    !/                  |      Lorenzo Mentaschi            |
    !/                  |                        FORTRAN 90 |
    !/                  | Last update :         14-Feb-2026 |
    !/                  +-----------------------------------+
    !/
    !/    14-Feb-2026 : Origination.                        ( version 7.xx )
    !/                                                       (L. Mentaschi)
    !/
    !  1. Purpose :
    !
    !     S_{ice} source term using UOST method.
    !
    !/ ------------------------------------------------------------------- /
    !
    !  2. Method :
    !
    !     From the ice concentration an isotropic transparency coefficient
    !     is estimated. The obstruction section is proportional to the 
    !     square root of ice concentration. The transparency coefficient
    !     alpha is 1-obstruction, and for simplicity beta==alpha.
    !     This means all energy is dissipated in the current cell.
    !
    !     The gamma parameter is computed as:
    !       gamma = (1 - beta) / beta
    !     And is capped at GAMMAUP = 200.
    !
    !     The diagonal term is:
    !       D(IK) = - CG(IK) / CELLSIZE * gamma
    !     And the source term is:
    !       S = D * A
    !
    !     Cell area is obtained from GSQRT (the metric tensor determinant)
    !     which properly handles all grid types and coordinate systems.
    !     Cell size is computed as the radius of an equivalent circle
    !     with the same area as the grid cell.
    !
    !  3. Parameters :
    !
    !     Parameter list
    !     ----------------------------------------------------------------
    !       A       R.A.  I   Action density spectrum (1-D)
    !       CG      R.A.  I   Group velocities.
    !       IX,IY   I.S.  I   Grid indices.
    !       S       R.A.  O   Source term (1-D version).
    !       D       R.A.  O   Diagonal term of derivative (1-D version).
    !     ----------------------------------------------------------------
    !
    !  4. Subroutines used :
    !
    !      Name      Type  Module   Description
    !     ----------------------------------------------------------------
    !      STRACE    Subr. W3SERVMD Subroutine tracing (!/S switch).
    !      PRT2DS    Subr. W3ARRYMD Print plot output (!/T1 switch).
    !      OUTMAT    Subr. W3ARRYMD Matrix output (!/T2 switch).
    !     ----------------------------------------------------------------
    !
    !  5. Called by :
    !
    !      Name      Type  Module   Description
    !     ----------------------------------------------------------------
    !      W3SRCE    Subr. W3SRCEMD Source term integration.
    !      W3EXPO    Subr.   N/A    ASCII Point output post-processor.
    !      W3EXNC    Subr.   N/A    NetCDF Point output post-processor.
    !      GXEXPO    Subr.   N/A    GrADS point output post-processor.
    !     ----------------------------------------------------------------
    !
    !  6. Error messages :
    !
    !     None.
    !
    !  7. Remarks :
    !
    !     If ice concentration is below threshold, no calculations are made.
    !
    !  8. Structure :
    !
    !     See source code.
    !
    !  9. Switches :
    !
    !     !/S   Enable subroutine tracing.
    !     !/T   Enable general test output.
    !     !/T0  2-D print plot of source term.
    !     !/T1  Print arrays.
    !
    ! 10. Source code :
    !
    !/ ------------------------------------------------------------------- /
    USE CONSTANTS, ONLY: TPI, PI
    USE W3ODATMD, ONLY: NDSE
    USE W3SERVMD, ONLY: EXTCDE
    USE W3GDATMD, ONLY: NK, NTH, NSPEC, SIG, MAPWN, GSQRT
    USE W3IDATMD, ONLY: ICEI, INFLAGS2
#ifdef W3_T
    USE W3ODATMD, ONLY: NDST
#endif
#ifdef W3_S
    USE W3SERVMD, ONLY: STRACE
#endif
#ifdef W3_T0
    USE W3ARRYMD, ONLY: PRT2DS
#endif
#ifdef W3_T1
    USE W3ARRYMD, ONLY: OUTMAT
#endif
    !
    IMPLICIT NONE
    !/
    !/ ------------------------------------------------------------------- /
    !/ Parameter list
    REAL, INTENT(IN)        :: CG(NK),   A(NSPEC)
    REAL, INTENT(OUT)       :: S(NSPEC), D(NSPEC)
    INTEGER, INTENT(IN)     :: IX, IY
    !/
    !/ ------------------------------------------------------------------- /
    !/ Local parameters
    !/
#ifdef W3_S
    INTEGER, SAVE           :: IENT = 0
#endif
#ifdef W3_T0
    INTEGER                 :: ITH
    REAL                    :: DOUT(NK,NTH)
#endif
    INTEGER                 :: IKTH, IK
    REAL                    :: D1D(NK)
    REAL                    :: ICEC, OBSTSECTION, BETA, CELLAREA, CELLSIZE
    REAL                    :: CGI, GAM
    REAL, PARAMETER         :: GAMMAUP = 200.0
    REAL, PARAMETER         :: ICETHR = 1.0E-6
    !/
    !/ ------------------------------------------------------------------- /
    !/
#ifdef W3_S
    CALL STRACE (IENT, 'W3SIC6')
#endif
    !
    ! 0.  Initializations ------------------------------------------------ *
    !
    D        = 0.0
    S        = 0.0
    D1D      = 0.0
    !
    ! Check if ice concentration field is available
    IF (.NOT.INFLAGS2(4)) THEN
      WRITE (NDSE,1001) 'ICE CONCENTRATION'
      CALL EXTCDE(2)
    ENDIF
    !
    ! Get ice concentration at this grid point
    ICEC = ICEI(IX,IY)
    ICEC = MAX(MIN(ICEC, 1.0), 0.0)
    !
    ! 1.  No ice or very low ice concentration --------------------------- /
    !
    IF ( ICEC <= ICETHR ) THEN
      D = 0.0
      S = 0.0
      RETURN
    END IF
    !
    ! 2.  Calculate ice dissipation -------------------------------------- /
    !
    ! 2.a Compute cell area ---------------------------------------------- /
    !
    ! Use GSQRT which is the area element (sqrt of metric tensor determinant)
    ! This works for all grid types and properly handles coordinate systems
    CELLAREA = ABS(GSQRT(IY,IX))
    !
    ! Compute cell size as radius of equivalent circle
    CELLSIZE = SQRT(CELLAREA / PI)
    !
#ifdef W3_T38
    WRITE (NDST,9000) ICEC, CELLAREA, CELLSIZE
#endif
    !
    ! 2.b Compute transparency coefficient ------------------------------- /
    !
    ! The total obstruction coefficient is given by sqrt(concentration)
    ! The transparency alpha is given by 1-obstruction
    ! Here we assume that beta==alpha
    OBSTSECTION = SQRT(ICEC)
    BETA = 1.0 - OBSTSECTION
    !
    ! Compute gamma parameter
    GAM = (1.0 - BETA) / BETA
    GAM = MIN(GAM, GAMMAUP)
    !
    ! 2.c Calculate diagonal term ---------------------------------------- /
    !
    DO IK=1, NK
      CGI = CG(IK)
      D1D(IK) = - CGI / CELLSIZE * GAM
    END DO
    !
    ! 2.d Fill diagonal matrix ------------------------------------------- /
    !
    DO IKTH=1, NSPEC
      D(IKTH) = D1D(MAPWN(IKTH))
    END DO
    !
    ! 2.e Calculate source term ------------------------------------------ /
    !
    S = D * A
    !
    ! ... Test output of arrays
    !
#ifdef W3_T0
    DO IK=1, NK
      DO ITH=1, NTH
        DOUT(IK,ITH) = D(ITH+(IK-1)*NTH)
      END DO
    END DO
    CALL PRT2DS (NDST, NK, NK, NTH, DOUT, SIG(1:), '  ', 1.,    &
         0.0, 0.001, 'Diag Sice6', ' ', 'NONAME')
#endif
    !
#ifdef W3_T1
    CALL OUTMAT (NDST, D, NTH, NTH, NK, 'diag Sice6')
#endif
    !
    ! Formats
    !
1001 FORMAT (/' *** WAVEWATCH III ERROR IN W3SIC6 : '/               &
         '     ',A,' REQUIRED BUT NOT SELECTED'/)
    !
#ifdef W3_T38
9000 FORMAT (' TEST W3SIC6 : ICEC, CELLAREA, CELLSIZE : ',3E10.3)
#endif
    !/
    !/ End of W3SIC6 ----------------------------------------------------- /
    !/
  END SUBROUTINE W3SIC6
  !/
  !/ End of module W3SIC6MD -------------------------------------------- /
  !/
END MODULE W3SIC6MD

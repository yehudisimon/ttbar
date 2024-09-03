C     Set of subroutines to project the Feynman diagrama amplitudes
C      AMPL onto the color flow space.

      SUBROUTINE ML5_2_COMPUTE_COLOR_FLOWS(HEL_MULT,DO_CUMULATIVE)
      IMPLICIT NONE
      LOGICAL DO_CUMULATIVE
      INTEGER HEL_MULT
      CALL ML5_2_REINITIALIZE_JAMPS()
      CALL ML5_2_DO_COMPUTE_COLOR_FLOWS(.FALSE.)
      ! HSS DEBUG
      CALL ML5_2_DO_COMPUTE_NEW_COLOR_FLOWS(.FALSE.)
      END SUBROUTINE

      SUBROUTINE ML5_2_COMPUTE_COLOR_FLOWS_DERIVED_QUANTITIES(HEL_MULT)
      IMPLICIT NONE
      INTEGER HEL_MULT
      CONTINUE
      END SUBROUTINE

      SUBROUTINE ML5_2_DEALLOCATE_COLOR_FLOWS()
      IMPLICIT NONE
      CALL ML5_2_DO_COMPUTE_COLOR_FLOWS(.TRUE.)
      ! HSS DEBUG
      CALL ML5_2_DO_COMPUTE_NEW_COLOR_FLOWS(.TRUE.)
      END SUBROUTINE

      SUBROUTINE ML5_2_DO_COMPUTE_COLOR_FLOWS(CLEANUP)
      IMPLICIT NONE
C     
C     CONSTANTS 
C     
      CHARACTER*512 PROC_PREFIX
      PARAMETER ( PROC_PREFIX='ML5_2_')
      CHARACTER*512 LOOPCOLORFLOWCOEFSNAME
      PARAMETER ( LOOPCOLORFLOWCOEFSNAME='LoopColorFlowCoefs.dat')
      CHARACTER*512 BORNCOLORFLOWCOEFSNAME
      PARAMETER ( BORNCOLORFLOWCOEFSNAME='BornColorFlowCoefs.dat')
      COMPLEX*16 IMAG1
      PARAMETER (IMAG1=(0.0D0,1.0D0))

      INTEGER NBORNAMPS
      PARAMETER (NBORNAMPS=3)
      INTEGER    NLOOPAMPS
      PARAMETER (NLOOPAMPS=120)
      INTEGER    NSQUAREDSO, NLOOPAMPSO
      PARAMETER (NSQUAREDSO=1, NLOOPAMPSO=2)
      INTEGER NLOOPFLOWS
      PARAMETER (NLOOPFLOWS=3)
      INTEGER NBORNFLOWS
      PARAMETER (NBORNFLOWS=2)
      INTEGER    NBORNAMPSO
      PARAMETER (NBORNAMPSO=NLOOPAMPSO)

C     
C     LOCAL VARIABLES 
C     
C     When this subroutine is called with CLEANUP=True, it deallocates
C      all its
C     arrays
      LOGICAL CLEANUP
C     Storing concatenated filenames
      CHARACTER*512 TMP
      CHARACTER*512 LOOPCOLORFLOWCOEFSN
      CHARACTER*512 BORNCOLORFLOWCOEFSN

      INTEGER I, J, K, SOINDEX, ARRAY_SIZE
      COMPLEX*16 PROJ_COEF


C     
C     FUNCTIONS
C     
      INTEGER ML5_2_ML5SOINDEX_FOR_BORN_AMP
      INTEGER ML5_2_ML5SOINDEX_FOR_LOOP_AMP
C     
C     GLOBAL VARIABLES
C     
      CHARACTER(512) MLPATH
      COMMON/MLPATH/MLPATH

      COMPLEX*16 AMP(NBORNAMPS)
      COMMON/ML5_2_AMPS/AMP
      COMPLEX*16 AMPL(3,NLOOPAMPS)
      COMMON/ML5_2_AMPL/AMPL
      COMPLEX*16 JAMPL(3,NLOOPFLOWS,NLOOPAMPSO)
      COMMON/ML5_2_JAMPL/JAMPL

      COMPLEX*16 JAMPB(NBORNFLOWS,NBORNAMPSO)
      COMMON/ML5_2_JAMPB/JAMPB

C     Now a more advanced data structure for storing the projection
C      coefficient
C     This is of course not F77 standard but widely supported by now.
      TYPE PROJCOEFFS
        INTEGER, DIMENSION(:), ALLOCATABLE :: NUM
        INTEGER, DIMENSION(:), ALLOCATABLE :: DENOM
        INTEGER, DIMENSION(:), ALLOCATABLE :: AMPID
      ENDTYPE PROJCOEFFS
      TYPE(PROJCOEFFS), DIMENSION(NLOOPFLOWS), SAVE ::
     $  LOOPCOLORPROJECTOR
      TYPE(PROJCOEFFS), DIMENSION(NBORNFLOWS), SAVE ::
     $  BORNCOLORPROJECTOR

C     ----------
C     BEGIN CODE
C     ----------

C     CleanUp duties, i.e. array deallocation
      IF (CLEANUP) THEN
        DO I=1,NLOOPFLOWS
          IF(ALLOCATED(LOOPCOLORPROJECTOR(I)%NUM))
     $      DEALLOCATE(LOOPCOLORPROJECTOR(I)%NUM)
          IF(ALLOCATED(LOOPCOLORPROJECTOR(I)%DENOM))
     $      DEALLOCATE(LOOPCOLORPROJECTOR(I)%DENOM)
          IF(ALLOCATED(LOOPCOLORPROJECTOR(I)%AMPID))
     $      DEALLOCATE(LOOPCOLORPROJECTOR(I)%AMPID)
        ENDDO
        DO I=1,NBORNFLOWS
          IF(ALLOCATED(BORNCOLORPROJECTOR(I)%NUM))
     $      DEALLOCATE(BORNCOLORPROJECTOR(I)%NUM)
          IF(ALLOCATED(BORNCOLORPROJECTOR(I)%DENOM))
     $      DEALLOCATE(BORNCOLORPROJECTOR(I)%DENOM)
          IF(ALLOCATED(BORNCOLORPROJECTOR(I)%AMPID))
     $      DEALLOCATE(BORNCOLORPROJECTOR(I)%AMPID)
        ENDDO
        RETURN
      ENDIF

C     Initialization; must allocate array from data files. All arrays
C      are allocated at once, so we only need to check if
C      LoopColorProjector(0)%Num is allocated and all other arrays
C      must share the same status.
      IF(.NOT.ALLOCATED(LOOPCOLORPROJECTOR(1)%NUM)) THEN
        CALL JOINPATH(MLPATH,PROC_PREFIX,TMP)
        CALL JOINPATH(TMP,LOOPCOLORFLOWCOEFSNAME,LOOPCOLORFLOWCOEFSN)
        CALL JOINPATH(TMP,BORNCOLORFLOWCOEFSNAME,BORNCOLORFLOWCOEFSN)

C       Initialize the LoopColorProjector
        OPEN(1, FILE=LOOPCOLORFLOWCOEFSN, ERR=201, STATUS='OLD',      
     $        ACTION='READ')
        DO I=1,NLOOPFLOWS
          READ(1,*,END=998) ARRAY_SIZE
          ALLOCATE(LOOPCOLORPROJECTOR(I)%NUM(ARRAY_SIZE))
          ALLOCATE(LOOPCOLORPROJECTOR(I)%DENOM(ARRAY_SIZE))
          ALLOCATE(LOOPCOLORPROJECTOR(I)%AMPID(ARRAY_SIZE))
          READ(1,*,END=998) (LOOPCOLORPROJECTOR(I)%NUM(J),J=1
     $     ,ARRAY_SIZE)
          READ(1,*,END=998) (LOOPCOLORPROJECTOR(I)%DENOM(J),J=1
     $     ,ARRAY_SIZE)
          READ(1,*,END=998) (LOOPCOLORPROJECTOR(I)%AMPID(J),J=1
     $     ,ARRAY_SIZE)
        ENDDO
        GOTO 203
 201    CONTINUE
        STOP 'Color projection coefficients could not be initialized'
     $   //' from file ML5_2_LoopColorFlowCoefs.dat.'
 203    CONTINUE
        CLOSE(1)

C       Initialize the BornColorProjector
        OPEN(1, FILE=BORNCOLORFLOWCOEFSN, ERR=301, STATUS='OLD',      
     $        ACTION='READ')
        DO I=1,NBORNFLOWS
          READ(1,*,END=998) ARRAY_SIZE
          ALLOCATE(BORNCOLORPROJECTOR(I)%NUM(ARRAY_SIZE))
          ALLOCATE(BORNCOLORPROJECTOR(I)%DENOM(ARRAY_SIZE))
          ALLOCATE(BORNCOLORPROJECTOR(I)%AMPID(ARRAY_SIZE))
          READ(1,*,END=998) (BORNCOLORPROJECTOR(I)%NUM(J),J=1
     $     ,ARRAY_SIZE)
          READ(1,*,END=998) (BORNCOLORPROJECTOR(I)%DENOM(J),J=1
     $     ,ARRAY_SIZE)
          READ(1,*,END=998) (BORNCOLORPROJECTOR(I)%AMPID(J),J=1
     $     ,ARRAY_SIZE)
        ENDDO
        GOTO 303
 301    CONTINUE
        STOP 'Color projection coefficients could not be initialized'
     $   //' from file ML5_2_BornColorFlowCoefs.dat.'
 303    CONTINUE
        CLOSE(1)

        GOTO 999
 998    CONTINUE
        STOP 'End of file reached. Should not have happened.'
 999    CONTINUE

      ENDIF

C     Start by projecting the Born amplitudes
      DO I=1,NBORNFLOWS
        DO J=1,SIZE(BORNCOLORPROJECTOR(I)%AMPID)
          SOINDEX = ML5_2_ML5SOINDEX_FOR_BORN_AMP(BORNCOLORPROJECTOR(I)
     $     %AMPID(J))
          PROJ_COEF=DCMPLX(BORNCOLORPROJECTOR(I)%NUM(J)
     $     /DBLE(ABS(BORNCOLORPROJECTOR(I)%DENOM(J))),0.0D0)
          IF(BORNCOLORPROJECTOR(I)%DENOM(J).LT.0) PROJ_COEF=PROJ_COEF
     $     *IMAG1
          JAMPB(I,SOINDEX) = JAMPB(I,SOINDEX) + PROJ_COEF
     $     *AMP(BORNCOLORPROJECTOR(I)%AMPID(J))
        ENDDO
      ENDDO

C     Projection of the loop amplitudes
      DO I=1,NLOOPFLOWS
        DO J=1,SIZE(LOOPCOLORPROJECTOR(I)%AMPID)
          SOINDEX = ML5_2_ML5SOINDEX_FOR_LOOP_AMP(LOOPCOLORPROJECTOR(I)
     $     %AMPID(J))
          PROJ_COEF=DCMPLX(LOOPCOLORPROJECTOR(I)%NUM(J)
     $     /DBLE(ABS(LOOPCOLORPROJECTOR(I)%DENOM(J))),0.0D0)
          IF(LOOPCOLORPROJECTOR(I)%DENOM(J).LT.0) PROJ_COEF=PROJ_COEF
     $     *IMAG1
          DO K=1,3
            JAMPL(K,I,SOINDEX) = JAMPL(K,I,SOINDEX) + PROJ_COEF*AMPL(K
     $       ,LOOPCOLORPROJECTOR(I)%AMPID(J))
          ENDDO
        ENDDO
      ENDDO

      END SUBROUTINE

C     This subroutine initializes the color flow matrix.
      SUBROUTINE ML5_2_INITIALIZE_FLOW_COLORMATRIX()
      IMPLICIT NONE
C     
C     CONSTANTS 
C     
      CHARACTER*512 PROC_PREFIX
      PARAMETER ( PROC_PREFIX='ML5_2_')
      CHARACTER*512 LOOPCOLORFLOWMATRIXNAME
      PARAMETER ( LOOPCOLORFLOWMATRIXNAME='LoopColorFlowMatrix.dat')
      CHARACTER*512 BORNCOLORFLOWMATRIXNAME
      PARAMETER ( BORNCOLORFLOWMATRIXNAME='BornColorFlowMatrix.dat')

      INTEGER NLOOPFLOWS
      PARAMETER (NLOOPFLOWS=3)
      INTEGER NBORNFLOWS
      PARAMETER (NBORNFLOWS=2)
C     
C     LOCAL VARIABLES 
C     
C     Storing concatenated filenames
      CHARACTER*512 TMP
      CHARACTER*512 LOOPCOLORFLOWMATRIXN
      CHARACTER*512 BORNCOLORFLOWMATRIXN
      INTEGER I, J
C     
C     GLOBAL VARIABLES
C     
      CHARACTER(512) MLPATH
      COMMON/MLPATH/MLPATH

C     Now a more advanced data structure for storing the projection
C      coefficient
C     This is of course not F77 standard but widely supported by now.
      TYPE COLORCOEFF
        SEQUENCE
        INTEGER :: NUM
        INTEGER :: DENOM
      ENDTYPE COLORCOEFF
      TYPE(COLORCOEFF), DIMENSION(NLOOPFLOWS,NBORNFLOWS) ::
     $  LOOPCOLORFLOWMATRIX
      TYPE(COLORCOEFF), DIMENSION(NBORNFLOWS,NBORNFLOWS) ::
     $  BORNCOLORFLOWMATRIX
      LOGICAL CMINITIALIZED
      DATA CMINITIALIZED/.FALSE./
      COMMON/ML5_2_FLOW_COLOR_MATRIX/LOOPCOLORFLOWMATRIX,
     $  BORNCOLORFLOWMATRIX, CMINITIALIZED

C     ----------
C     BEGIN CODE
C     ----------

C     Initialization
      IF(.NOT.CMINITIALIZED) THEN
        CMINITIALIZED = .TRUE.
        CALL JOINPATH(MLPATH,PROC_PREFIX,TMP)
        CALL JOINPATH(TMP,LOOPCOLORFLOWMATRIXNAME,LOOPCOLORFLOWMATRIXN)
        CALL JOINPATH(TMP,BORNCOLORFLOWMATRIXNAME,BORNCOLORFLOWMATRIXN)
C       Initialize the BornColorFlowMatrix
        OPEN(1, FILE=BORNCOLORFLOWMATRIXN, ERR=601, STATUS='OLD',     
     $         ACTION='READ')
        DO I=1,NBORNFLOWS
          READ(1,*,END=898) (BORNCOLORFLOWMATRIX(I,J)%NUM,J=1
     $     ,NBORNFLOWS)
          READ(1,*,END=898) (BORNCOLORFLOWMATRIX(I,J)%DENOM,J=1
     $     ,NBORNFLOWS)
        ENDDO
        GOTO 603
 601    CONTINUE
        STOP 'Color factors could not be initialized from file'
     $   //' ML5_2_BornColorFlowMatrix.dat.'
 603    CONTINUE
        CLOSE(1)
C       Initialize the LoopColorFlowMatrix
        OPEN(1, FILE=LOOPCOLORFLOWMATRIXN, ERR=701, STATUS='OLD',     
     $         ACTION='READ')
        DO I=1,NLOOPFLOWS
          READ(1,*,END=898) (LOOPCOLORFLOWMATRIX(I,J)%NUM,J=1
     $     ,NBORNFLOWS)
          READ(1,*,END=898) (LOOPCOLORFLOWMATRIX(I,J)%DENOM,J=1
     $     ,NBORNFLOWS)
        ENDDO
        GOTO 703
 701    CONTINUE
        STOP 'Color factors could not be initialized from file'
     $   //' ML5_2_BornColorFlowMatrix.dat.'
 703    CONTINUE
        CLOSE(1)

        GOTO 899
 898    CONTINUE
        STOP 'End of file reached. Should not have happened.'
 899    CONTINUE

      ENDIF
      END SUBROUTINE


C     This routine is used as a crosscheck only to make sure the loop
C      ME computed
C     from the JAMP is equal to the amplitude computed directly. It is
C      a consistency
C     check of the color projection and computation.
      SUBROUTINE ML5_2_COMPUTE_RES_FROM_JAMP(RES,HEL_MULT)
      IMPLICIT NONE
C     
C     CONSTANTS 
C     
      COMPLEX*16 IMAG1
      PARAMETER (IMAG1=(0.0D0,1.0D0))

      INTEGER    NSQUAREDSO, NLOOPAMPSO, NSO
      PARAMETER (NSQUAREDSO=1, NLOOPAMPSO=2, NSO=1)
      INTEGER NLOOPFLOWS
      PARAMETER (NLOOPFLOWS=3)
      INTEGER NBORNFLOWS
      PARAMETER (NBORNFLOWS=2)
C     We use the same list of amplitude split order for Born and Loop,
C      so that NBORNAMPSO is always equal to NLOOPAMPSO  
      INTEGER    NBORNAMPSO
      PARAMETER (NBORNAMPSO=NLOOPAMPSO)
C     
C     ARGUMENT 
C     
      REAL*8 RES(0:3,0:NSQUAREDSO)
C     HEL_MULT is the helicity multiplier which can be more than one
C      if several helicity configuration are mapped onto the one being
C      currently computed.
      INTEGER HEL_MULT
C     
C     LOCAL VARIABLES 
C     
      INTEGER I, J, M, N, K, ISQSO, DOUBLEFACT
      INTEGER ORDERS_A(NSO), ORDERS_B(NSO)
      COMPLEX*16 COLOR_COEF
      REAL*8 TEMP(3)
C     
C     FUNCTIONS
C     
C     This function belongs to the loop ME computation (loop_matrix.f)
C      and is prefixed with ML5 
      INTEGER ML5_2_ML5SQSOINDEX
C     This function belongs to the Born ME computation (born_matrix.f)
C      and is not prefixed at all
      INTEGER ML5_2_SQSOINDEX, ML5_2_SOINDEX_FOR_AMPORDERS
C     
C     GLOBAL VARIABLES
C     
      CHARACTER(512) MLPATH
      COMMON/MLPATH/MLPATH
      INTEGER SQSO_TARGET
      COMMON/ML5_2_SOCHOICE/SQSO_TARGET

      COMPLEX*16 JAMPL(3,NLOOPFLOWS,NLOOPAMPSO)
      COMMON/ML5_2_JAMPL/JAMPL
      COMPLEX*16 JAMPB(NBORNFLOWS,NBORNAMPSO)
      COMMON/ML5_2_JAMPB/JAMPB

      LOGICAL UVCT_REQ_SO_DONE,MP_UVCT_REQ_SO_DONE,CT_REQ_SO_DONE
     $ ,MP_CT_REQ_SO_DONE,LOOP_REQ_SO_DONE,MP_LOOP_REQ_SO_DONE
     $ ,CTCALL_REQ_SO_DONE,FILTER_SO
      COMMON/ML5_2_SO_REQS/UVCT_REQ_SO_DONE,MP_UVCT_REQ_SO_DONE
     $ ,CT_REQ_SO_DONE,MP_CT_REQ_SO_DONE,LOOP_REQ_SO_DONE
     $ ,MP_LOOP_REQ_SO_DONE,CTCALL_REQ_SO_DONE,FILTER_SO

C     Now a more advanced data structure for storing the projection
C      coefficient
C     This is of course not F77 standard but widely supported by now.
      TYPE COLORCOEFF
        SEQUENCE
        INTEGER :: NUM
        INTEGER :: DENOM
      ENDTYPE COLORCOEFF
      TYPE(COLORCOEFF), DIMENSION(NLOOPFLOWS,NBORNFLOWS) ::
     $  LOOPCOLORFLOWMATRIX
      TYPE(COLORCOEFF), DIMENSION(NBORNFLOWS,NBORNFLOWS) ::
     $  BORNCOLORFLOWMATRIX
      LOGICAL CMINITIALIZED
      COMMON/ML5_2_FLOW_COLOR_MATRIX/LOOPCOLORFLOWMATRIX,
     $  BORNCOLORFLOWMATRIX, CMINITIALIZED

C     ----------
C     BEGIN CODE
C     ----------

C     Initialization
      IF(.NOT.CMINITIALIZED) THEN
        CALL ML5_2_INITIALIZE_FLOW_COLORMATRIX()
      ENDIF

      DO I=0,NSQUAREDSO
        DO K=0,3
          RES(K,I)=0.0D0
        ENDDO
      ENDDO

C     Compute the Born ME from the Born color flow amplitudes (JAMP)
      DO I=1, NBORNFLOWS
        DO J=I, NBORNFLOWS
          COLOR_COEF=DCMPLX(BORNCOLORFLOWMATRIX(I,J)%NUM
     $     /DBLE(ABS(BORNCOLORFLOWMATRIX(I,J)%DENOM)),0.0D0)
          IF (BORNCOLORFLOWMATRIX(I,J)%DENOM.LT.0)
     $      COLOR_COEF=COLOR_COEF*IMAG1
          DO M=1,NBORNAMPSO
C           It may be that this AmpSO index does not receive
C            contribution by the Born amps (because we put the loop
C            and Born amplitude split orders in a common list)
            IF (ABS(JAMPB(I,M)).EQ.0.0D0) CYCLE
            DO N=1,NBORNAMPSO
              IF (ABS(JAMPB(J,N)).EQ.0.0D0) CYCLE
C             First fetch what orders the split order indices M, N
C              correspond to
              CALL ML5_2_ML5GET_ORDERS_FOR_AMPSOINDEX(M,ORDERS_A)
              CALL ML5_2_ML5GET_ORDERS_FOR_AMPSOINDEX(N,ORDERS_B)
C             Now figure out to which SQSOINDEX these orders together
C              correspond to *in the Born subroutine* (Notice the
C              absence of the ML5 prefix in the functions used then)
              ISQSO =
     $          ML5_2_SQSOINDEX(ML5_2_SOINDEX_FOR_AMPORDERS(ORDERS_A)
     $         ,ML5_2_SOINDEX_FOR_AMPORDERS(ORDERS_B))
              IF(J.NE.I) THEN
                DOUBLEFACT=2
              ELSE
                DOUBLEFACT=1
              ENDIF
              TEMP(1) = DOUBLEFACT*HEL_MULT*DBLE(COLOR_COEF*JAMPB(I,M)
     $         *DCONJG(JAMPB(J,N)))
              RES(0,ISQSO) = RES(0,ISQSO) + TEMP(1)
              IF((.NOT.FILTER_SO).OR.SQSO_TARGET.EQ.
     $         -1.OR.SQSO_TARGET.EQ.ISQSO) THEN
                RES(0,0) = RES(0,0) + TEMP(1)
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      ENDDO

C     Compute the Loop ME from the loop color flow amplitudes (JAMPL)
      DO I=1,NLOOPFLOWS
        DO J=1,NBORNFLOWS
          COLOR_COEF=DCMPLX(LOOPCOLORFLOWMATRIX(I,J)%NUM
     $     /DBLE(ABS(LOOPCOLORFLOWMATRIX(I,J)%DENOM)),0.0D0)
          IF (LOOPCOLORFLOWMATRIX(I,J)%DENOM.LT.0)
     $      COLOR_COEF=COLOR_COEF*IMAG1
          DO M=1,NLOOPAMPSO
C           It may be that this AmpSO index does not receive
C            contribution by the Loop amps (because we put the loop
C            and Born amplitude split orders in a common list)
            IF((ABS(JAMPL(1,I,M))+ABS(JAMPL(2,I,M))+ABS(JAMPL(3,I,M)))
     $       .EQ.0.0D0) CYCLE
            DO N=1,NBORNAMPSO
C             Same for contributions of split order index N by the
C              Born amps
              IF(ABS(JAMPB(J,N)).EQ.0.0D0) CYCLE
              ISQSO = ML5_2_ML5SQSOINDEX(M,N)
              DO K=1,3
                TEMP(K) = 2.0D0*HEL_MULT*DBLE(COLOR_COEF*JAMPL(K,I,M)
     $           *DCONJG(JAMPB(J,N)))
                RES(K,ISQSO) = RES(K,ISQSO) + TEMP(K)
              ENDDO
              IF((.NOT.FILTER_SO).OR.SQSO_TARGET.EQ.
     $         -1.OR.SQSO_TARGET.EQ.ISQSO) THEN
                DO K=1,3
                  RES(K,0) = RES(K,0) + TEMP(K)
                ENDDO
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      ENDDO

      END SUBROUTINE


C     This subroutine resets to 0 the common arrays JAMPL, JAMPB and
C      possibly JAMPL_FOR_AMP2 if used
      SUBROUTINE ML5_2_REINITIALIZE_JAMPS()
      IMPLICIT NONE
C     
C     CONSTANTS 
C     
      COMPLEX*16 CMPLXZERO
      PARAMETER (CMPLXZERO=(0.0D0,0.0D0))
      REAL*8 ZERO
      PARAMETER (ZERO=0.0D0)
      INTEGER    NLOOPAMPSO
      PARAMETER (NLOOPAMPSO=2)
      INTEGER NLOOPFLOWS
      PARAMETER (NLOOPFLOWS=3)
      INTEGER NBORNFLOWS
      PARAMETER (NBORNFLOWS=2)
C     We use the same list of amplitude split order for Born and Loop,
C      so that NBORNAMPSO is always equal to NLOOPAMPSO
      INTEGER    NBORNAMPSO
      PARAMETER (NBORNAMPSO=NLOOPAMPSO)
C     
C     LOCAL VARIABLES
C     
      INTEGER I,J,K,L
C     
C     GLOBAL VARIABLES
C     
      COMPLEX*16 JAMPL(3,NLOOPFLOWS,NLOOPAMPSO)
      COMMON/ML5_2_JAMPL/JAMPL
      COMPLEX*16 JAMPB(NBORNFLOWS,NBORNAMPSO)
      COMMON/ML5_2_JAMPB/JAMPB

C     ----------
C     BEGIN CODE
C     ----------

      DO I=1,NLOOPAMPSO
        DO J=1,NLOOPFLOWS
          DO K=1,3
            JAMPL(K,J,I)=CMPLXZERO
          ENDDO
        ENDDO
      ENDDO
      DO I=1,NBORNAMPSO
        DO J=1,NBORNFLOWS
          JAMPB(J,I)=CMPLXZERO
        ENDDO
      ENDDO

      END SUBROUTINE


      ! HSS DEBUG
      ! get <Ci|A(0)> and <Ci|A(1)>
      SUBROUTINE ML5_2_DO_COMPUTE_NEW_COLOR_FLOWS(DUMMY)
      IMPLICIT NONE
      LOGICAL DUMMY
      INTEGER I,J
      INCLUDE 'NCOLORBASIS.inc'
      INTEGER NLOOPAMPSO
      PARAMETER (NLOOPAMPSO=2)
      INTEGER NLOOPFLOWS
      PARAMETER (NLOOPFLOWS=3)
      INTEGER NBORNFLOWS
      PARAMETER (NBORNFLOWS=2)
      INTEGER NBORNAMPSO
      PARAMETER (NBORNAMPSO=NLOOPAMPSO)
      COMPLEX*16 ZERO
      PARAMETER (ZERO=DCMPLX(0d0,0d0))
      COMPLEX*16 JAMPL(3,NLOOPFLOWS,NLOOPAMPSO)
      COMMON/ML5_2_JAMPL/JAMPL
      COMPLEX*16 JAMPB(NBORNFLOWS,NBORNAMPSO)
      COMMON/ML5_2_JAMPB/JAMPB
      COMPLEX*16 JAMPL_OCB(3,NCOLORBASIS)
      COMMON/ML5_2_JAMPL_OCB/JAMPL_OCB
      COMPLEX*16 JAMPB_OCB(NCOLORBASIS)
      COMMON/ML5_2_JAMPB_OCB/JAMPB_OCB
      INCLUDE 'COLORROTATE.inc'

      ! check JAMPB, whether some of them are zero
      DO I=1,NCOLORBASIS
         JAMPB_OCB(I)=ZERO
         JAMPL_OCB(1,I)=ZERO
         JAMPL_OCB(2,I)=ZERO
         JAMPL_OCB(3,I)=ZERO
         DO J=1,NCOLORBASIS
            IF(J.LT.NCOLORBASIS)THEN
            ! WARNING WARNING WARNING WARNING WARNING
            ! check whether the Born JAMPB is zero for some color basis
            ! WARNING WARNING WARNING WARNING WARNING
            JAMPB_OCB(I)=JAMPB_OCB(I)+DBLE(ColorRotateNum(I,J))
     $           /DBLE(ColorRotateDen(I,J))*JAMPB(J,1)
            ENDIF
            JAMPL_OCB(1,I)=JAMPL_OCB(1,I)+DBLE(ColorRotateNum(I,J))
     $           /DBLE(ColorRotateDen(I,J))*JAMPL(1,J,1)
            JAMPL_OCB(2,I)=JAMPL_OCB(2,I)+DBLE(ColorRotateNum(I,J))
     $           /DBLE(ColorRotateDen(I,J))*JAMPL(2,J,1)
            JAMPL_OCB(3,I)=JAMPL_OCB(3,I)+DBLE(ColorRotateNum(I,J))
     $           /DBLE(ColorRotateDen(I,J))*JAMPL(3,J,1)
         ENDDO
      ENDDO
      RETURN
      END SUBROUTINE


      ! HSS DEBUG
      ! Born:     <ci|A(0)><A(0)|cj>/(<ci|ci><cj|cj>)
      ! One loop: (<ci|A(0)><A(1)|cj>+<ci|A(1)><A(0)|cj>)/(<ci|ci><cj|cj>)
      SUBROUTINE ML5_2_COMPUTE_HARDFUNCTION(RES,HEL_MULT)
      IMPLICIT NONE
      INCLUDE 'NCOLORBASIS.inc'
      COMPLEX*16 RES(0:3,NCOLORBASIS,NCOLORBASIS)
      INTEGER HEL_MULT
      INCLUDE 'COLORNORM.inc'
      COMPLEX*16 JAMPL_OCB(3,NCOLORBASIS)
      COMMON/ML5_2_JAMPL_OCB/JAMPL_OCB
      COMPLEX*16 JAMPB_OCB(NCOLORBASIS)
      COMMON/ML5_2_JAMPB_OCB/JAMPB_OCB
      INTEGER I,J,K
      REAL*8 TEMP
      COMPLEX*16 ZERO
      PARAMETER (ZERO=DCMPLX(0d0,0d0))

      DO I=1,NCOLORBASIS
         DO J=1,NCOLORBASIS
            DO K=0,3
               RES(K,I,J)=ZERO
            ENDDO
         ENDDO
      ENDDO

C     Compute hard function
      DO I=1,NCOLORBASIS
         DO J=1,NCOLORBASIS
            ! helicity multiply factor/<ci|ci>/<cj|cj>
            TEMP=DBLE(HEL_MULT)/(DBLE(CiCiNum(I))*DBLE(CiCiNum(J)))
            TEMP=TEMP*(DBLE(CiCiDen(I))*DBLE(CiCiDen(J)))
            RES(0,I,J)=TEMP*JAMPB_OCB(I)*DCONJG(JAMPB_OCB(J))
            RES(1,I,J)=JAMPB_OCB(I)*DCONJG(JAMPL_OCB(1,J))
     $           +JAMPL_OCB(1,I)*DCONJG(JAMPB_OCB(J))
            RES(1,I,J)=TEMP*RES(1,I,J)
            RES(2,I,J)=JAMPB_OCB(I)*DCONJG(JAMPL_OCB(2,J))
     $           +JAMPL_OCB(2,I)*DCONJG(JAMPB_OCB(J))
            RES(2,I,J)=TEMP*RES(2,I,J)
            RES(3,I,J)=JAMPB_OCB(I)*DCONJG(JAMPL_OCB(3,J))
     $           +JAMPL_OCB(3,I)*DCONJG(JAMPB_OCB(J))
            RES(3,I,J)=TEMP*RES(3,I,J)
         ENDDO
      ENDDO
      RETURN
      END SUBROUTINE



C     This subroutine resets all cumulative (for each helicity) arrays
C      in common block and deriving from the color flow computation.
C     This subroutine is called by loop_matrix.f at appropriate time
C      when these arrays must be reset.
C     An example of this is the JAMP2 which must be cumulatively
C      computed as loop_matrix.f loops over helicity configs but must
C      be resets when loop_matrix.f starts over with helicity one (for
C      the stability test for example).
      SUBROUTINE ML5_2_REINITIALIZE_CUMULATIVE_ARRAYS()
      IMPLICIT NONE
C     
C     CONSTANTS 
C     
      COMPLEX*16 CMPLXZERO
      PARAMETER (CMPLXZERO=(0.0D0,0.0D0))
      REAL*8 ZERO
      PARAMETER (ZERO=0.0D0)
      INTEGER    NLOOPDIAGRAMS
      PARAMETER (NLOOPDIAGRAMS=42)
      INTEGER NLOOPFLOWS
      PARAMETER (NLOOPFLOWS=3)
C     
C     LOCAL VARIABLES
C     
      INTEGER I,J,K
C     
C     GLOBAL VARIABLES
C     

C     ----------
C     BEGIN CODE
C     ----------


      CONTINUE

      END SUBROUTINE


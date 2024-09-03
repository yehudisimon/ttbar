      SUBROUTINE ML5_0_MP_COEF_CONSTRUCTION_1(P,NHEL,H,IC)
C     
      USE ML5_0_POLYNOMIAL_CONSTANTS
      IMPLICIT NONE
C     
C     CONSTANTS
C     
      INTEGER    NEXTERNAL
      PARAMETER (NEXTERNAL=4)
      INTEGER    NCOMB
      PARAMETER (NCOMB=16)

      INTEGER NBORNAMPS
      PARAMETER (NBORNAMPS=1)
      INTEGER    NLOOPS, NLOOPGROUPS, NCTAMPS
      PARAMETER (NLOOPS=10, NLOOPGROUPS=10, NCTAMPS=27)
      INTEGER    NLOOPAMPS
      PARAMETER (NLOOPAMPS=37)
      INTEGER    NWAVEFUNCS,NLOOPWAVEFUNCS
      PARAMETER (NWAVEFUNCS=6,NLOOPWAVEFUNCS=24)
      REAL*16     ZERO
      PARAMETER (ZERO=0.0E0_16)
      COMPLEX*32     IZERO
      PARAMETER (IZERO=CMPLX(0.0E0_16,0.0E0_16,KIND=16))
C     These are constants related to the split orders
      INTEGER    NSO, NSQUAREDSO, NAMPSO
      PARAMETER (NSO=1, NSQUAREDSO=1, NAMPSO=2)
C     
C     ARGUMENTS
C     
      REAL*16 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL), IC(NEXTERNAL)
      INTEGER H
C     
C     LOCAL VARIABLES
C     
      INTEGER I,J,K
      COMPLEX*32 COEFS(MAXLWFSIZE,0:VERTEXMAXCOEFS-1,MAXLWFSIZE)
C     
C     GLOBAL VARIABLES
C     
      INCLUDE 'mp_coupl_same_name.inc'

      INTEGER GOODHEL(NCOMB)
      LOGICAL GOODAMP(NSQUAREDSO,NLOOPGROUPS)
      COMMON/ML5_0_FILTERS/GOODAMP,GOODHEL

      INTEGER SQSO_TARGET
      COMMON/ML5_0_SOCHOICE/SQSO_TARGET

      LOGICAL UVCT_REQ_SO_DONE,MP_UVCT_REQ_SO_DONE,CT_REQ_SO_DONE
     $ ,MP_CT_REQ_SO_DONE,LOOP_REQ_SO_DONE,MP_LOOP_REQ_SO_DONE
     $ ,CTCALL_REQ_SO_DONE,FILTER_SO
      COMMON/ML5_0_SO_REQS/UVCT_REQ_SO_DONE,MP_UVCT_REQ_SO_DONE
     $ ,CT_REQ_SO_DONE,MP_CT_REQ_SO_DONE,LOOP_REQ_SO_DONE
     $ ,MP_LOOP_REQ_SO_DONE,CTCALL_REQ_SO_DONE,FILTER_SO

      COMPLEX*32 AMP(NBORNAMPS)
      COMMON/ML5_0_MP_AMPS/AMP
      COMPLEX*32 W(20,NWAVEFUNCS)
      COMMON/ML5_0_MP_W/W

      COMPLEX*32 WL(MAXLWFSIZE,0:LOOPMAXCOEFS-1,MAXLWFSIZE,
     $ -1:NLOOPWAVEFUNCS)
      COMPLEX*32 PL(0:3,-1:NLOOPWAVEFUNCS)
      COMMON/ML5_0_MP_WL/WL,PL

      COMPLEX*32 AMPL(3,NLOOPAMPS)
      COMMON/ML5_0_MP_AMPL/AMPL

C     
C     ----------
C     BEGIN CODE
C     ----------

C     The target squared split order contribution is already reached
C      if true.
      IF (FILTER_SO.AND.MP_LOOP_REQ_SO_DONE) THEN
        GOTO 1001
      ENDIF

C     Coefficient construction for loop diagram with ID 2
      CALL MP_VVV1L2P0_1(PL(0,0),W(1,5),GC_4,ZERO,ZERO,PL(0,1),COEFS)
      CALL MP_ML5_0_UPDATE_WL_0_1(WL(1,0,1,0),4,COEFS,4,4,WL(1,0,1,1))
      CALL MP_VVV1L2P0_1(PL(0,1),W(1,6),GC_4,ZERO,ZERO,PL(0,2),COEFS)
      CALL MP_ML5_0_UPDATE_WL_1_1(WL(1,0,1,1),4,COEFS,4,4,WL(1,0,1,2))
      CALL MP_ML5_0_CREATE_LOOP_COEFS(WL(1,0,1,2),2,4,1,2,1,28,H)
C     Coefficient construction for loop diagram with ID 3
      CALL MP_FFV1L3_1(PL(0,0),W(1,3),GC_5,MDL_MT,MDL_WT,PL(0,3),COEFS)
      CALL MP_ML5_0_UPDATE_WL_0_1(WL(1,0,1,0),4,COEFS,4,4,WL(1,0,1,3))
      CALL MP_FFV1L2P0_3(PL(0,3),W(1,4),GC_5,ZERO,ZERO,PL(0,4),COEFS)
      CALL MP_ML5_0_UPDATE_WL_1_0(WL(1,0,1,3),4,COEFS,4,4,WL(1,0,1,4))
      CALL MP_VVV1L2P0_1(PL(0,4),W(1,5),GC_4,ZERO,ZERO,PL(0,5),COEFS)
      CALL MP_ML5_0_UPDATE_WL_1_1(WL(1,0,1,4),4,COEFS,4,4,WL(1,0,1,5))
      CALL MP_ML5_0_CREATE_LOOP_COEFS(WL(1,0,1,5),2,4,2,1,1,29,H)
C     Coefficient construction for loop diagram with ID 4
      CALL MP_FFV1L1P0_3(PL(0,0),W(1,3),GC_5,ZERO,ZERO,PL(0,6),COEFS)
      CALL MP_ML5_0_UPDATE_WL_0_0(WL(1,0,1,0),4,COEFS,4,4,WL(1,0,1,6))
      CALL MP_FFV1L3_2(PL(0,6),W(1,4),GC_5,MDL_MT,MDL_WT,PL(0,7),COEFS)
      CALL MP_ML5_0_UPDATE_WL_0_1(WL(1,0,1,6),4,COEFS,4,4,WL(1,0,1,7))
      CALL MP_FFV1L1_2(PL(0,7),W(1,5),GC_5,MDL_MT,MDL_WT,PL(0,8),COEFS)
      CALL MP_ML5_0_UPDATE_WL_1_1(WL(1,0,1,7),4,COEFS,4,4,WL(1,0,1,8))
      CALL MP_ML5_0_CREATE_LOOP_COEFS(WL(1,0,1,8),2,4,3,1,1,30,H)
C     Coefficient construction for loop diagram with ID 5
      CALL MP_FFV1L2P0_3(PL(0,0),W(1,1),GC_5,ZERO,ZERO,PL(0,9),COEFS)
      CALL MP_ML5_0_UPDATE_WL_0_0(WL(1,0,1,0),4,COEFS,4,4,WL(1,0,1,9))
      CALL MP_FFV1L3_1(PL(0,9),W(1,2),GC_5,ZERO,ZERO,PL(0,10),COEFS)
      CALL MP_ML5_0_UPDATE_WL_0_1(WL(1,0,1,9),4,COEFS,4,4,WL(1,0,1,10))
      CALL MP_FFV1L2_1(PL(0,10),W(1,6),GC_5,ZERO,ZERO,PL(0,11),COEFS)
      CALL MP_ML5_0_UPDATE_WL_1_1(WL(1,0,1,10),4,COEFS,4,4,WL(1,0,1,11)
     $ )
      CALL MP_ML5_0_CREATE_LOOP_COEFS(WL(1,0,1,11),2,4,4,1,1,31,H)
C     Coefficient construction for loop diagram with ID 6
      CALL MP_FFV1L3_2(PL(0,0),W(1,1),GC_5,ZERO,ZERO,PL(0,12),COEFS)
      CALL MP_ML5_0_UPDATE_WL_0_1(WL(1,0,1,0),4,COEFS,4,4,WL(1,0,1,12))
      CALL MP_FFV1L1P0_3(PL(0,12),W(1,2),GC_5,ZERO,ZERO,PL(0,13),COEFS)
      CALL MP_ML5_0_UPDATE_WL_1_0(WL(1,0,1,12),4,COEFS,4,4,WL(1,0,1,13)
     $ )
      CALL MP_VVV1L2P0_1(PL(0,13),W(1,6),GC_4,ZERO,ZERO,PL(0,14),COEFS)
      CALL MP_ML5_0_UPDATE_WL_1_1(WL(1,0,1,13),4,COEFS,4,4,WL(1,0,1,14)
     $ )
      CALL MP_ML5_0_CREATE_LOOP_COEFS(WL(1,0,1,14),2,4,5,1,1,32,H)
C     Coefficient construction for loop diagram with ID 7
      CALL MP_FFV1L3_2(PL(0,13),W(1,4),GC_5,MDL_MT,MDL_WT,PL(0,15)
     $ ,COEFS)
      CALL MP_ML5_0_UPDATE_WL_1_1(WL(1,0,1,13),4,COEFS,4,4,WL(1,0,1,15)
     $ )
      CALL MP_FFV1L1P0_3(PL(0,15),W(1,3),GC_5,ZERO,ZERO,PL(0,16),COEFS)
      CALL MP_ML5_0_UPDATE_WL_2_0(WL(1,0,1,15),4,COEFS,4,4,WL(1,0,1,16)
     $ )
      CALL MP_ML5_0_CREATE_LOOP_COEFS(WL(1,0,1,16),2,4,6,1,1,33,H)
C     Coefficient construction for loop diagram with ID 8
      CALL MP_FFV1L3_1(PL(0,13),W(1,3),GC_5,MDL_MT,MDL_WT,PL(0,17)
     $ ,COEFS)
      CALL MP_ML5_0_UPDATE_WL_1_1(WL(1,0,1,13),4,COEFS,4,4,WL(1,0,1,17)
     $ )
      CALL MP_FFV1L2P0_3(PL(0,17),W(1,4),GC_5,ZERO,ZERO,PL(0,18),COEFS)
      CALL MP_ML5_0_UPDATE_WL_2_0(WL(1,0,1,17),4,COEFS,4,4,WL(1,0,1,18)
     $ )
      CALL MP_ML5_0_CREATE_LOOP_COEFS(WL(1,0,1,18),2,4,7,1,1,34,H)
C     Coefficient construction for loop diagram with ID 9
      CALL MP_GHGHGL2_1(PL(0,0),W(1,5),GC_4,ZERO,ZERO,PL(0,19),COEFS)
      CALL MP_ML5_0_UPDATE_WL_0_1(WL(1,0,1,0),1,COEFS,1,1,WL(1,0,1,19))
      CALL MP_GHGHGL2_1(PL(0,19),W(1,6),GC_4,ZERO,ZERO,PL(0,20),COEFS)
      CALL MP_ML5_0_UPDATE_WL_1_1(WL(1,0,1,19),1,COEFS,1,1,WL(1,0,1,20)
     $ )
      CALL MP_ML5_0_CREATE_LOOP_COEFS(WL(1,0,1,20),2,1,8,1,1,35,H)
C     Coefficient construction for loop diagram with ID 10
      CALL MP_FFV1L2_1(PL(0,0),W(1,5),GC_5,ZERO,ZERO,PL(0,21),COEFS)
      CALL MP_ML5_0_UPDATE_WL_0_1(WL(1,0,1,0),4,COEFS,4,4,WL(1,0,1,21))
      CALL MP_FFV1L2_1(PL(0,21),W(1,6),GC_5,ZERO,ZERO,PL(0,22),COEFS)
      CALL MP_ML5_0_UPDATE_WL_1_1(WL(1,0,1,21),4,COEFS,4,4,WL(1,0,1,22)
     $ )
      CALL MP_ML5_0_CREATE_LOOP_COEFS(WL(1,0,1,22),2,4,9,1,5,36,H)
C     Coefficient construction for loop diagram with ID 11
      CALL MP_FFV1L2_1(PL(0,0),W(1,5),GC_5,MDL_MT,MDL_WT,PL(0,23)
     $ ,COEFS)
      CALL MP_ML5_0_UPDATE_WL_0_1(WL(1,0,1,0),4,COEFS,4,4,WL(1,0,1,23))
      CALL MP_FFV1L2_1(PL(0,23),W(1,6),GC_5,MDL_MT,MDL_WT,PL(0,24)
     $ ,COEFS)
      CALL MP_ML5_0_UPDATE_WL_1_1(WL(1,0,1,23),4,COEFS,4,4,WL(1,0,1,24)
     $ )
      CALL MP_ML5_0_CREATE_LOOP_COEFS(WL(1,0,1,24),2,4,10,1,1,37,H)
C     At this point, all loop coefficients needed for (QCD=6), i.e. of
C      split order ID=1, are computed.
      IF(FILTER_SO.AND.SQSO_TARGET.EQ.1) GOTO 4000

      GOTO 1001
 4000 CONTINUE
      MP_LOOP_REQ_SO_DONE=.TRUE.
 1001 CONTINUE
      END


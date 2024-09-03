C     This File is Automatically generated by ALOHA 
C     The process calculated in this file is: 
C     P(3,2)
C     
      SUBROUTINE GHGHGL1_2(P1, V3, COUP, M2, W2, P2, COEFF)
      IMPLICIT NONE
      COMPLEX*16 CI
      PARAMETER (CI=(0D0,1D0))
      INCLUDE 'coef_specs.inc'
      COMPLEX*16 COEFF(MAXLWFSIZE,0:VERTEXMAXCOEFS-1,MAXLWFSIZE)
      REAL*8 W2
      COMPLEX*16 V3(*)
      COMPLEX*16 P2(0:3)
      COMPLEX*16 P1(0:3)
      COMPLEX*16 COUP
      COMPLEX*16 TMP1
      REAL*8 M2
      P2(0) = +P1(0)+V3(1)
      P2(1) = +P1(1)+V3(2)
      P2(2) = +P1(2)+V3(3)
      P2(3) = +P1(3)+V3(4)
      TMP1 = (V3(5)*P2(0)-V3(6)*P2(1)-V3(7)*P2(2)-V3(8)*P2(3))
      COEFF(1,0,1)= COUP*-CI * TMP1
      COEFF(1,1,1)= COUP*-CI * V3(5)
      COEFF(1,2,1)= COUP*CI * V3(6)
      COEFF(1,3,1)= COUP*CI * V3(7)
      COEFF(1,4,1)= COUP*CI * V3(8)
      END



PROGRAM test
  USE genps_nbody
  IMPLICIT NONE
  INTEGER::n
  REAL(KIND(1d0))::Q
  REAL(KIND(1d0)),DIMENSION(3)::mass
  REAL(KIND(1d0)),DIMENSION(6)::x
  REAL(KIND(1d0)),DIMENSION(5,0:3)::p
  REAL(KIND(1d0))::jac
  INTEGER::i
  n=2
  mass(1)=173.2d0
  mass(2)=173.2d0
  Q=1450d0
  x(1)=0.3d0
  x(2)=0.1d0
  x(3)=0.7d0
  CALL generate_mom(n,Q,mass,x,p,jac)
  DO i=1,n+2
     PRINT *, p(i,0:3)
  ENDDO
  PRINT *, "jac=",jac
  RETURN
END PROGRAM TEST

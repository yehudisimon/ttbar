      PROGRAM DRIVER
      IMPLICIT NONE
      double complex res(0:2)
      double precision s,t,mt2,mu2
      integer i,j
      mt2=173.3d0**2
      mu2=1d4
      s=12d0*mt2
      t=-0.3d0*mt2
      do i=1,3
         do j=1,3
            CALL gg2ttx_sHFmatrix(i,j,s,t,mt2,mu2,res)
c     PRINT *, "virtual weights:"
            PRINT *, i,j,res(0:2)
         enddo
      enddo
      END

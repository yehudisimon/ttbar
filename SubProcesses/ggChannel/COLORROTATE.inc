      ! rotate color basis
      INTEGER ColorRotateNum(NCOLORBASIS,NCOLORBASIS)
      INTEGER ColorRotateDen(NCOLORBASIS,NCOLORBASIS)
      DATA (ColorRotateNum(I,1),I=1,3) /4,6,10/
      DATA (ColorRotateDen(I,1),I=1,3) /1,1,3/
      DATA (ColorRotateNum(I,2),I=1,3) /4,-6,10/
      DATA (ColorRotateDen(I,2),I=1,3) /1,1,3/
      DATA (ColorRotateNum(I,3),I=1,3) /12,0,0/
      DATA (ColorRotateDen(I,3),I=1,3) /1,1,1/

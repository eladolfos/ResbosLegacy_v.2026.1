
C **************************************************************************    
C         
C  MODULE: 3GAUSS                                   Version 3.0
C  --------------       
C
C           adapted from CERN library by Peter Arnold
C         
C    THIS MODULE CONTAINS ROUTINES FOR ADAPTIVE AND NON-ADAPTIVE GAUSSIAN
C  INTEGRATION.  THEY INTEGRATE F FROM A TO B.
C
C    AGAUSS(F,A,B,N,ERR)
C       IF N > 0, WILL DIVIDE [A,B] INTO N INTERVALS AND PERFORM 16-POINT
C       GAUSSIAN INTEGRATION IN EACH.  IN THIS CASE, ERR IS IGNORED.
C       IF N = 0, CALLS AAGAUS INSTEAD.
C
C    AAGAUS(F,A,B,ERR)
C       ADAPTIVE GAUSSIAN INTEGRATION, REQURING RESULT TO BE WITHIN
C       RELATIVE ERROR SPECIFIED BY ERR.
C
C  THESE ROUTINES HAVE BEEN ADAPTED FROM THE CERN LIBRARY.
C  BGAUSS, BBGAUS, XGAUSS, AND XXGAUS ARE COPIES PROVIDED FOR 
C  3-DIMENSIONAL INTEGRATION.   
C
C    GAUSS8(F,A,B,N,ERR)
C      PERFORMS 8-PT. INTEGRATION IN N INTERVALS.  ERR IS IGNORED
C      AND N MUST BE POSITIVE.
C         
C **************************************************************************    
          


C **  THE FIRST 4 ENTRIES OF W AND X ARE THE WEIGHTS AND COORDINATES FOR
C **  8-POINT GAUSSIAN INTEGRATION.  THE REMAINING 8 ENTRIES ARE THOSE
C **  FOR 16-POINT GAUSSIAN INTEGRATION.
      
      BLOCK DATA GDAT

      IMPLICIT NONE
 
      DOUBLE PRECISION W,X
      COMMON /GAUDAT/ W(12),X(12)

      DATA W    
     */1.01228536290376E-01, 2.22381034453374E-01, 3.1370664587789 E-01,    
     * 3.6268378337836 E-01, 2.71524594117541E-02, 6.2253523938648 E-02,    
     * 9.5158511682493 E-02, 1.24628971255534E-01, 1.49595988816577E-01,    
     * 1.69156519395003E-01, 1.82603415044924E-01, 1.89450610455069E-01/    
   
      DATA X    
     */9.6028985649754 E-01, 7.9666647741363 E-01, 5.2553240991633 E-01,    
     * 1.83434642495650E-01, 9.8940093499165 E-01, 9.4457502307323 E-01,    
     * 8.6563120238783 E-01, 7.5540440835500 E-01, 6.1787624440264 E-01,    
     * 4.5801677765723 E-01, 2.8160355077926 E-01, 9.5012509837637 E-02/    
   
      END
  

      FUNCTION AGAUSS(F,A,B,N,RELERR)
        IMPLICIT NONE
	
        DOUBLE PRECISION AGAUSS,A,B,RELERR,XR,SS,XM,DX,F,AAGAUS
        INTEGER N,I,J
        EXTERNAL F
        DOUBLE PRECISION W,X
          COMMON /GAUDAT/ W(12),X(12)
      IF (N.LE.0) THEN
        AGAUSS = AAGAUS(F,A,B,RELERR)
        RETURN
      END IF
      XR = (B-A)/2/N
      SS = 0
      DO 20, I = 1,N
        XM = A + (I-0.5D0)*(B-A)/N
        DO 10, J = 5,12
          DX = XR*X(J)
          SS = SS + W(J)*(F(XM+DX)+F(XM-DX))
10      CONTINUE
20    CONTINUE
      AGAUSS = XR*SS
      RETURN
      END
               

      FUNCTION BGAUSS(F,A,B,N,RELERR)
        IMPLICIT NONE
        
        DOUBLE PRECISION BGAUSS,A,B,RELERR,XR,SS,XM,DX,F,BBGAUS
        INTEGER N,I,J
        EXTERNAL F
        DOUBLE PRECISION W,X
          COMMON /GAUDAT/ W(12),X(12)
      IF (N.LE.0) THEN
        BGAUSS = BBGAUS(F,A,B,RELERR)
        RETURN
      END IF
      XR = (B-A)/2/N
      SS = 0
      DO 20, I = 1,N
        XM = A + (I-0.5D0)*(B-A)/N
        DO 10, J = 5,12
          DX = XR*X(J)
          SS = SS + W(J)*(F(XM+DX)+F(XM-DX))
10      CONTINUE
20    CONTINUE
      BGAUSS = XR*SS
      RETURN
      END


      FUNCTION XGAUSS(F,A,B,N,RELERR)
        IMPLICIT NONE

        DOUBLE PRECISION XGAUSS,A,B,RELERR,XR,SS,XM,DX,F,XXGAUS
        INTEGER N,I,J
        EXTERNAL F
        DOUBLE PRECISION W,X
          COMMON /GAUDAT/ W(12),X(12)
      IF (N.LE.0) THEN
        XGAUSS = XXGAUS(F,A,B,RELERR)
        RETURN
      END IF
      XR = (B-A)/2/N
      SS = 0
      DO 20, I = 1,N
        XM = A + (I-0.5D0)*(B-A)/N
        DO 10, J = 5,12
          DX = XR*X(J)
          SS = SS + W(J)*(F(XM+DX)+F(XM-DX))
10      CONTINUE
20    CONTINUE
      XGAUSS = XR*SS
      RETURN
      END
               

      FUNCTION AAGAUS(F,A,B,EPS) 
        IMPLICIT NONE

        DOUBLE PRECISION AAGAUS,F,A,B,EPS,AA,BB,C1,C2,S8,U,S16
        DOUBLE PRECISION ERRLIM,ERRLM0,CONST
        INTEGER I
        LOGICAL FIRST,MFLAG,RFLAG   
        EXTERNAL F    
        DOUBLE PRECISION W,X
          COMMON /GAUDAT/ W(12),X(12) 
      AAGAUS=0.  
      IF(B.EQ.A) RETURN 
      CONST=0.005/(B-A)
      FIRST = .TRUE. 
      BB=A  
C   
C  COMPUTATIONAL LOOP.  
    1 AA=BB 
      BB=B  
    2    C1=0.5*(BB+AA) 
         C2=0.5*(BB-AA) 
         S8=0.  
         DO 3 I=1,4 
            U=C2*X(I)   
            S8=S8+W(I)*(F(C1+U)+F(C1-U))   
    3    CONTINUE   
         S8=C2*S8   
         S16=0. 
         DO 4 I=5,12    
            U=C2*X(I)   
            S16=S16+W(I)*(F(C1+U)+F(C1-U))  
    4    CONTINUE   
         S16=C2*S16
         ERRLIM = EPS*ABS(S16)
         IF (FIRST) THEN
C          ! SAVE PERMITTED ERROR BASED ON FIRST ESTIMATE.
           ERRLM0 = ERRLIM
           FIRST = .FALSE.
         ELSE
           ERRLIM = MAX(ERRLIM,ERRLM0*(BB-AA)/(B-A))
         END IF
         IF  (ABS(S16-S8) .LE. ERRLIM) GOTO 5   
         BB=C1  
         IF( 1.+ABS(CONST*C2) .NE. 1. ) GO TO 2 
      AAGAUS=0.  
      WRITE(*,6)  
      RETURN    
    5 AAGAUS=AAGAUS+S16   
      IF(BB.NE.B) GO TO 1
      RETURN    
C   
    6 FORMAT( 4X, 'FUNCTION AAGAUS ... TOO HIGH ACCURACY REQUIRED')  
      END   
                                            

      FUNCTION BBGAUS(F,A,B,EPS) 
        IMPLICIT NONE

        DOUBLE PRECISION BBGAUS,F,A,B,EPS,AA,BB,C1,C2,S8,U,S16
        DOUBLE PRECISION ERRLIM,ERRLM0,CONST
        LOGICAL FIRST,MFLAG,RFLAG   
        INTEGER I
        EXTERNAL F    
        DOUBLE PRECISION W,X
          COMMON /GAUDAT/ W(12),X(12) 
      BBGAUS=0.  
      IF(B.EQ.A) RETURN 
      CONST=0.005/(B-A)
      FIRST = .TRUE. 
      BB=A  
C   
C  COMPUTATIONAL LOOP.  
    1 AA=BB 
      BB=B  
    2    C1=0.5*(BB+AA) 
         C2=0.5*(BB-AA) 
         S8=0.  
         DO 3 I=1,4 
            U=C2*X(I)   
            S8=S8+W(I)*(F(C1+U)+F(C1-U))       
    3    CONTINUE   
         S8=C2*S8   
         S16=0. 
         DO 4 I=5,12    
            U=C2*X(I)   
            S16=S16+W(I)*(F(C1+U)+F(C1-U))  
    4    CONTINUE   
         S16=C2*S16
         ERRLIM = EPS*ABS(S16)
         IF (FIRST) THEN
C          ! SAVE PERMITTED ERROR BASED ON FIRST ESTIMATE.
           ERRLM0 = ERRLIM
           FIRST = .FALSE.
         ELSE
           ERRLIM = MAX(ERRLIM,ERRLM0*(BB-AA)/(B-A))
         END IF
         IF  (ABS(S16-S8) .LE. ERRLIM) GOTO 5   
         BB=C1  
         IF( 1.+ABS(CONST*C2) .NE. 1. ) GO TO 2 
      BBGAUS=0.  
      WRITE(*,6)  
      RETURN    
    5 BBGAUS=BBGAUS+S16   
      IF(BB.NE.B) GO TO 1   
      RETURN    
C   
    6 FORMAT( 4X, 'FUNCTION BBGAUS ... TOO HIGH ACCURACY REQUIRED')  
      END   

                
      FUNCTION XXGAUS(F,A,B,EPS) 
         IMPLICIT NONE

        DOUBLE PRECISION XXGAUS,F,A,B,EPS,AA,BB,C1,C2,S8,U,S16
        DOUBLE PRECISION ERRLIM,ERRLM0,CONST
        LOGICAL FIRST,MFLAG,RFLAG   
        INTEGER I
        EXTERNAL F    
        DOUBLE PRECISION W,X
          COMMON /GAUDAT/ W(12),X(12) 
      XXGAUS=0.  
      IF(B.EQ.A) RETURN 
      CONST=0.005/(B-A)
      FIRST = .TRUE. 
      BB=A  
C   
C  COMPUTATIONAL LOOP.  
    1 AA=BB 
      BB=B  
    2    C1=0.5*(BB+AA) 
         C2=0.5*(BB-AA) 
         S8=0.  
         DO 3 I=1,4 
            U=C2*X(I)   
            S8=S8+W(I)*(F(C1+U)+F(C1-U))   
    3    CONTINUE   
         S8=C2*S8   
         S16=0. 
         DO 4 I=5,12    
            U=C2*X(I)   
            S16=S16+W(I)*(F(C1+U)+F(C1-U))  
    4    CONTINUE   
         S16=C2*S16
         ERRLIM = EPS*ABS(S16)
         IF (FIRST) THEN
C          ! SAVE PERMITTED ERROR BASED ON FIRST ESTIMATE.
           ERRLM0 = ERRLIM
           FIRST = .FALSE.
         ELSE
           ERRLIM = MAX(ERRLIM,ERRLM0*(BB-AA)/(B-A))
         END IF
         IF  (ABS(S16-S8) .LE. ERRLIM) GOTO 5   
         BB=C1  
         IF( 1.+ABS(CONST*C2) .NE. 1. ) GO TO 2 
      XXGAUS=0.  
      WRITE(*,6)  
      RETURN    
    5 XXGAUS=XXGAUS+S16   
      IF(BB.NE.B) GO TO 1   
      RETURN    
C   
    6 FORMAT( 4X, 'FUNCTION XXGAUS ... TOO HIGH ACCURACY REQUIRED')  
      END   


C ***********************************************************
C *  ZGAUSS IS LIKE AGAUGG BUT TAKES A LAST ARGUMENT, A FUNCTION, WHICH
C *  IT PASSES ALONG TO THE INTEGRAND.


      FUNCTION ZGAUSS(F,A,B,N,RELERR,FARG)
        IMPLICIT NONE

        DOUBLE PRECISION ZGAUSS,A,B,RELERR,XR,SS,XM,DX,F,ZZGAUS,FARG
        INTEGER N,I,J
        EXTERNAL F,FARG
        DOUBLE PRECISION W,X
          COMMON /GAUDAT/ W(12),X(12)
      IF (N.LE.0) THEN
        ZGAUSS = ZZGAUS(F,A,B,RELERR,FARG)
        RETURN
      END IF
      XR = (B-A)/2/N
      SS = 0
      DO 20, I = 1,N
        XM = A + (I-0.5D0)*(B-A)/N
        DO 10, J = 5,12
          DX = XR*X(J)
          SS = SS + W(J)*(F(XM+DX,FARG)+F(XM-DX,FARG))
10      CONTINUE
20    CONTINUE
      ZGAUSS = XR*SS
      RETURN
      END
               

      FUNCTION ZZGAUS(F,A,B,EPS,FARG) 
        IMPLICIT NONE

        DOUBLE PRECISION ZZGAUS,F,A,B,EPS,AA,BB,C1,C2,S8,U,S16,FARG
        DOUBLE PRECISION ERRLIM,ERRLM0,CONST
        INTEGER I
        LOGICAL FIRST,MFLAG,RFLAG   
        EXTERNAL F,FARG
        DOUBLE PRECISION W,X
          COMMON /GAUDAT/ W(12),X(12)
      ZZGAUS=0.  
      IF(B.EQ.A) RETURN 
      CONST=0.005/(B-A)
      FIRST = .TRUE. 
      BB=A  
C   
C  COMPUTATIONAL LOOP.  
    1 AA=BB 
      BB=B  
    2    C1=0.5*(BB+AA) 
         C2=0.5*(BB-AA) 
         S8=0.  
         DO 3 I=1,4 
            U=C2*X(I)   
            S8=S8+W(I)*(F(C1+U,FARG)+F(C1-U,FARG))   
    3    CONTINUE   
         S8=C2*S8   
         S16=0. 
         DO 4 I=5,12    
            U=C2*X(I)   
            S16=S16+W(I)*(F(C1+U,FARG)+F(C1-U,FARG))  
    4    CONTINUE   
         S16=C2*S16
         ERRLIM = EPS*ABS(S16)
         IF (FIRST) THEN
C          ! SAVE PERMITTED ERROR BASED ON FIRST ESTIMATE.
           ERRLM0 = ERRLIM
           FIRST = .FALSE.
         ELSE
           ERRLIM = MAX(ERRLIM,ERRLM0*(BB-AA)/(B-A))
         END IF
         IF  (ABS(S16-S8) .LE. ERRLIM) GOTO 5   
         BB=C1  
         IF( 1.+ABS(CONST*C2) .NE. 1. ) GO TO 2 
      ZZGAUS=0.  
      WRITE(*,6)  
      RETURN    
    5 ZZGAUS=ZZGAUS+S16   
      IF(BB.NE.B) GO TO 1
      RETURN    
C   
    6 FORMAT( 4X, 'FUNCTION ZZGAUS ... TOO HIGH ACCURACY REQUIRED')  
      END   
                                            

      FUNCTION GAUSS8(F,A,B,N,RELERR)
        IMPLICIT NONE

CAB-BUG!        DOUBLE PRECISION AGAUSS,A,B,RELERR,XR,SS,XM,DX,F,AAGAUS

        DOUBLE PRECISION GAUSS8,A,B,RELERR,XR,SS,XM,DX,F,AAGAUS
        INTEGER N,I,J
        EXTERNAL F
        DOUBLE PRECISION W,X
          COMMON /GAUDAT/ W(12),X(12)
      IF (N.LE.0) PAUSE 'GAUSS8: N <=0.'
      XR = (B-A)/2/N
      SS = 0
      DO 20, I = 1,N
        XM = A + (I-0.5D0)*(B-A)/N
        DO 10, J = 1,4
          DX = XR*X(J)
          SS = SS + W(J)*(F(XM+DX)+F(XM-DX))
10      CONTINUE
20    CONTINUE
      GAUSS8 = XR*SS
      RETURN
      END

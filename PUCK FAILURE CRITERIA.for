C     Written by Karan Kodagali @ University of South Carolina at Cololumbia     
C     Based largely on the theories of Alfred Puck
C     Please send any bugs or errors to kodagali@email.sc.edu
      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATEV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
C
      INCLUDE 'ABA_PARAM.INC'
C
      INTEGER
     1 I,J
C
      CHARACTER*80 CMNAME
      CHARACTER*80 CPNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATEV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)
      DOUBLE PRECISION    
     1 E1,E2,E3,G12,G13,G23,V12,V13,V23,V21,V31,V32,S(NTENS),
     2 XT,XC,YT,YC,VF12,EF1,
     3 C(NTENS,NTENS),S21,FFT,FFC,MFT,MFC,NMP,NW1,G23O,
     4 MFLC,DFT,DFC,DMT,DMC,DMLC,STRANT(NTENS),CD(NTENS,NTENS),
     5 T21C,RVVA,ZERO,ONE,CF(NTENS,NTENS),MAXIM,E2O,G12O,
     6 CFULL(NTENS,NTENS),MFMV,DMG,COUNTER,XX,VAR,E1O,XTF
      INTEGER DEG,MAT
      PARAMETER (ZERO=0.D0, ONE=1.D0)
C
C     INITIALIZING VARIABLES---------------------------------------------
C
      E1 = PROPS(1)           !YOUNG'S MODULUS IN DIRECTION 1 (L)
      E2 = PROPS(2)           !YOUNG'S MODULUS IN DIRECTION 2 (T) 
      E3=E2               
      G12 = PROPS(3)          !SHEAR MODULUS IN 12 PLANE
      G13=G12                 !SHEAR MODULUS IN 13 PLANE
      V12=PROPS(4)            !POISSON RATIO IN 12
      V23=PROPS(5)            !POISSON RATIO IN 23
      V13=V12                 !POISSON RATIO IN 13 
      EF1=PROPS(6)            !MODULUS OF FIBER PARALLEL TO FIBER
	VF12 = PROPS(7)         !VOLUME FRACTION OF FIBER
      XT = PROPS(8)           !TENSILE STRENGTH PARALLEL TO FIBER
      XC = PROPS(9)           !COMPRESSIVE STRENGTH PARALLEL TO FIBER
      YT = PROPS(10)          !TENSILE STRENGTH PERPENDICULAR TO FIBER
      YC = PROPS(11)          !COMPRESSIVE STRENGTH PERPENDICULAR TO FIBER
      S21 = PROPS(12)         !IN PLANE SHEAR STRENGTH
      MAT = PROPS(13)         !MATERIAL TYPE FOR INCLINATION PARAMETERS
      DEG = PROPS(14)         !EWM/CSE
      G23 = E2/2/(1.+V23)     !SHEAR MODULUS IN 23 PLANE
      XTF=XT*EF1/E1           !EFFECTIVE TENSILE STRENGTH OF FIBER
C     Damage variables from previous time step
      DFT=STATEV(6)       
      DFC=STATEV(7)
      DMT=STATEV(8)
      DMC=STATEV(9)
C     Saving original stiffness values before degradation
      E2O=E2
      G12O=G12
      G23O=G23
      E1O=E1
C
      COUNTER=ZERO
      DMG=ONE
      NMP=1
      NW1=1
      V21=(E2/E1)*V12
      V31=(E3/E1)*V13
      V32=(E3/E2)*V23
C
C     STRAIN---------------------------------------------
C     
      DO I = 1, NTENS
         STRANT(I) = STRAN(I) + DSTRAN(I)
      END DO
C
C     DEGRADATION DUE TO PREVIOUS DAMAGE---------------------------------------------
C
      IF(DEG.EQ.1) THEN
          E1=(1-DFC)*(1-DFT)*E1O
          G12=(1-DFT)*(1-(DMC*0.77))*(1-(DMT*0.77))*G12O
          E2=(1-DFT)*(1-DMT)*E2O
          G13=G12
          G23=(1-DFT)*(1-(DMC*0.77))*(1-(DMT*0.77))*G23O
      END IF
C
C     CONSTITUTIVE RESPONSE AND STRESS---------------------------------------------
C
      CALL CONSTITUTIVE(CF,E1,E2,E3,G12,G13,G23,V12,V13,V23,V21,
     1 V31,V32,NTENS,NDI,NSHR,DFT,DFC,DMT,DMC,DEG)
      DO K1=1,NTENS
          S(K1)=0.0D0
          DO K2=1,NTENS
              S(K1)=S(K1)+CF(K2,K1)*STRANT(K2)
          ENDDO
      ENDDO
      STRESS=S
      DDSDDE=CF   !Updating jacobian and stress is carried out with degradation at the previous load step to ensure easier convergence
C                 !Thus a small time step is required for accurate results.
C     DAMAGE EVALUATION AND DEGRADATION---------------------------------------------
C    
      DO WHILE ( DMG.EQ.ONE )  
          CALL CONSTITUTIVE(CF,E1,E2,E3,G12,G13,G23,V12,V13,V23,V21,
     1    V31,V32,NTENS,NDI,NSHR,DFT,DFC,DMT,DMC,DEG)
          DO K1=1,NTENS
              S(K1)=0.0D0
              DO K2=1,NTENS
                  S(K1)=S(K1)+CF(K2,K1)*STRANT(K2)
              ENDDO
          ENDDO
          CALL THETAFP(S,S21,XTF,XC,YT,YC,THETA,NTENS,NMP,MAXIM,NDI,
     1     NSHR,MAT)
          CALL CFAILURE(S,V12,VF12,EF1,S21,XTF,XC,YT,YC,NDI,NSHR,
     1     MFT,FFT,FFC,MFC,DMG,NTENS,THETA,NMP,NW1,MAXIM,SIGN,MAT)
          IF(DFT.GT.0.9) FFT=1
          IF(DFC.GT.0.9) FFC=1
          IF(DMT.EQ.0.99) MFT=1
          IF(DMC.EQ.0.99) MFC=1          
          IF(FFT.LE.ONE.AND.FFC.LE.ONE.AND.MFT.LE.ONE.AND.MFC.LE.
     1    ONE) THEN 
             DMG=ZERO
          END IF
C     Degradation due to damage in step for CSE/EWM
          IF(DEG.EQ.1) THEN
              IF(FFT.GT.ONE) THEN
                  DFT=0.99
              END IF
              IF(FFC.GT.ONE) THEN
                  DFC=0.99
              END IF
              IF(MFT.GT.ONE) THEN
                  DMT=DMT+ONE/100
                  IF(DMT.GE.ONE) THEN
                        DMT=0.99
                  END IF
              END IF
              IF(MFC.GT.ONE) THEN
                  DMC=DMC+ONE/100
                  IF(DMC.GE.ONE) THEN
                        DMC=0.99
                  END IF
              END IF
              E1=(1-DFC)*(1-DFT)*E1O
              G12=(1-DFT)*(1-(DMC*0.7))*(1-(DMT*0.7))*G12O
              E2=(1-DFT)*(1-DMT)*E2O
              G13=G12
              G23=(1-DFT)*(1-(DMC*0.7))*(1-(DMT*0.7))*G23O
          ELSE IF(DEG.EQ.2) THEN
              IF(FFT.GT.ONE) THEN
                  DFT=0.99
              END IF
                IF(FFC.GT.ONE) THEN
                  DFC=0.99
              END IF
              IF(MFT.GT.ONE.AND.DMT.EQ.ZERO) THEN
                  DMT=0.99
              END IF
              IF(MFC.GT.ONE.AND.DMC.EQ.ZERO) THEN
                  DMC=0.99     
              END IF
          END IF
      END DO
C
C     SAVE STATE VARIABLES---------------------------------------------	  
C
      STATEV(1) = FFT
      STATEV(2) = FFC
      STATEV(3) = MFT
      STATEV(4) = MFC
      STATEV(5) = THETA
      STATEV(6) = DFT
      STATEV(7) = DFC
      STATEV(8) = DMT
      STATEV(9) = DMC
      STATEV(10) = MAX(DFT,DFC,DMT,DMC)
      STATEV(11) = SIGN
      STATEV(12) = E1
      STATEV(13) = FFT
      STATEV(14) = FFC
      STATEV(15) = MFT
      STATEV(16) = MFC
      RETURN
      END
C******************************************************************************
C CALCULATE THE CONSTITUTIVE RESPONSE******************************************
C******************************************************************************
      SUBROUTINE CONSTITUTIVE(CF,E1,E2,E3,G12,G13,G23,V12,V13,V23,V21,
     1 V31,V32,NTENS,NDI,NSHR,DFT,DFC,DMT,DMC,DEG)
      INCLUDE 'ABA_PARAM.INC'
      DOUBLE PRECISION
     1 E1,E2,G12,G23,V12,V13,V23,DMG,G13,V21,V31,V32,E3,DF,DM,DMT,DMC,
     2 CF(NTENS,NTENS),S(NTENS),STRANT(6),ATEMP,DELTA,SMT,SMC,DFT,DFC
      INTEGER NDI,NTENS,DEG,NSHR
      PARAMETER (ZERO=0.D0, ONE=1.D0)
      DO K1=1,NTENS
          DO K2=1,NTENS
              CF(K1,K2)=0.D0
          ENDDO
      ENDDO
      IF(NDI.EQ.3) THEN
          IF(DEG.EQ.1) THEN
              DELTA=1/(1-V12*V21-V23*V32-V13*V31-2*V21*V32*V13)
              CF(1,1) = E1*(1-V23*V32)*DELTA
	        CF(1,2) = E2*(V12+V32*V13)*DELTA
	        CF(1,3) = E1*(V31+V21*V32)*DELTA
	        CF(2,1) = CF(1,2)   
	        CF(2,2) = E2*(1-V13*V31)*DELTA     
	        CF(2,3) = E2*(V32+V12*V31)*DELTA      
              CF(3,1) = CF(1,3)    
	        CF(3,2) = CF(2,3)    
	        CF(3,3) = E3*(1-V12*V21)*DELTA      
	        CF(4,4)	= G12   
	        CF(5,5)	= G13   
	        CF(6,6)	= G23  
          ELSE
              SMT=0.9
              SMC=0.5
              DF=1-(1-DFT)*(1-DFC)
              DM=1-(1-DMT)*(1-DMMC)
              DELTA=1/(1-V12*V21-V23*V32-V13*V31-2*V21*V32*V13)
              CF(1,1) = (1-DF)*E1*(1-V23*V32)*DELTA
	        CF(1,2) = (1-DF)*(1-DM)*E1*(V21+V31*V23)*DELTA
	        CF(1,3) = (1-DF)*(1-DM)*E1*(V31+V21*V32)*DELTA
	        CF(2,1) = CF(1,2)   
	        CF(2,2) = (1-DF)*(1-DM)*E2*(1-V13*V31)*DELTA     
	        CF(2,3) = (1-DF)*(1-DM)*E2*(V32+V12*V31)*DELTA      
              CF(3,1) = CF(1,3)    
	        CF(3,2) = CF(2,3)    
	        CF(3,3) = (1-DF)*(1-DM)*E3*(1-V12*V21)*DELTA      
	        CF(4,4)	= (1-DF)*(1-SMT*DMT)*(1-SMC*DMC)*G12   
	        CF(5,5)	= (1-DF)*(1-SMT*DMT)*(1-SMC*DMC)*G13   
	        CF(6,6)	= (1-DF)*(1-SMT*DMT)*(1-SMC*DMC)*G23 
          END IF
      ELSE IF(NDI.EQ.2) THEN
          IF(DEG.EQ.1) THEN
              DELTA = 1-V12*V21
              CF(1,1) = E1/DELTA
              CF(2,2) = E2/DELTA
              CF(1,2) = V12*E2/DELTA
              CF(2,1) = CF(1,2)
              CF(3,3) = G12
              IF(NSHR.GT.1) THEN
                  CF(4,4) = G13
                  CF(5,5) = G23
              END IF
          ELSE
              DELTA = 1-V12*V21
              CF(1,1) = (1-DF)*E1/DELTA
              CF(2,2) = (1-DF)*(1-DM)*E2/DELTA
              CF(1,2) = (1-DF)*(1-DM)*V12*E2/DELTA
              CF(2,1) = (1-DF)*(1-DM)*CF(1,2)
              CF(3,3) = (1-DF)*(1-SMT*DMT)*(1-SMC*DMC)*G12
              IF(NSHR.GT.1) THEN
                  CF(4,4) = (1-DF)*(1-SMT*DMT)*(1-SMC*DMC)*G13
                  CF(5,5) = (1-DF)*(1-SMT*DMT)*(1-SMC*DMC)*G23
              END IF
          END IF
      END IF
      RETURN
      END
C******************************************************************************
C ANGLE OF FRACTURE PLANE******************************************************
C******************************************************************************
      SUBROUTINE THETAFP(S,S21,XT,XC,YT,YC,THETA,NTENS,NMP,MAXIM,NDI,
     1 NSHR,MAT)
      INCLUDE 'ABA_PARAM.INC'
      INTEGER NTENS,NDI,NSHR,MAT
      DOUBLE PRECISION
     1 P21T,P21C,P22C,S21,XT,XC,YT,YC,SIG13,SIG23,
     2 S(NTENS),FE,MAXIM,MAXT,SIG11,SIG22,SIG33,SIG12,
     3 RVVA,THETA,P22T,TAUNT,SIGN,TAUNL,PTR,PCR,COS2PSI,
     4 SIN2PSI,SMAX,SREF,SI,NMP
      PARAMETER (ZERO=0.D0, ONE=1.D0)
      IF(MAT.EQ.1) THEN
          P21T = 0.3
	    P21C = 0.25
      	P22C = 0.2
          P22T=P22C
      ELSE
          P21T = 0.35
	    P21C = 0.3
      	P22C = 0.25
          P22T=P22C
      END IF
      IF(NDI.EQ.3) THEN
          SIG11=S(1)
          SIG22=S(2)
          SIG33=S(3)
          SIG12=S(4)
          SIG13=S(5)
          SIG23=S(6)
      ELSE
          SIG11=S(1)
          SIG22=S(2)
          SIG12=S(3)
          IF(NSHR.GT.1) THEN
              SIG33=0
              SIG13=S(4)
              SIG23=S(5)
          ELSE
              SIG33=0
              SIG13=0
              SIG23=0
          END IF
      END IF
      SMAX=90
      SREF=30
      SI=0
      RVVA = (S21/(2*P21C))*(sqrt((1+2*P21C*YC/S21))-1)
      DO I = -90,90
          THETA = I
          SIGN=SIG22*(COS(THETA))**2+SIG33*(SIN(THETA))**2+2*SIG23*
     1    SIN(THETA)*COS(THETA)
          TAUNT=-SIG22*SIN(THETA)*COS(THETA)+SIG33*SIN(THETA)*COS(THETA)
     1    +SIG23*((COS(THETA))**2-(SIN(THETA))**2)
          TAUNL=SIG13*SIN(THETA)+SIG12*COS(THETA)
          COS2PSI=TAUNT**2/(TAUNT**2+TAUNL**2)
          SIN2PSI=TAUNL**2/(TAUNT**2+TAUNL**2)
          PTR=(P22T/RVVA)*COS2PSI+(P21T/S21)*SIN2PSI
          PCR=(P22C/RVVA)*COS2PSI+(P21C/S21)*SIN2PSI
          IF(SIGN.GE.ZERO) THEN
              FE =SQRT((((1/YT)-PTR)*SIGN)**2+(TAUNT/RVVA)**2+
     1        (TAUNL/S21)**2)+PTR*SIGN
          ELSE
              FE = SQRT((TAUNT/RVVA)**2+(TAUNL/S21)**2+(PCR*SIGN)**2)+
     1        PCR*SIGN
          END IF
          IF(FE.GT.MAXIM) THEN
              MAXIM=FE
              MAXT=THETA
          END IF
          IF(FE.GT.0.5) THEN
              SI=SI+(FE-0.5)
          END IF
      END DO
      THETA=MAXT
      IF(SI.LT.SREF) THEN
          NMP=1
      ELSE
          NMP=1-0.2*(SI-SREF)/(SMAX-SREF)
      END IF
      IF(NMP.LT.0.5) THEN
          NMP=0.5
      END IF
      RETURN
      END
C******************************************************************************
C PUCK FAILURE CRITERIA********************************************************
C******************************************************************************
      SUBROUTINE CFAILURE(S,V12,VF12,EF1,S211,XT,XC,YT1,YC1,NDI,NSHR,
     1 MFT,FFT,FFC,MFC,DMG,NTENS,THETA,NMP,NW1,MAXIM,SIGN,MAT)
      INCLUDE 'ABA_PARAM.INC'
      INTEGER NTENS
      DOUBLE PRECISION
     1 P21T,P21C,P22C,MSIG,VF12,EF1,S21,XT,XC,YT,YC,MAXIM,
     2 S(NTENS),STRANT(NTENS),MFT,FFT,FFC,MFC,DMG,V12,SIG23,
     3 T21C,RVVA,THETA,P22T,TAUNT,SIGN,TAUNL,PTR,PCR,COS2PSI,SIG13,
     4 SIN2PSI,NMP,M,SC,A,C,NW1,YT1,YC1,S211,SIG11,SIG22,SIG33,SIG12
      PARAMETER (ZERO=0.D0, ONE=1.D0)
      IF(NDI.EQ.3) THEN
          SIG11=S(1)
          SIG22=S(2)
          SIG33=S(3)
          SIG12=S(4)
          SIG13=S(5)
          SIG23=S(6)
      ELSE
          SIG11=S(1)
          SIG22=S(2)
          SIG12=S(3)
          IF(NSHR.GT.1) THEN
              SIG33=0
              SIG13=S(4)
              SIG23=S(5)
          ELSE
              SIG33=0
              SIG13=0
              SIG23=0
          END IF
      END IF
	IF(MAT.EQ.1) THEN
          P21T = 0.3
	    P21C = 0.25
      	P22C = 0.2
          P22T = P22C
          MSIG = 1.3
      ELSE
          P21T = 0.35
	    P21C = 0.3
      	P22C = 0.25
          P22T = P22C
          MSIG = 1.1
      END IF
c	RVVA = YC/(2*(1+P22C))    
C
C     FAILURE CRITERIA	  
C
!     FIBER TENSILE
      IF((SIG11-(V12-VF12*MSIG*E1/EF1)*(SIG22+SIG33)).GE.ZERO) THEN
          FFT = (SIG11-(V12-VF12*MSIG*E1/EF1)*(SIG22+SIG33))/XT
          FFC = ZERO
          IF(FFT.GE.ONE) THEN
		    DMG=ONE
          END IF
      ELSE
!     FIBER COMPRESSIVE
          FFC = ABS((SIG11-(V12-VF12*MSIG*E1/EF1)*(SIG22+SIG33)))/XC
          FFT = ZERO
          IF(FFC.GE.ONE) THEN
	 	    DMG=ONE
          END IF
      END IF
      M=0.5
      SC=0.5
      A=(1-SC)/SQRT(1-M**2)
      C=MAX(FFT,FFC)/MAXIM
      IF(C.GE.0.5.AND.C.LE.2) THEN
          NW1=C*(A*SQRT(C**2*(A**2-SC**2)+1)+SC)/((C*A)**2+1)
      ELSE
          NW1=1
      END IF
      YT=YT1/(NW1*NMP)
      YC=YC1/(NW1*NMP)
      S21=S211/(NW1*NMP)
      RVVA = (S21/(2*P21C))*(sqrt((1+2*P21C*YC/S21))-1)
      SIGN=SIG22*(COS(THETA))**2+SIG33*(SIN(THETA))**2+2*SIG23*SIN(THETA
     1 )*COS(THETA)
      TAUNT=-SIG22*SIN(THETA)*COS(THETA)+SIG33*SIN(THETA)*COS(THETA)+
     1 SIG23*((COS(THETA))**2-(SIN(THETA))**2)
      TAUNL=SIG13*SIN(THETA)+SIG12*COS(THETA)
      COS2PSI=TAUNT**2/(TAUNT**2+TAUNL**2)
      SIN2PSI=TAUNL**2/(TAUNT**2+TAUNL**2)
      PTR=(P22T/RVVA)*COS2PSI+(P21T/S21)*SIN2PSI
      PCR=(P22C/RVVA)*COS2PSI+(P21C/S21)*SIN2PSI
      SIGN=SIGN*NMP
      TAUNT=TAUNT*NMP
      TAUNL=TAUNL*NMP
!     MATRIX TENSILE
	IF(SIGN.GE.ZERO) THEN
          MFT =SQRT((((1/YT)-PTR)*SIGN)**2+(TAUNT/RVVA)**2+(TAUNL/S21)
     1    **2)+PTR*SIGN
          MFC=ZERO
      END IF
!     MATRIX  COMPRESSION
	IF(SIGN.LT.ZERO) THEN    
	    MFC = SQRT((TAUNT/RVVA)**2+(TAUNL/S21)**2+(PCR*SIGN)**2)+PCR*SIGN
          MFT=ZERO
      END IF
      MFT=MFT/(NW1*NMP)
      MFC=MFC/(NW1*NMP)
      IF(MFC.GE.ONE.OR.MFT.GE.ONE.OR.FFT.GE.ONE.OR.FFC.GE.ONE) THEN
          DMG=ONE
      ELSE
          DMG=ZERO
      END IF
      RETURN
      END
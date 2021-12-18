! Roger Bellido Peralta
! Modelització Molecular
! Pràctica final
! 10/12/21

	PROGRAM PRF
	 IMPLICIT NONE
	 DOUBLE PRECISION NU,SIGMA_LJ,EPS,DT,SIGMA_T,RHO,CUTOFF,L,V_DT,M
	 DOUBLE PRECISION AV,RHO_2
	 INTEGER N,DM,TS,I
	 COMMON/LJ/SIGMA_LJ,EPS,CUTOFF
	 COMMON/CONSTANTS/DT,RHO,L,RHO_2
	 COMMON/THERMO/NU,SIGMA_T
       REAL*4 TIME1,TIME2
       CHARACTER*30 DATE1,DATE2
	 PARAMETER(N=125,DM=3,TS=500000)
	 DIMENSION V_DT(4),RHO_2(4)
	
       CALL FDATE(DATE1)
       CALL CPU_TIME(TIME1)
       WRITE(*,*)"BEGINNING:",DATE1

	 OPEN(10,FILE='PRf_en_eu_10-3.dat')
	 OPEN(11,FILE='PRf_vel_eu_10-3.dat')

! First part
	 RHO=0.7D0
	 SIGMA_LJ=1D0
	 EPS=1D0
	 V_DT=(/1D-1,1D-2,1D-3,1D-4/)
	 NU=1D0
	 L=(DBLE(N)/RHO)**(1D0/DBLE(DM))
	 CUTOFF=L/3D0

	 DT=V_DT(3)
	 CALL APARTAT_1(N,DM,TS)
	



	 CLOSE(10)
	 CLOSE(11)

	 
       CALL CPU_TIME(TIME2)
       CALL FDATE(DATE2)
       WRITE(*,*)"ENDING...:",DATE2
       WRITE(*,*)"CPUTIME=",TIME2-TIME1,'s ', (TIME2-TIME1)/60d0,'m'

	
	END PROGRAM PRF


!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-------------------ZONA DE SUBRUTINES----------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

! Function to apply the Periodic Boundary Conditions. In this case, we
! consider a cage initialised between 0 and L. 
	FUNCTION PBC(D,L)
	 IMPLICIT NONE
	 DOUBLE PRECISION PBC,D,L
	 
	 IF ((D.LE.(L/2D0)) .AND. (D.GE.(-L/2D0))) THEN
	  PBC=D
	 ELSE IF (D.GT.(L/2D0)) THEN
	  PBC=D-L
	 ELSE IF (D.LT.(-L/2D0)) THEN
	  PBC=D+L
	 ENDIF

	RETURN
	END

! Subroutine that calculates gaussean distributed random numbers 
! providing it of uniform distributed ones.	
      SUBROUTINE BOXMULLER(SIGMA,MU,X1,X2,XOUT1,XOUT2)
       IMPLICIT NONE
       double precision PI, sigma, MU, x1, x2, xout1, xout2
       PI = 4d0*datan(1d0)
       
       XOUT1=sigma*dsqrt(-2d0*(dlog(1d0-x1)))*dcos(2d0*PI*x2)+MU
       XOUT2=sigma*dsqrt(-2d0*(dlog(1d0-x1)))*dsin(2d0*PI*x2)+MU
       
      RETURN
      END

! Subroutine that initialises a simple cubic lattice.
	SUBROUTINE INIT_SCC(N,DM,R)
	 IMPLICIT NONE
	 DOUBLE PRECISION L,A,RHO,R,DT,RHO_2
	 INTEGER N,I,J,K,DM,M,NPART
	 COMMON/CONSTANTS/DT,RHO,L,RHO_2
	 DIMENSION R(N,DM),RHO_2(4)

	 M=NINT(DBLE(N)**(1D0/DBLE(DM)))
	 A=L/DBLE(M)
	 NPART=1
	 DO I=1,M
	  DO J=1,M
	   DO K=1,M
	    R(NPART,1)=DBLE(A*I)
	    R(NPART,2)=DBLE(A*J)
	    R(NPART,3)=DBLE(A*K)
	    NPART=NPART+1
	   ENDDO
	  ENDDO
	 ENDDO

	 RETURN
	 END

	SUBROUTINE INIT_BIMODAL(N,DM,T,V)
	 IMPLICIT NONE
	 DOUBLE PRECISION V,SIGMA_LJ,EPS,CUTOFF
	 DOUBLE PRECISION X,KIN,T_INST,T
	 INTEGER N,DM,I,J
	 COMMON/LJ/SIGMA_LJ,EPS,CUTOFF
	 DIMENSION V(N,DM)

	 DO I=1,N
	  DO J=1,DM
	   CALL RANDOM_NUMBER(X)
	   IF (X.GT.0.5D0) THEN
	    V(I,J)=DSQRT(T)
	   ELSE
	    V(I,J)=-DSQRT(T)
	   ENDIF
	  ENDDO
	 ENDDO
	  KIN=0D0
	  DO J=1,N
	   KIN=KIN+0.5D0*(V(J,1)**2D0+V(J,2)**2D0+V(J,3)**2D0)
	  ENDDO	

	  T_INST=2D0/DBLE(3*N-3)*KIN
	  WRITE(*,*)'T_INST= ', T_INST

	 
	 
	RETURN
	END

! Subroutine that calculates distances between the particles of our
! system and the forces according to a Lennard-Jones potential. 
	SUBROUTINE EN_LJ(R,N,DM,L,F,POT,D)
	 IMPLICIT NONE
	 DOUBLE PRECISION R,L,CUTOFF,POT,F,SIGMA_LJ,EPS
	 DOUBLE PRECISION DX,DY,DZ,D,PBC,PROD
	 INTEGER I,J,N,DM
	 COMMON/LJ/SIGMA_LJ,EPS,CUTOFF
	 DIMENSION R(N,DM),F(N,DM)
	 
	 POT=0D0
	 F=0D0
	 PROD=1D0
	 DO I=1,N
	  DO J=I+1,N
	   DX=R(I,1)-R(J,1)
	   DY=R(I,2)-R(J,2)
	   DZ=R(I,3)-R(J,3)

	   DX=PBC(DX,L)
	   DY=PBC(DY,L)
	   DZ=PBC(DZ,L)

	   D=(DX**2D0+DY**2D0+DZ**2D0)**0.5D0

	   IF (D.LT.CUTOFF) THEN
	    F(I,1)=F(I,1)+(48D0/D**14D0-24D0/D**8D0)*DX
	    F(I,2)=F(I,2)+(48D0/D**14D0-24D0/D**8D0)*DY
	    F(I,3)=F(I,3)+(48D0/D**14D0-24D0/D**8D0)*DZ
	    
	    F(J,1)=F(J,1)-(48D0/D**14D0-24D0/D**8D0)*DX
	    F(J,2)=F(J,2)-(48D0/D**14D0-24D0/D**8D0)*DY
	    F(J,3)=F(J,3)-(48D0/D**14D0-24D0/D**8D0)*DZ
	    
	    POT=POT+4D0*(1D0/D**12D0-1D0/D**6D0)
     &       -4D0*(1D0/CUTOFF**12D0-1D0/CUTOFF**6D0)
         ENDIF
	  ENDDO
	 ENDDO
	RETURN
	END

! Integrator that uses the velocity of the current timestep to compute
! the positions. In each step, the calculation of the next velocity
! is needed. 
	SUBROUTINE VERLET_V(R,V,N,DM,TS,V_X,V_Y,V_Z,U,KE,E_TOT,T_INST,MOM)
	 IMPLICIT NONE
	 DOUBLE PRECISION T,R,F,ROLD,ROLDAUX,V,SIGMA,EPS,L,CUTOFF,POT
	 DOUBLE PRECISION DT,KIN,E_TOT,U,D,NU,SIGMA_R
	 DOUBLE PRECISION RHO,V_OUT,PBC,T_INST,KE,MOM,RHO_2
	 DOUBLE PRECISION V_X,V_Y,V_Z,R_F,V_F
	 INTEGER I,J,K,TS,N,DM,T_AN,W
	 COMMON/CONSTANTS/DT,RHO,L,RHO_2
	 DIMENSION R(N,DM),V(N,DM),F(N,DM),ROLD(N,DM),ROLDAUX(N,DM)
	 DIMENSION V_X(N,TS),V_Y(N,TS),V_Z(N,TS),KIN(N,TS)
	 DIMENSION E_TOT(TS),U(TS),T_INST(TS),KE(TS),MOM(TS)
	 DIMENSION V_OUT(N,DM),RHO_2(4)
	 DIMENSION R_F(N,DM),V_F(N,DM)
	 	
	 T=0D0
	 ROLD=R
	 MOM=0D0
	 DO I=1,TS
	  IF (MOD(I,100).EQ.0) THEN
	   WRITE(*,*)'I= ',I,'   TS= ', TS
	  ENDIF
	  T=T+DT
	  KE(I)=0D0
	  CALL EN_LJ(R,N,DM,L,F,POT,D)

	  ROLDAUX=R
	  DO J=1,N
	   R(J,1)=PBC(R(J,1)+V(J,1)*DT+F(J,1)*DT*DT/2D0,L)
	   R(J,2)=PBC(R(J,2)+V(J,2)*DT+F(J,2)*DT*DT/2D0,L)
	   R(J,3)=PBC(R(J,3)+V(J,3)*DT+F(J,3)*DT*DT/2D0,L)

	   V(J,1)=V(J,1)+F(J,1)*DT/2D0
	   V(J,2)=V(J,2)+F(J,2)*DT/2D0
	   V(J,3)=V(J,3)+F(J,3)*DT/2D0
	  ENDDO

	  CALL EN_LJ(R,N,DM,L,F,POT,D)

	  U(I)=POT
	  DO J=1,N
	   V(J,1)=V(J,1)+F(J,1)*DT/2D0
	   V(J,2)=V(J,2)+F(J,2)*DT/2D0
	   V(J,3)=V(J,3)+F(J,3)*DT/2D0
	  ENDDO
	  	  
	  DO J=1,N
	   KIN(J,I)=0.5D0*(V(J,1)**2D0+V(J,2)**2D0+V(J,3)**2D0)
	   KE(I)=KE(I)+KIN(J,I)
	  ENDDO	

	  T_INST(I)=2D0/DBLE(3*N-3)*KE(I)
	  ROLD=ROLDAUX
	  
	  DO K=1,N
	   V_X(K,I)=V(K,1)
	   V_Y(K,I)=V(K,2)
	   V_Z(K,I)=V(K,3)
	   
	   MOM(I)=MOM(I)+DSQRT(V(K,1)**2D0+V(K,2)**2D0+V(K,3)**2D0)/N
	  ENDDO
	 E_TOT(I)=KE(I)+U(I)

	 IF ((I.EQ.TS).AND.(W.EQ.1)) THEN
	  DO K=1,N
	  R_F(K,1)=R(K,1)
	  R_F(K,2)=R(K,2)
	  R_F(K,3)=R(K,3)
	  
	  V_F(K,1)=V(K,1)
	  V_F(K,2)=V(K,2)
	  V_F(K,3)=V(K,3)
	
	  ENDDO
	 ENDIF

	 ENDDO

	RETURN
	END

! Subroutine that computes the trajectory of a particle based on the 
! current position and velocity. 
	SUBROUTINE EULER(R,V,N,DM,TS,V_X,V_Y,V_Z,U,KE,E_TOT,T_INST,MOM)
	 IMPLICIT NONE
	 DOUBLE PRECISION T,R,F,ROLD,ROLDAUX,V,SIGMA,EPS,L,CUTOFF,POT
	 DOUBLE PRECISION TR_X,TR_Y,TR_Z,DT,KIN,E_TOT,U,D,PBC,RHO,RHO_2
	 DOUBLE PRECISION T_INST,KE,MOM,V_X,V_Y,V_Z
	 INTEGER I,J,K,TS,N,DM
	 COMMON/CONSTANTS/DT,RHO,L,RHO_2
	 DIMENSION R(N,DM),V(N,DM),F(N,DM),ROLD(N,DM),ROLDAUX(N,DM),U(TS)
	 DIMENSION TR_X(N,TS),TR_Y(N,TS),TR_Z(N,TS),KIN(N,TS),E_TOT(TS)
	 DIMENSION V_X(N,TS),V_Y(N,TS),V_Z(N,TS)
	 DIMENSION T_INST(TS),KE(TS),MOM(TS)
	 DIMENSION RHO_2(4)
	 	
	 T=0D0
	 ROLD=R
	 MOM=0D0
	 DO I=1,TS
	  IF (MOD(I,100).EQ.0) THEN
	   WRITE(*,*)'I= ',I,'   TS= ', TS
	  ENDIF
	  T=T+DT
	  KE(I)=0D0
	  CALL EN_LJ(R,N,DM,L,F,POT,D)
	  U(I)=POT
	  ROLDAUX=R
	  DO J=1,N
	   R(J,1)=PBC(R(J,1)+V(J,1)*DT+F(J,1)*DT*DT/2D0,L)
	   R(J,2)=PBC(R(J,2)+V(J,2)*DT+F(J,2)*DT*DT/2D0,L)
	   R(J,3)=PBC(R(J,3)+V(J,3)*DT+F(J,3)*DT*DT/2D0,L)

	   V(J,1)=V(J,1)+F(J,1)*DT/2D0
	   V(J,2)=V(J,2)+F(J,2)*DT/2D0
	   V(J,3)=V(J,3)+F(J,3)*DT/2D0
	  ENDDO
	  
	  DO J=1,N
	   KIN(J,I)=0.5D0*(V(J,1)**2D0+V(J,2)**2D0+V(J,3)**2D0)
	   KE(I)=KE(I)+KIN(J,I)
	  ENDDO	
	  T_INST(I)=2D0/DBLE(3*N-3)*KE(I)
	  ROLD=ROLDAUX
	  DO K=1,N
	   V_X(K,I)=V(K,1)
	   V_Y(K,I)=V(K,2)
	   V_Z(K,I)=V(K,3)
	   
	   MOM(I)=MOM(I)+DSQRT(V(K,1)**2D0+V(K,2)**2D0+V(K,3)**2D0)/N
	  ENDDO
	 E_TOT(I)=KE(I)+U(I)
	 ENDDO
	RETURN
	END
	

! Subroutine that computes the first part of the final exercise.
	SUBROUTINE APARTAT_1(N,DM,TS)
	 IMPLICIT NONE
	 DOUBLE PRECISION SIGMA_LJ,EPS,CUTOFF,DT,RHO,L,NU,SIGMA_T
	 DOUBLE PRECISION R,V,T,V_X,V_Y,V_Z,U,KE,E_TOT,T_INST,TEMP,MOM
	 DOUBLE PRECISION RHO_2
	 INTEGER N,DM,TS,I,T_AN	
	 COMMON/LJ/SIGMA_LJ,EPS,CUTOFF
	 COMMON/CONSTANTS/DT,RHO,L,RHO_2
	 COMMON/THERMO/NU,SIGMA_T
	 DIMENSION R(N,DM),V(N,DM)
	 DIMENSION V_X(N,TS),V_Y(N,TS),V_Z(N,TS)
	 DIMENSION E_TOT(TS),U(TS),T_INST(TS),KE(TS),MOM(TS)
	 DIMENSION RHO_2(4)

	 CALL INIT_SCC(N,DM,R)
	 TEMP=100D0
	 CALL INIT_BIMODAL(N,DM,TEMP,V)
	 WRITE(11,*)'# VELOCITATS INICIALS'
	  WRITE(11,101)V(:,1)
	  WRITE(11,101)V(:,2)
	  WRITE(11,101)V(:,3)
101	  FORMAT(1(F20.12))
	 WRITE(11,*)' '
	 WRITE(11,*)' '
	 T_AN=0
!	 CALL VERLET_V(R,V,N,DM,TS,V_X,V_Y,V_Z,U,KE,E_TOT,T_INST,MOM)
	    CALL EULER(R,V,N,DM,TS,V_X,V_Y,V_Z,U,KE,E_TOT,T_INST,MOM)
	 WRITE(11,*)'# VELOCITATS FINALS'
	  WRITE(11,101)V_X(:,TS)
	  WRITE(11,101)V_Y(:,TS)
	  WRITE(11,101)V_Z(:,TS)
	 T=0D0
	 WRITE(10,*)'T,EN_KIN,EN_POT,EN_TOT,T_INST,MOM_LIN'
	 DO I=1,TS
	  T=T+DT
	  WRITE(10,100)T,KE(I),U(I),E_TOT(I),T_INST(I),MOM(I)
100	  FORMAT(6(F20.12))	   
	 ENDDO
	 WRITE(*,*)'T_FINAL= ', T_INST(TS)

	RETURN
	END


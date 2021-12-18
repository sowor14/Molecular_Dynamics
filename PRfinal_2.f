! Roger Bellido Peralta
! Modelització Molecular
! Pràctica final
! 10/12/21

	PROGRAM PRF
	 IMPLICIT NONE
	 DOUBLE PRECISION NU,SIGMA_LJ,EPS,DT,SIGMA_T,RHO,CUTOFF,L,V_DT,M
	 DOUBLE PRECISION AV,RHO_2,KB
	 INTEGER N,DM,TS,I
	 COMMON/LJ/SIGMA_LJ,EPS,CUTOFF
	 COMMON/CONSTANTS/DT,RHO,L,RHO_2,M,AV,KB
	 COMMON/THERMO/NU,SIGMA_T
       REAL*4 TIME1,TIME2
       CHARACTER*30 DATE1,DATE2
	 PARAMETER(N=125,DM=3,TS=500000)
	 DIMENSION V_DT(4),RHO_2(4)
	
       CALL FDATE(DATE1)
       CALL CPU_TIME(TIME1)
       WRITE(*,*)"BEGINNING:",DATE1

	 OPEN(12,FILE='PRf_plots2.dat')
	 OPEN(13,FILE='PRf_Ar2.dat')

	 V_DT=(/1D-1,1D-2,1D-3,1D-4/)
	 DT=V_DT(4)
	
! Second part
	 EPS=0.998D0
	 SIGMA_LJ=3.4D0
	 M=40D0
	 AV=6.02214D23
	 KB=1.380649D-23
	 RHO_2=(/0.2D0,0.4D0,0.6D0,0.8D0/)
	 NU=1D0
	 
	 WRITE(13,*)'DENSITAT,<KE>,<U>,<E_TOT>,<T_INST>,<PRESS>'
	 DO I=3,4

	  RHO=RHO_2(I)
	  L=(DBLE(N)/RHO)**(1D0/DBLE(DM))
	  CUTOFF=L/3D0
	  CALL APARTAT_2(N,DM,TS)
	  WRITE(12,*)' '
	  WRITE(12,*)' '
	 ENDDO


	 CLOSE(12)
	 CLOSE(13)
	 
	 CALL STATS(TS)
	 

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
	 DOUBLE PRECISION L,A,RHO,R,DT,RHO_2,MASS,AV,KB
	 INTEGER N,I,J,K,DM,M,NPART
	 COMMON/CONSTANTS/DT,RHO,L,RHO_2,MASS,AV,KB
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

! Subroutine that calculates distances between the particles of our
! system and the forces according to a Lennard-Jones potential. 
	SUBROUTINE EN_LJ(R,N,DM,L,F,POT)
	 IMPLICIT NONE
	 DOUBLE PRECISION R,L,CUTOFF,POT,F,SIGMA_LJ,EPS
	 DOUBLE PRECISION DX,DY,DZ,D,PBC
	 INTEGER I,J,N,DM
	 COMMON/LJ/SIGMA_LJ,EPS,CUTOFF
	 DIMENSION R(N,DM),F(N,DM)
	 
	 POT=0D0
	 F=0D0

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


! Subroutine that simulates a thermal bath according to the Andersen
! algorithm.
	SUBROUTINE ANDERSEN_T(V,N,DM,V_OUT)
	 IMPLICIT NONE
	 DOUBLE PRECISION V,X1,X2,X3,X4,X5,V_OUT
	 DOUBLE PRECISION XOUT1,XOUT2,XOUT3,XOUT4
	 DOUBLE PRECISION SIGMA_T,NU,DT,RHO,L,RHO_2,M,AV,KB
	 INTEGER N,DM,I
	 COMMON/THERMO/NU,SIGMA_T
	 COMMON/CONSTANTS/DT,RHO,L,RHO_2,M,AV,KB
	 DIMENSION V(N,DM),V_OUT(N,DM),RHO_2(4)

	 DO I=1,N
	   CALL RANDOM_NUMBER(X1)
	   CALL RANDOM_NUMBER(X2)
	   CALL RANDOM_NUMBER(X3)
	   CALL RANDOM_NUMBER(X4)
	   CALL RANDOM_NUMBER(X5)
         CALL BOXMULLER(SIGMA_T, 0D0, X1, X2, XOUT1, XOUT2)
         CALL BOXMULLER(SIGMA_T, 0D0, X3, X4, XOUT3, XOUT4)
         IF ((X5.LT.(NU*DT))) THEN
          V_OUT(I,1)=XOUT1
          V_OUT(I,2)=XOUT2
          V_OUT(I,3)=XOUT3
         ELSE
          V_OUT(I,1)=V(I,1)
          V_OUT(I,2)=V(I,2)
          V_OUT(I,3)=V(I,3)
         ENDIF
	 ENDDO
	 
	RETURN
	END

! Integrator that uses the velocity of the current timestep to compute
! the positions. In each step, the calculation of the next velocity
! is needed. 
	SUBROUTINE VERLET_EV(R,V,N,DM,TS,R_F,V_F)                                
	 IMPLICIT NONE
	 DOUBLE PRECISION T,R,F,ROLD,ROLDAUX,V,L,CUTOFF,POT
	 DOUBLE PRECISION DT,KIN,E_TOT,U,D,NU,SIGMA_R
	 DOUBLE PRECISION RHO,V_OUT,PBC,T_INST,KE,MOM,RHO_2,KB
	 DOUBLE PRECISION R_F,V_F,M,AV
	 INTEGER I,J,K,TS,N,DM
	 COMMON/CONSTANTS/DT,RHO,L,RHO_2,M,AV,KB
	 DIMENSION R(N,DM),V(N,DM),F(N,DM),ROLD(N,DM),ROLDAUX(N,DM)
	 DIMENSION KIN(N,TS),E_TOT(TS),U(TS),T_INST(TS),KE(TS),MOM(TS)
	 DIMENSION V_OUT(N,DM),RHO_2(4)
	 DIMENSION R_F(N,DM),V_F(N,DM)


	 T=0D0
	 ROLD=R
	 DO I=1,TS
	  IF (MOD(I,100).EQ.0) THEN
	   WRITE(*,*)'I= ',I,'   TS= ', TS
	  ENDIF
	  T=T+DT
	  CALL EN_LJ(R,N,DM,L,F,POT)

	  ROLDAUX=R
	  DO J=1,N
	   R(J,1)=PBC(R(J,1)+V(J,1)*DT+F(J,1)*DT*DT/(2D0),L)
	   R(J,2)=PBC(R(J,2)+V(J,2)*DT+F(J,2)*DT*DT/(2D0),L)
	   R(J,3)=PBC(R(J,3)+V(J,3)*DT+F(J,3)*DT*DT/(2D0),L)

	   V(J,1)=V(J,1)+F(J,1)*DT/(2D0)
	   V(J,2)=V(J,2)+F(J,2)*DT/(2D0)
	   V(J,3)=V(J,3)+F(J,3)*DT/(2D0)
	  ENDDO

	  CALL EN_LJ(R,N,DM,L,F,POT)

	  DO J=1,N
	   V(J,1)=V(J,1)+F(J,1)*DT/(2D0)
	   V(J,2)=V(J,2)+F(J,2)*DT/(2D0)
	   V(J,3)=V(J,3)+F(J,3)*DT/(2D0)
	  ENDDO

	  CALL ANDERSEN_T(V,N,DM,V_OUT)
	  V=V_OUT
	  	  
	  ROLD=ROLDAUX

	  IF (I.EQ.TS) THEN
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

! Integrator that uses the velocity of the current timestep to compute
! the positions. In each step, the calculation of the next velocity
! is needed. 
	SUBROUTINE VERLET_V(R,V,N,DM,TS,U,KE,E_TOT,T_INST,PR)
	 IMPLICIT NONE
	 DOUBLE PRECISION T,R,F,ROLD,ROLDAUX,V,L,CUTOFF,POT
	 DOUBLE PRECISION DT,KIN,E_TOT,U,D,NU,SIGMA_R
	 DOUBLE PRECISION RHO,V_OUT,PBC,T_INST,KE,RHO_2,M,AV,KB
	 DOUBLE PRECISION DX,DY,DZ,FX,FY,FZ,PR
	 INTEGER I,J,K,TS,N,DM
	 COMMON/CONSTANTS/DT,RHO,L,RHO_2,M,AV,KB
	 DIMENSION R(N,DM),V(N,DM),F(N,DM),ROLD(N,DM),ROLDAUX(N,DM)
	 DIMENSION KIN(N,TS)
	 DIMENSION E_TOT(TS),U(TS),T_INST(TS),KE(TS),PR(TS)
	 DIMENSION V_OUT(N,DM),RHO_2(4)
	 	
	 T=0D0
	 ROLD=R
	 DO I=1,TS
	  IF (MOD(I,100).EQ.0) THEN
	   WRITE(*,*)'I= ',I,'   TS= ', TS
	  ENDIF
	  T=T+DT
	  KE(I)=0D0
	  PR(I)=0D0
	  CALL EN_LJ(R,N,DM,L,F,POT)

	  ROLDAUX=R
	  DO J=1,N
	   R(J,1)=PBC(R(J,1)+V(J,1)*DT+F(J,1)*DT*DT/(2D0),L)
	   R(J,2)=PBC(R(J,2)+V(J,2)*DT+F(J,2)*DT*DT/(2D0),L)
	   R(J,3)=PBC(R(J,3)+V(J,3)*DT+F(J,3)*DT*DT/(2D0),L)

	   V(J,1)=V(J,1)+F(J,1)*DT/(2D0)
	   V(J,2)=V(J,2)+F(J,2)*DT/(2D0)
	   V(J,3)=V(J,3)+F(J,3)*DT/(2D0)
	  ENDDO

	  CALL EN_LJ(R,N,DM,L,F,POT)

	  U(I)=POT
	  
	  DO J=1,N
	   V(J,1)=V(J,1)+F(J,1)*DT/(2D0)
	   V(J,2)=V(J,2)+F(J,2)*DT/(2D0)
	   V(J,3)=V(J,3)+F(J,3)*DT/(2D0)
	  ENDDO
	  DO J=1,N
	   DO K=J+1,N
	   DX=R(J,1)-R(K,1)
	   DY=R(J,2)-R(K,2)
	   DZ=R(J,3)-R(K,3)
	   
	   FX=F(J,1)-F(K,1)
	   FY=F(J,2)-F(K,2)
	   FZ=F(J,3)-F(K,3)

         PR(I)=PR(I)+1D0/(3D0*L**3D0)*(ABS(FX*DX)+ABS(FY*DY)+ABS(FZ*DZ))
	   ENDDO
	  ENDDO
	  
	  CALL ANDERSEN_T(V,N,DM,V_OUT)
	  V=V_OUT
	  	  
	  DO J=1,N
	   KIN(J,I)=0.5D0*(V(J,1)**2D0+V(J,2)**2D0+V(J,3)**2D0)
	   KE(I)=KE(I)+KIN(J,I)
	  ENDDO	

	  T_INST(I)=2D0/DBLE(3*N-3)*KE(I)
	  PR(I)=PR(I)+RHO*T_INST(I)
	  ROLD=ROLDAUX
	  
	 E_TOT(I)=KE(I)+U(I)

	 ENDDO

	RETURN
	END
	
! Subroutine that computes the second part of the final exercise.
	SUBROUTINE APARTAT_2(N,DM,TS)
	 IMPLICIT NONE
	 DOUBLE PRECISION SIGMA_LJ,EPS,CUTOFF,DT,RHO,L,NU,SIGMA_T
	 DOUBLE PRECISION R,V,T,V_X,V_Y,V_Z,U,KE,E_TOT,T_INST,TEMP,PR
	 DOUBLE PRECISION RHO_2,KE_ESP,U_ESP,TOT_ESP,TEMP_ESP,PR_ESP
	 DOUBLE PRECISION R_F,V_F,M,RHO_CONV,AV,KB,CONV_T,CONV_P
	 DOUBLE PRECISION T1,T2,CONV_TEMPS
	 INTEGER N,DM,TS,I,J,TS_EV,MARGE
	 COMMON/LJ/SIGMA_LJ,EPS,CUTOFF
	 COMMON/CONSTANTS/DT,RHO,L,RHO_2,M,AV,KB
	 COMMON/THERMO/NU,SIGMA_T
	 DIMENSION R(N,DM),V(N,DM)
	 DIMENSION E_TOT(TS),U(TS),T_INST(TS),KE(TS),PR(TS)
	 DIMENSION RHO_2(4)
	 DIMENSION R_F(N,DM),V_F(N,DM)


	 CALL INIT_SCC(N,DM,R)
	 V=0D0
	 T1=100D0
	 SIGMA_T=DSQRT(T1)
	 PRINT*,'INITIAL CONF, TEMP=',T1
	 TS_EV=10000
	 CALL VERLET_EV(R,V,N,DM,TS_EV,R_F,V_F)

	 T2=2D0
	 SIGMA_T=DSQRT(T2)
	 PRINT*,'MAIN SYSTEM, TEMP=',T2
	 CALL VERLET_V(R_F,V_F,N,DM,TS,U,KE,E_TOT,T_INST,PR)
	 
	 CONV_T=EPS*1D3/(AV*KB)
 	 RHO_CONV=RHO*M/(AV*(SIGMA_LJ*1D-8)**3D0)
 	 CONV_P=EPS*1D3/(AV*(SIGMA_LJ*1D-10)**3D0)
 	 CONV_TEMPS=DSQRT(M*1D-3/(EPS*1D3))*SIGMA_LJ*1D-10

	 T=0D0
	 DO I=1,TS
	  T=T+DT
	  WRITE(12,100)T*CONV_TEMPS,KE(I)*EPS,U(I)*EPS,E_TOT(I)*EPS
     &  ,T_INST(I)*CONV_T,PR(I)*CONV_P
100	  FORMAT(E20.12,4(F20.12),E20.12)
	 ENDDO

	 KE_ESP=0D0
	 U_ESP=0D0
	 TOT_ESP=0D0
	 TEMP_ESP=0D0
	 PR_ESP=0D0
	  
	 MARGE=200000
	 DO I=MARGE,TS
	  KE_ESP=KE_ESP+KE(I)*EPS/DBLE(TS-MARGE)
	  U_ESP=U_ESP+U(I)*EPS/DBLE(TS-MARGE)
	  TOT_ESP=TOT_ESP+E_TOT(I)*EPS/DBLE(TS-MARGE)
	  TEMP_ESP=TEMP_ESP+T_INST(I)*CONV_T/DBLE(TS-MARGE)
	  PR_ESP=PR_ESP+(PR(I)*CONV_P)/DBLE(TS-MARGE)
	 ENDDO
	 
	 WRITE(*,*)'T_FINAL= ', T_INST(TS)

	 WRITE(13,104)RHO_CONV,KE_ESP,U_ESP,TOT_ESP,TEMP_ESP,PR_ESP
104	 FORMAT(5(F20.12),E20.12)

	RETURN
	END
	


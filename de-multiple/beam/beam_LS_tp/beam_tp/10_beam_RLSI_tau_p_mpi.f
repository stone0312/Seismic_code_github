
*==================================================================*
*                                                                  *
*              Beam Forming 2D RLSI                                * 
*                                                                  *
*==================================================================*
*              ALL   RIGHTS   RESERVED                             *
*==================================================================*
*       Departmelt of Marine Geology and Geophysics , TongJi       *
*       University                                                 *
*==================================================================*
*       Remarks:                                                   *
*        Author:                                                   *
*                 Jiangtao Hu                                      *
*                 08, 2012                                         *
*==================================================================*
*       Summary:                                                   *
*             This program do a beam forming with Least square     *
!          inversion with Radon constraint to subtract multiple    *
!       note:input must be common shot gather                      *
*==================================================================*

      INCLUDE 'mpif.h'
      INTEGER MYID,NPROC,IERR
      INTEGER  status(MPI_STATUS_SIZE)
      PARAMETER(Lbyte=1)
      parameter(PAI2=2.0*3.1415926)

      INTEGER LT
      INTEGER NW1,NW2,NW3,NW4,NW
      INTEGER NP,NP_h,Nbeam,Nbeam_h
      REAL    DX,DT
      REAL    THETA_MAX, Vmin, Pmax
      REAL    DP,DW
      REAL    F1,F2,F3,F4

      CHARACTER*256 FN1,FN2,FN3

!==================================================================
       CALL MPI_INIT( ierr )
       CALL MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
       CALL MPI_COMM_SIZE( MPI_COMM_WORLD, NPROC, ierr )
!=================================================================
      CALL READPAR(FN1, FN2, FN3, LT, DT, DX, THETA_MAX, Vmin, 
     +   F1, F2, F3, F4,
     +   Nbeam_h, NP_h, NS_START, NS_END)

      Nbeam=2*Nbeam_h+1
      NP=2*NP_h+1
!      NSHOT=NS_END-NS_START+1

      Dw=PAI2*1000.0/(LT*DT)
      W1=PAI2*F1
      Nw1=W1/Dw+0.5
      W2=PAI2*F2
      Nw2=W2/Dw+0.5
      W3=PAI2*F3
      Nw3=W3/Dw+0.5
      W4=PAI2*F4
      Nw4=W4/Dw+0.5
      NW=NW4-NW1+1

      Pmax=sin(theta_max*PAI2/360.)/Vmin
      Dp=2*Pmax/(NP-1)

      print*,'Pmax,Dp=',Pmax,Dp
      print*,'dx,dt=',dx,dt
      print*,'nw1,nw4,dw=',nw1,nw4,dw

      PRINT*,'BEAM FORMING BEGIN!!!!'

      CALL BEAM_FORMING(FN1, FN2, FN3, LT, DT, DX, 
     +  NW1, NW2, NW3, NW4, NW,
     +  DP, DW, Nbeam_h, NP_h, Nbeam, NP,
     +  NS_START, NS_END, NSHOT,
     +  Lbyte,
     +  MYID, NPROC, ierr)

      PRINT*,'BEAM FORMING END!!!!!'
      call MPI_FINALIZE(ierr)

      END

!==================================================================

      SUBROUTINE BEAM_FORMING(FN1, FN2, FN3, LT, DT, DX, 
     +  NW1, NW2, NW3, NW4, NW,
     +  DP, DW, Nbeam_h, NP_h, Nbeam, NP,
     +  NS_START, NS_END, NSHOT,
     +  Lbyte,
     +  MYID, NPROC, ierr)

      INCLUDE 'mpif.h'
      INTEGER MYID, NPROC, ierr

      INTEGER LT,Lbyte
      INTEGER NW1, NW2, NW3, NW4, NW
      INTEGER Nbeam_h, NP_h, Nbeam, NP
      INTEGER NS_START, NS_END, NSHOT
      REAL    DX
      REAL    DP, DW, DT

      INTEGER,ALLOCATABLE :: INDEX_TABLE(:,:)
      REAL,ALLOCATABLE :: BUFF(:)
      
      CHARACTER*256 FN1, FN2, FN3

!=========definition of index table=====
!1:ishot
!3:irec
!4:ntrace
!=========================================

      CALL SEARCH_CORD(FN1, LT, Lbyte, NSHOT)

      PRINT*,'NSHOT=',NSHOT

      ALLOCATE(INDEX_TABLE(3,NSHOT),STAT=IERR)
      IF(IERR.NE.0.0)THEN
      PRINT*,'ALLOCATE ERR'
      STOP
      ENDIF

      CALL SEARCH_AND_INDEX(INDEX_TABLE, FN1, LT, Lbyte, NSHOT)

      PRINT*,'INDEX END'
      do ishot=1,3
        print*,INDEX_TABLE(1,ishot),INDEX_TABLE(2,ishot),
     +INDEX_TABLE(3,ishot)
      enddo

      OPEN(11,FILE=FN1,ACCESS='DIRECT',RECL=Lbyte*(LT+60))
      OPEN(22,FILE=FN2,ACCESS='DIRECT',RECL=Lbyte*(LT+60))
      OPEN(33,FILE=FN3,ACCESS='DIRECT',RECL=Lbyte*(LT+60),
     +    STATUS='REPLACE')

      IRECP=1

      PRINT*,'RLSI BEGIN!!!!!!!'

      DO 6666 ISHOT=1+MYID,NSHOT , NPROC
			  print*,ISHOT
        !ISHOT=1
        NX=INDEX_TABLE(3, ISHOT)

        Nbeam_ALL=NX

        L1 = 1
        L2 = L1+LT*NX                   !CO_ALL(LT,NX)
        L3 = L2+LT*NX                   !CO_ALL_M(LT,NX)
        L4 = L3+LT*Nbeam                !CO_LOCAL(LT,Nbeam)
        L5 = L4+LT*Nbeam                !CO_LOCAL_M(LT,Nbeam)
        L6 = L5+LT*NP                   !CO_TP_LOCAL(LT,NP)
        L7 = L6+LT*NX                   !CO_ALL_R(LT,NX)
        L8 = L7+LT*NP                   !CO_TP_SEMB(LT,NP)
        L9 = L8+LT*NP                   !CO_TP_SEMB_M(LT,NP)
        L10= L9+2*LT                    !CTRACE(LT)
        L11= L10+2*LT*Nbeam             !CCO_LOCAL(LT,Nbeam)
        L12= L11+2*LT*NP                !CCO_TP_LOCAL(LT,NP)
        L13= L12+2*Nbeam*NP             !A
        L14= L13+2*NP*Nbeam             !AH
        L15= L14+2*NP*NP                !AHA
        L16= L15+2*NP                   !F
        L17= L16+2*NP                   !X1
        L18= L17+2*NP                   !X2
        L19= L18+2*NP                   !R1
        L20= L19+2*NP                   !R2
        L21= L20+2*NP                   !P1
        L22= L21+2*NP                   !P2
        L23= L22+2*NP                   !AP
        L24= L23+NP                     !ENERGY
        L25= L24+60*NX                  !HEAD_ALL
        L26= L25+NP*LT                   !SMOOTHING SUBS
	                                      !+SEMBLING:CO_TP_SEMB_SM
        PRINT*,'MEMORY USED:',REAL(L26)/256.0/1024.0,'Mb'

        ALLOCATE(BUFF(L26),STAT=IERR)
        IF(IERR.NE.0.0)THEN
        PRINT*,'ALLOCATE ERR'
        STOP
        ENDIF

        CALL BEAMING(BUFF(L1), BUFF(L2), BUFF(L3), BUFF(L4), BUFF(L5),
     +         BUFF(L6), BUFF(L7), BUFF(L8), BUFF(L9), BUFF(L10),
     +         BUFF(L11), BUFF(L12), BUFF(L13), BUFF(L14), BUFF(L15),
     +         BUFF(L16), BUFF(L17), BUFF(L18), BUFF(L19), BUFF(L20),
     +         BUFF(L21), BUFF(L22), BUFF(L23), BUFF(L24), BUFF(L25),
     +         INDEX_TABLE,
     +         LT, DT, NX, DX,
     +         NW1, NW2, NW3, NW4, NW,
     +         DP, DW, Nbeam_h, NP_h, Nbeam, NP, Nbeam_ALL,
     +         ISHOT, NSHOT,
     +         IRECP)

      DEALLOCATE(BUFF)

6666  CONTINUE

      CLOSE(11)
      CLOSE(22)
      CLOSE(33)

      DEALLOCATE(INDEX_TABLE)

      END

!================================================================

      SUBROUTINE BEAMING(CO_ALL, CO_ALL_M, CO_LOCAL, CO_LOCAL_M, 
     +      CO_TP_LOCAL, CO_ALL_R, CO_TP_SEMB, CO_TP_SEMB_M,
     +      CTRACE, CCO_LOCAL, CCO_TP_LOCAL, A, AH, AHA, F, X1, X2,
     +      R1, R2, P1, P2, AP, ENERGY, HEAD_ALL, CO_TP_SEMB_SM, 
     +      INDEX_TABLE,
     +      LT, DT, NX, DX,
     +      NW1, NW2, NW3, NW4, NW,
     +      DP, DW, Nbeam_h, NP_h, Nbeam, NP, Nbeam_ALL,
     +      ISHOT, NSHOT,
     +      IRECP)

      include'fftw3.f'
      integer*8 plan_f,plan_b

      INTEGER LT, Lbyte
      INTEGER ISHOT, NSHOT
      INTEGER NW1, NW2, NW3, NW4, NW
      INTEGER Nbeam_h, NP_h, Nbeam, NP, Nbeam_ALL
      INTEGER NX
      INTEGER IRECP
      REAL    DT, DX, DP, DW

      REAL    CO_ALL(LT,NX), CO_LOCAL(LT,Nbeam)
      REAL    CO_ALL_M(LT,NX), CO_LOCAL_M(LT,Nbeam)
      REAL    CO_TP_LOCAL(LT,NP), CO_ALL_R(LT,NX)
      REAL    CO_TP_SEMB(LT,NP), CO_TP_SEMB_M(LT,NP)
			REAL    CO_TP_SEMB_SM(LT,NP)
      REAL    ENERGY(NP)

      INTEGER   INDEX_TABLE(3, NSHOT)
      DIMENSION HEAD_ALL(60, Nx)

      COMPLEX CTRACE(LT)
      COMPLEX CCO_LOCAL(LT,Nbeam), CCO_TP_LOCAL(LT,NP)
      COMPLEX A(Nbeam,NP),AH(NP,Nbeam),AHA(NP,NP)
      COMPLEX F(NP),X1(NP),X2(NP)
      COMPLEX P1(NP),P2(NP)
      COMPLEX R1(NP),R2(NP)
      COMPLEX AP(NP)

      call sfftw_plan_dft_1d(plan_f,lt,ctrace,ctrace,
     +  FFTW_BACKWARD,FFTW_MEASURE)
      call sfftw_plan_dft_1d(plan_b,lt,ctrace,ctrace,
     +  FFTW_FORWARD,FFTW_MEASURE)

!========input common shot gather====================

      IREC=INDEX_TABLE(2, ISHOT)
      NTRACE=INDEX_TABLE(3, ISHOT)

			DO IX=1,NX
        DO IT=1,LT
          CO_ALL_R(IT,IX)=0.0
        ENDDO
      ENDDO
!      print*,nx
      CALL INPUT_GATHER(CO_ALL,CO_ALL_M,HEAD_ALL,IREC,NX,LT)
!      open(55,file='CO_ALL',access='direct',recl=lt)
!      do ix=1,nx
!        write(55,rec=ix)(CO_ALL(it,ix),it=1,lt)
!      enddo
!      close(55)
!      open(55,file='CO_ALL_M',access='direct',recl=lt)
!      do ix=1,nx
!        write(55,rec=ix)(CO_ALL_M(it,ix),it=1,lt)
!      enddo
!      close(55)
!			stop

!      DO 6666 Ibeam=1,Nbeam_ALL

      NX_H=0
      NX_ALL=2*NX_H

!				DO 8888 icenter=5,5
        DO 8888 icenter=1,NX   !NX_H+1, Nbeam_ALL, NX_ALL
!        icenter=396
        DO IX=1,Nbeam
          DO IT=1,LT
            CO_LOCAL(IT,IX)=0.0
            CO_LOCAL_M(IT,IX)=0.0
          ENDDO
        ENDDO

        IRECCC=icenter-Nbeam_h
				DO IX=1,Nbeam
          IF(IRECCC.GE.1.AND.IRECCC.LE.NX)THEN
          DO IT=1,LT
             CO_LOCAL(IT,IX)=CO_ALL(IT,IRECCC)
             CO_LOCAL_M(IT,IX)=CO_ALL_M(IT,IRECCC)
          ENDDO
          ENDIF
          IRECCC=IRECCC+1
        ENDDO
!      open(55,file='CO_LOCAL',access='direct',recl=lt)
!      do ix=1,Nbeam
!        write(55,rec=ix)(CO_LOCAL(it,ix),it=1,lt)
!      enddo
!      close(55)
!      open(55,file='CO_LOCAL_M',access='direct',recl=lt)
!      do ix=1,Nbeam
!        write(55,rec=ix)(CO_LOCAL_M(it,ix),it=1,lt)
!      enddo
!      close(55)

        CALL RLSI_BEAM(CO_LOCAL, CO_TP_LOCAL,
     +      CTRACE, CCO_LOCAL, CCO_TP_LOCAL, A, AH, AHA, F, X1, X2,
     +      R1, R2, P1, P2, AP, ENERGY,
     +      LT, DT, DX,
     +      NW1, NW2, NW3, NW4, NW,
     +      DP, DW, Nbeam_h, NP_h, Nbeam, NP,
     +      plan_f,plan_b)

        CALL SEMB(CO_LOCAL, CO_TP_SEMB, 
     +      LT, DT, DX, Nbeam_h, Nbeam, 
     +      DP, NP_h, NP)
        CALL SEMB(CO_LOCAL_M, CO_TP_SEMB_M, 
     +      LT, DT, DX, Nbeam_h, Nbeam, 
     +      DP, NP_h, NP)

!inversion forward tau-p result
!			open(55,file='thre01_CO_TP_LOCAL_shot1_icenter200',
!     + access='direct',recl=lt)
!      do ix=1,NP
!        write(55,rec=ix)(CO_TP_LOCAL(it,ix),it=1,lt)
!      enddo
!      close(55)

!initial sembling cofficients 			
!			open(55,file='thre01_CO_TP_SEMB',
!     + access='direct',recl=lt)
!      do ix=1,NP
!        write(55,rec=ix)(CO_TP_SEMB(it,ix),it=1,lt)
!      enddo
!      close(55)
!			open(55,file='thre01_CO_TP_SEMB_M',
!     +access='direct',recl=lt)
!      do ix=1,NP
!        write(55,rec=ix)(CO_TP_SEMB_M(it,ix),it=1,lt)
!			enddo
!      close(55)			

      DO IP=1,NP
			 DO IT=1,LT
			   if(CO_TP_SEMB(IT,IP).GE.0.01)THEN
					 CO_TP_SEMB(IT,IP)=1.0
				 else 
					 CO_TP_SEMB(IT,IP)=0.0
			   ENDIF
		    ENDDO
			ENDDO
		  DO IP=1,NP
	        DO IT=1,LT
						if(CO_TP_SEMB_M(IT,IP).GE.0.01)THEN
		           CO_TP_SEMB_M(IT,IP)=1.0
						else
							 CO_TP_SEMB_M(IT,IP)=0.0
					  ENDIF	
		 	  	ENDDO
		    ENDDO
										 
      DO IP=1,NP
        DO IT=1,LT
          CO_TP_SEMB(IT,IP)=CO_TP_SEMB(IT,IP)-CO_TP_SEMB(IT,IP)
     +     *CO_TP_SEMB_M(IT,IP)
!          CO_TP_LOCAL(IT,IP)=CO_TP_LOCAL(IT,IP)*CO_TP_SEMB(IT,IP)
        ENDDO
      ENDDO
!initial substracted sembling cofficients
!			open(55,file='thre01_CO_TP_SEMB_INI_SUBSTRACTED_10',
!     + access='direct',recl=lt)
!      do ix=1,NP
!        write(55,rec=ix)(CO_TP_SEMB(it,ix),it=1,lt)
!      enddo
!      close(55)

!smoothing the initial substracted sembling cofficients
!      CALL SM_SUBS_SEM(CO_TP_SEMB_SM,CO_TP_SEMB,NP,LT)    		 	
!		  DO IP=1,NP 
!		    DO IT=1,LT 
!           CO_TP_LOCAL(IT,IP)=CO_TP_LOCAL(IT,IP)*CO_TP_SEMB_SM(IT,IP)
! 		 	  ENDDO  
!		  ENDDO   
!smoothing substracted sembling cofficients
!		  open(55,file='thre01_CO_TP_SUBS_SEMB_SM_shot_1_icenter200',
!     + access='direct',recl=lt)
!			do ix=1,NP
!					write(55,rec=ix)(CO_TP_SEMB_SM(it,ix),it=1,lt)
!			enddo
!			close(55)
													
!CO_TP_LOCAL:SUBSTRACTED FORWARD TAU-P RESULT TIME DOMAIN				
!  		open(55,file='thre01_CO_TP_LOCAL_SM_shot_1_icenter200',
!     + access='direct',recl=lt)
!      do ix=1,NP
!        write(55,rec=ix)(CO_TP_LOCAL(it,ix),it=1,lt)
!      enddo
!      close(55)
!			stop
!CCO_TP_LOCAL:SUBSTRACTED FORWARD TAU-P RESULT FREQUENCY DOMAIN
			
        CALL BACK_TRANSFORM(CO_TP_LOCAL, CO_LOCAL, 
     +       CCO_TP_LOCAL, CCO_LOCAL, A, CTRACE, 
     +       LT, DT, DX, 
     +       NW1, NW2, NW3, NW4, NW,
     +       DP, DW, Nbeam_h, NP_h, Nbeam, NP,
     +       plan_f,plan_b)

        ilocal=Nbeam_h+1-NX_H
				DO IX=icenter-NX_H, icenter+NX_H
          IF(IX.GE.1.AND.IX.LE.NX)THEN
          DO IT=1,LT
             CO_ALL_R(IT,IX)=CO_LOCAL(IT,ilocal)
          ENDDO
!          print*,ix,ilocal
          ilocal=ilocal+1
          ENDIF
        ENDDO
!      open(55,file='CO_LOCAL1',access='direct',recl=lt)
!      do ix=1,Nbeam
!        write(55,rec=ix)(CO_LOCAL(it,ix),it=1,lt)
!      enddo
!      close(55)
!      open(55,file='CO_ALL_R',access='direct',recl=lt)
!      do ix=1,11
!        write(55,rec=ix)(CO_ALL_R(it,ix+94),it=1,lt)
!      enddo
!      close(55)
!			stop

8888  CONTINUE
!6666  CONTINUE

!=======================================================
!      IREC=1
      IREC=INDEX_TABLE(2, ISHOT)
      DO IX=1,NX
        WRITE(33,REC=IREC)(HEAD_ALL(IH,IX),IH=1,60),
     +  (CO_ALL_R(IT,IX),IT=1,LT)
        IREC=IREC+1
      ENDDO

      call sfftw_destroy_plan(plan_f)
      call sfftw_destroy_plan(plan_b)


      END

!==============================================================

      SUBROUTINE BACK_TRANSFORM(CO_TP_LOCAL, CO_LOCAL, 
     +       CCO_TP_LOCAL, CCO_LOCAL, A, CTRACE, 
     +       LT, DT, DX,
     +       NW1, NW2, NW3, NW4, NW,
     +       DP, DW, Nbeam_h, NP_h, Nbeam, NP,
     +       plan_f,plan_b)

      include'fftw3.f'
      integer*8 plan_f,plan_b

      PARAMETER(PAI2=2.0*3.1415926)

      REAL    DT, DX, DP, DW
      INTEGER LT, NW1, NW2, NW3, NW4, NW
      INTEGER Nbeam_h, NP_h, Nbeam, NP

      REAL    CO_TP_LOCAL(LT,NP), CO_LOCAL(LT,Nbeam)
      COMPLEX CTRACE(LT), A(Nbeam, NP)
      COMPLEX CCO_TP_LOCAL(LT,NP), CCO_LOCAL(LT,Nbeam)
      COMPLEX sum1

!=====TRANSFORM TO FREQUENCY DOMAIN
      DO IP=1,NP

         DO IT=1,LT
           CTRACE(IT)=CMPLX(CO_TP_LOCAL(IT,IP),0.0)
         ENDDO

        CALL sfftw_execute(plan_b)
        DO IT=1,LT
           CTRACE(IT)=CTRACE(IT)/real(LT)
        ENDDO

        DO IT=1,LT
           CCO_TP_LOCAL(IT,IP)=CTRACE(IT)
        ENDDO

      ENDDO
         
!=====back transform
      
      do 2222 iw=Nw1+1,Nw4-1

        w=(iw-1)*dw

         do ip=1,np
           p=(ip-NP_h-1)*dp
           do ix=1,Nbeam
             x=(ix-Nbeam_h-1)*dx
             phase=w*p*x
             coe=1.0/sqrt(real(Nbeam))
             A(ix,ip)=coe*cmplx(cos(phase),sin(phase))
           enddo
         enddo

         do ix=1,Nbeam
           sum1=0.0
           do ip=1,np
             sum1=sum1+A(ix,ip)*CCO_TP_LOCAL(iw,ip)
           enddo
           CCO_LOCAL(iw,ix)=sum1
           IF(IW.GT.1.AND.IW.NE.LT+2-IW)THEN
             CCO_LOCAL(LT+2-IW,ix)=CONJG(CCO_LOCAL(iw,ix))
           ENDIF
         enddo

2222  continue

!=====back fft

      DO IX=1, Nbeam
         DO IT=1,LT
            CTRACE(IT)=CCO_LOCAL(IT, IX)
         ENDDO
				 call sfftw_execute(plan_f)
         DO IT=1,LT
            CO_LOCAL(IT, IX)=REAL(CTRACE(IT))
         ENDDO
      ENDDO

			END
!====================================

			SUBROUTINE SEMB(CO_LOCAL, CO_TP_SEMB, 
     +      LT, DT, DX, Nbeam_h, Nbeam, 
     +      DP, NP_h, NP)

      PARAMETER(threshold=0.1)
      INTEGER LT, Nbeam_h, Nbeam, NP_h, NP
      REAL    DT, DX, DP
      REAL    CO_LOCAL(LT, Nbeam), CO_TP_SEMB(LT, NP)


!====compute semb
			DO 111 IT=1,LT
        DO 222 IP=1,NP
         SUM1=0.0
         SUM2=0.0
         P=(IP-NP_h-1)*DP
         DO 333 IX=1, Nbeam
           X=(IX-Nbeam_h-1)*DX
           T=-X*P
           NTT=IT+1000.0*T/DT
           IF(NTT.GE.1.AND.NTT.LE.LT)THEN
             SUM1=SUM1+CO_LOCAL(NTT,IX)
             SUM2=SUM2+CO_LOCAL(NTT,IX)*CO_LOCAL(NTT,IX)
           ENDIF
333   CONTINUE
        CO_TP_SEMB(IT,IP)=(SUM1*SUM1)/(SUM2+0.001)
        CO_TP_SEMB(IT,IP)=CO_TP_SEMB(IT,IP)/REAL(Nbeam)
222   CONTINUE				 
111   CONTINUE

!=====filter

      rmax=0.0
      DO IP=1,NP
         DO IT=1,LT
           IF(CO_TP_SEMB(IT,IP).LT.threshold)CO_TP_SEMB(IT,IP)=0.0
           IF(CO_TP_SEMB(IT,IP).GT.rmax)rmax=CO_TP_SEMB(IT,IP)
         ENDDO
      ENDDO

      END
			
!=============================================================			
			SUBROUTINE SM_SUBS_SEM(CO_TP_SEMB_SM,CO_TP_SEMB,NP,LT)
      INTEGER NP,LT
      REAL    CO_TP_SEMB_SM(LT,NP),CO_TP_SEMB(LT,NP)	
	    DO IP=1,NP
	      DO IT=1,LT
				   CO_TP_SEMB_SM(IT,IP)=0.0
				ENDDO
			ENDDO	 
!=====smooth

			DO IP=1,NP
        DO IT=1,LT
          sum1=0.0
          do ipp=ip-4, ip+4
            if(ipp.ge.1.and.ipp.le.np)then
            do itt=it-4, it+4
              if(itt.ge.1.and.itt.le.lt)then
              sum1=sum1+CO_TP_SEMB(itt,ipp)
              endif
            enddo
            endif
          enddo
          CO_TP_SEMB_SM(IT,IP)=sum1/81.0
        ENDDO
      ENDDO

!=====normalization

      rmax1=0.0
      DO IP=1,NP
        DO IT=1,LT
         IF(CO_TP_SEMB(IT,IP).GT.rmax1)rmax1=CO_TP_SEMB(IT,IP)
        ENDDO
      ENDDO

      DO IP=1,NP
        DO IT=1,LT
          CO_TP_SEMB(IT,IP)=rmax*CO_TP_SEMB(IT,IP)/rmax1
        ENDDO
      ENDDO	

      END

!==============================================================

      SUBROUTINE RLSI_BEAM(CO_LOCAL, CO_TP_LOCAL,
     +      CTRACE, CCO_LOCAL, CCO_TP_LOCAL, A, AH, AHA, F, X1, X2,
     +      R1, R2, P1, P2, AP, ENERGY,
     +      LT, DT, DX,
     +      NW1, NW2, NW3, NW4, NW,
     +      DP, DW, Nbeam_h, NP_h, Nbeam, NP,
     +      plan_f,plan_b)

      include'fftw3.f'
      integer*8 plan_f,plan_b

      INTEGER LT, Lbyte
      INTEGER Ioffset
      INTEGER NW1, NW2, NW3, NW4, NW
      INTEGER Nbeam_h, NP_h, Nbeam, NP, Nbeam_ALL
      REAL    DT, DX, DP, DW
      REAL    Lamda

      REAL    CO_LOCAL(LT,Nbeam)
      REAL    CO_TP_LOCAL(LT,NP)
      REAL    ENERGY(NP)

      COMPLEX CTRACE(LT)
      COMPLEX CCO_LOCAL(LT,Nbeam), CCO_TP_LOCAL(LT,NP)
      COMPLEX A(Nbeam,NP),AH(NP,Nbeam),AHA(NP,NP)
      COMPLEX F(NP),X1(NP),X2(NP)
      COMPLEX P1(NP),P2(NP)
      COMPLEX R1(NP),R2(NP)
      COMPLEX AP(NP)

      COMPLEX SUM1,SUM2

      !transform time domain to frequency domain

      DO 1111 IX=1, Nbeam

        DO IT=1,LT
          CTRACE(IT)=CMPLX(CO_LOCAL(IT,IX),0.0)
        ENDDO

        CALL sfftw_execute(plan_b)
        DO IT=1,LT
         CTRACE(IT)=CTRACE(IT)/real(LT)
        ENDDO

        call HAMMING_WINDOW(CTRACE, NW1, NW2, NW3, NW4, LT)

        DO IT=1,LT
          CCO_LOCAL(IT,IX)=CTRACE(IT)
        ENDDO

1111  CONTINUE

      !get the max energy for constrain

      do 6666 iw=Nw1+1,Nw4-1
         w=(iw-1)*dw

         do ip=1,np
           p=(ip-NP_h-1)*dp
           do ix=1,Nbeam
             x=(ix-Nbeam_h-1)*dx
             phase=w*p*x
             coe=1.0/sqrt(real(Nbeam))
             A(ix,ip)=coe*cmplx(cos(phase),sin(phase))
           enddo
         enddo

         do ip=1,NP
           do ix=1,Nbeam
             AH(ip,ix)=CONJG(A(ix,ip))
           enddo
         enddo

         do ip=1,np
           sum1=0.0
           do ix=1,Nbeam
             sum1=sum1+AH(ip,ix)*CCO_LOCAL(iw,ix)
           enddo
           f(ip)=sum1
           energy(ip)=energy(ip)+real(f(ip)*CONJG(f(ip)))
         enddo

6666  CONTINUE

      Emax=0.0
      DO IP=1,NP
        IF(energy(IP).GT.Emax)Emax=energy(IP)
        energy(IP)=0.0
      ENDDO

      Lamda=0.1*Emax

      !RLSI Tau-P

      do 2222 iw=Nw1+1,Nw4-1
         w=(iw-1)*dw

         do ip=1,np
           p=(ip-NP_h-1)*dp
           do ix=1,Nbeam
             x=(ix-Nbeam_h-1)*dx
             phase=w*p*x
             coe=1.0/sqrt(real(Nbeam))
             A(ix,ip)=coe*cmplx(cos(phase),sin(phase))
           enddo
         enddo

         do ip=1,NP
           do ix=1,Nbeam
             AH(ip,ix)=CONJG(A(ix,ip))
           enddo
         enddo

         do ip=1,np
           do ip1=1,np
             sum1=0.0
             do ix=1,Nbeam
               sum1=sum1+AH(ip,ix)*A(ix,ip1)
             enddo
             AHA(ip,ip1)=sum1
           enddo
         enddo

         do ip=1,np
           sum1=0.0
           do ix=1,Nbeam
             sum1=sum1+AH(ip,ix)*CCO_LOCAL(iw,ix)
           enddo
           f(ip)=sum1
           energy(ip)=energy(ip)+real(f(ip)*CONJG(f(ip)))
         enddo

         do ip=1,np
!           IF(energy(ip).LT.1.0E-6)energy(ip)=1.0E-6
           AHA(ip,ip)=AHA(ip,ip)+0.1!Lamda/energy(ip)
           AHA(ip,ip)=cmplx(real(AHA(ip,ip)),0.0)
         enddo

         err=1.0E-15
         call complex_CG(AHA,x1,x2,f,P1,P2,R1,R2,AP,np,err)

         do ip=1,np
           CCO_TP_LOCAL(iw,ip)=x2(ip)
           IF(IW.GT.1.AND.IW.NE.LT+2-IW)THEN
             CCO_TP_LOCAL(LT+2-IW,ip)=CONJG(CCO_TP_LOCAL(iw,ip))
           ENDIF
         enddo

2222  continue

      !transform frequency domain to time domain

      do 3333 ip=1,np

        do it=1,lt
          ctrace(it)=CCO_TP_LOCAL(it,ip)
        enddo

        call sfftw_execute(plan_f)

        do it=1,lt
          CO_TP_LOCAL(it,ip)=real(ctrace(it))
        enddo

3333  continue

      END

!==============================================================

      SUBROUTINE INPUT_GATHER(CO_ALL,CO_ALL_M, HEAD_ALL,IREC,NX, LT)

      INTEGER IREC,NX

      REAL    CO_ALL(LT,NX), CO_ALL_M(LT,NX)
      DIMENSION HEAD_ALL(60, NX)

      DIMENSION HEAD(60)

      INTEGER*4 fldr, SX, SY, GX, GY
      EQUIVALENCE (HEAD(3), fldr)
      EQUIVALENCE (HEAD(19), SX), (HEAD(20), SY)
      EQUIVALENCE (HEAD(21), GX), (HEAD(22), GY)

      DO IX=1,NX
           
        READ(11,REC=IREC,ERR=5555)(HEAD(IH), IH=1, 60),
     +     (CO_ALL(IT,IX),IT=1,LT)

        DO IH=1,60
          HEAD_ALL(IH,IX)=HEAD(IH)
        ENDDO
        READ(22,REC=IREC,ERR=5555)(HEAD(IH), IH=1, 60),
     +     (CO_ALL_M(IT,IX),IT=1,LT)

        IREC=IREC+1

      ENDDO

5555  CONTINUE

      END

!================================================================

      SUBROUTINE SEARCH_AND_INDEX(INDEX_TABLE, FN1, LT, Lbyte, NSHOT)

      INTEGER LT,Lbyte
      INTEGER NSHOT
      INTEGER INDEX_TABLE(3,NSHOT)

      DIMENSION HEAD(60)

      INTEGER*4 fldr, SX, SY, GX, GY
      EQUIVALENCE (HEAD(3), fldr)
      EQUIVALENCE (HEAD(19), SX), (HEAD(20), SY)
      EQUIVALENCE (HEAD(21), GX), (HEAD(22), GY)

      CHARACTER*256 FN1

!=========definition of index table=====
!1:ishot
!2:irec
!3:ntrace
!=========================================

      OPEN(11,FILE=FN1,ACCESS='DIRECT',RECL=Lbyte*(60+LT))
      IREC=1

      DO 6666 ISHOT=1,NSHOT

      NTRACE=0

1001  CONTINUE

      READ(11,REC=IREC,ERR=5555)(HEAD(IH), IH=1, 60)

      IF(MOD(IREC,200000).EQ.0)PRINT*,'IREC=',IREC
      IF(IREC.EQ.1)ifldr=fldr

      IF(ifldr.NE.fldr)THEN
        ifldr=fldr
        goto 6666
      ELSE

        IF(NTRACE.EQ.0)INDEX_TABLE(2,ISHOT)=IREC
        INDEX_TABLE(1,ISHOT)=fldr
        NTRACE=NTRACE+1

        INDEX_TABLE(3, ISHOT)=NTRACE

      ENDIF

      IREC=IREC+1
      GOTO 1001

6666  CONTINUE

5555  CONTINUE

      CLOSE(11)

      END
!=================================================================

      SUBROUTINE SEARCH_CORD(FN1, LT, Lbyte, NSHOT)
!count no. of shots
      INTEGER LT,Lbyte
      INTEGER Noffset

      DIMENSION HEAD(60)

      INTEGER*4 fldr, SX, SY, GX, GY
      EQUIVALENCE (HEAD(3), fldr)
      EQUIVALENCE (HEAD(19), SX), (HEAD(20), SY)
      EQUIVALENCE (HEAD(21), GX), (HEAD(22), GY)

      CHARACTER*256 FN1

      OPEN(11,FILE=FN1,ACCESS='DIRECT',RECL=Lbyte*(60+LT))
      IREC=1
      NSHOT=0
			ifldr=0

1001  CONTINUE

      READ(11,REC=IREC,ERR=5555)(HEAD(IH), IH=1, 60)

      IF(MOD(IREC,200000).EQ.0)PRINT*,'IREC=',IREC

      IF(ifldr.NE.fldr)THEN
        NSHOT=NSHOT+1
				ifldr=fldr
      ENDIF

      IREC=IREC+1
      GOTO 1001

5555  CONTINUE

      CLOSE(11)

      END

!==================================================================

      SUBROUTINE READPAR(FN1, FN2, FN3, LT, DT, DX, THETA_MAX, Vmin, 
     +   F1, F2, F3, F4,
     +   Nbeam_h, NP_h, NS_START, NS_END)

      INTEGER LT
      INTEGER NP_h,Nbeam_h
      INTEGER NS_START, NS_END
      REAL    DX,DT
      REAL    THETA_MAX, Vmin
      REAL    F1,F2,F3,F4

      CHARACTER*256 FN1,FN2,FN3

      OPEN(11,FILE='10_beam_RLSI_tau_p_mpi.par')

      READ(11,'(A)')FN1
      WRITE(*,*) 'THE INITIAL SHOT GATHER FILENAME'
      WRITE(*,'(A)') FN1

      READ(11,'(A)')FN2
      WRITE(*,*) 'THE MULTIPLE MODEL SHOT GATHER FILENAME'
      WRITE(*,'(A)') FN2

      READ(11,'(A)')FN3
      WRITE(*,*) 'THE SUBTRACTED SHOT GATHER FILENAME'
      WRITE(*,'(A)') FN3

      READ(11,*) DX
      WRITE(*,*) 'THE TRACE INTERVAL'
      WRITE(*,*) DX

      READ(11,*) LT, DT
      WRITE(*,*) 'THE TRACE LENGTH (SAMPLE POINT NUMBER) AND DT(ms)'
      WRITE(*,*) LT, DT

      READ(11,*) F1,F2,F3,F4
      WRITE(*,*) 'THE FREQUENCY RANGE (Hz)'
      WRITE(*,*) F1,F2,F3,F4

      READ(11,*) THETA_MAX
      WRITE(*,*) 'THE MAX DIP ANGLE FOR BEAM FORMING'
      WRITE(*,*) THETA_MAX

      READ(11,*) Vmin
      WRITE(*,*) 'THE SURFACE VELOCITY'
      WRITE(*,*) Vmin

      READ(11,*) Nbeam_h
      WRITE(*,*) 'THE NUMBER OF TRACES IN HALF BEAM'
      WRITE(*,*) Nbeam_h

      READ(11,*) NP_h
      WRITE(*,*) 'THE HALF NUMBER OF RAY PARAMETER'
      WRITE(*,*) NP_h

      READ(11,*) NS_START, NS_END 
			WRITE(*,*) 'THE START AND END SHOT NUMBER FOR PROCESS' 
      WRITE(*,*) NS_START, NS_END

      CLOSE(11)

      END

*===================================================================

      SUBROUTINE HAMMING_WINDOW(CTRACE, NW1, NW2, NW3, NW4, LT)

      PARAMETER(PAI=3.14159165359)
      INTEGER  Nw1, NW2, Nw3, Nw4
      complex  CTRACE(LT)

      DO IW=1, LT/2+1

         IF(IW.GE.NW1.AND.IW.LE.NW2) THEN
           Hammingw=0.5+0.5*cos(PAI*(Iw-Nw1)/(NW2-Nw1)-PAI)
           CTRACE(Iw)=CTRACE(Iw)*Hammingw
           IF(IW.GT.1.AND.IW.NE.LT+2-IW)THEN
             CTRACE(LT+2-IW)=CTRACE(LT+2-IW)*Hammingw
           ENDIF
         ELSE IF(IW.GE.NW3.AND.IW.LE.NW4) THEN
           Hammingw=0.5+0.5*cos(pai*(Nw3-Iw)/(Nw4-Nw3))
           CTRACE(Iw)=CTRACE(Iw)*Hammingw
           IF(IW.GT.1.AND.IW.NE.LT+2-IW)THEN
             CTRACE(LT+2-IW)=CTRACE(LT+2-IW)*Hammingw
           ENDIF
         ELSE IF(IW.GT.Nw4.OR.IW.LT.Nw1) THEN
          CTRACE(IW)=0.0
           IF(IW.GT.1.AND.IW.NE.LT+2-IW)THEN
             CTRACE(LT+2-IW)=0.0
           ENDIF
         END IF

       END DO

       RETURN
       END

!===================================================

      subroutine complex_CG(A,x1,x2,f,P1,P2,R1,R2,AP,np,err)

      integer np
      real    err, alfa1, alfa2, beta1, beta2
      complex r1norm1,r2norm1,p1ap1

      complex A(np,np),f(np)
      complex X1(np),X2(np)
      complex P1(np),P2(np)
      complex R1(np),R2(np)
      complex AP(np)

      norm2min=1.0E-3

      !zero buff
      do ip=1,np
        X1(ip)=cmplx(0.0,0.0)
        X2(ip)=cmplx(0.0,0.0)
        R1(ip)=cmplx(0.0,0.0)
        R2(ip)=cmplx(0.0,0.0)
        P1(ip)=cmplx(0.0,0.0)
        P2(ip)=cmplx(0.0,0.0)
      enddo

      !set initial value
      do ip=1,np
        X1(ip)=cmplx(0.0,0.0)
        R1(ip)=f(ip)
        P1(ip)=R1(ip)
      enddo

      !begin CG process

      iteration=1

      do 6666 while(1)

        !compute alfa1
        r1norm1=0.0
        do ip=1,np
          r1norm1=r1norm1+R1(ip)*CONJG(R1(ip))  !R1**2
        enddo
        r1norm=real(r1norm1)

        do ip=1,np
          AP(ip)=cmplx(0.0,0.0)
        enddo
        do ip=1,np
          do ip1=1,np
            AP(ip)=AP(ip)+A(ip,ip1)*P1(ip1)   !A*P1
          enddo
        enddo

        P1AP1=0.0
        do ip=1,np
           P1AP1=P1AP1+P1(ip)*CONJG(AP(ip))      !(P1,AP)
        enddo
        p1ap=real(p1ap1)

        if(abs(P1AP).lt.norm2min)then
          print*,'err:P1AP=',P1AP
          iteration=iteration+1
          goto 8888
        endif
        alfa1=r1norm/P1AP      !r1norm/P1AP

        !compute X2
        do ip=1,np
          X2(ip)=X1(ip)+alfa1*P1(ip)   !X2=X1+alfa1*P1
        enddo

        !compute R2
        do ip=1,np
          R2(ip)=R1(ip)-alfa1*AP(ip)   !R2=R1-alfa1*AP
        enddo

        !compute Beta2
        r2norm1=0.0
        do ip=1,np
          r2norm1=r2norm1+R2(ip)*CONJG(R2(ip))   !R2**2
        enddo
        r2norm=real(r2norm1)

        if(r1norm.lt.norm2min)then
          print*,'err:r1norm=',r1norm
          iteration=iteration+1
          goto 8888
        endif
        beta2=r2norm/r1norm

        !compute P2
        do ip=1,np
          P2(ip)=R2(ip)+beta2*P1(ip)
        enddo

        !update
        do ip=1,np
          X1(ip)=X2(ip)
          R1(ip)=R2(ip)
          P1(ip)=P2(ip)
        enddo

      iteration=iteration+1
      !print*,'r2norm,r1norm',r2norm,r1norm
      if(iteration.ge.1000)goto 8888
      !if(r2norm.gt.r1norm)goto 8888
      if(sqrt(r2norm).le.err)goto 8888

6666  enddo

8888  continue

!      print*,'sqrt(r2norm),err=',sqrt(r2norm),err
!      print*,'iteration=',iteration

      end

!===============================================================



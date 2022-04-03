       SUBROUTINE LS_TAU_P(shot_all1,ctp,tp,lt,
     + dx,dt,vmin,np_h,nx_h,np,nx,theta_max,threshold,f1,f2,f3,f4)

!      EXTERNAL HAMMING_WINDOW
      EXTERNAL COMPLEX_CG
      INCLUDE'fftw3.f'

      integer*8 plan_f,plan_b
      parameter(PAI2=2.0*3.1415926)
      REAL pmax
      REAL dp
      INTEGER np,nx,lt,np_h,nx_h

      REAL dx,dt,vmin
      REAL f1,f2,f3,f4
      REAL w1,w2,w3,w4
      REAL EMAX,LAMDA
      real threshold,theta_max
      complex sum1,sum2
      integer ix,it,i

!      real shot_all(lt,nx)
!      real semb(lt,np),tp(lt,np)
!      REAL energy(np)
!      complex ctrace(lt)
!      complex cshot(lt,nx)
!      complex csemb(lt,np),ctp(lt,np)
!      complex A(nx,np),AH(np,nx),AHA(np,np)
!      complex f(np),x1(np),x2(np)
!      complex P1(np),P2(np)
!      complex R1(np),R2(np)
!      complex AP(np)

      real shot_all1(lt*nx)
      complex ctp(lt,np)
      real tp(lt,np)
       
!      real,allocatable :: shot_all1(:)
!      complex,allocatable :: ctp(:,:)
!      real,allocatable :: tp(:,:)
      real,allocatable :: shot_all(:,:)
      real,allocatable :: semb(:,:)
      real,allocatable :: energy(:)
      complex,allocatable :: ctrace(:)
      complex,allocatable :: cshot(:,:)
      complex,allocatable :: csemb(:,:)
      complex,allocatable :: A(:,:)
      complex,allocatable :: AH(:,:)
      complex,allocatable :: AHA(:,:)
      complex,allocatable :: f(:)
      complex,allocatable :: x1(:)
      complex,allocatable :: x2(:)
      complex,allocatable :: P1(:)
      complex,allocatable :: P2(:)
      complex,allocatable :: R1(:)
      complex,allocatable :: R2(:)
      complex,allocatable :: AP(:)

!      ALLOCATE(shot_all1(lt*nx),STAT=IERR)
!      IF(IERR.NE.0.0)THEN
!      PRINT*,'ALLOCATE ERR'
!      STOP
!      ENDIF

!      ALLOCATE(ctp(lt,np),STAT=IERR)
!      IF(IERR.NE.0.0)THEN
!      PRINT*,'ALLOCATE ERR'
!      STOP
!      ENDIF

!      ALLOCATE(tp(lt,np),STAT=IERR)
!      IF(IERR.NE.0.0)THEN
!      PRINT*,'ALLOCATE ERR'
!      STOP
!      ENDIF

      ALLOCATE(shot_all(lt,nx),STAT=IERR)
      IF(IERR.NE.0.0)THEN
      PRINT*,'ALLOCATE ERR'
      STOP
      ENDIF

      ALLOCATE(semb(lt,np),STAT=IERR)
      IF(IERR.NE.0.0)THEN
      PRINT*,'ALLOCATE ERR'
      STOP
      ENDIF

      ALLOCATE(energy(np),STAT=IERR)
      IF(IERR.NE.0.0)THEN
      PRINT*,'ALLOCATE ERR'
      STOP
      ENDIF

      ALLOCATE(ctrace(lt),STAT=IERR)
      IF(IERR.NE.0.0)THEN
      PRINT*,'ALLOCATE ERR'
      STOP
      ENDIF

      ALLOCATE(cshot(lt,nx),STAT=IERR)
      IF(IERR.NE.0.0)THEN
      PRINT*,'ALLOCATE ERR'
      STOP
      ENDIF

      ALLOCATE(csemb(lt,np),STAT=IERR)
      IF(IERR.NE.0.0)THEN
      PRINT*,'ALLOCATE ERR'
      STOP
      ENDIF
 
      ALLOCATE(A(nx,np),STAT=IERR)
      IF(IERR.NE.0.0)THEN
      PRINT*,'ALLOCATE ERR'
      STOP
      ENDIF

      ALLOCATE(AH(np,nx),STAT=IERR)
      IF(IERR.NE.0.0)THEN
      PRINT*,'ALLOCATE ERR'
      STOP
      ENDIF

      ALLOCATE(AHA(np,np),STAT=IERR)
      IF(IERR.NE.0.0)THEN
      PRINT*,'ALLOCATE ERR'
      STOP
      ENDIF

      ALLOCATE(f(np),STAT=IERR)
      IF(IERR.NE.0.0)THEN
      PRINT*,'ALLOCATE ERR'
      STOP
      ENDIF

      ALLOCATE(x1(np),STAT=IERR)
      IF(IERR.NE.0.0)THEN
      PRINT*,'ALLOCATE ERR'
      STOP
      ENDIF

      ALLOCATE(x2(np),STAT=IERR)
      IF(IERR.NE.0.0)THEN
      PRINT*,'ALLOCATE ERR'
      STOP
      ENDIF

      ALLOCATE(P1(np),STAT=IERR)
      IF(IERR.NE.0.0)THEN
      PRINT*,'ALLOCATE ERR'
      STOP
      ENDIF

      ALLOCATE(P2(np),STAT=IERR)
      IF(IERR.NE.0.0)THEN
      PRINT*,'ALLOCATE ERR'
      STOP
      ENDIF

      ALLOCATE(R1(np),STAT=IERR)
      IF(IERR.NE.0.0)THEN
      PRINT*,'ALLOCATE ERR'
      STOP
      ENDIF

      ALLOCATE(R2(np),STAT=IERR)
      IF(IERR.NE.0.0)THEN
      PRINT*,'ALLOCATE ERR'
      STOP
      ENDIF

      ALLOCATE(AP(np),STAT=IERR)
      IF(IERR.NE.0.0)THEN
      PRINT*,'ALLOCATE ERR'
      STOP
      ENDIF

      Pmax=sin(theta_max*PAI2/360.0)/Vmin
      Dp=2*Pmax/(np-1)
      dw=PAI2*(1.0/(lt*dt))
!      print*,'theta_max,pmax,dp,dw=',theta_max,pmax,dp,dw

      i=0
      do it=1,lt
        do ix=1,nx
          i=i+1
          shot_all(it,ix)=shot_all1(i)
        enddo
      enddo

!      f1=2.0
!      f2=10.0
!      f3=60.0
!      f4=66.0

      call sfftw_plan_dft_1d(plan_f,lt,ctrace,ctrace,
     +  FFTW_BACKWARD,FFTW_MEASURE)
      call sfftw_plan_dft_1d(plan_b,lt,ctrace,ctrace,
     +  FFTW_FORWARD,FFTW_MEASURE)

      Dw=PAI2/(lt*dt)
      W1=PAI2*F1
      Nw1=W1/Dw+0.5
      W2=PAI2*F2
      Nw2=W2/Dw+0.5
      W3=PAI2*F3
      Nw3=W3/Dw+0.5
      W4=PAI2*F4
      Nw4=W4/Dw+0.5
      NW=NW4-NW1+1

!      print*,'Pmax,Dp=',Pmax,Dp
!      print*,'dx,dt=',dx,dt
!      print*,'nw1,nw4,dw=',nw1,nw4,dw

      !transform time domain to frequency domain
      print*,'fft'


      do 1111 ix=1,nx

        do it=1,lt
          ctrace(it)=cmplx(shot_all(it,ix),0.0)
        enddo


        CALL sfftw_execute(plan_b)
        DO IT=1,LT
         CTRACE(IT)=CTRACE(IT)/real(LT)
        ENDDO

!        call HAMMING_WINDOW(CTRACE, NW1, NW2, NW3, NW4, LT)

        do it=1,lt
          cshot(it,ix)=ctrace(it)
        enddo

1111  continue
**************  get the max energy for constrain*********

      do 7777 iw=1,lt/2+1

         w=(iw-1)*dw

         do ip=1,np
           p=(ip-np_h)*dp
           do ix=1,nx

            x=(ix-nx_h)*dx

             phase=w*p*x
             COE=1.0/SQRT(REAL(Nx))
             A(IX,IP)=COE*CMPLX(COS(PHASE),SIN(PHASE))
           enddo
         enddo

         do ip=1,np
           do ix=1,nx
             AH(ip,ix)=CONJG(A(ix,ip))
           enddo
         enddo


         DO ip=1,np
           sum1=0.0
            DO ix=1,nx
              sum1=sum1+AH(ip,ix)*cshot(iw,ix)
            END DO
           f(ip)=sum1
           ENERGY(ip)=ENERGY(ip)+REAL(F(ip)*CONJG(F(ip)))
         END DO
7777   continue

       EMAX=0.0
       DO ip=1,np
        IF(ENERGY(IP).GT.EMAX)EMAX=ENERGY(IP)
          ENERGY(IP)=0.0
       END DO
       LAMDA=0.1*EMAX
*****************   rlsi tau-p transform   *******************
       DO 2222 IW=1,lt/2+1
         W=(IW-1)*DW

         DO IP=1,NP
           P=(ip-np_h)*dp
           DO IX=1,nx
             X=(IX-nx_h)*DX
             PHASE=W*P*X
             COE=1.0/SQRT(REAL(Nx))
             A(IX,IP)=COE*CMPLX(COS(PHASE),SIN(PHASE))
           END DO
         END DO

         DO IP=1,NP
           DO IX=1,Nx
             AH(IP,IX)=CONJG(A(IX,IP))
           END DO
         END DO

         DO IP=1,NP
           DO IP1=1, NP
             SUM1= 0.0
             DO IX=1, nx
               SUM1= SUM1+AH(IP,IX)*A(IX,IP1)
             END DO
             AHA(IP,IP1)=SUM1
           END DO
         END DO

         DO IP=1,NP
           SUM1=0.0
            DO IX=1,nx
              SUM1=SUM1+AH(IP,IX)*cshot(IW,IX)
            END DO
           F(IP)= SUM1
           ENERGY(IP)=ENERGY(IP)+REAL(F(IP)*CONJG(F(IP)))
         END DO

         DO IP=1,NP
           IF(ENERGY(IP).LT.1.0E-6)ENERGY(IP)=1.0E-6
           AHA(IP,IP)=AHA(IP,IP)+LAMDA/ENERGY(IP)
           AHA(ip,ip)=AHA(ip,ip)+0.1
           AHA(IP,IP)=CMPLX(REAL(AHA(IP,IP)),0.0)
         END DO

         ERR= 1.0E-15
         CALL COMPLEX_CG(AHA,X1,X2,F,P1,P2,R1,R2,AP,NP,ERR)

         DO IP=1,np
            ctp(IW,IP)= X2(IP)
           IF(IW.Gt.1.AND.IW.NE.LT+2-IW)
     +        ctp(LT+2-IW,ip)=CONJG(ctp(iw,ip))
         END DO
2222  CONTINUE

!======transform frequency domain to time domain
      print*,'ifft'
      do 3333 ip=1,np

        do it=1,lt
          ctrace(it)=ctp(it,ip)
        enddo

        call sfftw_execute(plan_f)

        do it=1,lt
          tp(it,ip)=real(ctrace(it))
        enddo

3333  continue


      call sfftw_destroy_plan(plan_f)
      call sfftw_destroy_plan(plan_b)

!      deallocate (shot_all1)
!      deallocate (ctp)
!      deallocate (tp)
      deallocate (shot_all)
      deallocate (semb)
      deallocate (energy)
      deallocate (ctrace)
      deallocate (cshot)
      deallocate (csemb)
      deallocate (A)
      deallocate (AH)
      deallocate (AHA)
      deallocate (f)
      deallocate (x1)
      deallocate (x2)
      deallocate (P1)
      deallocate (P2)
      deallocate (R1)
      deallocate (R2)
      deallocate (AP)  

      END SUBROUTINE

!=============complex_CG=============================

      subroutine complex_CG(A,x1,x2,f,P1,P2,R1,R2,AP,np,err)

      integer np
      real    err, alfa1, alfa2, beta1, beta2
      complex r1norm1,r2norm1,p1ap1

      complex A(np,np),f(np)
      complex X1(np),X2(np)
      complex P1(np),P2(np)
      complex R1(np),R2(np)
      complex AP(np)

!      complex,allocatable :: A(:,:)
!      complex,allocatable :: f(:)
!      complex,allocatable :: X1(:)
!      complex,allocatable :: X2(:)
!      complex,allocatable :: P1(:)
!      complex,allocatable :: P2(:)
!      complex,allocatable :: R1(:)
!      complex,allocatable :: R2(:)
!      complex,allocatable :: AP(:)

!      ALLOCATE(A(np,np),STAT=IERR)
!      IF(IERR.NE.0.0)THEN
!      PRINT*,'ALLOCATE ERR'
!      STOP
!      ENDIF     
      
!      ALLOCATE(f(np),STAT=IERR)
!      IF(IERR.NE.0.0)THEN
!      PRINT*,'ALLOCATE ERR'
!      STOP
!      ENDIF     

!      ALLOCATE(X1(np),STAT=IERR)
!      IF(IERR.NE.0.0)THEN
!      PRINT*,'ALLOCATE ERR'
!      STOP
!      ENDIF     

!      ALLOCATE(X2(np),STAT=IERR)
!      IF(IERR.NE.0.0)THEN
!      PRINT*,'ALLOCATE ERR'
!      STOP
!      ENDIF     

!      ALLOCATE(P1(np),STAT=IERR)
!      IF(IERR.NE.0.0)THEN
!      PRINT*,'ALLOCATE ERR'
!      STOP
!      ENDIF     

!      ALLOCATE(P2(np),STAT=IERR)
!      IF(IERR.NE.0.0)THEN
!      PRINT*,'ALLOCATE ERR'
!      STOP
!      ENDIF     

!      ALLOCATE(R1(np),STAT=IERR)
!      IF(IERR.NE.0.0)THEN
!      PRINT*,'ALLOCATE ERR'
!      STOP
!      ENDIF     

!     ALLOCATE(R2(np),STAT=IERR)
!      IF(IERR.NE.0.0)THEN
!      PRINT*,'ALLOCATE ERR'
!      STOP
!      ENDIF     

!      ALLOCATE(AP(np),STAT=IERR)
!      IF(IERR.NE.0.0)THEN
!      PRINT*,'ALLOCATE ERR'
!      STOP
!      ENDIF     

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

!      do 6666 while(1)
      do iteration = 1, 1000

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
!          print*,'err:P1AP=',P1AP
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
!          print*,'err:r1norm=',r1norm
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

!      iteration=iteration+1
      !print*,'r2norm,r1norm',r2norm,r1norm
      if(iteration.ge.1000)goto 8888
      !if(r2norm.gt.r1norm)goto 8888
      if(sqrt(r2norm).le.err)goto 8888

6666  enddo

8888  continue

!      print*,'sqrt(r2norm),err=',sqrt(r2norm),err
!      print*,'iteration=',iteration

!      deallocate (A)
!      deallocate (f)
!      deallocate (X1)
!      deallocate (X2)
!      deallocate (P1)
!      deallocate (P2)
!      deallocate (R1)
!      deallocate (R2)
!      deallocate (AP)


      END SUBROUTINE


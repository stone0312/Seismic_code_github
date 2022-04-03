      PROGRAM Aij_estimation
      IMPLICIT NONE
      REAL,PARAMETER::pi=3.1415926
      INTEGER lt
      INTEGER nx
      INTEGER np
      INTEGER nlocal
      INTEGER nx_expand
      INTEGER len_local
      INTEGER nx_local_h,np_h
      INTEGER dit,dip
      REAL    delt_tp_value
      REAL    vmin
      REAL    z_line
      REAL    dx
      REAL    dt
      REAL    theta_max
      
      CHARACTER*256 FN1
      
      REAL    pmax
      REAL    dp
      
      CALL READPAR(FN1,lt,nx,np,nlocal,nx_expand,
     +  len_local,np_h,nx_local_h,dit,dip,delt_tp_value,
     +  vmin,z_line,dx,dt,theta_max)
      
      pmax=sin(2*pi*theta_max/360.0)/vmin
      dp=2*pmax/(np-1)
      
      CALL COMPUTER_Aij_ALL(FN1,lt,nx,np,nlocal,nx_expand,
     +  len_local,np_h,nx_local_h,dit,dip,delt_tp_value,    
     +  vmin,z_line,dx,dt,dp)


      END PROGRAM
!****************************************************
!
!     cumputer the Aij_all
!
!****************************************************
      SUBROUTINE COMPUTER_Aij_ALL(FN1,lt,nx,np,nlocal,nx_expand,
     +  len_local,np_h,nx_local_h,dit,dip,delt_tp_value,
     +  vmin,z_line,dx,dt,dp)
      IMPLICIT NONE
      INTEGER lt 
      INTEGER nx
      INTEGER np
      INTEGER nlocal         !  the number of local_tau_p
      INTEGER nx_expand      !  the length of expanded pg
      INTEGER len_local      !  the length of local
      INTEGER np_h
      INTEGER nx_local_h,dit,dip
      REAL    dx
      REAL    dt
      REAL    dp
      REAL    vmin,z_line
      REAL delt_tp_value

      CHARACTER*256 FN1

      REAL pg(lt,nx)
      REAL pg_local(lt,len_local)   
      REAL pg_expand(lt,nx_expand)
      REAL tp(lt,np)
      COMPLEX ctp_local(lt,np)
      REAL tp_sum(6000,40000)
      REAL aij_medium(nx)
      REAL aij_all(nx,nx)     ! which we want

      INTEGER ilocal,it,ix,itrace,ip
      INTEGER len_local_h
!
      INTEGER ghost_number(nx)! the number of ghost for each receiver
!                               gather  2013.10.28.pm
      INTEGER::ghost_number0=0
	  
	  
      OPEN(55,FILE=FN1,ACTION='READ',FORM='BINARY')
      READ(55)pg
      CLOSE(55)
      
      len_local_h=NINT((len_local-1)/2.0)
      pg_expand(:lt,(len_local_h+1):(nx+len_local_h))=pg(:lt,1:nx)
      
      DO ilocal=1,nlocal
                
        pg_local(:lt,1:len_local)=
     +                         pg_expand(:lt,ilocal:ilocal+len_local-1)
      
	 
        CALL LS_TAU_P(pg_local,ctp_local,tp
     +       ,lt,len_local,np,dx,dt,vmin,np_h,nx_local_h)
              
!*************************************************88       
        IF(ilocal>=1.AND.ilocal<=len_local_h)THEN
      
        tp=tp*len_local/(len_local_h+ilocal)

        tp_sum(:lt,(ilocal-1)*np+1:(ilocal-1)*np+np)
     +   =tp(:lt,:np)
       print*,'text_local_tp'
       OPEN(61,FILE='text_local_tp.dat',
     +   FORM='BINARY',ACCESS='SEQUENTIAL')
       WRITE(61)tp
       CLOSE(61)  
       pause               
        
        CALL COMPUTER_LEAD_Aij(tp,lt,nx,np,delt_tp_value,
     + dit,dip,dp,dx,z_line,vmin,np_h,
     + ilocal,aij_medium,ghost_number0)
        aij_all(ilocal,:nx)=aij_medium
        ghost_number(ilocal)=ghost_number0
!        print*,'ilocal=',ilocal
!        print*,aij_medium
                
        END IF
        
!***************************************************
        IF(ilocal>=(len_local_h+1).AND.ilocal<=(np-len_local_h))THEN

        tp_sum(:lt,(ilocal-1)*np+1:(ilocal-1)*np+np)
     +   =tp(:lt,:np)

        CALL COMPUTER_LEAD_Aij(tp,lt,nx,np,delt_tp_value,
     +       dit,dip,dp,dx,z_line,vmin,np_h,ilocal,aij_medium,
     +       ghost_number0)
        aij_all(ilocal,:nx)=aij_medium
        ghost_number(ilocal)=ghost_number0
!        print*,'ilocal=',ilocal
!        print*,aij_medium
!        pause            

        END IF
        
!****************************************************8
        IF(ilocal>=(np-len_local_h+1).AND.ilocal<=(np))THEN

        tp=tp*len_local/(np+len_local_h+1-ilocal)

        tp_sum(:lt,(ilocal-1)*np+1:(ilocal-1)*np+np)
     +   =tp(:lt,:np)

        CALL COMPUTER_LEAD_Aij(tp,lt,nx,np,delt_tp_value,
     +       dit,dip,dp,dx,z_line,vmin,np_h,ilocal,aij_medium,
     +       ghost_number0)
        aij_all(ilocal,:nx)=aij_medium
        ghost_number(ilocal)=ghost_number0        
        END IF
      
      END DO
      
      print*,'output the aij_all:::'
      OPEN(61,FILE='aij_matrix_estimation_pg1.dat',
     +     FORM='BINARY',ACCESS='SEQUENTIAL')
      WRITE(61)aij_all
      CLOSE(61)

      OPEN(61,FILE='ghost_number_pg1.dat',
     +     FORM='BINARY',ACCESS='SEQUENTIAL')
      WRITE(61)ghost_number
      CLOSE(61)

      
!      DO ilocal=1,nlocal
!        PRINT*,'the reveiver =',ilocal
!        print*,aij_all(ilocal,:nx)
!        pause        
!      END DO
      
      END SUBROUTINE
!*******************************************************
!
!       parameter card
!
!**********************************************************
      SUBROUTINE READPAR(FN1,lt,nx,np,nlocal,nx_expand,
     +  len_local,np_h,nx_local_h,dit,dip,delt_tp_value,
     +  vmin,z_line,dx,dt,theta_max)     
      INTEGER lt
      INTEGER nx
      INTEGER np
      INTEGER nlocal
      INTEGER nx_expand
      INTEGER len_local
      INTEGER nx_local_h,np_h
      INTEGER dit,dip
      REAL    delt_tp_value
      REAL    vmin
      REAL    z_line
      REAL    dx
      REAL    dt
      REAL    theta_max
      CHARACTER*256 FN1

      
      OPEN(11,FILE='aij_estimation.par')
    
      READ(11,'(A)')FN1
      
      READ(11,*)lt

      READ(11,*)nx

      READ(11,*)np

      READ(11,*)nlocal

      READ(11,*)nx_expand

      READ(11,*)len_local

      READ(11,*)np_h

      READ(11,*)nx_local_h

      READ(11,*)dit

      READ(11,*)dip

      READ(11,*)delt_tp_value

      READ(11,*)vmin

      READ(11,*)z_line

      READ(11,*)dx

      READ(11,*)dt

      READ(11,*)theta_max

      CLOSE(11)

      END SUBROUTINE
!***********************************************************
!
!      cumputer the final aij for one receiver
!
!***********************************************************
      SUBROUTINE COMPUTER_LEAD_Aij(tp_value,lt,nx,np,delt_tp_value,
     +           dit,dip,dp,dx,z_line,vmin,np_h,target_trace,Aij,
     +           ghost_number0)
      IMPLICIT NONE
      INTEGER lt    !  the number of time_sampling
      INTEGER nx
      INTEGER GHOST_NUMBER0
      INTEGER target_trace
      INTEGER np    !  the number of ray-parameter p_sampling
      INTEGER np_h
      INTEGER dit   !
      INTEGER dip   !   determine the size of tp_block
      REAL    delt_tp_value
      REAL    dp
      REAL    dx
      REAL    z_line
      REAL    vmin
      REAL    Aij(nx)
      REAL    tp_value(lt,np)

      INTEGER it,ip
      INTEGER count_minimum1
      REAL    minimum1
      REAL    min_value(3)  !  minimum tp_value,it,ip

      aij=0.0

!   search  the minimum 1         
      minimum1=9999.0
      DO it=1,lt
        DO ip=1,np

          IF(minimum1>tp_value(it,ip))THEN
             minimum1=tp_value(it,ip)
             min_value(2)=it
             min_value(3)=ip
          END IF

         END DO
       END DO
       
       min_value(1)=minimum1
       
!    search the min_range 1      
       
      count_minimum1=0
      DO it=1,lt
        DO ip=1,np
         IF(tp_value(it,ip)>=min_value(1).AND.
     +     tp_value(it,ip)<=(min_value(1)+delt_tp_value))THEN
          count_minimum1=count_minimum1+1
         END IF
        END DO
      END DO
      
      CALL CIRCLE1(tp_value,lt,nx,np,dp,minimum1,delt_tp_value,
     +           dx,np_h,count_minimum1,min_value,dit,dip,
     =           target_trace,vmin,z_line,Aij,ghost_number0)

      
      END SUBROUTINE
!******************************************************************
!
!    cumputer the 1th ghost_tp,prapare for the 2nd,3rd,.....
!
!******************************************************************
      SUBROUTINE CIRCLE1(tp_value,lt,nx,np,dp,minimum1,delt_tp_value,
     +           dx,np_h,count_minimum1,min_value,dit,dip,
     =           target_trace,vmin,z_line,Aij,ghost_number0)
      IMPLICIT NONE
      INTEGER lt,ghost_number0
      INTEGER nx
      INTEGER np
      INTEGER np_h
      INTEGER dit
      INTEGER dip
      INTEGER target_trace
      REAL minimum1
      REAL minimum2
      REAL delt_tp_value
      REAL dx
      REAL dp
      REAL vmin
      REAL z_line
      INTEGER count_minimum1
      REAL Aij(nx)
      REAL min_value(3)
      REAL tp_value(lt,np)
      REAL min_range_value1(3,count_minimum1)
      REAL ghost_tp(3,lt)  !minimum tp-value
      REAL ghost_tp_max(3,lt)  ! maximum tp_value
      REAL ghost_tp_real(3,lt) ! the real tp_value
      REAL min_range_value1_1(3,count_minimum1)
      INTEGER it_medium
      REAL maximum
      INTEGER it,ip,icount
      INTEGER count_left,count_ghost
     
      ghost_tp_max=0.0
      ghost_tp_real=0.0
      ghost_tp=0.0
      ghost_tp(:3,1)=min_value
        
!      print*,'count_minimum1=',target_trace,count_minimum1
!      print*,'ghost_tp',ghost_tp(:3,1)
!      pause
      
      icount=0

      DO it=1,lt
        DO ip=1,np
         IF(tp_value(it,ip)>=minimum1.AND.
     +     tp_value(it,ip)<=(minimum1+delt_tp_value))THEN
           icount=icount+1
           min_range_value1(1,icount)=tp_value(it,ip)
           min_range_value1(2,icount)=it
           min_range_value1(3,icount)=ip
         END IF
        END DO
      END DO
           
!      print*,'icount=count_minimum1',icount,count_minimum1
      
!      DO it=1,count_minimum1
!        print*,it
!        print*,min_range_value1(:3,it)  
!      END DO
!      pause

      count_left=0
      min_range_value1_1=0.0
      DO it=1,count_minimum1
        IF(min_range_value1(2,it)>(min_value(2)+dit).OR.
     +     min_range_value1(2,it)<(min_value(2)-dit))THEN
          count_left=count_left+1
          min_range_value1_1(1,count_left)=min_range_value1(1,it)
          min_range_value1_1(2,count_left)=min_range_value1(2,it)
          min_range_value1_1(3,count_left)=min_range_value1(3,it)
        END IF
      END DO
           
!      print*,'count_left1',count_left
!      pause
      count_ghost=1
!      DO it=1,count_minimum1
!        print*,it
!        print*,min_range_value1_1(:3,it)
!      end do
!      pause
 
      DO WHILE(min_range_value1_1(1,1)/=0.0)
        
        count_ghost=count_ghost+1
        CALL CIRCLE2(min_range_value1_1,count_minimum1,count_left,
     +   delt_tp_value,dit,dip,ghost_tp(:3,count_ghost),ghost_number0)  
!        print*,'count_ghost',count_ghost
!        pause
      END DO

!       print*,'output the count_ghost'
!       print*,'count_ghost',count_ghost
!       print*,'ghost_information',ghost_tp(:3,:6)
!       pause           
!****************************************************************      
!        verify the real ghost_tp_real matrix
!
!***************************************************************
      it_medium=NINT((2*z_line/vmin)/0.002)

       DO icount=1,count_ghost
       maximum=-99999.99

          DO it=ghost_tp(2,icount)-it_medium,
     *          ghost_tp(2,icount)+it_medium
            DO ip=1,np
              IF(maximum<=tp_value(it,ip))THEN
                maximum=tp_value(it,ip)
                ghost_tp_max(2,icount)=it
                ghost_tp_max(3,icount)=ip
              END IF
            END DO
          END DO


       END DO

       DO icount=1,count_ghost
        IF(ghost_tp(2,icount)<ghost_tp_max(2,icount))THEN
          ghost_tp_real(1,icount)=ghost_tp(1,icount)
          ghost_tp_real(2,icount)=ghost_tp(2,icount)
          ghost_tp_real(3,icount)=ghost_tp(3,icount)
        ELSE
          ghost_tp_real(1,icount)=ghost_tp_max(1,icount)
          ghost_tp_real(2,icount)=ghost_tp_max(2,icount)
          ghost_tp_real(3,icount)=ghost_tp_max(3,icount)
        END IF
       END DO
!*************************************************************************

      CALL CUMPUTER_Aij(ghost_tp,nx,dp,np_h,z_line,vmin,Aij
     + ,dx,target_trace,lt,ghost_number0)
            
      END SUBROUTINE
!********************************************************************
!
!            involved in the circle 
!     cumputer the 3rd,4th,5th.....ghost for the target_trace
!
!********************************************************************
      SUBROUTINE CIRCLE2(min_range_value,count_minimum,count_left
     +   ,delt_tp_value,dit,dip,ghost_tp1)
      INTEGER count_minimum                    !the 1th range ensure
!                     size of matrix to avoid allocate
      INTEGER count_left,dit,dip               ! the length of this
!                                                range
      REAL    delt_tp_value
      REAL    min_range_value(3,count_minimum) ! the left range after
!                                                last election
      REAL    min_range_value1_1(3,count_minimum) 
      REAL    ghost_tp1(3)  ! the mini-value(3) of the min_range_value(3,)
      REAL    minimum
  
      INTEGER it,ip,icount    
      
      min_range_value1_1=0.0
      ghost_tp1=0.0
      minimum=9999.0
      
      DO it=1,count_minimum
        IF(minimum>min_range_value(1,it))THEN
          minimum=min_range_value(1,it)
          ghost_tp1(2)=min_range_value(2,it)
          ghost_tp1(3)=min_range_value(3,it)
        END IF
      END DO

      ghost_tp1(1)=minimum

!      print*,'ghost_tp1'
!      print*,ghost_tp1,count_left
!      pause

      count_left=0
      DO it=1,count_minimum
      IF(min_range_value(2,it)>(ghost_tp1(2)+dit).OR.
     +   min_range_value(2,it)<(ghost_tp1(2)-dit))THEN
          count_left=count_left+1
          min_range_value1_1(1,count_left)=min_range_value(1,it)
          min_range_value1_1(2,count_left)=min_range_value(2,it)
          min_range_value1_1(3,count_left)=min_range_value(3,it)
       END IF
       END DO

!       print*,'count_left',count_left
!       pause
!       print*,'min_range_value1_1'
!       do it=1,count_minimum
!         print*,it
!         print*,min_range_value1_1(:3,it)
!       end do
!       pause
       min_range_value=0.0
       min_range_value=min_range_value1_1
      

      END SUBROUTINE
!********************************************************************
!
!                   cumputer the Aij
!
!**********************************************************************
      SUBROUTINE CUMPUTER_Aij(min_tp_value,nx,dp,np_h,z_line,vmin,Aij
     + ,dx,target_trace,lt,ghost_number0)
      IMPLICIT NONE
      INTEGER lt
      INTEGER nx      !  the number of recevers
      REAL    dx      !  the sampling interval of recevers
      REAL    dp      !  the sampling interval of ray_parameter p
      INTEGER    np_h    !  the central np of local tau-p transform
      INTEGER target_trace !  the trace which we analyze
      REAL    z_line  !  the depth of streamer
      REAL    vmin    !  the water of velocity
      REAL    Aij(nx)   ! the parameter matrix which we want
      REAL min_tp_value(3,lt)  !  the searched ghost tp_value
      
      INTEGER  itrace,ix
      INTEGER  count_ghost
      REAL     minimum
      REAL     time_del(nx)         !  the optimal time_delay from p
      REAL     time_del_all(nx)  ! the delt_tij  matrix
      REAL     time_err(nx)      ! the difference

      INTEGER count_media,ghost_number0! computer the number of ghost
!     +                                   for each receiver gather 
      REAL    sum_errinverse
      REAL    errinverse_arry(nx)

      
      count_ghost=0 
!*******  cumputer the exect time_del(count_ghost)  ***************
      DO itrace=1,lt
         IF(min_tp_value(1,itrace)/=0.0)THEN  
         count_ghost=count_ghost+1
         time_del(count_ghost)=2*z_line/     
     +  (vmin*SQRT(1.0-vmin**2*
     +  ((-min_tp_value(3,itrace)+np_h)*dp)**2))

!         time_del(count_ghost)=(2*z_line*SQRT(1.0-vmin**2*
!     +   ((-min_tp_value(3,itrace)+np_h)*dp)**2))/vmin
       
         END IF
      END DO

      ghost_number0=count_ghost
!*****************************************************************8
!
!   cumpute the  time_del-all(nx) ,shot one side
! 
!********************************************************************8
        time_del_all=0.0
        DO ix=1,1
         IF(target_trace>=ix)THEN
         time_del_all(ix)=
     + SQRT((2*z_line)**2+((target_trace-ix)*dx)**2)/vmin
         END IF
        END DO
       
        DO ix=2,nx
         IF(target_trace>=ix)THEN      !  geng gai le 13.11.11
         time_del_all(ix)=
     + SQRT((2*z_line)**2+((target_trace-ix)*dx)**2)/vmin
         END IF
        END DO
!******************************************************************
!
!    compute the exact time_del with time_del_all to ensure the 2th
!    source
!
!******************************************************************              
      aij=0
      DO itrace=1,count_ghost
        minimum=999999.0

        DO ix=1,nx
         time_err(ix)=
     +      ABS(time_del_all(ix)-time_del(itrace))
        END DO

        DO ix=1,nx
          IF(minimum>time_err(ix))THEN
            minimum=time_err(ix)
          END IF
        END DO

        DO ix=1,nx
          IF(time_err(ix)==minimum)THEN
            Aij(ix)=aij(ix)+1.0
          END IF
        END DO

       END DO

       END SUBROUTINE
!**********************************************************************
!
!         complex ls tau-p transform
!
!*********************************************************************
      SUBROUTINE LS_TAU_P(shot_all,ctp,tp,lt,nx,np,
     + dx,dt,vmin,np_h,nx_h)
      EXTERNAL HAMMING_WINDOW
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
      real shot_all(lt,nx)
      real semb(lt,np),tp(lt,np)
      REAL energy(np)
      complex ctrace(lt)
      complex cshot(lt,nx)
      complex csemb(lt,np),ctp(lt,np)
      complex A(nx,np),AH(np,nx),AHA(np,np)
      complex f(np),x1(np),x2(np)
      complex P1(np),P2(np)
      complex R1(np),R2(np)
      complex AP(np)
      complex sum1,sum2

      REAL::threshold=0.5
      REAL::theta_max=80.0

      Pmax=sin(theta_max*PAI2/360.0)/Vmin
      Dp=2*Pmax/(np-1)
      dw=PAI2*(1.0/lt*dt)
      print*,'theta_max,pmax,dp,dw=',theta_max,pmax,dp,dw

  
      f1=2.0
      f2=10.0
      f3=60.0
      f4=66.0

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

      print*,'Pmax,Dp=',Pmax,Dp
      print*,'dx,dt=',dx,dt
      print*,'nw1,nw4,dw=',nw1,nw4,dw

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

      END SUBROUTINE
!****************************************************************8
!
!     complex backwrd tau-p transform
!
!*****************************************************************
      SUBROUTINE BACK_CTAU_P(ctp,local,lt,dt,dx,dp
     +     ,dw,nx_local_h,np_h,nx_local,np)
      INCLUDE'fftw3.f'
      INTEGER*8 plan_f,plan_b
      REAL,PARAMETER::pai2=2.0*3.1415926
      INTEGER lt
      INTEGER nx_local,np,nx_local_h,np_h
      REAL dt,dx,dp,dw
      REAL p,x,w
      REAL phase
      REAL local(lt,nx_local)
      COMPLEX ctp(lt,np)
      COMPLEX c_local(lt,nx_local)
      COMPLEX A(nx_local,np)
      COMPLEX ctrace(lt)
      COMPLEX sum1,sum2

      INTEGER iw,ix,ip,it

      call sfftw_plan_dft_1d(plan_f,lt,ctrace,ctrace,
     +  FFTW_BACKWARD,FFTW_MEASURE)
      call sfftw_plan_dft_1d(plan_b,lt,ctrace,ctrace,
     +  FFTW_FORWARD,FFTW_MEASURE)

      DO iw=1,lt/2+1
        w=(iw-1)*dw

        DO ip=1,np
          p=(ip-np_h)*dp
          DO ix=1,nx_local
            x=(ix-nx_local_h)*dx
            phase=w*p*x
            coe=1.0/sqrt(real(nx_local))
            A(ix,ip)=coe*CMPLX(cos(phase),sin(phase))
          END DO
         END DO

         DO ix=1,nx_local
           sum1=0.0
           DO ip=1,np
            sum1=sum1+A(ix,ip)*ctp(iw,ip)
           END DO
           c_local(iw,ix)=sum1
           IF(iw.gt.1.and.iwlt.ne.lt+2-iw)THEN
            c_local(lt+2-iw,ix)=CONJG(c_local(iw,ix))
           END IF
         END DO
      END DO

      DO ix=1,nx_local
        DO it=1,lt
          ctrace(it)=c_local(it,ix)
        END DO
        CALL  sfftw_execute(plan_f)
        DO it=1,lt
          local(it,ix)=REAL(CTRACE(it))
        END DO

      END DO

      call sfftw_destroy_plan(plan_f)
      call sfftw_destroy_plan(plan_b)

      END SUBROUTINE
!***********************************************************************
!
!    complex  CG
!
!**********************************************************************
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

      print*,'sqrt(r2norm),err=',sqrt(r2norm),err
      print*,'iteration=',iteration

      END SUBROUTINE



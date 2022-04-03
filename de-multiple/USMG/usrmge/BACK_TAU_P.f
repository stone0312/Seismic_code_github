      SUBROUTINE BACK_TRANSFORM(CO_TP_LOCAL, CO_LOCAL,
     +       CCO_TP_LOCAL, CCO_LOCAL, A, 
     +       LT, DT, DX,
     +       NW1, NW2, NW3, NW4, NW,
     +       DP, DW, Nbeam_h, NP_h, Nbeam, NP)
!     +       plan_f,plan_b)

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

2222   continue
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



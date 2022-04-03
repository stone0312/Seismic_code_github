       SUBROUTINE LS_TAU_P(shot_all,ctp,tp,lt,
     + dx,dt,vmin,np_h,nx_h,np,nx,theta_max,threshold,f1,f2,f3,f4)

      INCLUDE'fftw3.f'

      integer*8 plan_f,plan_b
      parameter(PAI2=2.0*3.1415926)
      REAL pmax,dw
      REAL dp
      INTEGER np,nx,lt,np_h,nx_h

      REAL dx,dt,vmin
      REAL f1,f2,f3,f4
      REAL w1,w2,w3,w4
      REAL EMAX,LAMDA
      real threshold,theta_max
      real shot_all(lt,nx)
      real tp(lt,np)
      complex ctp(lt,np)


      print*,'befor fftw'


      END SUBROUTINE


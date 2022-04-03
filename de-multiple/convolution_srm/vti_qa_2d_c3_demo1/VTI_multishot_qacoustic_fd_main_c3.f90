
        !Modified to conduct Apture_control during forward modelling
        !Modified to add su head to shot_gather
        !Modified for nanjing
        include 'VTI_multishot_qacoustic_fd_subroutine_c3.f90'

        program VTI_fd_qacoustic_modelling
        implicit none 
        character(len=256)::par_fn=&
        './vti_qa_2d_modelling_new.par'
        call multi_shot_modelling(par_fn)
        end program VTI_fd_qacoustic_modelling

        subroutine multi_shot_modelling(par_fn)
        use mpi
        use constant
        implicit none
        !Dummy variables
        !*Parameter_card_name*
        character(len=256)::par_fn
        !*Read from parameter_card*
        character(len=256)::vz_fn
        character(len=256)::epsilon_fn
        character(len=256)::delta_fn
        character(len=256)::shot_fn1
        character(len=256)::shot_fn2
        integer::nx_v,nz_v
        integer::nt
        integer::fd_order_explicit
        integer::order_pml_1st
        integer::order_pml_2nd
        integer::nshot
        real::offset_min,offset_max
        real::dx_v,dz_v
!        real::pml_thick!(m)

        real::x_bound_l,x_bound_r
        real::z_bound_u,z_bound_d!**zy

        real::sx0,sz0,rz0!(m)
        real::dsx!(m)
        real::f0!(Hz)
        real::dt!(s)
        !*Dummy variables(Out)*
        integer::nx_bound_l,nx_bound_r,nz_bound_u,nz_bound_d
        integer::nsx0,nsz0,nrz0
        integer::currshot_xmin,currshot_xmax
        integer::ndsx,nsx
        integer::currshot_no
        character(len=256)::currshot_name
        real,allocatable::vz(:,:)
        real,allocatable::VTI_epsilon(:,:)
        real,allocatable::VTI_delta(:,:)
        !Local variables
        !*Variables for MPI interface*
        integer::myid,nproc,ierr
        integer::status(MPI_STATUS_SIZE)  
        !*Other local variables*
        integer::i,j,k,err

        !Initiallization of MPI interface
        call MPI_init(ierr)
        call MPI_comm_rank(MPI_COMM_WORLD,myid,ierr)
        call MPI_comm_size(MPI_COMM_WORLD,nproc,ierr)
        call MPI_barrier(MPI_COMM_WORLD,ierr)

        !Read parameter_card
        call read_par(par_fn,vz_fn,epsilon_fn,delta_fn,shot_fn1,shot_fn2,&
                      nshot,offset_min,offset_max,sx0,sz0,rz0,dsx,nx_v,nz_v,&
                      dx_v,dz_v,nt,dt,f0,x_bound_l,x_bound_r,z_bound_u,z_bound_d,&
                      fd_order_explicit,order_pml_1st,order_pml_2nd)

        !Some test              
!        if(myid==0)then
!
!        write(*,'(A)')'vz_fn'
!        write(*,'(A)')vz_fn
!        write(*,'(A)')'epsilon_fn'
!        write(*,'(A)')epsilon_fn
!
!        write(*,'(A)')'delta_fn'
!        write(*,'(A)')delta_fn
!
!        write(*,'(A)')'shot_fn1,shot_fn2'
!        write(*,'(A)')shot_fn1,shot_fn2
!
!        write(*,*)'nshot'
!        write(*,*)nshot
!
!        write(*,*)'sx0,sz0,rz0'
!        write(*,*)sx0,sz0,rz0
!
!        write(*,*)'dsx'
!        write(*,'(F)')dsx
!
!
!        write(*,*)'nx_v,nz_v'
!        write(*,*)nx_v,nz_v
!        
!        write(*,*)'dx_v,dz_v'
!        write(*,*)dx_v,dz_v
!
!        write(*,*)'nt'
!        write(*,'(I)')nt
!
!        write(*,*)'dt'
!        write(*,'(F)')dt
!
!        write(*,*)'f0'
!        write(*,'(F)')f0
!
!        write(*,*)'pml_thick'
!        write(*,*)x_bound_l,x_bound_r,z_bound_u,z_bound_d
!        
!        write(*,*)'fd_order_explicit'
!        write(*,'(I)')fd_order_explicit
!
!        write(*,*)"order_pml_1st"
!        write(*,'(I)')order_pml_1st
!
!        write(*,*)'order_pml_2nd'
!        write(*,'(I)')order_pml_2nd
!      endif



        allocate(vz(nx_v,nz_v),STAT=err)
        allocate(VTI_delta(nx_v,nz_v),STAT=err)
        allocate(VTI_epsilon(nx_v,nz_v),STAT=err)
    
        !Initiallization
        vz=0.0
        VTI_epsilon=0.0
        VTI_delta=0.0
        !*Transfrom coordinates to grid_num for the conveniency of computation*!

              nx_bound_l=nint(x_bound_l/dx_v)
              nx_bound_r=nint(x_bound_r/dx_v)
        
              nz_bound_u=nint(z_bound_u/dz_v)
              nz_bound_d=nint(z_bound_d/dz_v)


!          write(*,*)'nx_bound_l,nx_bound_r,nz_bound_u,nz_bound_d'
!          write(*,*)nx_bound_l,nx_bound_r,nz_bound_u,nz_bound_d

        nsx0=nint(sx0/dx_v)+1!**zy                                    
        nsz0=nint(sz0/dz_v)+1!**zy                                    
        nrz0=nint(rz0/dz_v)+1!**zy
        !*Read parameters from files*!
        call read_data(vz_fn,epsilon_fn,delta_fn,&
                       vz,VTI_epsilon,VTI_delta,&
                       nx_v,nz_v,myid)
        

        !Multishot modelling
        if (myid==0)write(*,*)'==========2D_VTI_qusi_acoustic_modelling &
        began=========='


        do i=1+myid,nshot,nproc
          nsx=nsx0+nint((i-1)*dsx/dx_v)
          currshot_xmin=nsx+nint(offset_min/dx_v)
          currshot_xmax=nsx+nint(offset_max/dx_v)
          currshot_no=i


          write(currshot_name,'(I4)')i

          call single_shot_modelling(vz,VTI_epsilon,VTI_delta,&
               currshot_xmin,currshot_xmax,nsx,nsz0,nrz0,nx_v,nz_v,dx_v,dz_v,&
               nt,dt,f0,nx_bound_l,nx_bound_r,nz_bound_u,nz_bound_d,&
               fd_order_explicit,order_pml_1st,order_pml_2nd,&
               shot_fn1,shot_fn2,currshot_name,currshot_no,myid)
        enddo

        call MPI_barrier(MPI_COMM_WORLD,ierr)
        if (myid==0)then
          write(*,*)'Now merging shot_record files begins'
          call merge_shot_files(nshot,shot_fn1,shot_fn2,&
              currshot_xmin,currshot_xmax,nt)
        write(*,*)'==========2D_VTI_qusi_acoustic_modelling end&
        =========='
        endif
        
        call MPI_finalize(ierr)


        deallocate(vz,STAT=err)
        deallocate(VTI_delta,STAT=err)
        deallocate(VTI_epsilon,STAT=err)
        
        return
        end subroutine multi_shot_modelling


        subroutine single_shot_modelling(vz,VTI_epsilon,VTI_delta,&
                    currshot_xmin,currshot_xmax,nsx,nsz0,nrz0,nx_v,nz_v,dx,dz,&
                    nt,dt,f0,nx_bound_l,nx_bound_r,nz_bound_u,nz_bound_d,&
                    fd_order_explicit,order_pml_1st,order_pml_2nd,&
                    shot_fn1,shot_fn2,currshot_name,currshot_no,myid)
        
        use constant
        implicit none

        !Dummy variables
        real::vz(nx_v,nz_v)
        real::VTI_delta(nx_v,nz_v)
        real::VTI_epsilon(nx_v,nz_v)

        integer::nsx,nsz0,nrz0
        integer::currshot_xmin,currshot_xmax
        integer::nx_v,nz_v,nt,&
                 nx_bound_l,nx_bound_r,&
                 nz_bound_u,nz_bound_d
        integer::fd_order_explicit,order_pml_1st,order_pml_2nd
        integer::myid
        real::dx,dz,dt,f0
        integer::currshot_no
        character(len=256)::shot_fn1,shot_fn2,currshot_name
        !Local variables

        !*Buffers for finite_difference*
        real,allocatable::u1(:,:)
        real,allocatable::u2(:,:)
        real,allocatable::q1(:,:)
        real,allocatable::q2(:,:)
        real,allocatable::record(:,:)
        real,allocatable::fd_coe_explicit(:)
        real,allocatable::fd_coe_pml_2nd(:)
        real,allocatable::fd_coe_pml_1st(:)
        real,allocatable::currshot_vz(:,:)!**zy
        real,allocatable::currshot_epsilon(:,:)!**zy
        real,allocatable::currshot_delta(:,:)!**zy

        !**Buffers for PML wavefield**
      
        !***Buffers for p wave***
        real,allocatable::u1_1_x(:,:),u1_2_x(:,:),&
                          u2_1_x(:,:),u2_2_x(:,:),&
                          u2_tmp_1_x(:,:),u2_tmp_2_x(:,:),&
                          u3_1_x(:,:),u3_2_x(:,:),&
                          u1_1_z(:,:),u1_2_z(:,:),&
                          u2_1_z(:,:),u2_2_z(:,:),&
                          u2_tmp_1_z(:,:),u2_tmp_2_z(:,:),&
                          u3_1_z(:,:),u3_2_z(:,:)
                          
        !***Buffers  for q wave***
        real,allocatable::u1_1_x_q(:,:),u1_2_x_q(:,:),&
                          u2_1_x_q(:,:),u2_2_x_q(:,:),&
                          u2_tmp_1_x_q(:,:),u2_tmp_2_x_q(:,:),&
                          u3_1_x_q(:,:),u3_2_x_q(:,:),&
                          u1_1_z_q(:,:),u1_2_z_q(:,:),&
                          u2_1_z_q(:,:),u2_2_z_q(:,:),&
                          u2_tmp_1_z_q(:,:),u2_tmp_2_z_q(:,:),&
                          u3_1_z_q(:,:),u3_2_z_q(:,:)
    
        !*Variables wavefield extrapolation*
        integer::currshot_range,currshot_range_all,nz,&
                  currshot_nsx,currshot_nsz
        !*Other local variables*
        integer::i,i_ori,j,k,it,iit,err
        integer::nwt
        character(len=256)::snap_fn=&
        '/sldata2/tongji_data/hess_snap/snap'

        character(len=256)::currtime
        integer::counter=1


        !Computing currshot range
        currshot_range=(currshot_xmax-currshot_xmin)+1
        currshot_range_all=currshot_range+(nx_bound_l+nx_bound_r)-2  
        nz=nz_v+(nz_bound_u+nz_bound_d)-2!**zy
        !compute the length of wavelet in time dimension
!        nwt=nint(2.0*1.2/(f0*dt))+1
        nwt=nint(2.0*0.07/dt)+1


        !Allocate memories for buffers
        allocate(u1(currshot_range_all+10,nz+10),STAT=err)
        allocate(u2(currshot_range_all+10,nz+10),STAT=err)
        allocate(q1(currshot_range_all+10,nz+10),STAT=err)
        allocate(q2(currshot_range_all+10,nz+10),STAT=err)
        allocate(record(nt,currshot_range),STAT=err)!**new


        allocate(fd_coe_explicit(fd_order_explicit/2),STAT=err)
        allocate(fd_coe_pml_2nd(order_pml_2nd/2),STAT=err)
        allocate(fd_coe_pml_1st((order_pml_1st+1)/2),STAT=err)
        allocate(currshot_vz(currshot_range_all,nz))
        allocate(currshot_epsilon(currshot_range_all,nz))
        allocate(currshot_delta(currshot_range_all,nz))
        !*Buffers for PML wavefield*
        !**Buffers for p wave**
        !***X direction***
        allocate(u1_1_x(nx_bound_l+nx_bound_r,nz))
        allocate(u1_2_x(nx_bound_l+nx_bound_r,nz))
        allocate(u2_1_x(nx_bound_l+nx_bound_r,nz))
        allocate(u2_2_x(nx_bound_l+nx_bound_r,nz))
        allocate(u2_tmp_1_x(nx_bound_l+nx_bound_r,nz))
        allocate(u2_tmp_2_x(nx_bound_l+nx_bound_r,nz))
        allocate(u3_1_x(nx_bound_l+nx_bound_r,nz))
        allocate(u3_2_x(nx_bound_l+nx_bound_r,nz))
        !***Z direction***
        allocate(u1_1_z(currshot_range_all,nz_bound_u+nz_bound_d))
        allocate(u1_2_z(currshot_range_all,nz_bound_u+nz_bound_d))
        allocate(u2_1_z(currshot_range_all,nz_bound_u+nz_bound_d))
        allocate(u2_2_z(currshot_range_all,nz_bound_u+nz_bound_d))
        allocate(u2_tmp_1_z(currshot_range_all,nz_bound_u+nz_bound_d))
        allocate(u2_tmp_2_z(currshot_range_all,nz_bound_u+nz_bound_d))
        allocate(u3_1_z(currshot_range_all,nz_bound_u+nz_bound_d))
        allocate(u3_2_z(currshot_range_all,nz_bound_u+nz_bound_d))
        !**Buffers for q wave**
        
        !***X direction***
        allocate(u1_1_x_q(nx_bound_l+nx_bound_r,nz))
        allocate(u1_2_x_q(nx_bound_l+nx_bound_r,nz))
        allocate(u2_1_x_q(nx_bound_l+nx_bound_r,nz))
        allocate(u2_2_x_q(nx_bound_l+nx_bound_r,nz))
        allocate(u2_tmp_1_x_q(nx_bound_l+nx_bound_r,nz))
        allocate(u2_tmp_2_x_q(nx_bound_l+nx_bound_r,nz))
        allocate(u3_1_x_q(nx_bound_l+nx_bound_r,nz))
        allocate(u3_2_x_q(nx_bound_l+nx_bound_r,nz))
        !***Z direction***
        allocate(u1_1_z_q(currshot_range_all,nz_bound_u+nz_bound_d))
        allocate(u1_2_z_q(currshot_range_all,nz_bound_u+nz_bound_d))
        allocate(u2_1_z_q(currshot_range_all,nz_bound_u+nz_bound_d))
        allocate(u2_2_z_q(currshot_range_all,nz_bound_u+nz_bound_d))
        allocate(u2_tmp_1_z_q(currshot_range_all,nz_bound_u+nz_bound_d))
        allocate(u2_tmp_2_z_q(currshot_range_all,nz_bound_u+nz_bound_d))
        allocate(u3_1_z_q(currshot_range_all,nz_bound_u+nz_bound_d))
        allocate(u3_2_z_q(currshot_range_all,nz_bound_u+nz_bound_d))
    
        if (err.ne.0)then
          write(*,*)'Allocate memories fails,please check'
          stop
        endif
    
        !Initiallizations
        !*'Zeros'the buffer*
        u1=0.0
        u2=0.0
        q1=0.0
        q2=0.0
        record=0.0
          
        u1_1_x=0.0
        u1_2_x=0.0
        u2_1_x=0.0
        u2_2_x=0.0
        u2_tmp_1_x=0.0
        u2_tmp_2_x=0.0
        u3_1_x=0.0
        u3_2_x=0.0
    
        u1_1_z=0.0
        u1_2_z=0.0
        u2_1_z=0.0
        u2_2_z=0.0
        u2_tmp_1_z=0.0
        u2_tmp_2_z=0.0
        u3_1_z=0.0
        u3_2_z=0.0
    
    
        u1_1_x_q=0.0
        u1_2_x_q=0.0
        u2_1_x_q=0.0
        u2_2_x_q=0.0
        u2_tmp_1_x_q=0.0
        u2_tmp_2_x_q=0.0
        u3_1_x_q=0.0
        u3_2_x_q=0.0
    
        u1_1_z_q=0.0
        u1_2_z_q=0.0
        u2_1_z_q=0.0
        u2_2_z_q=0.0
        u2_tmp_1_z_q=0.0
        u2_tmp_2_z_q=0.0
        u3_1_z_q=0.0
        u3_2_z_q=0.0
    
    
        fd_coe_explicit=0.0
        fd_coe_pml_2nd=0.0
        fd_coe_pml_1st=0.0
        currshot_vz=0.0
        currshot_epsilon=0.0
        currshot_delta=0.0
    
!        open(unit=11,file=trim(adjustl(snap_fn))//'.bin',&
!        form='unformatted',access='direct',status='replace',recl=nz_v)
        !Obtain the coefficients for explicit FD extrapolation 
        call coefficient_2nd(fd_order_explicit,fd_coe_explicit)  
        call coefficient_2nd(order_pml_2nd,fd_coe_pml_2nd)
        call coefficient_1st(order_pml_1st,fd_coe_pml_1st)
    
        !******************************
        !    write(*,*)'test coe'      ! 
        !    write(*,*)fd_coe_explicit !
        !    write(*,*)fd_coe_pml_2nd  !
        !    write(*,*)fd_coe_pml_1st  !
        !    pause                     !
        !*****************************

        !Get currshot velocity&anisotropic parameters
        call get_currshot_parameters(vz,VTI_epsilon,VTI_delta,&
             currshot_vz,currshot_epsilon,currshot_delta,nx_v,nz_v,&
             currshot_range,currshot_range_all,currshot_xmin,currshot_xmax,&
             nz,nx_bound_l,nx_bound_r,nz_bound_u,nz_bound_d,nsx,nsz0,&
             currshot_nsx,currshot_nsz)

        !Test
        if (myid==0)then
          write(*,*)'currshot_range,currshot_range_all,currshot_nsx,nz,nsz0'
          write(*,*)currshot_range,currshot_range_all,currshot_nsx,nz,nsz0
        endif


        do j=1,nsz0+nz_bound_u+6
            currshot_epsilon(:,j)=0.0
            currshot_delta(:,j)=0.0
        enddo


    
        do it=1,nt+nint((nwt-1)/2.0)
          if(modulo(it,2)==1)then


            if(it<=nwt)then  
              u2(currshot_nsx,currshot_nsz)=&
              u2(currshot_nsx,currshot_nsz)+&
              exp(-(pi*f0*(it-(0.07)/dt)*dt)**2)*(1-2*(pi*f0*(it-(0.07)/dt)*dt)**2)
            endif


            call extrapolation_one_step_1(currshot_vz,currshot_epsilon,currshot_delta,&
                                           currshot_nsx,nsz0,nrz0,currshot_range,nz_v,currshot_range_all,&
                                          nz,dx,dz,nt,dt,f0,nx_bound_l,nx_bound_r,nz_bound_u,nz_bound_d,&
                                          fd_order_explicit,order_pml_1st,order_pml_2nd,&
                                          q2,u1,u2,fd_coe_explicit,&
                                          fd_coe_pml_2nd,fd_coe_pml_1st,&
                                          u1_1_x,u1_2_x,u2_1_x,u2_2_x,&
                                          u2_tmp_1_x,u2_tmp_2_x,u3_1_x,u3_2_x,&
                                          u1_1_z,u1_2_z,u2_1_z,u2_2_z,&
                                          u2_tmp_1_z,u2_tmp_2_z,u3_1_z,u3_2_z&
                                          )
    

            call extrapolation_one_step_2(currshot_vz,currshot_epsilon,currshot_delta,&
                                           currshot_nsx,nsz0,nrz0,currshot_range,nz_v,currshot_range_all,&
                                          nz,dx,dz,nt,dt,f0,nx_bound_l,nx_bound_r,nz_bound_u,nz_bound_d,&
                                          fd_order_explicit,order_pml_1st,order_pml_2nd,&
                                          u2,q1,q2,fd_coe_explicit,&
                                          fd_coe_pml_2nd,fd_coe_pml_1st,&
                                          u1_1_x_q,u1_2_x_q,u2_1_x_q,u2_2_x_q,&
                                          u2_tmp_1_x_q,u2_tmp_2_x_q,u3_1_x_q,u3_2_x_q,&
                                          u1_1_z_q,u1_2_z_q,u2_1_z_q,u2_2_z_q,&
                                          u2_tmp_1_z_q,u2_tmp_2_z_q,u3_1_z_q,u3_2_z_q,&
                                          )!Modified for correction in tx


!            if(it>=nint((nwt-1)/2))record(it,1:currshot_range)=&
!                        u1(nx_bound_l+1:currshot_range+nx_bound_l,nrz0+nz_bound_u)  
   
            iit=it-nint((nwt-1)/2.0)
            if(iit>=1.and.nrz0+nz_bound_u>=6)then
                record(iit,1:currshot_range)=&
                u1(nx_bound_l+1:currshot_range+nx_bound_l,nrz0+nz_bound_u)   
            else if(iit>1.and.nrz0+nz_bound_u<6)then
                record(iit,1:currshot_range)=&
                u1(nx_bound_l+1:currshot_range+nx_bound_l,6)   
             endif

        else

            if(it<=nwt)then  
              u1(currshot_nsx,currshot_nsz)=&
              u1(currshot_nsx,currshot_nsz)+&
              exp(-(pi*f0*(it-(0.07)/dt)*dt)**2)*(1-2*(pi*f0*(it-(0.07)/dt)*dt)**2)
            endif

            call  extrapolation_one_step_1(currshot_vz,currshot_epsilon,currshot_delta,&
                                           currshot_nsx,nsz0,nrz0,currshot_range,nz_v,currshot_range_all,&
                                          nz,dx,dz,nt,dt,f0,nx_bound_l,nx_bound_r,nz_bound_u,nz_bound_d,&
                                          fd_order_explicit,order_pml_1st,order_pml_2nd,&
                                          q1,u2,u1,fd_coe_explicit,&
                                          fd_coe_pml_2nd,fd_coe_pml_1st,&
                                          u1_2_x,u1_1_x,u2_2_x,u2_1_x,&
                                          u2_tmp_2_x,u2_tmp_1_x,u3_2_x,u3_1_x,&
                                          u1_2_z,u1_1_z,u2_2_z,u2_1_z,&
                                          u2_tmp_2_z,u2_tmp_1_z,u3_2_z,u3_1_z&
                                          )
    

            call extrapolation_one_step_2(currshot_vz,currshot_epsilon,currshot_delta,&
                                          currshot_nsx,nsz0,nrz0,currshot_range,nz_v,currshot_range_all,&
                                          nz,dx,dz,nt,dt,f0,nx_bound_l,nx_bound_r,nz_bound_u,nz_bound_d,&
                                          fd_order_explicit,order_pml_1st,order_pml_2nd,&
                                          u1,q2,q1,fd_coe_explicit,&
                                          fd_coe_pml_2nd,fd_coe_pml_1st,&
                                          u1_2_x_q,u1_1_x_q,u2_2_x_q,u2_1_x_q,&
                                          u2_tmp_2_x_q,u2_tmp_1_x_q,u3_2_x_q,u3_1_x_q,&
                                          u1_2_z_q,u1_1_z_q,u2_2_z_q,u2_1_z_q,&
                                          u2_tmp_2_z_q,u2_tmp_1_z_q,u3_2_z_q,u3_1_z_q,&
                                          )

    
               iit=it-nint((nwt-1)/2.0)
               if(iit>=1.and.nrz0+nz_bound_u>=6)then
                record(iit,1:currshot_range)=&
                u2(nx_bound_l+1:currshot_range+nx_bound_l,nrz0+nz_bound_u)   
               else if(iit>1.and.nrz0+nz_bound_u<6)then
                record(iit,1:currshot_range)=&
                u2(nx_bound_l+1:currshot_range+nx_bound_l,6)   
               endif
		  
		  
		  
          endif




          if (myid==0)then
          if(modulo(it,201)==0)then
              write(currtime,'(I4)')nint((it-1)*dt*1.0E3)
              write(*,*)'shot is',trim(currshot_name),'  myid=',myid,'nt=',nt,'it=',it
          endif

          endif
        enddo
    
    !*************************************************************************************  
      
        write(*,*)'shot  ',trim(adjustl(currshot_name)),&
        '  extrapolation finished,and now is writing disk'

        call write_currshot_disk(record,shot_fn1,currshot_name,currshot_no,&
                                 nsx,nsz0,currshot_range,currshot_xmin,&
                                 currshot_xmax,nt,dx,dz,dt)

        write(*,*)'shot  ',trim(adjustl(currshot_name)),'  is done'
    
        deallocate(u1,STAT=err)
        deallocate(u2,STAT=err)
        deallocate(fd_coe_explicit,STAT=err)
        deallocate(fd_coe_pml_2nd,STAT=err)
        deallocate(fd_coe_pml_1st,STAT=err)
        deallocate(q1,STAT=err)
        deallocate(q2,STAT=err)
        deallocate(record,STAT=err)
        deallocate(u1_1_x,STAT=err)
        deallocate(u1_2_x,STAT=err)
        deallocate(u2_1_x,STAT=err)
        deallocate(u2_2_x,STAT=err)
        deallocate(u2_tmp_1_x,STAT=err)
        deallocate(u2_tmp_2_x,STAT=err)
        deallocate(u3_1_x,STAT=err)
        deallocate(u3_2_x,STAT=err)
        deallocate(u1_1_z,STAT=err)
        deallocate(u1_2_z,STAT=err)
        deallocate(u2_1_z,STAT=err)
        deallocate(u2_2_z,STAT=err)
        deallocate(u2_tmp_1_z,STAT=err)
        deallocate(u2_tmp_2_z,STAT=err)
        deallocate(u3_1_z,STAT=err)
        deallocate(u3_2_z,STAT=err)
        deallocate(u1_1_x_q,STAT=err)
        deallocate(u1_2_x_q,STAT=err)
        deallocate(u2_1_x_q,STAT=err)
        deallocate(u2_2_x_q,STAT=err)
        deallocate(u2_tmp_1_x_q,STAT=err)
        deallocate(u2_tmp_2_x_q,STAT=err)
        deallocate(u3_1_x_q,STAT=err)
        deallocate(u3_2_x_q,STAT=err)
        deallocate(u1_1_z_q,STAT=err)
        deallocate(u1_2_z_q,STAT=err)
        deallocate(u2_1_z_q,STAT=err)
        deallocate(u2_2_z_q,STAT=err)
        deallocate(u2_tmp_1_z_q,STAT=err)
        deallocate(u2_tmp_2_z_q,STAT=err)
        deallocate(u3_1_z_q,STAT=err)
        deallocate(u3_2_z_q,STAT=err)
        deallocate(currshot_vz,STAT=err)
        deallocate(currshot_epsilon,STAT=err)
        deallocate(currshot_delta,STAT=err)
    
    
        if (err.ne.0)then
          write(*,*)'Deallocate memories fails,please check'
          stop
        endif

        return
        end subroutine single_shot_modelling


























         






        






























        

         


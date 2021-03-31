		!2D Finite-difference Time-domain modeling of P-SV wave propagation
		!Modified to use high order approximation for both spatial derivative
		!and free surface boundary conditions  

		include 'Iso_multishot_fd_module_zy_m2.f90'
		include 'Iso_multishot_fd_subroutine_fsbc_zy_m2_pri.f90'
        include 'Iso_multishot_fd_subroutine_fsbc_zy_m2_pub.f90'

        program Iso_fd_elastic_modelling
        implicit none 
        character(len=256)::par_fn=&
        './iso_2d_modelling_sin_test.par'
        call multi_shot_modelling(par_fn)
        end program Iso_fd_elastic_modelling


        subroutine multi_shot_modelling(par_fn)
        use mpi
        use constant
		use allocbuff
		use deallocbuff
        implicit none

        !Dummy variables
        character(len=256),intent(in)::par_fn

        !Local variables
        character(len=256)::molmod,vp_fn,vs_fn,rho_fn,&
						topo_fn,snapx_fn,snapz_fn,&
						shot_fn1,shot_fn2,shot_fn3,shot_fn4

        integer::csno,nx,nz,nt,nshot,rgau,&
        		myid,nproc,ierr,i,err,&
				status(MPI_STATUS_SIZE)  

        real::ofsmin,ofsmax,dx,dz,&
			ppxn,ppxp,ppzn,ppzp,sx0,sz0,rz0,&
        	dsx,f0,dt,dtsnap,decay

        real,allocatable::vp(:,:),vs(:,:),&
						rho(:,:),topo(:)


		!Initiallization of MPI interface
        call MPI_init(ierr)
        call MPI_comm_rank(MPI_COMM_WORLD,myid,ierr)
        call MPI_comm_size(MPI_COMM_WORLD,nproc,ierr)
        call MPI_barrier(MPI_COMM_WORLD,ierr)


        !Read parameter_card
        call read_par(par_fn,molmod,vp_fn,vs_fn,rho_fn,topo_fn,&
				snapx_fn,snapz_fn,shot_fn1,shot_fn2,shot_fn3,shot_fn4,&
				nshot,ofsmin,ofsmax,sx0,sz0,rz0,dsx,nx,nz,dx,dz,nt,dt,&
				dtsnap,f0,decay,ppxn,ppxp,ppzn,ppzp,rgau)


		!Allocate buffers for common data
		call mallocbuff(err,nx,nz,topo,vp,vs,rho)

        !Read parameters from files!
        call read_data(vp_fn,vs_fn,rho_fn,topo_fn,vp,vs,rho,topo,nx,nz,myid)


        !Multishot modelling
        if (myid==0)&
		write(*,*)'************************************************'
		write(*,*)'*  2D isotropic elastic-wave modelling begins  *'
		write(*,*)'************************************************'

        do i=1+myid,nshot,nproc


          call single_shot_modelling(molmod,snapx_fn,snapz_fn,&
			   	shot_fn1,shot_fn2,shot_fn3,shot_fn4,i,rgau,nx,nz,nt,&
				dx,dz,dt,dtsnap,f0,decay,ofsmin,ofsmax,sx0,sz0,rz0,dsx,&
				ppxn,ppxp,ppzn,ppzp,rho,topo,vp,vs,myid)



        enddo

        call MPI_barrier(MPI_COMM_WORLD,ierr)

        if (myid==0)then

        	write(*,*)'Now merging shot_record files begins'

          	call merge_shot_files(nshot,shot_fn1,shot_fn2,&
              shot_fn3,shot_fn4,ofsmin,ofsmax,dx,nt)

			write(*,*)'**********************************************'
			write(*,*)'*  2D isotropic elastic-wave modelling ends  *'
			write(*,*)'**********************************************'

        endif
       

        call MPI_finalize(ierr)

      	call mdeallocbuff(err,nx,nz,topo,vp,vs,rho) 


        return
        end subroutine multi_shot_modelling


!******************************************************************************************
         subroutine single_shot_modelling(molmod,snapx_fn,snapz_fn,&
			   		shot_fn1,shot_fn2,shot_fn3,shot_fn4,csno,rgau,nx,nz,nt,&
					dx,dz,dt,dtsnap,f0,decay,ofsmin,ofsmax,sx0,sz0,rz0,dsx,&
					ppxn,ppxp,ppzn,ppzp,rho,topo,vp,vs,myid)

        use constant
		use allocbuff
		use deallocbuff

        implicit none

        !Dummy variables
        character(len=256),intent(in)::molmod,snapx_fn,snapz_fn,&
							shot_fn1,shot_fn2,shot_fn3,shot_fn4

        integer,intent(in)::nx,nz,nt,ppxn,ppxp,ppzn,ppzp,&
        				csno,myid,rgau

        real,intent(in)::dx,dz,dt,dtsnap,f0,decay,ofsmin,&
					ofsmax,sx0,sz0,rz0,dsx,vp(nx,nz),vs(nx,nz),&
       				rho(nx,nz),topo(nx)

        !Local variables
		character(len=256)::currt,currtsnap,csname
        integer::i,k,it,iit,err,nwt,nxe,nze,nxea,nzea,&
                isx0,isz0,isx,isz,irz0,isxe,isze,ixmine,&
				ixmaxe,idecay,idtsnap,npxn,npxp,npzn,npzp

        integer,allocatable::itopoe(:)
        real,allocatable::vx(:,:),vz(:,:),strxx(:,:),&
					strzz(:,:),strxz(:,:),vxx(:,:),&
					vxz(:,:),vzx(:,:),vzz(:,:),&
        			strxxx(:,:),strxxz(:,:),strzzx(:,:),&
					strzzz(:,:),strxzx(:,:),strxzz(:,:),&
        			der1(:,:),der2(:,:),vpe(:,:),&
        			vse(:,:),rhoe(:,:),record1(:,:),&
        			record2(:,:),topoe(:),mu(:,:),&
 					lamb(:,:),musg(:,:),musgef(:,:),&
 					lambsg(:,:),lambmusg(:,:),busg1(:,:),&
 					busg2(:,:),busg1ef(:,:),busg2ef(:,:),&
					spgxn(:),spgxp(:),spgzn(:),&
					spgzp(:),spg1xn(:),spg2xp(:),&
					spg1zn(:),spg2zp(:),gc2d(:,:)


		!Precompute some parameters for allocating buffers and
		!other purposes
		call spcalpa(csname,csno,rgau,ixmine,ixmaxe,isx0,isz0,irz0,&
				isz,isx,idecay,idtsnap,nz,npxn,npxp,npzn,npzp,nxe,nze,&
				nxea,nzea,ppxn,ppxp,ppzn,ppzp,sx0,sz0,rz0,dx,dz,dsx,&
				ofsmin,ofsmax,decay,dt,dtsnap)

		!Allocate buffers for data used in single-shot modelling
		call sallocbuff(err,nt,nxe,nxea,nzea,npxn,npxp,npzn,npzp,rgau,&
				itopoe,vx,vz,strxx,strzz,strxz,vxx,vxz,vzx,vzz,strxxx,&
				strxxz,strzzx,strzzz,strxzx,strxzz,der1,der2,vpe,vse,&
				rhoe,record1,record2,topoe,mu,lamb,musg,musgef,lambsg,&
				lambmusg,busg1,busg2,busg1ef,busg2ef,spgxn,spgxp,spgzn,&
				spgzp,spg1xn,spg2xp,spg1zn,spg2zp,gc2d)


		!Get currshot velocity&anisotropic parameters
        call get_currshot_parameters(vp,vs,rho,topo,vpe,vse,rhoe,topoe,&
				itopoe,gc2d,nx,nz,nxe,nze,nxea,nzea,ixmine,ixmaxe,npxn,&
				npxp,npzn,npzp,rgau,isx,isz,isxe,isze,dx,dz)


		!Print some basic information on the screen
!		write(*,*)'nxe,nxea,nze,nzea,isz,isx,isze,isxe'
!		write(*,*)nxe,nxea,nze,nzea,isz,isx,isze,isxe


		!Build staggered grids
 		call substagger(vpe,vse,rhoe,mu,lamb,musg,musgef,lambsg,&
				lambmusg,busg1,busg2,busg1ef,busg2ef,nxea,nzea,itopoe)


		!Building absorbing coefficients
		call csponge(spgxn,spgxp,spgzn,spgzp,spg1xn,spg2xp,spg1zn,spg2zp,&
				npxn,npxp,npzn,npzp,dx,dz)

		select case(molmod)


		case ('fsbc')

            call isoewm2dfsbc(csname,snapx_fn,snapz_fn,nt,idecay,idtsnap,&
					record1,record2,myid,irz0,rgau,gc2d,isxe,isze,nxe,nze,&
					nxea,nzea,npxn,npxp,npzn,npzp,f0,dx,dz,dt,vx,vz,strxx,&
					strzz,strxz,vxx,vxz,vzx,vzz,strxxx,strxxz,strzzx,strzzz,&
					strxzx,strxzz,der1,der2,spgxn,spgxp,spgzp,spg1xn,spg2xp,&
					spg2zp,vpe,vse,rhoe,itopoe,lamb,lambsg,lambmusg,mu,musg,&
					musgef,busg1,busg2,busg1ef,busg2ef)

		case ('abc')

   	         call isoewm2dabc(csname,snapx_fn,snapz_fn,nt,idecay,idtsnap,&
			 		record1,record2,myid,irz0,rgau,gc2d,isxe,isze,nxe,nze,&
					nxea,nzea,npxn,npxp,npzn,npzp,f0,dx,dz,dt,vx,vz,strxx,&
					strzz,strxz,vxx,vxz,vzx,vzz,strxxx,strxxz,strzzx,strzzz,&
					strxzx,strxzz,der1,der2,spgxn,spgxp,spgzn,spgzp,spg1xn,spg2xp,&
					spg1zn,spg2zp,vpe,vse,rhoe,itopoe,lamb,lambsg,lambmusg,mu,musg,&
					musgef,busg1,busg2,busg1ef,busg2ef)


		case default

			write(*,*)'Wrong modelling mode flag, please check'
			stop

		end select 

    
    !*************************************************************************************  
      
        write(*,*)'shot  ',trim(adjustl(csname)),&
        '  extrapolation finished,and now is writing disk'

        call write_currshot_disk(record1,shot_fn1,record2,shot_fn3,&
								csname,csno,isx,isz,nxe,ixmine,ixmaxe,&
								nt,dx,dz,dt)


        write(*,*)'shot  ',trim(adjustl(csname)),'  is done'
   

		call sdeallocbuff(err,itopoe,vx,vz,strxx,strzz,strxz,vxx,&
					vxz,vzx,vzz,strxxx,strxxz,strzzx,strzzz,strxzx,&
					strxzz,der1,der2,vpe,vse,rhoe,record1,record2,&
					topoe,mu,lamb,musg,musgef,lambsg,lambmusg,busg1,&
 					busg2,busg1ef,busg2ef,spgxn,spgxp,spgzn,spgzp,&
					spg1xn,spg2xp,spg1zn,spg2zp,gc2d)


        return
        end subroutine single_shot_modelling


























         






        






























        

         


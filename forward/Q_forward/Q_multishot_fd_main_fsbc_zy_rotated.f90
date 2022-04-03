		!2D Finite-difference Time-domain modeling of P-SV wave propagation
		!Modified to use high order approximation for both spatial derivative
		!and free surface boundary conditions  
		!Modified to a rotated-staggered version
		!Modified to add a separation method
		!Separating function caluxabs into calder and caluxabs
		!Removing useless functions in subroutines(fsbc concerned)
		!Adding the bond transform to adapt TTI case
		!Reducing computational cost by introducing more buffers
		!for specfic derivatives
		!Modified for a qp case
		!Modified to modelling visco_elastic wave propagation&&
		!using rotated staggered grid

        include './Q_multishot_fd_subroutine_fsbc_zy_rotated.f90'
        program Iso_fd_elastic_modelling
        implicit none 
        character(len=256)::par_fn=&
        './isoq_2d_modelling_simple_hess.par'
        call multi_shot_modelling(par_fn)
        end program Iso_fd_elastic_modelling

        subroutine multi_shot_modelling(par_fn)
        use mpi
        use constant
        implicit none
        !Dummy variables
        !*Parameter_card_name*
        character(len=256)::par_fn
        !*Read from parameter_card*
        character(len=256)::vp_fn
        character(len=256)::vs_fn
        character(len=256)::rho_fn
        character(len=256)::qfp_fn
        character(len=256)::qfs_fn
        character(len=256)::shot_fn1
        character(len=256)::shot_fn2
        character(len=256)::shot_fn3
        character(len=256)::shot_fn4
        integer::nx,nz
        integer::nt
        integer::fdo1
        integer::fdo2
        integer::nshot
        real::offset_min,offset_max
        real::dx,dz
        real::ppxn,ppxp
        real::ppzn,ppzp!**zy
        real::sx0,sz0,rz0!(m)
        real::dsx!(m)
        real::f0!(Hz)
        real::dt!(s)
        !*Dummy variables(Out)*
        integer::npxn,npxp,npzn,npzp
        integer::isx0,nsz0,nrz0
        integer::ixmine,ixmaxe
        integer::isx,nsz
        integer::csno
        character(len=256)::csname
        real,allocatable::vp(:,:)
        real,allocatable::vs(:,:)
        real,allocatable::rho(:,:)
        real,allocatable::qfp(:,:)!quality factor
        real,allocatable::qfs(:,:)!quality factor
        real,allocatable::topo(:)

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
        call read_par(par_fn,vp_fn,vs_fn,rho_fn,qfp_fn,qfs_fn,&
					shot_fn1,shot_fn2,shot_fn3,shot_fn4,&
					nshot,offset_min,offset_max,sx0,sz0,rz0,&
					dsx,nx,nz,dx,dz,nt,dt,f0,ppxn,ppxp,ppzn,ppzp,&
					fdo2,fdo1)

        allocate(vp(nz,nx),STAT=err)
        allocate(vs(nz,nx),STAT=err)
        allocate(rho(nz,nx),STAT=err)
        allocate(qfp(nz,nx),STAT=err)
        allocate(qfs(nz,nx),STAT=err)
        allocate(topo(nx),STAT=err)

        !Initiallization
        vp=0.0
        vs=0.0
        rho=0.0
		qfp=0.0
		qfs=0.0
		topo=0.0


        !*Transfrom coordinates to grid_num for the conveniency of computation*!
        npxn=nint(ppxn/dx)+1
        npxp=nint(ppxp/dx)+1
        npzn=nint(ppzn/dz)+1
        npzp=nint(ppzp/dz)+1

        isx0=nint(sx0/dx)+1!**zy      
        nsz0=nint(sz0/dz)+1!**zy pay attenuation                              
        nrz0=nint(rz0/dz)+1!**zy


        !*Read parameters from files*!
        call read_data(vp_fn,vs_fn,rho_fn,qfp_fn,qfs_fn,&
				vp,vs,rho,qfp,qfs,nx,nz,myid)

        !Multishot modelling
        if (myid==0)write(*,*)'==========2D_Iso_qusi_acoustic_modelling &
        began=========='


        do i=1+myid,nshot,nproc

          isx=isx0+nint((i-1)*dsx/dx)
		  nsz=nsz0
          ixmine=isx+nint(offset_min/dx)
          ixmaxe=isx+nint(offset_max/dx)
          csno=i

          write(csname,'(I4)')i

          call single_shot_modelling(vp,vs,rho,qfp,qfs,topo,&
               			ixmine,ixmaxe,isx,nsz0,nrz0,nx,nz,dx,dz,&
               			nt,dt,f0,npxn,npxp,npzn,npzp,fdo2,fdo1,&
			   			shot_fn1,shot_fn2,shot_fn3,shot_fn4,&
			   			csname,csno,myid)


        enddo

        call MPI_barrier(MPI_COMM_WORLD,ierr)

        if (myid==0)then
          write(*,*)'Now merging shot_record files begins'
          call merge_shot_files(nshot,shot_fn1,shot_fn2,&
              shot_fn3,shot_fn4,ixmine,ixmaxe,nt)

        write(*,*)'==========2D_Iso_acoustic_modelling end&
        =========='
        endif
        
        call MPI_finalize(ierr)


        deallocate(vp,STAT=err)
        deallocate(vs,STAT=err)
        deallocate(rho,STAT=err)
        deallocate(qfp,STAT=ierr)
        deallocate(qfs,STAT=ierr)
        deallocate(topo,STAT=ierr)


        return
        end subroutine multi_shot_modelling

!******************************************************************************************
		subroutine  single_shot_modelling(vp,vs,rho,qfp,qfs,topo,&
								ixmine,ixmaxe,isx,isz,nrz0,nx,nz,dx,dz,&
								nt,dt,f0,npxn,npxp,npzn,npzp,fdo2,fdo1,&
			   					shot_fn1,shot_fn2,shot_fn3,shot_fn4,&
								csname,csno,myid)


        use constant
        implicit none
        !Dummy variables
        character(len=256)::shot_fn1,shot_fn2,&
							shot_fn3,shot_fn4,&
							csname
        integer::ixmine,ixmaxe,&
        		nx,nz,nt,npxn,npxp,npzn,npzp,&
        		fdo2,fdo1,csno,myid

        real::dx,dz,dt,f0,vp(nz,nx),vs(nz,nx),&
 			 rho(nz,nx),qfp(nz,nx),qfs(nz,nx),topo(nx)

        !Local variables
        !*Buffers for finite_difference*
        real,allocatable::vx(:,:)
        real,allocatable::vz(:,:)
        real,allocatable::strxx(:,:)
        real,allocatable::strzz(:,:)
        real,allocatable::strxz(:,:)
        real,allocatable::rxx(:,:)
        real,allocatable::rzz(:,:)
        real,allocatable::rxz(:,:)
        real,allocatable::vxx(:,:)
        real,allocatable::vxz(:,:)
        real,allocatable::vzx(:,:)
        real,allocatable::vzz(:,:)

        real,allocatable::strxxx(:,:)
        real,allocatable::strxxz(:,:)
        real,allocatable::strxxr(:,:)

        real,allocatable::strzzx(:,:)
        real,allocatable::strzzz(:,:)
        real,allocatable::strzzr(:,:)

        real,allocatable::strxzx(:,:)
        real,allocatable::strxzz(:,:)
        real,allocatable::strxzr(:,:)

        real,allocatable::rxxx(:,:)
        real,allocatable::rxxz(:,:)
        real,allocatable::rxxr(:,:)

        real,allocatable::rzzx(:,:)
        real,allocatable::rzzz(:,:)
        real,allocatable::rzzr(:,:)

        real,allocatable::rxzx(:,:)
        real,allocatable::rxzz(:,:)
        real,allocatable::rxzr(:,:)


        real,allocatable::der1(:,:)
       	real,allocatable::der2(:,:)
       	real,allocatable::der3(:,:)
       	real,allocatable::der4(:,:)
        real,allocatable::vpe(:,:)
        real,allocatable::vse(:,:)
        real,allocatable::rhoe(:,:)
        real,allocatable::qfpe(:,:)
        real,allocatable::qfse(:,:)
 		real,allocatable::busge(:,:)!1/rho
 		real,allocatable::lambe(:,:)
 		real,allocatable::mue(:,:)
 		real,allocatable::lambmue(:,:)
 		real,allocatable::taue(:,:)
 		real,allocatable::taupstine(:,:)
 		real,allocatable::tausstine(:,:)
 		real,allocatable::taustsse(:,:)
 		real,allocatable::rtaustsse(:,:)
 		real,allocatable::taupstinsse(:,:)
 		real,allocatable::tausstinsse(:,:)
        real,allocatable::record1(:,:)
        real,allocatable::record2(:,:)
        real,allocatable::topoe(:)
        integer,allocatable::itopoe(:)
		real,allocatable::spgxn(:)
		real,allocatable::spgxp(:)
		real,allocatable::spgzn(:)
		real,allocatable::spgzp(:)
		real,allocatable::spg1xn(:)
		real,allocatable::spg2xp(:)
		real,allocatable::spg1zn(:)
		real,allocatable::spg2zp(:)


        !*Variables wavefield extrapolation*
        integer::nxe,nze,nxea,nzea,&
                 isx,isz,nrz0,isxe,isze
        !*Other local variables*
        integer::i,i_ori,j,k,it,iit,err
        integer::nwt,decay
		character(len=256)::currt,currtsnap
		character(len=256)::snapx_fn='./snap/iso_vx_elastic_'
		character(len=256)::snapz_fn='./snap/iso_vz_elastic_'

        !Computing currshot range
        nxe=(ixmaxe-ixmine)+1
        nxea=nxe+(npxn+npxp)
		nze=nz
        nzea=nz+(npzn+npzp)! nz effective all
        !compute the length of wavelet in time dimension
        nwt=2*nint(1.2/(f0*dt))
		decay=nint(nwt/2.0)

        !Allocate memories for buffers
        allocate(vx(-4:nzea+5,-4:nxea+5),STAT=err)
        allocate(vz(-4:nzea+5,-4:nxea+5),STAT=err)
        allocate(strxx(-4:nzea+5,-4:nxea+5),STAT=err)
        allocate(strzz(-4:nzea+5,-4:nxea+5),STAT=err)
        allocate(strxz(-4:nzea+5,-4:nxea+5),STAT=err)
        allocate(rxx(nzea,nxea),STAT=err)
        allocate(rzz(nzea,nxea),STAT=err)
        allocate(rxz(nzea,nxea),STAT=err)
        allocate(vxx(0:nzea,nxea),STAT=err)
        allocate(vxz(0:nzea,nxea),STAT=err)
        allocate(vzx(0:nzea,nxea),STAT=err)
        allocate(vzz(0:nzea,nxea),STAT=err)
        allocate(strxxx(0:nzea,nxea),STAT=err)
        allocate(strxxz(0:nzea,nxea),STAT=err)
        allocate(strxxr(0:nzea,nxea),STAT=err)
        allocate(strzzx(0:nzea,nxea),STAT=err)
        allocate(strzzz(0:nzea,nxea),STAT=err)
        allocate(strzzr(0:nzea,nxea),STAT=err)
        allocate(strxzx(0:nzea,nxea),STAT=err)
        allocate(strxzz(0:nzea,nxea),STAT=err)
        allocate(strxzr(0:nzea,nxea),STAT=err)
        allocate(rxxx(0:nzea,nxea),STAT=err)
        allocate(rxxz(0:nzea,nxea),STAT=err)
        allocate(rxxr(0:nzea,nxea),STAT=err)
        allocate(rzzx(0:nzea,nxea),STAT=err)
        allocate(rzzz(0:nzea,nxea),STAT=err)
        allocate(rzzr(0:nzea,nxea),STAT=err)
        allocate(rxzx(0:nzea,nxea),STAT=err)
        allocate(rxzz(0:nzea,nxea),STAT=err)
        allocate(rxzr(0:nzea,nxea),STAT=err)
        allocate(der1(nzea,nxea),STAT=err)
        allocate(der2(nzea,nxea),STAT=err)
        allocate(der3(nzea,nxea),STAT=err)
        allocate(der4(nzea,nxea),STAT=err)
        allocate(vpe(nzea,nxea),STAT=err)
        allocate(vse(nzea,nxea),STAT=err)
        allocate(rhoe(nzea,nxea),STAT=err)
        allocate(qfpe(nzea,nxea),STAT=err)
        allocate(qfse(nzea,nxea),STAT=err)
		allocate(busge(nzea,nxea),STAT=err)
        allocate(lambe(nzea,nxea),STAT=err)
        allocate(mue(nzea,nxea),STAT=err)
        allocate(lambmue(nzea,nxea),STAT=err)
		allocate(taue(nzea,nxea),STAT=err)
		allocate(taupstine(nzea,nxea),STAT=err)
		allocate(tausstine(nzea,nxea),STAT=err)
		allocate(taustsse(nzea,nxea),STAT=err)
		allocate(rtaustsse(nzea,nxea),STAT=err)
		allocate(taupstinsse(nzea,nxea),STAT=err)
		allocate(tausstinsse(nzea,nxea),STAT=err)
        allocate(topoe(nxea),STAT=err)
        allocate(itopoe(nxea),STAT=err)
        allocate(record1(nt,nxe),STAT=err)
        allocate(record2(nt,nxe),STAT=err)
		allocate(spgxn(npxn+1),STAT=err)
		allocate(spgxp(npxp+1),STAT=err)
		allocate(spgzn(npzn+1),STAT=err)
		allocate(spgzp(npzp+1),STAT=err)
		allocate(spg1xn(npxn+1),STAT=err)
		allocate(spg2xp(npxp+1),STAT=err)
		allocate(spg1zn(npzn+1),STAT=err)
		allocate(spg2zp(npzp+1),STAT=err)

        !Initiallizations
        !*'Zeros'the buffer*
        vx=0.0
        vz=0.0
        strxx=0.0
        strzz=0.0
        strxz=0.0
        rxx=0.0
        rzz=0.0
        rxz=0.0
        vxx=0.0
        vxz=0.0
        vzx=0.0
        vzz=0.0
        strxxx=0.0
        strxxz=0.0
        strxxr=0.0
        strzzx=0.0
        strzzz=0.0
        strzzr=0.0
        strxzx=0.0
        strxzz=0.0
        strxzr=0.0
        rxxx=0.0
        rxxz=0.0
        rxxr=0.0
        rzzx=0.0
        rzzz=0.0
        rzzr=0.0
        rxzx=0.0
        rxzz=0.0
        rxzr=0.0
        der1=0.0
        der2=0.0
        der3=0.0
        der4=0.0
        vpe=0.0
        vse=0.0
        rhoe=0.0
        qfpe=0.0
        qfse=0.0
		busge=0.0
        lambe=0.0
        mue=0.0
        lambmue=0.0
		taue=0.0
		taupstine=0.0
		tausstine=0.0
		taustsse=0.0
		rtaustsse=0.0
		taupstinsse=0.0
		tausstinsse=0.0
        topoe=0.0
        itopoe=0.0
        record1=0.0
        record2=0.0
		spgxn=0.0
		spgxp=0.0
		spgzn=0.0
		spgzp=0.0
		spg1xn=0.0
		spg2xp=0.0
		spg1zn=0.0
		spg2zp=0.0

		!Get currshot velocity&anisotropic parameters
        call get_currshot_parameters(vp,vs,rho,qfp,qfs,topo,&
			 vpe,vse,rhoe,qfpe,qfse,topoe,itopoe,nx,nz,nxe,nze,&
			 nxea,nzea,ixmine,ixmaxe,npxn,npxp,npzn,npzp,&
			 isx,isz,isxe,isze,dx,dz)

!		rhoe=1.0
		qfpe=100000.0
		qfse=100000.0


		!Build staggered grids
 		call substagger(nxea,nzea,rhoe,busge)

		!Pre-calculate parameters in single shot
		call pcalpas(f0,vpe,vse,rhoe,qfpe,qfse,lambe,mue,lambmue,&
			taue,taupstine,tausstine,taustsse,rtaustsse,taupstinsse,&
			tausstinsse,nxea,nzea,dx,dz)

		!Print some basic information on the screen
		write(*,*)'nxe,nxea,nze,nzea,isz,isx,isze,isxe'
		write(*,*)nxe,nxea,nze,nzea,isz,isx,isze,isxe

		!Building absorbing coefficients
		call csponge(spgxn,spgxp,spgzn,spgzp,spg1xn,spg2xp,&
				spg1zn,spg2zp,npxn,npxp,npzn,npzp,dx,dz)


        do it=1,nt+decay
			!Modified from here
   	         call extrapolation_one_step_10nd_absbc(&
			 		it,isxe,isze,nxe,nze,nxea,nzea,npxn,npxp,npzn,npzp,&
					f0,dx,dz,dt,vx,vz,strxx,strzz,strxz,rxx,rzz,rxz,&
					vxx,vxz,vzx,vzz,strxxx,strxxz,strxxr,strzzx,strzzz,strzzr,&
					strxzx,strxzz,strxzr,rxxx,rxxz,rxxr,rzzx,rzzz,rzzr,&
					rxzx,rxzz,rxzr,der1,der2,der3,der4,spgxn,spgxp,spgzn,spgzp,&
					spg1xn,spg2xp,spg1zn,spg2zp,vpe,vse,rhoe,qfpe,qfse,&
					lambe,mue,lambmue,taue,taupstine,tausstine,taustsse,&
					rtaustsse,taupstinsse,tausstinsse,itopoe,busge)


			iit=it-decay
			if(iit>=1)then
				!This should be modified
				do i=1,nxe
					record1(iit,i)=strzz(itopoe(i+npxn)+nrz0,i+npxn) 
					record2(iit,i)=strxx(itopoe(i+npxn)+nrz0,i+npxn) 
				enddo
			endif


	 		if (myid==0)then
	          	if(modulo(it,50)==0)then
	              write(*,*)'shot is',trim(csname),'  myid=',myid,'nt=',nt,'it=',it
	
					write(currt,'(I6)')it

					write(currtsnap,*)trim(adjustl(snapx_fn)),trim(adjustl(currt))
					open(unit=8,file=currtsnap,form='unformatted',access='direct',&
					status='replace',recl=nze*nxe)
					write(8,rec=1)((strxx(k,i),k=npzn+1,nze+npzn),i=1+npxn,nxe+npxn)
					close(unit=8)

					write(currtsnap,*)trim(adjustl(snapz_fn)),trim(adjustl(currt))
					open(unit=8,file=currtsnap,form='unformatted',access='direct',&
					status='replace',recl=nze*nxe)
					write(8,rec=1)((strzz(k,i),k=npzn+1,nze+npzn),i=1+npxn,nxe+npxn)
					close(unit=8)

	          	endif
	 		endif


	enddo
    
    !*************************************************************************************  
      
        write(*,*)'shot  ',trim(adjustl(csname)),&
        '  extrapolation finished,and now is writing disk'

        call write_currshot_disk(record1,shot_fn1,record2,shot_fn3,&
								csname,csno,&!**zy
                               	isx,isz,nxe,ixmine,ixmaxe,nt,dx,dz,dt)


        write(*,*)'shot  ',trim(adjustl(csname)),'  is done'

        deallocate(vx,STAT=err)
        deallocate(vz,STAT=err)
        deallocate(strxx,STAT=err)
        deallocate(strzz,STAT=err)
        deallocate(strxz,STAT=err)
        deallocate(rxx,STAT=err)
        deallocate(rzz,STAT=err)
        deallocate(rxz,STAT=err)
        deallocate(vxx,STAT=err)
        deallocate(vxz,STAT=err)
        deallocate(vzx,STAT=err)
        deallocate(vzz,STAT=err)
        deallocate(strxxx,STAT=err)
        deallocate(strxxz,STAT=err)
        deallocate(strxxr,STAT=err)
        deallocate(strzzx,STAT=err)
        deallocate(strzzz,STAT=err)
        deallocate(strzzr,STAT=err)
        deallocate(strxzx,STAT=err)
        deallocate(strxzz,STAT=err)
        deallocate(strxzr,STAT=err)
        deallocate(rxxx,STAT=err)
        deallocate(rxxz,STAT=err)
        deallocate(rxxr,STAT=err)
        deallocate(rzzx,STAT=err)
        deallocate(rzzz,STAT=err)
        deallocate(rzzr,STAT=err)
        deallocate(rxzx,STAT=err)
        deallocate(rxzz,STAT=err)
        deallocate(rxzr,STAT=err)
        deallocate(der1,STAT=err)
        deallocate(der2,STAT=err)
        deallocate(der3,STAT=err)
        deallocate(der4,STAT=err)
        deallocate(vpe,STAT=err)
        deallocate(vse,STAT=err)
        deallocate(rhoe,STAT=err)
        deallocate(qfpe,STAT=err)
        deallocate(qfse,STAT=err)
		deallocate(busge,STAT=err)
        deallocate(lambe,STAT=err)
        deallocate(mue,STAT=err)
        deallocate(lambmue,STAT=err)
		deallocate(taue,STAT=err)
		deallocate(taupstine,STAT=err)
		deallocate(tausstine,STAT=err)
		deallocate(taustsse,STAT=err)
		deallocate(rtaustsse,STAT=err)
		deallocate(taupstinsse,STAT=err)
		deallocate(tausstinsse,STAT=err)
        deallocate(topoe,STAT=err)
        deallocate(itopoe,STAT=err)
        deallocate(record1,STAT=err)
        deallocate(record2,STAT=err)
		deallocate(spgxn,STAT=err)
		deallocate(spgxp,STAT=err)
		deallocate(spgzn,STAT=err)
		deallocate(spgzp,STAT=err)
		deallocate(spg1xn,STAT=err)
		deallocate(spg2xp,STAT=err)
		deallocate(spg1zn,STAT=err)
		deallocate(spg2zp,STAT=err)

        return
        end subroutine single_shot_modelling




         






        






























        

         


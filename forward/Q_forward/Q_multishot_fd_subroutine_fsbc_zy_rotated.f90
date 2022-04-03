!*******constant module and externel library*************
        module constant
        implicit none 
        real,parameter::pi=3.1415926
        real,parameter::R=1.0E-3

		
		real,parameter::a20=1.0

		real,parameter::a40=9.0/8.0
		real,parameter::a41=-1.0/24.0


		real,parameter::a60=75.0/64.0
		real,parameter::a61=-25.0/384.0
		real,parameter::a62=3.0/640.0

		real,parameter::a80=1225.0/1024.0
		real,parameter::a81=-245.0/3072.0
		real,parameter::a82=49.0/5120.0
		real,parameter::a83=-5.0/7168.0


		real,parameter::a100=19845.0/16384.0
		real,parameter::a101=-735.0/8192.0
		real,parameter::a102=567.0/40960.0
		real,parameter::a103=-405.0/229376.0
		real,parameter::a104=35.0/294912.0

        end module constant


        module header_module
        type :: segy
         integer(kind = 4) :: tracl
         integer(kind = 4) :: tracr
         integer(kind = 4) :: fldr
         integer(kind = 4) :: tracf
         integer(kind = 4) :: ep
         integer(kind = 4) :: cdp
         integer(kind = 4) :: cdpt
         integer(kind = 2) :: trid
         integer(kind = 2) :: nvs
         integer(kind = 2) :: nhs
         integer(kind = 2) :: duse
         integer(kind = 4) :: offset
         integer(kind = 4) :: gelev
         integer(kind = 4) :: selev
         integer(kind = 4) :: sdepth
         integer(kind = 4) :: gdel
         integer(kind = 4) :: sedl
         integer(kind = 4) :: swdep
         integer(kind = 4) :: gwdep
         integer(kind = 2) :: scalel
         integer(kind = 2) :: scalco
         integer(kind = 4) :: sx
         integer(kind = 4) :: sy
         integer(kind = 4) :: gx
         integer(kind = 4) :: gy
         integer(kind = 2) :: counit
         integer(kind = 2) :: wevel
         integer(kind = 2) :: swevel
         integer(kind = 2) :: sut
         integer(kind = 2) :: gut
         integer(kind = 2) :: sstat
         integer(kind = 2) :: gstat
         integer(kind = 2) :: tstat
         integer(kind = 2) :: laga
         integer(kind = 2) :: lagb
         integer(kind = 2) :: delrt
         integer(kind = 2) :: muts
         integer(kind = 2) :: mute
         integer(kind = 2) :: ns
         integer(kind = 2) :: dt
         integer(kind = 2) :: gain
         integer(kind = 2) :: igc
         integer(kind = 2) :: igi
         integer(kind = 2) :: corr
         integer(kind = 2) :: sfs
         integer(kind = 2) :: sfe
         integer(kind = 2) :: slen
         integer(kind = 2) :: styp
         integer(kind = 2) :: stas
         integer(kind = 2) :: stae
         integer(kind = 2) :: tatyp
         integer(kind = 2) :: afilf
         integer(kind = 2) :: afils
         integer(kind = 2) :: nofilf
         integer(kind = 2) :: nofils
         integer(kind = 2) :: lcf
         integer(kind = 2) :: hcf
         integer(kind = 2) :: lcs
         integer(kind = 2) :: hcs
         integer(kind = 2) :: year
         integer(kind = 2) :: day
         integer(kind = 2) :: hour
         integer(kind = 2) :: minute
         integer(kind = 2) :: sec
         integer(kind = 2) :: timbas
         integer(kind = 2) :: trwf
         integer(kind = 2) :: grnors
         integer(kind = 2) :: grnofr
         integer(kind = 2) :: grnlof
         integer(kind = 2) :: gaps
         integer(kind = 2) :: otrav
         real(kind = 4) :: d1
         real(kind = 4) :: f1
         real(kind = 4) :: d2
         real(kind = 4) :: f2
         real(kind = 4) :: ungpow
         real(kind = 4) :: unscale
         integer(kind = 4) :: ntr
         integer(kind = 2) :: mark
         integer(kind = 2) :: shortpad
         integer(kind = 2) :: unass(14)
         end type
        end module


!***********************************************************************************
        subroutine read_par(par_fn,vp_fn,vs_fn,rho_fn,qfp_fn,qfs_fn,&
					shot_fn1,shot_fn2,shot_fn3,shot_fn4,&
					nshot,offset_min,offset_max,sx0,sz0,rz0,&
					dsx,nx,nz,dx,dz,nt,dt,f0,ppxn,ppxp,ppzn,ppzp,&
					fdo2,fdo1)

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
        integer::fdo2
        integer::fdo1
        integer::nshot
        real::offset_min,offset_max
        real::dx,dz
!        real::pml_thick!(m)
        real::ppxn,ppxp,ppzn,ppzp!(m)
        real::sx0,sz0,rz0!(m)
        real::dsx!(m)
        real::f0!(Hz)
        real::dt!(s)
        !*variables of MPI interface*
        integer::myid
      
        !Local variables
        integer::i,j,ierr
        character(len=256)::par_jumpper


        open(10,file=par_fn,action='read',status='old',form='formatted',&
        access='sequential')
        
        read(10,'(A)')par_jumpper
        read(10,'(A)')vp_fn

        read(10,'(A)')par_jumpper
        read(10,'(A)')vs_fn

        read(10,'(A)')par_jumpper
        read(10,'(A)')rho_fn

        read(10,'(A)')par_jumpper
        read(10,'(A)')qfp_fn

        read(10,'(A)')par_jumpper
        read(10,'(A)')qfs_fn

        read(10,'(A)')par_jumpper
        read(10,'(A)')shot_fn1

        read(10,'(A)')par_jumpper
        read(10,'(A)')shot_fn2!**zy

        read(10,'(A)')par_jumpper
        read(10,'(A)')shot_fn3!**zy

        read(10,'(A)')par_jumpper
        read(10,'(A)')shot_fn4!**zy

        read(10,'(A)')par_jumpper
        read(10,*)nshot

        read(10,'(A)')par_jumpper
        read(10,*)offset_min,offset_max

        read(10,'(A)')par_jumpper
        read(10,*)sx0,sz0,rz0

        read(10,'(A)')par_jumpper
        read(10,'(F)')dsx

        read(10,'(A)')par_jumpper
        read(10,*)nx,nz

        read(10,'(A)')par_jumpper
        read(10,*)dx,dz

        read(10,'(A)')par_jumpper
        read(10,'(I)')nt

        read(10,'(A)')par_jumpper
        read(10,'(F)')dt

        read(10,'(A)')par_jumpper
        read(10,*)f0

        read(10,'(A)')par_jumpper
        read(10,*)ppxn,ppxp

        read(10,'(A)')par_jumpper
        read(10,*)ppzn,ppzp!**zy

        read(10,'(A)')par_jumpper
        read(10,'(I)')fdo2

        read(10,'(A)')par_jumpper
        read(10,'(I)')fdo1

        close(10)

		return
        end subroutine read_par


!********subroutine*******************************************
        subroutine read_data(vp_fn,vs_fn,rho_fn,qfp_fn,qfs_fn,&
					vp,vs,rho,qfp,qfs,nx,nz,myid)

        use constant
        implicit none
        !Dummy variables
        character(len=256)::vp_fn,vs_fn,rho_fn,&
							qfp_fn,qfs_fn,topo_fn

        integer::nx,nz,myid
        real::vp(nz,nx),vs(nz,nx),rho(nz,nx),&
			qfp(nz,nx),qfs(nz,nx)


        !Local variables
        integer::i,k

        !Read data from files

        open(10,file=vp_fn,action='read',&
        form='unformatted',access='direct',status='old',recl=nz*nx)
		read(10,rec=1)((vp(k,i),k=1,nz),i=1,nx)
        close (10)

        open(10,file=vs_fn,action='read',&
        form='unformatted',access='direct',status='old',recl=nz*nx)
		read(10,rec=1)((vs(k,i),k=1,nz),i=1,nx)
        close (10)

        open(10,file=rho_fn,action='read',&
        form='unformatted',access='direct',status='old',recl=nz*nx)
		read(10,rec=1)((rho(k,i),k=1,nz),i=1,nx)
        close (10)

        open(10,file=qfp_fn,action='read',&
        form='unformatted',access='direct',status='old',recl=nz*nx)
		read(10,rec=1)((qfp(k,i),k=1,nz),i=1,nx)
        close (10)


        open(10,file=qfs_fn,action='read',&
        form='unformatted',access='direct',status='old',recl=nz*nx)
		read(10,rec=1)((qfs(k,i),k=1,nz),i=1,nx)
        close (10)

        return
        end subroutine read_data


!====================================================================================
        subroutine get_currshot_parameters(vp,vs,rho,qfp,qfs,topo,&
			 		vpe,vse,rhoe,qfpe,qfse,topoe,itopoe,nx,nz,nxe,nze,&
					nxea,nzea,ixmine,ixmaxe,npxn,npxp,npzn,npzp,&
					isx,isz,isxe,isze,dx,dz)

        use constant

        implicit none
        !Dummy variables
        integer,intent(in out)::isxe,isze,itopoe(nxea)
        integer,intent(in)::nx,nz,nxe,nze,nxea,nzea,&
							npxn,npxp,npzn,npzp,ixmine,&
							ixmaxe,isx,isz

		real,intent(in)::dx,dz
        real,intent(in)::vp(nz,nx),vs(nz,nx),rho(nz,nx),&
						qfp(nz,nx),qfs(nz,nx),topo(nx)


        real,intent(in out)::vpe(nzea,nxea),vse(nzea,nxea),&
							rhoe(nzea,nxea),qfpe(nzea,nxea),&
							qfse(nzea,nxea),topoe(nxea)

        !Local  variables
        integer::i,k,ka,i_ori

        !Assigning&Padding the  parameters
        do i=1,nxe
          i_ori=i+ixmine
          if(i_ori<1)i_ori=1
          if(i_ori>nx)i_ori=nx

           vpe(npzn+1:nze+npzn,i+npxn)=vp(1:nze,i_ori)
           vse(npzn+1:nze+npzn,i+npxn)=vs(1:nze,i_ori)
           rhoe(npzn+1:nze+npzn,i+npxn)=rho(1:nze,i_ori)
           qfpe(npzn+1:nze+npzn,i+npxn)=qfp(1:nze,i_ori)
           qfse(npzn+1:nze+npzn,i+npxn)=qfs(1:nze,i_ori)
           topoe(i+npxn)=topo(i_ori)
        enddo


		!Using a simpler way to pa velocity
		do k=1,npzn
			do i=1,nxea
				vpe(k,i)=vpe(npzn+1,i)
				vse(k,i)=vse(npzn+1,i)
				rhoe(k,i)=rhoe(npzn+1,i)
				qfpe(k,i)=qfpe(npzn+1,i)
				qfse(k,i)=qfse(npzn+1,i)
			enddo
		enddo


		do k=npzn+nze+1,nzea
			do i=1,nxea
				vpe(k,i)=vpe(npzn+nz,i)
				vse(k,i)=vse(npzn+nz,i)
				rhoe(k,i)=rhoe(npzn+nz,i)
				qfpe(k,i)=qfpe(npzn+nz,i)
				qfse(k,i)=qfse(npzn+nz,i)
			enddo
		enddo


		do k=1,nzea
			do i=1,npxn
				vpe(k,i)=vpe(k,npxn+1)
				vse(k,i)=vse(k,npxn+1)
				rhoe(k,i)=rhoe(k,npxn+1)
				qfpe(k,i)=qfpe(k,npxn+1)
				qfse(k,i)=qfse(k,npxn+1)
				topoe(i)=topoe(npxn+1)
			enddo
		enddo


		do k=1,nzea
			do i=npxn+nxe+1,nxea
				vpe(k,i)=vpe(k,npxn+nxe)
				vse(k,i)=vse(k,npxn+nxe)
				rhoe(k,i)=rhoe(k,npxn+nxe)
				qfpe(k,i)=qfpe(k,npxn+nxe)
				qfse(k,i)=qfse(k,npxn+nxe)
				topoe(i)=topoe(npxn+nxe)
			enddo
		enddo


		!Pay attenuation to the defination of itopoe
		do i=1,nxea
			itopoe(i)=nint(topoe(i)/dz)+npzn+1
		enddo


        !Compute currshot source_position
        isxe=isx-ixmine+npxn+1!**zy
        isze=isz+itopoe(isxe)


        open(10,file='vpe_test',action='write',&
        form='unformatted',access='direct',status='replace',recl=nxe*nze)
		 write(10,rec=1)((vpe(k,i),k=npzn+1,npzn+nze),i=npxn+1,npxn+nxe)
        close (10)

        open(10,file='vse_test',action='write',&
        form='unformatted',access='direct',status='replace',recl=nxe*nze)
		 write(10,rec=1)((vse(k,i),k=npzn+1,npzn+nze),i=npxn+1,npxn+nxe)
        close (10)


        open(10,file='rhoe_test',action='write',&
        form='unformatted',access='direct',status='replace',recl=nxe*nze)
		 write(10,rec=1)((rhoe(k,i),k=npzn+1,npzn+nze),i=npxn+1,npxn+nxe)
        close (10)

!        open(10,file='epie_test',action='write',&
!        form='unformatted',access='direct',status='replace',recl=nxe*nze)
!		 write(10,rec=1)((epie(k,i),k=npzn+1,npzn+nze),i=npxn+1,npxn+nxe)
!        close (10)


        open(10,file='qfpe_test',action='write',&
        form='unformatted',access='direct',status='replace',recl=nxe*nze)
		 write(10,rec=1)((qfpe(k,i),k=npzn+1,npzn+nze),i=npxn+1,npxn+nxe)
        close (10)


        open(10,file='qfse_test',action='write',&
        form='unformatted',access='direct',status='replace',recl=nxe*nze)
		 write(10,rec=1)((qfse(k,i),k=npzn+1,npzn+nze),i=npxn+1,npxn+nxe)
        close (10)


        open(10,file='topoe_test',action='write',&
        form='unformatted',access='direct',status='replace',recl=nxe)
		 write(10,rec=1)(topoe(i),i=npxn+1,npxn+nxe)
        close (10)


        return
        end subroutine get_currshot_parameters


!===========================================================================
!**************Modfied by zy to adjust rotated staggered grid***************

	 subroutine substagger(nxea,nzea,rhoe,busg1)
	 implicit none
	 !Dummy variables
	 integer::nxea,nzea
	 real::rhoe(nzea,nxea),busg1(nzea,nxea)
	 !Local variables
	 integer::i,k

	! #########################################################################
	! Busg1
	 do i=1,nxea-1
	 	do k=1,nzea-1
	 		rhoe(k,i)=0.25*(rhoe(k,i)+rhoe(k+1,i)+rhoe(k+1,i+1)+rhoe(k,i+1))
	 	end do	
	 end do


	 do i=1,nxea
	 	rhoe(nzea,i)=rhoe(nzea-1,i)
	 end do
	
	 do k=1,nzea
	 	rhoe(k,nxea)=rhoe(k,nxea-1)
	 end do	

	
	 ! convert to 1/rho
	 do i=1,nxea
	 	do k=1,nzea
	 		if (rhoe(k,i).eq.0.) stop 'FIND A DENSITY EQUAL TO 0'
	 			busg1(k,i)=1./rhoe(k,i)
	 	end do
	 end do

	 return
	 end subroutine substagger


!****************************************************************************
	subroutine pcalpas(f0,vpe,vse,rhoe,qfpe,qfse,lambe,mue,lambmue,&
			taue,taupstine,tausstine,taustsse,rtaustsse,taupstinsse,&
			tausstinsse,nxea,nzea,dx,dz)

	use constant
	implicit none
	!Dummy variables
	integer::nxea,nzea
	real::f0,dx,dz,vpe(nzea,nxea),vse(nzea,nxea),&
		rhoe(nzea,nxea),qfpe(nzea,nxea),&
		qfse(nzea,nxea),lambe(nzea,nxea),&
		mue(nzea,nxea),lambmue(nzea,nxea),&
		taue(nzea,nxea),taupstine(nzea,nxea),&
		tausstine(nzea,nxea),taustsse(nzea,nxea),&
		rtaustsse(nzea,nxea),taupstinsse(nzea,nxea),&
		tausstinsse(nzea,nxea)
	!Local variables
	integer::i,k
	real::tau0,taup,taus,tp,fup,flow,nup,nlow,dup,dlow


!$OMP PARALLEL PRIVATE(taup,taus,nlow,nup,dlow,dup)
!$OMP DO 
	 do i=1,nxea
	 	do k=1,nzea

			!Pay attenuation to the expression of mue
			mue(k,i)=rhoe(k,i)*vse(k,i)**2
			lambe(k,i)=rhoe(k,i)*vpe(k,i)**2-2.0*mue(k,i)
			lambmue(k,i)=rhoe(k,i)*vpe(k,i)**2

			!Using Carcione's approach
!			tau0=1.0/f0	
!			taupstine(k,i)=tau0/qfpe(k,i)*(sqrt(qfpe(k,i)**2+1)+1)
!			tausstine(k,i)=tau0/qfse(k,i)*(sqrt(qfse(k,i)**2+1)+1)
!			taustsse(k,i)=tau0/qfpe(k,i)*(sqrt(qfpe(k,i)**2+1)-1)


			!Using Blanch's approach
			flow=2*2*pi
			fup=2.5*f0*2*pi

			qfpe(k,i)=qfpe(k,i)*(18.0/20.0)
			qfse(k,i)=qfse(k,i)*(18.0/20.0)

			taustsse(k,i)=1.0/((fup-flow)/2.0)

			nlow=log(taustsse(k,i)**2*flow**2+1)
			nup=log(taustsse(k,i)**2*fup**2+1)

			dlow=atan(flow*taustsse(k,i))-(flow*taustsse(k,i))/(1+taustsse(k,i)**2*flow**2)
			dup=atan(fup*taustsse(k,i))-(fup*taustsse(k,i))/(1+taustsse(k,i)**2*fup**2)


			if(dup-dlow/=0.0)then
				taup=1.0/qfpe(k,i)*(nup-nlow)/(dup-dlow)
				taus=1.0/qfse(k,i)*(nup-nlow)/(dup-dlow)
			endif

			taupstine(k,i)=taup*taustsse(k,i)+taustsse(k,i)
			tausstine(k,i)=taus*taustsse(k,i)+taustsse(k,i)


!			write(*,*)f0,1/f0,taustsse(1,1)*1.0E3,taupstine(1,1)*1.0E3
!			pause

			!Other variables used in this program
			rtaustsse(k,i)=1.0/taustsse(k,i)
			taupstinsse(k,i)=taupstine(k,i)/taustsse(k,i)
			tausstinsse(k,i)=tausstine(k,i)/taustsse(k,i)

!			write(*,*) taupstinsse(k,i), tausstinsse(k,i)

		enddo
	enddo
!$OMP END DO
!$OMP END PARALLEL
	return
	end subroutine pcalpas


!******************************************************************************
	subroutine csponge(spgxn,spgxp,spgzn,spgzp,spg1xn,spg2xp,spg1zn,spg2zp,&
				npxn,npxp,npzn,npzp,dx,dz)

	use constant
	implicit none
	!Dummy variables
	integer::npxn,npxp,npzn,npzp
	real::dx,dz,&
		spgxn(npxn+1),spgxp(npxp+1),&
		spgzn(npzn+1),spgzp(npzp+1),&
		spg1xn(npxn+1),spg2xp(npxp+1),&
		spg1zn(npzn+1),spg2zp(npzp+1)


	!Local variables
	integer::i,j,k
	real::facxn,facxp,faczn,faczp,dxro,dzro


	facxn=acos(R)/(float(npxn)*dx)
	facxp=acos(R)/(float(npxp)*dx)
	faczn=acos(R)/(float(npzn)*dz)
	faczp=acos(R)/(float(npzp)*dz)


	do i=1,npxn+1
		spgxn(i)=cos(facxn*float(npxn-i+1)*dx)
		spg1xn(i)=cos(facxn*(float(npxn-i+1)*dx-dx/2))
		if (i.eq.npxn+1) spg1xn(i)=1.
	end do


	do i=1,npzn+1
		spgzn(i)=cos(faczn*float(npzn-i+1)*dz)
		spg1zn(i)=cos(faczn*(float(npzn-i+1)*dz-dz/2))
		if(i.eq.npzn+1)spg1zn(i)=1.
	enddo


	do i=1,npxp+1
		spgxp(i)=cos(facxp*float(npxp-i+1)*dx)
		spg2xp(i)=cos(facxp*(float(npxp-i+1)*dx+dx/2))
	end do

	do i=1,npzp+1
		spgzp(i)=cos(faczp*float(npzp-i+1)*dz)
		spg2zp(i)=cos(faczp*(float(npzp-i+1)*dz+dz/2))
	end do

	return
	end 


!****************************************************************************
	subroutine calpmlbclrtb(npxn,npxp,npzn,npzp,nx,nz,dx,dz,&
							itopo,spgxn,spgxp,spgzn,spgzp,ux,uz)

	implicit none
	!Dummy variables
	integer::nx,nz,npxn,npxp,npzn,npzp,itopo(nx)
	real::dx,dz,spgxn(npxn+1),spgxp(npxp+1),&
			spgzn(npzn+1),spgzp(npzp+1),&
        	ux(0:nz,nx),uz(0:nz,nx)

	!Local variables
	integer::i,k

	!Left
	do i=1,npxn+1
		do k=itopo(i)-npzn,nz
			ux(k,i)=ux(k,i)*spgxn(i)
		end do
	end do


	!Right
	do i=1,npxp+1
		do k=itopo(i)-npzn,nz
			ux(k,nx-i+1)=ux(k,nx-i+1)*spgxp(i)
		end do
	end do


	!Top
	do i=1,nx
		do k=itopo(i)-npzn,itopo(i)
			uz(k,i)=uz(k,i)*spgzn(k-(itopo(i)-npzn)+1)
		end do
	end do


	!Bottom
	do i=1,nx
        do k=1,npzp+1
        	uz(nz-k+1,i)=uz(nz-k+1,i)*spgzp(k)
        end do
	end do

	return
	end subroutine calpmlbclrtb


!*****************************************************************************
	subroutine calpmlbclrb(npxn,npxp,npzn,npzp,nxe,nze,nx,nz,dx,dz,itopo,&
				spgxn,spgxp,spgzp,ux,uz)
	implicit none
	!Dummy variables
	integer::nxe,nze,nx,nz,npxn,npxp,npzn,npzp,itopo(nx)
	real::dx,dz,spgxn(npxn+1),spgxp(npxp+1),spgzp(npzp+1),&
        	ux(0:nz,nx),uz(0:nz,nx)

	!Local variables
	integer::i,k

	!Left
	do i=1,npxn+1
		do k=itopo(i),nz
			ux(k,i)=ux(k,i)*spgxn(i)
		end do
	end do

	!Right
	do i=1,npxp+1
		do k=itopo(i),nz
			ux(k,nx-i+1)=ux(k,nx-i+1)*spgxp(i)
		end do
	end do

	!Bottom
	do i=1,nx
       do k=1,npzp+1
        	uz(nz-k+1,i)=uz(nz-k+1,i)*spgzp(k)
        end do
	end do

	return
	end subroutine calpmlbclrb


!****************************************************************************
	subroutine extrapolation_one_step_10nd_absbc(&
 			it,isxe,isze,nxe,nze,nx,nz,npxn,npxp,npzn,npzp,&
			f0,dx,dz,dt,vx,vz,strxx,strzz,strxz,rxx,rzz,rxz,&
			vxx,vxz,vzx,vzz,strxxx,strxxz,strxxr,strzzx,strzzz,strzzr,&
			strxzx,strxzz,strxzr,rxxx,rxxz,rxxr,rzzx,rzzz,rzzr,&
			rxzx,rxzz,rxzr,der1,der2,der3,der4,spgxn,spgxp,spgzn,spgzp,&
			spg1xn,spg2xp,spg1zn,spg2zp,vp,vs,rho,qfp,qfs,&
			lamb,mu,lambmu,tau,taupstin,tausstin,taustss,&
			rtaustss,taupstinss,tausstinss,itopo,busg)

	use constant
	implicit none
	!Dummy variables
	integer::it,isxe,isze,nxe,nze,nx,nz,npxn,npxp,npzn,npzp
    integer::itopo(nx)
	real::f0,dx,dz,dt
	real::vx(-4:nz+5,-4:nx+5),vz(-4:nz+5,-4:nx+5),&
        strxx(-4:nz+5,-4:nx+5),strzz(-4:nz+5,-4:nx+5),&
        strxz(-4:nz+5,-4:nx+5),&
        rxx(nz,nx),rxz(nz,nx),rzz(nz,nx),&
        vxx(0:nz,nx),vxz(0:nz,nx),&
        vzx(0:nz,nx),vzz(0:nz,nx),&
        strxxx(0:nz,nx),strxxz(0:nz,nx),&
		strxxr(0:nz,nx),&
        strzzx(0:nz,nx),strzzz(0:nz,nx),&
		strzzr(0:nz,nx),&
        strxzx(0:nz,nx),strxzz(0:nz,nx),&
		strxzr(0:nz,nx),&
        rxxx(0:nz,nx),rxxz(0:nz,nx),&
		rxxr(0:nz,nx),&
        rzzx(0:nz,nx),rzzz(0:nz,nx),&
		rzzr(0:nz,nx),&
        rxzx(0:nz,nx),rxzz(0:nz,nx),&
		rxzr(0:nz,nx),&
        der1(nz,nx),der2(nz,nx),&
        der3(nz,nx),der4(nz,nx),&
        spgxn(npxn+1),&
		spgxp(npxp+1),spgzn(npzn+1),&
		spgzp(npzp+1),spg1xn(npxn+1),&
		spg2xp(npxp+1),spg1zn(npzn+1),&
		spg2zp(npzp+1),&
		busg(nz,nx),vp(nz,nx),&
		vs(nz,nx),rho(nz,nx),&
		qfp(nz,nx),qfs(nz,nx),&
		lamb(nz,nx),mu(nz,nx),&
		lambmu(nz,nx),tau(nz,nx),&
		taupstin(nz,nx),tausstin(nz,nx),&
		taustss(nz,nx),rtaustss(nz,nx),&
		taupstinss(nz,nx),tausstinss(nz,nx)

		!Local variables
		integer::i,k,idir,ii,oop,xflag,zflag
		real::dts,tsmax,wavelet



!******************Updating velocity component******************
		!Updating vx
		xflag=0
		call calderxabs(xflag,nx,nz,npzn,dx,dt,strxx,der1,itopo)
		call calpvabs(nx,nz,npzn,dt,der1,vxx,busg,itopo)
		zflag=0
		call calderzabs(zflag,nx,nz,npzn,dz,dt,strxz,der2,itopo)
		call calpvabs(nx,nz,npzn,dt,der2,vxz,busg,itopo)

		call calpmlbclrtb(npxn,npxp,npzn,npzp,nx,nz,dx,dz,&
					itopo,spgxn,spgxp,spgzn,spgzp,vxx,vxz)
		call calsumvel(nx,nz,vxx,vxz,vx)


		!Updating vz
		xflag=0
		call calderxabs(xflag,nx,nz,npzn,dx,dt,strxz,der1,itopo)
		call calpvabs(nx,nz,npzn,dt,der1,vzx,busg,itopo)
		zflag=0
		call calderzabs(zflag,nx,nz,npzn,dz,dt,strzz,der2,itopo)
		call calpvabs(nx,nz,npzn,dt,der2,vzz,busg,itopo)

		call calpmlbclrtb(npxn,npxp,npzn,npzp,nx,nz,dx,dz,&
				itopo,spgxn,spgxp,spgzn,spgzp,vzx,vzz)
		call calsumvel(nx,nz,vzx,vzz,vz)


!******************Updating stress component******************
		!der(odd) implies  x direciton
		!der(even) implies  z direciton
		xflag=1
		zflag=1
		call calderxabs(xflag,nx,nz,npzn,dx,dt,vx,der1,itopo)
		call calderxabs(xflag,nx,nz,npzn,dx,dt,vz,der2,itopo)
		call calderzabs(zflag,nx,nz,npzn,dz,dt,vx,der3,itopo)
		call calderzabs(zflag,nx,nz,npzn,dz,dt,vz,der4,itopo)


		!Updating strxx
		call calpstrxx(nx,nz,npzn,dt,der1,der4,rxx,strxxx,&
				strxxz,strxxr,mu,lambmu,taupstinss,tausstinss,itopo)
		call calpmlbclrtb(npxn,npxp,npzn,npzp,nx,nz,dx,dz,&
				itopo,spg1xn,spg2xp,spg1zn,spg2zp,strxxx,strxxz)
		call addsrc(nx,nz,isxe,isze,it,f0,dx,dz,dt,strxxx)
		call calsumstr(nx,nz,strxxx,strxxz,strxxr,strxx)

		!Updating strzz
		call calpstrzz(nx,nz,npzn,dt,der1,der4,rzz,strzzx,&
				strzzz,strzzr,mu,lambmu,taupstinss,tausstinss,itopo)
		call calpmlbclrtb(npxn,npxp,npzn,npzp,nx,nz,dx,dz,&
				itopo,spg1xn,spg2xp,spg1zn,spg2zp,strzzx,strzzz)
		call addsrc(nx,nz,isxe,isze,it,f0,dx,dz,dt,strzzx)
		call calsumstr(nx,nz,strzzx,strzzz,strzzr,strzz)


		!Updating strxz
		call calpstrxz(nx,nz,npzn,dt,der2,der3,rxz,strxzx,&
				strxzz,strxzr,mu,lambmu,taupstinss,tausstinss,itopo)
		call calpmlbclrtb(npxn,npxp,npzn,npzp,nx,nz,dx,dz,&
				itopo,spg1xn,spg2xp,spg1zn,spg2zp,strxzx,strxzz)
		call calsumstr(nx,nz,strxzx,strxzz,strxzr,strxz)


!****************Updating memory variables component***************
		!Updating rxx
		call calprxx(nx,nz,npzn,dt,der1,der4,rxx,rxxx,rxxz,rxxr,&
				mu,lambmu,rtaustss,taupstinss,tausstinss,itopo)
		call calpmlbclrtb(npxn,npxp,npzn,npzp,nx,nz,dx,dz,&
				itopo,spg1xn,spg2xp,spg1zn,spg2zp,rxxx,rxxz)
		call calsummr(nx,nz,rxxx,rxxz,rxxr,rxx)

		!Updating rzz
		call calprzz(nx,nz,npzn,dt,der1,der4,rzz,rzzx,rzzz,rzzr,&
				mu,lambmu,rtaustss,taupstinss,tausstinss,itopo)
		call calpmlbclrtb(npxn,npxp,npzn,npzp,nx,nz,dx,dz,&
				itopo,spg1xn,spg2xp,spg1zn,spg2zp,rzzx,rzzz)
		call calsummr(nx,nz,rzzx,rzzz,rzzr,rzz)


		!Updating rxz
		call calprxz(nx,nz,npzn,dt,der2,der3,rxz,rxzx,rxzz,rxzr,&
				mu,lambmu,rtaustss,taupstinss,tausstinss,itopo)
		call calpmlbclrtb(npxn,npxp,npzn,npzp,nx,nz,dx,dz,&
				itopo,spg1xn,spg2xp,spg1zn,spg2zp,rxzx,rxzz)
		call calsummr(nx,nz,rxzx,rxzz,rxzr,rxz)


	return
	end subroutine extrapolation_one_step_10nd_absbc


!*****************************************************************************
	subroutine calderxabs(xflag,nx,nz,npzn,dx,dt,u,der1,itopo)
	use constant
	implicit none
	!Dummy variables

	integer::xflag,nx,nz,npzn,itopo(nx)

	real::dx,dt,&
		u(-4:nz+5,-4:nx+5),&
    	der1(nz,nx)

	!Local variables
	integer::k,i

	if(xflag==0)then
!$OMP PARALLEL PRIVATE(k) 
!$OMP DO 
		do i=1,nx
			do k=itopo(i)-npzn,nz
				der1(k,i)=0.5/dx*&
				(&
					a100*(u(k,i)-u(k-1,i-1)+u(k-1,i)-u(k,i-1))+&
					a101*(u(k+1,i+1)-u(k-2,i-2)+u(k-2,i+1)-u(k+1,i-2))+&
					a102*(u(k+2,i+2)-u(k-3,i-3)+u(k-3,i+2)-u(k+2,i-3))+&
					a103*(u(k+3,i+3)-u(k-4,i-4)+u(k-4,i+3)-u(k+3,i-4))+&
					a104*(u(k+4,i+4)-u(k-5,i-5)+u(k-5,i+4)-u(k+4,i-5))&
				)
			enddo
		enddo
!$OMP END DO
!$OMP END PARALLEL


	else if(xflag==1)then
!$OMP PARALLEL PRIVATE(k) 
!$OMP DO 
		do i=1,nx
			do k=itopo(i)-npzn,nz
				der1(k,i)=0.5/dx*&
				(&

					a100*(u(k+1,i+1)-u(k,i)+u(k,i+1)-u(k+1,i))+&
					a101*(u(k+2,i+2)-u(k-1,i-1)+u(k-1,i+2)-u(k+2,i-1))+&
					a102*(u(k+3,i+3)-u(k-2,i-2)+u(k-2,i+3)-u(k+3,i-2))+&
					a103*(u(k+4,i+4)-u(k-3,i-3)+u(k-3,i+4)-u(k+4,i-3))+&
					a104*(u(k+5,i+5)-u(k-4,i-4)+u(k-4,i+5)-u(k+5,i-4))&
				)
			enddo
		enddo

!$OMP END DO
!$OMP END PARALLEL

	else
		write(*,*)'Wrong xflag value,please check'
		stop
	endif

	return
	end subroutine calderxabs

!**************************************************************************
	subroutine calderzabs(zflag,nx,nz,npzn,dz,dt,u,der2,itopo)
	use constant
	implicit none
	!Dummy variables
	integer::zflag,nx,nz,npzn,&
			itopo(nx)
	real::	dz,dt,&
			u(-4:nz+5,-4:nx+5),&
    		der2(nz,nx)

	!Local variables
	integer::k,i

	if(zflag==0)then


!$OMP PARALLEL PRIVATE(k) 
!$OMP DO 
		do i=1,nx
			do k=itopo(i)-npzn,nz
				der2(k,i)=0.5/dz*&
				(&
				a100*(u(k,i)-u(k-1,i-1)-u(k-1,i)+u(k,i-1))+&
				a101*(u(k+1,i+1)-u(k-2,i-2)-u(k-2,i+1)+u(k+1,i-2))+&
				a102*(u(k+2,i+2)-u(k-3,i-3)-u(k-3,i+2)+u(k+2,i-3))+&
				a103*(u(k+3,i+3)-u(k-4,i-4)-u(k-4,i+3)+u(k+3,i-4))+&
				a104*(u(k+4,i+4)-u(k-5,i-5)-u(k-5,i+4)+u(k+4,i-5))&
				)
			enddo
		enddo
!$OMP END DO
!$OMP END PARALLEL



	else if(zflag==1)then
!$OMP PARALLEL PRIVATE(k) 
!$OMP DO 
		do i=1,nx
			do k=itopo(i)-npzn,nz
				der2(k,i)=0.5/dz*&
				(&
				a100*(u(k+1,i+1)-u(k,i)-u(k,i+1)+u(k+1,i))+&
				a101*(u(k+2,i+2)-u(k-1,i-1)-u(k-1,i+2)+u(k+2,i-1))+&
				a102*(u(k+3,i+3)-u(k-2,i-2)-u(k-2,i+3)+u(k+3,i-2))+&
				a103*(u(k+4,i+4)-u(k-3,i-3)-u(k-3,i+4)+u(k+4,i-3))+&
				a104*(u(k+5,i+5)-u(k-4,i-4)-u(k-4,i+5)+u(k+5,i-4))&
				)
			enddo
		enddo

!$OMP END DO
!$OMP END PARALLEL
	else

		write(*,*)'Wrong zflag value,please check'
		stop

	endif

	return
	end subroutine calderzabs

!*******************************************************************
	subroutine calpvabs(nx,nz,npzn,dt,der1,u,pa,itopo)
	use constant
	implicit none
	!Dummy variables
	integer::nx,nz,npzn,itopo(nx)
	real::dt,&
    	u(0:nz,nx),pa(nz,nx),&
		der1(nz,nx)
	!Local variables
	integer::k,i

!$OMP PARALLEL PRIVATE(k) 
!$OMP DO 
		do i=1,nx
			do k=itopo(i)-npzn,nz
				u(k,i)=u(k,i)+dt*pa(k,i)*der1(k,i)
			enddo
		enddo
!$OMP END DO
!$OMP END PARALLEL
	return
	end subroutine calpvabs

!********************************************************************
	subroutine calpstrxx(nx,nz,npzn,dt,derx,derz,r,strx,strz,strr,&
				mu,lambmu,taupstinss,tausstinss,itopo)

	implicit none
	!Dummy variables
	integer::nx,nz,npzn,&
			itopo(nx)

	real::	dt,r(nz,nx),&
			strx(0:nz,nx),&
			strz(0:nz,nx),&
			strr(0:nz,nx),&
			derx(nz,nx),&
			derz(nz,nx),&
			taupstinss(nz,nx),&
			tausstinss(nz,nx),&
			mu(nz,nx),&
			lambmu(nz,nx)


	!Local variables
	integer::k,i

	!$OMP PARALLEL PRIVATE(k) 
	!$OMP DO 
		do i=1,nx
			do k=itopo(i)-npzn,nz

				strx(k,i)=strx(k,i)+dt*(lambmu(k,i)*taupstinss(k,i)*derx(k,i))

				strz(k,i)=strz(k,i)+dt*((lambmu(k,i)*taupstinss(k,i)-&
					2.0*mu(k,i)*tausstinss(k,i))*derz(k,i))
	
				strr(k,i)=strr(k,i)+dt*r(k,i)


!				strx(k,i)=strx(k,i)+dt*(lambmu(k,i)*derx(k,i))
!				strz(k,i)=strz(k,i)+dt*((lambmu(k,i)-2.0*mu(k,i))*derz(k,i))

			enddo
		enddo
	
	!$OMP END DO
	!$OMP END PARALLEL
	return
	end subroutine calpstrxx

!*********************************************************************
	subroutine calpstrzz(nx,nz,npzn,dt,derx,derz,r,strx,strz,strr,&
				mu,lambmu,taupstinss,tausstinss,itopo)

	implicit none
	!Dummy variables
	integer::nx,nz,npzn,&
			itopo(nx)

	real::	dt,r(nz,nx),&
			strx(0:nz,nx),&
			strz(0:nz,nx),&
			strr(0:nz,nx),&
			derx(nz,nx),&
			derz(nz,nx),&
			taupstinss(nz,nx),&
			tausstinss(nz,nx),&
			mu(nz,nx),&
			lambmu(nz,nx)


	!Local variables
	integer::k,i

	!$OMP PARALLEL PRIVATE(k) 
	!$OMP DO 
		do i=1,nx
			do k=itopo(i)-npzn,nz

				strx(k,i)=strx(k,i)+dt*((lambmu(k,i)*taupstinss(k,i)-&
					2.0*mu(k,i)*tausstinss(k,i))*derx(k,i))

				strz(k,i)=strz(k,i)+dt*(lambmu(k,i)*taupstinss(k,i)*derz(k,i))
				strr(k,i)=strr(k,i)+dt*r(k,i)


!				strx(k,i)=strx(k,i)+dt*((lambmu(k,i)-2.0*mu(k,i))*derx(k,i))
!				strz(k,i)=strz(k,i)+dt*(lambmu(k,i)*derz(k,i))

			enddo
		enddo
	
	!$OMP END DO
	!$OMP END PARALLEL
	return
	end subroutine calpstrzz

!*********************************************************************
	subroutine calpstrxz(nx,nz,npzn,dt,derx,derz,r,strx,strz,strr,&
				mu,lambmu,taupstinss,tausstinss,itopo)

	implicit none
	!Dummy variables
	integer::nx,nz,npzn,&
			itopo(nx)

	real::	dt,r(nz,nx),&
			strx(0:nz,nx),&
			strz(0:nz,nx),&
			strr(0:nz,nx),&
			derx(nz,nx),&
			derz(nz,nx),&
			taupstinss(nz,nx),&
			tausstinss(nz,nx),&
			mu(nz,nx),&
			lambmu(nz,nx)

	!Local variables
	integer::k,i

	!$OMP PARALLEL PRIVATE(k) 
	!$OMP DO 
		do i=1,nx
			do k=itopo(i)-npzn,nz

				strx(k,i)=strx(k,i)+dt*(mu(k,i)*tausstinss(k,i)*derx(k,i))

				strz(k,i)=strz(k,i)+dt*(mu(k,i)*tausstinss(k,i)*derz(k,i))

				strr(k,i)=strr(k,i)+dt*r(k,i)


!				strx(k,i)=strx(k,i)+dt*(mu(k,i)*derx(k,i))
!				strz(k,i)=strz(k,i)+dt*(mu(k,i)*derz(k,i))


			enddo
		enddo
	
	!$OMP END DO
	!$OMP END PARALLEL
	return
	end subroutine calpstrxz

!******************************************************************
	subroutine calprxx(nx,nz,npzn,dt,derx,derz,r,rx,rz,rr,&
				mu,lambmu,rtaustss,taupstinss,tausstinss,itopo)

	implicit none
	!Dummy variables
	integer::nx,nz,npzn,&
			itopo(nx)
	real::	dt,r(nz,nx),&
			rx(0:nz,nx),&
			rz(0:nz,nx),&
			rr(0:nz,nx),&
			derx(nz,nx),&
			derz(nz,nx),&
			mu(nz,nx),&
			lambmu(nz,nx),&
			rtaustss(nz,nx),&
			taupstinss(nz,nx),&
			tausstinss(nz,nx)


	!Local variables
	integer::k,i


	!$OMP PARALLEL PRIVATE(k) 
	!$OMP DO 
		do i=1,nx
			do k=itopo(i)-npzn,nz

				rx(k,i)=rx(k,i)+dt*(-1.0*rtaustss(k,i))*(&
				lambmu(k,i)*(taupstinss(k,i)-1)*derx(k,i)&
				)

				rz(k,i)=rz(k,i)+dt*(-1.0*rtaustss(k,i))*(&
				(lambmu(k,i)*(taupstinss(k,i)-1)-2.0*mu(k,i)*(tausstinss(k,i)-1))*derz(k,i)&
				)

				rr(k,i)=rr(k,i)+dt*(-1.0*rtaustss(k,i))*r(k,i)
			enddo
		enddo
	
	!$OMP END DO
	!$OMP END PARALLEL
	return
	end subroutine calprxx


!******************************************************************
	subroutine calprzz(nx,nz,npzn,dt,derx,derz,r,rx,rz,rr,&
				mu,lambmu,rtaustss,taupstinss,tausstinss,itopo)

	implicit none
	!Dummy variables
	integer::nx,nz,npzn,&
			itopo(nx)
	real::	dt,r(nz,nx),&
			rx(0:nz,nx),&
			rz(0:nz,nx),&
			rr(0:nz,nx),&
			derx(nz,nx),&
			derz(nz,nx),&
			mu(nz,nx),&
			lambmu(nz,nx),&
			rtaustss(nz,nx),&
			taupstinss(nz,nx),&
			tausstinss(nz,nx)

	!Local variables
	integer::k,i


	!$OMP PARALLEL PRIVATE(k) 
	!$OMP DO 
		do i=1,nx
			do k=itopo(i)-npzn,nz

				rz(k,i)=rz(k,i)+dt*(-1.0*rtaustss(k,i))*(&
				lambmu(k,i)*(taupstinss(k,i)-1)*derz(k,i)&
				)

				rx(k,i)=rx(k,i)+dt*(-1.0*rtaustss(k,i))*(&
				(lambmu(k,i)*(taupstinss(k,i)-1)-2.0*mu(k,i)*(tausstinss(k,i)-1))*derx(k,i)&
				)

				rr(k,i)=rr(k,i)+dt*(-1.0*rtaustss(k,i))*r(k,i)
			enddo
		enddo
	
	!$OMP END DO
	!$OMP END PARALLEL
	return
	end subroutine calprzz

!****************************************************************
	subroutine calprxz(nx,nz,npzn,dt,derx,derz,r,rx,rz,rr,&
				mu,lambmu,rtaustss,taupstinss,tausstinss,itopo)

	implicit none
	!Dummy variables
	integer::nx,nz,npzn,&
			itopo(nx)
	real::	dt,r(nz,nx),&
			rx(0:nz,nx),&
			rz(0:nz,nx),&
			rr(0:nz,nx),&
			derx(nz,nx),&
			derz(nz,nx),&
			mu(nz,nx),&
			lambmu(nz,nx),&
			rtaustss(nz,nx),&
			taupstinss(nz,nx),&
			tausstinss(nz,nx)

	!Local variables
	integer::k,i


	!$OMP PARALLEL PRIVATE(k) 
	!$OMP DO 
		do i=1,nx
			do k=itopo(i)-npzn,nz

				rx(k,i)=rx(k,i)+dt*(-1.0*rtaustss(k,i))*(mu(k,i)*(tausstinss(k,i)-1)*derx(k,i))

				rz(k,i)=rz(k,i)+dt*(-1.0*rtaustss(k,i))*(mu(k,i)*(tausstinss(k,i)-1)*derz(k,i))

				rr(k,i)=rr(k,i)+dt*(-1.0*rtaustss(k,i))*r(k,i)

			enddo
		enddo
	
	!$OMP END DO
	!$OMP END PARALLEL
	return
	end subroutine calprxz


!******************************************************************
	subroutine calsumvel(nx,nz,vx,vz,v)
	implicit none
	!Dummy variables
	integer::nx,nz
	real::	v(-4:nz+5,-4:nx+5),&
    		vx(0:nz,nx),vz(0:nz,nx)
	!Local variables
	integer::k,i
	!$OMP PARALLEL 
	!$OMP DO 
	do i=1,nx
		do k=1,nz
			v(k,i)=vx(k,i)+vz(k,i)
		enddo
	enddo
	!$OMP END DO
	!$OMP END PARALLEL
	return
	end subroutine calsumvel

!******************************************************************
	subroutine calsumstr(nx,nz,strx,strz,strr,str)
	implicit none
	!Dummy variables
	integer::nx,nz
	real::	str(-4:nz+5,-4:nx+5),&
    		strx(0:nz,nx),strz(0:nz,nx),&
			strr(0:nz,nx)
	!Local variables
	integer::k,i

	!$OMP PARALLEL 
	!$OMP DO 
	do i=1,nx
		do k=1,nz
			str(k,i)=strx(k,i)+strz(k,i)+strr(k,i)
		enddo
	enddo
	!$OMP END DO
	!$OMP END PARALLEL
	return
	end subroutine calsumstr
!******************************************************************
	subroutine calsummr(nx,nz,rx,rz,rr,r)
	implicit none
	!Dummy variables
	integer::nx,nz
	real::	r(nz,nx),&
    		rx(0:nz,nx),rz(0:nz,nx),&
			rr(0:nz,nx)
	!Local variables
	integer::k,i


	!$OMP PARALLEL 
	!$OMP DO 
	do i=1,nx
		do k=1,nz
			r(k,i)=rx(k,i)+rz(k,i)+rr(k,i)
		enddo
	enddo
	!$OMP END DO
	!$OMP END PARALLEL
	return
	end subroutine calsummr

!***********************************************************************************
	subroutine addsrc(nx,nz,isxe,isze,it,f0,dx,dz,dt,u)
	use constant

	implicit none
	!Dummy variables
	integer::nx,nz,isxe,isze,it
	real::f0,dx,dz,dt
	real::u(0:nz,nx)

	!Local variables
	integer::i,k,hw1,hw2,idx1,idx2,m,n,ierr
	integer::n1,n2
	real,parameter::sigma=1.5
	real,parameter::alpha=0.5
	real::wavelet,atten
    real,allocatable::g2c(:,:)

!	n1=50
!	n2=50

	n1=nz
	n2=nx

!	write(*,*)'isxe,isze inside'
!	write(*,*)isxe,isze
!	pause

	hw1=nint((n1+1)/2.0)
	hw2=nint((n2+1)/2.0)


	allocate(g2c(n1,n2),STAT=ierr)
	g2c=0.0
	idx1=1
	idx2=1

	wavelet=exp(-(pi*f0*(it*dt-1.2/f0))**2)*(1-2*(pi*f0*(it*dt-1.2/f0))**2)

	call get_gauss_filter_coe(g2c,sigma,n1,n2)

	do n=isxe-(hw2-1),isxe+(hw2-1)
		idx1=1
		do m=isze-(hw1-1),isze+(hw1-1)
			if(n>=1.and.n<=nx.and.m>=1.and.m<=nz.and.idx2<=n2.and.idx1<=n1)then
				u(m,n)=wavelet*g2c(idx1,idx2)*dt/(dz**2)+u(m,n)
			endif
			idx1=idx1+1 
		enddo
		idx2=idx2+1
	enddo
	deallocate(g2c)
	return
	end subroutine addsrc




!======================================================================*
      subroutine get_gauss_filter_coe(g2c,sigma,n1,n2)
      implicit none
      !Dummy variables
      integer::n1,n2
	  real::sigma,g2c(n1,n2)
      !Local variables
      integer::i,j,k,l,ixcen,izcen,ierr
      real,parameter::pi=3.1415926
      real::dis,sigma2,sum_tmp
      sum_tmp=0.0
      sigma2=sigma**2
      ixcen=(n1+1)/2
      izcen=(n2+1)/2
      do j=1,n2
          do i=1,n1
	          dis=abs(j-izcen)**2+abs(i-ixcen)**2
              g2c(i,j)=1/(2*pi*sigma2)*exp(-dis/(2*sigma2))
          enddo
      enddo

      do j=1,n2
	      do i=1,n1
		    sum_tmp=sum_tmp+g2c(i,j)
	      enddo
      enddo

      !pay attenuation
      g2c=g2c/sum_tmp


!        open(10,file='gc_test',action='write',&
!       form='unformatted',access='direct',status='replace',recl=n1*n2)
!		 write(10,rec=1)((g2c(k,i),k=1,n1),i=1,n2)
!       close (10)

!      sum_tmp=0.0
!      do j=1,n2
!	      do i=1,n1
!		    sum_tmp=sum_tmp+g2c(i,j)
!			write(*,*)'i=.j=,g2c(i,j)=',i,j,g2c(i,j)
!	      enddo
!		enddo
!
!	  write(*,*)'sum=',sum_tmp

      return
	 end subroutine get_gauss_filter_coe



!**********************************************************************************
	subroutine  subsum(u1,u2,u,nx,nz)
	implicit none
	!Dummy variables
	integer::nx,nz
	real:: u1(0:nz,nx),u2(0:nz,nx),u(0:nz+5,-4:nx+5)
	!Local variables
	integer::i,k

	do i=1,nx
		do k=1,nz
			u(k,i)=u1(k,i)+u2(k,i)
		enddo
	enddo

	return
	end subroutine subsum



!*********************************************************************************

!	subroutine get_record()
!
!	end subroutine get_record


!************************************************************************************
        subroutine write_currshot_disk(record,shot_fn1,record_acc,shot_fn3,&
								currshot_name,currshot_no,&
                                 nsx,nsz0,currshot_range,currshot_xmin,&
                                 currshot_xmax,nt,dx_v,dz_v,dt)

        use header_module               
        implicit none
        !Dummy variables
        integer::currshot_no,nsx,nsz0,&
                 currshot_range,currshot_xmin,&
                 currshot_xmax,nt
        real::record(nt,currshot_range)
        real::record_acc(nt,currshot_range)
        character(len=256)::shot_fn1
        character(len=256)::shot_fn3
        character(len=256)::currshot_name
        real::dx_v,dz_v,dt
        !Local variables
        type(segy)::su_head
        real::offx_tmp,offy_tmp
        integer::i,j,k,irec

        !Initiallization of the su_head
        su_head=segy(0,0,0,0,0,0,0,0,0,0,0,&
                    0,0,0,0,0,0,0,0,0,0,&
                    0,0,0,0,0,0,0,0,0,0,&
                    0,0,0,0,0,0,0,0,0,0,&
                    0,0,0,0,0,0,0,0,0,0,&
                    0,0,0,0,0,0,0,0,0,0,&
                    0,0,0,0,0,0,0,0,0,0,&
                    0,0,0,0,0,0,0,0,0,0)

        su_head%fldr=currshot_no!shot_no
!        su_head%tracl=currshot_range!number of traces per shot in inline direction
        su_head%trid=1!data_type
        su_head%sx=(nsx-1)*dx_v
        su_head%ns=nt
        su_head%dt=dt*1.0E6!us SEG standard
!        su_head%d1=dx_v!Self_defined
!        su_head%f1=0.0
!        su_head%d2=dz_v


		!Pressure
        open(unit=8,file=trim(adjustl(shot_fn1))//trim(adjustl(currshot_name))//'.su',&
        form='unformatted',access='direct',status='replace',recl=(nt+60))


		!Acceleration
        open(unit=9,file=trim(adjustl(shot_fn3))//trim(adjustl(currshot_name))//'.su',&
        form='unformatted',access='direct',status='replace',recl=(nt+60))


        do i=1,currshot_range

            su_head%gx=(currshot_xmin-1+i-1)*dx_v
			su_head%offset=su_head%gx-su_head%sx

!			su_head%cdp=ceiling(nint((su_head%gx-su_head%sx)/(dx_v)))
			su_head%cdp=nint((su_head%gx+su_head%sx)/(dx_v))+1


            write(8,rec=i)su_head,(record(k,i),k=1,nt)!**zy
            write(9,rec=i)su_head,(record_acc(k,i),k=1,nt)!**zy

        enddo
        close(unit=8)
		close(unit=9)

        return
        end subroutine write_currshot_disk

!****************************************************************************************

        subroutine merge_shot_files(nshot,shot_fn1,shotall_fn1,&
									shot_fn2,shotall_fn2,currshot_xmin,currshot_xmax,nt)

         implicit none
          !Dummy variables
         character(len=256)::shot_fn1,shotall_fn1,shot_fn2,shotall_fn2
         integer::nshot,currshot_xmin,currshot_xmax,nt
         !Local variables
         real,allocatable::record(:,:)
         integer::currshot_ran_x,&
                   is,i,j,k,ierr
         character(len=256)::currshot_noc
         character(len=256)::currshot_name
  
        !Initiallization
        currshot_ran_x=(currshot_xmax-currshot_xmin)+1
        allocate(record(nt+60,currshot_ran_x),STAT=ierr)
        record=0.0



		!================merge No.1 shot fuiles=====================
        open(11,file=trim(adjustl(shotall_fn1))//'.su',access='direct',form='unformatted',&
              action='write',status='replace',recl=currshot_ran_x*(nt+60))  

        do is=1,nshot
          record=0.0
          write(currshot_noc,'(I4)')is
          !Need some test here
          write(currshot_name,*)trim(adjustl(shot_fn1))//trim(adjustl(currshot_noc))//'.su'
          write(*,*)trim(currshot_name),' is bening merged'
          
          open(10,file=currshot_name,access='direct',form='unformatted',&
               action='read',recl=currshot_ran_x*(nt+60))

          !Try to read in this way for better efficiency
           read(10,rec=1)((record(k,i),k=1,nt+60),i=1,currshot_ran_x)
          close (10,status='delete')
!          close (10)

          write(11,rec=is)((record(k,i),k=1,nt+60),i=1,currshot_ran_x)
        enddo

        close(11)


		!================merge No.2 shot fuiles=====================

        open(11,file=trim(adjustl(shotall_fn2))//'.su',access='direct',form='unformatted',&
              action='write',status='replace',recl=currshot_ran_x*(nt+60))  

        do is=1,nshot
          record=0.0
          write(currshot_noc,'(I4)')is
          !Need some test here
          write(currshot_name,*)trim(adjustl(shot_fn2))//trim(adjustl(currshot_noc))//'.su'
          write(*,*)trim(currshot_name),' is bening merged'
          
          open(10,file=currshot_name,access='direct',form='unformatted',&
               action='read',recl=currshot_ran_x*(nt+60))

          !Try to read in this way for better efficiency
           read(10,rec=1)((record(k,i),k=1,nt+60),i=1,currshot_ran_x)
          close (10,status='delete')
!          close (10)

          write(11,rec=is)((record(k,i),k=1,nt+60),i=1,currshot_ran_x)
        enddo

        close(11)


        deallocate(record,STAT=ierr)
        return
        end subroutine merge_shot_files
!=========================================================================

    subroutine coefficient_2nd(order_even,coefficients)
    !use IMSL
    implicit none
     !Dummy variables
     integer::order_even
     real::coefficients(order_even/2)
    !routine variables
    real,allocatable::a(:,:),b(:)
    real::i,j
    real::fact
     integer::k1
     k1=order_even/2
    allocate(a(k1,k1))
    allocate(b(k1))
     !Initiallization
     a=0.0
     b=0.0
    !Computation
    do i=1,k1
      do j=1,k1
        fact=1.0
        a(i,j)=(i**(2*j))/factorical((2*j),fact)
      enddo
    enddo
    !The operation of calling can not be omitted
    call inv()
    deallocate(a)
    deallocate(b)
    return
     contains
       real function factorical(l,fact)
       implicit none
       real::fact
       real::m,l
       do m=1,l
         fact=fact*m
       enddo
      factorical=fact
      !write(*,*)fact
      return
      end function factorical
           
      subroutine inv()
      implicit none
      integer::i1,j1
      real::ftmp1,ftmp5
      real,allocatable::inv_a(:,:),temp(:),ftmp2(:,:),ftmp3(:,:),ftmp4(:),p1(:)
      allocate(p1(k1))
      allocate(temp(k1))
      allocate(ftmp2(k1,k1))
      allocate(ftmp3(k1,k1))
      allocate(inv_a(k1,k1))
      allocate(ftmp4(k1))
      !Initiallization
      p1=0.0
      temp=0.0
      ftmp2=0.0
      ftmp3=0.0
      ftmp4=0.0
      inv_a=0.0
      !Computation
      do i=1,k1
         do j=1,k1
           if (i==j)then
             inv_a(i,j)=1
           endif
         enddo
      enddo
      do i=1,k1
         temp=a(i,:)
        temp(i)=temp(i)-1
        ftmp1=sum(temp*inv_a(:,i))+1
        ftmp4=inv_a(:,i)
        do j1=1,k1
          do i1=1,k1
              ftmp2(j1,i1)=ftmp4(j1)*temp(i1)
            enddo
        enddo
         ftmp3=matmul(ftmp2,inv_a)
         do j=1,k1
           do j1=1,k1
            ftmp3(j,j1)=ftmp3(j,j1)/ftmp1
          enddo
        enddo
        inv_a=inv_a-ftmp3
      enddo
      do i=1,k1
        b=0
        b(i)=0.5
        p1=matmul(inv_a,b)
        coefficients(i)=p1(1)
      enddo
      deallocate(p1)
       deallocate(temp)
       deallocate(ftmp2)
       deallocate(ftmp3)
       deallocate(inv_a)
       deallocate(ftmp4)
       end subroutine inv
    end subroutine  coefficient_2nd


    subroutine coefficient_1st(order_odd,coefficients)
    !use IMSL
    implicit none
    !Dummy variables
    integer::order_odd
    real::coefficients((order_odd+1)/2)
    !routine variables
    real,allocatable::a(:,:),b(:)
    real::i,j
    real::fact
    integer::k1_1
    k1_1=((order_odd+1)/2)!Pay attenuation and need more thicking
    allocate(a(k1_1,k1_1))
    allocate(b(k1_1))
    !Initiallization
    a=0.0
    b=0.0
    !Computation
    do i=1,k1_1
      do j=1,k1_1
        fact=1.0
        a(i,j)=(i**(2*j-1))/factorical((2*j-1),fact)
      enddo
    enddo
    !The operation of calling can not be omitted
     call inv_1st()
     
     deallocate(a)
     deallocate(b)
     return
     contains
       real function factorical(l,fact)
       implicit none
       real::fact
       real::m,l
       do m=1,l
         fact=fact*m
       enddo
      factorical=fact
      return
      end function factorical
           
      subroutine inv_1st()
      implicit none
      integer::i1,j1
      real::ftmp1,ftmp5
      real,allocatable::inv_a(:,:),temp(:),ftmp2(:,:),ftmp3(:,:),ftmp4(:),p1(:)
      allocate(p1(k1_1))
      allocate(temp(k1_1))
      allocate(ftmp2(k1_1,k1_1))
      allocate(ftmp3(k1_1,k1_1))
      allocate(inv_a(k1_1,k1_1))
      allocate(ftmp4(k1_1))
      !Initiallization
      p1=0.0
      temp=0.0
      ftmp2=0.0
      ftmp3=0.0
      inv_a=0.0
      ftmp4=0.0
      do i=1,k1_1
         do j=1,k1_1
           if (i==j)then
             inv_a(i,j)=1
           endif
         enddo
      enddo
     do i=1,k1_1
       temp=a(i,:)
      temp(i)=temp(i)-1
      ftmp1=sum(temp*inv_a(:,i))+1
      ftmp4=inv_a(:,i)
      do j1=1,k1_1
        do i1=1,k1_1
          ftmp2(j1,i1)=ftmp4(j1)*temp(i1)
        enddo
      enddo
       ftmp3=matmul(ftmp2,inv_a)
       do j=1,k1_1
         do j1=1,k1_1
          ftmp3(j,j1)=ftmp3(j,j1)/ftmp1
        enddo
      enddo
      inv_a=inv_a-ftmp3
     enddo
    do i=1,k1_1
      b=0
      b(i)=0.5
      p1=matmul(inv_a,b)
      coefficients(i)=p1(1)!Pay attenuation
    enddo
    !=======Print some basic information on the screen====
!    write(*,*)'=====The order of time derivatives is======'
!    write(*,*)'          4'
!    write(*,*)'=====The order of 1st spatial derivatives is==='
!    write(*,*)k_1
!    write(*,*)'=====The difference coefficients for 1st spatial derivatives are========'
!    write(*,*)p_1

    deallocate(p1)
    deallocate(temp)
    deallocate(ftmp2)
    deallocate(ftmp3)
    deallocate(inv_a)
    deallocate(ftmp4)
     end subroutine inv_1st
    end subroutine  coefficient_1st

        

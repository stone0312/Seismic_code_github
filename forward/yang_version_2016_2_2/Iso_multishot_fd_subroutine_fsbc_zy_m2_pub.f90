		!*****************************************************************
		subroutine get_record(it,npxn,nrz0,nxe,nx,nz,nt,itopo,vx,record1)
		implicit none
		!Dummy variables
		integer,intent(in)::it,npxn,nrz0,nxe,nx,nz,nt
		integer,intent(in)::itopo(nx)
		real,intent(in)::vx(-4:nz+5,-4:nx+5)
		real,intent(inout)::record1(nt,nxe)
		!Local variables
		integer::i
	
		do i=1,nxe
			record1(it,i)=vx(itopo(i+npxn)+nrz0,i+npxn)
		enddo
		return
		end subroutine get_record


		!***************************************************************
		subroutine write_currtsnap(currtfnpart1,it,nxe,nze,nx,nz,&
					npxn,npzn,u)

		implicit none
		!Dummy variables
		character(len=256),intent(inout)::currtfnpart1
		integer,intent(in)::it,nxe,nze,nx,nz,npxn,npzn
		real,intent(in)::u(-4:nz+5,-4:nx+5)

		!Local variables
		character(len=256)::currtfnall,currt
		integer::i,k

		write(currt,'(I6)')it

		write(currtfnall,*)trim(adjustl(currtfnpart1)),trim(adjustl(currt))

		open(unit=8,file=currtfnall,form='unformatted',access='direct',&
		status='replace',recl=nze*nxe)

		write(8,rec=1)((u(k,i),k=npzn+1,nze+npzn),i=1+npxn,nxe+npxn)
		close(unit=8)

		return
		end subroutine write_currtsnap


		!****************************************************************
        subroutine read_par(par_fn,molmod,vp_fn,vs_fn,rho_fn,topo_fn,&
					snapx_fn,snapz_fn,shot_fn1,shot_fn2,shot_fn3,shot_fn4,&
					nshot,ofsmin,ofsmax,sx0,sz0,rz0,dsx,nx,nz,dx,dz,nt,dt,&
					dtsnap,f0,decay,ppxn,ppxp,ppzn,ppzp,rgau)

        implicit none

        !Dummy variables
        !*Parameter_card_name*
        character(len=256)::par_fn
        !*Read from parameter_card*
        character(len=256)::molmod
        character(len=256)::vp_fn
        character(len=256)::vs_fn
        character(len=256)::rho_fn
        character(len=256)::topo_fn
        character(len=256)::snapx_fn
        character(len=256)::snapz_fn
        character(len=256)::shot_fn1
        character(len=256)::shot_fn2
        character(len=256)::shot_fn3
        character(len=256)::shot_fn4
        integer::nx,nz
        integer::nt
        integer::nshot
		integer::rgau
        real::ofsmin,ofsmax
        real::dx,dz
!        real::pml_thick!(m)
        real::ppxn,ppxp,ppzn,ppzp!(m)
        real::sx0,sz0,rz0!(m)
        real::dsx!(m)
        real::dt,dtsnap!(s)
        real::f0!(Hz)
        real::decay!(s)

        !*variables of MPI interface*
        integer::myid
      
        !Local variables
        integer::i,j,ierr
        character(len=256)::par_jumpper

        open(10,file=par_fn,action='read',status='old',form='formatted',&
        access='sequential')

        read(10,'(A)')par_jumpper
        read(10,'(A)')molmod
        read(10,'(A)')par_jumpper
        read(10,'(A)')vp_fn
        read(10,'(A)')par_jumpper
        read(10,'(A)')vs_fn
        read(10,'(A)')par_jumpper
        read(10,'(A)')rho_fn
        read(10,'(A)')par_jumpper
        read(10,'(A)')topo_fn
        read(10,'(A)')par_jumpper
        read(10,'(A)')snapx_fn
        read(10,'(A)')par_jumpper
        read(10,'(A)')snapz_fn
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
        read(10,*)ofsmin,ofsmax
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
        read(10,*)dt,dtsnap
        read(10,'(A)')par_jumpper
        read(10,'(F)')f0
        read(10,'(A)')par_jumpper
        read(10,'(F)')decay
        read(10,'(A)')par_jumpper
        read(10,*)ppxn,ppxp
        read(10,'(A)')par_jumpper
        read(10,*)ppzn,ppzp!**zy
        read(10,'(A)')par_jumpper
        read(10,'(I)')rgau

        close(10)
		return
        end subroutine read_par

!********subroutine*******************************************
        subroutine read_data(vp_fn,vs_fn,rho_fn,topo_fn,&
							vp,vs,rho,topo,nx,nz,myid)

        use constant
        implicit none
        !Dummy variables
        character(len=256)::vp_fn,vs_fn,rho_fn,topo_fn
        integer::nx,nz,myid
        real::vp(nz,nx),vs(nz,nx),rho(nz,nx),topo(nx)

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

        open(10,file=topo_fn,action='read',&
        form='unformatted',access='direct',status='old',recl=nx)
		read(10,rec=1)(topo(i),i=1,nx)
        close (10)

        return
        end subroutine read_data


!====================================================================================
        subroutine get_currshot_parameters(vp,vs,rho,topo,vpe,vse,rhoe,topoe,&
						itopoe,gc2d,nx,nz,nxe,nze,nxea,nzea,ixmine,ixmaxe,npxn,&
						npxp,npzn,npzp,rgau,isx,isz,isxe,isze,dx,dz)

        use constant

        implicit none
        !Dummy variables
        integer,intent(in out)::isxe,isze,itopoe(nxea)
        integer,intent(in)::nx,nz,nxe,nze,nxea,nzea,&
							npxn,npxp,npzn,npzp,ixmine,&
							ixmaxe,isx,isz,rgau

		
		real,intent(in)::dx,dz
        real,intent(in)::vp(nz,nx),vs(nz,nx),rho(nz,nx),&
						topo(nx)

        real,intent(in out)::vpe(nzea,nxea),vse(nzea,nxea),&
							rhoe(nzea,nxea),topoe(nxea),&
							gc2d(2*rgau-1,2*rgau-1)



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
           topoe(i+npxn)=topo(i_ori)

        enddo


		!Using a simpler way to pa velocity
		do k=1,npzn
			do i=1,nxea
				vpe(k,i)=vpe(npzn+1,i)
				vse(k,i)=vse(npzn+1,i)
				rhoe(k,i)=rhoe(npzn+1,i)
			enddo
		enddo


		do k=npzn+nze+1,nzea
			do i=1,nxea
				vpe(k,i)=vpe(npzn+nz,i)
				vse(k,i)=vse(npzn+nz,i)
				rhoe(k,i)=rhoe(npzn+nz,i)
			enddo
		enddo


		do k=1,nzea
			do i=1,npxn
				vpe(k,i)=vpe(k,npxn+1)
				vse(k,i)=vse(k,npxn+1)
				rhoe(k,i)=rhoe(k,npxn+1)
				topoe(i)=topoe(npxn+1)
			enddo
		enddo



		do k=1,nzea
			do i=npxn+nxe+1,nxea
				vpe(k,i)=vpe(k,npxn+nxe)
				vse(k,i)=vse(k,npxn+nxe)
				rhoe(k,i)=rhoe(k,npxn+nxe)
				topoe(i)=topoe(npxn+nxe)
			enddo
		enddo


		!Pay attenuation to the defination of itopoe
		do i=1,nxea
			itopoe(i)=nint(topoe(i)/dz)+npzn
		enddo


        !Compute currshot source_position
        isxe=isx-ixmine+npxn+1!**zy
        isze=isz+itopoe(isxe)


		!Get 2D gauss coefficietns
		call get_gauss_filter_coe(gc2d,2.0,2*rgau-1,2*rgau-1)


!        open(10,file='vpe_test',action='write',&
!        form='unformatted',access='direct',status='replace',recl=nxe*nze)
!		 write(10,rec=1)((vpe(k,i),k=npzn+1,npzn+nze),i=npxn+1,npxn+nxe)
!        close (10)


!       open(10,file='vpe_test',action='write',&
!	    form='unformatted',access='direct',status='replace',recl=nxea*nzea)
!		write(10,rec=1)((vpe(k,i),k=1,nzea),i=1,nxea)
!       close (10)

!        open(10,file='vse_test',action='write',&
!        form='unformatted',access='direct',status='replace',recl=nxe*nze)
!		 write(10,rec=1)((vse(k,i),k=npzn+1,npzn+nze),i=npxn+1,npxn+nxe)
!        close (10)
!
!
!        open(10,file='rhoe_test',action='write',&
!        form='unformatted',access='direct',status='replace',recl=nxe*nze)
!		 write(10,rec=1)((rhoe(k,i),k=npzn+1,npzn+nze),i=npxn+1,npxn+nxe)
!        close (10)
!
!
!        open(10,file='topoe_test',action='write',&
!        form='unformatted',access='direct',status='replace',recl=nxe)
!		 write(10,rec=1)(topoe(i),i=npxn+1,npxn+nxe)
!        close (10)

        return
        end subroutine get_currshot_parameters


!======================================================================================
! SUBROUTINE SUBSTAGGER: BUILD STAGGERED GRIDS
!
! I: VP, VS, RHO, MU, LAMB, N1, N2
! O: MUSG,MUSGE,LAMBSG,LAMBMUSG,BUSG1,BUSG2,BUSG1E,BUSG2E
!=======================================================================================

	 subroutine substagger(vp,vs,rho,mu,lamb,musg,musge,lambsg,&
	 						lambmusg,busg1,busg2,busg1e,busg2e,nx,nz,itopoe)
	
	 implicit none
	 !Dummy variables
	 integer::nx,nz,itopoe(nx)
	 real::vp(nz,nx),vs(nz,nx),rho(nz,nx),&
			mu(nz,nx),lamb(nz,nx),musg(nz,nx),&
			lambsg(nz,nx),musge(nz,nx),lambmusg(nz,nx),&
			busg1(nz,nx),busg2(nz,nx),busg1e(nz,nx),&
			busg2e(nz,nx)

	 !Local variables
	 integer::i,k

	! #########################################################################
	! Lame parameters
	 do i=1,nx
	 	do k=1,nz
	 		mu(k,i)=rho(k,i)*vs(k,i)*vs(k,i)
	 		lamb(k,i)=rho(k,i)*vp(k,i)*vp(k,i)-2.*mu(k,i)
	 	end do
	 end do

	! #########################################################################

	! Busg1
	 do i=1,nx-1
	 	do k=1,nz-1
	 		busg1(k,i)=0.25*(rho(k,i)+rho(k+1,i)+rho(k+1,i+1)+rho(k,i+1))
	 	end do	
	 end do


	 do i=1,nx
	 	busg1(nz,i)=busg1(nz-1,i)
	 end do
	
	
	 do k=1,nz
	 	busg1(k,nx)=busg1(k,nx-1)
	 end do	

	! convert to 1/rho
	
	 do i=1,nx
	 	do k=1,nz
	 		if (busg1(k,i).eq.0.) stop 'FIND A DENSITY EQUAL TO 0'
	 			busg1(k,i)=1./busg1(k,i)
	 	end do
	 end do
	
	
	! Effective parameters
	 do i=1,nx
	 	do k=2,nz
	 		busg1e(k,i)=0.5*(busg1(k,i)+busg1(k-1,i))
	 	end do
	 end do
	
	 do i=1,nx
	 	busg1e(1,i)=busg1(1,i)
	 end do

	! #########################################################################
	! RHOSG2
	
	 do i=1,nx
	 	do k=1,nz
	 		busg2(k,i)=rho(k,i)
	 	end do
	 end do


	! convert to 1/rho
	 do i=1,nx
	 	do k=1,nz
	 		if (busg2(k,i).eq.0.) stop 'FIND A DENSITY EQUAL TO 0'
	 		busg2(k,i)=1./busg2(k,i)
	 	end do
	 end do
	
	! Effective parameters
	 do i=1,nx-1
	 	do k=1,nz
	 		busg2e(k,i)=0.5*(busg2(k,i)+busg2(k,i+1))
	 	end do
	 end do

	 do k=1,nz
	 	busg2e(k,nx)=busg2(k,nx)
	 end do
	



	! #########################################################################
	! Musg
	 do i=1,nx
	 	do k=1,nz-1
	 		musg(k,i)=0.5*(mu(k+1,i)+mu(k,i))
	 	end do
	 end do
	
	 do i=1,nx
	 	musg(nz,i)=mu(nz,i)
	 end do


	! Effective parameters
	!Modified by zy
	 do i=1,nx-1
	 	do k=2,nz

			 if(musg(k,i)/=0.0.and.musg(k-1,i)/=0.0.and.musg(k,i+1)/=0.0.and.musg(k-1,i+1)/=0.0)then
	 			musge(k,i)=4./(1./musg(k,i)+1./musg(k,i+1)+1./musg(k-1,i+1)+1./musg(k-1,i))

			else
				musge(k,i)=0.0

			endif

	 	end do
	 end do


	 do i=1,nx
	 	musge(1,i)=musg(1,i)
	 end do
	
	 do k=1,nz
	 	musge(k,nx)=musg(k,nx)
	 end do


	! #########################################################################
	! Lambsg
	 do k=1,nz
	 	do i=1,nx-1
	 		lambsg(k,i)=0.5*(lamb(k,i+1)+lamb(k,i))
	 	end do
	 end do
	
	 do k=1,nz
	 	lambsg(k,nx)=lamb(k,nx)
	 end do



	! #########################################################################
	! Lambmusg
	 do i=1,nx-1
	 	do k=1,nz
	 		lambmusg(k,i)=0.5*(lamb(k,i)+2*mu(k,i)+lamb(k,i+1)+2*mu(k,i+1))
	 	end do
	 end do
	
	 do k=1,nz
	 	lambmusg(k,nx)=lamb(k,nx)+2*mu(k,nx)
	 end do

	!Refining parameters to adjust PML
	do i=1,nx
		do k=itopoe(i),2,-1
			vp(k-1,i)=vp(k,i)
			vs(k-1,i)=vs(k,i)
			rho(k-1,i)=rho(k,i)
			lamb(k-1,i)=lamb(k,i)
			lambsg(k-1,i)=lambsg(k,i)
			lambmusg(k-1,i)=lambmusg(k,i)
			mu(k-1,i)=mu(k,i)
			musg(k-1,i)=musg(k,i)
			musge(k-1,i)=musge(k,i)
			busg1(k-1,i)=busg1(k,i)
			busg2(k-1,i)=busg2(k,i)
			busg1e(k-1,i)=busg1e(k,i)
			busg2e(k-1,i)=busg2e(k,i)
		enddo
	enddo


!        open(10,file='vpe_test',action='write',&
!        form='unformatted',access='direct',status='replace',recl=nxe*nze)
!		 write(10,rec=1)((vpe(,i),k=npzn+1,npzn+nze),i=npxn+1,npxn+nxe)
!        close (10)
!
!
!        open(10,file='vse_test',action='write',&
!        form='unformatted',access='direct',status='replace',recl=nxe*nze)
!		 write(10,rec=1)((vse(k,i),k=npzn+1,npzn+nze),i=npxn+1,npxn+nxe)
!        close (10)
!
!
!        open(10,file='rhoe_test',action='write',&
!        form='unformatted',access='direct',status='replace',recl=nxe*nze)
!		 write(10,rec=1)((rhoe(k,i),k=npzn+1,npzn+nze),i=npxn+1,npxn+nxe)
!        close (10)


!        open(10,file='vpe_test_new',action='write',&
!        form='unformatted',access='direct',status='replace',recl=nx*nz)
!		 write(10,rec=1)((vp(k,i),k=1,nz),i=1,nx)
!        close (10)

!
!        open(10,file='vse_test_new',action='write',&
!        form='unformatted',access='direct',status='replace',recl=nx*nz)
!		 write(10,rec=1)((vs(k,i),k=1,nz),i=1,nx)
!        close (10)
!
!
!        open(10,file='rho_test_new',action='write',&
!        form='unformatted',access='direct',status='replace',recl=nx*nz)
!		 write(10,rec=1)((rho(k,i),k=1,nz),i=1,nx)
!        close (10)


        open(10,file='musge_test',action='write',&
        form='unformatted',access='direct',status='replace',recl=nx*nz)
		 write(10,rec=1)((musge(k,i),k=1,nz),i=1,nx)
        close (10)


	 return
	 end subroutine substagger

	subroutine modipapml(nx,nz,rgau,itopoe,gc2d,vp,vs,rho,lamb,lambsg,&
				lambmusg,mu,musg,musge,busg1,busg2,busg1e,busg2e)

	implicit none
	!Dummy variables
	integer::nx,nz,rgau,itopoe(nx)
	real::vp(nz,nx),vs(nz,nx),rho(nz,nx),&
		mu(nz,nx),lamb(nz,nx),musg(nz,nx),&
		lambsg(nz,nx),musge(nz,nx),lambmusg(nz,nx),&
		busg1(nz,nx),busg2(nz,nx),busg1e(nz,nx),&
		busg2e(nz,nx),gc2d(2*rgau-1,2*rgau-1)

	!Local variables
	integer::i,k,idx1,idx2,m,n
	real::sutmp
	
!	write(*,*)gc2d
!	write(*,*)itopoe

	call gauss_filter_2d(nx,nz,rgau,itopoe,gc2d,vp)
	call gauss_filter_2d(nx,nz,rgau,itopoe,gc2d,vs)
	call gauss_filter_2d(nx,nz,rgau,itopoe,gc2d,rho)
	call gauss_filter_2d(nx,nz,rgau,itopoe,gc2d,mu)
	call gauss_filter_2d(nx,nz,rgau,itopoe,gc2d,lamb)
	call gauss_filter_2d(nx,nz,rgau,itopoe,gc2d,musg)
	call gauss_filter_2d(nx,nz,rgau,itopoe,gc2d,lambsg)
	call gauss_filter_2d(nx,nz,rgau,itopoe,gc2d,musge)
	call gauss_filter_2d(nx,nz,rgau,itopoe,gc2d,lambmusg)
	call gauss_filter_2d(nx,nz,rgau,itopoe,gc2d,busg1)
	call gauss_filter_2d(nx,nz,rgau,itopoe,gc2d,busg2)
	call gauss_filter_2d(nx,nz,rgau,itopoe,gc2d,busg1e)
	call gauss_filter_2d(nx,nz,rgau,itopoe,gc2d,busg2e)


!	open(10,file='vpe_test_after_sm',action='write',&
!	form='unformatted',access='direct',status='replace',recl=nx*nz)
!	write(10,rec=1)((vp(k,i),k=1,nz),i=1,nx)
!	close (10)

	return
	end subroutine modipapml


	subroutine gauss_filter_2d(nx,nz,rgau,itopoe,gc2d,u)

	implicit none
	!Dummy variables
	integer::nx,nz,rgau,itopoe(nx)
	real::gc2d(2*rgau-1,2*rgau-1),u(nz,nx)
	!Local variables
	integer::i,k,idx1,idx2,m,n
	real::sutmp


	!$OMP PARALLEL PRIVATE(idx1,idx2,sutmp)
	!$OMP DO 

	do i=1,nx
		do k=1,itopoe(i)

			idx1=1
			idx2=1
			sutmp=0.0

			do n=i-(rgau-1),i+(rgau-1)
				idx1=1
				do m=k-(rgau-1),k+(rgau-1)
					if(n>=1.and.n<=nx.and.m>=1.and.m<=nz&
					.and.idx2<=2*rgau-1.and.idx1<=2*rgau-1)then
						sutmp=sutmp+u(m,n)*gc2d(idx1,idx2)
					endif
					idx1=idx1+1 
				enddo
				idx2=idx2+1
			enddo

			u(k,i)=sutmp

		enddo
	enddo

	!$OMP END DO
	!$OMP END PARALLEL

	return
	end subroutine gauss_filter_2d

!*********************************************************************
	subroutine susmooth2(n1,n2,r1,r2,u)
	implicit none
	!Dummy variables
	integer::n1,n2
	real::r1,r2,u(n1,n2)

	!Local variables
	integer::i,k,nmax,ierr
	real,allocatable::tripd_d(:),&
	tripd_e(:),tripd_f(:)
	real::r1r,r2r

	nmax=max(n1,n2)
!	write(*,*)nmax

	allocate(tripd_d(nmax),STAT=ierr)
	allocate(tripd_e(nmax),STAT=ierr)
	allocate(tripd_f(nmax),STAT=ierr)

	tripd_d=0.0
	tripd_e=0.0
	tripd_f=0.0

!	write(*,*)r1,r2
!	pause

	r1r=r1**2*0.25
	r2r=r2**2*0.25

!	write(*,*)'here1'
	call smooth2(n1,n2,nmax,r1r,r2r,tripd_d,tripd_e,tripd_f,u)


	deallocate(tripd_d)
	deallocate(tripd_e)
	deallocate(tripd_f)

	return
	end subroutine susmooth2


	!**************************************************************
	subroutine smooth2(n1,n2,nmax,r1,r2,tripd_d,tripd_e,tripd_f,vel)
	implicit none
	!Dummy variables
	integer::n1,n2,nmax
	real::r1,r2,vel(0:n1-1,0:n2-1),&
		tripd_d(0:nmax-1),&
		tripd_e(0:nmax-1),&
		tripd_f(0:nmax-1)

	!Local variables
	integer::ix,iz

	!Slow dimension
!!$OMP PARALLEL
!!$OMP DO 
	do iz=0,n1-1
		do ix=0,n2-2
			
			tripd_d(ix)=1.0+r2*2.0
			tripd_e(ix)	=-r2
			tripd_f(ix)=vel(iz,ix)
		enddo
			
		tripd_d(0)=tripd_d(0)-r2
		tripd_d(n2-1)=1.0+r2
		tripd_f(n2-1)=vel(iz,n2-1)
		call tripd(n2,nmax,tripd_d,tripd_e,tripd_f)
		do ix=0,n2-1
			vel(iz,ix)=tripd_f(ix)
		enddo

	enddo
!!$OMP END DO
!!$OMP END PARALLEL

!!$OMP PARALLEL
!!$OMP DO 
	!Fast dimension	
	do ix=0,n2-1
		do iz=0,n1-3

			tripd_d(iz)=1.0+r1*2.0
			tripd_e(iz)=-r1
			tripd_f(iz)=vel(iz+1,ix)

		enddo

		tripd_f(0)=tripd_f(0)+r1*vel(0,ix)
		tripd_d(n1-2)=1.0+r1
		tripd_f(n1-2)=vel(n1-1,ix)

		call  tripd(n1-1,nmax,tripd_d,tripd_e,tripd_f)

		do iz=0,n1-2
			vel(iz+1,ix)=tripd_f(iz)
		enddo

	enddo
!!$OMP END DO
!!$OMP END PARALLEL


	return
	end subroutine smooth2


	!*****************************************
	subroutine tripd(n,nmax,d,e,b)
	implicit none
	!Dummy variables
	integer::n,nmax
	!*Trying to incorporating with C*
	real::d(0:nmax-1),e(0:nmax-1),b(0:nmax-1)
	!Local variables
	integer::k
	real::temp

	!*decomposition*!

!	!$OMP PARALLEL
!	!$OMP DO 
	do k=1,n-1
		temp=e(k-1)
		e(k-1)=temp/d(k-1)
		d(k)=d(k)-temp*e(k-1)
	enddo
!	!$OMP END DO
!	!$OMP END PARALLEL

	!*substitution*!
!	!$OMP PARALLEL
!	!$OMP DO 
	do k=1,n-1
		b(k)=b(k)-e(k-1)*b(k-1)
	enddo
!	!$OMP END DO
!	!$OMP END PARALLEL

	b(n-1)=b(n-1)/d(n-1)

!	!$OMP PARALLEL
!	!$OMP DO 
	do k=n-1,1,-1
		b(k-1)=b(k-1)/d(k-1)-e(k-1)*b(k)
	enddo
!	!$OMP END DO
!	!$OMP END PARALLEL

	return
	end subroutine tripd

!**********************************************************************************
! SUBROUTINE CSPONGE: BUILD PML FUNCTIONS FOR THE DIFFERENT STAGGERED GRIDS      
! I: NSP, H
! SPG, SPG1, SPG2
! Modified by zy
!**********************************************************************************
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
	real::facxn,facxp,faczn,faczp



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


!********************************************************************************
 subroutine calfsbcv(vx,vz,lambsg,musgef,lambmusg,npxn,&
 					nxe,nze,nx,nz,itopo)
 implicit none

 !Dummy variables
 integer::npxn,nxe,nze,nx,nz,itopo(nx)
 real::vx(0:nz,1:nx),vz(0:nz,1:nx),lambsg(nz,nx),&
! real::vx(-4:nz+5,-4:nx+5),vz(-4:nz+5,-4:nx+5),lambsg(nz,nx),&
 		musgef(nz,nx),lambmusg(nz,nx)
 !Local variables
 integer::ix,izh


! TAUzz=0 ---> compute vz above free surface 
!(this is required to compute Tauxx at free surface)


 	do ix=1,nx-1
 		izh=itopo(ix)
 		vz(izh-1,ix)=vz(izh,ix)+(lambsg(izh,ix)/lambmusg(izh,ix))*(vx(izh,ix+1)-vx(izh,ix))
 	end do

! TAUxz=0 <---> TAUxz(izh-1)=-TAUxz(izh)
! Useless, because we compute Tauxz half a grid point below topo with 2nd-order accurate
! operator for vertical derivative. In this case, we only need vx at free surface
! In addition, computing this will cause stability problem

!	 do ix=2,nx
!		izh=itopo(ix)
!	  	vx(izh-1,ix)=vx(izh+1,ix)+vz(izh-1,ix)-vz(izh-1,ix-1)+vz(izh,ix)-vz(izh,ix-1)
!	end do
!	vx(itopo(1),1)=vx(itopo(2),2)

 return
 end


!****************************************************************************
! SUBROUTINE SUBFREESTAU2ZZ: Set tau_zz=0 at free surface
!****************************************************************************
	subroutine calfsbcstrzz(tauzz,npxn,nxe,nze,nx,nz,itopo)
	implicit none
	!Dummy variables
	integer::npxn,nxe,nze,nx,nz,itopo(nx)
	real::tauzz(0:nz,nx)
	!Local variables
	integer::ix,iz,izh


! Useless, because we compute vz half a grid point below topo with 2nd-order accurate
! operator for vertical derivative. In this case, we only need tauzz at free surface
!		tauzz(izh-1,ix)=-tauzz(izh+1,ix)


	!Normal stress
	do ix=1,nx
		izh=itopo(ix)
		do iz=0,izh
			tauzz(iz,ix)=0.0
		enddo
	 end do

	return
	end




!*****************************************************************************
! Set tau_xz=0 at free surface by image theory
!*****************************************************************************
	subroutine calfsbcstrxz(tauxz,npxn,nxe,nze,nx,nz,itopo)
	implicit none
	!Dummy variables
	integer::npxn,nxe,nze,nx,nz,itopo(nx)
	real::tauxz(0:nz,nx)
	!Local variables
	integer::ix,izh

	do ix=1,nx
		izh=itopo(ix)
		tauxz(izh-1,ix)=-tauxz(izh,ix)
	end do

	return
	end

!*********************************************************************************
	subroutine calux(xflag,nx,nz,dx,dt,u,ux,pa,itopo)
	use constant
	implicit none
	!Dummy variables

	integer::xflag,nx,nz,&
			itopo(nx)

	real::dx,dt,&
		u(-4:nz+5,-4:nx+5),&
    	ux(0:nz,nx),pa(nz,nx)


	!Local variables
	integer::k,i



	if(xflag==0)then
		
!$OMP PARALLEL PRIVATE(k) 
!$OMP DO 
		do i=1,nx
			do k=itopo(i),nz
				ux(k,i)=ux(k,i)+(dt/dx)*pa(k,i)*&
				(&
					a100*(u(k,i)-u(k,i-1))+&
					a101*(u(k,i+1)-u(k,i-2))+&
					a102*(u(k,i+2)-u(k,i-3))+&
					a103*(u(k,i+3)-u(k,i-4))+&
					a104*(u(k,i+4)-u(k,i-5))&
				)
			enddo
		enddo

!$OMP END DO
!$OMP END PARALLEL

	else if(xflag==1)then

!$OMP PARALLEL PRIVATE(k) 
!$OMP DO 
		do i=1,nx
			do k=itopo(i),nz
				ux(k,i)=ux(k,i)+(dt/dx)*pa(k,i)*&
				(&
					a100*(u(k,i+1)-u(k,i))+&
					a101*(u(k,i+2)-u(k,i-1))+&
					a102*(u(k,i+3)-u(k,i-2))+&
					a103*(u(k,i+4)-u(k,i-3))+&
					a104*(u(k,i+5)-u(k,i-4))&
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
	end subroutine calux

!******************************************************************************
	subroutine caluxabs(xflag,nx,nz,npzn,dx,dt,u,ux,pa,itopo)
	use constant
	implicit none
	!Dummy variables

	integer::xflag,nx,nz,npzn,itopo(nx)

	real::dx,dt,&
		u(-4:nz+5,-4:nx+5),&
    	ux(0:nz,nx),pa(nz,nx)
	!Local variables
	integer::k,i


	if(xflag==0)then
		
!$OMP PARALLEL PRIVATE(k) 
!$OMP DO 
		do i=1,nx
			do k=itopo(i)-npzn,nz
				ux(k,i)=ux(k,i)+(dt/dx)*pa(k,i)*&
				(&
					a100*(u(k,i)-u(k,i-1))+&
					a101*(u(k,i+1)-u(k,i-2))+&
					a102*(u(k,i+2)-u(k,i-3))+&
					a103*(u(k,i+3)-u(k,i-4))+&
					a104*(u(k,i+4)-u(k,i-5))&
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
				ux(k,i)=ux(k,i)+(dt/dx)*pa(k,i)*&
				(&
					a100*(u(k,i+1)-u(k,i))+&
					a101*(u(k,i+2)-u(k,i-1))+&
					a102*(u(k,i+3)-u(k,i-2))+&
					a103*(u(k,i+4)-u(k,i-3))+&
					a104*(u(k,i+5)-u(k,i-4))&
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
	end subroutine caluxabs

!******************************************************************************
	subroutine caluz(zflag,nx,nz,dz,dt,u,uz,pa,itopo)
	use constant
	implicit none
	!Dummy variables
	integer::zflag,nx,nz,&
			itopo(nx)
	real::	dz,dt,&
			u(-4:nz+5,-4:nx+5),&
    		uz(0:nz,nx),pa(nz,nx)

	!Local variables
	integer::k,i

	if(zflag==0)then

!$OMP PARALLEL PRIVATE(k) 
!$OMP DO 
		do i=1,nx
			do k=itopo(i)+5,nz
				uz(k,i)=uz(k,i)+(dt/dz)*pa(k,i)*&
				(&
				a100*(u(k,i)-u(k-1,i))+&
				a101*(u(k+1,i)-u(k-2,i))+&
				a102*(u(k+2,i)-u(k-3,i))+&
				a103*(u(k+3,i)-u(k-4,i))+&
				a104*(u(k+4,i)-u(k-5,i))&
				)
			enddo
		enddo
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL PRIVATE(k) 
!$OMP DO 
		do i=1,nx
			k=itopo(i)+4
			uz(k,i)=uz(k,i)+(dt/dz)*pa(k,i)*&
			(&
			a80*(u(k,i)-u(k-1,i))+&
			a81*(u(k+1,i)-u(k-2,i))+&
			a82*(u(k+2,i)-u(k-3,i))+&
			a83*(u(k+3,i)-u(k-4,i))&
			)
		enddo
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL PRIVATE(k) 
!$OMP DO 
		do i=1,nx
			k=itopo(i)+3
			uz(k,i)=uz(k,i)+(dt/dz)*pa(k,i)*&
			(&
			a60*(u(k,i)-u(k-1,i))+&
			a61*(u(k+1,i)-u(k-2,i))+&
			a62*(u(k+2,i)-u(k-3,i))&
			)
		enddo
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL PRIVATE(k) 
!$OMP DO 
		do i=1,nx
			k=itopo(i)+2
			uz(k,i)=uz(k,i)+(dt/dz)*pa(k,i)*&
			(&
			a40*(u(k,i)-u(k-1,i))+&
			a41*(u(k+1,i)-u(k-2,i))&
			)
		enddo
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL PRIVATE(k) 
!$OMP DO 
		do i=1,nx
			k=itopo(i)+1
			uz(k,i)=uz(k,i)+(dt/dz)*pa(k,i)*&
			(&
			a20*(u(k,i)-u(k-1,i))&
			)
		enddo
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL PRIVATE(k) 
!$OMP DO 
		do i=1,nx
			k=itopo(i)
			uz(k,i)=uz(k,i)+(dt/dz)*pa(k,i)*&
			(&
			a20*(u(k,i)-u(k-1,i))&
			)
		enddo
!$OMP END DO
!$OMP END PARALLEL

	else if(zflag==1)then

!$OMP PARALLEL PRIVATE(k) 
!$OMP DO 
		do i=1,nx
			do k=itopo(i)+4,nz
				uz(k,i)=uz(k,i)+(dt/dz)*pa(k,i)*&
				(&
				a100*(u(k+1,i)-u(k,i))+&
				a101*(u(k+2,i)-u(k-1,i))+&
				a102*(u(k+3,i)-u(k-2,i))+&
				a103*(u(k+4,i)-u(k-3,i))+&
				a104*(u(k+5,i)-u(k-4,i))&
				)
			enddo
		enddo

!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL PRIVATE(k) 
!$OMP DO 
		do i=1,nx
			k=itopo(i)+3
			uz(k,i)=uz(k,i)+(dt/dz)*pa(k,i)*&
			(&
			a80*(u(k+1,i)-u(k,i))+&
			a81*(u(k+2,i)-u(k-1,i))+&
			a82*(u(k+3,i)-u(k-2,i))+&
			a83*(u(k+4,i)-u(k-3,i))&
			)
		enddo
!$OMP END DO
!$OMP END PARALLEL


!$OMP PARALLEL PRIVATE(k) 
!$OMP DO 
		do i=1,nx
			k=itopo(i)+2
			uz(k,i)=uz(k,i)+(dt/dz)*pa(k,i)*&
			(&
			a60*(u(k+1,i)-u(k,i))+&
			a61*(u(k+2,i)-u(k-1,i))+&
			a62*(u(k+3,i)-u(k-2,i))&
			)
		enddo
!$OMP END DO
!$OMP END PARALLEL


!$OMP PARALLEL PRIVATE(k) 
!$OMP DO 
		do i=1,nx
			k=itopo(i)+1
			uz(k,i)=uz(k,i)+(dt/dz)*pa(k,i)*&
			(&
			a40*(u(k+1,i)-u(k,i))+&
			a41*(u(k+2,i)-u(k-1,i))&
			)
		enddo
!$OMP END DO
!$OMP END PARALLEL


!$OMP PARALLEL PRIVATE(k) 
!$OMP DO 
		do i=1,nx
			k=itopo(i)
			uz(k,i)=uz(k,i)+(dt/dz)*pa(k,i)*&
			(&
			a20*(u(k+1,i)-u(k,i))&
			)
		enddo
!$OMP END DO
!$OMP END PARALLEL

	else

		write(*,*)'Wrong zflag value,please check'
		stop

	endif

	return
	end subroutine caluz


!**************************************************************
	subroutine caluzabs(zflag,nx,nz,npzn,dz,dt,u,uz,pa,itopo)
	use constant
	implicit none
	!Dummy variables
	integer::zflag,nx,nz,npzn,&
			itopo(nx)
	real::	dz,dt,&
			u(-4:nz+5,-4:nx+5),&
    		uz(0:nz,nx),pa(nz,nx)

	!Local variables
	integer::k,i

	if(zflag==0)then

!$OMP PARALLEL PRIVATE(k) 
!$OMP DO 
		do i=1,nx
			do k=itopo(i)-npzn,nz
				uz(k,i)=uz(k,i)+(dt/dz)*pa(k,i)*&
				(&
				a100*(u(k,i)-u(k-1,i))+&
				a101*(u(k+1,i)-u(k-2,i))+&
				a102*(u(k+2,i)-u(k-3,i))+&
				a103*(u(k+3,i)-u(k-4,i))+&
				a104*(u(k+4,i)-u(k-5,i))&
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
				uz(k,i)=uz(k,i)+(dt/dz)*pa(k,i)*&
				(&
				a100*(u(k+1,i)-u(k,i))+&
				a101*(u(k+2,i)-u(k-1,i))+&
				a102*(u(k+3,i)-u(k-2,i))+&
				a103*(u(k+4,i)-u(k-3,i))+&
				a104*(u(k+5,i)-u(k-4,i))&
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
	end subroutine caluzabs

!*********************************************************************************
	subroutine  calu(nx,nz,ux,uz,u)
	implicit none
	!Dummy variables
	integer::nx,nz
	real::	u(-4:nz+5,-4:nx+5),&
    		ux(0:nz,nx),uz(0:nz,nx)
	!Local variables
	integer::k,i
	
	do i=1,nx
		do k=1,nz
			u(k,i)=ux(k,i)+uz(k,i)
		enddo
	enddo

	return
	end subroutine calu

!***********************************************************************************
	subroutine addsrc(nx,nz,isxe,isze,it,idecay,rgau,f0,dx,dz,dt,gc2d,u)
	use constant

	implicit none
	!Dummy variables
	integer::nx,nz,isxe,isze,it,idecay,rgau
	real::f0,dx,dz,dt
	real::u(0:nz,nx)
	real::gc2d(2*rgau-1,2*rgau-1)

	!Local variables
	integer::i,k,hw1,hw2,idx1,idx2,m,n,ierr
	integer::n1,n2
	real::wavelet,atten

	idx1=1
	idx2=1

	wavelet=exp(-(pi*f0*(it*dt-idecay*dt))**2)*&
			(1-2*(pi*f0*(it*dt-idecay*dt))**2)


!	write(*,*)rgau

	do n=isxe-(rgau-1),isxe+(rgau-1)
		idx1=1
		do m=isze-(rgau-1),isze+(rgau-1)
			if(n>=1.and.n<=nx.and.m>=1.and.m<=nz&
				.and.idx2<=2*rgau-1.and.idx1<=2*rgau-1)then
				u(m,n)=wavelet*gc2d(idx1,idx2)*dt/(dz**2)+u(m,n)
			endif
			idx1=idx1+1 
		enddo
		idx2=idx2+1
	enddo

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

            write(8,rec=i)su_head,(record(k,i),k=1,nt)!**zy
            write(9,rec=i)su_head,(record_acc(k,i),k=1,nt)!**zy


        enddo
        close(unit=8)
		close(unit=9)



        return
        end subroutine write_currshot_disk

!****************************************************************************************
        subroutine merge_shot_files(nshot,shot_fn1,shotall_fn1,&
					shot_fn2,shotall_fn2,ofsmin,ofsmax,dx,nt)

         implicit none
          !Dummy variables
         character(len=256)::shot_fn1,shotall_fn1,shot_fn2,shotall_fn2
         integer::nshot,nt
		 real::ofsmin,ofsmax,dx
         !Local variables
         real,allocatable::record(:,:)
         integer::currshot_ran_x,&
                   is,i,j,k,ierr
         character(len=256)::currshot_noc
         character(len=256)::currshot_name
  
        !Initiallization


		currshot_ran_x=nint(ofsmax/dx)-nint(ofsmin/dx)+1

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

        

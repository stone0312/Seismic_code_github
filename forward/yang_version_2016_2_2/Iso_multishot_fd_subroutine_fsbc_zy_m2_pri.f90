!*************************************************************************************
		subroutine spcalpa(csname,csno,rgau,ixmine,ixmaxe,isx0,isz0,irz0,&
					isz,isx,idecay,idtsnap,nz,npxn,npxp,npzn,npzp,nxe,nze,&
					nxea,nzea,ppxn,ppxp,ppzn,ppzp,sx0,sz0,rz0,dx,dz,dsx,&
					ofsmin,ofsmax,decay,dt,dtsnap)


		use constant
		implicit none
		!Dummy variables
		character(len=256)::csname

		integer::csno,ixmine,ixmaxe,isx0,isz0,irz0,&
				isx,isz,idecay,idtsnap,&
				nz,npxn,npxp,npzn,npzp,&
				nxe,nze,nxea,nzea,rgau

		real::ppxn,ppxp,ppzn,ppzp,sx0,sz0,rz0,&
			dx,dz,dsx,ofsmin,ofsmax,decay,dt,&
			dtsnap



		!Computing the number of
		!grids with pml abc
        npxn=nint(ppxn/dx)+1
        npxp=nint(ppxp/dx)+1
        npzn=nint(ppzn/dz)+1
        npzp=nint(ppzp/dz)+1

		!Computing idtsnap,idecay
		idtsnap=nint(dtsnap/dt)
		idecay=nint(decay/dt)

        !Computing currshot range
		isx0=nint(sx0/dx)+1!**zy                                    
        isz0=nint(sz0/dz)+1!**zy                                    
        irz0=nint(rz0/dz)+1!**zy

		isx=isx0+nint((csno-1)*dsx/dx)
		isz=isz0
		ixmine=isx+nint(ofsmin/dx)
		ixmaxe=isx+nint(ofsmax/dx)

        nxe=(ixmaxe-ixmine)+1
        nxea=nxe+(npxn+npxp)
		nze=nz
        nzea=nz+(npzn+npzp)

		write(csname,'(I4)')csno

		return
		end subroutine spcalpa


	!**************************************************************************************************
	subroutine isoewm2dfsbc(csname,snapx_fn,snapz_fn,nt,idecay,idtsnap,&
				record1,record2,myid,irz0,rgau,gc2d,isxe,isze,nxe,nze,&
				nx,nz,npxn,npxp,npzn,npzp,f0,dx,dz,dt,vx,vz,strxx,&
				strzz,strxz,vxx,vxz,vzx,vzz,strxxx,strxxz,strzzx,strzzz,&
				strxzx,strxzz,der1,der2,spgxn,spgxp,spgzp,spg1xn,spg2xp,&
				spg2zp,vpe,vse,rhoe,itopoe,lamb,lambsg,lambmusg,mu,musg,&
				musgef,busg1,busg2,busg1ef,busg2ef)


	use constant
	implicit none
	!Dummy variables

	character(len=256)::snapx_fn,snapz_fn,csname
	integer::nt,isxe,isze,nxe,nze,nx,nz,&
			npxn,npxp,npzn,npzp,idecay,idtsnap,&
			myid,irz0,rgau
    integer::itopoe(nx)
	real::f0,dx,dz,dt
	real::&
		gc2d(2*rgau-1,2*rgau-1),&
		record1(nt,nxe),&
		record2(nt,nxe),&
		vx(-4:nz+5,-4:nx+5),&
 		vz(-4:nz+5,-4:nx+5),&
        strxx(-4:nz+5,-4:nx+5),&
        strzz(-4:nz+5,-4:nx+5),&
        strxz(-4:nz+5,-4:nx+5),&
        vxx(0:nz,nx),&
        vxz(0:nz,nx),&
        vzx(0:nz,nx),&
        vzz(0:nz,nx),&
        strxxx(0:nz,nx),&
        strxxz(0:nz,nx),&
        strzzx(0:nz,nx),&
        strzzz(0:nz,nx),&
        strxzx(0:nz,nx),&
        strxzz(0:nz,nx),&
        der1(nz,nx),&
        der2(nz,nx),&
		spgxn(npxn+1),&
		spgxp(npxp+1),&
		spgzp(npzp+1),&
		spg1xn(npxn+1),&
		spg2xp(npxp+1),&
		spg2zp(npzp+1),&
        vpe(nz,nx),&
        vse(nz,nx),&
        rhoe(nz,nx),&
		lamb(nz,nx),&
		lambsg(nz,nx),&
		lambmusg(nz,nx),&
		mu(nz,nx),&
		musg(nz,nx),&
		musgef(nz,nx),&
		busg1(nz,nx),&
		busg2(nz,nx),&
		busg1ef(nz,nx),&
		busg2ef(nz,nx)

		!Local variables
		integer::it,iit,xflag,zflag
	

		do it=1,nt+idecay

			!******************Updating velocity component******************
			!Updating vx
			xflag=0
			call calux(xflag,nx,nz,dx,dt,strxx,vxx,busg2ef,itopoe)
			zflag=0
			call caluz(zflag,nx,nz,dz,dt,strxz,vxz,busg2ef,itopoe)
			call calpmlbclrb(npxn,npxp,npzn,npzp,nxe,nze,nx,nz,&
					dx,dz,itopoe,spgxn,spgxp,spgzp,vxx,vxz)
	
			call calu(nx,nz,vxx,vxz,vx)
			!Updating vz
			xflag=1
			call calux(xflag,nx,nz,dx,dt,strxz,vzx,busg1ef,itopoe)
			zflag=1
			call caluz(zflag,nx,nz,dz,dt,strzz,vzz,busg1ef,itopoe)
			call calpmlbclrb(npxn,npxp,npzn,npzp,nxe,nze,nx,nz,&
					dx,dz,itopoe,spg1xn,spg2xp,spg2zp,vzx,vzz)
	
			call calu(nx,nz,vzx,vzz,vz)
	

			!Implementing free-surface boundary condition on 
			!velocity components
			call calfsbcv(vxx,vzx,lambsg,musgef,lambmusg,npxn,&
						nxe,nze,nx,nz,itopoe)

			call calfsbcv(vxz,vzz,lambsg,musgef,lambmusg,npxn,&
						nxe,nze,nx,nz,itopoe)


			!******************Updating stress component******************
	
			!Updating strxx
			xflag=1
			call calux(xflag,nx,nz,dx,dt,vx,strxxx,lambmusg,itopoe)
			zflag=0
			call caluz(zflag,nx,nz,dz,dt,vz,strxxz,lambsg,itopoe)
			call calpmlbclrb(npxn,npxp,npzn,npzp,nxe,nze,nx,nz,&
					dx,dz,itopoe,spg1xn,spg2xp,spgzp,strxxx,strxxz)
	
	
			!if source  at free surface,don't increment strxxx
			if(isze>itopoe(isxe))then
	!			write(*,*)here
				call addsrc(nx,nz,isxe,isze,it,idecay,rgau,f0,dx,dz,dt,gc2d,strxxx)
			endif
			call calu(nx,nz,strxxx,strxxz,strxx)

			!Updating strzz
			xflag=1
			call calux(xflag,nx,nz,dx,dt,vx,strzzx,lambsg,itopoe)
			zflag=0
			call caluz(zflag,nx,nz,dz,dt,vz,strzzz,lambmusg,itopoe)
			call calpmlbclrb(npxn,npxp,npzn,npzp,nxe,nze,nx,nz,&
					dx,dz,itopoe,spg1xn,spg2xp,spgzp,strzzx,strzzz)
	
			call calfsbcstrzz(strzzz,npxn,nxe,nze,nx,nz,itopoe)
			call calfsbcstrzz(strzzx,npxn,nxe,nze,nx,nz,itopoe)
	
			call addsrc(nx,nz,isxe,isze,it,idecay,rgau,f0,dx,dz,dt,gc2d,strzzx)
			call calu(nx,nz,strzzx,strzzz,strzz)
	
			!Updating strxz
			xflag=0
			call calux(xflag,nx,nz,dx,dt,vz,strxzx,musgef,itopoe)
			zflag=1
			call caluz(zflag,nx,nz,dz,dt,vx,strxzz,musgef,itopoe)
			call calpmlbclrb(npxn,npxp,npzn,npzp,nxe,nze,nx,nz,&
					dx,dz,itopoe,spgxn,spgxp,spg2zp,strxzx,strxzz)
	
			call calfsbcstrxz(strxzx,npxn,nxe,nze,nx,nz,itopoe)
			call calfsbcstrxz(strxzz,npxn,nxe,nze,nx,nz,itopoe)
	
			call calu(nx,nz,strxzx,strxzz,strxz)

			!For common use
			iit=it-idecay
			if(iit>=1)then
			!For comparison with analytical programs
!			iit=it
!			if(iit>=1.and.iit<=nt)then

				call get_record(iit,npxn,irz0,nxe,nx,nz,nt,&
								itopoe,vx,record1)

				call get_record(iit,npxn,irz0,nxe,nx,nz,nt,&
								itopoe,vz,record2)
			endif

			if(myid==0)then
				if(modulo(it,idtsnap)==0)then		
					write(*,*)'shot is',trim(csname),'myid=',myid,&
					'nt=',nt,'it=',it

					call write_currtsnap(snapx_fn,it,nxe,nze,nx,nz,&
						npxn,npzn,vx)

					call write_currtsnap(snapz_fn,it,nxe,nze,nx,nz,&
						npxn,npzn,vz)

				endif
			endif

		enddo

  	return

	end subroutine isoewm2dfsbc


!***************************************************************************************************
	subroutine isoewm2dabc(csname,snapx_fn,snapz_fn,nt,idecay,idtsnap,&
				record1,record2,myid,irz0,rgau,gc2d,isxe,isze,nxe,nze,&
				nx,nz,npxn,npxp,npzn,npzp,f0,dx,dz,dt,vx,vz,strxx,strzz,&
				strxz,vxx,vxz,vzx,vzz,strxxx,strxxz,strzzx,strzzz,&
				strxzx,strxzz,der1,der2,spgxn,spgxp,spgzn,spgzp,spg1xn,spg2xp,&
				spg1zn,spg2zp,vpe,vse,rhoe,itopoe,lamb,lambsg,lambmusg,mu,musg,&
				musgef,busg1,busg2,busg1ef,busg2ef)


	use constant
	implicit none
	!Dummy variables
	character(len=256)::snapx_fn,snapz_fn,csname
	integer::nt,isxe,isze,nxe,nze,nx,nz,&
			npxn,npxp,npzn,npzp,idecay,idtsnap,&
			myid,irz0,rgau

    integer::itopoe(nx)
	real::f0,dx,dz,dt
	real::&
		gc2d(2*rgau-1,2*rgau-1),&
		record1(nt,nxe),&
		record2(nt,nxe),&
		vx(-4:nz+5,-4:nx+5),&
 		vz(-4:nz+5,-4:nx+5),&
        strxx(-4:nz+5,-4:nx+5),&
        strzz(-4:nz+5,-4:nx+5),&
        strxz(-4:nz+5,-4:nx+5),&
        vxx(0:nz,nx),&
        vxz(0:nz,nx),&
        vzx(0:nz,nx),&
        vzz(0:nz,nx),&
        strxxx(0:nz,nx),&
        strxxz(0:nz,nx),&
        strzzx(0:nz,nx),&
        strzzz(0:nz,nx),&
        strxzx(0:nz,nx),&
        strxzz(0:nz,nx),&
        der1(nz,nx),&
        der2(nz,nx),&
		spgxn(npxn+1),&
		spgxp(npxp+1),&
		spgzn(npzn+1),&
		spgzp(npzp+1),&
		spg1xn(npxn+1),&
		spg2xp(npxp+1),&
		spg1zn(npzn+1),&
		spg2zp(npzp+1),&
        vpe(nz,nx),&
        vse(nz,nx),&
        rhoe(nz,nx),&
		lamb(nz,nx),&
		lambsg(nz,nx),&
		lambmusg(nz,nx),&
		mu(nz,nx),&
		musg(nz,nx),&
		musgef(nz,nx),&
		busg1(nz,nx),&
		busg2(nz,nx),&
		busg1ef(nz,nx),&
		busg2ef(nz,nx)


		!Local variables
		integer::it,iit,xflag,zflag


		call modipapml(nx,nz,rgau,itopoe,gc2d,vpe,vse,rhoe,lamb,lambsg,&
				lambmusg,mu,musg,musgef,busg1,busg2,busg1ef,busg2ef)


		do it=1,nt+idecay

			!******************Updating velocity component******************
	
			!Updating vx
			xflag=0
			call caluxabs(xflag,nx,nz,npzn,dx,dt,strxx,vxx,busg2ef,itopoe)
			zflag=0
			call caluzabs(zflag,nx,nz,npzn,dz,dt,strxz,vxz,busg2ef,itopoe)
			call calpmlbclrtb(npxn,npxp,npzn,npzp,nx,nz,dx,dz,&
						itopoe,spgxn,spgxp,spgzn,spgzp,vxx,vxz)
			call calu(nx,nz,vxx,vxz,vx)
	
	
			!Updating vz
			xflag=1
			call caluxabs(xflag,nx,nz,npzn,dx,dt,strxz,vzx,busg1ef,itopoe)
			zflag=1
			call caluzabs(zflag,nx,nz,npzn,dz,dt,strzz,vzz,busg1ef,itopoe)
			call calpmlbclrtb(npxn,npxp,npzn,npzp,nx,nz,dx,dz,&
					itopoe,spg1xn,spg2xp,spg1zn,spg2zp,vzx,vzz)
			call calu(nx,nz,vzx,vzz,vz)
	
			!******************Updating stress component******************
	
			!Updating strxx
			xflag=1
			call caluxabs(xflag,nx,nz,npzn,dx,dt,vx,strxxx,lambmusg,itopoe)
			zflag=0
			call caluzabs(zflag,nx,nz,npzn,dz,dt,vz,strxxz,lambsg,itopoe)
			call calpmlbclrtb(npxn,npxp,npzn,npzp,nx,nz,dx,dz,&
				itopoe,spg1xn,spg2xp,spgzn,spgzp,strxxx,strxxz)
			call addsrc(nx,nz,isxe,isze,it,idecay,rgau,f0,dx,dz,dt,gc2d,strxxx)
			call calu(nx,nz,strxxx,strxxz,strxx)
	
	
			!Updating strzz
			xflag=1
			call caluxabs(xflag,nx,nz,npzn,dx,dt,vx,strzzx,lambsg,itopoe)
			zflag=0
			call caluzabs(zflag,nx,nz,npzn,dz,dt,vz,strzzz,lambmusg,itopoe)
			call calpmlbclrtb(npxn,npxp,npzn,npzp,nx,nz,dx,dz,&
				itopoe,spg1xn,spg2xp,spgzn,spgzp,strzzx,strzzz)
			call addsrc(nx,nz,isxe,isze,it,idecay,rgau,f0,dx,dz,dt,gc2d,strzzx)
			call calu(nx,nz,strzzx,strzzz,strzz)
	
	
			!Updating strxz
			xflag=0
			call caluxabs(xflag,nx,nz,npzn,dx,dt,vz,strxzx,musgef,itopoe)
			zflag=1
			call caluzabs(zflag,nx,nz,npzn,dz,dt,vx,strxzz,musgef,itopoe)
			call calpmlbclrtb(npxn,npxp,npzn,npzp,nx,nz,dx,dz,&
				itopoe,spgxn,spgxp,spg1zn,spg2zp,strxzx,strxzz)
			call calu(nx,nz,strxzx,strxzz,strxz)

			iit=it-idecay
			if(iit>=1)then
				call get_record(iit,npxn,irz0,nxe,nx,nz,nt,&
								itopoe,vx,record1)
				call get_record(iit,npxn,irz0,nxe,nx,nz,nt,&
								itopoe,vz,record2)
			endif

			if(myid==0)then
				if(modulo(it,idtsnap)==0)then		
					write(*,*)'shot is',trim(csname),'myid=',myid,&
					'nt=',nt,'it=',it

					call write_currtsnap(snapx_fn,it,nxe,nze,nx,nz,&
						npxn,npzn,vx)

					call write_currtsnap(snapz_fn,it,nxe,nze,nx,nz,&
						npxn,npzn,vz)

				endif
			endif

		enddo

	return
	end subroutine isoewm2dabc



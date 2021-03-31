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
        end module header_module



		module allocbuff
			contains

			subroutine mallocbuff(err,nx,nz,topo,vp,vs,rho)
			implicit none
			!Dummy variables
			integer,intent(in)::nx,nz
			integer,intent(inout)::err
			real,allocatable::topo(:),rho(:,:),vp(:,:),vs(:,:)

			!Allocate memories for buffers
        	allocate(vp(nz,nx),STAT=err)
        	allocate(vs(nz,nx),STAT=err)
        	allocate(rho(nz,nx),STAT=err)
        	allocate(topo(nx),STAT=err)

        	!Initiallization
        	vp=0.0
        	vs=0.0
        	rho=0.0
        	topo=0.0

			return
			end subroutine mallocbuff


			subroutine sallocbuff(err,nt,nxe,nxea,nzea,npxn,npxp,npzn,npzp,rgau,&
						itopoe,vx,vz,strxx,strzz,strxz,vxx,vxz,vzx,vzz,strxxx,&
						strxxz,strzzx,strzzz,strxzx,strxzz,der1,der2,vpe,vse,&
						rhoe,record1,record2,topoe,mu,lamb,musg,musgef,lambsg,&
						lambmusg,busg1,busg2,busg1ef,busg2ef,spgxn,spgxp,spgzn,&
						spgzp,spg1xn,spg2xp,spg1zn,spg2zp,gc2d)

			implicit none
			!Dummy variables
			integer,intent(in)::nt,nxe,nxea,nzea,npxn,npxp,npzn,npzp,&
								rgau
			integer,intent(inout)::err
	        integer,allocatable::itopoe(:)
	        real,allocatable::vx(:,:),vz(:,:),strxx(:,:),&
						strzz(:,:),strxz(:,:),vxx(:,:),&
						vxz(:,:),vzx(:,:),vzz(:,:),&
						strxxx(:,:),strxxz(:,:),strzzx(:,:),&
	        			strzzz(:,:),strxzx(:,:),strxzz(:,:),&
						der1(:,:),der2(:,:),vpe(:,:),&
						vse(:,:),rhoe(:,:),record1(:,:),&
	        			record2(:,:),topoe(:),mu(:,:),lamb(:,:),&
						musg(:,:),musgef(:,:),&
						lambsg(:,:),lambmusg(:,:),busg1(:,:),&
	 					busg2(:,:),busg1ef(:,:),busg2ef(:,:),&
						spgxn(:),spgxp(:),spgzn(:),&
						spgzp(:),spg1xn(:),spg2xp(:),&
						spg1zn(:),spg2zp(:),gc2d(:,:)


			!Allocate memories for buffers
	        allocate(vx(-4:nzea+5,-4:nxea+5),STAT=err)
	        allocate(vz(-4:nzea+5,-4:nxea+5),STAT=err)
	        allocate(strxx(-4:nzea+5,-4:nxea+5),STAT=err)
	        allocate(strzz(-4:nzea+5,-4:nxea+5),STAT=err)
	        allocate(strxz(-4:nzea+5,-4:nxea+5),STAT=err)
	        allocate(vxx(0:nzea,nxea),STAT=err)
	        allocate(vxz(0:nzea,nxea),STAT=err)
	        allocate(vzx(0:nzea,nxea),STAT=err)
	        allocate(vzz(0:nzea,nxea),STAT=err)
	        allocate(strxxx(0:nzea,nxea),STAT=err)
	        allocate(strxxz(0:nzea,nxea),STAT=err)
	        allocate(strzzx(0:nzea,nxea),STAT=err)
	        allocate(strzzz(0:nzea,nxea),STAT=err)
	        allocate(strxzx(0:nzea,nxea),STAT=err)
	        allocate(strxzz(0:nzea,nxea),STAT=err)
	        allocate(der1(nzea,nxea),STAT=err)
	        allocate(der2(nzea,nxea),STAT=err)
	        allocate(vpe(nzea,nxea),STAT=err)
	        allocate(vse(nzea,nxea),STAT=err)
	        allocate(rhoe(nzea,nxea),STAT=err)
	        allocate(itopoe(nxea),STAT=err)
			allocate(mu(nzea,nxea),STAT=err)
			allocate(lamb(nzea,nxea),STAT=err)
			allocate(musg(nzea,nxea),STAT=err)
			allocate(musgef(nzea,nxea),STAT=err)
			allocate(lambsg(nzea,nxea),STAT=err)
			allocate(lambmusg(nzea,nxea),STAT=err)
			allocate(busg1(nzea,nxea),STAT=err)
			allocate(busg2(nzea,nxea),STAT=err)
			allocate(busg1ef(nzea,nxea),STAT=err)
			allocate(busg2ef(nzea,nxea),STAT=err)
			allocate(record1(nt,nxe),STAT=err)
			allocate(record2(nt,nxe),STAT=err)
			allocate(topoe(nxea),STAT=err)
			allocate(spgxn(npxn+1),STAT=err)
			allocate(spgxp(npxp+1),STAT=err)
			allocate(spgzn(npzn+1),STAT=err)
			allocate(spgzp(npzp+1),STAT=err)
			allocate(spg1xn(npxn+1),STAT=err)
			allocate(spg2xp(npxp+1),STAT=err)
			allocate(spg1zn(npzn+1),STAT=err)
			allocate(spg2zp(npzp+1),STAT=err)
			allocate(gc2d(2*rgau-1,2*rgau-1),STAT=err)

			!Zeros the buffers
	        vx=0.0
	        vz=0.0
	        strxx=0.0
	        strzz=0.0
	        strxz=0.0
	        vxx=0.0
	        vxz=0.0
	        vzx=0.0
	        vzz=0.0
	        strxxx=0.0
	        strxxz=0.0
	        strzzx=0.0
	        strzzz=0.0
	        strxzx=0.0
	        strxzz=0.0
	        der1=0.0
	        der2=0.0
	        vpe=0.0
	        vse=0.0
	        rhoe=0.0
	        itopoe=0
			mu=0.0
			lamb=0.0
			musg=0.0
			musgef=0.0
			lambsg=0.0
			lambmusg=0.0
			busg1=0.0
			busg2=0.0
			busg1ef=0.0
			busg2ef=0.0
			record1=0.0
			record2=0.0
			topoe=0.0
			spgxn=0.0
			spgxp=0.0
			spgzn=0.0
			spgzp=0.0
			spg1xn=0.0
			spg2xp=0.0
			spg1zn=0.0
			spg2zp=0.0
			gc2d=0.0

			return
			end subroutine  sallocbuff

		end module allocbuff


		module deallocbuff
			contains

			subroutine mdeallocbuff(err,nx,nz,topo,vp,vs,rho)
			implicit none
			!Dummy variables
			integer,intent(in)::nx,nz
			integer,intent(inout)::err
			real,allocatable::topo(:),rho(:,:),vp(:,:),vs(:,:)
			!Dellocate memories for buffers
        	deallocate(vp,STAT=err)
        	deallocate(vs,STAT=err)
        	deallocate(rho,STAT=err)
        	deallocate(topo,STAT=err)

			return
			end subroutine mdeallocbuff


			subroutine sdeallocbuff(err,itopoe,&
	        			vx,vz,strxx,strzz,strxz,vxx,&
						vxz,vzx,vzz,strxxx,strxxz,strzzx,&
	        			strzzz,strxzx,strxzz,der1,der2,vpe,&
						vse,rhoe,record1,record2,topoe,mu,&
						lamb,musg,musgef,lambsg,lambmusg,busg1,&
	 					busg2,busg1ef,busg2ef,spgxn,spgxp,spgzn,&
						spgzp,spg1xn,spg2xp,spg1zn,spg2zp,gc2d)

			implicit none
			!Dummy variables
			integer,intent(inout)::err
			integer,allocatable::itopoe(:)
	        real,allocatable::vx(:,:),vz(:,:),strxx(:,:),&
						strzz(:,:),strxz(:,:),vxx(:,:),&
						vxz(:,:),vzx(:,:),vzz(:,:),&
						strxxx(:,:),strxxz(:,:),strzzx(:,:),&
	        			strzzz(:,:),strxzx(:,:),strxzz(:,:),&
						der1(:,:),der2(:,:),vpe(:,:),&
						vse(:,:),rhoe(:,:),record1(:,:),&
	        			record2(:,:),topoe(:),mu(:,:),lamb(:,:),&
						musg(:,:),musgef(:,:),&
						lambsg(:,:),lambmusg(:,:),busg1(:,:),&
	 					busg2(:,:),busg1ef(:,:),busg2ef(:,:),&
						spgxn(:),spgxp(:),spgzn(:),&
						spgzp(:),spg1xn(:),spg2xp(:),&
						spg1zn(:),spg2zp(:),gc2d(:,:)

	        deallocate(vx,STAT=err)
	        deallocate(vz,STAT=err)
	        deallocate(strxx,STAT=err)
	        deallocate(strzz,STAT=err)
	        deallocate(strxz,STAT=err)
	        deallocate(vxx,STAT=err)
	        deallocate(vxz,STAT=err)
	        deallocate(vzx,STAT=err)
	        deallocate(vzz,STAT=err)
	        deallocate(strxxx,STAT=err)
	        deallocate(strxxz,STAT=err)
	        deallocate(strzzx,STAT=err)
	        deallocate(strzzz,STAT=err)
	        deallocate(strxzx,STAT=err)
	        deallocate(strxzz,STAT=err)
	        deallocate(der1,STAT=err)
	        deallocate(der2,STAT=err)
	        deallocate(vpe,STAT=err)
	        deallocate(vse,STAT=err)
	        deallocate(rhoe,STAT=err)
	        deallocate(itopoe,STAT=err)
			deallocate(mu,STAT=err)
			deallocate(lamb,STAT=err)
			deallocate(musg,STAT=err)
			deallocate(musgef,STAT=err)
			deallocate(lambsg,STAT=err)
			deallocate(lambmusg,STAT=err)
			deallocate(busg1,STAT=err)
			deallocate(busg2,STAT=err)
			deallocate(busg1ef,STAT=err)
			deallocate(busg2ef,STAT=err)
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
			deallocate(gc2d,STAT=err)

			return
			end subroutine  sdeallocbuff

		end module deallocbuff











































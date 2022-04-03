#include "iostream.h"
#include "math.h"
#include "fstream.h"
#include "stdlib.h"
#include <complex.h>
#include "alloc.c"
#include "fftw3.h"
#define pai 3.14159265

int main()
{
   char fn1[256],fn2[256],fn3[256],fn4[256];
   int nx,ns,nr,lt,ltt,nsr,zmax,zmin,bou,nxb,wid,nxx;
   float dt,dsx,drx,fmin,fmax,dx,dz,v,dip_max;
   int itmin,itmax,shot,len;

   int icount,icount1,nt,abs1;
   float st_beg,st_end;   
   float k,k1,mir_k,theta,fai;


   int is,ir,it,iz,is1,ir1,is2,ir2,itmp,ix,tmp,tmp1;
   
    ifstream swq;
    swq.open("mirror_green_function_slope_single_survey.par");
    if(!swq)
      {
          cout<<"cannot open mirror_green_function_single_survey.par"<<endl;
          abort();
      }
    swq>>fn1>>fn2>>fn3>>fn4>>nx>>dx>>dz>>ns>>nr>>lt>>dt>>dsx>>drx>>fmin>>fmax>>v>>bou>>shot>>dip_max>>wid>>len>>nt>>st_beg>>st_end;
    swq.close();

    cout<<"fna of water layer depth is==== "<<fn1<<endl;
    cout<<"fna of real parts of green's function for positive frequency is==== "<<fn2<<endl;
    cout<<"fna of imaginary parts of green's function for positive frequency is==== "<<fn3<<endl;
    cout<<"fna of final green's function is==== "<<fn4<<endl;
    cout<<"No. of vel model width is==== "<<nx<<endl;
    cout<<"lateral interval of vel model is==== "<<dx<<endl;
    cout<<"vertical interval of vel model is==== "<<dz<<endl;
    cout<<"No. of shots is==== "<<ns<<endl;
    cout<<"No. of traces per shot is==== "<<nr<<endl;
    cout<<"No. of samples is==== "<<lt<<endl;
    cout<<"temperal interal (ms) is==== "<<dt<<endl;
    cout<<"spatial interval of sourses is==== "<<dsx<<endl;
    cout<<"spatial interval of receivers is==== "<<dsx<<endl;
    cout<<"minimum frequency is==== "<<fmin<<endl;
    cout<<"maximum frequency is==== "<<fmax<<endl;
    cout<<"velocity of water is==== "<<v<<endl;
    cout<<"half of width of boundary is==== "<<bou<<endl;
    cout<<"No. of shot to be predicted is===="<<shot<<endl;
    cout<<"maximum dip is===="<<dip_max<<endl;
    cout<<"width of taper for absorbing  is===="<<wid<<endl;
    cout<<"length of taper for gibbs effects of frequency is===="<<len<<endl;

    ltt=2*lt; 
    nxb=2*nx-1+2*bou;
    nsr=nx*nx;

    cout<<nxb<<endl;

//ordering the depth of water layer to find the maximum depth
    float *dep;
    dep=alloc1float(2*nx-1);

    float *dep1;
    dep1=alloc1float(2*nx);

    float *mir_dep;
    mir_dep=alloc1float(2*nx-1);
 
   int *receiver;
    receiver=alloc1int(2*nx-1);

   int *mir_receiver;
    mir_receiver=alloc1int(2*nx-1);

    ifstream swq1;
    swq1.open(fn1,ios::binary);
    if(!swq1)
       {
          cout<<"cannot open "<<fn1<<endl;
          abort();
       } 
    for(is=0;is<2*nx-1;is++)
       swq1.read((char*)&dep[is],sizeof(dep[is]));
    swq1.close();

    for(is=0;is<2*nx-1;is++)
       cout<<is+1<<"  "<<dep[is]<<endl;
//    return 0;


    for(is=0;is<2*nx-1;is++)
       dep1[is]=dep[is];
    dep1[2*nx-1]=dep[2*nx-2];

    complex<float> *s;
    s=alloc1complex(nxb);
    for(is=0;is<nxb;is++)
       s[is]=(0.0,0.0);

    complex<float> *s_k;
    s_k=alloc1complex(nxb);
 
    float *omega;
    omega=alloc1float(ltt);
    float *kx;
    kx=alloc1float(nxb);

    complex <float> *g_tmp;
    g_tmp=alloc1complex(nxb);

    complex <float> *g_k; 
    g_k=alloc1complex(nxb);

    complex <float> **g;
    g=alloc2complex(nx,nx);

    complex <float> *gf;
    gf=alloc1complex(ltt);

    float *gt;
    gt=alloc1float(ltt);

    complex<float> gtmp1;

//define the plans and arrays when fftw is complemented
    fftwf_complex *in1,*out1,*in2,*out2,*in3,*out3;
    fftwf_plan p1,p2,p3;
    in1=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * nxb);
    out1=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * nxb);
    in2=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * nxb);
    out2=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * nxb);
    in3=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ltt);
    out3=(fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ltt);
 
    p1=fftwf_plan_dft_1d(nxb,in1,out1,FFTW_FORWARD,FFTW_ESTIMATE);
    p2=fftwf_plan_dft_1d(nxb,in2,out2,FFTW_BACKWARD,FFTW_ESTIMATE);
    p3=fftwf_plan_dft_1d(ltt,in3,out3,FFTW_BACKWARD,FFTW_ESTIMATE);


//calculate the omega and kx
   for(it=0;it<ltt/2+1;it++)
      omega[it]=2*pai*it*1000/(dt*ltt);
   for(it=ltt/2+1;it<ltt;it++)
      omega[it]=2*pai*(-1000/(2*dt)+(it-ltt/2)*1000/(dt*ltt));

   itmin=int(fmin*dt*ltt/1000);
   itmax=int(fmax*dt*ltt/1000)+1;

   cout<<"totally "<<(itmax-itmin)<<" frequency slices need to be calculated"<<endl;

   float *ham;
   ham=alloc1float(itmax);
   for(it=0;it<itmax;it++)
      ham[it]=1.0;
   for(it=itmin;it<(itmin+len);it++)
      ham[it]=sqrt(sin(pai*(it-itmin+1)/(2*(len))));
   for(it=itmax-len;it<itmax;it++)
      ham[it]=ham[itmin+itmax-it-1];

   for(is=0;is<nxb/2+1;is++)
      kx[is]=2*pai*float(is)/float(dx*nxb); 
   for(is=nxb/2+1;is<nxb;is++)
      kx[is]=-2*pai*1.0/float(2*dx)+2*pai*float(is-nxb/2)/float(dx*nxb);
 
//read the original data

    ofstream swq2;
    swq2.open(fn2,ios::binary);
    if(!swq2)
       {
          cout<<"cannot open "<<fn2<<endl;
          abort();
       }
    ofstream swq3;
    swq3.open(fn3,ios::binary);
    if(!swq3)
       {
          cout<<"cannot open "<<fn3<<endl;
          abort();
       }

    cout<<"Prediction Starts..."<<endl;
	
    for(it=0;it<ltt/2+1;it++)    
        {
          if(it<itmin)
            {
               for(is=0;is<nx;is++)
                 for(ir=0;ir<nx;ir++) 
                 {
                    g[is][ir]=(0.0,0.0);  
                    swq2.write((char*)&g[is][ir].real(),sizeof(g[is][ir].real()));  
                    swq3.write((char*)&g[is][ir].imag(),sizeof(g[is][ir].imag()));                
                 }
               cout<<it<<"frequency slices have been finished!"<<endl;
            }      
          else if(it>itmax&&it<ltt/2+1)
            {
               for(is=0;is<nx;is++)
                 for(ir=0;ir<nx;ir++)
                 {
                    g[is][ir]=(0.0,0.0);
                    swq2.write((char*)&g[is][ir].real(),sizeof(g[is][ir].real()));
                    swq3.write((char*)&g[is][ir].imag(),sizeof(g[is][ir].imag())); 
                 }
               cout<<it<<"frequency slices have been finished!"<<endl;
            }


	  else
            {

//         cout<<"2222"<<endl; 
 
             if((it+1)%5==0)
                cout<<it<<"frequency slices have been finished!"<<endl;

             for(is=0;is<nx;is++)
               for(ir=0;ir<nx;ir++)
                 g[is][ir]=(0.0,0.0);

//calculating Green's functions of each source (or equally receiver) in kx-w domain;    
           for(is2=nx-1;is2<2*nx-1;is2++)
            {
     
//              cout<<it<<"   "<<is2<<"   2222"<<endl; 
     
              for(ix=0;ix<2*nx-1;ix++)
               {

                 k=(dep1[ix+1]-dep1[ix])/dx;
                 theta=atan(k);

                 k1=dep1[ix]/((ix-is2)*dx);
                 fai=atan(k1);

                 mir_k=tan(2*theta-fai);

                 if(-(float(dep1[ix])-float(mir_k*ix*dx))/mir_k/dx>=0.0)
                   receiver[ix]=int(-(float(dep1[ix])-float(mir_k*ix*dx))/mir_k/dx+0.5);
                 else
                   receiver[ix]=int(-(float(dep1[ix])-float(mir_k*ix*dx))/mir_k/dx-0.5);

                 mir_dep[ix]=fabs((sqrt(dep1[ix]*dep1[ix]+(ix*dx-receiver[ix]*dx)*(ix*dx-receiver[ix]*dx))+sqrt(dep1[ix]*dep1[ix]+(ix-is2)*dx*(ix-is2)*dx))*sin(fai));

//                 mir_dep[ix]=2*fabs((sqrt(dep1[ix]*dep1[ix]+(ix*dx-receiver[ix]*dx)*(ix*dx-receiver[ix]*dx)))*sin(theta-fai));

                if(cos(fai)>=0)
                  mir_receiver[ix]=is2+int(fabs((sqrt(dep1[ix]*dep1[ix]+(ix*dx-receiver[ix]*dx)*(ix*dx-receiver[ix]*dx))+sqrt(dep1[ix]*dep1[ix]+(ix-is2)*dx*(ix-is2)*dx))*cos(fai)/dx)+0.5);
                else
                  mir_receiver[ix]=is2+int(fabs((sqrt(dep1[ix]*dep1[ix]+(ix*dx-receiver[ix]*dx)*(ix*dx-receiver[ix]*dx))+sqrt(dep1[ix]*dep1[ix]+(ix-is2)*dx*(ix-is2)*dx))*cos(fai)/dx-0.5));


//                cout<<"Reflection Point==== "<<ix<<", Ishot==== "<<is2<<" ,Receiver==== "<<receiver[ix]<<" ,Mir_Depth==== "<<mir_dep[ix]<<" ,Mir_Receiver==== "<<mir_receiver[ix]<<"   2theta-fai===="<<2*theta-fai<<endl; 


               }

//            cout<<is2<<"   2222"<<endl; 

//              return 0;

              for(is=0;is<nxb;is++)
                 s[is]=(0.0,0.0);
              s[bou+is2].real()=1.0*ham[it];
              s[bou+is2].imag()=0.0;

              if((in1==NULL)||(out1==NULL))
                cout<<"memory insufficient"<<endl;
              else
                {
                  for(is=0;is<nxb;is++)
                    {
                       in1[is][0]=s[is].real();
                       in1[is][1]=s[is].imag();
                    }
                }

              fftwf_execute(p1);

              for(is=0;is<nxb;is++)
                {
                   s_k[is].real()=out1[is][0];
                   s_k[is].imag()=out1[is][1];
                }
         

//               cout<<"2222222===="<<endl; 


//forward extrapolation
               for(ir=0;ir<2*nx-1;ir++) 
                {

                 tmp=receiver[ir];
                 tmp1=bou+mir_receiver[ir];

//               cout<<ir<<"  tmp==="<<tmp<<endl; 


                 if(tmp>=is2-nx+1&&tmp<=is2)
                 {//1111

//                cout<<ir<<"====2222222===="<<endl; 

                  for(is=0;is<nxb;is++)
                   {
                       if(pow(omega[it]/v,2)-pow(kx[is],2)<0.0)
                        {
                           g_k[is].real()=0.0;
                           g_k[is].imag()=0.0;
                        }

                       else if(pow(omega[it]/v,2)-pow(kx[is],2)>pow(cos(2*pai*dip_max/360.0)*omega[it]/v,2))
                        {
                          gtmp1.real()=0.0;               //phase-shift extrapolating operator
                          gtmp1.imag()=-mir_dep[ir]*sqrt(pow(omega[it]/v,2)-pow(kx[is],2));
                          g_k[is]=s_k[is]*exp(gtmp1);
                        }
                       else
                        {
                           g_k[is].real()=0.0;
                           g_k[is].imag()=0.0;
                        }
                   }


//               cout<<ir<<"====2222222===="<<endl; 


//final Green's functions in kx-w domain for a certain (kx,w) pair;
                 if((in2==NULL)||(out2==NULL))
                 cout<<"memory insufficient"<<endl;
                 else
                   {   
                      for(is=0;is<nxb;is++)
                        {   
                           in2[is][0]=g_k[is].real();
                           in2[is][1]=g_k[is].imag();
                        }   
                   }   
                 fftwf_execute(p2);

                 for(is=0;is<nxb;is++)
                   {   
                      g_tmp[is].real()=out2[is][0]/nxb;
                      g_tmp[is].imag()=out2[is][1]/nxb;
                   }  

//                cout<<ir<<"====2222222===="<<is2-tmp<<"  "<<tmp1<<endl; 

                  if(tmp1>-1&&tmp1<nxb)
                    g[is2-nx+1][is2-tmp]+=g_tmp[tmp1];  
                  else
                    g[is2-nx+1][is2-tmp]+=(0.0,0.0);

//                cout<<ir<<"====2222222===="<<endl; 

                 }//1111  
                else
                  g[is2-nx+1][is2-tmp]+=(0.0,0.0);


                }   


//               cout<<"====2222===="<<endl; 


             }//end of is2 loop, .ie. Green's function of all shots for a certain frequency 
           
            for(is=0;is<nx;is++)
              for(ir=0;ir<nx;ir++)
                 {
                    swq2.write((char*)&g[is][ir].real(),sizeof(g[is][ir].real()));
                    swq3.write((char*)&g[is][ir].imag(),sizeof(g[is][ir].imag()));
                 } 

            }//end of else,.ie. a certain frequency
         }//end of all positive frequency

        swq2.close();
        swq3.close();

        cout<<"Positive frequcies prediction Done!"<<endl;

       ifstream swq22;
       swq22.open(fn2,ios::binary);
       if(!swq22)
         {
           cout<<"cannot open "<<fn2<<endl;
           abort();
         }

       ifstream swq33;
       swq33.open(fn3,ios::binary);
       if(!swq33)
         {
           cout<<"cannot open "<<fn3<<endl;
           abort();
         }

       ofstream swq4;
       swq4.open(fn4,ios::binary);
       if(!swq4)
         {
             cout<<"cannot open "<<fn4<<endl;
             abort();
         }

//start to compute the predicticted model of the negative frequency
        cout<<"Negative frequcies prediction Done!"<<endl;
        cout<<"IFFT Starts..."<<endl;
        
        for(is=0;is<nsr;is++)
          {
           swq22.seekg(is*4,ios::beg); 
           swq22.read((char*)&gf[0].real(),sizeof(gf[0].real())); 
           swq33.seekg(is*4,ios::beg);  
           swq33.read((char*)&gf[0].imag(),sizeof(gf[0].imag()));

           for(it=1;it<ltt/2+1;it++)
             {
               swq22.seekg((nsr-1)*4,ios::cur);
               swq33.seekg((nsr-1)*4,ios::cur);
               swq22.read((char*)&gf[it].real(),sizeof(gf[it].real()));
               swq33.read((char*)&gf[it].imag(),sizeof(gf[it].imag()));
             } 
           for(it=ltt/2+1;it<ltt;it++)
               {
                 gf[it].real()=gf[ltt-it].real();
                 gf[it].imag()=-gf[ltt-it].imag();
               } 
              
           if((in3==NULL)||(out3==NULL))
              cout<<"memory insufficient"<<endl;
           else
              {
                for(it=0;it<ltt;it++)
                  {
                     in3[it][0]=gf[it].real();
                     in3[it][1]=gf[it].imag();
                  }
              }

              fftwf_execute(p3);

              for(it=0;it<ltt;it++)
                 gt[it]=out3[it][0]/ltt;

              for(it=0;it<lt;it++)
                swq4.write((char*)&gt[it],sizeof(gt[it]));
   
          }

       
        swq22.close(); 
        swq33.close();
        swq4.close();

        cout<<"disk written."<<endl;

        cout<<"all done!"<<endl;

		
        return 0;
 
   
}







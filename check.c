#include "global_var.h"
#include "global_ex.h"
#include "ex_func.h"
#include "math.h"
#include <stdio.h>

void check(){
  // v_pp at 100 au
  double a_p,r=100,r_in=1.0,vr, d_v,v_Brown,d_vr,d_vt,d_vz,v_turb,tgrowth,St,hdust,alpha,m_peb;
  int i,Nsize=300;
  FILE *fp,*fp0,*fp1,*fp2,*fp3,*fp4,*fp5,*fp6,*fp7,*fp8;
  
  fp=fopen("size_check.txt","w");
  fp0=fopen("vpp_check.txt","w");
  fp1=fopen("vBrown_check.txt","w");
  fp2=fopen("dvr_check.txt","w");
  fp3=fopen("dvt_check.txt","w");
  fp4=fopen("dvz_check.txt","w");
  fp5=fopen("vturb_check.txt","w");
  fp6=fopen("tgrowth100au.txt","w");
  fp7=fopen("tgrowth_inner.txt","w");
  fp8=fopen("vr_inner.txt","w");




  for(i=0;i<Nsize;i++){
    a_p=1e-5*exp(i*1.0/Nsize*log(1e3/1e-5));
    d_v=v_pp(r,a_p,0);
    v_Brown=v_pp(r,a_p,1);
    d_vr=v_pp(r,a_p,2);
    d_vt=v_pp(r,a_p,3);
    d_vz=v_pp(r,a_p,4);
    v_turb=v_pp(r,a_p,5);
    m_peb=4.*M_PI*rho_peb*a_p*a_p*a_p/3.;
    St=stokes(100.,a_p);
    alpha=alpha_func(100.);
    hdust=height(100.)/sqrt(1+St*(1+2*St)\
        /alpha/(1+St));
    tgrowth=3*m_peb*hdust/2/sqrt(M_PI)/\
            a_p/a_p/Sigma_gas(100.)/dust_gas/d_v;

    fprintf(fp,"%e\n",a_p);
    fprintf(fp0,"%e\n",d_v);
    fprintf(fp1,"%e\n",v_Brown);
    fprintf(fp2,"%e\n",d_vr);
    fprintf(fp3,"%e\n",d_vt);
    fprintf(fp4,"%e\n",d_vz);
    fprintf(fp5,"%e\n",v_turb);
    fprintf(fp6,"%e\n",tgrowth);
    
    d_v=v_pp(r_in,a_p,0);
    St=stokes(r_in,a_p);
    alpha=alpha_func(r_in);
    hdust=height(r_in)/sqrt(1+St*(1+2*St)\
        /alpha/(1+St));
    tgrowth=3*m_peb*hdust/2/sqrt(M_PI)/\
            a_p/a_p/Sigma_gas(r_in)/dust_gas/d_v;
    vr=v_r(r_in,St);
    fprintf(fp7,"%e\n",tgrowth);
    fprintf(fp8,"%e\n",vr);


  }
  fclose(fp);
  fclose(fp0);
  fclose(fp1);
  fclose(fp2);
  fclose(fp3);
  fclose(fp4);
  fclose(fp5);
  fclose(fp6);
  fclose(fp7);
  fclose(fp8);


  

}

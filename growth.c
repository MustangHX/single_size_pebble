#include "global_var.h"
#include "global_ex.h"
#include "ex_func.h"
#include "math.h"
#include <stdio.h>
void upwind_size(double dt0){
  double dt=dt0*TUNIT;
  int i;
  double umax=0.,u,u1,u2,qx,r,r1,r2,a_p; //r1 is i-1,r2 is i+1
  double qtemp[ring_num]={0.};
  double srcfac=1.0,drtfac=1.;
  for(i=0;i<ring_num-1;i++){
    r=dust[i].r*LUNIT;
    r2=dust[i+1].r*LUNIT;
    if (fabs(v_r(dust[i].r,dust[i].a_p)/(r2-r))>umax){
      umax=fabs(v_r(dust[i].r,dust[i].a_p))/(r2-r);
    }
  }
  
  for(i=0;i<ring_num;i++){
    qtemp[i]=dust[i].m_peb;//*dust[i].sigma;
  }
  for(i=1;i<ring_num-1;i++){
    r=dust[i].r*LUNIT;
    r1=dust[i-1].r*LUNIT;
    r2=dust[i+1].r*LUNIT;
    a_p=dust[i].a_p;

    u=v_r(r/LUNIT,a_p);
    if (u>0.){
      qx=u*(qtemp[i]-qtemp[i-1])/(r-r1);
    }
    else{
      qx=u*(qtemp[i+1]-qtemp[i])/(r2-r);
    }
    dust[i].m_peb=(qtemp[i]-drtfac*dt*qx);
    //a_p=dust[i].a_p;
    //printf("a_p=%e\t m_peb=%e\n",dust[i].a_p,dust[i].m_peb);

    double St,alpha;
    a_p=cbrt(3.*dust[i].m_peb/4./M_PI/rho_peb);
    St=stokes(r/LUNIT,a_p);
    dust[i].h=disk[i].h/sqrt(1.+St*(1+2.*St)/alpha/(1.+St));

    //St=stokes(r/LUNIT,a_p);
    alpha=alpha_func(r/LUNIT);
    // growth souce term
    double d_v;
    /*double d_vr,d_vt;
    d_vr=v_r(r/LUNIT,St)-v_r(r/LUNIT,ratio_st*St);
    d_vt=v_t(r/LUNIT,St)-v_t(r/LUNIT,ratio_st*St);

    double d_v=sqrt(d_vr*d_vr+d_vt*d_vt);
    double v_Brown=sqrt(8*qtemp[i]*disk[i].temp*k_B/M_PI/qtemp[i]/qtemp[i]);
    double v_turb,Re_turb,D_turb;
    D_turb=alpha*disk[i].cs*disk[i].h;
    Re_turb=D_turb/disk[i].visc_mol;
    if(St<1/sqrt(Re_turb)) v_turb=sqrt(alpha)*pow(Re_turb,0.25)*St*(1-ratio_st);
    else v_turb=sqrt(3.*alpha)*disk[i].cs*sqrt(St);
    if(St>0.1) v_turb=0.;
    d_v=sqrt(d_v*d_v+v_Brown*v_Brown+v_turb*v_turb);*/
    d_v=v_pp(r/LUNIT,a_p,0);
    double rho_dust=dust[i].sigma/dust[i].h/sqrt(2.*M_PI);
    double srcterm=srcfac*2.*sqrt(M_PI)*a_p*a_p*d_v*dust[i].sigma/dust[i].h;//4*M_PI*a_p*a_p*d_v*rho_dust/srcfac;
    //printf("rhodust=%e\t src=%e\n",dust[i].a_p,d_vr);
    
    double grow_fac=1.0;
    if(log(v_frag/d_v)/log(5.)<1.) grow_fac=log(v_frag/d_v)/log(5.);
    //printf("grow_fac=%e\t, d_v=%e\n",grow_fac,d_v);
    dust[i].m_peb+=srcterm*dt;
    dust[i].a_p=cbrt(3*dust[i].m_peb/4/M_PI/rho_peb);
    dust[i].St=stokes(r/LUNIT,a_p);
    dust[i].h=disk[i].h/sqrt(1+St*(1+2*St)/alpha/(1+St)); 
    dust[i].d_v=d_v;
    dust[i].vr=v_r(dust[i].r,St);
  }
  // boundary
  //inner
  r=dust[0].r*LUNIT;
  r2=dust[1].r*LUNIT;
  a_p=dust[0].a_p;
  u=v_r(r/LUNIT,a_p);
  qx=u*(qtemp[1]-0.*qtemp[0])/(r2-r);
  double St,alpha,d_v,d_vr,d_vt,rho_dust,srcterm,grow_fac;
  St=stokes(r/LUNIT,a_p);
  alpha=alpha_func(r/LUNIT);
  dust[0].m_peb=(qtemp[0]-dt*qx);

  /*d_vr=v_r(r/LUNIT,St)-v_r(r/LUNIT,ratio_st*St);
  d_vt=v_t(r/LUNIT,St)-v_t(r/LUNIT,ratio_st*St);

  d_v=sqrt(d_vr*d_vr+d_vt*d_vt);
  double v_Brown=sqrt(8*qtemp[0]*disk[0].temp*k_B/M_PI/qtemp[0]/qtemp[0]);
  double v_turb,Re_turb,D_turb;
  D_turb=alpha*disk[0].cs*disk[0].h;
  Re_turb=D_turb/disk[0].visc_mol;
  if(St<1/sqrt(Re_turb)) v_turb=sqrt(alpha)*pow(Re_turb,0.25)*St*(1-ratio_st);
  else v_turb=sqrt(3.*alpha)*disk[0].cs*sqrt(St);
  if(St>0.1) v_turb=0.;
  d_v=sqrt(d_v*d_v+v_Brown*v_Brown+v_turb*v_turb);*/
    d_v=v_pp(r/LUNIT,a_p,0);
  rho_dust=dust[0].sigma/dust[0].h/sqrt(2.*M_PI);
  srcterm=srcfac*2.*sqrt(M_PI)*a_p*a_p*d_v*dust[i].sigma/dust[i].h;//4*M_PI*a_p*a_p*d_v*rho_dust/srcfac;

  grow_fac=1.0;
  if(log(v_frag/d_v)/log(5.)<1.) grow_fac=log(v_frag/d_v)/log(5.);
    //printf("grow_fac=%e\t, d_v=%e\n",grow_fac,d_v);
  dust[0].m_peb+=srcterm*dt;
  dust[0].a_p=cbrt(3*dust[0].m_peb/4/M_PI/rho_peb);
  dust[0].St=stokes(r/LUNIT,a_p);
  dust[0].h=disk[0].h/sqrt(1+St*(1+2*St)/alpha/(1+St));
  dust[0].d_v=d_v;
  dust[0].vr=v_r(dust[0].r,St);
  //
  //outer
  
  r=dust[ring_num-1].r*LUNIT;
  r1=dust[ring_num-2].r*LUNIT;
  a_p=dust[ring_num-1].a_p;
  u=v_r(r/LUNIT,a_p);
  qx=u*(0.-qtemp[ring_num-1])/(r1-r);
  St=stokes(r/LUNIT,a_p);
  alpha=alpha_func(r/LUNIT);
 /* 
  d_vr=v_r(r/LUNIT,St)-v_r(r/LUNIT,ratio_st*St);
  d_vt=v_t(r/LUNIT,St)-v_t(r/LUNIT,ratio_st*St);

  d_v=sqrt(d_vr*d_vr+d_vt*d_vt);
  v_Brown=sqrt(8*qtemp[ring_num-1]*disk[ring_num-1].temp*k_B/M_PI/qtemp[ring_num-1]/qtemp[ring_num-1]);
  D_turb=alpha*disk[ring_num-1].cs*disk[ring_num-1].h;
  Re_turb=D_turb/disk[ring_num-1].visc_mol;
  if(St<1/sqrt(Re_turb)) v_turb=sqrt(alpha)*pow(Re_turb,0.25)*St*(1-ratio_st);
  else v_turb=sqrt(3.*alpha)*disk[ring_num-1].cs*sqrt(St);
  if(St>0.1) v_turb=0.;
  d_v=sqrt(d_v*d_v+v_Brown*v_Brown+v_turb*v_turb);
  */
  d_v=v_pp(r/LUNIT,a_p,0);
  rho_dust=dust[ring_num-1].sigma/dust[ring_num-1].h/sqrt(2.*M_PI);
  srcterm=srcfac*2.*sqrt(M_PI)*a_p*a_p*d_v*dust[i].sigma/dust[i].h;//4*M_PI*a_p*a_p*d_v*rho_dust/srcfac;
  //dust[ring_num-1].m_peb=(qtemp[ring_num-1]-dt*qx);

  grow_fac=1.0;
  if(log(v_frag/d_v)/log(5.)<1.) grow_fac=log(v_frag/d_v)/log(5.);
    //printf("grow_fac=%e\t, d_v=%e\n",grow_fac,d_v);
  dust[ring_num-1].m_peb+=srcterm*dt;
  dust[ring_num-1].a_p=cbrt(3*dust[ring_num-1].m_peb/4/M_PI/rho_peb);
  dust[ring_num-1].St=stokes(r/LUNIT,a_p);
  dust[ring_num-1].h=disk[ring_num-1].h/sqrt(1+St*(1+2*St)/alpha/(1+St));
  dust[ring_num-1].d_v=d_v;

}
void grow_two_pop(double dt0, double tot_time){
  double dt=dt0*TUNIT;
  double tot_t=tot_time*TUNIT;
  double tau_grow, a_frag, a_drift, a_df, v_f,\
    cs, alpha,f_m, St, St_df, a0=a_min, a1, a_p, r;
  double srcterm, srcfac=1.0, rho_dust, d_v;
  int i;
  for (i=0;i<ring_num;i++){
    r=dust[i].r;
    alpha=alpha_func(r);
    cs=disk[i].cs;
    v_f=v_frag;
    tau_grow=1./(dust[i].sigma/disk[i].sigma*w_K(r));
    a_frag=FF*(2./3./M_PI)*(disk[i].sigma/rho_peb/alpha)*(v_f*v_f/cs/cs);
    a_drift=FD*(2*dust[i].sigma/M_PI/rho_peb)*\
            (v_K(r)*v_K(r)/cs/cs)/fabs(k_P_func(r));
    St_df=v_f*v_K(r)/k_P_func(r)/cs/cs/(1.-ratio_st);
    a_df=FF*St_df*2.*disk[i].sigma/M_PI/rho_peb;

    if(a_drift < a_frag && a_drift < a_df){
      f_m= FMD;
      a1 = a_drift;
    }
    else if (a_frag < a_drift && a_frag < a_df){
      f_m = FMF;
      a1  = a_frag;
    }
    else{
      f_m = FMF;
      a1  = a_df;
    }
    a_p=dust[i].a_p;
    d_v=v_pp(r,a_p,0);
    rho_dust=dust[0].sigma/dust[0].h/sqrt(2.*M_PI);
    srcterm=srcfac*2.*sqrt(M_PI)*a_p*a_p*d_v*dust[i].sigma/dust[i].h;
    dust[i].m_peb+=srcterm*dt;
    a_p=cbrt(3*dust[i].m_peb/4/M_PI/rho_peb);

    //a_p=a0*exp(tot_t/tau_grow);
    if (a_p>a1) a_p=a1;
    dust[i].a_p=a_p;
    dust[i].f_m=f_m;
    dust[i].a_frag=a_frag;
    dust[i].a_drift=a_drift;
    dust[i].a_df=a_df;
    St=stokes(dust[i].r,a_p);
    dust[i].St=St;
    dust[i].h=disk[i].h/sqrt(1+St*(1+2*St)/alpha/(1+St));
    dust[i].vr=v_r(dust[i].r,dust[i].St);
    
  }

}

#include "global_var.h"
#include "global_ex.h"
#include "ex_func.h"
#include "math.h"
#include <stdio.h>

void opa_init(){// initialize opacity profile
  int i,nr=30;
  double rad[nr],opac1d[nr];
  for(i=nr-1;i>=0;i--){
  if (i>nr-2)  opa=0.1;
  else opa=opac1d[i+1];
//  rad[i]=(r_min*0.9)*exp(i*1.0/nr*log(R_OUT*1.5/(r_min*0.9)));
  rad[i]=rmin+(rmax-rmin)*i*1.0/nr;
  opac1d[i]=opa_iter(rad[i],opa);
  printf("OPA_SAMPLE=%g\t%g\n",rad[i],opac1d[i]);
  opa_line.x[i]=rad[i];
  opa_line.y[i]=opac1d[i];
  }
  opa_line.point_num=nr;
  opa_line.begin_k2=0.0;
  opa_line.end_k2=0.0;

  line1(p_opa_line);
}

void init(double tot_t){
  tot_t*=TUNIT;
  int i;
  for(i=0;i<ring_num+1;i++){
    //dust[i].r=rmin+(rmax-rmin)*i*1.0/ring_num;
    dust[i].rf=rmin*exp(i*1.0/ring_num*log(rmax/rmin));
  }
  if (!(mdot<2e-10 && alpha>8e-4)){
    ITER=1;
    printf("TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT");
    opa_init();
  }
  ITER=0;
  for(i=0;i<ring_num;i++){
    dust[i].r=(dust[i].rf+dust[i+1].rf)/2.;
    disk[i].sigma=Sigma_gas(dust[i].r);
    disk[i].h=height(dust[i].r);
    disk[i].temp=temperature(dust[i].r);
    disk[i].rho=disk[i].sigma/disk[i].h/sqrt(2*M_PI);//density(dust[i].r);
    disk[i].visc_mol=viscosity(dust[i].r);
    disk[i].cs=sound_sp(dust[i].r);
    disk[i].yeta=yeta(dust[i].r);
    disk[i].yetavk=yeta(dust[i].r)*v_K(dust[i].r);

    printf("r=%e\t sigma_gas=%e\n",dust[i].r,\
        Sigma_gas(dust[i].r));
    if(1 || ( dust[i].r>0. && dust[i].r<199.)) \
      dust[i].sigma=Sigma_gas(dust[i].r)*dust_gas;
    dust[i].a_p=a_min;//*pow(dust[i].r,-1.57);
    double tau_grow=1./(dust[i].sigma/disk[i].sigma*w_K(dust[i].r));
    dust[i].a_gr=a_min*exp(tot_t/tau_grow);
    dust[i].m_peb=4.*M_PI*dust[i].a_p*dust[i].a_p\
                  *dust[i].a_p*rho_peb/3.;
    dust[i].Nd=dust[i].sigma/dust[i].m_peb;
    double St,alpha;
    St=stokes(dust[i].r,dust[i].a_p);
    alpha=alpha_func(dust[i].r);
    dust[i].h=disk[i].h/sqrt(1+St*(1+2*St)/alpha/(1+St));
    dust[i].St=St;
    dust[i].vr=v_r(dust[i].r,St);
    dust[i].St0=stokes(dust[i].r,a_min);
    dust[i].vr=v_r(dust[i].r,dust[i].St0);
    dust[i].f_m=0.0;
    
  }
  
}

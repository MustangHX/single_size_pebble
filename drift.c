#include "global_var.h"
#include "global_ex.h"
#include "ex_func.h"
#include "math.h"
#include <stdio.h>
double stokes(double r,double a_p){
 double t_stop;
 double rho_gas=density(r);
 double lambda=mean_path(r);
 double v_thermal=sqrt(8*k_B*temperature(r)/M_PI/mu/m_p);
 //printf("rhogas=%e\t meanpath=%e\t vthermal=%e\n",rho_gas,lambda,v_thermal);
 if ( a_p<2.25*lambda){
   t_stop=rho_peb*a_p/rho_gas/v_thermal;
 }
 else{
   t_stop = 4.0*rho_peb*a_p*a_p/9.0/rho_gas/v_thermal/lambda;
 }

 return w_K(r)*t_stop;
}

double v_r(double r, double a_p){
  double St,vk;
  St=stokes(r,a_p);
  vk=v_K(r);
 // printf("St=%e\t yeta=%e\t vk=%e\n",St,yeta(r),vk);
  return 2.*St*yeta(r)*vk/(1+St*St);
}

double v_t(double r, double a_p){//peb azimutal v respect to gas
  double St,vr;
  St=stokes(r,a_p);
  vr=v_r(r,a_p);
  return 0.5*St*vr;
}



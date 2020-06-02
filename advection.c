#include "global_var.h"
#include "global_ex.h"
#include "ex_func.h"
#include "math.h"
#include <stdio.h>
double upwind(double cfl){
  int i;
  double dt,umax=0.,u,u1,u2,qx,r,r1,r2,a_p,a_p1,a_p2; //r1 is i-1,r2 is i+1
  double qtemp[ring_num]={0.};
  for(i=0;i<ring_num-1;i++){
    //printf("vr=%e\t at r=%e\n",v_r(dust[i].r,dust[i].a_p),dust[i].r);
    r=dust[i].r*LUNIT;
    r2=dust[i+1].r*LUNIT;
    if (fabs(v_r(dust[i].r,dust[i].a_p)/(r2-r))>umax){
      //printf("vr=%e\t at r=%e\n",v_r(dust[i].r,dust[i].a_p),dust[i].r);
      umax=fabs(v_r(dust[i].r,dust[i].a_p))/(r2-r);
    }
  }
  dt=cfl/umax;
  for(i=0;i<ring_num;i++){
    qtemp[i]=dust[i].sigma;
  }
  for(i=1;i<ring_num-1;i++){
    r=dust[i].r*LUNIT;
    r1=dust[i-1].r*LUNIT;
    r2=dust[i+1].r*LUNIT;
    a_p=dust[i].a_p;
    a_p1=dust[i-1].a_p;
    a_p2=dust[i+1].a_p;
    u=v_r(r/LUNIT,a_p);
    if (u>0.){
      qx=(r*u*dust[i].sigma-r1*v_r(r1/LUNIT,a_p1)\
         *dust[i-1].sigma)/(r-r1);
    }
    else{
      qx=(r2*v_r(r2/LUNIT,a_p2)*dust[i+1].sigma\
         -r*u*dust[i].sigma)/(r2-r);
    }
    dust[i].sigma=qtemp[i]-dt*qx/r;
  }
  return dt/TUNIT;
}

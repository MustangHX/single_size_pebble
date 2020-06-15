#include "global_var.h"
#include "global_ex.h"
#include "ex_func.h"
#include "math.h"
#include <stdio.h>
double simple_advection(double cfl){
  int i;
  double dt,umax=0.,u,u1,u2,qx,r,r1,r2,St,St1,St2,a_p,a_p1,a_p2; //r1 is i-1,r2 is i+1
  double qtemp[ring_num]={0.};
  for(i=0;i<ring_num-1;i++){
    r=dust[i].r*LUNIT;
    r2=dust[i+1].r*LUNIT;
    if (fabs(v_r(dust[i].r,dust[i].St)/(r2-r))>umax){
      umax=fabs(v_r(dust[i].r,dust[i].St))/(r2-r);
    }
  }
  dt=cfl/umax;
  //dt=1e1*TUNIT;
  for(i=0;i<ring_num;i++){
    qtemp[i]=dust[i].r*LUNIT*dust[i].sigma;
  }
  for(i=1;i<ring_num-1;i++){
    r=dust[i].r*LUNIT;
    r1=dust[i-1].r*LUNIT;
    r2=dust[i+1].r*LUNIT;
    a_p=dust[i].a_p;
    St=dust[i].St;
    St1=dust[i-1].St;
    St2=dust[i+1].St;
    a_p1=dust[i-1].a_p;
    a_p2=dust[i+1].a_p;
    u=v_r(r/LUNIT,St);
    u1=v_r(r1/LUNIT,St1);
    u2=v_r(r2/LUNIT,St2);

    if (u>0.){
      
    }
  }

  return dt/TUNIT;
}



double upwind(double cfl){
  int i;
  double dt,umax=0.,u,u1,u2,qx,r,r1,r2,St,St1,St2,a_p,a_p1,a_p2; //r1 is i-1,r2 is i+1
  double qtemp[ring_num]={0.};
  double vr_fac=0.2;
  for(i=0;i<ring_num-1;i++){
    //printf("vr=%e\t at r=%e\n",v_r(dust[i].r,dust[i].a_p),dust[i].r);
    r=dust[i].r*LUNIT;
    r2=dust[i+1].r*LUNIT;
    if (fabs(v_r(dust[i].r,dust[i].St)/(r2-r))>umax){
      //printf("vr=%e\t at r=%e\n",v_r(dust[i].r,dust[i].a_p),dust[i].r);
      umax=fabs(v_r(dust[i].r,dust[i].St))/(r2-r);
    }
  }
  //dt=cfl/umax;
  dt=dt_fix*TUNIT*vr_fac;
  for(i=0;i<ring_num;i++){
    qtemp[i]=dust[i].r*LUNIT*dust[i].sigma;
  }
  for(i=1;i<ring_num-1;i++){
    r=dust[i].r*LUNIT;
    r1=dust[i-1].r*LUNIT;
    r2=dust[i+1].r*LUNIT;
    a_p=dust[i].a_p;
    St=dust[i].St;
    St1=dust[i-1].St;
    St2=dust[i+1].St;
    a_p1=dust[i-1].a_p;
    a_p2=dust[i+1].a_p;
    u=v_r(r/LUNIT,St);
    u1=v_r(r1/LUNIT,St1);
    u2=v_r(r2/LUNIT,St2);

    if (u1>0. && u>=0. && u2>=0.){
  //    qx=(r*u*dust[i].sigma-r1*v_r(r1/LUNIT,a_p1)\
         *dust[i-1].sigma)/(r-r1);
      qx=(u*qtemp[i]-u1*qtemp[i-1])/(r-r1);
      //printf("outward\n");
    }
    else if (u1<0. && u<0. && u2<0.){
  //    qx=(r2*v_r(r2/LUNIT,a_p2)*dust[i+1].sigma\
         -r*u*dust[i].sigma)/(r2-r);
      qx=(u2*qtemp[i+1]-u*qtemp[i])/(r2-r);
    }
    else if (u1>0. && u>0. && u2<0.){
      qx=(u*qtemp[i]-u1*qtemp[i-1])/(r-r1)+(u2*qtemp[i+1])/(r2-r);
    }
    else if (u1>0. && u<0. && u2<0.){
      qx=(-1.*u*qtemp[i]-u1*qtemp[i-1])/(r-r1)+(u2*qtemp[i+1])/(r2-r);
    }
    else if (u1<0. && u<0. && u2>0.){
      qx=(-1.*u*qtemp[i])/(r-r1);
    }
    else if (u1<0. && u>0. && u2>0.){
      qx=(u*qtemp[i])/(r-r1);
    }

    if((qtemp[i]-dt*qx)/r>sigdust_floor){
      dust[i].sigma=(qtemp[i]-dt*qx)/r;
    }
    else dust[i].sigma=sigdust_floor;
  }
  //boundary 
  r=dust[0].r*LUNIT;
  r2=dust[1].r*LUNIT;
  St=dust[0].St;
  St2=dust[1].St;
  u=v_r(r/LUNIT,St);
  u2=v_r(r2/LUNIT,St2);
  qx=(u2*qtemp[1]-0.*u*qtemp[0])/(r2-r);
  if((qtemp[0]-dt*qx)/r>sigdust_floor){
    dust[0].sigma=(qtemp[0]-dt*qx)/r;
  }
  else dust[0].sigma=sigdust_floor;

  r=dust[ring_num-1].r*LUNIT;
  r1=dust[ring_num-2].r*LUNIT;
  St=dust[ring_num-1].St;
  u=v_r(r/LUNIT,St);
  qx=u*(0.-qtemp[ring_num-1])/(r-r1);
  if((qtemp[ring_num-1]-dt*qx)/r>sigdust_floor){
    dust[ring_num-1].sigma=(qtemp[ring_num-1]-dt*qx)/r;
  }
  else dust[ring_num-1].sigma=sigdust_floor;


  return dt/TUNIT/vr_fac;
}


double van_Leer(double cfl){

  return dt/TUNIT/vr_fac;
}

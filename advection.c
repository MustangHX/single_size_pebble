#include "global_var.h"
#include "global_ex.h"
#include "ex_func.h"
#include "math.h"
#include <stdio.h>
double simple_advection(double cfl){
  int i;
  double dt,umax=0.,u,u1,u2,qx,r,r1,r2,St,St1,St2,Stf,Stf1,a_p,a_p01,a_p1; //r01 is i-1,r1 is i+1
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
    a_p01=dust[i-1].a_p;
    a_p1=dust[i+1].a_p;
    u=v_r(r/LUNIT,St);
    u1=v_r(r1/LUNIT,St1);
    u2=v_r(r2/LUNIT,St2);

    if (u>0.){
      
    }
  }

  return dt/TUNIT;
}



double upwind(double cfl, double dt0){
  int i;
  double dt,umax=0.,u,u01,u1,u2, uf, uf1, qx,r,r01,r1,r2,\
rf,rf1,rf01,St,St1,St2, St01,Stf,Stf1,a_p,a_p01,a_p1, a_p2; //r01 is i-1,r1 is i+1
  double qtemp[ring_num]={0.};
  double vr_fac=1.0;
  for(i=0;i<ring_num-1;i++){
    //printf("vr=%e\t at r=%e\n",v_r(dust[i].r,dust[i].a_p),dust[i].r);
    rf=dust[i].rf*LUNIT;
    rf1=dust[i+1].rf*LUNIT;
    if (fabs(v_r(dust[i].r,dust[i].St)/(rf1-rf))>umax){
      //printf("vr=%e\t at r=%e\n",v_r(dust[i].r,dust[i].a_p),dust[i].r);
      umax=fabs(v_r(dust[i].r,dust[i].St))/(rf1-rf);
    }
  }
  //dt=cfl/umax;
  dt=dt0*TUNIT*vr_fac;
  for(i=0;i<ring_num;i++){
    qtemp[i]=dust[i].r*LUNIT*dust[i].sigma;
  }
  for(i=1;i<ring_num-1;i++){
    r=dust[i].r*LUNIT;
    r01=dust[i-1].r*LUNIT;
    r1=dust[i+1].r*LUNIT;
    rf=dust[i].rf*LUNIT;
    rf1=dust[i+1].rf*LUNIT;

    a_p=dust[i].a_p;
    a_p01=dust[i-1].a_p;
    a_p1=dust[i+1].a_p;
    St=dust[i].St;
    St01=dust[i-1].St;
    St1=dust[i+1].St;
    St1=stokes(r1/LUNIT,a_p);
    //if(fabs(St1-dust[i+1].St)/St1>0.1) printf("St1_0=%e\tSt1=%e\tr1=%e\ta_p=%e\n",dust[i+1].St,St1,r1/LUNIT,a_p);
    Stf=stokes(rf/LUNIT,(a_p+a_p01)/2.);
    Stf1=stokes(rf1/LUNIT,(a_p+a_p1)/2.);
    u=v_r(r/LUNIT,St);
    uf=v_r(rf/LUNIT,Stf);
    u1=v_r(r1/LUNIT,St1);
    uf1=v_r(rf1/LUNIT,Stf1);


//    r01=dust[i-1].r*LUNIT;
//    r1=dust[i+1].r*LUNIT;
//    a_p=dust[i].a_p;
//    St=dust[i].St;
//    St01=dust[i-1].St;
//    St1=dust[i+1].St;
//    a_p01=dust[i-1].a_p;
//    a_p1=dust[i+1].a_p;
//    u=v_r(r/LUNIT,St);
//    u01=v_r(r01/LUNIT,St01);
//    u1=v_r(r1/LUNIT,St1);

    if(TWO_POP>0){
     uf=uf*dust[i].f_m+dust[i].vr0*(1.-dust[i].f_m);
     //printf("i=%d\tu=%e\f_m=t%e\n",i,u,dust[i].f_m);
     //u01=u01*dust[i-1].f_m+dust[i-1].vr0*(1.-dust[i-1].f_m);
     uf1=uf1*dust[i+1].f_m+dust[i+1].vr0*(1.-dust[i+1].f_m);
    }

    if (uf>=0. && uf1>=0.){
  //    qx=(r*u*dust[i].sigma-r1*v_r(r1/LUNIT,a_p1)\
         *dust[i-1].sigma)/(r-r1);
      qx=(uf*qtemp[i-1]-uf1*qtemp[i])/(rf1-rf);
      //printf("outward\n");
    }
    else if (uf<0. && uf1<0.){
  //    qx=(r2*v_r(r2/LUNIT,a_p2)*dust[i+1].sigma\
         -r*u*dust[i].sigma)/(r2-r);
      qx=(uf*qtemp[i]-uf1*qtemp[i+1])/(rf1-rf);
      //qx=(u*qtemp[i]-u1*qtemp[i+1])/(r1-r);
//      qx=-1.*(u1*qtemp[i+1]-u*qtemp[i])/(r1-r);

    //printf("qx=%e\t%e\t%e\t%e\t%e\n",qx,u,qtemp[i],rf1/LUNIT,rf/LUNIT);

      //if(qx<0 && r/LUNIT>2. && r/LUNIT < 3. ) printf("qx r=%e\t",r/LUNIT);
    }
    else if (uf > 0.0 && uf1 < 0.0){
      //qx=(u*qtemp[i]-u1*qtemp[i-1])/(r-r1)+(u2*qtemp[i+1])/(r2-r);
      qx=(uf*qtemp[i-1]-uf1*qtemp[i+1])/(rf1-rf);
    }
    else {//(uf<0.0 && uf1 > 0.0){
      //qx=(-1.*u*qtemp[i]-u1*qtemp[i-1])/(r-r1)+(u2*qtemp[i+1])/(r2-r);
      qx=(uf*qtemp[i]-uf1*qtemp[i])/(rf1-rf);
    }
    //else if (u1<0. && u<0. && u2>0.){
     // qx=(-1.*u*qtemp[i])/(r-r1);
   // }
    //else if (u1<0. && u>0. && u2>0.){
     // qx=(u*qtemp[i])/(r-r1);
   // }

    if((qtemp[i]+dt*qx)/r>sigdust_floor){
      dust[i].sigma=(qtemp[i]+dt*qx)/r;
    }
    else dust[i].sigma=sigdust_floor;
  }
  //boundary 
  r=dust[0].r*LUNIT;
  r1=dust[1].r*LUNIT;
  rf=dust[0].rf*LUNIT;
  rf1=dust[1].rf*LUNIT;
  St=dust[0].St;
  St01=dust[1].St;
  u=v_r(r/LUNIT,St);
  u1=v_r(r1/LUNIT,St01);
  
  a_p=dust[0].a_p;
  a_p01=dust[0].a_p;
  a_p1=dust[1].a_p;

  Stf=stokes(rf/LUNIT,(a_p+a_p01)/2.);
  Stf1=stokes(rf1/LUNIT,(a_p+a_p1)/2.);


  uf=v_r(rf/LUNIT,Stf);
  uf1=v_r(rf1/LUNIT,Stf1);
  if(TWO_POP>0){
    uf=uf*dust[0].f_m+dust[0].vr0*(1.-dust[0].f_m);
    uf1=uf1*dust[1].f_m+dust[1].vr0*(1.-dust[1].f_m);
  }

          
  qx=(0.*uf*qtemp[0]-uf1*qtemp[1])/(rf1-rf);
  if((qtemp[0]+dt*qx)/r>sigdust_floor){
    dust[0].sigma=(qtemp[0]+dt*qx)/r;
  }
  else dust[0].sigma=sigdust_floor;

  r=dust[ring_num-1].r*LUNIT;
  r01=dust[ring_num-2].r*LUNIT;
  rf=dust[ring_num-1].rf*LUNIT;
  rf1=dust[ring_num].rf*LUNIT;
  St=dust[ring_num-1].St;

  a_p=dust[ring_num-1].a_p;
  a_p01=dust[ring_num-2].a_p;
  a_p1=dust[ring_num-1].a_p;

  Stf=stokes(rf/LUNIT,(a_p+a_p01)/2.);
  Stf1=stokes(rf1/LUNIT,(a_p+a_p1)/2.);

  uf=v_r(rf/LUNIT,Stf);
  u=v_r(r/LUNIT,St);


  if(TWO_POP>0){
    uf=uf*dust[ring_num-1].f_m+dust[ring_num-1].vr0*(1.-dust[ring_num-1].f_m);
  }
  qx=uf*(qtemp[ring_num-1]-0.)/(rf1-rf);
  if((qtemp[ring_num-1]+dt*qx)/r>sigdust_floor){
    dust[ring_num-1].sigma=(qtemp[ring_num-1]+dt*qx)/r;
    //printf("qx=%e\t%e\t%e\t%e\t%e\n",qx,u,qtemp[ring_num-1],rf1/LUNIT,rf/LUNIT);
  }
  else dust[ring_num-1].sigma=sigdust_floor;


  return dt/TUNIT/vr_fac;
}


double van_Leer(double cfl){
  double dt, vr_fac;
  return dt/TUNIT/vr_fac;
}

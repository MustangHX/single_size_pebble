#include <stdio.h>
#include <math.h>
#include "global_var.h"
#include "ex_func.h"
#include "global_ex.h"

void diffusion(double dt){
  dt=dt*TUNIT;
  int i;
  double alpha,alpha01,alpha1,qx,r,r1,r2,r01,r02,a_p,a_p01,a_p1,St,St01,St1; //r01 is i-1,r1 is i+1
  double qtemp[ring_num]={0.}, q[ring_num]={0.};
  double diffco=0.01,diffco01,diffco1,Fdiff01,Fdiff1;
  for(i=0;i<ring_num;i++){
    qtemp[i]=dust[i].sigma/disk[i].sigma;
  }
  for(i=2;i<ring_num-2;i++){
    r=dust[i].r*LUNIT;
    r01=dust[i-1].r*LUNIT;
    r02=dust[i-2].r*LUNIT;
    r1=dust[i+1].r*LUNIT;
    r2=dust[i+2].r*LUNIT;

    a_p=dust[i].a_p;
    St=stokes(r/LUNIT,a_p);
    alpha=alpha_func(r/LUNIT);
    diffco=alpha*disk[i].cs*disk[i].h/(1+St*St);
    a_p01=dust[i-1].a_p;
    St01=stokes(r01/LUNIT,a_p01);
    alpha01=alpha_func(r01/LUNIT);
    diffco01=alpha*disk[i-1].cs*disk[i-1].h/(1+St01*St01);
    a_p1=dust[i+1].a_p;
    St1=stokes(r1/LUNIT,a_p1);
    alpha1=alpha_func(r1/LUNIT);
    diffco1=alpha1*disk[i+1].cs*disk[i+1].h/(1+St1*St1);



    Fdiff01=diffco01*disk[i-1].sigma*(qtemp[i-2]-qtemp[i])\
           /(r02-r);
    Fdiff1=diffco1*disk[i+1].sigma*(qtemp[i]-qtemp[i+2])\
           /(r-r2);
    //printf("diff %e\n",1./r*dt*diffco*2./(r2-r1)*((qtemp[i+1]\
                     -qtemp[i])/(r2-r)*r2*disk[i-1].sigma-\
                  (qtemp[i]-qtemp[i-1])/(r-r1)*r1*disk[i-1].sigma));
    q[i]=qtemp[i]*disk[i].sigma+1./r*dt\
         /(r01-r1)*(r01*Fdiff01-r1*Fdiff1);
    //q[i]=qtemp[i]+1./r*dt*diffco*2./(r2-r1)*((qtemp[i+1]\
         -qtemp[i])/(r2-r)*r2-\
        (qtemp[i]-qtemp[i-1])/(r-r1)*r1);
    q[i]=q[i]/disk[i].sigma;
  }
  q[0]=q[1];
  q[ring_num-1]=q[ring_num-2];
 
  for(i=0;i<ring_num;i++){
    dust[i].sigma=q[i]*disk[i].sigma;
  }
}



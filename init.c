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

void init(){
  int i;
  for(i=0;i<ring_num;i++){
    dust[i].r=rmin+(rmax-rmin)*i*1.0/ring_num;
  }
  if (!(mdot<2e-10 && alpha>8e-4)){
    ITER=1;
    printf("TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT");
    opa_init();
  }
  ITER=0;
  for(i=0;i<ring_num;i++){
    disk[i].sigma=Sigma_gas(dust[i].r);
    printf("r=%e\t sigma_gas=%e\n",dust[i].r,Sigma_gas(dust[i].r));
    dust[i].sigma=Sigma_gas(dust[i].r)*dust_gas;
    dust[i].a_p=1.0;
  }
  
}

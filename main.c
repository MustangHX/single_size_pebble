#include "global_var.h"
#include "global_ex.h"
#include "ex_func.h"
#include "math.h"
#include <stdio.h>
int main(){
  int i;
  double tsum=0.,cfl=0.3,dt=0.;
  char outdustsig[256], outdustsize[256];
  FILE *fp1, *fp2;

  mdot=mdot_init;
  init();

  while(tsum<tlim){
    mdot=mdot_init;
    dt=upwind(cfl);
    printf("time=%e\t dt=%e\n",tsum,dt);
    if( ((int)(tsum))%((int)(outp_step))==0){
      printf("out put\n");
      sprintf(outdustsig,"dust_sigma%d.txt",(int)tsum);
      sprintf(outdustsize,"dust_size%d.txt",(int)tsum);
      fp1=fopen(outdustsig,"w");
      fp2=fopen(outdustsize,"w");
      for(i=0;i<ring_num;i++){
        fprintf(fp1,"%e\t%e\n",dust[i].r,dust[i].sigma);
        fprintf(fp2,"%e\t%e\n",dust[i].r,dust[i].a_p);
      }
      fclose(fp1);
      fclose(fp2);
    }
    
    tsum+=dt;

  }
  return 0;
}

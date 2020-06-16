#include "global_var.h"
#include "global_ex.h"
#include "ex_func.h"
#include "math.h"
#include <stdio.h>
int main(){
  int i;
  double tsum=0.,cfl=0.3,dt=dt_fix,outtime=0., tot_mass=0.;
  char outdustsig[256], outdustsize[256], \
    outdustvr[256], outdustst[256], outmass[256],\
    dust_frag[256],dust_drift[256], dust_df[256];
  FILE *fp1, *fp2, *fp3, *fp4, *fp5, *fp6, *fp7, *fp8;

  mdot=mdot_init;
  init();
  check();
  fp1=fopen("rad.txt","w");
  fp2=fopen("gassig.txt","w");
  fp3=fopen("gasrho.txt","w");
  fp4=fopen("gashei.txt","w");
  fp5=fopen("gastemp.txt","w");
  fp6=fopen("yeta.txt","w");
  fp7=fopen("yetavk.txt","w");
  fp8=fopen("cs.txt","w");

  for(i=0;i<ring_num;i++){
    fprintf(fp1,"%e\n",dust[i].r);
    fprintf(fp2,"%e\n",disk[i].sigma);
    fprintf(fp3,"%e\n",disk[i].rho);
    fprintf(fp4,"%e\n",disk[i].h);
    fprintf(fp5,"%e\n",disk[i].temp);
    fprintf(fp6,"%e\n",disk[i].yeta);
    fprintf(fp7,"%e\n",disk[i].yetavk);
    fprintf(fp8,"%e\n",disk[i].cs);
  }
  fclose(fp1);
  fclose(fp2);
  fclose(fp3);
  fclose(fp4);
  fclose(fp5);
  fclose(fp6);
  fclose(fp7);
  fclose(fp8);


  fp5=fopen("mass_check.txt","w");
  for(i=0;i<ring_num;i++){
    tot_mass+=dust[i].sigma*2*M_PI*(dust[i+1].rf*dust[i+1].rf\
        -dust[i].rf*dust[i].rf)*LUNIT*LUNIT/m_earth;
  }
  fprintf(fp5,"%e\t%e\n",0.,tot_mass);

  while(tsum<tlim){
    mdot=mdot_init;
    dt=upwind(cfl);
    if (TWO_POP>0){
    grow_two_pop(dt,tsum);
    }
    else upwind_size(dt);
    //diffusion(dt);
    //if( ((int)(tsum))%((int)(outp_step))==0){
    if( tsum>=outtime){
      printf("time=%e\t dt=%e\n",tsum,dt);
      printf("out put\n");
      sprintf(outdustsig,"dust_sigma%d.txt",(int)outtime);
      sprintf(outdustsize,"dust_size%d.txt",(int)outtime);
      sprintf(outdustvr,"dust_vr%d.txt",(int)outtime);
      sprintf(outdustst,"dust_st%d.txt",(int)outtime);
      sprintf(dust_frag,"dust_frag%d.txt",(int)outtime);
      sprintf(dust_drift,"dust_drift%d.txt",(int)outtime);
      sprintf(dust_df,"dust_df%d.txt",(int)outtime);
      
      
      
      fp1=fopen(outdustsig,"w");
      fp2=fopen(outdustsize,"w");
      fp3=fopen(outdustvr,"w");
      fp4=fopen(outdustst,"w");

      fp6=fopen(dust_drift,"w");
      fp7=fopen(dust_frag,"w");
      fp8=fopen(dust_df,"w");



      for(i=0;i<ring_num;i++){
        fprintf(fp1,"%e\n",dust[i].sigma);
        fprintf(fp2,"%e\n",dust[i].a_p);
        fprintf(fp3,"%e\n",dust[i].vr);
        fprintf(fp4,"%e\n",dust[i].St);
        fprintf(fp6,"%e\n",dust[i].a_drift);
        fprintf(fp7,"%e\n",dust[i].a_frag);
        fprintf(fp8,"%e\n",dust[i].a_df);




      }
      tot_mass=0.;
      fp5=fopen("mass_check.txt","a+");
      for(i=0;i<ring_num;i++){
        tot_mass+=dust[i].sigma*2*M_PI*(dust[i+1].rf*dust[i+1].rf\
            -dust[i].rf*dust[i].rf)*LUNIT*LUNIT/m_earth;
      }
      fprintf(fp5,"%e\t%e\n",tsum,tot_mass);
      printf("%e\t%e\n",tsum,tot_mass);


      fclose(fp1);
      fclose(fp2);
      fclose(fp3);
      fclose(fp4);
      fclose(fp5);
      fclose(fp6);
      fclose(fp7);
      fclose(fp8);

      outtime+=outp_step;

    }
    
    tsum+=dt;

  }
  return 0;
}

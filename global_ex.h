#include <math.h>

#define rho_peb 3.0
//#define alpha 0.001
#define gamma0 1.4
#define m_star 1.0
//#define M_PI 3.14159265358979323846264338327950288
#define AU_km 1.49597871e8
#define outp_step 100.0  //output time step in yr
#define init_step 1.0 //inial time step
#define LUNIT 1.49597871e13
#define MSUN 1.9891e33
#define MUNIT 1.9891e33
#define mdot_init 1e-9
#define alpha_init 1e-4
#define temp0 172.043610986 //temp at 1 AU
//#define r_tran 300.365  //transition from active to passive
#define TUNIT 3.15569e7
#define k_P 2.55
#define k_B 1.38e-16
#define sig_sb 5.6704e-5
//#define opa 10.0
#define r_star 0.01395141 //T_Tauri with 3 solar radius
#define mu 2.33
#define m_p 1.660538921e-24
#define G 6.674e-8
#define peb_low_lim 1e-40
#define COAG_SW 0 // coagulation must be turned on to have fragmentation
#define FRAG_SW 0
#define DIFF_SW 1 //diffusion


#define tlim 10000. //in yr
#define peb_num 50
#define v_coag_max 500.0 //cm/s
#define v_tran_width 500.0 // cm/s
#define v_frag_min0  1500.0 //cm/s
#define v_frag_min1  2000.0 //cm/s
#define v_frag_max1  2000.0 //cm/s
#define v_frag_max0  2500.0 //cm/s
#define frag_slope -0.2
#define interaction_min_radius 0.2 // in AU, min r with peb-peb interaction
#define size_ring 0.25
#define size_min 0.01
#define size_min_inj 1e-2
#define size_min_init 1e-2
#define outp_time 1
#define NUM_LIM 100
#define peb_size_num 90
#define peb_size_lim 0.3 //in cm
#define size_step 0.05
#define size_slope 5.0
#define size_slope_pow 1.0
#define size_slope_birn 0.0
#define dust_gas 0.01
#define peb_dust 0.01
#define peb_dust_inj 0.0
#define rmax 100.0
#define rmin 10.0
#define ring_num 90
#define SINEALPHA 1
#define MDOT_INT 1 //0 for mdot=1e-8, 1 for mdot=1e-9
#define VISCOSITYRATIO 10
#define RTRAN 50.0 
#define DRTRAN 10.0
#define PEB_IMPOSE 0

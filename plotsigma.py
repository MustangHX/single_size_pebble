#!bin/python

import sys
sys.path.append("/Users/xiaohu/work/py_package")
from readcol import *
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import matplotlib.cm as cm
from matplotlib.colors import LogNorm
LUNIT=1.49597871e13
TUNIT=3.15569e7
m_earth=5.97219e27

TWO_POP=False

name=sys.argv[1]
num=sys.argv[2]

fig=plt.figure(figsize=(7,6))
gs=gridspec.GridSpec(3,2)#,width_ratios=[1,1])

r=readcol(name+"/rad.txt")

sig0=readcol(name+"/dust_sigma0.txt")
sig=readcol(name+"/dust_sigma"+num+".txt")

size0=readcol(name+"/dust_size0.txt")
size=readcol(name+"/dust_size"+num+".txt")

if TWO_POP:
  a_drift=readcol(name+"/dust_drift"+num+".txt")
  a_frag=readcol(name+"/dust_frag"+num+".txt")
  a_df=readcol(name+"/dust_df"+num+".txt")



vr0=readcol(name+"/dust_vr0.txt")
vr=readcol(name+"/dust_vr"+num+".txt")

St0=readcol(name+"/dust_st0.txt")
St=readcol(name+"/dust_st"+num+".txt")

ax1=plt.subplot(gs[0,0])
ax1.plot(r,sig0)
ax1.plot(r,sig)
ax1.set_yscale('log')
ax1.set_xscale('log')
ax1.set_ylabel('Sig_dust')
ax1.set_ylim(bottom=1e-5)
#ax1.set_xlim(4,100)
  
ax2=plt.subplot(gs[0,1])
ax2.plot(r,size0)
ax2.plot(r,size)

if TWO_POP:
 ax2.plot(r,a_drift,label="drift")
 ax2.plot(r,a_frag,label="frag")
 ax2.plot(r,a_df,label="df")
 ax2.legend()


ax2.plot(r,20*r**-1.57,ls='--',lw=1)
ax2.set_yscale('log')
ax2.set_xscale('log')
ax2.set_ylabel('a_p (cm)')
#ax2.set_xlim(4,100)
  
ax3=plt.subplot(gs[1,0])
ax3.plot(r,vr0)
ax3.plot(r,vr)
ax3.set_ylabel('vr')
ax3.set_xscale('log')
#ax3.set_xlim(4,100)
  

ax4=plt.subplot(gs[1,1])
#ax4.plot(r,-vr0*sig0*2*np.pi*r*LUNIT*TUNIT*1e6/m_earth)
#ax4.plot(r,-vr*sig*2*np.pi*r*LUNIT*TUNIT*1e6/m_earth)
ax4.plot(r,St0)
ax4.plot(r,St)
ax4.set_ylabel('St')
ax4.set_yscale('log')
ax4.set_xscale('log')
#ax4.set_xlim(4,100)
gassig=readcol(name+"/gassig.txt")
ax5=plt.subplot(gs[2,:])
ax5.plot(r,sig0/gassig)
ax5.plot(r,sig/gassig)
ax5.set_yscale('log')
ax5.set_xscale('log')
plt.show()

size_check=readcol(name+"/size_check.txt")
vpp_check=readcol(name+"/vpp_check.txt");
vBrown_check=readcol(name+"/vBrown_check.txt");
dvr_check=readcol(name+"/dvr_check.txt");
dvt_check=readcol(name+"/dvt_check.txt");
dvz_check=readcol(name+"/dvz_check.txt");
vturb_check=readcol(name+"/vturb_check.txt");


fig,ax=plt.subplots()

ax.plot(size_check,vpp_check,label='vpp',c='k')
ax.plot(size_check,vBrown_check,label='vBrown',ls='--',c='grey')
ax.plot(size_check,dvr_check,label='dvr',ls='--',c='g',alpha=0.5)
ax.plot(size_check,dvt_check,label='dvt',ls='--',c='orange',alpha=0.5)
ax.plot(size_check,dvz_check,label='dvz',c='b')
ax.plot(size_check,vturb_check,label='vturb',c='r')

ax.set_xlim(1e-5,1e3)
ax.set_ylim(1e-2,1e4)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('a_p (cm)')
ax.set_ylabel('v_pp (cm/s)')
plt.legend()
plt.show()

tgrowth100au=readcol(name+"/tgrowth100au.txt");
tgrowth1au=readcol(name+"/tgrowth_inner.txt");

fig,ax=plt.subplots()

ax.plot(size_check,tgrowth100au/TUNIT,label='100au',c='b')
ax.plot(size_check,tgrowth1au/TUNIT,label='inner',c='r')

ax.set_xlim(1e-5,1000)
#ax.set_ylim(1e0,1e6)
ax.set_xscale('log')
ax.set_yscale('log')

plt.legend()
plt.show()

"""
tgrowth1au=readcol(name+"/tgrowth1au.txt");
tgrowth2au=readcol(name+"/tgrowth2au.txt");
tgrowth3au=readcol(name+"/tgrowth3au.txt");
tgrowth4au=readcol(name+"/tgrowth4au.txt");
tgrowth5au=readcol(name+"/tgrowth5au.txt");






fig,ax=plt.subplots()
ax.plot(size_check,tgrowth1au/TUNIT,label='1au')
ax.plot(size_check,tgrowth2au/TUNIT,label='2au')
ax.plot(size_check,tgrowth3au/TUNIT,label='3au')
ax.plot(size_check,tgrowth4au/TUNIT,label='4au')
ax.plot(size_check,tgrowth5au/TUNIT,label='5au')


ax.set_xlim(1e-5,1000)
#ax.set_ylim(1e0,1e6)
ax.set_xscale('log')
ax.set_yscale('log')

plt.legend()
plt.show()

vr1au=readcol(name+"/vr1au.txt");

fig,ax=plt.subplots()
ax.plot(size_check,-vr1au,label='1au')


#ax.set_xlim(1e-5,1000)
#ax.set_ylim(1e0,1e6)
ax.set_xscale('log')
ax.set_yscale('log')

plt.legend()
plt.show()
"""
gassig=readcol(name+"/gassig.txt")
gashei=readcol(name+"/gashei.txt")
gasrho=readcol(name+"/gasrho.txt")
yeta=readcol(name+"/yeta.txt")
yetavk=readcol(name+"/yetavk.txt")
cs=readcol(name+"/cs.txt")



fig=plt.figure(figsize=(6,7))
gs=gridspec.GridSpec(6,1)

ax1=plt.subplot(gs[0,0])
ax1.plot(r,gassig)
ax1.plot(r,1700*r**-1.5,ls='--')
ax1.set_yscale('log')
ax1.set_xscale('log')
ax1.set_ylabel("gassig")


ax2=plt.subplot(gs[1,0])
ax2.plot(r,gashei/LUNIT)
ax2.plot(r,0.026*r**1.25,ls='--')
ax2.set_yscale('log')
ax2.set_xscale('log')
ax2.set_ylabel("h_g")

ax3=plt.subplot(gs[2,0])
ax3.plot(r,gasrho)
ax3.plot(r,1.7e-9*r**(-11./4.),ls='--')
ax3.set_yscale('log')
ax3.set_xscale('log')
ax3.set_ylabel("gasrho")

ax4=plt.subplot(gs[3,0])
ax4.plot(r,yeta)
ax4.plot(r,2.2e-3*r**0.5,ls='--')
ax4.set_yscale('log')
ax4.set_xscale('log')
ax4.set_ylabel("yeta")

ax5=plt.subplot(gs[4,0])
ax5.plot(r,yetavk)
ax5.plot(r,6600*r**0,ls='--')
#ax5.set_yscale('log')
ax5.set_xscale('log')
ax5.set_ylabel("yetavk")

ax6=plt.subplot(gs[5,0])
ax6.plot(r,cs)
ax6.plot(r,7.8e4*r**-0.25,ls='--')
ax6.set_yscale('log')
ax6.set_xscale('log')
ax6.set_ylabel("cs")

plt.show()

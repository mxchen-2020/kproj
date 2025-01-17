import numpy as np
import matplotlib.pyplot as plt

data=open('bs_projected.dat','r')

info=data.readline().split()
NKPOINT = int(info[0])
NSPIN   = int(info[1])
E_fermi = float(info[2])

info=data.readline().split()
KSPLIT=int(info[0])+1

KSPLIT_POS = np.empty([KSPLIT], dtype=float)
info=data.readline().split()
KSPLIT_POS[0]=0
for i in range(1,KSPLIT):
    KSPLIT_POS[i]= float(info[i-1])

#-------init----

KPOINT_DIR = np.empty([NKPOINT,3], dtype=float)
info = data.readline().split()
NBAND=int(info[1])

KPOINT = np.empty([NKPOINT], dtype=float)

KPOINT_DIR[0,0]=info[3]
KPOINT_DIR[0,1]=info[4]
KPOINT_DIR[0,2]=info[5]
KPOINT[0]=float(info[6])

EIGENV = np.empty([NSPIN,NKPOINT,NBAND], dtype=float)
WEIGHT = np.empty([NSPIN,NKPOINT,NBAND], dtype=float)

for iband in range(NBAND):
    info = data.readline().split()
    EIGENV[0,0,iband] = float(info[0])
    WEIGHT[0,0,iband] = float(info[1])

#-------star loop----

for ikpoint in range(1,NKPOINT):
    info = data.readline().split()
    KPOINT_DIR[ikpoint,0]=info[3]
    KPOINT_DIR[ikpoint,1]=info[4]
    KPOINT_DIR[ikpoint,2]=info[5]
    KPOINT[ikpoint]=float(info[6])
    for iband in range(NBAND):
        info = data.readline().split()
        EIGENV[0,ikpoint,iband] = float(info[0])
        WEIGHT[0,ikpoint,iband] = float(info[1])

if NSPIN == 2 :
    for ikpoint in range(1,NKPOINT):
        info = data.readline().split()
        for iband in range(NBAND):
           info = data.readline().split()
           EIGENV[1,ikpoint,iband] =float(info[0])
           WEIGHT[1,ikpoint,iband] =float(info[1])


EIGENV_PLOT = np.empty([NSPIN*NKPOINT*NBAND], dtype=float)
WEIGHT_PLOT = np.empty([NSPIN*NKPOINT*NBAND], dtype=float)
KPOINT_PLOT = np.empty([NSPIN*NKPOINT*NBAND], dtype=float)

n=0
ispin=0
for ispin in range(NSPIN):
    ikpoint=0
    for ikpoint in range(NKPOINT):
        iband=0
        for iband in range(NBAND):
            KPOINT_PLOT[n]=KPOINT[ikpoint]
            EIGENV_PLOT[n]=EIGENV[ispin][ikpoint][iband]
            WEIGHT_PLOT[n]=WEIGHT[ispin][ikpoint][iband]
            n=n+1
#-------------plot parameter------------------------------------------
E_fermi=-4.2304
EIGENV_PLOT=EIGENV_PLOT-E_fermi
plt.figure(figsize=(6,6))               #picture-size：x:y
plt.xlim(KPOINT[0], KPOINT[NKPOINT-1])    #x-axis(KPOINT RANGE)
plt.ylim( -20, 10)                         #y-axis(ENERGY RANGE)
plt.axhline(y=0,ls="--",c="black")   #E-Fermi
factor = 50                               #Control the size of the point           ↓
plt.scatter(KPOINT_PLOT,EIGENV_PLOT,alpha=WEIGHT_PLOT,c=WEIGHT_PLOT,cmap="Reds",s=WEIGHT_PLOT*factor,marker='o',vmin=0,vmax=1)
for i in range(2):                   #High SYMMETRY POINT
    plt.axvline(x=KSPLIT_POS[i],ls="--",c='blue',alpha=0.5)
plt.xticks([KPOINT[0],KPOINT[60],KPOINT[90]],[r'$\Gamma$',r'K1',r'M'])
plt.ylabel('Energy(ev)',size=20) 

plt.colorbar()
plt.show()

    
#已知的问题：Python版本小于3.8时,plot.scatter中的alpha参数会报错
#该参数是将投影权重表示为数据点的透明度

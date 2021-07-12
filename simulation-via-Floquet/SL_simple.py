import numpy as np
import nwlib as my
import matplotlib.pyplot as plt
import matplotlib
import Phase_Method as Pha
import Lagrange_Method as Lag
import os

pi = np.pi

P = 3e-2 # power

dt = 1e-4 # step width

Delta = 0 # ω - Ω

# Stuart-Landau
# dx/dt = x - self.a *y - (x - self.b *y) * (x**2 + y**2)
# dy/dt = self.a *x + y - (self.b *x + y) * (x**2 + y**2)

a = 11
b = 1
SL = my.SL(a, b)

omega = a - b
Omega = omega - Delta

Tsimu = 500 # simulation time
Tsimu_num = int(Tsimu / dt)

T = 2*pi/omega
Tnum = int(T / dt)

Te = 2*pi/Omega # External force term
Tenum = int(omega / Omega * Tnum) 

Devision = 40 # Number of phase measurements per cycle
PhiCount = int(Tnum/Omega*omega/Devision) # phase measurement interval

simple_SL = Lag.SL(b, omega, P, Delta)
mu, nu = simple_SL.Calc_mu_nu() # calculate lagrange multipliers

Initial = 1/4 # initial phase (2π*Initial)

x = SL.lc_theta(Initial*2*pi) 

SL_theta = np.empty((2,int(Tsimu_num/PhiCount))) # save Φ(t)
SL_X = np.empty((2,2*Tenum)) # save state X of last two cycle

PP = 0 # phase array position

for tt in range(Tsimu_num):
    # External Force Phase
    EFP = Omega * tt * dt
    if(tt%PhiCount==0):
        X_phase = SL.phase(x)
        Phi = Pha.Trans_PI_SL(X_phase - EFP) # Φ(t) [-π, π]
        SL_theta[:,PP:PP+1] = np.array([[tt*dt],[Phi]]) # record Φ(t)
        PP = PP + 1
        
    q = omega/2/nu * (-np.array([[-np.cos(EFP)+b*np.sin(EFP)],[-np.sin(EFP)-b*np.cos(EFP)]]) + mu * np.array([[-np.sin(EFP)-b*np.cos(EFP)],[np.cos(EFP)-b*np.sin(EFP)]]))
    k1 = dt*SL.dif_per(x, q)
    k2 = dt*SL.dif_per(x+0.5*k1, q)
    k3 = dt*SL.dif_per(x+0.5*k2, q)
    k4 = dt*SL.dif_per(x+k3, q)
    x += (k1+2*k2+2*k3+k4)/6
    #x += SL.dif_per(x, q) * dt

print("final phi = ", Phi)
plt.plot(SL_theta[0,:],SL_theta[1,:]) 
plt.grid()

# save state X after conversation
for tt2 in range(2*Tenum):
    tt = Tsimu_num + tt2
    EFP = tt*dt*Omega
    X_phase = SL.phase(x)
    SL_X[:,tt2:tt2+1] = x
    q = omega/2/nu * (-np.array([[-np.cos(EFP)+b*np.sin(EFP)],[-np.sin(EFP)-b*np.cos(EFP)]]) + mu * np.array([[-np.sin(EFP)-b*np.cos(EFP)],[np.cos(EFP)-b*np.sin(EFP)]]))
    x += SL.dif_per(x, q) * dt


datapath = "data/SL/simple/P{}Delta{}Initial{}/".format(P,Delta,Initial)
os.makedirs(datapath, exist_ok=True) # make folder

np.savetxt(datapath + 'mu_nu.txt', np.array([mu,nu]), delimiter = ',')
SL_theta.dump(datapath + 'theta.dat')
SL_X.dump(datapath + 'X.dat')


#%%
## plot ##########################################################

# initial parameter 

plt.rcParams['xtick.direction'] = 'in'#x軸の目盛線
plt.rcParams['ytick.direction'] = 'in'#y軸の目盛線

plt.rcParams["xtick.top"] = True            # 上部に目盛り線を描くかどうか
plt.rcParams["xtick.bottom"] = True         # 下部に目盛り線を描くかどうか
plt.rcParams["ytick.left"] = True           # 左部に目盛り線を描くかどうか
plt.rcParams["ytick.right"] = True          # 右部に目盛り線を描くかどうか

plt.rcParams["axes.linewidth"] = 2.0 
plt.rcParams["xtick.major.width"] = 2.0     # x軸主目盛り線の線幅
plt.rcParams["ytick.major.width"] = 2.0     # y軸主目盛り線の線幅

plt.rcParams["xtick.labelsize"] = 18        # 目盛りのフォントサイズ
plt.rcParams["ytick.labelsize"] = 18        # 目盛りのフォントサイズ

Delta = 0
Initial = 1/4
dt = 1e-4 # step width

figpath = "fig/SL/simple/"
os.makedirs(figpath, exist_ok=True) # make folder


############
## theta ###
############
filename = "theta.pdf"
filepath = figpath + filename

plt.figure(1)

P = 1e-3
datapath = "data/SL/simple/P{}Delta{}Initial{}/".format(P,Delta,Initial)
with open(datapath + 'theta.dat','rb') as f:
    SL0 = np.load(f, allow_pickle=True)
    
P = 3e-2
datapath = "data/SL/simple/P{}Delta{}Initial{}/".format(P,Delta,Initial)
with open(datapath + 'theta.dat','rb') as f:
    SL1 = np.load(f, allow_pickle=True)

P = 1e0
datapath = "data/SL/simple/P{}Delta{}Initial{}/".format(P,Delta,Initial)
with open(datapath + 'theta.dat','rb') as f:
    SL2 = np.load(f, allow_pickle=True)

plt.plot(SL0[0,:], SL0[1,:], label = r"$P=0.001$", color='b')
plt.plot(SL1[0,:], SL1[1,:], label = r"$P=0.03$", color='orange')
plt.plot(SL2[0,:], SL2[1,:], label = r"$P=1.0$", color = 'g')

Zeroline = np.array([[0,150],[0,0]])
plt.plot(Zeroline[0,:],Zeroline[1,:], color='k', alpha=0.8, linestyle='--')

plt.legend(fontsize = 20)

plt.xlabel('$\it{t}$', size = 20)
plt.ylabel('$\it{\phi}$', size=20)
plt.xlim(0,150)
plt.ylim(-0.35,1.65)

plt.savefig(filepath, bbox_inches="tight", pad_inches=0.0, transparent=True)

############
## X #######
############
filename = "X.pdf"
filepath = figpath + filename

plt.figure(2, figsize=(5, 5))

# when a=2 b=1, limitcycle is a cycle with a radius of 1 and the origin at the center
limitcycle = np.empty((2,100))
for tt in range(100):
    limitcycle[0,tt] = np.cos(0.1*tt)
    limitcycle[1,tt] = np.sin(0.1*tt)

P = 1e-3
datapath = "data/SL/simple/P{}Delta{}Initial{}/".format(P,Delta,Initial)
with open(datapath + 'X.dat','rb') as f:
    SL0 = np.load(f, allow_pickle=True)
    
P = 3e-2
datapath = "data/SL/simple/P{}Delta{}Initial{}/".format(P,Delta,Initial)
with open(datapath + 'X.dat','rb') as f:
    SL1 = np.load(f, allow_pickle=True)

P = 1e0
datapath = "data/SL/simple/P{}Delta{}Initial{}/".format(P,Delta,Initial)
with open(datapath + 'X.dat','rb') as f:
    SL2 = np.load(f, allow_pickle=True)

plt.xlabel(r"$x$", size = 20)
plt.ylabel(r"$y$", size = 20)

plt.plot(SL0[0,:], SL0[1,:], label = r"$P=0.001$", color='b')
plt.plot(SL1[0,:], SL1[1,:], label = r"$P=0.03$", color = 'orange')
plt.plot(SL2[0,:], SL2[1,:], label = r"$P=1.0$", color = 'g')
plt.plot(limitcycle[0,:],limitcycle[1,:], alpha = 0.6, label = "limit cycle", linestyle = '--', color ='r')

plt.yticks([-1, 0, 1])
plt.legend(fontsize=14)

plt.savefig(filepath, bbox_inches="tight", pad_inches=0.0, transparent=True)
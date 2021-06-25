#%% FitzHugh nagumo optimal inputs by Zlotnik(2013)
import numpy as np
import matplotlib.pyplot as plt
import Phase_Method as Pha
import Lagrange_Method as Lag
import os

P = 1e-4 # power

Delta = -0.005 # ω - Ω

Omega = omega - Delta
Tenum = omega / Omega * Tnum # period with external force

# calculate lagrange multipliers
print("calculate lagrange multipliers ...")
simple_2D = Lag.simple(Tnum, v0_, v0_dif, omega, P, Delta)
mu, nu = simple_2D.Calc_mu_nu()
print("nu = ", nu, " mu = ", mu)
#%%
Initial = 1/4 # initial phase (2π*Initial)
X_phase = int(Tnum*Initial)
x = np.copy(X0_[:,X_phase:X_phase+1])

Tsimu = 410 # simulation time
Tsimu_num = int(Tsimu/dt)
Devision = 40 # Number of phase measurements per cycle
PhiCount = int(Tnum/Omega*omega/Devision) # phase measurement interval

FHN_theta = np.empty((2,int(Tsimu_num/PhiCount))) # save Φ(t)
FHN_Average = np.array([[0],[X_phase/Tnum*2*pi]]) # save [Φ(t)]_t
FHN_X = np.empty((2,2*int(Tenum))) # save state X of last two cycle

FHN_input = np.empty((2,int(Tenum))) # save input

SX = Tsimu_num - 2*int(Tenum) # start of state X measurement
PhaseAverage = 0

PP = 0 # phase array position
AP = 0 # Average array position

XP = 0 # X array position

for tt in range(Tsimu_num):
    #External Force Phase
    EFP = Omega * tt / omega
    EFP = EFP - int(EFP / Tnum) * Tnum
    IDP = EFP - int(EFP) # internally dividing point
    
    if(tt%PhiCount==0):
        X_phase = Pha.Calc_phase_directory(x, FHN.dif, X0_, Tnum, dt) # measure phase [0, Tnum-1]
        Phi = Pha.Trans_PI(X_phase - EFP, Tnum) # Φ(t) [-π, π]
        FHN_theta[:,PP:PP+1] = np.array([[tt*dt],[Phi]]) # save Φ(t)
        PP = PP + 1
        PhaseAverage += Phi
        AP += 1 
    
    # Linear interpolation
    v0 = IDP * v0_[:,int(EFP):int(EFP)+1] + (1 - IDP) * np.array([v0_[:,(int(EFP)+1)%Tnum]]).T
    v0dif = IDP * v0_dif[:,int(EFP):int(EFP)+1] + (1 - IDP) * np.array([v0_dif[:,(int(EFP)+1)%Tnum]]).T     
    
    # input
    q = 1/2/nu * (-1/omega * v0dif + mu * v0)
    
    # Runge Kutta
    k1 = dt*FHN.dif_per(x, q)
    k2 = dt*FHN.dif_per(x+0.5*k1, q)
    k3 = dt*FHN.dif_per(x+0.5*k2, q)
    k4 = dt*FHN.dif_per(x+k3, q)
    x += (k1+2*k2+2*k3+k4)/6
    
    # save state X
    if tt >= SX:
        FHN_X[:,XP:XP+1] = x
        XP += 1
    
    # save phase Average [Φ(t)]_t
    if AP == Devision:
        FHN_Average = np.append(FHN_Average, np.array([[(tt-int(Tenum/2))*dt],[PhaseAverage / Devision]]), axis = 1)
        PhaseAverage = 0 # average reset
        AP = 0
    
################################################################
X_phase = Pha.Calc_phase_directory(x, FHN.dif, X0_, Tnum, dt)
print("final phi = ", Pha.Trans_PI(X_phase - EFP, Tnum))
plt.plot(FHN_theta[0,:],FHN_theta[1,:]) 
plt.grid()

# save input
for tt in range(int(Tenum)):
    #External Force Phase
    EFP = Omega * tt / omega
    EFP = EFP - int(EFP / Tnum) * Tnum
    IDP = EFP - int(EFP) # internally dividing point
   
    # Linear interpolation
    v0 = IDP * v0_[:,int(EFP):int(EFP)+1] + (1 - IDP) * np.array([v0_[:,(int(EFP)+1)%Tnum]]).T
    v0dif = IDP * v0_dif[:,int(EFP):int(EFP)+1] + (1 - IDP) * np.array([v0_dif[:,(int(EFP)+1)%Tnum]]).T     
    
    # input
    q = 1/2/nu * (-1/omega * v0dif + mu * v0)
    
    # save 
    FHN_input[:,tt:tt+1] = q
     
#%%
Gamma = simple_2D.Calc_Gamma(mu, nu)

datapath = "data/FHN/simple/P{}Delta{}Initial{}/".format(P,Delta,Initial)
os.makedirs(datapath, exist_ok=True) # make folder

np.savetxt(datapath + 'mu_nu.txt', np.array([mu,nu]), delimiter = ',')
FHN_theta.dump(datapath + 'theta.dat')
FHN_Average.dump(datapath + 'Average.dat')
FHN_X.dump(datapath + 'X.dat')
FHN_input.dump(datapath + "input.dat")
Gamma.dump(datapath + 'Gamma.dat')

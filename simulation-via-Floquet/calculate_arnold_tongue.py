#%% Van der Pol (Calculate arnold tongue) #############

# Note that the calculation time is very long.

import matplotlib.pyplot as plt
import Phase_Method as Pha
import Lagrange_Method as Lag
import numpy as np
from scipy.optimize import brenth
import os
import matplotlib
from matplotlib.colors import Normalize

def simple_Average(P, Delta):
    Omega = omega - Delta
    Tenum = omega / Omega * Tnum # period with external force

    simple_2D = Lag.simple(Tnum, v0_, v0_dif, omega, P, Delta)
    mu, nu = simple_2D.Calc_mu_nu_arnold()
    
    if mu == 0 and nu == 0:
        PhaseAverage = 1e4 # not synchronize
        return PhaseAverage
        
    X_phase = 0 # initial phase
    x = np.copy(X0_[:,X_phase:X_phase+1])
    Tsimu = int(10/np.sqrt(P)) # Adjust according to P
    Tsimu_num = int(Tsimu/dt)
    
    Division = 10
    PhiCount = int(Tenum/Division)
    SA = Tsimu_num - int(Tenum) # start of theta measurement
    ti = 0 
    PhaseAverage = 0

    
    for tt in range(Tsimu_num):
        #External Force Phase
        EFP = Omega * tt / omega
        EFP = EFP - int(EFP / Tnum) * Tnum
        IDP = EFP - int(EFP) # internally dividing point
    
        # Linear interpolation
        v0 = IDP * v0_[:,int(EFP):int(EFP)+1] + (1 - IDP) * np.array([v0_[:,(int(EFP)+1)%Tnum]]).T
        v0dif = IDP * v0_dif[:,int(EFP):int(EFP)+1] + (1 - IDP) * np.array([v0_dif[:,(int(EFP)+1)%Tnum]]).T     
    
        # input
        q = 1/2/nu * (-1/omega * v0dif + mu * v0)
    
        # Runge Kutta
        k1 = dt*VAN.dif_per(x, q)
        k2 = dt*VAN.dif_per(x+0.5*k1, q)
        k3 = dt*VAN.dif_per(x+0.5*k2, q)
        k4 = dt*VAN.dif_per(x+k3, q)
        x += (k1+2*k2+2*k3+k4)/6
    
        if tt >= SA:
            ti += 1
            if(ti%PhiCount==0):
                X_phase = Pha.Calc_phase_directory(x, VAN.dif, X0_, Tnum, dt, rotations = 5)
                Phi = Pha.Trans_PI(X_phase - EFP, Tnum) # Φ(t) [-π, π]
                PhaseAverage += Phi
    PhaseAverage = PhaseAverage/Division
    
    return PhaseAverage

def feedback_Average(P, Delta, alpha):
    Omega = omega - Delta
    Tenum = omega / Omega * Tnum # period with external force

    simple_2D = Lag.simple(Tnum, v0_, v0_dif, omega, P, Delta)
    mu, nu = simple_2D.Calc_mu_nu_arnold()
    
    if mu == 0 and nu == 0:
        PhaseAverage = 1e4 # not synchronize
        return PhaseAverage
        
    X_phase = 0 # initial phase
    x = np.copy(X0_[:,X_phase:X_phase+1])
    Tsimu = int(10/np.sqrt(P)) # Adjust according to P
    Tsimu_num = int(Tsimu/dt)
    
    Division = 10
    PhiCount = int(Tenum/Division)
    SA = Tsimu_num - int(Tenum) # start of theta measurement
    ti = 0 
    PhaseAverage = 0

    for tt in range(Tsimu_num):
        #External Force Phase
        EFP = Omega * tt / omega
        EFP = EFP - int(EFP / Tnum) * Tnum
        IDP = EFP - int(EFP) # internally dividing point
        
        y, X_phase = Pha.Calc_phase_via_Floquet(x, X0_, v0_, X_phase, Tnum)
        
        # Linear interpolation
        v0 = IDP * v0_[:,int(EFP):int(EFP)+1] + (1 - IDP) * np.array([v0_[:,(int(EFP)+1)%Tnum]]).T
        v0dif = IDP * v0_dif[:,int(EFP):int(EFP)+1] + (1 - IDP) * np.array([v0_dif[:,(int(EFP)+1)%Tnum]]).T     
    
        # input
        q = 1/2/nu * (-1/omega * v0dif + mu * v0) - alpha * y
    
        # Runge Kutta
        k1 = dt*VAN.dif_per(x, q)
        k2 = dt*VAN.dif_per(x+0.5*k1, q)
        k3 = dt*VAN.dif_per(x+0.5*k2, q)
        k4 = dt*VAN.dif_per(x+k3, q)
        x += (k1+2*k2+2*k3+k4)/6
    
        if tt >= SA:
            ti += 1
            if(ti%PhiCount==0):
                X_phase = Pha.Calc_phase_directory(x, VAN.dif, X0_, Tnum, dt, rotations = 5) # measure phase [0, Tnum-1]
                Phi = Pha.Trans_PI(X_phase - EFP, Tnum) # Φ(t) [-π, π]
                PhaseAverage += Phi
            
    PhaseAverage = PhaseAverage/Division
    return PhaseAverage

def penalty_Average(P, Delta, k):
    Omega = omega - Delta
    Tenum = omega / Omega * Tnum # period with external force

    #  Verification for existence of mu and nu
    simple_2D = Lag.simple(Tnum, v0_, v0_dif, omega, P, Delta)
    mu, nu = simple_2D.Calc_mu_nu_arnold()
    
    if mu == 0 and nu == 0:
        PhaseAverage = 1e4 # not synchronize
        return PhaseAverage
    
    penalty_2D = Lag.penalty_2D(k, Tnum, v0_, v0_dif, v1_, omega, P, Delta)

    nu = brenth(penalty_2D.Calc_nu, 1e-5, 1e3)
    mu = penalty_2D.Calc_mu(nu)
        
    X_phase = 0 # initial phase
    x = np.copy(X0_[:,X_phase:X_phase+1])
    Tsimu = int(10/np.sqrt(P)) # Adjust according to P
    Tsimu_num = int(Tsimu/dt)
    
    Division = 10
    PhiCount = int(Tenum/Division)
    SA = Tsimu_num - int(Tenum) # start of theta measurement
    ti = 0 
    PhaseAverage = 0

    for tt in range(Tsimu_num):
        #External Force Phase
        EFP = Omega * tt / omega
        EFP = EFP - int(EFP / Tnum) * Tnum
        IDP = EFP - int(EFP) # internally dividing point
        
        # Linear interpolation
        v0 = IDP * v0_[:,int(EFP):int(EFP)+1] + (1 - IDP) * np.array([v0_[:,(int(EFP)+1)%Tnum]]).T
        v0dif = IDP * v0_dif[:,int(EFP):int(EFP)+1] + (1 - IDP) * np.array([v0_dif[:,(int(EFP)+1)%Tnum]]).T     
        v1x = IDP * v1_[0,int(EFP)] + (1 - IDP) * v1_[0,(int(EFP)+1)%Tnum]
        v1y = IDP * v1_[1,int(EFP)] + (1 - IDP) * v1_[1,(int(EFP)+1)%Tnum]

        InverseMatrix = 1/(nu*nu + k *nu * (v1x*v1x + v1y*v1y))*np.array([[nu+k*v1y*v1y, -k*v1x*v1y],[-k*v1x*v1y, nu+k*v1x*v1x]])
    
        #input
        q = 1/2*np.dot(InverseMatrix, -v0dif/omega + mu * v0)
    
        # Runge Kutta
        k1 = dt*VAN.dif_per(x, q)
        k2 = dt*VAN.dif_per(x+0.5*k1, q)
        k3 = dt*VAN.dif_per(x+0.5*k2, q)
        k4 = dt*VAN.dif_per(x+k3, q)
        x += (k1+2*k2+2*k3+k4)/6
    
        if tt >= SA:
            ti += 1
            if(ti%PhiCount==0):
                X_phase = Pha.Calc_phase_directory(x, VAN.dif, X0_, Tnum, dt, rotations = 5) # measure phase [0, Tnum-1]
                Phi = Pha.Trans_PI(X_phase - EFP, Tnum) # Φ(t) [-π, π]
                PhaseAverage += Phi
            
    PhaseAverage = PhaseAverage/Division
    return PhaseAverage

#%% simple
P_number = 11
D_number = 21

arnold_x = np.empty((D_number, P_number))
arnold_y = np.empty((D_number, P_number))
arnold_Average = np.empty((D_number, P_number))
for n in range(P_number):
    P = 1e-4*(100**(1/(P_number-1)))**n*(timescale**2)
    for m in range(D_number):
        Delta = 1e-3 * (m - D_number//2) * (10**(1/(P_number-1)))**n*timescale
        Average = simple_Average(P, Delta)
        print("P = ", P, " Delta = ", Delta, " Average = ", Average)
        arnold_x[m, n] = Delta
        arnold_y[m, n] = P
        arnold_Average[m, n] = Average

maximum = 2e-1
minimum = -2e-1
fig, ax = plt.subplots(ncols=1, figsize=(6,4))
v = np.linspace(minimum, maximum, 200, endpoint=True)
ax.contourf(arnold_x, arnold_y, arnold_Average, v, norm=Normalize(minimum,maximum))
ax.set_xlim(-1, 1)
ax.set_xlabel(r"$\Delta$", fontsize=16)

# save data
datapath = "data/VAN/simple/arnold/"
os.makedirs(datapath, exist_ok=True) # make folder
arnold_x.dump(datapath + 'X.dat')
arnold_y.dump(datapath + 'Y.dat')
arnold_Average.dump(datapath + 'Average.dat')

# save fig
figpath = "fig/VAN/arnold/"
filename = "simple.pdf"
filepath = figpath + filename

os.makedirs(figpath, exist_ok=True) # make folder
plt.savefig(filepath, bbox_inches="tight", pad_inches=0.0, transparent=True)

#%% feedback

P_number = 11
D_number = 21

arnold_x = np.empty((D_number, P_number))
arnold_y = np.empty((D_number, P_number))
arnold_Average = np.empty((D_number, P_number))

alpha = 50 # feedback gain

for n in range(P_number):
    P = 1e-4*(100**(1/(P_number-1)))**n*(timescale**2)
    for m in range(D_number):
        Delta = 1e-3 * (m - D_number//2) * (10**(1/(P_number-1)))**n*timescale
        Average = feedback_Average(P, Delta, alpha)
        print("P = ", P, " Delta = ", Delta, " Average = ", Average)
        arnold_x[m, n] = Delta
        arnold_y[m, n] = P
        arnold_Average[m, n] = Average

maximum = 2e-1
minimum = -2e-1
fig, ax = plt.subplots(ncols=1, figsize=(6,4))
v = np.linspace(minimum, maximum, 200, endpoint=True)
ax.contourf(arnold_x, arnold_y, arnold_Average, v, norm=Normalize(minimum,maximum))
ax.set_xlim(-1, 1)
ax.set_xlabel(r"$\Delta$", fontsize=16)

# save data
datapath = "data/VAN/feedback/arnold/"
os.makedirs(datapath, exist_ok=True) # make folder
arnold_x.dump(datapath + 'X.dat')
arnold_y.dump(datapath + 'Y.dat')
arnold_Average.dump(datapath + 'Average.dat')

# save fig
figpath = "fig/VAN/arnold/"
filename = "feedback.pdf"
filepath = figpath + filename

os.makedirs(figpath, exist_ok=True) # make folder
plt.savefig(filepath, bbox_inches="tight", pad_inches=0.0, transparent=True)


#%% penalty

P_number = 11
D_number = 21
k = 2e1 #penalty weight

arnold_x = np.empty((D_number, P_number))
arnold_y = np.empty((D_number, P_number))
arnold_Average = np.empty((D_number, P_number))

for n in range(P_number):
    P = 1e-4*(100**(1/(P_number-1)))**n*(timescale**2)
    for m in range(D_number):
        Delta = 1e-3 * (m - D_number//2) * (10**(1/(P_number-1)))**n*timescale
        Average = penalty_Average(P, Delta, k)
        print("P = ", P, " Delta = ", Delta, " Average = ", Average)
        arnold_x[m, n] = Delta
        arnold_y[m, n] = P
        arnold_Average[m, n] = Average

maximum = 2e-1
minimum = -2e-1
fig, ax = plt.subplots(ncols=1, figsize=(6,4))
v = np.linspace(minimum, maximum, 200, endpoint=True)
ax.contourf(arnold_x, arnold_y, arnold_Average, v, norm=Normalize(minimum,maximum))
ax.set_xlim(-1, 1)
ax.set_xlabel(r"$\Delta$", fontsize=16)

datapath = "data/VAN/penalty/arnold/"
os.makedirs(datapath, exist_ok=True) # make folder
arnold_x.dump(datapath + 'X.dat')
arnold_y.dump(datapath + 'Y.dat')
arnold_Average.dump(datapath + 'Average.dat')

# save fig
figpath = "fig/VAN/arnold/"
filename = "penalty.pdf"
filepath = figpath + filename

os.makedirs(figpath, exist_ok=True) # make folder
plt.savefig(filepath, bbox_inches="tight", pad_inches=0.0, transparent=True)


#%% save colorbar

maximum = 2e-1
minimum = -2e-1

norm = matplotlib.colors.Normalize(vmin=minimum, vmax=maximum)     

fig, ax = plt.subplots(figsize=(0.25,5))
cmap = plt.get_cmap("Wistia")
cbar = matplotlib.colorbar.ColorbarBase(
    ax=ax,
    norm=norm,
    orientation="vertical",
)
cbar.set_ticks([-0.2,-0.1,0,0.1,0.2])
cbar.set_ticklabels(["-0.2","-0.1","0","0.1","0.2"])
font_size = 25 # Adjust as appropriate.
cbar.ax.tick_params(labelsize=font_size)
cbar.set_label(r"$\it{[\Phi(t)]_t}$", size=30)

# save fig
figpath = "fig/VAN/arnold/"
filename = "colorbar.pdf"
filepath = figpath + filename

os.makedirs(figpath, exist_ok=True) # make folder
plt.savefig(filepath, bbox_inches="tight", pad_inches=0.0, transparent=True)
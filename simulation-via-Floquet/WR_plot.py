#%% draw and save simulating data for Willamowski-Rossker

# initial parameter 

import numpy as np
import matplotlib.pyplot as plt 
from mpl_toolkits.mplot3d import Axes3D
import matplotlib
import os

plt.rcParams['xtick.direction'] = 'in' # draw a tick on x axis
plt.rcParams['ytick.direction'] = 'in' # draw a tick on y axis

plt.rcParams["xtick.top"] = True            # draw a tick line at the top
plt.rcParams["xtick.bottom"] = True         # draw a tick line at the bottom
plt.rcParams["ytick.left"] = True           # draw tick marks on the left
plt.rcParams["ytick.right"] = True          # draw tick marks on the right

plt.rcParams["axes.linewidth"] = 2.0 
plt.rcParams["xtick.major.width"] = 2.0     # Line width of x-axis major tick line
plt.rcParams["ytick.major.width"] = 2.0     # Line width of the y-axis major tick line

plt.rcParams["xtick.labelsize"] = 18        # Font size of x-ticks
plt.rcParams["ytick.labelsize"] = 18        # Font size of y-ticks

pi = np.pi

#%%
############################################
## feedback ################################
############################################

# initial set

P = 1e1
alpha = 1e3
Delta = -0.3
Initial = 1/4
dt = 2.5e-6 # step width

############
## X #######
############

matplotlib.rcParams['text.usetex'] = True
figpath = "fig/WR/feedback/"
os.makedirs(figpath, exist_ok=True) # make folder

filename = "X.pdf"
filepath = figpath + filename

# limitcycle
datapath = "data/WR/Floquet/"
with open(datapath + 'X0_.dat','rb') as f:
    X0_ = np.load(f, allow_pickle=True)

# simple
datapath = "data/WR/simple/P{}Delta{}Initial{}/".format(P,Delta,Initial)
with open(datapath + 'X.dat','rb') as f:
    WR_simple = np.load(f, allow_pickle=True)

# feedback
datapath = "data/WR/feedback/P{}Delta{}Initial{}alpha{}/".format(P,Delta,Initial,alpha)
with open(datapath + 'X.dat','rb') as f:
    WR_feedback = np.load(f, allow_pickle=True)

fig = plt.figure(figsize=(10,8))
ax = fig.add_subplot(111, projection='3d')

ax.set_xlabel(r"$x$", size = 30)
ax.set_ylabel(r"$y$", size = 30)
ax.set_zlabel(r"$z$", size = 30)

matplotlib.rcParams['text.usetex'] = False

plt.plot(X0_[0,:],X0_[1,:],X0_[2,:], label='limit cycle', linestyle='--', color='r')
plt.plot(WR_simple[0,:],WR_simple[1,:],WR_simple[2,:], label='without feedback', color='orange')
plt.plot(WR_feedback[0,:],WR_feedback[1,:],WR_feedback[2,:], label='with feedback', color='b', alpha=0.5)

ax.legend(fontsize=20)
ax.grid(False)
ax.set_xticks([0, 20, 40, 60])
ax.set_yticks([0, 20, 40, 60])
ax.set_zticks([0, 40, 80, 120])
ax.set_xlim(0,max(WR_simple[0,:])+1)
ax.set_ylim(0,max(WR_simple[1,:])+1)
ax.set_zlim(min(WR_simple[2,:])-1,max(WR_simple[2,:])+1)
plt.savefig(filepath, bbox_inches="tight", pad_inches=0.0, transparent=True)


#######################
## theta and Average###
#######################
filename = "theta_Average.pdf"
filepath = figpath + filename

plt.figure(2)

# simple(theta)
datapath = "data/WR/simple/P{}Delta{}Initial{}/".format(P,Delta,Initial)
with open(datapath + 'theta.dat','rb') as f:
    WR_simple = np.load(f, allow_pickle=True)

# feedback(theta)
datapath = "data/WR/feedback/P{}Delta{}Initial{}alpha{}/".format(P,Delta,Initial,alpha)
with open(datapath + 'theta.dat','rb') as f:
    WR_feedback = np.load(f, allow_pickle=True)

# simple(Average)
datapath = "data/WR/simple/P{}Delta{}Initial{}/".format(P,Delta,Initial)
with open(datapath + 'Average.dat','rb') as f:
    WR_simple_a = np.load(f, allow_pickle=True)

# feedback(Average)
datapath = "data/WR/feedback/P{}Delta{}Initial{}alpha{}/".format(P,Delta,Initial,alpha)
with open(datapath + 'Average.dat','rb') as f:
    WR_feedback_a = np.load(f, allow_pickle=True)

Zeroline = np.array([[0,50],[0,0]])

plt.plot(WR_simple[0,:],WR_simple[1,:], linestyle = '--', color='orange', alpha=0.5)
plt.plot(WR_feedback[0,:],WR_feedback[1,:], linestyle = '--', color='b', alpha=0.4)
plt.plot(WR_simple_a[0,:],WR_simple_a[1,:], linestyle = '-', label='without feedback', linewidth=2.0, color='orange')
plt.plot(WR_feedback_a[0,:],WR_feedback_a[1,:], linestyle = '-', label='with feedback', linewidth=2.0, color='b', alpha=1)
plt.plot(Zeroline[0,:],Zeroline[1,:], color='k', alpha=0.5, linestyle='--')

plt.legend(fontsize = 20)
plt.xlim(0,20)
plt.ylim(-0.2,1.7)
plt.xticks([0, 5, 10, 15, 20])
plt.xlabel('$\it{t}$', size = 20)
plt.ylabel('$\it{\phi}$', size = 20)
plt.savefig(filepath, bbox_inches="tight", pad_inches=0.0, transparent=True)


#########################
## theta and Average   ##
## Average for 2 terms ##
#########################
filename = "theta_Average_2term.pdf"
filepath = figpath + filename

plt.figure(3)

# simple(theta)
datapath = "data/WR/simple/P{}Delta{}Initial{}/".format(P,Delta,Initial)
with open(datapath + 'theta.dat','rb') as f:
    WR_simple = np.load(f, allow_pickle=True)

# feedback(theta)
datapath = "data/WR/feedback/P{}Delta{}Initial{}alpha{}/".format(P,Delta,Initial,alpha)
with open(datapath + 'theta.dat','rb') as f:
    WR_feedback = np.load(f, allow_pickle=True)

# simple(Average)
datapath = "data/WR/simple/P{}Delta{}Initial{}/".format(P,Delta,Initial)
with open(datapath + 'Average2term.dat','rb') as f:
    WR_simple_a = np.load(f, allow_pickle=True)

# feedback(Average)
datapath = "data/WR/feedback/P{}Delta{}Initial{}alpha{}/".format(P,Delta,Initial,alpha)
with open(datapath + 'Average2term.dat','rb') as f:
    WR_feedback_a = np.load(f, allow_pickle=True)

Zeroline = np.array([[0,50],[0,0]])

plt.plot(WR_simple[0,:],WR_simple[1,:], linestyle = '--', color='orange', alpha=0.5)
plt.plot(WR_feedback[0,:],WR_feedback[1,:], linestyle = '--', color='b', alpha=0.4)
plt.plot(WR_simple_a[0,:],WR_simple_a[1,:], linestyle = '-', label='without feedback', linewidth=2.0, color='orange')
plt.plot(WR_feedback_a[0,:],WR_feedback_a[1,:], linestyle = '-', label='with feedback', linewidth=2.0, color='b', alpha=1)
plt.plot(Zeroline[0,:],Zeroline[1,:], color='k', alpha=0.5, linestyle='--')

plt.legend(fontsize = 20)
plt.xlim(0,20)
plt.ylim(-0.2,1.7)
plt.xticks([0, 5, 10, 15, 20])
plt.xlabel('$\it{t}$', size = 20)
plt.ylabel('$\it{\phi}$', size = 20)

plt.savefig(filepath, bbox_inches="tight", pad_inches=0.0, transparent=True)

#%% 
############################################
## penalty ################################
############################################

#initial set

P = 1e1
k = 1e-1
Delta = -0.3
Initial = 1/4
dt = 2.5e-6 # step width

############
## X #######
############

matplotlib.rcParams['text.usetex'] = True
figpath = "fig/WR/penalty/"
os.makedirs(figpath, exist_ok=True) # make folder

filename = "X.pdf"
filepath = figpath + filename

# limitcycle
datapath = "data/WR/Floquet/"
with open(datapath + 'X0_.dat','rb') as f:
    X0_ = np.load(f, allow_pickle=True)

# simple
datapath = "data/WR/simple/P{}Delta{}Initial{}/".format(P,Delta,Initial)
with open(datapath + 'X.dat','rb') as f:
    WR_simple = np.load(f, allow_pickle=True)

# penalty
datapath = "data/WR/penalty/P{}Delta{}Initial{}k{}/".format(P,Delta,Initial,k)
with open(datapath + 'X.dat','rb') as f:
    WR_penalty = np.load(f, allow_pickle=True)

fig = plt.figure(figsize=(10.0, 8.0))
ax = fig.add_subplot(111, projection='3d')

ax.set_xlabel(r"$x$", size = 30)
ax.set_ylabel(r"$y$", size = 30)
ax.set_zlabel(r"$z$", size = 30)

matplotlib.rcParams['text.usetex'] = False

plt.plot(X0_[0,:],X0_[1,:],X0_[2,:], label='limit cycle', linestyle='--', color='r')
plt.plot(WR_simple[0,:],WR_simple[1,:], WR_simple[2,:], label='without penalty', color='orange')
plt.plot(WR_penalty[0,:],WR_penalty[1,:], WR_penalty[2,:], label='with penalty', color='b', alpha=0.5)

ax.legend(fontsize=20)
ax.grid(False)
ax.set_xticks([0, 20, 40, 60])
ax.set_yticks([0, 20, 40, 60])
ax.set_zticks([0, 40, 80, 120])
ax.set_xlim(0,max(WR_simple[0,:])+1)
ax.set_ylim(0,max(WR_simple[1,:])+1)
ax.set_zlim(min(WR_simple[2,:])-1,max(WR_simple[2,:])+1)
plt.savefig(filepath, bbox_inches="tight", pad_inches=0.0, transparent=True)


#######################
## theta and Average ##
#######################
filename = "theta_Average.pdf"
filepath = figpath + filename

plt.figure(2)

# simple(theta)
datapath = "data/WR/simple/P{}Delta{}Initial{}/".format(P,Delta,Initial)
with open(datapath + 'theta.dat','rb') as f:
    WR_simple = np.load(f, allow_pickle=True)

# penalty(theta)
datapath = "data/WR/penalty/P{}Delta{}Initial{}k{}/".format(P,Delta,Initial,k)
with open(datapath + 'theta.dat','rb') as f:
    WR_penalty = np.load(f, allow_pickle=True)

# simple(Average)
datapath = "data/WR/simple/P{}Delta{}Initial{}/".format(P,Delta,Initial)
with open(datapath + 'Average.dat','rb') as f:
    WR_simple_a = np.load(f, allow_pickle=True)

# penalty(Average)
datapath = "data/WR/penalty/P{}Delta{}Initial{}k{}/".format(P,Delta,Initial,k)
with open(datapath + 'Average.dat','rb') as f:
    WR_penalty_a = np.load(f, allow_pickle=True)

Zeroline = np.array([[0,50],[0,0]])

plt.plot(WR_simple[0,:],WR_simple[1,:], linestyle = '--', color='orange', alpha=0.5)
plt.plot(WR_penalty[0,:],WR_penalty[1,:], linestyle = '--', color='b', alpha=0.4)
plt.plot(WR_simple_a[0,:],WR_simple_a[1,:], linestyle = '-', label='without penalty', linewidth=2.0, color='orange')
plt.plot(WR_penalty_a[0,:],WR_penalty_a[1,:], linestyle = '-', label='with penalty', linewidth=2.0, color='b', alpha=1)
plt.plot(Zeroline[0,:],Zeroline[1,:], color='k', alpha=0.5, linestyle='--')

plt.legend(fontsize=20)
plt.xlim(0,20)
plt.ylim(-0.2,1.7)
plt.xticks([0, 5, 10, 15, 20])
plt.xlabel('$\it{t}$', size = 20)
plt.ylabel('$\it{\phi}$', size=20)
plt.savefig(filepath, bbox_inches="tight", pad_inches=0.0, transparent=True)

#########################
## theta and Average   ##
## Average for 2 terms ##
#########################
filename = "theta_Average_2term.pdf"
filepath = figpath + filename

plt.figure(3)

# simple(theta)
datapath = "data/WR/simple/P{}Delta{}Initial{}/".format(P,Delta,Initial)
with open(datapath + 'theta.dat','rb') as f:
    WR_simple = np.load(f, allow_pickle=True)

# penalty(theta)
datapath = "data/WR/penalty/P{}Delta{}Initial{}k{}/".format(P,Delta,Initial,k)
with open(datapath + 'theta.dat','rb') as f:
    WR_penalty = np.load(f, allow_pickle=True)

# simple(Average)
datapath = "data/WR/simple/P{}Delta{}Initial{}/".format(P,Delta,Initial)
with open(datapath + 'Average2term.dat','rb') as f:
    WR_simple_a = np.load(f, allow_pickle=True)

# penalty(Average)
datapath = "data/WR/penalty/P{}Delta{}Initial{}k{}/".format(P,Delta,Initial,k)
with open(datapath + 'Average2term.dat','rb') as f:
    WR_penalty_a = np.load(f, allow_pickle=True)

Zeroline = np.array([[0,50],[0,0]])

plt.plot(WR_simple[0,:],WR_simple[1,:], linestyle = '--', color='orange', alpha=0.5)
plt.plot(WR_penalty[0,:],WR_penalty[1,:], linestyle = '--', color='b', alpha=0.4)
plt.plot(WR_simple_a[0,:],WR_simple_a[1,:], linestyle = '-', label='without penalty', linewidth=2.0, color='orange')
plt.plot(WR_penalty_a[0,:],WR_penalty_a[1,:], linestyle = '-', label='with penalty', linewidth=2.0, color='b', alpha=1)
plt.plot(Zeroline[0,:],Zeroline[1,:], color='k', alpha=0.5, linestyle='--')

plt.legend(fontsize = 20)
matplotlib.rcParams['text.usetex'] = True

plt.xlim(0,20)
plt.ylim(-0.2,1.7)
plt.xticks([0, 5, 10, 15, 20])
plt.xlabel(r'$\it{t}$', size = 20)
plt.ylabel(r'$\it{\phi}$', size=20)
matplotlib.rcParams['text.usetex'] = False

plt.savefig(filepath, bbox_inches="tight", pad_inches=0.0, transparent=True)

##############
###  Gamma ###
##############

figpath = "fig/WR/penalty/"
os.makedirs(figpath, exist_ok=True) # make folder

filename = "Gamma.pdf"
filepath = figpath + filename

plt.figure(4)

# simple
datapath = "data/WR/simple/P{}Delta{}Initial{}/".format(P,Delta,Initial)
with open(datapath + 'Gamma.dat','rb') as f:
    WR_simple = np.load(f, allow_pickle=True)

# penalty
datapath = "data/WR/penalty/P{}Delta{}Initial{}k{}/".format(P,Delta,Initial,k)
with open(datapath + 'Gamma.dat','rb') as f:
    WR_penalty = np.load(f, allow_pickle=True)

plt.xlabel('$\it{\phi}$', size=20)
plt.ylabel('$\it{\Delta + \Gamma(\phi)}$', size=20)

Zeroline = np.array([[-pi,pi],[0,0]])
Zeroline2 = np.array([[0,0],[-1,1]])

plt.plot(WR_simple[0,:],WR_simple[1,:], label='without penalty', color='orange')
plt.plot(WR_penalty[0,:],WR_penalty[1,:], label='with penalty', color='b', alpha=0.5)
plt.plot(Zeroline[0,:],Zeroline[1,:], color='k', alpha=1.0, linestyle='--')
plt.plot(Zeroline2[0,:],Zeroline2[1,:], color='k', alpha=1.0, linestyle='--')

plt.xlim(-pi,pi)
plt.ylim(-0.90,0.90)
plt.xticks([-pi, -pi/2, 0, pi/2, pi], ["-π", "-π/2", "0", "π/2", "π"])
plt.yticks([-0.75, -0.50, -0.25, 0, 0.25, 0.50, 0.75])

plt.legend()
plt.savefig(filepath, bbox_inches="tight", pad_inches=0.0, transparent=True)


###########
## input ##
###########
filename = "input.pdf"
filepath = figpath + filename

plt.figure(5)

# simple(theta)
datapath = "data/WR/simple/P{}Delta{}Initial{}/".format(P,Delta,Initial)
with open(datapath + 'input.dat','rb') as f:
    WR_simple = np.load(f, allow_pickle=True)

# penalty(theta)
datapath = "data/WR/penalty/P{}Delta{}Initial{}k{}/".format(P,Delta,Initial,k)
with open(datapath + 'input.dat','rb') as f:
    WR_penalty = np.load(f, allow_pickle=True)

Te = WR_simple.shape[1]*dt
Time_ax = np.linspace(0,WR_simple.shape[1]-1, WR_simple.shape[1])*dt

plt.plot(Time_ax, WR_simple[0,:], label='without penalty (qx)', linestyle = '--', color='orange', alpha=1)
plt.plot(Time_ax, WR_simple[1,:], label='without penalty (qy)', linestyle = '--', color='green', alpha=1)
plt.plot(Time_ax, WR_simple[2,:], label='without penalty (qz)', linestyle = '--', color='r', alpha=1)
plt.plot(Time_ax, WR_penalty[0,:], label='with penalty (qx)', linestyle = '-', color='b', alpha=0.7)
plt.plot(Time_ax, WR_penalty[1,:], label='with penalty (qy)', linestyle = '-', color='c', alpha=0.7)
plt.plot(Time_ax, WR_penalty[2,:], label='with penalty (qz)', linestyle = '-', color='m', alpha=0.7)


plt.legend(fontsize=8.5)
plt.xlim(0,Te)
plt.xticks([0, Te/4, Te/2, 3*Te/4, Te], [ "0", "Te/4", "Te/2", "3Te/4", "Te"])
#plt.ylim(-0.3,1.7)
plt.xlabel('$\it{t}$', size = 20)
plt.ylabel('$\it{q}$', size = 20)
plt.savefig(filepath, bbox_inches="tight", pad_inches=0.0, transparent=True)






#%% draw and save simulating data for van der Pol

# initial parameter

import numpy as np
import matplotlib.pyplot as plt 
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

P = 1e0
alpha = 50
Delta = -0.5
Initial = 1/4
dt = 2.5e-5 # step width

############
## X #######
############

figpath = "fig/VAN/feedback/"
os.makedirs(figpath, exist_ok=True) # make folder

filename = "X.pdf"
filepath = figpath + filename

plt.figure(0)

# limitcycle
datapath = "data/VAN/Floquet/"
with open(datapath + 'X0_.dat','rb') as f:
    X0_ = np.load(f, allow_pickle=True)

# simple
datapath = "data/VAN/simple/P{}Delta{}Initial{}/".format(P,Delta,Initial)
with open(datapath + 'X.dat','rb') as f:
    VAN_simple = np.load(f, allow_pickle=True)

# feedback
datapath = "data/VAN/feedback/P{}Delta{}Initial{}alpha{}/".format(P,Delta,Initial,alpha)
with open(datapath + 'X.dat','rb') as f:
    VAN_feedback = np.load(f, allow_pickle=True)

plt.xticks([-1, 0, 1])
plt.yticks([-1, 0, 1])
plt.xlabel(r"$x$", size = 20)
plt.ylabel(r"$y$", size = 20)

plt.plot(X0_[0,:],X0_[1,:], label='limit cycle', linestyle='--', color='r')
plt.plot(VAN_simple[0,:],VAN_simple[1,:], label='without feedback', color='orange')
plt.plot(VAN_feedback[0,:],VAN_feedback[1,:], label='with feedback', color='b', alpha=0.5)

plt.legend(fontsize=14)
plt.savefig(filepath, bbox_inches="tight", pad_inches=0.0, transparent=True)


#######################
## theta and Average###
#######################
filename = "theta_Average.pdf"
filepath = figpath + filename

plt.figure(1)

# simple(theta)
datapath = "data/VAN/simple/P{}Delta{}Initial{}/".format(P,Delta,Initial)
with open(datapath + 'theta.dat','rb') as f:
    VAN_simple = np.load(f, allow_pickle=True)

# feedback(theta)
datapath = "data/VAN/feedback/P{}Delta{}Initial{}alpha{}/".format(P,Delta,Initial,alpha)
with open(datapath + 'theta.dat','rb') as f:
    VAN_feedback = np.load(f, allow_pickle=True)

# simple(Average)
datapath = "data/VAN/simple/P{}Delta{}Initial{}/".format(P,Delta,Initial)
with open(datapath + 'Average.dat','rb') as f:
    VAN_simple_a = np.load(f, allow_pickle=True)

# feedback(Average)
datapath = "data/VAN/feedback/P{}Delta{}Initial{}alpha{}/".format(P,Delta,Initial,alpha)
with open(datapath + 'Average.dat','rb') as f:
    VAN_feedback_a = np.load(f, allow_pickle=True)

Zeroline = np.array([[0,100],[0,0]])

plt.plot(VAN_simple[0,:],VAN_simple[1,:], linestyle = '--', color='orange', alpha=0.5)
plt.plot(VAN_feedback[0,:],VAN_feedback[1,:], linestyle = '--', color='b', alpha=0.4)
plt.plot(VAN_simple_a[0,:],VAN_simple_a[1,:], linestyle = '-', label='without feedback', linewidth=2.0, color='orange')
plt.plot(VAN_feedback_a[0,:],VAN_feedback_a[1,:], linestyle = '-', label='with feedback', linewidth=2.0, color='b', alpha=1)
plt.plot(Zeroline[0,:],Zeroline[1,:], color='k', alpha=0.5, linestyle='--')

plt.legend(fontsize=20)
plt.xlim(0,20)
plt.ylim(-0.3,1.7)
plt.xlabel('$\it{t}$', size = 20)
plt.ylabel('$\it{\phi}$', size = 20)
plt.savefig(filepath, bbox_inches="tight", pad_inches=0.0, transparent=True)
 
#%%
############################################
## penalty ################################
############################################

# initial set

P = 1e0
k = 20
Delta = -0.5
Initial = 1/4
dt = 2.5e-5 # step width

figpath = "fig/VAN/penalty/"
os.makedirs(figpath, exist_ok=True) # make folder

############
## X #######
############

filename = "X.pdf"
filepath = figpath + filename

plt.figure(0)

# limitcycle
datapath = "data/VAN/Floquet/"
with open(datapath + 'X0_.dat','rb') as f:
    X0_ = np.load(f, allow_pickle=True)

# simple
datapath = "data/VAN/simple/P{}Delta{}Initial{}/".format(P,Delta,Initial)
with open(datapath + 'X.dat','rb') as f:
    VAN_simple = np.load(f, allow_pickle=True)

# penalty
datapath = "data/VAN/penalty/P{}Delta{}Initial{}k{}/".format(P,Delta,Initial,k)
with open(datapath + 'X.dat','rb') as f:
    VAN_penalty = np.load(f, allow_pickle=True)

plt.xticks([-1, 0, 1])
plt.yticks([-1, 0, 1])
plt.xlabel(r"$x$", size = 20)
plt.ylabel(r"$y$", size = 20)

plt.plot(X0_[0,:],X0_[1,:], label='limit cycle', linestyle='--', color='r')
plt.plot(VAN_simple[0,:],VAN_simple[1,:], label='without penalty', color='orange')
plt.plot(VAN_penalty[0,:],VAN_penalty[1,:], label='with penalty', color='b', alpha=0.5)

plt.legend(fontsize=14)
plt.savefig(filepath, bbox_inches="tight", pad_inches=0.0, transparent=True)


#######################
## theta and Average ##
#######################
filename = "theta_Average.pdf"
filepath = figpath + filename

plt.figure(1)

# simple(theta)
datapath = "data/VAN/simple/P{}Delta{}Initial{}/".format(P,Delta,Initial)
with open(datapath + 'theta.dat','rb') as f:
    VAN_simple = np.load(f, allow_pickle=True)

# penalty(theta)
datapath = "data/VAN/penalty/P{}Delta{}Initial{}k{}/".format(P,Delta,Initial,k)
with open(datapath + 'theta.dat','rb') as f:
    VAN_penalty = np.load(f, allow_pickle=True)

# simple(Average)
datapath = "data/VAN/simple/P{}Delta{}Initial{}/".format(P,Delta,Initial)
with open(datapath + 'Average.dat','rb') as f:
    VAN_simple_a = np.load(f, allow_pickle=True)

# penalty(Average)
datapath = "data/VAN/penalty/P{}Delta{}Initial{}k{}/".format(P,Delta,Initial,k)
with open(datapath + 'Average.dat','rb') as f:
    VAN_penalty_a = np.load(f, allow_pickle=True)

Zeroline = np.array([[0,100],[0,0]])

plt.plot(VAN_simple[0,:],VAN_simple[1,:], linestyle = '--', color='orange', alpha=0.5)
plt.plot(VAN_penalty[0,:],VAN_penalty[1,:], linestyle = '--', color='b', alpha=0.4)
plt.plot(VAN_simple_a[0,:],VAN_simple_a[1,:], linestyle = '-', label='without penalty', linewidth=2.0, color='orange')
plt.plot(VAN_penalty_a[0,:],VAN_penalty_a[1,:], linestyle = '-', label='with penalty', linewidth=2.0, color='b', alpha=1)
plt.plot(Zeroline[0,:],Zeroline[1,:], color='k', alpha=0.5, linestyle='--')

plt.legend(fontsize=20)
plt.xlim(0,20)
plt.ylim(-0.3,1.7)
plt.xlabel('$\it{t}$', size = 20)
plt.ylabel('$\it{\phi}$', size = 20)
plt.savefig(filepath, bbox_inches="tight", pad_inches=0.0, transparent=True)


##############
###  Gamma ###
##############


filename = "Gamma.pdf"
filepath = figpath + filename

plt.figure(2)

# simple
datapath = "data/VAN/simple/P{}Delta{}Initial{}/".format(P,Delta,Initial)
with open(datapath + 'Gamma.dat','rb') as f:
    VAN_simple = np.load(f, allow_pickle=True)

# penalty
datapath = "data/VAN/penalty/P{}Delta{}Initial{}k{}/".format(P,Delta,Initial,k)
with open(datapath + 'Gamma.dat','rb') as f:
    VAN_penalty = np.load(f, allow_pickle=True)

plt.xlabel('$\it{\phi}$', size=20)
plt.ylabel('$\it{\Delta + \Gamma(\phi)}$', size=20)

Zeroline = np.array([[-pi,pi],[0,0]])
Zeroline2 = np.array([[0,0],[-1,1]])

plt.plot(VAN_simple[0,:],VAN_simple[1,:], label='without penalty', color='orange')
plt.plot(VAN_penalty[0,:],VAN_penalty[1,:], label='with penalty', color='b', alpha=0.5)
plt.plot(Zeroline[0,:],Zeroline[1,:], color='k', alpha=1.0, linestyle='--')
plt.plot(Zeroline2[0,:],Zeroline2[1,:], color='k', alpha=1.0, linestyle='--')

plt.xlim(-pi,pi)
plt.ylim(-1.5,0.5)
plt.xticks([-pi, -pi/2, 0, pi/2, pi], ["-π", "-π/2", "0", "π/2", "π"])

plt.legend(fontsize=15)
plt.savefig(filepath, bbox_inches="tight", pad_inches=0.0, transparent=True)


###########
## input ##
###########
filename = "input.pdf"
filepath = figpath + filename

plt.figure(3)

# simple(theta)
datapath = "data/VAN/simple/P{}Delta{}Initial{}/".format(P,Delta,Initial)
with open(datapath + 'input.dat','rb') as f:
    VAN_simple = np.load(f, allow_pickle=True)

# penalty(theta)
datapath = "data/VAN/penalty/P{}Delta{}Initial{}k{}/".format(P,Delta,Initial,k)
with open(datapath + 'input.dat','rb') as f:
    VAN_penalty = np.load(f, allow_pickle=True)

Te = VAN_simple.shape[1]*dt
Time_ax = np.linspace(0,VAN_simple.shape[1]-1, VAN_simple.shape[1])*dt

plt.plot(Time_ax, VAN_simple[0,:], label='without penalty (qx)', linestyle = '--', color='orange', alpha=1)
plt.plot(Time_ax, VAN_simple[1,:], label='without penalty (qy)', linestyle = '--', color='green', alpha=1)
plt.plot(Time_ax, VAN_penalty[0,:], label='with penalty (qx)', linestyle = '-', color='b', alpha=0.7)
plt.plot(Time_ax, VAN_penalty[1,:], label='with penalty (qy)', linestyle = '-', color='c', alpha=0.7)

plt.legend()
plt.xlim(0,Te)
plt.xticks([0, Te/4, Te/2, 3*Te/4, Te], [ "0", "Te/4", "Te/2", "3Te/4", "Te"])
#plt.ylim(-0.3,1.7)
plt.xlabel('$\it{t}$', size = 20)
plt.ylabel('$\it{q}$', size = 20)
plt.savefig(filepath, bbox_inches="tight", pad_inches=0.0, transparent=True)
 

#%%
############################################
## feedback and penalty ####################
############################################

# initial set

alpha = 50
k = 10
Delta = 0
Initial = 1/4

############
## X #######
############

figpath = "fig/VAN/feedback_penalty/"
os.makedirs(figpath, exist_ok=True) # make folder

filename = "X.pdf"
filepath = figpath + filename

plt.figure(0)

# limitcycle
datapath = "data/VAN/Floquet/"
with open(datapath + 'X0_.dat','rb') as f:
    X0_ = np.load(f, allow_pickle=True)

P = 1e-2
# simple(weak)
datapath = "data/VAN/simple/P{}Delta{}Initial{}/".format(P,Delta,Initial)
with open(datapath + 'X.dat','rb') as f:
    VAN_simple = np.load(f, allow_pickle=True)

P = 1e0
# simple(strong)
datapath = "data/VAN/simple/P{}Delta{}Initial{}/".format(P,Delta,Initial)
with open(datapath + 'X.dat','rb') as f:
    VAN_simple2 = np.load(f, allow_pickle=True)
    
# feedback
datapath = "data/VAN/feedback/P{}Delta{}Initial{}alpha{}/".format(P,Delta,Initial,alpha)
with open(datapath + 'X.dat','rb') as f:
    VAN_feedback = np.load(f, allow_pickle=True)

# penalty
datapath = "data/VAN/penalty/P{}Delta{}Initial{}k{}/".format(P,Delta,Initial,k)
with open(datapath + 'X.dat','rb') as f:
    VAN_penalty = np.load(f, allow_pickle=True)

plt.xlabel(r"$x$", size = 20)
plt.ylabel(r"$y$", size = 20)


plt.plot(VAN_simple[0,:],VAN_simple[1,:], label="simple(P=0.01)", color='orange')
plt.plot(VAN_simple2[0,:],VAN_simple2[1,:], label="simple(P=1.0)", color='c')
plt.plot(VAN_feedback[0,:],VAN_feedback[1,:], label='with feedback(P=1.0)', color='darkblue')
plt.plot(VAN_penalty[0,:],VAN_penalty[1,:], color='lawngreen', alpha = 0.7, label="with penalty(P=1.0)")
plt.plot(X0_[0,:],X0_[1,:], label='limit cycle', linestyle='--', color='r')

plt.xticks([-1,0,1])
plt.legend(fontsize=12)
plt.savefig(filepath, bbox_inches="tight", pad_inches=0.0, transparent=True)


##############
##  Average ##
##############

filename = "Average.pdf"
filepath = figpath + filename

plt.figure(1)

P = 1e-2
# simple(weak)
datapath = "data/VAN/simple/P{}Delta{}Initial{}/".format(P,Delta,Initial)
with open(datapath + 'Average.dat','rb') as f:
    VAN_simple = np.load(f, allow_pickle=True)

P = 1e0
# simple(strong)
datapath = "data/VAN/simple/P{}Delta{}Initial{}/".format(P,Delta,Initial)
with open(datapath + 'Average.dat','rb') as f:
    VAN_simple2 = np.load(f, allow_pickle=True)
    
# feedback
datapath = "data/VAN/feedback/P{}Delta{}Initial{}alpha{}/".format(P,Delta,Initial,alpha)
with open(datapath + 'Average.dat','rb') as f:
    VAN_feedback = np.load(f, allow_pickle=True)

# penalty
datapath = "data/VAN/penalty/P{}Delta{}Initial{}k{}/".format(P,Delta,Initial,k)
with open(datapath + 'Average.dat','rb') as f:
    VAN_penalty = np.load(f, allow_pickle=True)

plt.xlabel(r"$x$", size = 20)
plt.ylabel(r"$y$", size = 20)

plt.plot(VAN_simple[0,:],VAN_simple[1,:], label="simple(P=0.01)", color='orange')
plt.plot(VAN_simple2[0,:],VAN_simple2[1,:], label="simple(P=1.0)", color='c')
plt.plot(VAN_feedback[0,:],VAN_feedback[1,:], label='with feedback(P=1.0)', color='darkblue')
plt.plot(VAN_penalty[0,:],VAN_penalty[1,:], color='lawngreen', alpha = 0.7, label="with penalty(P=1.0)")

Zeroline = np.array([[0,200],[0,0]])

plt.plot(Zeroline[0,:],Zeroline[1,:], color='k', alpha=0.5, linestyle='--')

plt.legend(fontsize=15)
plt.xlim(0,60)
plt.ylim(-0.3,1.7)
plt.xlabel('$\it{t}$', size = 20)
plt.ylabel('$[\it{\phi}(t)]_t$', size = 20)
plt.savefig(filepath, bbox_inches="tight", pad_inches=0.0, transparent=True)

########################
##  theta and Average ##
########################

filename = "theta_Average.pdf"
filepath = figpath + filename

plt.figure(2)

P = 1e-2
# simple(weak)
datapath = "data/VAN/simple/P{}Delta{}Initial{}/".format(P,Delta,Initial)
with open(datapath + 'theta.dat','rb') as f:
    VAN_simple = np.load(f, allow_pickle=True)
    
datapath = "data/VAN/simple/P{}Delta{}Initial{}/".format(P,Delta,Initial)
with open(datapath + 'Average.dat','rb') as f:
    VAN_simple_a = np.load(f, allow_pickle=True)

P = 1e0
# simple(strong)
datapath = "data/VAN/simple/P{}Delta{}Initial{}/".format(P,Delta,Initial)
with open(datapath + 'theta.dat','rb') as f:
    VAN_simple2 = np.load(f, allow_pickle=True)
    
datapath = "data/VAN/simple/P{}Delta{}Initial{}/".format(P,Delta,Initial)
with open(datapath + 'Average.dat','rb') as f:
    VAN_simple2_a = np.load(f, allow_pickle=True)
    
# feedback
datapath = "data/VAN/feedback/P{}Delta{}Initial{}alpha{}/".format(P,Delta,Initial,alpha)
with open(datapath + 'theta.dat','rb') as f:
    VAN_feedback = np.load(f, allow_pickle=True)
    
datapath = "data/VAN/feedback/P{}Delta{}Initial{}alpha{}/".format(P,Delta,Initial,alpha)
with open(datapath + 'Average.dat','rb') as f:
    VAN_feedback_a = np.load(f, allow_pickle=True)

# penalty
datapath = "data/VAN/penalty/P{}Delta{}Initial{}k{}/".format(P,Delta,Initial,k)
with open(datapath + 'theta.dat','rb') as f:
    VAN_penalty = np.load(f, allow_pickle=True)
    
datapath = "data/VAN/penalty/P{}Delta{}Initial{}k{}/".format(P,Delta,Initial,k)
with open(datapath + 'Average.dat','rb') as f:
    VAN_penalty_a = np.load(f, allow_pickle=True)

plt.xlabel(r"$x$", size = 20)
plt.ylabel(r"$y$", size = 20)

plt.plot(VAN_simple[0,:],VAN_simple[1,:], linestyle="--", alpha = 0.4, color='orange')
plt.plot(VAN_simple2[0,:],VAN_simple2[1,:], linestyle="--", alpha = 0.4, color='c')
plt.plot(VAN_feedback[0,:],VAN_feedback[1,:], linestyle="--", alpha = 0.4, color='darkblue')
plt.plot(VAN_penalty[0,:],VAN_penalty[1,:], color='lawngreen', linestyle="--", alpha = 0.4)
plt.plot(VAN_simple_a[0,:],VAN_simple_a[1,:], label="simple(P=0.01)", color='orange')
plt.plot(VAN_simple2_a[0,:],VAN_simple2_a[1,:], label="simple(P=1.0)", color='c')
plt.plot(VAN_feedback_a[0,:],VAN_feedback_a[1,:], label='with feedback(P=1.0)', color='darkblue')
plt.plot(VAN_penalty_a[0,:],VAN_penalty_a[1,:], color='lawngreen', label="with penalty(P=1.0)")

Zeroline = np.array([[0,200],[0,0]])

plt.plot(Zeroline[0,:],Zeroline[1,:], color='k', alpha=0.5, linestyle='--')

plt.legend(fontsize=15)
plt.xlim(0,60)
plt.ylim(-0.3,1.7)
plt.xlabel('$\it{t}$', size = 20)
plt.ylabel('$\it{\phi}$', size = 20)
plt.savefig(filepath, bbox_inches="tight", pad_inches=0.0, transparent=True)

#%%
############################################
## tangent  ################################
############################################

# initial set

P = 1e0
Delta = -0.1
Initial = 1/4
dt = 2.5e-5 # step width

figpath = "fig/VAN/tangent/"
os.makedirs(figpath, exist_ok=True) # make folder

############
## X #######
############

filename = "X.pdf"
filepath = figpath + filename

plt.figure(0)

# limitcycle
datapath = "data/VAN/Floquet/"
with open(datapath + 'X0_.dat','rb') as f:
    X0_ = np.load(f, allow_pickle=True)

# simple
datapath = "data/VAN/simple/P{}Delta{}Initial{}/".format(P,Delta,Initial)
with open(datapath + 'X.dat','rb') as f:
    VAN_simple = np.load(f, allow_pickle=True)

# penalty
k = 0.5
datapath = "data/VAN/penalty/P{}Delta{}Initial{}k{}/".format(P,Delta,Initial,k)
with open(datapath + 'X.dat','rb') as f:
    VAN_penalty = np.load(f, allow_pickle=True)
    
k = 1
datapath = "data/VAN/penalty/P{}Delta{}Initial{}k{}/".format(P,Delta,Initial,k)
with open(datapath + 'X.dat','rb') as f:
    VAN_penalty2 = np.load(f, allow_pickle=True)

k = 5
datapath = "data/VAN/penalty/P{}Delta{}Initial{}k{}/".format(P,Delta,Initial,k)
with open(datapath + 'X.dat','rb') as f:
    VAN_penalty3 = np.load(f, allow_pickle=True)

#tangent
datapath = "data/VAN/tangent/P{}Delta{}Initial{}/".format(P,Delta,Initial)
with open(datapath + 'X.dat','rb') as f:
    VAN_tangent = np.load(f, allow_pickle=True)

plt.xticks([-1, 0, 1])

plt.xlabel(r"$x$", size = 20)
plt.ylabel(r"$y$", size = 20)

plt.plot(X0_[0,:],X0_[1,:], label='limit cycle', linestyle='--', color='r')
plt.plot(VAN_simple[0,:],VAN_simple[1,:], label='without penalty', color='orange')
plt.plot(VAN_penalty[0,:],VAN_penalty[1,:], color='brown', label='with penalty k=0.5', alpha = 0.8)
plt.plot(VAN_penalty2[0,:],VAN_penalty2[1,:], color='m', label='with penalty k=1', alpha = 0.6)
plt.plot(VAN_penalty3[0,:],VAN_penalty3[1,:], color='b', label='with penalty k=5')
plt.plot(VAN_tangent[0,:],VAN_tangent[1,:], label='tangent only', color = 'g', alpha=0.6)

plt.legend(fontsize=12)
plt.savefig(filepath, bbox_inches="tight", pad_inches=0.0, transparent=True)


###########
## theta ##
###########
filename = "theta.pdf"
filepath = figpath + filename

plt.figure(1)

# simple
datapath = "data/VAN/simple/P{}Delta{}Initial{}/".format(P,Delta,Initial)
with open(datapath + 'theta.dat','rb') as f:
    VAN_simple = np.load(f, allow_pickle=True)

# penalty
k = 0.5
datapath = "data/VAN/penalty/P{}Delta{}Initial{}k{}/".format(P,Delta,Initial,k)
with open(datapath + 'theta.dat','rb') as f:
    VAN_penalty = np.load(f, allow_pickle=True)
    
k = 1
datapath = "data/VAN/penalty/P{}Delta{}Initial{}k{}/".format(P,Delta,Initial,k)
with open(datapath + 'theta.dat','rb') as f:
    VAN_penalty2 = np.load(f, allow_pickle=True)

k = 5
datapath = "data/VAN/penalty/P{}Delta{}Initial{}k{}/".format(P,Delta,Initial,k)
with open(datapath + 'theta.dat','rb') as f:
    VAN_penalty3 = np.load(f, allow_pickle=True)

#tangent
datapath = "data/VAN/tangent/P{}Delta{}Initial{}/".format(P,Delta,Initial)
with open(datapath + 'theta.dat','rb') as f:
    VAN_tangent = np.load(f, allow_pickle=True)


Zeroline = np.array([[0,400],[0,0]])
plt.plot(VAN_simple[0,:],VAN_simple[1,:], label='without penalty', color='orange')
plt.plot(VAN_penalty[0,:],VAN_penalty[1,:], color='brown', label='with penalty k=0.5', alpha = 0.8)
plt.plot(VAN_penalty2[0,:],VAN_penalty2[1,:], color='m', label='with penalty k=1', alpha = 0.6)
plt.plot(VAN_penalty3[0,:],VAN_penalty3[1,:], color='b', label='with penalty k=5')
plt.plot(VAN_tangent[0,:],VAN_tangent[1,:], label='tangent only', color = 'g', alpha=0.6)
plt.plot(Zeroline[0,:],Zeroline[1,:], color='k', alpha=0.5, linestyle='--')

plt.legend()
plt.xlim(0,40)
plt.ylim(-0.3,1.7)
plt.xlabel('$\it{t}$', size = 20)
plt.ylabel('$\it{\phi}$', size = 20)
plt.savefig(filepath, bbox_inches="tight", pad_inches=0.0, transparent=True)

#######################
## theta and Average ##
#######################
filename = "theta_Average.pdf"
filepath = figpath + filename

plt.figure(2)

# simple
datapath = "data/VAN/simple/P{}Delta{}Initial{}/".format(P,Delta,Initial)
with open(datapath + 'theta.dat','rb') as f:
    VAN_simple = np.load(f, allow_pickle=True)
    
datapath = "data/VAN/simple/P{}Delta{}Initial{}/".format(P,Delta,Initial)
with open(datapath + 'Average.dat','rb') as f:
    VAN_simple_a = np.load(f, allow_pickle=True)

# penalty
k = 0.5
datapath = "data/VAN/penalty/P{}Delta{}Initial{}k{}/".format(P,Delta,Initial,k)
with open(datapath + 'theta.dat','rb') as f:
    VAN_penalty = np.load(f, allow_pickle=True)

datapath = "data/VAN/penalty/P{}Delta{}Initial{}k{}/".format(P,Delta,Initial,k)
with open(datapath + 'Average.dat','rb') as f:
    VAN_penalty_a = np.load(f, allow_pickle=True)
    
k = 1
datapath = "data/VAN/penalty/P{}Delta{}Initial{}k{}/".format(P,Delta,Initial,k)
with open(datapath + 'theta.dat','rb') as f:
    VAN_penalty2 = np.load(f, allow_pickle=True)
    
datapath = "data/VAN/penalty/P{}Delta{}Initial{}k{}/".format(P,Delta,Initial,k)
with open(datapath + 'Average.dat','rb') as f:
    VAN_penalty2_a = np.load(f, allow_pickle=True)

k = 5
datapath = "data/VAN/penalty/P{}Delta{}Initial{}k{}/".format(P,Delta,Initial,k)
with open(datapath + 'theta.dat','rb') as f:
    VAN_penalty3 = np.load(f, allow_pickle=True)
    
datapath = "data/VAN/penalty/P{}Delta{}Initial{}k{}/".format(P,Delta,Initial,k)
with open(datapath + 'Average.dat','rb') as f:
    VAN_penalty3_a = np.load(f, allow_pickle=True)

#tangent
datapath = "data/VAN/tangent/P{}Delta{}Initial{}/".format(P,Delta,Initial)
with open(datapath + 'theta.dat','rb') as f:
    VAN_tangent = np.load(f, allow_pickle=True)

datapath = "data/VAN/tangent/P{}Delta{}Initial{}/".format(P,Delta,Initial)
with open(datapath + 'Average.dat','rb') as f:
    VAN_tangent_a = np.load(f, allow_pickle=True)


Zeroline = np.array([[0,400],[0,0]])
plt.plot(VAN_simple[0,:],VAN_simple[1,:], linestyle = '--', color='orange', alpha = 0.4)
plt.plot(VAN_penalty[0,:],VAN_penalty[1,:], color='brown',linestyle = '--', alpha = 0.4)
plt.plot(VAN_penalty2[0,:],VAN_penalty2[1,:], color='m',linestyle = '--', alpha = 0.4)
plt.plot(VAN_penalty3[0,:],VAN_penalty3[1,:], color='b',linestyle = '--', alpha = 0.4)
plt.plot(VAN_tangent[0,:],VAN_tangent[1,:], linestyle = '--', color = 'g', alpha=0.4)
plt.plot(VAN_simple_a[0,:],VAN_simple_a[1,:], label='without penalty', linestyle = '-', color='orange', linewidth=2.0)
plt.plot(VAN_penalty_a[0,:],VAN_penalty_a[1,:], color='brown',linestyle = '-', label='with penalty k=0.5', linewidth=2.0)
plt.plot(VAN_penalty2_a[0,:],VAN_penalty2_a[1,:], color='m',linestyle = '-', label='with penalty k=1', linewidth=2.0)
plt.plot(VAN_penalty3_a[0,:],VAN_penalty3_a[1,:], color='b',linestyle = '-', label='with penalty k=5', linewidth=2.0)
plt.plot(VAN_tangent_a[0,:],VAN_tangent_a[1,:], label='tangent only',linestyle = '-', color = 'g', linewidth=2.0)




plt.plot(Zeroline[0,:],Zeroline[1,:], color='k', alpha=0.5, linestyle='--')


plt.legend(fontsize=14)
plt.xlim(0,40)
plt.ylim(-0.2,1.7)
plt.xlabel('$\it{t}$', size = 20)
plt.ylabel('$\it{\phi}$', size = 20)
plt.savefig(filepath, bbox_inches="tight", pad_inches=0.0, transparent=True)

##############
###  Gamma ###
##############

filename = "Gamma.pdf"
filepath = figpath + filename

plt.figure(3)

# simple
datapath = "data/VAN/simple/P{}Delta{}Initial{}/".format(P,Delta,Initial)
with open(datapath + 'Gamma.dat','rb') as f:
    VAN_simple = np.load(f, allow_pickle=True)

# penalty
k = 0.5
datapath = "data/VAN/penalty/P{}Delta{}Initial{}k{}/".format(P,Delta,Initial,k)
with open(datapath + 'Gamma.dat','rb') as f:
    VAN_penalty = np.load(f, allow_pickle=True)
    
k = 1
datapath = "data/VAN/penalty/P{}Delta{}Initial{}k{}/".format(P,Delta,Initial,k)
with open(datapath + 'Gamma.dat','rb') as f:
    VAN_penalty2 = np.load(f, allow_pickle=True)

k = 5
datapath = "data/VAN/penalty/P{}Delta{}Initial{}k{}/".format(P,Delta,Initial,k)
with open(datapath + 'Gamma.dat','rb') as f:
    VAN_penalty3 = np.load(f, allow_pickle=True)

#tangent
datapath = "data/VAN/tangent/P{}Delta{}Initial{}/".format(P,Delta,Initial)
with open(datapath + 'Gamma.dat','rb') as f:
    VAN_tangent = np.load(f, allow_pickle=True)


Zeroline = np.array([[-pi,pi],[0,0]])
Zeroline2 = np.array([[0,0],[-1,1]])

plt.plot(VAN_simple[0,:],VAN_simple[1,:], label='without penalty', color='orange')
plt.plot(VAN_penalty[0,:],VAN_penalty[1,:], color='brown', label='with penalty k=0.5', alpha = 0.8)
plt.plot(VAN_penalty2[0,:],VAN_penalty2[1,:], color='m', label='with penalty k=1', alpha = 0.6)
plt.plot(VAN_penalty3[0,:],VAN_penalty3[1,:], color='b', label='with penalty k=5')
plt.plot(VAN_tangent[0,:],VAN_tangent[1,:], label='tangent only', color = 'g', alpha=0.6)
plt.plot(Zeroline[0,:],Zeroline[1,:], color='k', alpha=1.0, linestyle='--')
plt.plot(Zeroline2[0,:],Zeroline2[1,:], color='k', alpha=1.0, linestyle='--')

plt.xlim(-pi,pi)
plt.ylim(-1.0,0.8)
plt.xticks([-pi, -pi/2, 0, pi/2, pi], ["-π", "-π/2", "0", "π/2", "π"])

plt.xlabel('$\it{\phi}$', size=20)
plt.ylabel('$\it{\Delta + \Gamma(\phi)}$', size=20)

plt.legend(fontsize=12)
plt.savefig(filepath, bbox_inches="tight", pad_inches=0.0, transparent=True)

###########
## input ##
###########
filename = "inputX.pdf"
filepath = figpath + filename

plt.figure(4)

# simple
datapath = "data/VAN/simple/P{}Delta{}Initial{}/".format(P,Delta,Initial)
with open(datapath + 'input.dat','rb') as f:
    VAN_simple = np.load(f, allow_pickle=True)

# penalty
k = 0.5
datapath = "data/VAN/penalty/P{}Delta{}Initial{}k{}/".format(P,Delta,Initial,k)
with open(datapath + 'input.dat','rb') as f:
    VAN_penalty = np.load(f, allow_pickle=True)
    
k = 1
datapath = "data/VAN/penalty/P{}Delta{}Initial{}k{}/".format(P,Delta,Initial,k)
with open(datapath + 'input.dat','rb') as f:
    VAN_penalty2 = np.load(f, allow_pickle=True)

k = 5
datapath = "data/VAN/penalty/P{}Delta{}Initial{}k{}/".format(P,Delta,Initial,k)
with open(datapath + 'input.dat','rb') as f:
    VAN_penalty3 = np.load(f, allow_pickle=True)

#tangent
datapath = "data/VAN/tangent/P{}Delta{}Initial{}/".format(P,Delta,Initial)
with open(datapath + 'input.dat','rb') as f:
    VAN_tangent = np.load(f, allow_pickle=True)

Te = VAN_simple.shape[1]*dt
Time_ax = np.linspace(0,VAN_simple.shape[1]-1, VAN_simple.shape[1])*dt

plt.plot(Time_ax, VAN_simple[0,:], label='without penalty', color='orange')
plt.plot(Time_ax, VAN_penalty[0,:], color='brown', label='with penalty k=0.5', alpha = 0.8)
plt.plot(Time_ax, VAN_penalty2[0,:], color='m', label='with penalty k=1', alpha = 0.6)
plt.plot(Time_ax, VAN_penalty3[0,:], color='b', label='with penalty k=5')
plt.plot(Time_ax, VAN_tangent[0,:], label='tangent only', color = 'g', alpha=0.6)


plt.legend(bbox_to_anchor=(0.6, 0.41))
plt.xlim(0,Te)
plt.xticks([0, Te/4, Te/2, 3*Te/4, Te], [ "0", "Te/4", "Te/2", "3Te/4", "Te"])
#plt.ylim(-0.3,1.7)
plt.xlabel('$\it{t}$', size = 20)
plt.ylabel('$\it{q_x}$', size = 20)
plt.savefig(filepath, bbox_inches="tight", pad_inches=0.0, transparent=True)
 

filename = "inputY.pdf"
filepath = figpath + filename

plt.figure(5)

# simple
datapath = "data/VAN/simple/P{}Delta{}Initial{}/".format(P,Delta,Initial)
with open(datapath + 'input.dat','rb') as f:
    VAN_simple = np.load(f, allow_pickle=True)

# penalty
k = 0.5
datapath = "data/VAN/penalty/P{}Delta{}Initial{}k{}/".format(P,Delta,Initial,k)
with open(datapath + 'input.dat','rb') as f:
    VAN_penalty = np.load(f, allow_pickle=True)
    
k = 1
datapath = "data/VAN/penalty/P{}Delta{}Initial{}k{}/".format(P,Delta,Initial,k)
with open(datapath + 'input.dat','rb') as f:
    VAN_penalty2 = np.load(f, allow_pickle=True)

k = 5
datapath = "data/VAN/penalty/P{}Delta{}Initial{}k{}/".format(P,Delta,Initial,k)
with open(datapath + 'input.dat','rb') as f:
    VAN_penalty3 = np.load(f, allow_pickle=True)

#tangent
datapath = "data/VAN/tangent/P{}Delta{}Initial{}/".format(P,Delta,Initial)
with open(datapath + 'input.dat','rb') as f:
    VAN_tangent = np.load(f, allow_pickle=True)

Te = VAN_simple.shape[1]*dt
Time_ax = np.linspace(0,VAN_simple.shape[1]-1, VAN_simple.shape[1])*dt

plt.plot(Time_ax, VAN_simple[1,:], label='without penalty', color='orange')
plt.plot(Time_ax, VAN_penalty[1,:], color='brown', label='with penalty k=0.5', alpha = 0.8)
plt.plot(Time_ax, VAN_penalty2[1,:], color='m', label='with penalty k=1', alpha = 0.6)
plt.plot(Time_ax, VAN_penalty3[1,:], color='b', label='with penalty k=5')
plt.plot(Time_ax, VAN_tangent[1,:], label='tangent only', color = 'g', alpha=0.6)


plt.legend()
plt.xlim(0,Te)
plt.xticks([0, Te/4, Te/2, 3*Te/4, Te], [ "0", "Te/4", "Te/2", "3Te/4", "Te"])
#plt.ylim(-0.3,1.7)
plt.xlabel('$\it{t}$', size = 20)
plt.ylabel('$\it{q_y}$', size = 20)
plt.savefig(filepath, bbox_inches="tight", pad_inches=0.0, transparent=True)
 

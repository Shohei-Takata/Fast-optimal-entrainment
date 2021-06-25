import numpy as np
import nwlib as my
import matplotlib.pyplot as plt
import os

pi = np.pi

# simulation parameter
tmax = 2000 # running time for convergence
dt = 2.5e-3 # step width

ND = my.ND(tmax, dt)

# FHN der Pol
#di/dt = x - x*x*x/3.0 - y + self.x0
#dv/dt = (x - self.b*y + self.y0) / self.tau
tau = 5
a = -1.0
b = 0.8
x0 = 0.5
y0 = 0.7
FHN = my.FHN(a, b, tau, x0, y0)

# inital point
Xstart = np.ones((2,1))

print("calculating X to convergence...", end="")

# run to convergence
Xstart = ND.evol_to_convergence(FHN.dif, Xstart)
print("converged X = ", Xstart.T)

# search period

# The phase is 0 when the y falls below y_basis = 0.0
Tnum, Xstart = ND.find_Tnum_Y(FHN.dif, Xstart, 0.0)
T = dt * Tnum
omega = 2*pi/T

print("period T = ", T)
print("frequency omega = ", omega)

print("calculate limit cycle X0_ ...", end="")

# calculate floquet vector
Floquet = my.Floquet(FHN.dif, FHN.Jacobian, Tnum, T, omega, dt)

# limit cycle X0_ and the first floquet vector u0
X0_, u0_ = Floquet.Calc_X0_u0(Xstart)

print("done")
print("calculate v0_ ...", end="")

# the first left floquet vector v0_ and v0_dif = d(v0)/dt
# For convergence from any initial point of v0(0), we evol v0 for number_of_rotations = 2 times.
v0_, v0_dif = Floquet.Calc_v0(X0_, 10)

print("done")
print("calculate u1_ ...", end="")

# the second right floquet vector u1_ and floqet eigenvalue lambda1
# For convergence from any initial point of u0(0), we evol u0 for number_of_rotations = 2 times.
lambda1, u1_ = Floquet.Calc_u1(X0_, u0_, v0_, 2)

print("done")
print("calculate v1_ ...", end="")
# the second left floquet vector v1_
# For convergence from any initial point of v1(0), we evol v1 for number_of_rotations = 2 times.
v1_ = Floquet.Calc_v1(lambda1, X0_, u0_, v0_, u1_, 2)

print("done")
print("")
print("calculating floquet vector has finished")

#%%
###### plot #######################################################
## verificate production of floquet vector
fig = plt.figure(0, figsize=(8.0, 6.0))
ax1 = fig.add_subplot(2,2,1)
ax2 = fig.add_subplot(2,2,2)
ax3 = fig.add_subplot(2,2,3)
ax4 = fig.add_subplot(2,2,4)

Time_ax = np.linspace(0,Tnum-1, Tnum)*dt

fig.suptitle("Verification of floquet product")
# <u0,v0>
ax1.plot(Time_ax, Floquet.product(u0_, v0_))
ax1.set_title('<u0,v0>')
ax1.grid()
# <u0,v1>
ax2.plot(Time_ax, Floquet.product(u0_, v1_))
ax2.set_title('<u0,v1>')
ax2.grid()
# <u0,v2>
ax3.plot(Time_ax, Floquet.product(u1_, v0_))
ax3.set_title('<u1,v0>')
ax3.grid()
# <u0,v0>
ax4.plot(Time_ax, Floquet.product(u1_, v1_))
ax4.set_title('<u1,v1>')
ax4.grid()

## preparation for drawing and saving floquet vector #####################################

figpath = "fig/FHN/Floquet/"
os.makedirs(figpath, exist_ok=True) # make folder

datapath = "data/FHN/Floquet/"
os.makedirs(datapath, exist_ok=True) # make folder

## u0_ 
filename = "u0_.pdf"
filepath = figpath + filename

plt.figure(1)
plt.plot(Time_ax, u0_[0,:], label = r"$u_{0x}$", color = 'r')
plt.plot(Time_ax, u0_[1,:], linestyle = '--', label = r"$u_{0y}$", color = 'g')
plt.xlabel(r"$\theta$")
plt.ylabel(r"$u_0$", size = 20)

plt.xticks([0, T/4, T/2, 3*T/4, T], [ "0", "π/2", "π", "3π/2", "2π"])
#plt.xticks([0, T/4, T/2, 3*T/4, T], [ "0", "T/4", "T/2", "3T/4", "T"])
plt.xlim(0,T)

plt.legend(fontsize=17)
plt.savefig(filepath, bbox_inches="tight", pad_inches=0.0, transparent=True)

## u1_
filename = "u1_.pdf"
filepath = figpath + filename

plt.figure(2)
plt.plot(Time_ax, u1_[0,:], label = r"$u_{1x}$", color = 'b')
plt.plot(Time_ax, u1_[1,:], linestyle = '--', label = r"$u_{1y}$", color = 'orange')
plt.xlabel(r"$\theta$")
plt.ylabel(r"$u_1$", size = 20)

plt.xticks([0, T/4, T/2, 3*T/4, T], [ "0", "π/2", "π", "3π/2", "2π"])
plt.xlim(0,T)

plt.legend(fontsize=17)
plt.savefig(filepath, bbox_inches="tight", pad_inches=0.0, transparent=True)

## v0_
filename = "v0_.pdf"
filepath = figpath + filename

plt.figure(3)
plt.plot(Time_ax, v0_[0,:], label = r"$v_{0x}$", color = 'r')
plt.plot(Time_ax, v0_[1,:], linestyle = '--', label = r"$v_{0y}$", color = 'g')
plt.xlabel(r"$\theta$")
plt.ylabel(r"$v_0$", size = 20)

plt.xticks([0, T/4, T/2, 3*T/4, T], [ "0", "π/2", "π", "3π/2", "2π"])
plt.xlim(0,T)

plt.legend(fontsize=17)
plt.savefig(filepath, bbox_inches="tight", pad_inches=0.0, transparent=True)

## v1_
filename = "v1_.pdf"
filepath = figpath + filename

plt.figure(4)
plt.plot(Time_ax, v1_[0,:], label = r"$v_{1x}$", color = 'b')
plt.plot(Time_ax, v1_[1,:], linestyle = '--', label = r"$v_{1y}$", color = 'orange')
plt.xlabel(r"$\theta$")
plt.ylabel(r"$v_1$", size = 20)

plt.xticks([0, T/4, T/2, 3*T/4, T], [ "0", "π/2", "π", "3π/2", "2π"])
plt.xlim(0,T)

plt.legend(fontsize=17)
plt.savefig(filepath, bbox_inches="tight", pad_inches=0.0, transparent=True)

## save data
X0_.dump(datapath + 'X0_.dat')
u0_.dump(datapath + 'u0_.dat')
u1_.dump(datapath + 'u1_.dat')
v0_.dump(datapath + 'v0_.dat')
v1_.dump(datapath + 'v1_.dat')


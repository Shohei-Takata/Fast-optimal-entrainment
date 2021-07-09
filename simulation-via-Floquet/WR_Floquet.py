#%% Calculate Floquet(Willamowski-Rossler model)
import numpy as np
import nwlib as my
import matplotlib.pyplot as plt 
import os

pi = np.pi

# simulation parameter
tmax = 100 # running time for convergence
dt = 2.5e-6 # step width

ND = my.ND(tmax, dt)

# Willamowski-Rossler
#dx/dt = x*(b1 - d1*x - y - z) 
#dy/dt = y*(b2 - d2*y - x) 
#dz/dt = z*(x - d3)
b1 = 80
b2 = 20
d1 = 0.16
d2 = 0.13
d3 = 16
WR = my.WR(b1, b2, d1, d2, d3)

# inital point
Xstart = np.ones((3,1))

print("calculating X to convergence...", end="")

# run to convergence
Xstart = ND.evol_to_convergence(WR.dif, Xstart)
print("converged X = ", Xstart.T)

# search period

# The phase is 0 when the y falls below y_basis = 0.0
Tnum, Xstart = ND.find_Tnum_X(WR.dif, Xstart, x_basis=15)
T = dt * Tnum
omega = 2*pi/T

print("period T = ", T)
print("frequency omega = ", omega)

print("calculate limit cycle X0_ ...", end="")

# calculate floquet vector
Floquet = my.Floquet(WR.dif, WR.Jacobian, Tnum, T, omega, dt)

# limit cycle X0_ and the first floquet vector u0
X0_, u0_ = Floquet.Calc_X0_u0(Xstart)

print("done")
print("calculate v0_ ...", end="")

# the first left floquet vector v0_ and v0_dif = d(v0)/dt
# For convergence from any initial point of v0(0), we evol v0 for rotations = 40 times.
v0_, v0_dif = Floquet.Calc_v0(X0_, 40)

print("done")
print("calculate u1_u2_ ...", end="")

# Conjugated floquet vectors u1_ and u2_

# First, we find initial Floquet eigenvectors u1(0), u2(0) 
# and Floquet eigenvalues lambda1, lambda2 via Monodromy matrix. 
M = Floquet.monodromy(X0_) # monodromy matrix
EigenValue , EigenVector = np.linalg.eig(M)
u1 = np.copy(EigenVector[:,1:2]) # u1(0)
u2 = np.copy(EigenVector[:,2:3]) # u2(0)

lambda1 = np.log(EigenValue[1]) / T # Floquet eigenvalue
lambda2 = np.log(EigenValue[2]) / T

u1_, u2_ = Floquet.Calc_u1u2(lambda1, lambda2, u1, u2, X0_, u0_, v0_)

print("done")
print("Floquet eigenvalue lambda1 = ", lambda1)

print("calculate v1_v2_ ...", end="")
# Conjugated floquet vectors v1_ and v2_

# First, we find initial Floquet eigenvectors v1(0), v2(0) 
# and Floquet eigenvalues lambda1, lambda2 via Monodromy matrix. 

EigenValue , EigenVector = np.linalg.eig(M.T) 

v1 = np.copy(EigenVector[:,2:3]) # v1(0)
v2 = np.copy(EigenVector[:,1:2]) 

v1 = v1 / np.dot(np.conjugate(u1).T,v1) # normalization for <u1, v1> = 1
v2 = v2 / np.dot(np.conjugate(u2).T,v2)

v1_, v2_ = Floquet.Calc_v1v2(lambda1, lambda2, v1, v2, X0_, u0_, v0_)

print("done")
print("")
print("calculating floquet vector has finished")

#%%
###### plot #######################################################
## verificate production of floquet vector
fig = plt.figure(0, figsize=(18, 12))
ax = fig.subplots(3, 3)
Time_ax = np.linspace(0,Tnum-1, Tnum)*dt

# 左上 <u0,v0>
ax[0, 0].plot(Time_ax, my.inner_product(u0_, v0_, Tnum))
ax[0, 0].set_title('<u0,v0>')
ax[0, 0].grid()
# 中央上 <u0,v1>
ax[0, 1].plot(Time_ax, my.inner_product(u0_, v1_, Tnum))
ax[0, 1].set_title('<u0,v1>')
ax[0, 1].grid()
# 右上 <u0,v2>
ax[0, 2].plot(Time_ax, my.inner_product(u0_, v2_, Tnum))
ax[0, 2].set_title('<u0,v2>')
ax[0, 2].grid()
# 左中央 <u0,v0>
ax[1, 0].plot(Time_ax, my.inner_product(u1_, v0_, Tnum))
ax[1, 0].set_title('<u1,v0>')
ax[1, 0].grid()
# 中央 <u0,v1>
ax[1, 1].plot(Time_ax, my.inner_product(u1_, v1_, Tnum))
ax[1, 1].set_title('<u1,v1>')
ax[1, 1].grid()
# 右中央 <u0,v2>
ax[1, 2].plot(Time_ax, my.inner_product(u1_, v2_, Tnum))
ax[1, 2].set_title('<u1,v2>')
ax[1, 2].grid()
# 左下 <u2,v0>
ax[2, 0].plot(Time_ax, my.inner_product(u2_, v0_, Tnum))
ax[2, 0].set_title('<u2,v0>')
ax[2, 0].grid()
# 中央上 <u2,v1>
ax[2, 1].plot(Time_ax, my.inner_product(u2_, v1_, Tnum))
ax[2, 1].set_title('<u2,v1>')
ax[2, 1].grid()
# 右上 <u2,v2>
ax[2, 2].plot(Time_ax, my.inner_product(u2_, v2_, Tnum))
ax[2, 2].set_title('<u2,v2>')
ax[2, 2].grid()

## preparation for drawing and saving floquet vector #####################################

figpath = "fig/WR/Floquet/"
os.makedirs(figpath, exist_ok=True) # make folder

datapath = "data/WR/Floquet/"
os.makedirs(datapath, exist_ok=True) # make folder

## u0_ 
filename = "u0_.pdf"
filepath = figpath + filename

plt.figure(1)
plt.plot(Time_ax, u0_[0,:], label = r"$u_{0x}$", color = 'r')
plt.plot(Time_ax, u0_[1,:], linestyle = '--', label = r"$u_{0y}$", color = 'g')
plt.plot(Time_ax, u0_[2,:], linestyle = 'dotted', label = r"$u_{0z}$", color = 'b')

plt.xlabel(r"$\theta$", size = 20)
plt.ylabel(r"$u_0$", size = 22)
#plt.xticks([0, T/4, T/2, 3*T/4, T], [ "0", "T/4", "T/2", "3T/4", "T"])
plt.xticks([0, T/4, T/2, 3*T/4, T], [ "0", "π/2", "π", "3π/2", "2π"])
plt.yticks([-100,0,100,200])
plt.xlim(0,T)
plt.ylim(-120,220)
plt.legend(fontsize=17)
plt.savefig(filepath, bbox_inches="tight", pad_inches=0.0, transparent=True)

## u1_(real)
filename = "Re[u1_].pdf"
filepath = figpath + filename

plt.figure(2)
plt.plot(Time_ax, u1_[0,:].real, label = r"Re$~u_{1x}$", color = 'r')
plt.plot(Time_ax, u1_[1,:].real, linestyle = '--', label = r"Re$~u_{1y}$", color = 'g')
plt.plot(Time_ax, u1_[2,:].real, linestyle = 'dotted', label = r"Re$~u_{1z}$", color = 'b')

plt.xlabel(r"$\theta$", size = 20)
plt.ylabel(r"Re$~u_1$", size = 20)
#plt.xticks([0, T/4, T/2, 3*T/4, T], [ "0", "T/4", "T/2", "3T/4", "T"])
plt.xticks([0, T/4, T/2, 3*T/4, T], [ "0", "π/2", "π", "3π/2", "2π"])
plt.xlim(0,T)
plt.ylim(-3,5)
plt.legend(fontsize=17)
plt.savefig(filepath, bbox_inches="tight", pad_inches=0.0, transparent=True)

## u1_(imag)
filename = "Im[u1_].pdf"
filepath = figpath + filename

plt.figure(3)
plt.plot(Time_ax, u1_[0,:].imag, label = r"Im$~u_{1x}$", color = 'r')
plt.plot(Time_ax, u1_[1,:].imag, linestyle = '--', label = r"Im$~u_{1y}$", color = 'g')
plt.plot(Time_ax, u1_[2,:].imag, linestyle = 'dotted', label = r"Im$~u_{1z}$", color = 'b')

plt.xlabel(r"$\theta$", size = 20)
plt.ylabel(r"Im$~u_1$", size = 20)
#plt.xticks([0, T/4, T/2, 3*T/4, T], [ "0", "T/4", "T/2", "3T/4", "T"])
plt.xticks([0, T/4, T/2, 3*T/4, T], [ "0", "π/2", "π", "3π/2", "2π"])
plt.xlim(0,T)
plt.ylim(-11, 7)
plt.legend(fontsize=15)
plt.savefig(filepath, bbox_inches="tight", pad_inches=0.0, transparent=True)

## v0_
filename = "v0_.pdf"
filepath = figpath + filename

plt.figure(4)
plt.plot(Time_ax, v0_[0,:], label = r"$v_{0x}$", color = 'r')
plt.plot(Time_ax, v0_[1,:], linestyle = '--', label = r"$v_{0y}$", color = 'g')
plt.plot(Time_ax, v0_[2,:], linestyle = 'dotted', label = r"$v_{0z}$", color = 'b')

plt.xlabel(r"$\theta$", size = 20)
plt.ylabel(r"$v_0$", size = 20)

plt.xticks([0, T/4, T/2, 3*T/4, T], [ "0", "π/2", "π", "3π/2", "2π"])
plt.xlim(0,T)
plt.ylim(-1.2,0.6)
plt.legend(fontsize=17)
plt.savefig(filepath, bbox_inches="tight", pad_inches=0.0, transparent=True)

## v1_(real)
filename = "Re[v1_].pdf"
filepath = figpath + filename

plt.figure(5)
plt.plot(Time_ax, v1_[0,:].real, label = r"Re$~v_{1x}$", color = 'r')
plt.plot(Time_ax, v1_[1,:].real, linestyle = '--', label = r"Re$~v_{1y}$", color = 'g')
plt.plot(Time_ax, v1_[2,:].real, linestyle = 'dotted', label = r"Re$~v_{1z}$", color = 'b')
plt.xlabel(r"$\theta$", size = 20)
plt.ylabel(r"Re$~v_1$", size = 20)
#plt.xticks([0, T/4, T/2, 3*T/4, T], [ "0", "T/4", "T/2", "3T/4", "T"])
plt.xticks([0, T/4, T/2, 3*T/4, T], [ "0", "π/2", "π", "3π/2", "2π"])
plt.xlim(0,T)
plt.ylim(-55, 46)
plt.legend(fontsize=15)
plt.savefig(filepath, bbox_inches="tight", pad_inches=0.0, transparent=True)

## v1_(imag)
filename = "Im[v1_].pdf"
filepath = figpath + filename

plt.figure(6)
plt.plot(Time_ax, v1_[0,:].imag, label = r"Im$~v_{1x}$", color = 'r')
plt.plot(Time_ax, v1_[1,:].imag, linestyle = '--', label = r"Im$~v_{1y}$", color = 'g')
plt.plot(Time_ax, v1_[2,:].imag, linestyle = 'dotted', label = r"Im$~v_{1z}$", color = 'b')
plt.xlabel(r"$\theta$", size = 20)
plt.ylabel(r"Im$~v_1$", size = 20)
#plt.xticks([0, T/4, T/2, 3*T/4, T], [ "0", "T/4", "T/2", "3T/4", "T"])
plt.xticks([0, T/4, T/2, 3*T/4, T], [ "0", "π/2", "π", "3π/2", "2π"])
plt.yticks([0,-10,-20])
plt.xlim(0,T)
plt.ylim(-27, 7)
plt.legend(fontsize=17)
plt.savefig(filepath, bbox_inches="tight", pad_inches=0.0, transparent=True)


## save data
X0_.dump(datapath + 'X0_.dat')
u0_.dump(datapath + 'u0_.dat')
u1_.dump(datapath + 'u1_.dat')
u2_.dump(datapath + 'u2_.dat')
v0_.dump(datapath + 'v0_.dat')
v1_.dump(datapath + 'v1_.dat')
v2_.dump(datapath + 'v2_.dat')

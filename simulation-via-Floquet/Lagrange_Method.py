import numpy as np
import matplotlib.pyplot as plt
import sys

##### Lagrange Method ###########################
pi = np.pi

class SL:
    def __init__(self, b, omega, P, Delta):
        self.b = b
        self.omega = omega
        self.P = P
        self.Delta = Delta
    
    def Calc_mu_nu(self):
        nu = self.omega**2*(1+self.b**2)/2 * np.sqrt(1/(self.omega**2*self.P*(1+self.b**2)-self.Delta**2))
        mu = -2*nu*self.Delta/(self.omega**2*(1+self.b**2))
        return mu, nu


class simple:
    def __init__(self, Tnum, v0_, v0_dif, omega, P, Delta):
        self.Tnum = Tnum
        self.v0_ = v0_
        self.v0_dif = v0_dif
        self.omega = omega
        self.P = P
        self.Delta = Delta
    
    # Z = v0_  I1 = v1_
    # Z'(θ) = 1/ω dZ/dt
    def Calc_mu_nu(self):
        Z2 = 0 # <Z, Z>
        Zdif2 = 0 # <Z', Z'>
        for tt in range(self.Tnum):
            v0 = self.v0_[:,tt:tt+1]
            v0dif = self.v0_dif[:,tt:tt+1]
            Z2 += np.dot(v0.T, v0) # <Z, Z>
            Zdif2 += np.dot(v0dif.T, v0dif) # <dZ/dt, dZ/dt>
        Z2 = Z2 / self.Tnum
        Zdif2 = Zdif2 / self.Tnum / self.omega / self.omega # <Z', Z'>
        if Zdif2 / (self.P - self.Delta**2 / Z2) < 0:
            print("Lagrange Error")
            print("cannot take the square root")
            sys.exit()
        nu = 1/2*np.sqrt( Zdif2 / (self.P - self.Delta**2 / Z2))
        mu = (- 2 * nu * self.Delta ) / Z2
        return mu.item(), nu.item()
    
    def Calc_mu_nu_arnold(self):
        Z2 = 0 # <Z, Z>
        Zdif2 = 0 # <Z', Z'>
        for tt in range(self.Tnum):
            v0 = self.v0_[:,tt:tt+1]
            v0dif = self.v0_dif[:,tt:tt+1]
            Z2 += np.dot(v0.T, v0) # <Z, Z>
            Zdif2 += np.dot(v0dif.T, v0dif) # <dZ/dt, dZ/dt>
        Z2 = Z2 / self.Tnum
        Zdif2 = Zdif2 / self.Tnum / self.omega / self.omega # <Z', Z'>
        if Zdif2 / (self.P - self.Delta**2 / Z2) < 0:
            return 0, 0 # not synchronization
        nu = 1/2*np.sqrt( Zdif2 / (self.P - self.Delta**2 / Z2))
        mu = (- 2 * nu * self.Delta ) / Z2
        return mu.item(), nu.item()
    
    def Calc_Gamma(self, mu, nu):
        # Γ(Φ) = ∫Z(Φ+ωt)・q(t)dt
        Division = 100 # Divide [-π, π] into 100 pieces.
        GammaPhi = np.empty((self.v0_.shape[0],Division+1))
        for tp in range(Division+1):
            Pnum = tp * int(self.Tnum/Division) - int(self.Tnum/2) # Φ(num) [-Tnum/2, Tnum/2]
            if(Pnum < 0):
                te = Pnum + self.Tnum # Φ = Pnum and Φ = Pnum + Tnum is equivalent for Z(Φ+ωt). 
            else:
                te = Pnum
            # Calculate Γ(Φ)
            Gamma = 0
            for tt in range(self.Tnum):
                # Z0(Φ+ωt)
                v0phi = np.array([self.v0_[:,(tt+te)%self.Tnum]]).T
                # q(ωt)
                v0 = self.v0_[:,tt:tt+1]
                v0dif = self.v0_dif[:,tt:tt+1]
                q_ = 1/2/nu * (-1/self.omega * v0dif + mu * v0)
                # Γ(Φ) = ∫Z0(Φ+ωt)・q(ωt)dt
                Gamma += np.dot(q_.T, v0phi)
            # save
            GammaPhi[0,tp] = Pnum/self.Tnum*2*pi
            GammaPhi[1,tp] = self.Delta + Gamma/self.Tnum
        # plot
        plt.plot(GammaPhi[0,:], GammaPhi[1,:])
        plt.xlabel("Φ")
        plt.ylabel("Δ+Γ(Φ)")
        plt.xlim(-pi,pi)
        plt.grid()
        return GammaPhi
    
class penalty_2D:
    def __init__(self, k, Tnum, v0_, v0_dif, v1_, omega, P, Delta):
        self.k = k
        self.Tnum = Tnum
        self.v0_ = v0_
        self.v0_dif = v0_dif
        self.v1_ = v1_
        self.omega = omega
        self.P = P
        self.Delta = Delta
    
    # Z = v0_  I1 = v1_
    # Z'(θ) = 1/ω dZ/dt
    def Calc_nu(self, nu):
        mu_above = 0
        mu_under = 0
        for tt in range(self.Tnum):
            v0 = self.v0_[:,tt:tt+1]
            v0dif = self.v0_dif[:,tt:tt+1]
            v1x = self.v1_[0,tt]
            v1y = self.v1_[1,tt]
            InverseMatrix = 1/(nu*nu + self.k *nu * (v1x*v1x + v1y*v1y))*np.array([[nu+self.k*v1y*v1y, -self.k*v1x*v1y],[-self.k*v1x*v1y, nu+self.k*v1x*v1x]])
            mu_above += np.dot(np.dot(v0.T, InverseMatrix),v0dif)/self.omega # numerator of mu, v'(θ) = 1/ω v0dif
            mu_under += np.dot(np.dot(v0.T, InverseMatrix),v0) # denominator of mu
        mu = (mu_above/self.Tnum - 2*self.Delta) / (mu_under/self.Tnum)
        Qsum = 0
        for tt in range(self.Tnum):
            v0 = self.v0_[:,tt:tt+1]
            v0dif = self.v0_dif[:,tt:tt+1]
            v1x = self.v1_[0,tt]
            v1y = self.v1_[1,tt]
            InverseMatrix = 1/(nu*nu + self.k *nu * (v1x*v1x + v1y*v1y))*np.array([[nu+self.k*v1y*v1y, -self.k*v1x*v1y],[-self.k*v1x*v1y, nu+self.k*v1x*v1x]])
            q_ = 1/2*np.dot(InverseMatrix, -v0dif/self.omega + mu * v0)
            Qsum += np.dot(q_.T,q_)
        Qsum = Qsum / self.Tnum
        return Qsum.item() - self.P
    
    def Calc_mu(self, nu):
        mu_above = 0
        mu_under = 0
        for tt in range(self.Tnum):
            v0 = self.v0_[:,tt:tt+1]
            v0dif = self.v0_dif[:,tt:tt+1]
            v1x = self.v1_[0,tt]
            v1y = self.v1_[1,tt]
            InverseMatrix = 1/(nu*nu + self.k *nu * (v1x*v1x + v1y*v1y))*np.array([[nu+self.k*v1y*v1y, -self.k*v1x*v1y],[-self.k*v1x*v1y, nu+self.k*v1x*v1x]])
            mu_above += np.dot(np.dot(v0.T, InverseMatrix),v0dif)/self.omega # numerator of mu, v'(θ) = 1/ω v0dif
            mu_under += np.dot(np.dot(v0.T, InverseMatrix),v0) # denominator of mu
        mu = (mu_above/self.Tnum - 2*self.Delta) / (mu_under/self.Tnum)
        return mu.item()
    
    def Calc_Gamma(self, mu, nu):
        # Γ(Φ)
        Division = 100 # Divide [-π, π] into 100 pieces.
        GammaPhi = np.empty((2,Division+1))
        for tp in range(Division+1):
            Pnum = tp * int(self.Tnum/Division) - int(self.Tnum/2) # Φ(num) [-Tnum/2, Tnum/2]
            if(Pnum < 0):
                te = Pnum + self.Tnum # Φ = Pnum and Φ = Pnum + Tnum is equivalent for Z(Φ+ωt). 
            else:
                te = Pnum
            # Calculate Γ(Φ)
            Gamma = 0
            for tt in range(self.Tnum):
                # Z0(Φ+ωt)
                v0phi = np.array([self.v0_[:,(tt+te)%self.Tnum]]).T
                # q(ωt)
                v0 = self.v0_[:,tt:tt+1]
                v0dif = self.v0_dif[:,tt:tt+1]
                v1x = self.v1_[0,tt]
                v1y = self.v1_[1,tt]
                InverseMatrix = 1/(nu*nu + self.k *nu * (v1x*v1x + v1y*v1y))*np.array([[nu+self.k*v1y*v1y, -self.k*v1x*v1y],[-self.k*v1x*v1y, nu+self.k*v1x*v1x]])
                q_ = 1/2*np.dot(InverseMatrix, -v0dif/self.omega + mu * v0)
                # Γ(Φ) = ∫Z0(Φ+ωt)・q(ωt)dt
                Gamma += np.dot(q_.T, v0phi)
            # save
            GammaPhi[0,tp] = Pnum/self.Tnum*2*pi
            GammaPhi[1,tp] = self.Delta + Gamma/self.Tnum
        # plot
        plt.plot(GammaPhi[0,:], GammaPhi[1,:])
        plt.xlabel("Φ")
        plt.ylabel("Δ+Γ(Φ)")
        plt.xlim(-pi,pi)
        plt.grid()
        return GammaPhi

class penalty_3D_conj:
    def __init__(self, k, Tnum, v0_, v0_dif, v1_, omega, P, Delta):
        self.k = k
        self.Tnum = Tnum
        self.v0_ = v0_
        self.v0_dif = v0_dif
        self.v1_ = v1_
        self.omega = omega
        self.P = P
        self.Delta = Delta
    
    # Z = v0_  I1 = v1_
    # Z'(θ) = 1/ω dZ/dt
    def Calc_nu(self, nu):
        mu_above = 0
        mu_under = 0
        for tt in range(self.Tnum):
            v0 = self.v0_[:,tt:tt+1]
            v0dif = self.v0_dif[:,tt:tt+1]
            v1x = self.v1_[0,tt]
            v1y = self.v1_[1,tt]
            v1z = self.v1_[2,tt]
            # Calculate inverse matrix 
            # a_ij is component of (νE + k v1 v1† + k v2 v2†)
            a11 = (nu+self.k*(v1x*v1x.conjugate()+v1x.conjugate()*v1x)).real # For calculation stability, we eliminate imaginary part. 
            a12 = self.k*(v1x*v1y.conjugate() + v1y*v1x.conjugate()).real
            a13 = self.k*(v1x*v1z.conjugate() + v1x.conjugate()*v1z).real
            a21 = a12
            a22 = nu+self.k*(v1y*v1y.conjugate()+v1y.conjugate()*v1y).real
            a23 = self.k*(v1y*v1z.conjugate() + v1y.conjugate()*v1z).real
            a31 = a13
            a32 = a23
            a33 = nu+self.k*(v1z*v1z.conjugate()+v1z.conjugate()*v1z).real
            InverseMatrix = 1/(a11*a22*a33+a12*a23*a31+a13*a21*a32-a13*a22*a31-a12*a21*a33-a11*a23*a32) * np.array([[a22*a33-a23*a32,-(a12*a33-a13*a32),a12*a23-a13*a22],[-(a21*a33-a23*a31), a11*a33-a13*a31, -(a11*a23-a13*a21)],[a21*a32-a22*a31, -(a11*a32-a12*a31), (a11*a22-a21*a12)]])
            mu_above += np.dot(np.dot(v0.T, InverseMatrix),v0dif)/self.omega # numerator of mu, v'(θ) = 1/ω v0dif
            mu_under += np.dot(np.dot(v0.T, InverseMatrix),v0) # denominator of mu
        mu = (mu_above/self.Tnum - 2*self.Delta) / (mu_under/self.Tnum)
        Qsum = 0 # <q, q>
        for tt in range(self.Tnum):
            v0 = self.v0_[:,tt:tt+1]
            v0dif = self.v0_dif[:,tt:tt+1]
            v1x = self.v1_[0,tt]
            v1y = self.v1_[1,tt]
            v1z = self.v1_[2,tt]
            # Calculate inverse matrix 
            # a_ij is component of (νE + k v1 v1† + k v2 v2†)
            a11 = (nu+self.k*(v1x*v1x.conjugate()+v1x.conjugate()*v1x)).real # For calculation stability, we eliminate imaginary part. 
            a12 = self.k*(v1x*v1y.conjugate() + v1y*v1x.conjugate()).real
            a13 = self.k*(v1x*v1z.conjugate() + v1x.conjugate()*v1z).real
            a21 = a12
            a22 = nu+self.k*(v1y*v1y.conjugate()+v1y.conjugate()*v1y).real
            a23 = self.k*(v1y*v1z.conjugate() + v1y.conjugate()*v1z).real
            a31 = a13
            a32 = a23
            a33 = nu+self.k*(v1z*v1z.conjugate()+v1z.conjugate()*v1z).real
            InverseMatrix = 1/(a11*a22*a33+a12*a23*a31+a13*a21*a32-a13*a22*a31-a12*a21*a33-a11*a23*a32) * np.array([[a22*a33-a23*a32,-(a12*a33-a13*a32),a12*a23-a13*a22],[-(a21*a33-a23*a31), a11*a33-a13*a31, -(a11*a23-a13*a21)],[a21*a32-a22*a31, -(a11*a32-a12*a31), (a11*a22-a21*a12)]])
            q_ = 1/2*np.dot(InverseMatrix, -v0dif/self.omega + mu * v0)
            Qsum += np.dot(q_.T,q_)
        Qsum = Qsum / self.Tnum
        return Qsum.item() - self.P
    
    def Calc_mu(self, nu):
        mu_above = 0
        mu_under = 0
        for tt in range(self.Tnum):
            v0 = self.v0_[:,tt:tt+1]
            v0dif = self.v0_dif[:,tt:tt+1]
            v1x = self.v1_[0,tt]
            v1y = self.v1_[1,tt]
            v1z = self.v1_[2,tt]
            # Calculate inverse matrix 
            # a_ij is component of (νE + k v1 v1† + k v2 v2†)
            a11 = (nu+self.k*(v1x*v1x.conjugate()+v1x.conjugate()*v1x)).real # For calculation stability, we eliminate imaginary part. 
            a12 = self.k*(v1x*v1y.conjugate() + v1y*v1x.conjugate()).real
            a13 = self.k*(v1x*v1z.conjugate() + v1x.conjugate()*v1z).real
            a21 = a12
            a22 = nu+self.k*(v1y*v1y.conjugate()+v1y.conjugate()*v1y).real
            a23 = self.k*(v1y*v1z.conjugate() + v1y.conjugate()*v1z).real
            a31 = a13
            a32 = a23
            a33 = nu+self.k*(v1z*v1z.conjugate()+v1z.conjugate()*v1z).real
            InverseMatrix = 1/(a11*a22*a33+a12*a23*a31+a13*a21*a32-a13*a22*a31-a12*a21*a33-a11*a23*a32) * np.array([[a22*a33-a23*a32,-(a12*a33-a13*a32),a12*a23-a13*a22],[-(a21*a33-a23*a31), a11*a33-a13*a31, -(a11*a23-a13*a21)],[a21*a32-a22*a31, -(a11*a32-a12*a31), (a11*a22-a21*a12)]])
            mu_above += np.dot(np.dot(v0.T, InverseMatrix),v0dif)/self.omega # numerator of mu, v'(θ) = 1/ω v0dif
            mu_under += np.dot(np.dot(v0.T, InverseMatrix),v0) # denominator of mu
        mu = (mu_above/self.Tnum - 2*self.Delta) / (mu_under/self.Tnum)
        return mu.item()
    
    def Calc_Gamma(self, mu, nu):
        # Γ(Φ)
        Division = 100 # Divide [-π, π] into 100 pieces.
        GammaPhi = np.empty((3,Division+1))
        for tp in range(Division+1):
            Pnum = tp * int(self.Tnum/Division) - int(self.Tnum/2) # Φ(num) [-Tnum/2, Tnum/2]
            if(Pnum < 0):
                te = Pnum + self.Tnum # Φ = Pnum and Φ = Pnum + Tnum is equivalent for Z(Φ+ωt). 
            else:
                te = Pnum
            # Calculate Γ(Φ)
            Gamma = 0
            for tt in range(self.Tnum):
                # Z0(Φ+ωt)
                v0phi = np.array([self.v0_[:,(tt+te)%self.Tnum]]).T
                # q(ωt)
                v0 = self.v0_[:,tt:tt+1]
                v0dif = self.v0_dif[:,tt:tt+1]
                v1x = self.v1_[0,tt]
                v1y = self.v1_[1,tt]
                v1z = self.v1_[2,tt]
                # Calculate inverse matrix 
                # a_ij is component of (νE + k v1 v1† + k v2 v2†)
                a11 = (nu+self.k*(v1x*v1x.conjugate()+v1x.conjugate()*v1x)).real # For calculation stability, we eliminate imaginary part. 
                a12 = self.k*(v1x*v1y.conjugate() + v1y*v1x.conjugate()).real
                a13 = self.k*(v1x*v1z.conjugate() + v1x.conjugate()*v1z).real
                a21 = a12
                a22 = nu+self.k*(v1y*v1y.conjugate()+v1y.conjugate()*v1y).real
                a23 = self.k*(v1y*v1z.conjugate() + v1y.conjugate()*v1z).real
                a31 = a13
                a32 = a23
                a33 = nu+self.k*(v1z*v1z.conjugate()+v1z.conjugate()*v1z).real
                InverseMatrix = 1/(a11*a22*a33+a12*a23*a31+a13*a21*a32-a13*a22*a31-a12*a21*a33-a11*a23*a32) * np.array([[a22*a33-a23*a32,-(a12*a33-a13*a32),a12*a23-a13*a22],[-(a21*a33-a23*a31), a11*a33-a13*a31, -(a11*a23-a13*a21)],[a21*a32-a22*a31, -(a11*a32-a12*a31), (a11*a22-a21*a12)]])
            
                q_ = 1/2*np.dot(InverseMatrix, -v0dif/self.omega + mu * v0)
                # Γ(Φ) = ∫Z0(Φ+ωt)・q(ωt)dt
                Gamma += np.dot(q_.T, v0phi)
            # save
            GammaPhi[0,tp] = Pnum/self.Tnum*2*pi
            GammaPhi[1,tp] = self.Delta + Gamma/self.Tnum
        # plot
        plt.plot(GammaPhi[0,:], GammaPhi[1,:])
        plt.xlabel("Φ")
        plt.ylabel("Δ+Γ(Φ)")
        plt.xlim(-pi,pi)
        plt.grid()
        return GammaPhi

class tangent:
    def __init__(self, Tnum, u0_, v0_, v0_dif, omega, P, Delta):
        self.Tnum = Tnum
        self.u0_ = u0_
        self.v0_ = v0_
        self.v0_dif = v0_dif
        self.omega = omega
        self.P = P
        self.Delta = Delta
    
    # Z = v0_  I1 = v1_
    # Z'(θ) = 1/ω dZ/d
    def Calc_mu_nu(self):
        Zdif_u0_divide_u0_u0 = 0 # [<Z',u0>/|u0|^2]_t
        Zdif_u0_2_divide_u0_u0 = 0 # [<Z',u0>^2/|u0|^2]_t
        divide_u0_u0 = 0 # [1/|u0|^2]_t
        for tt in range(self.Tnum):
            u0 = self.u0_[:,tt:tt+1]
            v0dif = self.v0_dif[:,tt:tt+1]
            Zdif_u0_divide_u0_u0 += np.dot(v0dif.T, u0)/np.dot(u0.T, u0)
            Zdif_u0_2_divide_u0_u0 += np.dot(v0dif.T, u0)**2/np.dot(u0.T, u0)
            divide_u0_u0 += 1/np.dot(u0.T, u0)
        Zdif_u0_divide_u0_u0 = Zdif_u0_divide_u0_u0/self.Tnum/self.omega
        Zdif_u0_2_divide_u0_u0 = Zdif_u0_2_divide_u0_u0/self.Tnum/self.omega/self.omega
        divide_u0_u0 = divide_u0_u0/self.Tnum

        #lambdaを求める
        nu = 1/2*np.sqrt( (-Zdif_u0_divide_u0_u0**2 + Zdif_u0_2_divide_u0_u0*divide_u0_u0) / (self.P*divide_u0_u0 - self.Delta**2))
        #muを求める
        mu = (Zdif_u0_divide_u0_u0-2*nu*self.Delta)/(divide_u0_u0)
        return mu.item(), nu.item()
    
    def Calc_Gamma(self, mu, nu):
        # Γ(Φ) = ∫Z(Φ+ωt)・q(t)dt
        Division = 100 # Divide [-π, π] into 100 pieces.
        GammaPhi = np.empty((self.v0_.shape[0],Division+1))
        for tp in range(Division+1):
            Pnum = tp * int(self.Tnum/Division) - int(self.Tnum/2) # Φ(num) [-Tnum/2, Tnum/2]
            if(Pnum < 0):
                te = Pnum + self.Tnum # Φ = Pnum and Φ = Pnum + Tnum is equivalent for Z(Φ+ωt). 
            else:
                te = Pnum
            # Calculate Γ(Φ)
            Gamma = 0
            for tt in range(self.Tnum):
                # Z0(Φ+ωt)
                v0phi = np.array([self.v0_[:,(tt+te)%self.Tnum]]).T
                # q(ωt)
                u0 = self.u0_[:,tt:tt+1]
                v0dif = self.v0_dif[:,tt:tt+1]
                q_ = (-1/self.omega * np.dot(v0dif.T, u0) + mu)/(2*nu*np.dot(u0.T, u0)) * u0
    
                # Γ(Φ) = ∫Z0(Φ+ωt)・q(ωt)dt
                Gamma += np.dot(q_.T, v0phi)
            # save
            GammaPhi[0,tp] = Pnum/self.Tnum*2*pi
            GammaPhi[1,tp] = self.Delta + Gamma/self.Tnum
        # plot
        plt.plot(GammaPhi[0,:], GammaPhi[1,:])
        plt.xlabel("Φ")
        plt.ylabel("Δ+Γ(Φ)")
        plt.xlim(-pi,pi)
        plt.grid()
        return GammaPhi
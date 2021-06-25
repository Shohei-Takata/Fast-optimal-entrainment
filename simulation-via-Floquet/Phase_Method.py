import numpy as np
import sys

pi = np.pi

#############################################
#############################################
def Trans_PI(theta_num, Tnum): # Transform phase into [-π, π]
    theta = (theta_num % Tnum) / Tnum * 2*pi
    if(theta > pi):
        theta -= 2*pi
    return theta

def Trans_PI_SL(theta):
    theta = theta % (2*pi)
    if(theta > pi):
        theta -= 2*pi
    return theta

# We solve the equation, <x - ((1-beta)X0a + beta X0b), (1-beta)v0a + beta v0b> = 0.
def Linear_interpolation(x,v0a,v0b,X0a,X0b):
    # We solve the equation, a beta^2 + b beta + c = 0
    a = np.dot((X0b - X0a).T, (v0a - v0b)) 
    b = np.dot(x.T, (v0a - v0b)) + np.dot(v0b.T, (X0b - X0a)) - np.dot(X0b.T, (v0a - v0b))
    c = np.dot((x - X0b).T, v0b)
    if(b*b-4*a*c<0):
        beta = -b/2/a # For convenience, we define such. 
        # Note that phase is not accurate when the under error is displayed. 
        print("Phase Error (linear interpolation)")
    else:
        beta = ( - b - np.sqrt(b*b - 4*a*c) ) / 2 / a
    return beta

# Calculate phase via Floquet vector
# We define x = X(r, θ). 
# We use this equation, <x - X0(θ), v0(θ)> = 0. 
# θ' > θ → <x - X0(θ'), v0(θ')> is negative. 
# We find accurate phase by Linear interpolation. 

# This method is available when r << 1, x - X0(θ) << 1. 
# So, we can apply this method only to amplitude-feedback method. 
def Calc_phase_via_Floquet(x, X0_, v0_, X_phase, Tnum):
    # We assume that X0(θ) is in [X_phase-Interval, X_phase+Interval]. 
    # Here, X_phase is phase before 1 step. 
    Interval = int(Tnum/50) # search X0(θ) for "Interval" steps
    y = x - X0_[:,X_phase:X_phase+1]
    R = np.dot(y.T, v0_[:,X_phase:X_phase+1]) # R = <x - X0(θ'), v0(θ')>
    # First, we assume that X0(θ) is in [Phase_Temp-1, Phase_temp]. 
    if R >= 0: # The present phase goes ahead of 1 step before.
        for ts in range(Interval):
            Phase_Temp = (X_phase + ts) % Tnum
            y = x - X0_[:,Phase_Temp:Phase_Temp+1]
            R = np.dot(y.T, v0_[:,Phase_Temp:Phase_Temp+1])
            #Next, we find accurate phase by Linear interpolation. 
            if (R < 0): # θ in [Phase_Temp-1, Phase_Temp]
                beta = Linear_interpolation(x, v0_[:,Phase_Temp:Phase_Temp+1], np.array([v0_[:,(Phase_Temp-1)%Tnum]]).T, X0_[:,Phase_Temp:Phase_Temp+1], np.array([X0_[:,(Phase_Temp-1)%Tnum]]).T)
                #phase = Phase_Temp + beta
                x0 = (1-beta) * np.array([X0_[:,(Phase_Temp-1)%Tnum]]).T + beta * X0_[:,Phase_Temp:Phase_Temp+1]
                y = x - x0
                break
    else:# The present phase goes behind 1 step before.
        for ts in range(Interval):
            Phase_Temp = (X_phase - ts) % Tnum
            y = x - X0_[:,Phase_Temp:Phase_Temp+1]
            R = np.dot(y.T, v0_[:,Phase_Temp:Phase_Temp+1])
            #Next, we find accurate phase by Linear interpolation. 
            if (R >= 0): # θ in [Phase_Temp-1, Phase_Temp]
                beta = Linear_interpolation(x, v0_[:,Phase_Temp:Phase_Temp+1], np.array([v0_[:,(Phase_Temp-1)%Tnum]]).T, X0_[:,Phase_Temp:Phase_Temp+1], np.array([X0_[:,(Phase_Temp-1)%Tnum]]).T)
                #phase = Phase_Temp + beta
                x0 = (1-beta) * np.array([X0_[:,(Phase_Temp-1)%Tnum]]).T + beta * X0_[:,Phase_Temp:Phase_Temp+1]
                y = x - x0
                break
    if ts == Interval - 1:
        print("Phase Error (Interval over)")
        sys.exit()
    return y, Phase_Temp

# The phase is measured after the state is relaxed on the limit cycle.
def Calc_phase_directory(x, F, X0_, Tnum, dt, rotations = 3):
    X_Temp = np.copy(x)
    # we evol x for rotations times
    for tt in range(rotations*Tnum):
        #k1 = dt*F(X_Temp)
        #k2 = dt*F(X_Temp+0.5*k1)
        #k3 = dt*F(X_Temp+0.5*k2)
        #k4 = dt*F(X_Temp+k3)
        #X += (k1+2*k2+2*k3+k4)/6
        X_Temp += dt*F(X_Temp)
    
    # Find recent contacts for X_Temp
    Distance = 10000
    for tt in range(Tnum):
        Y = X_Temp - X0_[:,tt:tt+1]
        Distance_Temp = np.linalg.norm(Y)
        if Distance > Distance_Temp:
            X_phase = tt
            Distance = Distance_Temp
    return X_phase
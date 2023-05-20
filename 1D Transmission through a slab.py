# This file is the main program of Problem 1. Run this file to obtain 'Ez_vac_omega.txt' and 'Ez_slab_omega.txt' to plot transmission spectrum T in file 'Plot numerical transmission spectrum'.

import matplotlib.pyplot as plt
import numpy as np

# Initial parameters
mu0 = 1
c0 = 1
lambda0 = 8e-7
lambdaL = 6e-7
lambdaU = 1e-6
omega0 = 2*np.pi*c0/lambda0
sigma = 2*lambda0/omega0/(lambdaU-lambdaL)

dx = 1e-8
nx = 130

nt = int(1.5e4)
dt = dx

# Define sine wave envelope Gaussian pulse
def Source(time):
    tn = time*dt
    fsrc = np.exp(-(tn-4*sigma)**2/sigma**2)*np.sin(omega0*tn)
    return fsrc


ez = np.zeros(nx)
hy = np.zeros(nx)

x = np.arange(nx)

# Seting epsilon_0 and epsilon_r
epsilon = np.ones(nx)
# epsilon[30:100] = 10

# Initiate plotting motion graphics 
# fig = plt.figure()
# plt.ion()

# Initiate recording Ez at 120 space step
ez120 = []

for i in range(nt):
    
    # 1D FDTD loop
    hy[0:nx-1] = hy[0:nx-1] + dt*(ez[1:nx]-ez[0:nx-1])/(mu0*dx)
    # TFSF boundary for Hy
    hy[9] = hy[9] - dt/mu0/dx*Source(i)
 
    # Absorber condtion at the left
    ez[0] = ez[1]
    
    ez[1:nx] = ez[1:nx] + dt*(hy[1:nx]-hy[0:nx-1])/(epsilon[1:nx]*dx)

    # TFSF boundary for Ez
    ez[10] = ez[10] + dt/epsilon[10]/dx*Source(i+1)

    # Record Ez(t)
    ez120.append(ez[120])

    # Plot motion graphics
#     if i%10 == 0:
#         plt.clf()
        
#         ax1 = fig.add_subplot(211)
#         plt.title('ez time_step:'+str(i))
#         ax1.set_ylim([-1, 1])
#         ax1.plot(x,ez)

#         ax2 = fig.add_subplot(212)
#         plt.title('hy time_step:'+str(i))
#         ax2.set_ylim([-1, 1])
#         ax2.plot(x,hy)
#         plt.tight_layout() 

#         plt.pause(0.01)

# plt.ioff()
# plt.show()

# Plot Ez_slab(t) or Ez_vac(t)
plt.figure()
ez120 = np.array(ez120)
t = np.linspace(0,nt,nt)
plt.plot(t,ez120)
plt.xlabel('timestep')
# plt.ylabel('Ez_vac(t)')
plt.ylabel('Ez_slab(t)')

# plt.title('timestep vs Ez_vac(t)')
plt.title('timestep vs Ez_slab(t)')
plt.show()



from scipy.fftpack import fft,ifft

# Take FFT to obtain Ez(omega)
ezpmega = fft(ez120)
p2 = np.abs(ezpmega)   
p1 = p2[:int(nt/2)]
p1 = 2*p1/nt
Fs = 1/dt
f = np.arange(int(nt/2))*Fs/nt;
f[0] = f[1]
p1 = 2*p1/nt

# Transform frequency f into wavelengths l
l = 1/f # wavelengths l
a = np.where((6e-7<l)&(l<10e-7)) # Selecting l between lambdaU and lambdaL
l_re = l[a]

Ez_omega = p1[a] # Selecting Ez_omega corresponding to wavelength
Ez_omega = Ez_omega[::-1] # Reverse Ez_omega due to frequency switching to wavelength

# Save Ez(omega) data to plot numerical transmission spectrum T
# np.savetxt('Ez_slab_omega.txt',Ez_omega)
np.savetxt('Ez_vac_omega.txt',Ez_omega)

# Plot wavelengths l vs Ez_omega
plt.figure()
plt.plot(l_re,Ez_omega) 
plt.xlabel('wavelengths /m')
# plt.ylabel('Ez_vac(omega)')
plt.ylabel('Ez_slab(omega)')
# plt.title('wavelengths vs Ez_vac(omega)')
plt.title('wavelengths vs Ez_slab(omega)')
plt.show()





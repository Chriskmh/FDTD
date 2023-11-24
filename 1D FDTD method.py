import matplotlib.pyplot as plt
import numpy as np

# define three different current source function
def DC_Gaussian(t):
    t0 = 40
    tao = 12
    fsrc = np.exp(-0.5*np.power((t-t0)/tao,2))
    return fsrc

def CW_Sinusoid(t,dt,f0):
    fsrc = np.sin(2*np.pi*f0*t*dt)
    return fsrc

def Sine_Envelope_Gaussian(t,dt,f0):

# define Sine_Envelope_Gaussian by W_Sinusoid*DC_Gaussian
    fsrc = CW_Sinusoid(t,dt,f0)*DC_Gaussian(t) 
    return fsrc

# define constant 
epsilon0 = 8.85418782e-12
mu0 = 1.25663706e-6
c0 = 1/np.power(epsilon0*mu0,1/2)
f0 = 1e9
lambda0 = c0/f0
dx = lambda0/10
nx = 500

nt = 1000
dt = dx/c0

# initialize E field , H field and source 
ez = np.zeros(nx+1)
hy = np.zeros(nx)
jz = np.zeros(nx+1)

xe=np.linspace(0,(nx+1)*dx,nx+1)
xh=np.linspace(0,nx*dx,nx)

# initialize motion plot

# plt.figure()
# plt.ion()

ez_300 = []

for i in range(nt):
    
    # select source function
    
    # jz[250] = DC_Gaussian(i)
    # jz[250] = CW_Sinusoid(i,dt,f0)
    jz[250] = Sine_Envelope_Gaussian(i,dt,f0)

# Updating E and H
    ez[1:nx] = ez[1:nx] + dt*(hy[1:nx]-hy[0:nx-1])/(epsilon0*dx) - dt*jz[1:nx]
    hy[0:nx] = hy[0:nx] + dt*(ez[1:nx+1]-ez[0:nx])/(mu0*dx)

# PEC boundary Condition
    # ez[0] = 0
    # ez[(nx+1)-1] = 0

# Absorbing boundary Condition
    ez[0] = ez[1]
    ez[(nx+1)-1] = ez[(nx+1)-2]
   
# Recording E-field in i = 300 for item 2.3
    ez_300.append(ez[300])

# Ploting motion plot
    # if i%20 == 0:
    #     plt.clf()
        
    #     plt.subplot(2,1,1)
    #     plt.title('ez time_step:'+str(i))
    #     plt.ylim([-6e-11, 6e-11])
    #     plt.plot(xe,ez)

    #     plt.subplot(2,1,2)
    #     plt.title('hy time_step:'+str(i))
    #     plt.ylim([-2e-13, 2e-13])
    #     plt.plot(xh,hy)
    #     plt.tight_layout() 

    #     plt.pause(0.01)

# Plot for time steps: t = 100, 275, 500, 700
    if i == 100:
        plt.figure()
        plt.subplot(2,2,1)
        plt.title('Ez t=100')
        plt.ylim([-6e-11, 6e-11])
        plt.plot(xe,ez)

        plt.subplot(2,2,3)
        plt.title('Hy t=100')
        plt.ylim([-2e-13, 2e-13])
        plt.plot(xh,hy)


    if i ==275:
        plt.subplot(2,2,2)
        plt.title('Ez t=275')
        plt.ylim([-6e-11, 6e-11])
        plt.plot(xe,ez)

        plt.subplot(2,2,4)
        plt.title('Hy t=275')
        plt.ylim([-2e-13, 2e-13])
        plt.plot(xh,hy)
        plt.tight_layout() 

    if i ==500:
        plt.figure()
        plt.subplot(2,2,1)
        plt.title('Ez t=500')
        plt.ylim([-6e-11, 6e-11])
        plt.plot(xe,ez)

        plt.subplot(2,2,3)
        plt.title('Hy t=500')
        plt.ylim([-2e-13, 2e-13])
        plt.plot(xh,hy)

    if i ==700:
        plt.subplot(2,2,2)
        plt.title('Ez t=700')
        plt.ylim([-6e-11, 6e-11])
        plt.plot(xe,ez)

        plt.subplot(2,2,4)
        plt.title('Hy t=700')
        plt.ylim([-2e-13, 2e-13])
        plt.plot(xh,hy)
        plt.tight_layout() 

# Perform FFT for itwm 2.3
from scipy.fftpack import fft,ifft

ez_300 = np.array(ez_300)
freq = fft(ez_300)
bisidefreq = np.abs(freq)  
L = len(ez)
sisidefreq = bisidefreq[:int(L/2)]
f = np.arange(int(L/2))*f0/L;

plt.figure()
plt.plot(f,2*sisidefreq/L) 
plt.title('i = 300 in  frequency domain')
plt.ioff()
plt.show()

print("GPG_key")



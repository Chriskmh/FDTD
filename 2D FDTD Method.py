import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

# Initial parameters
mu0 = 1
c0 = 1

# Two group of wavelengths 
lambda0 = 8e-7
lambdaL = 4e-7
lambdaU = 12e-7

# lambda0 = 7e-7
# lambdaL = 7.15e-7
# lambdaU = 7.3e-7

# Define sine wave envelope Gaussian pulse
omega0 = 2*np.pi*c0/lambda0
sigma = 2*lambda0/omega0/(lambdaU-lambdaL)

def Source(time):
    tn = time*dt
    fsrc = np.exp(-(tn-4*sigma)**2/sigma**2)*np.sin(omega0*tn)
    return fsrc

# Design space step
dx = lambda0/80 # dx = 10nm
SizeX = 80 # SizeX = 80
dy = lambda0/80 # dy = 10nm
SizeY = 80 # SizeY = 80

# Design time step
MaxTime = int(1.5e4)
dt = dx/(2**0.5)/c0 

# Define convinient parameters
Cx = np.ones((SizeX,SizeY))*dt/dx
Cy = np.ones((SizeX,SizeY))*dt/dy

# Initiate Ez，Hx，Hy
Hx = np.zeros((SizeX,SizeY ))
Hy = np.zeros((SizeX ,SizeY))
Ez = np.zeros((SizeX,SizeY))

# Initiate recording Ez at (0,200)
Ez_record = []

# Initiate plotting coordinate
x = np.arange(SizeX)
y = np.arange(SizeY)
Y,X = np.meshgrid(y,x)


# Initiate plotting motion graphics 
# plt.figure()
# plt.ion

for time in range(MaxTime):

    # 2D FDTD loop
    Hx[0:SizeX,0:SizeY-1] = Hx[0:SizeX,0:SizeY-1] - Cy[0:SizeX,0:SizeY-1]*(Ez[0:SizeX,0+1:SizeY-1+1] - Ez[0:SizeX,0:SizeY-1])

    Hy[0:SizeX-1,0:SizeY] = Hy[0:SizeX-1,0:SizeY] + Cx[0:SizeX-1,0:SizeY]*(Ez[0+1:SizeX-1+1,0:SizeY] - Ez[0:SizeX-1,0:SizeY])

    Ez[1:SizeX-1,1:SizeY-1] = Ez[1:SizeX-1,1:SizeY-1] + Cx[1:SizeX-1,1:SizeY-1]*((Hy[1:SizeX-1,1:SizeY-1] - Hy[0:SizeX-2,1:SizeY-1]) - (Hx[1:SizeX-1,1:SizeY-1] - Hx[1:SizeX-1,0:SizeY-2]))

    # Generate source 
    Ez[int(SizeX/2),int(SizeY/4*3)] = Ez[int(SizeX/2),int(SizeY/4*3)] + dt*Source(time)

    # PEC Boundary conditions
    Ez[0:0,0:SizeY] = 0
    Ez[SizeX:SizeX,0:SizeY] = 0

    Ez[0:SizeX,0:0] = 0
    Ez[0:SizeX,SizeY:SizeY] = 0
    
    Hx[0:SizeX,0:0] = 0
    Hx[0:SizeX,SizeY-1:SizeY-1] = 0

    Hy[0:0,0:SizeY] = 0
    Hy[SizeX-1:SizeX-1,0:SizeY] = 0

    # Recording Ez at (0,200)
    Ez_record.append(Ez[int(SizeX/2),int(SizeY/4)])


    # Plot motion graphics
    # if time%50 == 0:
    #     print(time)
    #     plt.clf()
    #     plt.contour(X,Y,Ez)
    #     plt.colorbar()
    #     plt.title(str(time))
    #     plt.pause(0.1)

    # Record Ez,Hx,Hy at the end of the simulation
    if time == MaxTime-1:
        Ez_end = Ez
        Hy_end = Hy
        Hx_end = Hx

# Plotting motion graphics showing
# plt.ioff()
# plt.show()

Ez_record = np.array(Ez_record)

# Plot Ez(t)
def Plot_Ez_t():

    # Ez_record = np.array(Ez_record)
    t = np.linspace(0,MaxTime,MaxTime)

    plt.figure()
    plt.plot(t,Ez_record)
    plt.title('Ez(t) vs timestep')
    plt.xlabel('timestep')
    plt.ylabel('Ez(t)')
    plt.show()

# Plot Ez(omega) vs wavelengths
def Plot_Ez_omega():
    from scipy.fftpack import fft,ifft

    # Do fft
    Ezpmega = fft(Ez_record)
    p2 = np.abs(Ezpmega)   
    p1 = p2[:int(MaxTime/2)]
    p1 = 2*p1/MaxTime
    Fs = 1/dt
    f = np.arange(int(MaxTime/2))*Fs/MaxTime # frequency f
    f[0] = f[1]
    p1 = 2*p1/MaxTime

    # Transform frequency f into wavelengths l
    l = 1/f  # wavelengths l

    a = np.where((4e-7<l)&(l<12e-7)) # Selecting l between lambdaU and lambdaL
    l_re = l[a]

    Ez_omega = p1[a] # Selecting Ez_omega corresponding to wavelength
    Ez_omega = Ez_omega[::-1] # Reverse Ez_omega due to frequency switching to wavelength
     
    # Plot wavelengths l vs Ez_omega
    plt.figure()
    plt.plot(l_re,Ez_omega)

    # Select wave modes 
    b = np.where(Ez_omega>0.75e-14)
    print(l_re[b])
    
    # Mark wave modes
    plt.scatter(l_re[b],Ez_omega[b],c = 'r')
    plt.annotate('(1,1)',xy = (l_re[b][0],Ez_omega[b][0]))
    plt.annotate('(0,2)',xy = (l_re[b][1],Ez_omega[b][1]))
    plt.annotate('(1,3)',xy = (l_re[b][2],Ez_omega[b][2]))
    plt.annotate('(2,3)',xy = (l_re[b][3],Ez_omega[b][3]))

    plt.title('wavelengths vs Ez(omega)')
    plt.xlabel('wavelengths /m')
    plt.ylabel('Ez(omega)')

    plt.show()

# Plot Ez,Hx,Hy at the end of the simulation
def Plot_2d():

    # Plot Ez
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111,projection = '3d')
    ax1.plot_surface(X,Y,Ez_end,cmap = 'rainbow')
    ax1.set_xlabel('X /10nm')
    ax1.set_ylabel('Y /10nm')
    ax1.set_zlabel('Ez')
    plt.title('Ez at the end')

    # Plot Hy
    fig2 = plt.figure()
    ax2 = fig2.add_subplot(111,projection = '3d')
    ax2.plot_surface(X,Y,Hy_end,cmap = 'rainbow')
    ax2.set_xlabel('X /10nm')
    ax2.set_ylabel('Y /10nm')
    ax2.set_zlabel('Hy')
    plt.title('Hy at the end')

    # Plot Hx
    fig3 = plt.figure()
    ax3 = fig3.add_subplot(111,projection = '3d')
    ax3.plot_surface(X,Y,Hx_end,cmap = 'rainbow')
    ax3.set_xlabel('X /10nm')
    ax3.set_ylabel('Y /10nm')
    ax3.set_zlabel('Hx')
    plt.title('Hx at the end')

    plt.show()


# Plot_Ez_t()
Plot_Ez_omega()
# Plot_2d()

print('done')
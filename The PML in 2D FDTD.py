import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D


# Initial parameters
mu0 = 1
c0 = 1
episilon0 = 1

lambda0 = 10e-7
lambdaL = 5e-7
lambdaU = 15e-7

# Define sine wave envelope Gaussian pulse
omega0 = 2*np.pi*c0/lambda0
sigma = 2*lambda0/omega0/(lambdaU-lambdaL)


Lx = 3e-6
Ly = 3e-6

# Design space step
dx = 2e-8 # dx = 20nm
SizeX = int(Lx/dx) # SizeX = 150
dy = 2e-8 # dy = 20nm
SizeY = int(Ly/dy) # SizeY = 150

# Design time step
MaxTime = int(5e3)
dt = dx/(2**0.5)/c0 

def Source(time):
    tn = time*dt
    fsrc = np.exp(-(tn-4*sigma)**2/sigma**2)*np.sin(omega0*(tn-4*sigma))
    return fsrc

# Initiate Ez，Hx，Hy, Dz, Bx, By
# Initiate  Bx_pre, By_pre, Dz_pre to rememeber Dz, Bx, By from last timestep
Hx = np.zeros((SizeX,SizeY ))
Hy = np.zeros((SizeX ,SizeY))
Bx = np.zeros((SizeX,SizeY ))
By = np.zeros((SizeX ,SizeY))
Ez = np.zeros((SizeX,SizeY))
Dz = np.zeros((SizeX,SizeY))
Bx_pre = np.zeros((SizeX,SizeY))
By_pre = np.zeros((SizeX,SizeY))
Dz_pre = np.zeros((SizeX,SizeY))

# Initiate plotting coordinate
x = np.arange(SizeX)
y = np.arange(SizeY)
Y,X = np.meshgrid(y,x)

# Initiate PML
Lpml = 15 # PMl length
R_err = 1e-6 # R_err is the desired reflection
eta = (mu0*episilon0)**0.5

# Define Sigmap_p
def Sigmap(xnum):
    Sigma_max = -4*np.log(R_err)/(2*eta*Lpml*dx) 
    y = Sigma_max*(xnum/Lpml)**3
    return y

# Apply Sigmap_p to get sigma_x
SigmaX = np.zeros((SizeX,SizeY))
for i in range(Lpml):
    ipml = Lpml+1-i
    SigmaX[i,0:SizeY] = Sigmap(ipml)
    SigmaX[SizeX-i-1,0:SizeY] = SigmaX[i,0:SizeY]

# Apply Sigmap_p to get sigma_y
SigmaY = np.zeros((SizeX,SizeY))
for i in range(Lpml):
    ipml = Lpml+1-i
    SigmaY[0:SizeX,i] = Sigmap(ipml)
    SigmaY[0:SizeX,SizeY-i-1] = SigmaY[0:SizeX,i]

# Define convenient parameters from sigma_x and sigma_y
Cy_neg = 1-(dt*SigmaY/2)
Cy_pos = 1+(dt*SigmaY/2)
Cx_neg = 1-(dt*SigmaX/2)
Cx_pos = 1+(dt*SigmaX/2)

# Initiate recording Ez
Ez_record = []
Ez_230 = np.zeros((SizeX,SizeY))

# Initiate plotting motion graphics 
# plt.figure()
# plt.ion

for time in range(MaxTime):
    # 2D FDTD loop with pml
    Bx[0:SizeX,0:SizeY-1] = Cy_neg[0:SizeX,0:SizeY-1]/Cy_pos[0:SizeX,0:SizeY-1]*Bx[0:SizeX,0:SizeY-1] - dt/dy/mu0/Cy_pos[0:SizeX,0:SizeY-1]*(Ez[0:SizeX,0+1:SizeY-1+1] - Ez[0:SizeX,0:SizeY-1])

    By[0:SizeX-1,0:SizeY] = Cx_neg[0:SizeX-1,0:SizeY]/Cx_pos[0:SizeX-1,0:SizeY]*By[0:SizeX-1,0:SizeY] + dt/dx/mu0/Cx_pos[0:SizeX-1,0:SizeY]*(Ez[0+1:SizeX-1+1,0:SizeY] - Ez[0:SizeX-1,0:SizeY])

    Hx[0:SizeX,0:SizeY-1] = Hx[0:SizeX,0:SizeY-1] + Cx_pos[0:SizeX,0:SizeY-1]*Bx[0:SizeX,0:SizeY-1] - Cx_neg[0:SizeX,0:SizeY-1]*Bx_pre[0:SizeX,0:SizeY-1]

    Hy[0:SizeX-1,0:SizeY] = Hy[0:SizeX-1,0:SizeY] + Cy_pos[0:SizeX-1,0:SizeY]*By[0:SizeX-1,0:SizeY] - Cy_neg[0:SizeX-1,0:SizeY]*By_pre[0:SizeX-1,0:SizeY]
    
    Dz[1:SizeX-1,1:SizeY-1] = Cx_neg[1:SizeX-1,1:SizeY-1]/Cx_pos[1:SizeX-1,1:SizeY-1]*Dz[1:SizeX-1,1:SizeY-1] + dt/episilon0/Cx_pos[1:SizeX-1,1:SizeY-1]*((Hy[1:SizeX-1,1:SizeY-1]-Hy[0:SizeX-2,1:SizeY-1])/dy-(Hx[1:SizeX-1,1:SizeY-1] - Hx[1:SizeX-1,0:SizeY-2])/dx)

    Ez[1:SizeX-1,1:SizeY-1] = Cy_neg[1:SizeX-1,1:SizeY-1]/Cy_pos[1:SizeX-1,1:SizeY-1]*Ez[1:SizeX-1,1:SizeY-1] + 1/Cy_pos[1:SizeX-1,1:SizeY-1]*(Dz[1:SizeX-1,1:SizeY-1]-Dz_pre[1:SizeX-1,1:SizeY-1])

    # Generate source 
    Ez[int(SizeX/2),int(SizeY/2)] = Ez[int(SizeX/2),int(SizeY/2)] + dt*Source(time)

    # Record from last timestep
    Bx_pre[0:SizeX,0:SizeY] = Bx[0:SizeX,0:SizeY]
    By_pre[0:SizeX,0:SizeY] = By[0:SizeX,0:SizeY]
    Dz_pre[0:SizeX,0:SizeY] = Dz[0:SizeX,0:SizeY]

    # Recording Ez at (125,75)
    Ez_record.append(Ez[125,75])

    # Record Ez at the time timestep = 230
    if time == 230:
        Ez_230[:,:] = Ez[:,:]

    # Plot motion graphics
    # if time%10 == 0:
    #     print(time)
    #     plt.clf()
    #     plt.contourf(X,Y,Ez)
    #     plt.colorbar()
    #     plt.title(str(time))
    #     plt.pause(0.1)


# Plotting motion graphics showing
# plt.ioff()
# plt.show()

Ez_record = np.abs(np.array(Ez_record))

# Plot Ez at (125,75) with t 
def Plot_Ez_t():

    # Ez_record = np.array(Ez_record)
    t = np.linspace(0,MaxTime,MaxTime)

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.plot(t,Ez_record)
    ax1.set_xscale("log")
    plt.title('|Ez(t)| vs timestep')
    plt.xlabel('timestep')
    plt.ylabel('|Ez(t)|')


# Plot Ez(omega) vs wavelengths
Plot_Ez_t()

# Plot Ez at timestep = 230
def plot_Ez_230():
    fig = plt.figure()
    ax2 = fig.add_subplot(111, projection='3d')
    ax2.plot_surface(X, Y, Ez_230,cmap='coolwarm')
    # ax2.set_zlim([0,1.75e-10])
    ax2.set_xlabel('X')
    ax2.set_ylabel('Y')
    ax2.set_zlabel('Ez')
    ax2.set_title('Ez at time timestep = 230')
    plt.show()

plot_Ez_230()
print('done')
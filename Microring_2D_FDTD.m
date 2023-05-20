%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2D Finite-Difference Time-Domain (FDTD)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% I refer to your program and modify my program. Fdtd and PML code are
% still my version. I usc 450 spacesteps and 5e4 timestep. This program runs for about 20 minutes on my laptop.
% You do not have to run the simulation again. I have recorded the data in
% 'Probe.mat'. Run 'Data_process.m' to see the data. 

clc, clear, close all;

%% Initial parameters
tic;

mu0 = 1;
c0 = 1;
episilon0 = 1;
episilon_waveguide = 12;
lambda = [1555:1:1555];
length_lambda = length(lambda);

Probe_input_cell = cell(1,1);
Probe_drop_cell = cell(1,1);
Probe_through_cell = cell(1,1);

for lambda_index = 1:1:length_lambda
lambda0 = lambda(lambda_index); % center wavelength in nanometers
disp(['lambda0 = ',num2str(lambda0)]);

lambdaL = 1550;  % left wavelength in nanometers
lambdaU = 1558; % right wavelength in nanometers

% Here I refer to your program and reduce the number of space step to SizeX = 450, SizeY = 450, and increase Lx and Ly a little bit to place PML well.
% I think in your program 2200 space step is too slow and I do not think we need such accuracy. 440 steps is enough.
dx = 50;
dy = 50;
Lx = 22500;
Ly = 22500;
SizeX = round(Lx/dx); % Space step SizeX = 450
SizeY = round(Ly/dy); % Space step SizeY = 450

% I am not sure that how long the program need to run, so I run it in 5e4 timestep. 
% In such long time it is enough to reduce the Ez at drop port to quite small. 
MaxTime = round(7000);
dt = dx/(2^0.5)/c0;

%% Define geometry of Microring
% Although the geometry program is slight different with yours, they are actually the same.

episilon = ones(SizeX,SizeY);

waveguide_width = 400/dx; % waveguide bus width in nanometers (400 nm)
waveguide_ring_width = 400/dx; % waveguide ring width in nanometers (400 nm)
waveguide_outer_radius = 10000/dx; % outer radius in nanometers (10 um)
waveguide_inner_radius = waveguide_outer_radius-waveguide_ring_width; % inner radius in nanometers
waveguide_gap = 100/dx; % gap in nanometers

waveguide_ring_location = SizeX/2;

% Define ring waveguide
for i = 1:1:SizeY
    for j = 1:1:SizeX
        if (i-waveguide_ring_location)^2+(j-waveguide_ring_location)^2 < waveguide_outer_radius^2
            episilon(i,j) = episilon_waveguide;
        end
    end
end

for i = 1:1:SizeY
    for j = 1:1:SizeX
        if (i-waveguide_ring_location)^2+(j-waveguide_ring_location)^2 < waveguide_inner_radius^2
            episilon(i,j) = 1;
        end
    end
end

% Define top waveguide
top_waveguide_location = waveguide_ring_location + waveguide_gap + waveguide_outer_radius;
episilon(top_waveguide_location:top_waveguide_location + waveguide_width,:) = episilon_waveguide;

% Define bottom waveguide
bottom_waveguide_location = waveguide_ring_location - waveguide_gap - waveguide_outer_radius;
episilon(bottom_waveguide_location - waveguide_width:bottom_waveguide_location,:) = episilon_waveguide;

% Show the waveguide geometry
% figure(1);
% contourf(episilon);
% colorbar;
% titlestring = ('Waveguide Simulation Space Material (Relative Permittivity)');
% title(titlestring,'color','k');
% xlabel('x in space step');
% ylabel('y in space step');
% axis equal;
%% Define Source sine wave envelope Gaussian pulse
omega0 = 2*pi*c0/lambda0;
sigma = 2*lambda0/omega0/(lambdaU-lambdaL);
% Source = @(time)exp(-(time*dt-4*sigma)^2/sigma^2)*sin(omega0*(time*dt-4*sigma));
Source = @(time)sin(omega0*(time*dt-4*sigma));

% Define source location
Source_top_y = top_waveguide_location + waveguide_width - 1;
Source_bottom_y = top_waveguide_location + 1;
Source_x = 15;

%% Initiate field vectors
Hx = zeros(SizeX,SizeY);
Hy = zeros(SizeX,SizeY);
Bx = zeros(SizeX,SizeY);
By = zeros(SizeX,SizeY);
Ez = zeros(SizeX,SizeY);
Dz = zeros(SizeX,SizeY);
Bx_pre = zeros(SizeX,SizeY);
By_pre = zeros(SizeX,SizeY);
Dz_pre = zeros(SizeX,SizeY);

% Ploting array
x = [1:1:SizeX];
y = [1:1:SizeY];
% [X,Y] = meshgrid(x,y);

%% Initiate PML
Lpml = 10; % PMl length
R_err = 1e-6; % R_err is the desired reflection
eta = (mu0/episilon0)^0.5;
Sigma_max = -4*log(R_err)/(2*eta*Lpml*dx);

Sigmap = @(xnum)Sigma_max*(xnum/Lpml)^3;

% PML on y direction
SigmaY = zeros(SizeX,SizeY);
for i = 1:1:Lpml
    ipml = Lpml+1-i;
    SigmaY(1:SizeX,i) = Sigmap(ipml);
    SigmaY(1:SizeX,SizeY-i+1) = SigmaY(1:SizeX,i);
end

% PML on x direction
SigmaX = zeros(SizeX,SizeY);
for i = 1:1:Lpml
    ipml = Lpml+1-i;
    SigmaX(i,1:SizeY) = Sigmap(ipml);
    SigmaX(SizeX-i+1,1:SizeY) = SigmaX(i,1:SizeY);
end

%% Some other paratemers
% Define Update Equation Coefficients
Cy_neg = 1-(dt*SigmaY/2);
Cy_pos = 1+(dt*SigmaY/2);
Cx_neg = 1-(dt*SigmaX/2);
Cx_pos = 1+(dt*SigmaX/2);

%% Define probe 
% Probe location. Probes are located at input port, through port, and drop
% port.
Probe_input_bottom_y = Source_bottom_y;
Probe_input_top_y = Source_top_y;
Probe_input_x = Source_x + 3;

Probe_through_bottom_y = Source_bottom_y;
Probe_through_top_y = Source_top_y;
Probe_through_x = SizeX - Probe_input_x;

Probe_drop_bottom_y = SizeY - Source_top_y;
Probe_drop_top_y = SizeY - Source_bottom_y;
Probe_drop_x = Probe_input_x;

% Define probe array
length_probe = Probe_input_top_y  - Probe_input_bottom_y + 1;
Probe_drop = zeros(length_probe,MaxTime);
Probe_input = zeros(length_probe,MaxTime);
Probe_through = zeros(length_probe,MaxTime);

%% Plot the location of the waveguide, probe and source

plot_Source_probe = episilon;
plot_Source_probe(Source_bottom_y:Source_top_y,Source_x) = 5;
plot_Source_probe(Probe_input_bottom_y:Probe_input_top_y,Probe_input_x) = 8;
plot_Source_probe(Probe_through_bottom_y:Probe_through_top_y,Probe_through_x) = 8;
plot_Source_probe(Probe_drop_bottom_y:Probe_drop_top_y,Probe_drop_x) = 8;

figure(2);
contourf(plot_Source_probe);
title('Source and probe location');
text(Source_x,Source_bottom_y-3,'Source','Color','r');
text(Probe_input_x+3,(Probe_input_bottom_y+Probe_input_top_y)/2,'Probe input','Color','r');
text(Probe_through_x+3,(Probe_through_bottom_y+Probe_through_top_y)/2,'Probe through','Color','r');
text(Probe_drop_x+3,(Probe_drop_bottom_y+Probe_drop_top_y)/2,'Probe drop','Color','r');
xlabel('x');
ylabel('y');
axis equal;


%% 2D FDTD loop
% figure;
for time = 1:1:MaxTime
    Bx(1:SizeX,1:SizeY-1) = Cy_neg(1:SizeX,1:SizeY-1)./Cy_pos(1:SizeX,1:SizeY-1).*Bx(1:SizeX,1:SizeY-1) - dt/dy/mu0./Cy_pos(1:SizeX,1:SizeY-1).*(Ez(1:SizeX,2:SizeY) - Ez(1:SizeX,1:SizeY-1));
    By(1:SizeX-1,1:SizeY) = Cx_neg(1:SizeX-1,1:SizeY)./Cx_pos(1:SizeX-1,1:SizeY).*By(1:SizeX-1,1:SizeY) + dt/dx/mu0./Cx_pos(1:SizeX-1,1:SizeY).*(Ez(2:SizeX,1:SizeY) - Ez(1:SizeX-1,1:SizeY));
    Hx(1:SizeX,1:SizeY-1) = Hx(1:SizeX,1:SizeY-1) + Cx_pos(1:SizeX,1:SizeY-1).*Bx(1:SizeX,1:SizeY-1) - Cx_neg(1:SizeX,1:SizeY-1).*Bx_pre(1:SizeX,1:SizeY-1);
    Hy(1:SizeX-1,1:SizeY) = Hy(1:SizeX-1,1:SizeY) + Cy_pos(1:SizeX-1,1:SizeY).*By(1:SizeX-1,1:SizeY) - Cy_neg(1:SizeX-1,1:SizeY).*By_pre(1:SizeX-1,1:SizeY);
    Dz(2:SizeX-1,2:SizeY-1) = Cx_neg(2:SizeX-1,2:SizeY-1)./Cx_pos(2:SizeX-1,2:SizeY-1).*Dz(2:SizeX-1,2:SizeY-1) + dt./episilon(2:SizeX-1,2:SizeY-1)./Cx_pos(2:SizeX-1,2:SizeY-1).*((Hy(2:SizeX-1,2:SizeY-1) - Hy(1:SizeX-2,2:SizeY-1))./dy - (Hx(2:SizeX-1,2:SizeY-1) - Hx(2:SizeX-1,1:SizeY-2))./dx);
    Ez(2:SizeX-1,2:SizeY-1) = Cy_neg(2:SizeX-1,2:SizeY-1)./Cy_pos(2:SizeX-1,2:SizeY-1).*Ez(2:SizeX-1,2:SizeY-1) + 1./Cy_pos(2:SizeX-1,2:SizeY-1).*(Dz(2:SizeX-1,2:SizeY-1) - Dz_pre(2:SizeX-1,2:SizeY-1));
    
    % Source
    if time<500
    Ez(Source_bottom_y:Source_top_y,Source_x) = Ez(Source_bottom_y:Source_top_y,Source_x) + dt.*Source(time);
    end

    % Bx,By,Dz in last timestep
    Bx_pre(1:SizeX,1:SizeY) = Bx(1:SizeX,1:SizeY);
    By_pre(1:SizeX,1:SizeY) = By(1:SizeX,1:SizeY);
    Dz_pre(1:SizeX,1:SizeY) = Dz(1:SizeX,1:SizeY);

%     Probe Ez 
    Probe_input(1:length_probe,time) = Ez(Probe_input_bottom_y:Probe_input_top_y,Probe_input_x);
    Probe_through(1:length_probe,time) = Ez(Probe_through_bottom_y:Probe_through_top_y,Probe_through_x);
    Probe_drop(1:length_probe,time) = Ez(Probe_drop_bottom_y:Probe_drop_top_y,Probe_drop_x);

    if mod(time,1000) == 0
        disp(time);
    end

    % I the dynamic graph here the y axis of imagesc is rerverse, but the contourf is correct.
    % I could not figure out the reason.

%     dynamic graph 
%     if mod(time,100) == 0
%     imagesc(x,y',Ez,[0,100]);
% %     contourf(x,y,Ez);
% %     colorbar;
%     titlestr = ['2D FDTD for Ez',num2str(time)];
%     title(titlestr);
%     xlabel('x in space step');
%     ylabel('y in space step');
%     drawnow;
%     end
end

%% Data saving
% Data is save at 'Probe.mat', and processed in 'Data_process.m'
Probe_input_cell(1,lambda_index) = {Probe_input};
Probe_drop_cell(1,lambda_index) = {Probe_drop};
Probe_through_cell(1,lambda_index) = {Probe_through};
end
save('Probe2.mat','Probe_input_cell',"Probe_drop_cell","Probe_through_cell");

toc;
disp('done');
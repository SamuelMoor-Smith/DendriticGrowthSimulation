% ------------------------
% Define Parameters and Setup Simulation for Gold Memristor Simulation
% ------------------------
clear global

% Note: Wouldn't modify the !!! terms unless for good reason
% Commented out some of the parameters used for temperature evolution as I did not find this very important to model

% ------------------------
% Random Number Generator Settings
% ------------------------
% s = rng; % Get the current random number generator settings
% disp(s.Seed);
rng(2)

% ------------------------
% Basic Simulation Parameters
% ------------------------
params.i = 1; % Simulation Number

% ------------------------
% Well-understood Physical Constants
% ------------------------
% params.kB = 1.38e-11; % !!! Boltzmann constant [J/K]
params.T_coeff = 10000; % !!! Temperature coefficient [K^-1] -- this is a guess, not sure what to use here

% ------------------------
% Simulation Length and Particle Count
% ------------------------
params.n = 500; % Number of particles
n = params.n; % For convenience

simulation_length = 6; % in seconds
simulation_steps = 2000;
params.tspan=linspace(0,simulation_length,simulation_steps);

% ------------------------
% Initial position box and material size
% ------------------------
params.Lx = 100; % Size of material
params.Ly = 100; % Size of material

Lx = params.Lx; % For convenience
Ly = params.Ly; % For convenience

params.electrode_width = Lx/2; % Width for the gold deposits on either side
params.electrode_height = Ly - 0.2*Ly; % Height for gold deposits on either side

% ------------------------
% Heat related parameters
% ------------------------

params.To = 300; % !!! Initial Temperature
params.Q = 1; % Thermal energy [J]
% Jasons code?
params.k = 420/Lx; % !!! Heat transfer coefficient [W/m^2K] from https://www.spiraxsarco.com/learn-about-steam/steam-engineering-principles-and-heat-transfer/heat-transfer
params.CT = 0.24*5e-9; % !!! Heat capacity = specific heat capacity*mass [J/K] http://www2.ucdsb.on.ca/tiss/stretton/database/Specific_Heat_Capacity_Table.html -- check for gold

% ------------------------
% Pinning Force Parameters
% ------------------------
params.m = 250; % Number of sites
m = params.m; % For convenience

% ------------------------
% Pinning Site Locations and Forces
% ------------------------

params.wpa_repulse = 3500;
% params.wpa_attract = 500;
% wp_attract=(-(rand(1,0.5*m)).*params.wpa_attract - params.wpa_attract/2);
wp_repulse=((rand(1,m)).*params.wpa_repulse + 0.25*params.wpa_repulse);
% params.w_pin=rand(m,1)*5000-500;
params.w_pin= wp_repulse;% [wp_attract wp_repulse];
params.x_pin=rand(m,1).*0.8*Lx + 0.1*Lx;
% params.y_pin=(rand(m,1).*params.electrode_height)-0.5*params.electrode_height;
params.y_pin=rand(m,1).*0.8*Ly - 0.4*Ly;
% params.x_pin_raw = linspace(0, Lx, 16);
% params.y_pin_raw = linspace(-Ly/2, Ly/2, 16);

% params.x_pin = zeros(m,1);
% params.y_pin = zeros(m,1);
% for i=1:16
%     for j=1:16
%         params.x_pin((i-1)*16 + j) = params.x_pin_raw(i);
%         params.y_pin((i-1)*16 + j) = params.y_pin_raw(j);
%     end
% end


% params.w_pin = 100; % Pinning potential amplitude
params.R_pin = rand(m,1).*10+1; % 5; % Pinning potential distance

params.w_pin2 = 2000;

% ------------------------
% Applied Electric Field Parameters
% ------------------------
params.V = 1; % Applied Voltage
params.alpha = rand(n,1).*1e5+0.2e5;

% ------------------------
% Drag Force Parameters
% ------------------------
params.eta = 1; % Viscousity
params.Cd = 50; % 1.8e4; % Drag Coefficient

% ------------------------
% Interfacial Potential Parameters
% ------------------------
params.wI = 2e4; % Interfacial potential amplitude only with clusters
params.RI = 25; % Interfacial potential distance only with clusters


% ------------------------
% Residual Stress Parameters
% ------------------------
params.E=323.5e9; % Young's modulus of AlN in Pa = N/m^2
params.nu=0.23; % Poisson's ratio of AlN, not sure if this should be used or be constant
params.kC=8.99e9; % Coulomb's constant
params.q=1.602e-19; % Electric charge
params.r_Ag=126e-12; % Silver ion radius

% ------------------------
% Current Calculation Parameters
% ------------------------
params.num_e = 100; % Number of electrons to simulate
params.Rt = 1;                         
params.lambda = 1;
params.steps = 250;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%% RUN DENDRITIC SIMULATION %%%%%%%%%%%%%%%

addpath('Dendritic_Evolution');
runSimulation(params);
% macromodel_ode_2d_V2_updated(params);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

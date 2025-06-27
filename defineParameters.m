% ------------------------
% Define Parameters and Setup Simulation for Gold Memristor Simulation
% ------------------------

% Note: Wouldn't modify the !!! terms unless for good reason
% Commented out some of the parameters used for temperature evolution as I did not find this very important to model

% ------------------------
% Random Number Generator Settings
% ------------------------
% s = rng; % Get the current random number generator settings
% disp(s.Seed);
rng(42)

% ------------------------
% Basic Simulation Parameters
% ------------------------
params.i = 1; % Simulation Number
params.fin = 1; % Something to do with stopping the simulation at a point

% ------------------------
% Well-understood Physical Constants
% ------------------------
params.kB = 1.38e-23; % !!! Boltzmann constant [J/K]

% ------------------------
% Simulation Length and Particle Count
% ------------------------ß
params.n = 400; % Number of particles
n = params.n; % For convenience

simulation_length = 3e-3; % 10 millisecond
params.tspan=linspace(0,simulation_length,500);

% ------------------------
% Initial position box and material size
% ------------------------
params.L = 50e-6; % Size of material
L = params.L; % For convenience

params.xoWidth = 100e-6; % Width for initial position box
params.yoHeight = 10e-6; % Height for initial position box

% ------------------------
% Heat related parameters
% ------------------------

params.To = 300; % !!! Initial Temperature
params.Q = 1; % Thermal energy [J]
% Jasons code?
params.k = 420/L; % !!! Heat transfer coefficient [W/m^2K] from https://www.spiraxsarco.com/learn-about-steam/steam-engineering-principles-and-heat-transfer/heat-transfer
params.CT = 0.24*5e-9; % !!! Heat capacity = specific heat capacity*mass [J/K] http://www2.ucdsb.on.ca/tiss/stretton/database/Specific_Heat_Capacity_Table.html

% ------------------------
% Pinning Force Parameters
% ------------------------
params.m = 20; % Number of sites
m = params.m; % For convenience

params.wpa_repulse = 0;
params.wpa_attract = 300;
params.Rp = 8e-6; % Pinning potential distance

% ------------------------
% Pinning Site Locations and Forces
% ------------------------
wpplus=((rand(1,m)).*params.wpa_attract);
wpminus=(-(rand(1,0)).*params.wpa_repulse);
params.wp=[wpplus wpminus];
params.xp=rand(m,1).*L;
params.yp=(rand(m,1).*10e-6)-5e-6;

% ------------------------
% Applied Electric Field Parameters
% ------------------------
params.V = 1; % Applied Voltage
params.alpha = rand(n,1).*8e-3;

% ------------------------
% Drag Force Parameters
% ------------------------
params.eta = rand(n,1).*4 + 1; % Viscousity
params.Cd = 1.8e4; % Drag Coefficient

% ------------------------
% Interfacial Potential Parameters
% ------------------------
params.wI = 4e3; % Interfacial potential amplitude
params.RI = 7e-6; % Interfacial potential distance

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
params.num_e = 1000; % Number of electrons to simulate
params.Rt = 1;                         
params.lambda = 5e-6;
params.steps = 1000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%% RUN DENDRITIC SIMULATION %%%%%%%%%%%%%%%

addpath('Dendritic_Evolution');
runSimulation(params);
% macromodel_ode_2d_V2_updated(params);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

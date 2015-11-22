mrstModule add mimetic

clear;
close all;

%% Configuration
% Length of domain in meters
L = 400;

% Depth (vertical dimension) of domain in meters
D = 100;

% Number of cells in x, z and t dimensions
nx = 100;
nz = 100;
nt = 10;

% Size of time step in seconds
dt = 0.01;

% Gravitational acceleration in m/s^2
g = 9.81;

% The wave amplitude in meters (rough approximation)
A = 10;

% Mass of foreign object hitting water in kg
m = 1.404 * 10^11;

% Velocity of foreign object hitting water in m/s
v0 = 45;

% Kinetic energy of foreign object as it hits the water, in J
E_k = 0.5 * m * v0^2;

% Dimension of impact zone, in meters (should be small compared to L)
d = 5;

%% Boundary and initial conditions
%eta0 = @(x) 0.01 * sin(5 * pi * x);
eta0 = @(x) zeros(size(x));
h =   @(x) zeros(size(x));  % 0.05 * cos(3 * pi * x);
%h = @(x) (x > 0.5) .* (1.5 * x - 0.75);

X = L * linspace(0, 1, nx);
phi_amplitude = E_k / L .* (1 / erf(L / (sqrt(2) * d)));
phi_top = - phi_amplitude * exp(-(X - (L / 2)).^2/(d^2 / 2));
%phi_top = zeros(size(X));

%% Scales
% Scales for time and potential
T = L / sqrt(g * D);
P = A * L * sqrt(g / D);

% Set up scaled versions of variables
h_scaled = @(x) h(x) / D;
eta0_scaled = @(x) eta0(x) / A;
dt_scaled = dt / T;
phi_top_scaled = phi_top / P;

%% Setup simulation
sim = simulation(nx, nz, g, eta0_scaled, h_scaled);

%% Compute surface eta
eta = solve_wave(sim, phi_top_scaled, dt_scaled, nt, 'ShowProgress');

%% Inverse scale
% Recover variables in original scale
eta = A * eta;

x = L * linspace(0, 1, nx + 1);
t = dt * (0 : nt);

%% Visualize
eta_max = max(max(eta));
eta_min = min(min(eta));
eta_dist = abs(eta_max - eta_min);
figure
fprintf('Press enter to start plotting...\n');
while true
    pause;
    for i = 1:length(t)
        newplot;
        plot(x, eta(:, i));
        %         ylim([-D, A]);
        ylim([eta_min - eta_dist, eta_max + eta_dist]);
        pauselength = 50 / length(t);
        pause(min(pauselength, 0.2));
    end
    fprintf('Done plotting. Press enter to plot again or ctrl+c to abort. \n');
end
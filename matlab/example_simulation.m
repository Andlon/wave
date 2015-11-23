mrstModule add mimetic

clear;
close all;

%% Configuration
% Length of domain in meters
L = 1000;

% Depth (vertical dimension) of domain in meters
D = 100;

% Width of fjord, in meters
W = 500;

% Number of cells in x, z and t dimensions
nx = 100;
nz = 25;
nt = 500;

% Size of time step in seconds
dt = 0.01;

% Mass of foreign object hitting water in kg
m = 1.404 * 10^1;

% Velocity of foreign object hitting water in m/s
v0 = 45;

% Kinetic energy of foreign object as it hits the water, in J
E_k = 0.5 * m * v0^2;

% Dimension of impact zone, in meters (should be small compared to L)
d = 50;

%% Boundary and initial conditions
eta0 = @(x) zeros(size(x));
h =   @(x) zeros(size(x));

X = L * linspace(0, 1, nx);
phi_amplitude = 8 * E_k / (sqrt(2) * pi * d * erf(L / (sqrt(2) * d)));
phi_top = phi_amplitude * exp(-(X - (L / 2)).^2/(d^2 / 2));
% phi_top = - exp(-(X - (L / 2)).^2/(2 * d^2));
% phi_top = - 0.2 * ( (X > (L / 2)) .* (-X + L) + (X <= (L/2)) .* X ); %zeros(size(X));

%% Setup simulation
sim = simulation(L, D, nx, nz, eta0, h);
sim.g = 9.81;

%% Compute surface eta
eta = solve_wave(sim, phi_top, dt, nt, 'ShowProgress');

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
        ylim([eta_min - eta_dist, eta_max + eta_dist]);
        pause(dt);
    end
    fprintf('Done plotting. Press enter to plot again or ctrl+c to abort. \n');
end
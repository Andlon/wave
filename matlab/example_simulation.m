mrstModule add mimetic

clear;
close all;

%% Configuration
% Length of domain in meters
L = 12000;

% Depth (vertical dimension) of domain in meters
D = 400;

% Width of fjord, in meters
W = 500;

% Number of cells in x, z and t dimensions
nx = 75;
nz = 15;
nt = 200;

% Size of time step in seconds
dt = 1;

% Mass of foreign object hitting water in kg
m = 1 * 10^5;

% Velocity of foreign object hitting water in m/s
v0 = 45;

% Kinetic energy of foreign object as it hits the water, in J
E_k = 0.5 * m * v0^2;

% Fraction of kinetic energy transferred from the rock to the water
conservation_factor = 0.05;

% Dimension of impact zone, in meters (should be small compared to L)
d = 500;

%% Boundary and initial conditions
eta0 = @(x) zeros(size(x)); %exp(-(X - (L / 2)).^2/(2 * d^2)); zeros(size(x));  %-0.01 * exp(-(x - (L / 2)).^2/(2 * d^2));
h =   @(x) zeros(size(x));

X = L * linspace(0, 1, nx);
phi_amplitude = 8 * E_k * conservation_factor / (sqrt(2) * pi * d * W * erf(L / (sqrt(2) * d)));
phi_top = phi_amplitude * exp(-(X - (L / 2)).^2/(d^2 / 2));
% phi_top = 80 * exp(-(X - (L / 2)).^2/(2 * d^2));
% phi_top = - 0.2 * ( (X > (L / 2)) .* (-X + L) + (X <= (L/2)) .* X ); %zeros(size(X));

%% Setup simulation
sim = simulation(L, D, nx, nz, eta0, h);

%% Compute surface eta
[eta] = solve_wave(sim, phi_top, dt, nt, 'ShowProgress');

x = L * linspace(0, 1, nx + 1);
t = dt * (0 : nt);

%% Visualize
eta_max = max(max(eta));
eta_min = min(min(eta));
eta_dist = abs(eta_max - eta_min);

figure;
fprintf('Press enter to start plotting...\n');
while true
    pause;
    for i = 1:10:length(t)
        area(x, eta(:, i), -D);
        ylim([- eta_dist, eta_dist]);
        title(sprintf('t = %f', i * dt));
        
%         pause(min(dt, 0.1));
        pause
    end
    fprintf('Done plotting. Press enter to plot again or ctrl+c to abort. \n');
end
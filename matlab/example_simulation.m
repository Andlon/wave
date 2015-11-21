mrstModule add mimetic

clear;
close all;

%% Configuration
nx = 50;
nz = 50;
nt = 100;
dt = 0.005;
g = 9.81;

%% Boundary and initial conditions
eta0 = @(x) zeros(size(x)); %  0.1 * sin(5 * pi * x);
h =   @(x) zeros(size(x));  % 0.05 * cos(3 * pi * x);

X = linspace(0, 1, nx);
phi_top = exp(-(X-0.5).^2) - 0.8;
%phi_top = X;

%% Setup simulation
sim = simulation(nx, nz, g, eta0, h);

%% Compute surface eta
eta = solve_wave(sim, phi_top, dt, nt);

x = linspace(0, 1, nx + 1);
t = dt * (0 : nt);

%% Visualize
figure
for i = 1:length(t)
    newplot;
    plot(x, eta(:, i));
    ylim([-0.005, 0.005])
    pause(0.1)
end

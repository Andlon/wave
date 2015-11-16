mrstModule add mimetic

%% Configuration
nx = 50;
nz = 25;
nt = 500;
dt = 0.01;

dx = 1 / nx;

%% Boundary and initial conditions
eta_initial = @(x) 0.1 * sin(5 * pi*x);
h =   @(x) 0.05 * cos(3 * pi * x);

X = linspace(0, 1, nx);
%phi_top = 25 * exp(-25*(X-0.5).^2);
phi_top = X;

%% Grid construction
[ G, top, bottom, left, right ] = setup_grid(eta_initial, h, nx, nz);

%% Compute surface eta
eta = solve_wave(G, phi_top, top, left, right, h, nx, nz, dt, nt, 9.81 / 10);

x = linspace(0, 1, nx + 1);
t = dt * (0 : nt);
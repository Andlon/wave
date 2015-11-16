mrstModule add mimetic

%% Configuration
nx = 50;
nz = 50;
dt = 0.01;

dx = 1 / nx;

%% Boundary and initial conditions
eta_initial = @(x) 0.1 * sin(2 * pi*x);
h =   @(x) 0.05 * cos(3 * pi * x);

X = linspace(0, 1, nx);
% phi_top = 25 * exp(-25*(X-0.5).^2);
phi_top = X;

%% Grid construction
[ G, top, bottom, left, right ] = setup_grid(eta_initial, h, nx, nz);

% %% (Temporary) Velocity field plot
% gradients = cell_gradients(G, v, 1:G.cells.num);
% 
% x = G.cells.centroids(:, 1);
% y = G.cells.centroids(:, 2);
% ux = gradients(:, 1);
% uy = gradients(:, 2);
% quiver(x, y, ux, uy);

%% Compute surface eta
eta = solve_wave(G, phi_top, top, left, right, nx, dt, 10, 9.81 / 1000);

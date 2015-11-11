function [ G, top_faces, bottom_faces ] = setup_grid( eta, h, nx, nz )
% SETUP_GRID Set up grid for wave computation.
%
%   [ G, top_faces, bottom_faces ] = SETUP_GRID(eta, h, nx, nz) produces
%   a deformed Cartesian grid of (non-deformed) size [0, 1]x[0, 1] with
%   nx nodes along the x-axis and nz points along the z-axis. The grid
%   is deformed acording to the functions eta(x) and h(x), where eta
%   is a shape function for the surface, while h is a shape function
%   for the seabed. Note that h(x) is already shifted by 1 for convenience.
%
%   Example usage:
%
%       eta = @(x) 0.1 * sin(2 * pi*x);
%       h =   @(x) 0.05 * cos(3 * pi * x);
%       [G, top_faces, bottom_faces] = setup_grid(eta, h, 100, 100);
%       plotGrid(G);

require mimetic;

G = cartGrid([nx, nz], [1, 1]);
G = computeGeometry(G);

% Identify top and bottom faces. Note that direct comparison works here
% because 32-bit integers can be exactly represented as doubles,
% and so it will only work as long as the cartesian grid is instantiated
% as [0, 1]x[0, 1] or similar.
top_faces = (1 : G.faces.num)';
top_faces = top_faces(G.faces.centroids(:, 2) == 1);
bottom_faces = (1 : G.faces.num)';
bottom_faces = bottom_faces(G.faces.centroids(:, 2) == 0);

% Deform the grid according to supplied h and eta
h_shifted = @(x) ones(numel(x), 1) + h(x);
x = G.nodes.coords;
xx = zeros(size(x));
x1 = x(:, 1);
xx(:, 1) = x1;
xx(:, 2) = -h_shifted(x1) + (eta(x1) + h_shifted(x1)).*x(:, 2);
G.nodes.coords = xx;

end


function [ G, surface_faces, surface_nodes ] = setup_grid(L, D, eta, h, nx, nz )
% SETUP_GRID Set up grid for wave computation.
%
%   [ G, top, bottom, left, right ] = SETUP_GRID(L, D, eta, h, nx, nz) produces
%   a deformed Cartesian grid of (non-deformed) size [0, L]x[-D, 0] with
%   nx nodes along the x-axis and nz points along the z-axis. The grid
%   is deformed according to the functions eta(x) and h(x), where eta
%   is a shape function for the surface, while h is a shape function
%   for the seabed. Here we define a 'shape function' such that it only
%   describes the deviation from the zero surface, i.e. deviation from the
%   case when h = 0. The depth is thus described by D + h(x).
%   
%   The output G represents the geometry, while top, bottom, left and right
%   contain the indices of the top, bottom, left and right faces,
%   respectively.
%

require mimetic;

G = cartGrid([nx, nz], [L, D]);
G = computeGeometry(G);

all_faces = 1 : G.faces.num;
all_nodes = 1 : G.nodes.num;
surface_faces = all_faces(G.faces.centroids(:, 2) == D);
surface_nodes = all_nodes(G.nodes.coords(:, 2) == D);

% Deform the grid according to supplied h and eta
h_shifted = @(x) D * ones(size(x)) + h(x);
x = G.nodes.coords;
xx = zeros(size(x));
x1 = x(:, 1);
xx(:, 1) = x1;
xx(:, 2) = -h_shifted(x1) + (eta(x1) + h_shifted(x1)) .* x(:, 2) / D;
G.nodes.coords = xx;
G = computeGeometry(G);

end


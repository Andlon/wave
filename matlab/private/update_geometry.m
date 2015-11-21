function [ G, top, bottom, left, right ] = update_geometry( G, surface_nodes, surface_faces, surface_shape, h, nx, nz )

% Here we have new values for the z-position of the face centroids.
% We want to determine the z-positions of the nodes from these, however we
% have one problem: There are N + 1 nodes but only N centroids, so
% the node positions cannot be uniquely determined. 
% We can however impose some constraints and solve an optimization problem.
% In particular, we need to make sure that the new nodes result in the same
% face centroid positions as surface_shape. Denoting z_i the z-position of 
% node i and b(i) the z-position of face centroid i, 
% this gives the following linear system:
% 
%   z(i) + z(i + 1) = 2 * b(i) for i = 1, ..., N
% 
% At this point, z(N + 1) is free and the preceding gives a linear system
% of the form Az = b where A is N x (N + 1), so the system is
% underdetermined and has infinitely many solutions. We assume our surface
% to be a fairly smooth function, so it makes sense to select solutions so
% that the node z-values are close to each other (hence we avoid jagged
% functions), i.e. for any two nodes i and i + 1, we want to minimize
% (- z_i + z(i + 1))^2. In addition, we want to approximately estimate
% z(N + 1) so that it is close to the linear approximation given by the two
% preceding nodes z(i) and z(i - 1). Thus we also want to minimize
% (- z(N - 1) + 2 z(N) - z(N + 1))^2. The preceding gives a constrained
% least squares optimization problem in the form
% 
% min 1/2 ||Cz||^2 such that Ax = b
%
% which is what we solve for in the following. 

% where n(i) is the z-value of node i, c(i) is the z-value of the centroid
% of face i, and n(N + 1) is known.

N = numel(surface_faces);
A = spdiags([ones(N + 1, 1), ones(N + 1, 1)], [0, 1], N, N + 1);
A = full(A);
b = 2 * surface_shape;

C = spdiags([-ones(N + 1, 1), ones(N + 1, 1)], [0, 1], N + 1, N + 1);
C(end, end - 2:end) = [ -1, 2, -1 ];
C = full(C);

% Disable output from the optimization toolbox to avoid spam. Enable if
% debugging.
options = optimoptions('lsqlin', 'Display', 'off');
x0 = [ surface_shape; surface_shape(end) ];
warning('off', 'optimlib:lsqlin:LinConstraints');
[Z,~,~, EXITFLAG] = lsqlin(C, zeros(N + 1, 1), [], [], ...
    A, b, [], [], x0, options);
warning('on', 'optimlib:lsqlin:LinConstraints');

assert(EXITFLAG == 1);

dx = 1 / nx;
nodefunc = @(x) Z(round(1 + x / dx));
[G, top, bottom, left, right] = setup_grid(nodefunc, h, nx, nz);

end


function [ G, top, bottom, left, right ] = update_geometry( G, surface_nodes, surface_faces, surface_shape, h, nx, nz )

% Here we have new values for the z-position of the face centroids.
% We want to determine the z-positions of the nodes from these, however we
% have one problem: There are N + 1 nodes but only N centroids, so
% the node positions cannot be uniquely determined. Make the assumption
% that the _last_ node moves up or down in the same amount as the last
% centroid. In this case, we can build a system of equations such that for
% each centroid i, we have
%   n(i + 1) + n(i) = 2 * c(i) for i = 1, ..., N - 1
% and
%   n(N) = 2 c(i) - n(N + 1)
% where n(i) is the z-value of node i, c(i) is the z-value of the centroid
% of face i, and n(N + 1) is known.

N = numel(surface_faces);

last_face = surface_faces(end);
last_face_delta = surface_shape(end) - G.faces.centroids(last_face, 2);
last_node = surface_nodes(end);
last_node_z = G.nodes.coords(last_node, 2);
last_node_z = last_node_z + last_face_delta;

A = spdiags([ones(N, 1), ones(N, 1)], [0, 1], N, N);
b = 2 * surface_shape;
b(end) = b(end) - last_node_z;

nodepos = A \ b;
nodepos(end + 1) = last_node_z;
dx = 1 / nx;
nodefunc = @(x) nodepos(min(N, round(1 + x / dx)));

[G, top, bottom, left, right] = setup_grid(nodefunc, h, nx, nz);

%face_centroid_diff = top_faces_shape - G.faces.centroids(top_faces, 2);



% % So far we have the temporal difference in face centroids along the top.
% % We need to obtain the temporal difference in the node points.
% face_x = G.faces.centroids(top_faces, 1);
% node_x = G.nodes.coords(top_nodes, 1);
% node_diff = interp1(face_x, face_centroid_diff, node_x, 'linear', 'extrap');
% 
% % Compute new values for the nodes
% node_y = G.nodes.coords(top_nodes, 2) + node_diff;
% eta = @(x) interp1(node_x, node_y, x);
% [G, top, bottom, left, right] = setup_grid(eta, h, nx, nz);

end


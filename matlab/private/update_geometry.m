function [ G, top, bottom, left, right ] = update_geometry( G, top_nodes, top_faces, top_faces_shape, h, nx, nz )

face_centroid_diff = top_faces_shape - G.faces.centroids(top_faces, 2);

% So far we have the temporal difference in face centroids along the top.
% We need to obtain the temporal difference in the node points.
face_x = G.faces.centroids(top_faces, 1);
node_x = G.nodes.coords(top_nodes, 1);
node_diff = interp1(face_x, face_centroid_diff, node_x, 'linear', 'extrap');

% Compute new values for the nodes
node_y = G.nodes.coords(top_nodes, 2) + node_diff;
eta = @(x) interp1(node_x, node_y, x);
[G, top, bottom, left, right] = setup_grid(eta, h, nx, nz);

end


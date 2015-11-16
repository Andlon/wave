function eta = solve_wave( G, phi0_top, top_faces, ...
    left_faces, right_faces, nx, dt, nt, g)
% SOLVE_WAVE(...)

top_cells = boundary_cells(G, top_faces);
top_nodes = sort_boundary_nodes(G, face_nodes(G, top_faces));

phi_top_faces = phi0_top;
dx = 1/nx;

eta = zeros(nx + 1, nt + 1);
eta(:, 1) = G.nodes.coords(top_nodes, 2);

% Solve for initial phi in the entirety of the domain and compute gradients
[phi, v] = solve_laplace(G, phi_top_faces, ...
    top_faces, left_faces, right_faces);
grad_phi_top = cell_gradients(G, v, top_cells);
phi_top = phi(top_cells);

surface_shape = G.faces.centroids(top_faces, 2);

for n = 1:nt
    surface_shape = next_surface_shape(dx, dt, surface_shape, grad_phi_top);
    
    % TODO: Update geometry! We now have values for the centroids
    % of the top faces for the new domain. Interpolate the difference from
    % the old domain to attain new values at the surface nodes
    
    % eta(:, n + 1) = G.nodes.coords(top_nodes, 2);
    
    % At this point we have the shape of our new domain, but
    % only potential data for our old domain. The best we can do is
    % probably to interpolate our old data to fit our new domain.
    % TODO: Implement interpolation onto new domain
    
    phi_top = next_surface_potential(g, dt, surface_shape, ...
        phi_top, grad_phi_top);
    [phi, v] = solve_laplace(G, phi_top, ...
        top_faces, left_faces, right_faces);
    
    grad_phi_top = cell_gradients(G, v, top_cells);
    phi_top = phi(top_cells);
end

end

function top_nodes = sort_boundary_nodes(G, top_nodes)
% Sorts nodes on boundary along their x-axis coordinates
xcoords = G.nodes.coords(top_nodes, 1);
[~, I] = sort(xcoords);
top_nodes = top_nodes(I);
end

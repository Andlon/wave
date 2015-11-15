function eta = solve_wave( G, phi0_top, top_faces, ...
    left_faces, right_faces, bottom_faces, nx, dt, nt, g)
% SOLVE_WAVE(...)

top_nodes = face_nodes(G, top_faces);
phi_top_faces = phi0_top;
dx = 1/(nx + 1);
eta = zeros(nx + 1, nt + 1);
eta(:, 1) = G.nodes.coords(top_nodes, 2);

% Solve for initial phi in the entirety of the domain and compute gradients
[phi, v] = solve_laplace(G, phi_top_faces, ...
    top_faces, left_faces, right_faces);
grad_phi_top = node_gradients(G, v, top_nodes);
phi_top = node_potentials(G, phi, top_nodes);

for n = 1:nt
    eta(:, n + 1) = next_surface_shape(dx, dt, eta(:, n), grad_phi_top);
    
    % At this point we have eta(the shape of our new domain, but
    % only potential data for our old domain. The best we can do is
    % probably to interpolate our old data to fit our new domain.
    % TODO: Implement interpolation onto new domain
    
    % TODO: Update geometry!
    
    phi_top = next_surface_potential(g, dt, eta(:, n + 1), ...
        phi_top, grad_phi_top);
    [phi, v] = solve_laplace(G, phi_top_faces, ...
        top_faces, left_faces, right_faces);
    
    % MAJOR HACK! Pretend that the top face values are equal to the values
    % for the node to its left. Fix asap!
    phi_top_faces = phi_top(1:end-1);
    
    grad_phi_top = node_gradients(G, v, top_nodes);
    phi_top = node_potentials(G, phi, top_nodes);
end

end

function eta = solve_wave( sim, surface_phi0, dt, nt)
% SOLVE_WAVE(...)

dx = 1/sim.Nx;

eta = zeros(sim.Nx + 1, nt + 1);
eta(:, 1) = sim.eta();

% Solve for initial phi in the entirety of the domain and compute gradients
[phi, v] = solve_laplace(sim.grid, surface_phi0, sim.top_faces);

surface_shape = sim.surface_shape();
surface_phi_grad = sim.surface_potential_gradient(v);
surface_phi = sim.surface_potential(phi);

for n = 1:nt
    surface_shape = next_surface_shape(dx, dt, surface_shape, surface_phi_grad);    
    sim.update_surface(surface_shape);
    eta(:, n + 1) = sim.eta();
    
    % At this point we have the shape of our new domain, but
    % only potential data for our old domain. For now we just use the data
    % for our old domain, but this will inevitably introduce an error.
    % Unforuntately we're not likely to be able to extrapolate well either,
    % so it's unclear if we can do any better.
    
    surface_phi = next_surface_potential(sim.g, dt, surface_shape, ...
        surface_phi, surface_phi_grad);
    [phi, v] = solve_laplace(sim.grid, surface_phi, sim.top_faces);
    
    surface_phi_grad = sim.surface_potential_gradient(v);
    surface_phi = sim.surface_potential(phi);
end

end
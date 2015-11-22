function eta = solve_wave( sim, surface_phi0, dt, nt, varargin)
% SOLVE_WAVE(...)
options = process_options(varargin);
dx = 1/sim.Nx;

eta = zeros(sim.Nx + 1, nt + 1);
eta(:, 1) = sim.eta();

% Solve for initial phi in the entirety of the domain and compute gradients
[phi, v] = solve_laplace(sim.grid, surface_phi0, sim.top_faces);

surface_shape = sim.surface_shape();
surface_phi_grad = sim.surface_potential_gradient(v);
surface_phi = sim.surface_potential(phi);

if options.ShowProgress
   w_handle = waitbar(0, 'Computing solution...');
end

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

    % Our values for phi and grad phi are computed at the surface nodes, 
    % but our Laplace solver requires Dirichlet conditions on the faces.
    % Approximate face potentials from these nodal values.
    surface_face_phi = face_potential(surface_phi);       
    [phi, v] = solve_laplace(sim.grid, surface_face_phi, sim.top_faces);
    
    surface_phi_grad = sim.surface_potential_gradient(v);
    surface_phi = sim.surface_potential(phi);
    
    if options.ShowProgress
       waitbar(n / nt, w_handle)
    end
end

if options.ShowProgress
   close(w_handle); 
end

end

function potential = face_potential(node_potential)
% Average neighboring nodes
potential = (node_potential(1:end-1) + node_potential(2:end)) / 2;
end

function opt = process_options(options)
opt.ShowProgress = false;
i = 1;
while i <= length(options)
    switch options{i}
        case 'ShowProgress'
            opt.ShowProgress = true;
            i = i + 1;        
        otherwise
            i = i + 1;
    end    
end
end
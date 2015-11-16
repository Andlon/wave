function eta = next_surface_shape( dx, dt, eta, grad_phi )

% Differentiate eta. For now simply we use forward differences, and
% backward difference for the last element.
% TODO: Better ways to approximate this?
eta_x = [ diff(eta); eta(end) - eta(end-1) ] / dx;
eta = eta - dt * (grad_phi(:, 1) .* eta_x - grad_phi(:, 2));

end


function eta = next_surface_shape( dx, dt, eta, grad_phi )

% Differentiate eta. For now simply we use forward differences, and
% backward difference for the last element.
% TODO: Better ways to approximate this?
eta_x = [ diff(eta); eta(end) - eta(end-1) ] / dx;
A = [ eta_x'; - ones(1, numel(eta_x)) ];
dots = trace(grad_phi * A);
eta = eta - dt * dots;

end


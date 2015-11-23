function eta = next_surface_shape( dx, dt, eta, grad_phi )

eta_x = differentiate(eta, dx);
eta = eta - dt * (grad_phi(:, 1) .* eta_x - grad_phi(:, 2));

end

function Y = differentiate(X, dx)
centraldiff = X(3:end) - X(1:end-2) / (2 * dx);
boundary_forwarddiff = (X(2) - X(1)) / dx;
boundary_backwarddiff = (X(end) - X(end - 1)) / dx;
Y = [boundary_forwarddiff; centraldiff; boundary_backwarddiff];
end


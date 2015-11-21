function phi_top = next_surface_potential(g, dt, eta, phi_top, phi_grad_top)
norm_sqr_phi_grad_top = sum(phi_grad_top .^ 2, 2);
phi_top = phi_top - dt * (0.5 * norm_sqr_phi_grad_top - g * eta);
end


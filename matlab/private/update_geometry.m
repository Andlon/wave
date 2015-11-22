function [ G, top, bottom, left, right ] = update_geometry( surface_shape, h, nx, nz )
dx = 1 / nx;
nodefunc = @(x) surface_shape(round(1 + x / dx));
[G, top, bottom, left, right] = setup_grid(nodefunc, h, nx, nz);
end


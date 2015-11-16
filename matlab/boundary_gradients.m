function grad = boundary_gradients( G, v, boundary_faces)
neighbors = G.faces.neighbors(boundary_faces, :);
sorted_neighbors = sort(neighbors, 2);
cells = sorted_neighbors(:, 2);

% Additional layer of debugging to catch errors. Comment out for 
% increased performance. Checks that all faces are actually on the
% boundary by asserting that each face is only connected to one cell.
assert(~all(sorted_neighbors(:, 1)), 'faces must be on boundary!');

grad = cell_gradients(G, v, cells);
end


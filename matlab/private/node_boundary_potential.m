function potential = node_boundary_potential( G, phi, nodes )
N = numel(nodes);
potential = zeros(N, 1);

% Reconstruct the mapping from G.faces.nodes to global faces
global_faces = rldecode(1:G.faces.num, diff(G.faces.nodePos), 2) .';

% For each node
for i = 1:numel(nodes)
    node = nodes(i);
    
    local_face_is_node_neighbor = G.faces.nodes == node;
    faces = global_faces(local_face_is_node_neighbor);
    faces = faces(is_boundary_face(G, faces));
    cells = unique(boundary_cells(G, faces));
    
    cellpotentials = phi(cells);
    
    % Average the gradients of each neighboring cells
    potential(i) = sum(cellpotentials) / numel(cellpotentials);
end


end

function is_boundary = is_boundary_face(G, faces)
sorted_neighbors = sort(G.faces.neighbors(faces, :), 2);
is_boundary = sorted_neighbors(:, 1) == 0;
end




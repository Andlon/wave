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
    cells = boundary_cells(G, faces);
    
    cellpotentials = phi(cells);
    
    % Average the gradients of each neighboring cells
    potential(i) = sum(cellpotentials) / numel(cellpotentials);
end


end


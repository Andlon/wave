function grad = node_boundary_gradients( G, v, nodes)
N = numel(nodes);
grad = zeros(N, 2);
% Reconstruct the mapping from G.faces.nodes to global faces
global_faces = rldecode(1:G.faces.num, diff(G.faces.nodePos), 2) .';

for i = 1:numel(nodes)
    node = nodes(i);
    
    local_face_is_node_neighbor = G.faces.nodes == node;
    faces = global_faces(local_face_is_node_neighbor);
    cells = boundary_cells(G, faces);
    
    cellgrads = cell_gradients(G, v, cells);
    
    % Average the gradients of each neighboring cells
    grad(i, :) = sum(cellgrads, 1) / size(cellgrads, 1);
end

end


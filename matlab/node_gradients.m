function grad = node_gradients( G, v, nodes )
N = numel(nodes);
grad = zeros(N, 2);
for i = 1:N
   node = nodes(i);
   cells = node_cells(G, node);
   cellgrad = cell_gradients(G, v, cells);
   num_cellgrads = size(cellgrad, 1);
   
   % Approximate gradient at node by averaging all neighboring cells
   grad(i, :) = sum(cellgrad, 1) / num_cellgrads;
end

end

function cells = node_cells(G, node)
faces = node_faces(G, node);
cells = face_cells(G, faces);
end

function faces = node_faces(G, node)
% Note: The following way to obtain the faces associated with a node
% is an extremely inefficient hack, since it scans the entire G.faces.nodes
% matrix for every lookup, but it does the job. Fix if
% it becomes a bottleneck (i.e. precompute a mapping?)
global_faces = rldecode(1:G.faces.num, diff(G.faces.nodePos), 2) .';
faces_contain_node = G.faces.nodes == node;
faces = global_faces(faces_contain_node);

end

function cells = face_cells(G, faces)
cells = nonzeros(unique(G.faces.neighbors(faces, :)));
end


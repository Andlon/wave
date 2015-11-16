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
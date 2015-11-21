function nodes = face_nodes( G, faces )
% FACE_NODES Retrieve all nodes that are connected to the given faces.

% This is highly inefficient. Fix if bottleneck.
nodes = [];
for i = 1:numel(faces)
    face = faces(i);
    nodes = [
        nodes;
        G.faces.nodes(G.faces.nodePos(face) : G.faces.nodePos(face+1)-1, :)
        ];    
end

nodes = unique(nodes);

end


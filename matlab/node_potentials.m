function phi = node_potentials( G, celldata, nodes )
N = numel(nodes);
phi = zeros(N, 1);
for i = 1:N
   node = nodes(i);
   cells = node_cells(G, node);
   
   % Approximate potential at node by averaging all neighboring cells
   phi(i, :) = mean(celldata(cells));
end


end


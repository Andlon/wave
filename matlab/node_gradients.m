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

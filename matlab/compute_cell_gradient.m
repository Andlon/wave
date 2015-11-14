function grad = compute_cell_gradient( G, v, cell )
dimension = G.griddim;
half_faces = G.cells.facePos(cell) : G.cells.facePos(cell + 1) - 1;
faces = G.cells.faces(half_faces);
weighted_face_normals = G.faces.normals(faces, :);
areas = G.faces.areas(faces, :);
face_neighbors = G.faces.neighbors(faces, :);

% Compute the half face normals. These are +- face normals, so we only
% need to determine the sign (relative to face normals) to obtain
% the half face normals.
% The signs are computed as follows:
% dir = cell in first column of face_neighbors => 1, otherwise => -1
half_face_directions = 2 * (face_neighbors(:, 1) == cell) - 1;
D = repmat(half_face_directions, 1, dimension);
A = repmat(areas, 1, dimension);
half_face_normals = D .* weighted_face_normals ./ A;

% TODO: Determine if the following assumption is correct:
% The normals of the half-faces are defined such that they face 
% from the faces into the interior of their respective cells.
% In this case, we must flip the signs of our normals:
half_face_normals = - half_face_normals;

projected_gradients = v(half_faces);
grad = half_face_normals \ projected_gradients;

end


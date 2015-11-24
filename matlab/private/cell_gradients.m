function gradients = cell_gradients( G, v, cells)
gradients = zeros(numel(cells), 2);
for c = 1:numel(cells)
    cell = cells(c);
    gradients(c, :) = compute_cell_gradient(G, v, cell);
end
end

function grad = compute_cell_gradient( G, v, cell )
dimension = G.griddim;
half_faces = G.cells.facePos(cell) : G.cells.facePos(cell + 1) - 1;
faces = G.cells.faces(half_faces, :);

% Determine normals of North, West and East faces
faceN = faces(faces(:, 2) == 4, 1);
faceW = faces(faces(:, 2) == 1, 1);
faceE = faces(faces(:, 2) == 2, 1);
hfaceN = half_faces(faces(:, 1) == faceN);
hfaceW = half_faces(faces(:, 1) == faceW);
hfaceE = half_faces(faces(:, 1) == faceE);
areaN = G.faces.areas(faceN);
areaW = G.faces.areas(faceW);
areaE = G.faces.areas(faceE);

faceN_normal = G.faces.normals(faceN, :) ./ repmat(areaN, 1, dimension);
faceW_normal = G.faces.normals(faceW, :) ./ repmat(areaW, 1, dimension);
faceE_normal = G.faces.normals(faceE, :) ./ repmat(areaE, 1, dimension);

% Determine directions of half face normals
face_neighbors = G.faces.neighbors([faceN, faceW, faceE], :);

% Compute the half face normals. These are +- face normals, so we only
% need to determine the sign (relative to face normals) to obtain
% the half face normals.
% The signs are computed as follows:
% dir = cell in first column of face_neighbors => 1, otherwise => -1
half_face_directions = 2 * (face_neighbors(:, 1) == cell) - 1;

hfaceN_normal = half_face_directions(1) * faceN_normal;
hfaceW_normal = half_face_directions(2) * faceW_normal;
hfaceE_normal = half_face_directions(3) * faceE_normal;

% assert(dot(hfaceW_normal, hfaceE_normal) < 0);

gradN = (v(hfaceN) / areaN) * hfaceN_normal;
gradW = (v(hfaceW) / areaW) * hfaceW_normal;
gradE = (v(hfaceE) / areaE) * hfaceE_normal;

% Orthogonalize W and E relative to N (but not relative to each other)
u1 = gradN;
projection = @(v) u1 .* dot(v, u1) ./ dot(u1, u1);
orthW = gradW - projection(gradW);
orthE = gradE - projection(gradE);

% Take the average of the W and E vectors to create an orthogonal basis
% consisting of u1 and u2. This way we take equal contributions from W and
% E, but retain all the information from N.
u2 = (orthW + orthE) / 2;

if dot(u1, u1) > 0
    grad = u1 + u2;
else
    grad = (gradN + gradW + gradE) / 3;
end

end


function [ phi, v ] = solve_laplace( G, top_phi, top_faces)
%SOLVE_LAPLACE Solve the Laplace equation on the supplied geometry.
% TODO: Write more documentation
%
% Note: The code here is mostly copied from the supplied example code.

nc = G.cells.num;
nf = G.faces.num;
half_faces = G.cells.faces(:, 1);
nhf = numel(half_faces);

% Need to decide on boundary conditions here.
%dirich_faces = [ top, left, right ];
%dirich_p = [ top_phi, zeros(1, numel(left)), zeros(1, numel(right)) ];
dirich_faces = top_faces;
dirich_p = top_phi;

% Compute the mimetic scalar product
rock.perm = ones(G.cells.num, 1);
s   = computeMimeticIP(G, rock);
BI = s.BI;

% Assemble the matrices C and D.
cellNo = rldecode(1:G.cells.num, diff(G.cells.facePos), 2) .';
C = - sparse(1:numel(cellNo), cellNo, 1);
D = sparse(1:numel(cellNo), double(half_faces), 1, numel(cellNo), G.faces.num);

% We identify the flux unknown corresponding to half faces where the Dirichlet
% condition is applied.
is_dirich_faces = false(nf, 1);
is_dirich_faces(dirich_faces) = true;

is_ext_faces  = (G.faces.neighbors(:, 1) == 0) | (G.faces.neighbors(:, 2) == 0) ;
is_int_faces = ~is_ext_faces;
nif = nnz(is_int_faces);
is_neumann_faces = is_ext_faces & ~is_dirich_faces;
is_neumann_half_faces = is_neumann_faces(half_faces);
neuman_half_faces = 1 : nhf;
neuman_half_faces = neuman_half_faces(is_neumann_half_faces);
nnhf = nnz(is_neumann_half_faces);

% We construct the right-hand side corresponding to the source term coming from the
% Dirichlet boundary conditons.
dirich_pii_rhs = zeros(nf, 1); % called pii instead of pi
dirich_pii_rhs(dirich_faces) = dirich_p; 
dirich_rhs = -D*dirich_pii_rhs;

% Reduce the system to the unknown variables.
D = D(:, is_int_faces);
N = - sparse(neuman_half_faces, 1 : nnhf , 1, nhf, nnhf);

% Compute the right-hand side using the dirichlet conditions
rhs = [dirich_rhs; zeros(nc + nif + nnhf, 1)];

% Schur reduction
R = [[-C'; -D'; -N']*BI, eye(nc + nif + nnhf)];
A = [[C, D, N]; zeros(nc + nif + nnhf)];
Q = R*A;
rhs = R*rhs;

% Solve the system.
sol = Q\rhs;

% Recover phi, pi
phi = sol([true(nc, 1); false(nif, 1); false(nnhf, 1)]);
pi = sol([false(nc, 1); true(nif, 1); false(nnhf, 1)]);
phi_neum = sol([false(nc, 1); false(nif, 1); true(nnhf, 1)]);

% Compute v
v = BI * (dirich_rhs - C * phi - D * pi - N * phi_neum);

end


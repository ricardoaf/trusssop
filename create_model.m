function mdl = create_model (x, conn, supp, disp, force, E, A, yield)

%% Default values for bar strength
if nargin<8 || isempty(yield), yield = 1e+100*ones(1,size(conn,1)); end

%% Number of dimensions, bars, DOFs
ndim = size(supp, 2)-1; % number of dimensions
mdl.ne = size(conn, 1); %number of bars
mdl.ndof = length(x); % number of dofs

%% Nodal coords and element connectivities
mdl.x = x(:); % (ndof x 1)
node2dofs = @(n) ndim*n*ones(1,ndim) - (ndim-1:-1:0);
mdl.conn = [node2dofs(conn(:,1)) node2dofs(conn(:,2))]; % (ne x 2ndim)

%% Prescribed displacements, external forces
mdl.u = disp(:); % (ndof x 1)
mdl.f = force(:); % (ndof x 1)

%% Free DOFs
supp2dofs = @(s) nonzeros(node2dofs(s(1)).*s(2:end))';
supp_dofs = [];
for s = supp', supp_dofs = [supp_dofs supp2dofs(s')]; end
mdl.free = setdiff(1:mdl.ndof, supp_dofs); % (1 x nfree)

%% Young modulus, areas, strengths
mdl.E = E(:);
mdl.a = A(:);
mdl.strength = abs(yield(:));

%% Calc bar lengths
mdl.Le = zeros(mdl.ne,1);
for e = 1:mdl.ne
    % element nodes coords and lengths
    c = mdl.conn(e,:); xe = mdl.x(c)'; n = length(xe);
    xe_i = xe(1:n/2); xe_j = xe(n/2+1:end);
    mdl.Le(e) = sqrt(sum((xe_j-xe_i).^2));
end
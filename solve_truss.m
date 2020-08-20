function [displacements, stress, internal_force] = solve_truss(mdl, a)
%% global stiffness matrix
K = zeros(mdl.ndof);

% init transformation matrices
L = []; L{mdl.ne} = [];

% element loop
for e = 1:mdl.ne
    % element nodes coords
    c = mdl.conn(e,:); xe = mdl.x(c)'; n = length(xe);
    xe_i = xe(1:n/2); xe_j = xe(n/2+1:end); dx = xe_j-xe_i;
    
    % internal-to-local transformation (L)
    L{e} = zeros(2, n);
    L{e}(1,1:n/2) = dx./mdl.Le(e);
    L{e}(2,n/2+1:end) = dx./mdl.Le(e);
    
    % local stiffness matrix
    k = mdl.E(e)/mdl.Le(e) * [1 -1; -1 1];
    k0 = L{e}' * k * L{e};
    
    % assembling global matrix
    K(c,c) = K(c,c) + a(e)*k0;
end

%% partition method
free = mdl.free;
fix = setdiff(1:mdl.ndof, free);

% solve unknown displacements and reaction forces
mdl.u(free) = K(free,free)\(mdl.f(free) - K(free,fix)*mdl.u(fix));
mdl.f(fix) = K(fix,free)*mdl.u(free) + K(fix,fix)*mdl.u(fix);

%% calculate element strain, stress and internal force
stress = zeros(mdl.ne,1);
strain = zeros(mdl.ne,1);
internal_force = zeros(mdl.ne,1);

% element loop
for e = 1:mdl.ne
    
    % element nodes displacements and elongation
    c = mdl.conn(e,:); ue = mdl.u(c);
    uint = L{e}*ue; du = uint(2)-uint(1);
    
    % calc element strain and stress
    strain(e) = du/mdl.Le(e);
    stress(e) = mdl.E(e)*strain(e);
    
    % calc internal force
    internal_force(e) = stress(e) * a(e);
end

%% set results
displacements = mdl.u;

%%
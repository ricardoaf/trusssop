function [x, fmin, g] = truss_sop (input)
% truss-based size optimization
% volume minimization with strength and local-buckling constraints

% read (max stiff) truss from file
mdl = input();

% setup
f = @(a) objfcn (a, mdl); % objective function
x0 = mdl.a; % initial guess (design variable)
lb = 1e-4*mdl.a; ub = 1e+4*mdl.a; % lower/upper bounds
c = @(a) nl_con(a, mdl); % strength and local-buckling constraints

% optimization options
opt = optimoptions( 'fmincon', ...
                    'Algorithm', 'interior-point', ...
                    'Display', 'iter', ...
                    'SpecifyObjectiveGradient', true, ...
                    'StepTolerance', eps, ...
                    'ConstraintTolerance', eps, ...
                    'OptimalityTolerance', eps, ...
                    'MaxIterations', 10000, ...
                    'MaxFunctionEvaluations', 200000);

% call fmincon and set outputs
tic; x = fmincon(f, x0, [], [], [], [], lb, ub, c, opt); toc
fmin = f(x); g = c(x);

%% OBJECTIVE FUNCTION
function [f, g] = objfcn (a, mdl)

% minimize volume: f = volume
f = mdl.Le'*a;

% g = df/da
g = mdl.Le;

%% CONSTRAINTS
function [c, ceq] = nl_con(a, mdl)

% nonlinear inequalities
c = [G1(a, mdl); G2(a, mdl)]; 

% nonlinear equalities
ceq = [];

%% CONSTRAINT G1: yield condition
function c = G1(a, mdl)

% solve truss and find stress on bars
[~, sig] = solve_truss(mdl, a);

% yield condition for tensile bars
c = sig - mdl.strength; 

%% CONSTRAINT G2: local-buckling
function c = G2(a, mdl)

% solve truss and find internal force on bars
[~, ~, fi] = solve_truss(mdl, a); 

% area moment of inertia (full circular section)
I = 1/(4*pi)*(a.^2);

% critical load
pcr = pi*pi*mdl.E.*I./(mdl.Le.^2);

% buckling condition for compression bars
c = -fi - pcr;
%%
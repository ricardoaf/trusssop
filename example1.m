function mdl = example1
% spatial truss

%% input data
% nodal coord, connectivities and supports
nodal_coords = [0 0 6, 3 4 3, 0 0 0, 6 0 0, 6 0 6]; % [ndof]
bar_connectivities = [1 2; 3 2; 4 2; 5 2]; % [nbars x 2]
supports = [1, 1 1 1; 3, 1 1 1; 4, 1 1 1; 5, 1 1 1]; % [nsupp x 1+ndim]

% prescribed displacements and external forces
displacements = [0 0 0, 0 0 0, 0 0 0, 0 0 0, 0 0 0]; % [ndof]
forces = [0 0 0, 30e+3 -100e+3 0, 0 0 0, 0 0 0, 0 0 0]; % [ndof]

% material properties (for each bar)
young_modulus = 2.05*108*1e+3 * [1 1 1 1]; % [nbars]
area = 2.4e+3 * [1 1 1 1]; % [nbars]
yield = 1e+300 * [1 1 1 1]; % [nbars]

%% fill model object
mdl = create_model (nodal_coords, ...
                    bar_connectivities, ...
                    supports, ...
                    displacements, ...
                    forces, ...
                    young_modulus, ...
                    area, ...
                    yield);

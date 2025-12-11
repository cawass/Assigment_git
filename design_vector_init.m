function [x0, lb, ub] = design_vector_init()
% design_vector_init
% Loads initial values and bounds from CSV files.

    %% --------------------- Required Variable Order ---------------------
    orderedNames = { ...
        'Mcr','hcr','LambdaLEin','cr_in','s_out','perc_LE','perc_TE', ...
        'Au1','Au2','Au3','Au4','Au5', ...
        'Al1','Al2','Al3','Al4','Al5' ...
    };

    %% -------------------------- Load CSVs ------------------------------
    T0 = readtable('design_initial.csv');     % Initial values
    Tb = readtable('design_bounds.csv');      % Lower and upper bounds

    %% ---------------------- Allocate vectors ---------------------------
    n = length(orderedNames);
    x0 = zeros(n,1);
    lb = zeros(n,1);
    ub = zeros(n,1);

    %% ---------------------- Build vectors ------------------------------
    for i = 1:n
        var = orderedNames{i};

        % Initial value
        row0 = strcmp(T0.Name, var);
        if ~any(row0)
            error('Initial value for variable "%s" not found in design_initial.csv.', var);
        end
        x0(i) = T0.Value(row0);

        % Bounds
        rowb = strcmp(Tb.Name, var);
        if ~any(rowb)
            error('Bounds for variable "%s" not found in design_bounds.csv.', var);
        end
        lb(i) = Tb.LB(rowb);
        ub(i) = Tb.UB(rowb);
    end
end

function [x0, lb, ub] = design_vector_init()
% design_vector_init
% Initial design vector and bounds (Table 1.2 of your report).

% Initial values
Mcr0        = 0.7877;
hcr0        = 10668;      % [m]
LambdaLEin0 = 26.3;       % [deg]
cr_in0      = 5.996;      % [m]
s_out0      = 9.101;      % [m]
perc_LE0    = 17.5;       % [% of chord]
perc_TE0    = 57.5;       % [% of chord]

% CST coefficients (RAE 2822 fit)
Au0 = [ 0.10  0.20  0.10  0.30  0.20 ];
Al0 = [-0.10 -0.20 -0.13 -0.20 -0.10];  % sign chosen consistent with other lower values

% Assemble initial vector
x0 = [Mcr0, hcr0, LambdaLEin0, cr_in0, s_out0, perc_LE0, perc_TE0, Au0, Al0].';

% Lower bounds
lb = [ ...
    0.709,   ...  % Mcr
    9619,   ...   % hcr [m]
    0.01,   ...   % LambdaLEin [deg]
    5.11,   ...   % cr_in [m]
    7.48,   ...   % s_out [m]
    15.0,   ...   % %LE
    55.0,   ...   % %TE
    0.05, 0.10, 0.05, 0.15, 0.10, ...         % Au1..Au5
   -0.15,-0.30,-0.185,-0.30,-0.15  ...        % Al1..Al5
    ].';

% Upper bounds
ub = [ ...
    0.866,   ...  % Mcr
    11735,  ...   % hcr [m]
    37.5,   ...   % LambdaLEin [deg]
    6.10,   ...   % cr_in [m]
    9.88,   ...   % s_out [m]
    20.0,   ...   % %LE
    60.0,   ...   % %TE
    0.15, 0.30, 0.15, 0.45, 0.30, ...         % Au1..Au5
   -0.05,-0.10,-0.075,-0.10,-0.05 ...         % Al1..Al5
    ].';

end

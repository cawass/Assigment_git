function const = constants_ae4205()
% constants_ae4205
% Returns a struct with all constant / reference data used in the MDO
% formulation for the Fokker 100 assignment.

% Geometry
const.S          = 93.50;    % Wing area [m^2]
const.sin        = 4.939;    % Inboard span [m]
const.lambda_ref = 0.235;    % Reference taper ratio [-]

% Atmosphere / gas properties
const.gamma = 1.4;           % Ratio of specific heats
const.Rgas  = 287.058;       % Gas constant [J/(kg*K)]
const.T0    = 288.15;        % Sea-level ISA temperature [K]
const.Lapse = 0.0065;        % Temperature lapse rate [K/m]
const.p0    = 101325;        % Sea-level pressure [Pa]
const.g0    = 9.80665;       % Gravity [m/s^2]

% Reference cruise condition
const.Vcr_ref = 233.56;      % Reference cruise speed [m/s]
const.hcr_ref = 10668;       % Reference cruise altitude [m]

% Aerodynamic “aircraft minus wing” drag increment (to be set by you)
const.CD_AminusW = 0.0;      % Placeholder [-]

% Fuel / propulsion
const.CT_baseline = 1.8639e-4;   % Baseline SFC C_T [N/(N*s)]
const.ftank       = 0.93;        % Tank volume factor [-]
const.rho_f       = 0.81715e3;   % Fuel density [kg/m^3]

% Weights
const.W_f_ref     = 10731;   % Reference fuel weight [kg]
const.WA_minus_W  = 0.0;     % "Aircraft minus wing" weight [kg] (fill from assignment)
const.W_TO_max_ref = 45810;  % Reference MTOW [kg]

% Wing loading constraint
const.WS_ref = NaN;          % Reference wing loading [kg/m^2] – set from assignment

% Structural / material constants (used later in EMWET, etc.)
const.E_al    = 70e3;        % Aluminium Young's modulus [N/mm^2]
const.sigma_y = 295;         % Yield stress [N/mm^2]
const.rho_al  = 2800;        % Aluminium density [kg/m^3]

% Other
const.n_max   = 2.5;         % Load factor
const.rib_pitch = 0.5;       % Rib pitch [m]

end

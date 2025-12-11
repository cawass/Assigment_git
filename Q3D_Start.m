clear; close all; clc;

%% ============================================================
% 1. LOAD CONSTANTS + DESIGN VECTOR
% ============================================================
[x0, lb, ub] = design_vector_init();
const = constants_ae4205();

Mcr = x0(1);
hcr = x0(2);
LambdaLEin_c = x0(3);
cr_in_c      = x0(4);
s_out_c      = x0(5);
perc_LE_c    = x0(6)/100;
perc_TE_c    = x0(7)/100;
Au = x0(8:12)';
Al = x0(13:17)';

% Inboard span
if isfield(const,'s_in')
    s_in = const.s_in;
else
    s_in = const.sin;
end

% Sweep difference from reference
dLambda_LE = 5.9;

%% ============================================================
% 2. CREATE WING OBJECT
% ============================================================
W = Wing(Mcr, hcr, LambdaLEin_c, cr_in_c, s_in, s_out_c, ...
         perc_LE_c, perc_TE_c, Au, Al, dLambda_LE);

%% ============================================================
% 3. GEOMETRY: PLANFORM + AIRFOIL SECTIONS
% ============================================================
AC.Wing.Geom = W.make_AVL_geom();          % [x_LE y_LE z_LE chord twist]
AC.Wing.inc  = 0;                          % incidence

[AC.Wing.Airfoils, AC.Wing.eta] = W.make_airfoil_table_2();

disp('AC.Wing.Geom = ');     disp(AC.Wing.Geom);
disp('AC.Wing.Airfoils = '); disp(AC.Wing.Airfoils);
disp('AC.Wing.eta = ');      disp(AC.Wing.eta);

%% ============================================================
% 4. VISCOUS / INVISCID
% ============================================================
AC.Visc  = 0;
AC.Aero.MaxIterIndex = 150;

%% ============================================================
% 5. FLIGHT CONDITIONS (CALCULATED FROM INPUTS)
% ============================================================
% Atmosphere (ISA) using constants file
if const.Lapse ~= 0
    T = const.T0 - const.Lapse * hcr;
    P = const.p0 * (T/const.T0)^(const.g0/(const.Rgas*const.Lapse));
else
    T = const.T0;
    P = const.p0 * exp(-const.g0*hcr/(const.Rgas*T));
end

rho = P / (const.Rgas * T);
a   = sqrt(const.gamma * const.Rgas * T);

V = Mcr * a;               % true airspeed
MAC = const.MAC;          % 3.80 m
mu = 1.789e-5;            % viscosity

Re = rho * V * MAC / mu;

% Required CL from design weight
W_required = const.W_des * const.g0;
S = const.S;
CL_required = (2 * W_required) / (rho * V^2 * S);

% Fill AC struct
AC.Aero.V   = V;
AC.Aero.rho = rho;
AC.Aero.alt = hcr;
AC.Aero.Re  = Re;
AC.Aero.M   = Mcr;
AC.Aero.CL  = CL_required;

fprintf('\n===== FLIGHT CONDITION =====\n');
fprintf('Altitude       %.1f m\n', hcr);
fprintf('Temperature    %.2f K\n', T);
fprintf('Density        %.4f kg/m³\n', rho);
fprintf('Speed of sound %.2f m/s\n', a);
fprintf('Airspeed       %.2f m/s\n', V);
fprintf('Reynolds       %.3e\n', Re);
fprintf('Required CL    %.4f\n\n', CL_required);

%% ============================================================
% 6. RUN Q3D SOLVER
% ============================================================
tic;
Res = Q3D_solver(AC);
toc;

disp('Q3D Solver Output Fields:');
disp(fieldnames(Res.Wing))


plot_wing_3D(W, Res, AC)


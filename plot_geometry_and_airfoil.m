clear; close all; clc;

%% ============================================================
% 1. LOAD CONSTANTS + DESIGN VECTOR + BOUNDS
% ============================================================
const = constants_ae4205();                 % Already loads geometry constants
[x0, lb, ub] = design_vector_init();        % Loads initial values + bounds

% Variable ordering inside design vector:
% x = [Mcr, hcr, LambdaLEin, cr_in, s_out, perc_LE, perc_TE, Au1..Au5, Al1..Al5]
LambdaLEin = x0(3);
cr_in      = x0(4);
s_out      = x0(5);
perc_LE    = x0(6);
perc_TE    = x0(7);

% Extract CST coefficients
Au      = x0(8:12)';
Al      = x0(13:17)';
Au_min  = lb(8:12)';
Al_min  = lb(13:17)';
Au_max  = ub(8:12)';
Al_max  = ub(13:17)';

%% ============================================================
% 2. WING PLANFORM PLOT (USING CONSTANTS FROM FUNCTION)
% ============================================================

% Geometry from constants_ae4205()
S      = const.S;
sin    = const.sin;         % inboard span
lambda = const.lambda_ref;  % taper ratio

% Outboard tip chord from reference taper
c_root = cr_in;
c_tip  = lambda * c_root;

% Leading edge geometry
LE_root     = 0;
LE_kink     = sin * tand(LambdaLEin);

% Outboard leading edge
LE_tip      = LE_kink + s_out * tand(LambdaLEin);

% Trailing edges
TE_root = LE_root + c_root;
TE_kink = LE_kink + (lambda * c_root);     % consistent with taper
TE_tip  = LE_tip  + c_tip;

figure; hold on; grid on; axis equal;
title('Wing Planform from Loaded Constants & Initial Values');
xlabel('x [m]'); ylabel('y [m]');

% Inboard polygon
patch([LE_root TE_root TE_kink LE_kink], ...
      [0       0       sin     sin], ...
      'b', 'FaceAlpha',0.3);

% Outboard polygon
patch([LE_kink LE_kink + c_tip  LE_tip + c_tip  LE_tip], ...
      [sin     sin              sin+s_out      sin+s_out], ...
      'r', 'FaceAlpha',0.3);

legend('Inboard Section','Outboard Section');
hold off;

%% ============================================================
% 3. AIRFOIL PLOTS (BASELINE, LOWER BOUND, UPPER BOUND)
% ============================================================

X = linspace(0,1,200)';

[Xtu,    Xtl]    = D_airfoil2(Au,     Al,     X);
[Xtu_min,Xtl_min]= D_airfoil2(Au_min, Al_min, X);
[Xtu_max,Xtl_max]= D_airfoil2(Au_max, Al_max, X);

figure; hold on; axis equal; grid on;
title('CST Airfoil: Baseline, Lower Bound, Upper Bound');
xlabel('x'); ylabel('z');

% Baseline
plot(Xtu(:,1),    Xtu(:,2),    'b','LineWidth',1.4);
plot(Xtl(:,1),    Xtl(:,2),    'b','LineWidth',1.4);

% Lower bound
plot(Xtu_min(:,1),Xtu_min(:,2),'r--','LineWidth',1.2);
plot(Xtl_min(:,1),Xtl_min(:,2),'r--','LineWidth',1.2);

% Upper bound
plot(Xtu_max(:,1),Xtu_max(:,2),'g--','LineWidth',1.2);
plot(Xtl_max(:,1),Xtl_max(:,2),'g--','LineWidth',1.2);

% RAE2822 reference (same as in assignment)
rae = load('rae2822.dat');
plot(rae(:,1), rae(:,2), 'k', 'LineWidth',1.2);

legend('Baseline upper','Baseline lower', ...
       'Lower bound upper','Lower bound lower', ...
       'Upper bound upper','Upper bound lower', ...
       'RAE2822 Reference');

hold off;


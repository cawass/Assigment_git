clear; close all; clc;

%% ============================================================
% 1. LOAD CONSTANTS + DESIGN VECTOR + BOUNDS
% ============================================================
const = constants_ae4205();
[x0, lb, ub] = design_vector_init();

% Design vector ordering:
% x = [Mcr, hcr, LambdaLEin, cr_in, s_out, perc_LE, perc_TE, Au1..Au5, Al1..Al5]
LambdaLEin_c = x0(3);          % inboard LE sweep [deg] (design variable)
cr_in_c      = x0(4);          % root chord [m]
s_out_c      = x0(5);          % outboard span [m]
perc_LE_c    = x0(6) / 100;    % [% → fraction]
perc_TE_c    = x0(7) / 100;

Au      = x0(8:12)';           % baseline CST upper
Al      = x0(13:17)';          % baseline CST lower
Au_min  = lb(8:12)';
Al_min  = lb(13:17)';
Au_max  = ub(8:12)';
Al_max  = ub(13:17)';

% Get inboard span (support legacy name "sin")
if isfield(const,'s_in')
    s_in_const = const.s_in;
else
    s_in_const = const.sin;  % fallback
end

%% ============================================================
% 2. REFERENCE GEOMETRY (from drawing) & ANGLE DIFFERENCE
% ============================================================
% Reference / original wing (Fokker drawing, all in meters)
cr_in_ref   = 5.996;           % root chord
s_in_ref    = 4.939;           % inboard span
s_out_ref   = 14.040 - 4.939;  % outboard span
LambdaLEin_ref  = 26.3;        % inboard LE sweep [deg]
LambdaLEout_ref = 20.4;        % outboard LE sweep [deg]
LambdaTEin_ref  = 4.2;         % inboard TE sweep [deg]
LambdaTEout_ref = 4.2;         % outboard TE sweep [deg]
perc_LE_ref = 17.5/100;
perc_TE_ref = 57.5/100;

% Difference between inboard and outboard LE sweep (kept constant)
dLambda_LE = LambdaLEin_ref - LambdaLEout_ref;   % = 5.9 deg

% Common TE sweeps for all design cases (as in assignment)
LambdaTEin  = 4.2;
LambdaTEout = 4.2;

% MATLAB default colors (same as your airfoil plot)
blue  = [0 0.4470 0.7410];
red   = [0.8500 0.3250 0.0980];
green = [0.4660 0.6740 0.1880];
black = [0 0 0];

%% ============================================================
% 3. DEFINE CASES: REFERENCE, CURRENT, LOWER, UPPER
% ============================================================
cases = struct([]);

% 1) Reference/original wing
cases(1).name        = 'Reference';
cases(1).cr_in       = cr_in_ref;
cases(1).s_in        = s_in_ref;
cases(1).s_out       = s_out_ref;
cases(1).LambdaLEin  = LambdaLEin_ref;
cases(1).LambdaLEout = LambdaLEout_ref;
cases(1).LambdaTEin  = LambdaTEin_ref;
cases(1).LambdaTEout = LambdaTEout_ref;
cases(1).perc_LE     = perc_LE_ref;
cases(1).perc_TE     = perc_TE_ref;
cases(1).color       = black;
cases(1).style       = '-';
cases(1).isCurrent   = false;

% 2) Current design (from x0)
cases(2).name        = 'Current';
cases(2).cr_in       = cr_in_c;
cases(2).s_in        = s_in_const;
cases(2).s_out       = s_out_c;
cases(2).LambdaLEin  = LambdaLEin_c;
cases(2).LambdaLEout = LambdaLEin_c - dLambda_LE;
cases(2).LambdaTEin  = LambdaTEin;
cases(2).LambdaTEout = LambdaTEout;
cases(2).perc_LE     = perc_LE_c;
cases(2).perc_TE     = perc_TE_c;
cases(2).color       = blue;
cases(2).style       = '- -';
cases(2).isCurrent   = true;

% 3) Lower-bound design (using lb)
LambdaLEin_lo = lb(3);
cr_in_lo      = lb(4);
s_out_lo      = lb(5);
perc_LE_lo    = lb(6)/100;
perc_TE_lo    = lb(7)/100;

cases(3).name        = 'Lower bound';
cases(3).cr_in       = cr_in_lo;
cases(3).s_in        = s_in_const;
cases(3).s_out       = s_out_lo;
cases(3).LambdaLEin  = LambdaLEin_lo;
cases(3).LambdaLEout = LambdaLEin_lo - dLambda_LE;
cases(3).LambdaTEin  = LambdaTEin;
cases(3).LambdaTEout = LambdaTEout;
cases(3).perc_LE     = perc_LE_lo;
cases(3).perc_TE     = perc_TE_lo;
cases(3).color       = red;
cases(3).style       = '--';
cases(3).isCurrent   = false;

% 4) Upper-bound design (using ub)
LambdaLEin_hi = ub(3);
cr_in_hi      = ub(4);
s_out_hi      = ub(5);
perc_LE_hi    = ub(6)/100;
perc_TE_hi    = ub(7)/100;

cases(4).name        = 'Upper bound';
cases(4).cr_in       = cr_in_hi;
cases(4).s_in        = s_in_const;
cases(4).s_out       = s_out_hi;
cases(4).LambdaLEin  = LambdaLEin_hi;
cases(4).LambdaLEout = LambdaLEin_hi - dLambda_LE;
cases(4).LambdaTEin  = LambdaTEin;
cases(4).LambdaTEout = LambdaTEout;
cases(4).perc_LE     = perc_LE_hi;
cases(4).perc_TE     = perc_TE_hi;
cases(4).color       = green;
cases(4).style       = '--';
cases(4).isCurrent   = false;

%% ============================================================
% 4. PRECOMPUTE GEOMETRY FOR ALL CASES
% ============================================================
nCases = numel(cases);

% First case defines struct layout
geoms(1) = compute_planform_case(cases(1));
max_y = geoms(1).y_total;
min_x = min([geoms(1).LE0 geoms(1).LE1 geoms(1).LE2 ...
             geoms(1).TE0 geoms(1).TE1 geoms(1).TE2]);
max_x = max([geoms(1).LE0 geoms(1).LE1 geoms(1).LE2 ...
             geoms(1).TE0 geoms(1).TE1 geoms(1).TE2]);

for k = 2:nCases
    geoms(k) = compute_planform_case(cases(k));
    max_y = max(max_y, geoms(k).y_total);

    min_x = min(min_x, min([geoms(k).LE0 geoms(k).LE1 geoms(k).LE2 ...
                            geoms(k).TE0 geoms(k).TE1 geoms(k).TE2]));
    max_x = max(max_x, max([geoms(k).LE0 geoms(k).LE1 geoms(k).LE2 ...
                            geoms(k).TE0 geoms(k).TE1 geoms(k).TE2]));
end

%% ============================================================
% 5. FIGURE WITH PLANFORM (LEFT) + AIRFOIL (RIGHT)
% ============================================================
figure('Position',[100 100 900 500]);

%% ------------------------------------------------------------
% LEFT: PLANFORM WITH ALL CASES + FUEL TANK FOR CURRENT
% ------------------------------------------------------------
subplot(1,2,1); hold on; grid on; axis equal;
title('Wing Planform: Reference, Design, Lower & Upper Bounds');
xlabel('x [m]'); ylabel('y [m]');

% Dummy handles for legend
href  = plot(NaN,NaN,'-','Color',black,'LineWidth',1.8);
hcur  = plot(NaN,NaN,'-','Color',blue, 'LineWidth',1.8);
hlow  = plot(NaN,NaN,'--','Color',red,  'LineWidth',1.8);
hup   = plot(NaN,NaN,'--','Color',green,'LineWidth',1.8);
htank = patch(NaN,NaN,'y','FaceAlpha',0.3,'EdgeColor','k');

for k = 1:nCases
    p = cases(k);
    g = geoms(k);

    % Closed outline polygon: LE0->LE1->LE2->TE2->TE1->TE0->LE0
    X_outline = [g.LE0 g.LE1 g.LE2 g.TE2 g.TE1 g.TE0 g.LE0];
    Y_outline = [0      p.s_in p.s_in+p.s_out p.s_in+p.s_out p.s_in 0 0];

    plot(X_outline, Y_outline, ...
         'Color', p.color, 'LineStyle', p.style, 'LineWidth', 1.8);

    % Spar lines (LE spar and TE spar)
    y_span = linspace(0, g.y_total, 100);
    x_LE_spar = g.LE_spar_fun(y_span);
    x_TE_spar = g.TE_spar_fun(y_span);

    plot(x_LE_spar, y_span, ':', 'Color', p.color, 'LineWidth', 1.0);
    plot(x_TE_spar, y_span, ':', 'Color', p.color, 'LineWidth', 1.0);

    % Fuel tank only for current design

        b_semi_cur = g.y_total;
        s_tank = 0.85 * b_semi_cur;   % 85% of semi-span
        y_t0 = 0;
        y_t1 = s_tank;

        tank_LE_0 = g.LE_spar_fun(y_t0);
        tank_TE_0 = g.TE_spar_fun(y_t0);
        tank_LE_1 = g.LE_spar_fun(y_t1);
        tank_TE_1 = g.TE_spar_fun(y_t1);

        patch([tank_LE_0 tank_TE_0 tank_TE_1 tank_LE_1], ...
              [y_t0      y_t0      y_t1      y_t1     ], ...
              p.color, ...                 % fill color = case color
              'FaceAlpha', 0.02, ...       % transparency
              'EdgeColor', p.color);       % outline = case color
        
            
end

legend([href hcur hlow hup htank], ...
       {'Reference wing','Current design','Lower bound','Upper bound','Fuel tank'}, ...
       'Location','northoutside');

xlim([min_x-0.5, max_x+0.5]);
ylim([0, max_y + 0.5]);
hold off;

%% ------------------------------------------------------------
% RIGHT: CST AIRFOIL + BOUNDS (COLOR SCHEME MATCHING PLANFORM)
% ------------------------------------------------------------
subplot(1,2,2); hold on; grid on; axis equal;
title('CST Airfoil: Baseline, Lower & Upper Bounds');
xlabel('x'); ylabel('z');

X = linspace(0,1,200)';

[Xtu,    Xtl]     = D_airfoil2(Au,     Al,     X);
[Xtu_min,Xtl_min] = D_airfoil2(Au_min, Al_min, X);
[Xtu_max,Xtl_max] = D_airfoil2(Au_max, Al_max, X);

% Baseline (solid blue)
plot(Xtu(:,1), Xtu(:,2), 'Color', blue,  'LineWidth',1.4);
plot(Xtl(:,1), Xtl(:,2), 'Color', blue,  'LineWidth',1.4);

% Lower bound (dashed red)
plot(Xtu_min(:,1), Xtu_min(:,2), '--', 'Color', red, 'LineWidth',1.2);
plot(Xtl_min(:,1), Xtl_min(:,2), '--', 'Color', red, 'LineWidth',1.2);

% Upper bound (dashed green)
plot(Xtu_max(:,1), Xtu_max(:,2), '--', 'Color', green, 'LineWidth',1.2);
plot(Xtl_max(:,1), Xtl_max(:,2), '--', 'Color', green, 'LineWidth',1.2);

% Reference RAE2822 airfoil (solid black)
rae = load('rae2822.dat');
plot(rae(:,1), rae(:,2), 'Color', black, 'LineWidth',1.1);

legend({'Baseline upper','Baseline lower', ...
        'Lower bound upper','Lower bound lower', ...
        'Upper bound upper','Upper bound lower', ...
        'RAE2822 reference'}, ...
        'Location','southoutside');

hold off;

%% ============================================================
% LOCAL FUNCTION: compute planar wing geometry for one case
% ============================================================
function geom = compute_planform_case(p)
    % Span stations
    y0 = 0;
    y1 = p.s_in;
    y2 = p.s_in + p.s_out;

    % Root LE/TE
    LE0 = 0;
    TE0 = p.cr_in;

    % Kink LE/TE
    LE1 = LE0 + tand(p.LambdaLEin)  * (y1 - y0);
    TE1 = TE0 + tand(p.LambdaTEin)  * (y1 - y0);

    % Tip LE/TE
    LE2 = LE1 + tand(p.LambdaLEout) * (y2 - y1);
    TE2 = TE1 + tand(p.LambdaTEout) * (y2 - y1);

    geom.s_in    = p.s_in;
    geom.s_out   = p.s_out;
    geom.y_total = y2;

    geom.LE0 = LE0; geom.LE1 = LE1; geom.LE2 = LE2;
    geom.TE0 = TE0; geom.TE1 = TE1; geom.TE2 = TE2;

    % Piecewise LE/TE along span
    geom.LE_fun = @(y) (y <= y1) .* (LE0 + tand(p.LambdaLEin)  .* y) + ...
                       (y >  y1) .* (LE1 + tand(p.LambdaLEout) .* (y - y1));

    geom.TE_fun = @(y) (y <= y1) .* (TE0 + tand(p.LambdaTEin)  .* y) + ...
                       (y >  y1) .* (TE1 + tand(p.LambdaTEout) .* (y - y1));

    % Local chord
    geom.c_fun = @(y) geom.TE_fun(y) - geom.LE_fun(y);

    % Spars (as percentages of chord)
    geom.LE_spar_fun = @(y) geom.LE_fun(y) + p.perc_LE * geom.c_fun(y);
    geom.TE_spar_fun = @(y) geom.LE_fun(y) + p.perc_TE * geom.c_fun(y);
end

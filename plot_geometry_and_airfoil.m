clear; close all; clc;

%% ============================================================
% 1. LOAD CONSTANTS + DESIGN VECTOR + BOUNDS
% ============================================================
const = constants_ae4205();
[x0, lb, ub] = design_vector_init();

% Extract design vector components
Mcr = x0(1);
hcr = x0(2);
LambdaLEin_c = x0(3);
cr_in_c      = x0(4);
s_out_c      = x0(5);
perc_LE_c    = x0(6) / 100;
perc_TE_c    = x0(7) / 100;

% CST coefficients
Au     = x0(8:12)';
Al     = x0(13:17)';
Au_min = lb(8:12)';
Al_min = lb(13:17)';
Au_max = ub(8:12)';
Al_max = ub(13:17)';

% Inboard span
if isfield(const,'s_in')
    s_in_const = const.s_in;
else
    s_in_const = const.sin;
end

%% ============================================================
% 2. REFERENCE GEOMETRY
% ============================================================
cr_in_ref   = 5.996;
s_in_ref    = 4.939;
s_out_ref   = 14.040 - 4.939;

LambdaLEin_ref  = 26.3;
LambdaLEout_ref = 20.4;
dLambda_LE      = LambdaLEin_ref - LambdaLEout_ref;

perc_LE_ref = 17.5/100;
perc_TE_ref = 57.5/100;

%% ============================================================
% 3. CREATE WINGS FOR ALL CASES
% ============================================================
cases = struct([]);

cases(1).name  = 'Reference';
cases(1).wing  = Wing(Mcr, hcr, LambdaLEin_ref, cr_in_ref, s_in_ref, s_out_ref, ...
                      perc_LE_ref, perc_TE_ref, Au, Al, dLambda_LE);
cases(1).color = [0 0 0];
cases(1).style = '-';

cases(2).name  = 'Current';
cases(2).wing  = Wing(Mcr, hcr, LambdaLEin_c, cr_in_c, s_in_const, s_out_c, ...
                      perc_LE_c, perc_TE_c, Au, Al, dLambda_LE);
cases(2).color = [0 0.4470 0.7410];
cases(2).style = '-';

% Lower bound
cases(3).name  = 'Lower bound';
cases(3).wing  = Wing(Mcr, hcr, lb(3), lb(4), s_in_const, lb(5), ...
                      lb(6)/100, lb(7)/100, Au_min, Al_min, dLambda_LE);
cases(3).color = [0.8500 0.3250 0.0980];
cases(3).style = '--';

% Upper bound
cases(4).name  = 'Upper bound';
cases(4).wing  = Wing(Mcr, hcr, ub(3), ub(4), s_in_const, ub(5), ...
                      ub(6)/100, ub(7)/100, Au_max, Al_max, dLambda_LE);
cases(4).color = [0.4660 0.6740 0.1880];
cases(4).style = '--';

nCases = numel(cases);

%% ============================================================
% 4. COMPUTE GEOMETRY (SAFE PREALLOCATION)
% ============================================================
% First case defines struct layout
g1 = cases(1).wing.compute_planform();

% Preallocate struct array to avoid MATLAB assignment error
geoms = repmat(g1, nCases, 1);

% Fill
for k = 1:nCases
    geoms(k) = cases(k).wing.compute_planform();
end

max_y = max([geoms.y_total]);

all_x = [];
for k = 1:nCases
    g = geoms(k);
    all_x = [all_x g.LE0 g.LE1 g.LE2 g.TE0 g.TE1 g.TE2];
end

min_x = min(all_x);
max_x = max(all_x);

%% ============================================================
% 5. FIGURE: WING GEOMETRY + AIRFOIL
% ============================================================
figure('Position',[100 100 900 500]);

%% LEFT -- PLANFORM
subplot(1,2,1); hold on; grid on; axis equal;
xlabel('x [m]'); ylabel('y [m]');
title('Wing Planform: Reference, Design, Bounds');

% Dummy legend handles
href = plot(NaN,NaN,'-','Color',[0 0 0],'LineWidth',1.8);
hcur = plot(NaN,NaN,'-','Color',[0 0.4470 0.7410],'LineWidth',1.8);
hlow = plot(NaN,NaN,'--','Color',[0.8500 0.3250 0.0980],'LineWidth',1.8);
hup  = plot(NaN,NaN,'--','Color',[0.4660 0.6740 0.1880],'LineWidth',1.8);
htank= patch(NaN,NaN,'y','FaceAlpha',0.3,'EdgeColor','k');

for k = 1:nCases
    p = cases(k);
    g = geoms(k);

    y1 = p.wing.s_in;

    X_outline = [g.LE0 g.LE1 g.LE2 g.TE2 g.TE1 g.TE0 g.LE0];
    Y_outline = [0 y1 g.y_total g.y_total y1 0 0];

    plot(X_outline, Y_outline, ...
         'Color', p.color, 'LineStyle', p.style, 'LineWidth', 1.8);

    y_span = linspace(0, g.y_total, 200);
    plot(g.LE_spar_fun(y_span), y_span, ':', 'Color', p.color);
    plot(g.TE_spar_fun(y_span), y_span, ':', 'Color', p.color);

    % Fuel tank shading on current wing
    if strcmp(p.name,'Current')
        y_t1 = 0.85 * g.y_total;
        patch([g.LE_spar_fun(0) g.TE_spar_fun(0) g.TE_spar_fun(y_t1) g.LE_spar_fun(y_t1)], ...
              [0 0 y_t1 y_t1], p.color, 'FaceAlpha',0.05, 'EdgeColor',p.color);
    end
end

legend([href hcur hlow hup htank], ...
    {'Reference','Current','Lower bound','Upper bound','Fuel tank'}, ...
    'Location','northoutside');

xlim([min_x-0.5, max_x+0.5]);
ylim([0, max_y + 0.5]);

%% RIGHT -- AIRFOIL
subplot(1,2,2); hold on; grid on; axis equal;
title('CST Airfoil: Baseline + Bounds');
xlabel('x'); ylabel('z');

X = linspace(0,1,200)';

[Xtu,Xtl]     = D_airfoil2(Au,     Al,     X);
[Xtu_min,Xtl_min] = D_airfoil2(Au_min, Al_min, X);
[Xtu_max,Xtl_max] = D_airfoil2(Au_max, Al_max, X);

plot(Xtu(:,1),Xtu(:,2),'Color',[0 0.4470 0.7410],'LineWidth',1.4);
plot(Xtl(:,1),Xtl(:,2),'Color',[0 0.4470 0.7410],'LineWidth',1.4);

plot(Xtu_min(:,1),Xtu_min(:,2),'--','Color',[0.8500 0.3250 0.0980],'LineWidth',1.2);
plot(Xtl_min(:,1),Xtl_min(:,2),'--','Color',[0.8500 0.3250 0.0980],'LineWidth',1.2);

plot(Xtu_max(:,1),Xtu_max(:,2),'--','Color',[0.4660 0.6740 0.1880],'LineWidth',1.2);
plot(Xtl_max(:,1),Xtl_max(:,2),'--','Color',[0.4660 0.6740 0.1880],'LineWidth',1.2);

rae = load('rae2822.dat');
plot(rae(:,1),rae(:,2),'k','LineWidth',1.1);

legend('Baseline upper','Baseline lower', ...
       'Lower upper','Lower lower', ...
       'Upper upper','Upper lower', ...
       'RAE2822');

hold off;

function plot_wing_3D(W, Res, AC)

g = W.compute_planform();

Ny = 300;
y_geo = linspace(0, g.y_total, Ny);
LE = g.LE_fun(y_geo);
TE = g.TE_fun(y_geo);

%% Aerodynamic distributions
yq  = Res.Wing.Yst(:);
cl  = Res.Wing.cl(:);
cdi = Res.Wing.cdi(:);

rho = AC.Aero.rho;
V   = AC.Aero.V;
q   = 0.5 * rho * V^2;

c_fun = g.c_fun;
chord_q = c_fun(yq);

Lift_Npm = q .* chord_q .* cl;
Drag_Npm = q .* chord_q .* cdi;

L_scale = g.y_total * 0.12;
D_scale = max(TE-LE) * 0.20;

Lift_vis = Lift_Npm / max(Lift_Npm) * L_scale;
Drag_vis = Drag_Npm / max(Drag_Npm) * D_scale;

%% 3D figure
figure('Position',[100 100 1200 700]);
hold on; grid on; axis equal;

view(35,22);
title('3D Wing: Planform + Lift (blue) + Drag (red)');

xlabel('x [m]');
ylabel('y [m]');
zlabel('Lift ↑      Drag →');

%% Wing planform
fill3([LE TE fliplr(TE) fliplr(LE)], ...
      [y_geo y_geo fliplr(y_geo) fliplr(y_geo)], ...
      zeros(1,4*Ny), ...
      [0.90 0.90 0.90], ...
      'FaceAlpha',0.7, 'EdgeColor','k');

plot3(g.LE_spar_fun(y_geo), y_geo, zeros(size(y_geo)), 'k--');
plot3(g.TE_spar_fun(y_geo), y_geo, zeros(size(y_geo)), 'k--');

%% Arrays for curve endpoints
Lift_tip_x = zeros(size(yq));
Lift_tip_y = zeros(size(yq));
Lift_tip_z = zeros(size(yq));

Drag_tip_x = zeros(size(yq));
Drag_tip_y = zeros(size(yq));
Drag_tip_z = zeros(size(yq));

%% Lift (vertical bars)
for i = 1:length(yq)
    xc = g.LE_fun(yq(i)) + 0.25 * chord_q(i);

    Lift_tip_x(i) = xc;
    Lift_tip_y(i) = yq(i);
    Lift_tip_z(i) = Lift_vis(i);

    plot3([xc xc], [yq(i) yq(i)], [0 Lift_vis(i)], ...
          'Color',[0 0.3 1], 'LineWidth',2);
end

%% Drag (horizontal bars)
for i = 1:length(yq)
    xc = g.LE_fun(yq(i)) + 0.55 * chord_q(i);

    Drag_tip_x(i) = xc + Drag_vis(i);
    Drag_tip_y(i) = yq(i);
    Drag_tip_z(i) = 0;

    plot3([xc xc + Drag_vis(i)], [yq(i) yq(i)], [0 0], ...
          'Color',[1 0 0], 'LineWidth',2);
end

%% CONNECT the tips of lift and drag bars
plot3(Lift_tip_x, Lift_tip_y, Lift_tip_z, ...
      'Color',[0 0.3 1], 'LineWidth',2.5);

plot3(Drag_tip_x, Drag_tip_y, Drag_tip_z, ...
      'Color',[1 0 0], 'LineWidth',2.5);

%% Legend
h_plan = plot3(NaN,NaN,NaN,'k','LineWidth',1.5);
h_LE   = plot3(NaN,NaN,NaN,'k--','LineWidth',1.5);
h_TE   = plot3(NaN,NaN,NaN,'k--','LineWidth',1.5);
h_L    = plot3(NaN,NaN,NaN,'Color',[0 0.3 1],'LineWidth',3);
h_D    = plot3(NaN,NaN,NaN,'Color',[1 0 0],'LineWidth',3);

legend([h_plan h_LE h_TE h_L h_D], ...
       {'Planform','LE Spar','TE Spar','Lift','Drag'}, ...
       'Location','best');

%% Axis limits
xlim([min(LE)-1, max(TE)+D_scale+1]);
ylim([0, g.y_total+1]);
zlim([0, L_scale+1]);

hold off;
end

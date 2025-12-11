classdef Wing
    % WING  Parametric wing object holding planform and CST airfoil data.

    properties
        % --- DESIGN VARIABLES ---
        Mcr
        hcr

        LambdaLEin
        LambdaLEout
        LambdaTEin
        LambdaTEout

        cr_in
        s_in
        s_out

        perc_LE
        perc_TE

        Au   % CST upper coefficients
        Al   % CST lower coefficients
    end

    methods
        % ============================================================
        % Constructor
        % ============================================================
        function obj = Wing(Mcr, hcr, LambdaLEin, cr_in, s_in, s_out, ...
                            perc_LE, perc_TE, Au, Al, dLambda_LE)

            obj.Mcr = Mcr;
            obj.hcr = hcr;

            obj.LambdaLEin  = LambdaLEin;
            obj.LambdaLEout = LambdaLEin - dLambda_LE;

            % Assignment default trailing-edge sweep
            obj.LambdaTEin  = 4.2;
            obj.LambdaTEout = 4.2;

            obj.cr_in = cr_in;
            obj.s_in  = s_in;
            obj.s_out = s_out;

            obj.perc_LE = perc_LE;
            obj.perc_TE = perc_TE;

            obj.Au = Au(:);
            obj.Al = Al(:);
        end

        % ============================================================
        % Compute planform geometry
        % ============================================================
        function geom = compute_planform(obj)

            p = obj;

            y0 = 0;
            y1 = p.s_in;
            y2 = p.s_in + p.s_out;

            % Root LE / TE
            LE0 = 0;
            TE0 = p.cr_in;

            % Kink LE / TE
            LE1 = LE0 + tand(p.LambdaLEin)  * (y1 - y0);
            TE1 = TE0 + tand(p.LambdaTEin)  * (y1 - y0);

            % Tip LE / TE
            LE2 = LE1 + tand(p.LambdaLEout) * (y2 - y1);
            TE2 = TE1 + tand(p.LambdaTEout) * (y2 - y1);

            geom.s_in    = p.s_in;
            geom.s_out   = p.s_out;
            geom.y_total = y2;

            geom.LE0 = LE0; geom.LE1 = LE1; geom.LE2 = LE2;
            geom.TE0 = TE0; geom.TE1 = TE1; geom.TE2 = TE2;

            % Spanwise LE/TE functions
            geom.LE_fun = @(y) (y <= y1).* (LE0 + tand(p.LambdaLEin).*y) + ...
                               (y >  y1).* (LE1 + tand(p.LambdaLEout).*(y - y1));

            geom.TE_fun = @(y) (y <= y1).* (TE0 + tand(p.LambdaTEin).*y) + ...
                               (y >  y1).* (TE1 + tand(p.LambdaTEout).*(y - y1));

            % Local chord
            geom.c_fun = @(y) geom.TE_fun(y) - geom.LE_fun(y);

            % Spars
            geom.LE_spar_fun = @(y) geom.LE_fun(y) + p.perc_LE * geom.c_fun(y);
            geom.TE_spar_fun = @(y) geom.LE_fun(y) + p.perc_TE * geom.c_fun(y);
        end

        % ============================================================
        % AVL / Q3D geometry table: [x_LE  y_LE  z_LE  chord  twist]
        % ============================================================
        function AVL = make_AVL_geom(obj)

            geom = obj.compute_planform();

            y0 = 0;
            y1 = obj.s_in;
            y2 = obj.s_in + obj.s_out;

            LE0 = geom.LE_fun(y0);
            LE1 = geom.LE_fun(y1);
            LE2 = geom.LE_fun(y2);

            c0 = geom.c_fun(y0);
            c1 = geom.c_fun(y1);
            c2 = geom.c_fun(y2);

            twist0 = 0;  % [deg], set if you include twist distribution
            twist1 = 0;
            twist2 = 0;

            AVL = [ ...
                LE0   y0   0   c0   twist0;   % root
                LE1   y1   0   c1   twist1;   % kink
                LE2   y2   0   c2   twist2];  % tip
        end

        % ============================================================
        % Airfoil CST table for root and tip (2 sections)
        % ============================================================
        function [AF, eta] = make_airfoil_table_2(obj)
            % Returns:
            %   AF  = [Au1..Au5 Al1..Al5;
            %          Au1..Au5 Al1..Al5]
            %   eta = [0; 1]

            row = [obj.Au(:).' obj.Al(:).'];
            AF  = [row;
                   row];

            eta = [0; 1];
        end

        % (optional) CST → coordinates if you ever need it
        function [Xtu, Xtl] = airfoil(obj, X)
            [Xtu, Xtl] = D_airfoil2(obj.Au, obj.Al, X);
        end

    end
end
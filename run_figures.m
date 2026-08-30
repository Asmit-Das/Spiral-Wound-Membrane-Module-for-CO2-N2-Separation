function run_figures()
    close all;
    clc;

    fprintf('Generating figures...\n');

    fig2();
    fig3();
    fig8b();

    fprintf('Done.\n');
end

function fig2()
    TC = [57 67 77];
    pA = logspace(-1, 2, 400);

    colors = [
        0.15 0.30 0.75
        0.10 0.45 0.85
        0.35 0.60 0.95
    ];

    figure( ...
        'Color','w', ...
        'Position',[100 100 1000 450]);

    subplot(1,2,1);
    hold on;
    box on;

    for j = 1:numel(TC)

        P = ftmParams(TC(j));

        PA = permeanceCO2(pA, P);

        plot( ...
            pA, PA, ...
            '-', ...
            'LineWidth',2, ...
            'Color',colors(j,:));
        PAstar = permeanceCO2(P.pAst, P);

        plot( ...
            P.pAst, PAstar, ...
            'o', ...
            'MarkerSize',7, ...
            'MarkerFaceColor','w', ...
            'Color',colors(j,:), ...
            'LineWidth',1.5);
        text( ...
            0.12, ...
            permeanceCO2(0.12,P)*0.93, ...
            sprintf('%d ^\\circC',TC(j)), ...
            'Color',colors(j,:), ...
            'FontWeight','bold');
    end

    set(gca, ...
        'XScale','log', ...
        'YScale','log', ...
        'FontSize',11);

    xlabel('CO_2 partial pressure (kPa)');
    ylabel('CO_2 permeance (GPU)');

    xlim([0.1 100]);
    ylim([500 12000]);

    title('Fig 2a — Eq. 9 with Table 1');

    legend( ...
        {'Model','p_A^*'}, ...
        'Location','southwest', ...
        'Box','off');


   
    C = ftmConstants();

    Tsw = linspace(55,90,200);

    P = ftmParams(Tsw);

    invT = 1000 ./ (Tsw + C.T0);

    subplot(1,2,2);

    hold on;
    box on;

    % Left axis
    yyaxis left;

    semilogy( ...
        invT, ...
        P.PA0, ...
        '-', ...
        'LineWidth',2);

    semilogy( ...
        invT, ...
        P.PN2, ...
        '-', ...
        'LineWidth',2);

    ylabel('Parameter value');

    ylim([1 1e4]);


    % Right axis
    yyaxis right;

    semilogy( ...
        invT, ...
        P.pAst, ...
        '-', ...
        'LineWidth',2);

    semilogy( ...
        invT, ...
        P.eta, ...
        '-', ...
        'LineWidth',2);

    ylabel('Parameter value');

    ylim([0.1 100]);

    xlabel('1000/T (K^{-1})');

    title('Fig 2b — Arrhenius relationships');

    legend( ...
        {'P_{A,0}/\ell (GPU)', ...
         'P_{N_2}/\ell (GPU)', ...
         'p_A^* (kPa)', ...
         '\eta (-)'}, ...
        'Location','east', ...
        'Box','off');

    set(gca,'FontSize',11);
end


function fig3()

    C = ftmConstants();

    % Membrane parameters at 67 degC
    P = ftmParams(67);

    % Feed composition:
    % [CO2 N2 H2O O2]
    xF = [0.148 0.749 0.067 0.036];

    % Feed molar flow
    nF = 1.0;

    % Pressures [Pa]
    ph = 1.5 * C.atm;
    pl = 0.2 * C.atm;

    % Target CO2 recoveries
    recoveries = [0.30 0.70 0.90];

    % Colors
    colors = [
        0.15 0.30 0.75
        0.10 0.55 0.30
        0.80 0.25 0.15
    ];

    figure( ...
        'Color','w', ...
        'Position',[100 100 700 500]);

    Aguess = 200;

 
    for j = 1:numel(recoveries)

        target = recoveries(j);

        runModule = @(A) moduleCrossflow( ...
            nF, ...
            xF, ...
            ph, ...
            pl, ...
            A, ...
            P, ...
            struct('N',400));

        [Aguess,R] = solveArea( ...
            runModule, ...
            target, ...
            @(r) r.recovery, ...
            Aguess);

        z = R.A ./ R.A(end);

       
        yyaxis left;
        hold on;
        box on;

        plot( ...
            z, ...
            R.pA, ...
            '-', ...
            'LineWidth',2, ...
            'Color',colors(j,:));

        text( ...
            1.01, ...
            R.pA(end), ...
            sprintf('%d%%',round(100*target)), ...
            'Color',colors(j,:), ...
            'FontWeight','bold');


        %% CO2 permeance
        yyaxis right;
        hold on;

        plot( ...
            z, ...
            R.PCO2, ...
            '--', ...
            'LineWidth',2, ...
            'Color',colors(j,:));

    end

    %% Formatting
    yyaxis left;

    ylabel('CO_2 partial pressure (kPa)');
    ylim([0 50]);

    yyaxis right;

    ylabel('CO_2 permeance (GPU)');
    ylim([0 2500]);

    xlabel('Normalized distance from feed inlet');

    xlim([0 1.08]);

    title('Fig 3 — Profiles at 30 / 70 / 90% recovery');

    set(gca,'FontSize',11);
end



function fig8b()


    figure( ...
        'Color','w', ...
        'Position',[100 100 500 420]);

    hold on;
    box on;

    pA = logspace(0,2,300);

    TC = [57 67 77 87];

    colors = [
        0.05 0.15 0.45
        0.10 0.40 0.80
        0.40 0.65 0.95
        0.70 0.85 0.95
    ];

    for j = 1:numel(TC)

        P = ftmParams(TC(j));

        CO2permeance = permeanceCO2(pA,P);

        selectivity = CO2permeance ./ P.PN2;

        plot( ...
            pA, ...
            selectivity, ...
            '-', ...
            'LineWidth',2, ...
            'Color',colors(j,:));
    end

    set(gca, ...
        'XScale','log', ...
        'FontSize',11);

    xlabel('CO_2 partial pressure (kPa)');
    ylabel('CO_2/N_2 selectivity');

    xlim([1 100]);
    ylim([50 350]);

    title('Fig 8b — Selectivity');

    legend( ...
        {'57 ^\circC', ...
         '67 ^\circC', ...
         '77 ^\circC', ...
         '87 ^\circC'}, ...
        'Location','northeast', ...
        'Box','off');
end



function P = ftmParams(T_C,cal)

    if nargin < 2 || isempty(cal)
        cal = ftmCalibration('none');
    end

    C = ftmConstants();

    T_K = T_C + C.T0;

    RT = C.R .* T_K;

    P.T_C = T_C;
    P.T_K = T_K;

    % CO2 low-pressure permeance
    P.PA0 = ...
        cal.PA0 .* ...
        1.10e9 .* ...
        exp(-38.42 ./ RT);

    % Nonlinear parameter
    P.eta = ...
        cal.eta .* ...
        1.57e-2 .* ...
        exp(9.52 ./ RT);

    % Characteristic pressure
    P.pAst = ...
        cal.pAst .* ...
        2.66e5 .* ...
        exp(-29.30 ./ RT);

    % N2 permeance
    P.PN2 = ...
        cal.PN2 .* ...
        4.38e10 .* ...
        exp(-62.71 ./ RT);

    % Other species
    P.PH2O = P.PA0;

    P.PO2 = zeros(size(P.PN2));
end



function cal = ftmCalibration(mode)


    if nargin < 1 || isempty(mode)
        mode = 'none';
    end

   
    cal.PA0  = 1;
    cal.eta  = 1;
    cal.pAst = 1;
    cal.PN2  = 1;

    switch lower(mode)

        case 'none'
            % No calibration

        case 'paper67'

            % Raw parameters at 67 degC
            Praw = ftmParams(67);

            % Match p_A* = 7.5 kPa
            cal.pAst = 7.5 / Praw.pAst;

            % Match N2 permeance = 8.8 GPU
            cal.PN2 = 8.8 / Praw.PN2;

            % CO2 calibration
            pst = 7.5;

            f1 = sqrt(1 + pst/74.1) - 1;
            f2 = sqrt(1 + pst/3.9) - 1;

            M = [
                1 f1
                1 f2
            ];

            v = M \ [1464;1918];

            cal.PA0 = v(1) / Praw.PA0;

            cal.eta = ...
                (v(2)/v(1)) / Praw.eta;

        otherwise

            error( ...
                'ftmCalibration:mode', ...
                'Unknown calibration mode "%s".', ...
                mode);
    end
end



function C = ftmConstants()

    C.R = 8.314462618e-3;      % kPa*m^3/(mol*K)

    C.R_SI = 8.314462618;      % J/(mol*K)
    C.Vm_STP = 22414;

    % Pressure conversions
    C.atm  = 101325;            % Pa
    C.cmHg = 1333.2237;         % Pa
    C.psi  = 6894.757;          % Pa

    % GPU conversion
    C.GPU = ...
        1e-6 * ...
        (1/C.Vm_STP) * ...
        1e4 / C.cmHg;

    % Celsius -> Kelvin
    C.T0 = 273.15;
end



function PA = permeanceCO2(pA,P)

    pA = max(pA,eps);

    PA = P.PA0 .* ...
        (1 + ...
        P.eta .* ...
        (sqrt(1 + P.pAst ./ pA) - 1));
end



function Pi = permeanceVector(x,ph_kPa,P)

    pA = ph_kPa .* x(1);

    Pi = [
        permeanceCO2(pA,P)
        P.PN2
        P.PH2O
        P.PO2
    ];

    Pi = Pi(:).';
end



function [y,ok] = crossflowPermeate(x,Pi,ph,pl,k)


    if nargin < 5
        k = 1;
    end

    x  = x(:).';
    Pi = Pi(:).';

    n = numel(x);

    % Pressure ratio
    r = pl/ph;

    % Relative permeance
    a = Pi ./ Pi(k);

    % Active species
    active = ...
        (a > 0) & ...
        (x > 1e-14);

    if ~any(active)

        y = x;
        ok = false;

        return;
    end

    aa = a(active);
    xx = x(active);

    amin = min(aa);

    ucrit = r*(1-amin);

    % Nonlinear equation
    g = @(v) ...
        sum( ...
        aa .* xx ./ ...
        (r.*(aa-amin)+v)) - 1;

    % Lower bound
    vlo = max(x(k)-ucrit,1e-200);

    if g(vlo) <= 0

        y = zeros(1,n);

        y(active) = aa .* xx;

        total = sum(y);

        if total > 0
            y = y/total;
        else
            y = x;
        end

        ok = false;

        return;
    end

    % Upper bound
    vhi = max(10*vlo,1);

    guard = 0;

    while g(vhi) > 0 && guard < 400

        vhi = 10*vhi;

        guard = guard + 1;
    end

    if g(vhi) > 0

        y = x;
        ok = false;

        return;
    end

    % Bisection in log-space
    lo = log(vlo);
    hi = log(vhi);

    for it = 1:200

        mid = 0.5*(lo+hi);

        gm = g(exp(mid));

        if gm > 0
            lo = mid;
        else
            hi = mid;
        end

        if abs(hi-lo) < 1e-13
            break;
        end
    end

    v = exp(0.5*(lo+hi));

    % Permeate composition
    y = zeros(1,n);

    y(active) = ...
        aa .* xx ./ ...
        (r.*(aa-amin)+v);

    total = sum(y);

    if total > 0

        y = y/total;

        ok = true;

    else

        y = x;

        ok = false;
    end
end



function R = moduleCrossflow( ...
    nF,xF,ph_Pa,pl_Pa,A0,P,opt)


    if nargin < 7 || isempty(opt)
        opt = struct();
    end

    if ~isfield(opt,'N')
        opt.N = 600;
    end

    C = ftmConstants();

    N = opt.N;

    dA = A0/N;

    % Preallocate
    nh = zeros(N+1,1);

    X = zeros(N+1,4);

    Y = zeros(N+1,4);

    pA = zeros(N+1,1);

    PCO2 = zeros(N+1,1);

    % Membrane area
    A = (0:N)'*dA;


    z = [nF,xF(:).'];

    for j = 1:N+1

        [dz,yLocal,Pi] = ...
            deriv(z,ph_Pa,pl_Pa,P,C);

        nh(j) = z(1);

        X(j,:) = z(2:5);

        Y(j,:) = yLocal;

        % CO2 partial pressure
        pA(j) = ...
            ph_Pa*z(2)/1000;

        % CO2 permeance
        PCO2(j) = Pi(1);

        if j == N+1
            break;
        end

        z = rk4( ...
            z, ...
            dA, ...
            ph_Pa, ...
            pl_Pa, ...
            P, ...
            C);
    end

    % Results
    R.A = A;

    R.nh = nh;

    R.x = X;

    R.y_local = Y;

    R.pA = pA;

    R.PCO2 = PCO2;

    % Retentate
    R.nRet = nh(end);

    R.xRet = X(end,:);

    % Permeate
    R.nPerm = nF-R.nRet;

    if R.nPerm > 0

        R.yPerm = ...
            (nF*xF(:).' - ...
             R.nRet*R.xRet) / ...
            R.nPerm;

    else

        R.yPerm = zeros(1,4);
    end

    % Recovery
    R.recovery = ...
        1 - ...
        (R.nRet*R.xRet(1)) / ...
        (nF*xF(1));

    % Stage cut
    R.stagecut = R.nPerm/nF;
end



function [dz,yLocal,Pi] = deriv(z,ph,pl,P,C)

    nh = max(z(1),1e-12);

    % Composition
    x = max(z(2:5),0);

    sx = sum(x);

    if sx <= 0
        error( ...
            'deriv:InvalidComposition', ...
            'Gas composition became zero.');
    end

    x = x/sx;

    % Permeances [GPU]
    Pi_GPU = ...
        permeanceVector(x,ph/1000,P);

    Pi = Pi_GPU;

    % Convert GPU
    Pi_SI = Pi_GPU*C.GPU;

    % Local permeate composition
    [yLocal,ok] = ...
        crossflowPermeate(x,Pi_GPU,ph,pl);

    if ~ok
        % Do not stop calculation, but make the issue visible.
        warning( ...
            'deriv:PermeateSolver', ...
            'Permeate composition solver did not fully converge.');
    end

    % Partial-pressure driving force
    drivingForce = ...
        ph.*x - pl.*yLocal;

    % No negative flux
    drivingForce = max(drivingForce,0);

    % Species flux
    J = Pi_SI .* drivingForce;

    % Total flux
    Jtot = sum(J);

    % Differential equations
    dz = zeros(1,5);

    % Total molar flow
    dz(1) = -Jtot;

    % Component balances
    dz(2:5) = ...
        -(J-x*Jtot)/nh;
end



function zNext = rk4(z,h,ph,pl,P,C)


    k1 = deriv( ...
        z, ...
        ph,pl,P,C);

    k2 = deriv( ...
        z + 0.5*h*k1, ...
        ph,pl,P,C);

    k3 = deriv( ...
        z + 0.5*h*k2, ...
        ph,pl,P,C);

    k4 = deriv( ...
        z + h*k3, ...
        ph,pl,P,C);

    zNext = ...
        z + ...
        h/6 .* ...
        (k1+2*k2+2*k3+k4);

    % Numerical check
    if any(~isfinite(zNext))

        error( ...
            'rk4:NonFiniteState', ...
            'RK4 integration produced a non-finite state.');
    end

    % Prevent zero/negative total flow
    zNext(1) = ...
        max(zNext(1),1e-12);

    % Enforce physical composition
    zNext(2:5) = ...
        max(zNext(2:5),0);

    sx = sum(zNext(2:5));

    if sx <= 0

        error( ...
            'rk4:InvalidComposition', ...
            'Composition became invalid during integration.');
    end

    zNext(2:5) = ...
        zNext(2:5)/sx;
end



function [A0,R,info] = solveArea( ...
    runModule,target,getMetric,Aguess,opt)

    if nargin < 5 || isempty(opt)
        opt = struct();
    end

    if ~isfield(opt,'tol')
        opt.tol = 1e-7;
    end

    if ~isfield(opt,'maxit')
        opt.maxit = 60;
    end

    if Aguess <= 0

        error( ...
            'solveArea:InvalidGuess', ...
            'Initial area guess must be positive.');
    end

    % Work in log(area)
    residual = @(logA) ...
        getMetric(runModule(exp(logA)))-target;

    logA0 = log(Aguess);

    f0 = residual(logA0);

    % Already solved
    if abs(f0) < opt.tol

        A0 = Aguess;

        R = runModule(A0);

        info.iterations = 0;

        info.residual = f0;

        return;
    end

    % Initial bracket step
    step = 0.3;

    logA1 = logA0+step;

    f1 = residual(logA1);

    % Try opposite direction if needed
    if f0*f1 > 0 && abs(f1) > abs(f0)

        step = -step;

        logA1 = logA0+step;

        f1 = residual(logA1);
    end

    % Bracket search
    bracketIter = 0;

    while f0*f1 > 0

        logA0 = logA1;

        f0 = f1;

        step = 1.6*step;

        logA1 = logA1+step;

        f1 = residual(logA1);

        if ~isfinite(f1)

            error( ...
                'solveArea:NonFinite', ...
                'Module returned a non-finite metric.');
        end

        bracketIter = bracketIter+1;

        if bracketIter > 45

            error( ...
                'solveArea:Bracket', ...
                'Could not bracket the requested target.');
        end
    end

    % Bisection
    lo = logA0;
    hi = logA1;

    flo = f0;

    for it = 1:opt.maxit

        mid = 0.5*(lo+hi);

        fm = residual(mid);

        if ~isfinite(fm)

            error( ...
                'solveArea:NonFinite', ...
                'Non-finite residual encountered.');
        end

        if flo*fm <= 0

            hi = mid;

        else

            lo = mid;

            flo = fm;
        end

        if abs(fm) < opt.tol || ...
           abs(hi-lo) < 1e-9

            break;
        end
    end

    % Final area
    A0 = exp(0.5*(lo+hi));

    R = runModule(A0);

    info.iterations = ...
        it+bracketIter;

    info.residual = ...
        getMetric(R)-target;
end
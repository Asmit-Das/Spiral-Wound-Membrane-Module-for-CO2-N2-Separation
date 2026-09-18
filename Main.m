
T_C = 67;
T_K = T_C + 273.15;

P_feed = 1.5 * 101.325;       % kPa
P_perm = 0.2 * 101.325;       % kPa

xCO2_in = 0.15;


P0_CO2 =  arrhenius(1.1e9, 38.42, T_K);
eta    = arrhenius(1.57e-2, -9.52, T_K);
pstar  = arrhenius(2.66e5, 29.30, T_K);
P_N2   = arrhenius(4.38e10, 62.71, T_K);

GPU = 3.3464e-10;

P0_CO2_SI = P0_CO2 * GPU;
P_N2_SI   = P_N2   * GPU;

P_feed_Pa = P_feed * 1000;
P_perm_Pa = P_perm * 1000;

recovery = [0.30 0.70 0.90];
z = linspace(0,1,300);

pCO2_all = zeros(length(recovery),length(z));
perm_all = zeros(length(recovery),length(z));

for r = 1:length(recovery)

    target = recovery(r);

    % lambda = A_total / n_feed
    fun = @(lambda) recovery_error(lambda,target, P_feed_Pa,P_perm_Pa, P0_CO2_SI,eta,pstar,P_N2_SI, xCO2_in);


    % Find lambda

    lambda = fzero(fun,[0.001 100]);


    fprintf('Dimensionless area = %.6f\n',lambda);

    % Integrate feed-side equations
    options = odeset('RelTol',1e-8,'AbsTol',1e-10);


    [~,Y] = ode45(@(z,Y) feedODE(z,Y,lambda,P_feed_Pa,P_perm_Pa,P0_CO2_SI,eta,pstar,P_N2_SI),z,[1 xCO2_in],options);
    n = Y(:,1);
    xCO2 = Y(:,2);


    % CO2 partial pressure

    pCO2 = xCO2 * P_feed;

    % CO2 permeance

    pCO2_safe = max(pCO2,1e-12);

    permeance = P0_CO2 .* ( 1 + eta .* (sqrt(1 + pstar ./ pCO2_safe) - 1 ));

    pCO2_all(r,:) = pCO2;

    perm_all(r,:) = permeance;

    % Check recovery

    calculated_recovery = 1 - (n(end)*xCO2(end)) / (n(1)*xCO2(1));


    fprintf('Calculated recovery = %.4f %%\n',...
        calculated_recovery*100);

    fprintf('Outlet CO2 partial pressure = %.4f kPa\n',...
        pCO2(end));

end

% 9. FIGURE

figure('Color','w');

yyaxis left

plot(z,pCO2_all(1,:),...
    'LineWidth',2);

hold on;

plot(z,pCO2_all(2,:),...
    'LineWidth',2);

plot(z,pCO2_all(3,:),...
    'LineWidth',2);

ylabel('CO_2 partial pressure (kPa)');

ylim([0 50]);

yticks(0:10:50);


% RIGHT AXIS
yyaxis right

plot(z,perm_all(1,:),...
    'LineWidth',2);

plot(z,perm_all(2,:),...
    'LineWidth',2);

plot(z,perm_all(3,:),...
    'LineWidth',2);

ylabel('CO_2 permeance (GPU)');

ylim([0 2500]);

yticks(0:500:2500);


% X AXIS

xlabel('Normalized distance from feed inlet');

xlim([0 1]);

xticks(0:0.2:1);

grid on;
box on;

title('Effect of Carrier Saturation on Module Performance');


% LABELS

yyaxis left

text(0.84,pCO2_all(1,end)+1,...
    '30%','FontWeight','bold');

text(0.84,pCO2_all(2,end)+1,...
    '70%','FontWeight','bold');

text(0.84,pCO2_all(3,end)+1,...
    '90%','FontWeight','bold');

% FEED-SIDE DIFFERENTIAL EQUATIONS
function dYdz = feedODE( ...
                ~,Y,lambda,...
                P_feed,P_perm,...
                P0_CO2,eta,pstar,P_N2)

    % State variables

    n = Y(1);
    x = Y(2);


    % N2 composition

    xN2 = 1-x;

    % Permeate side approximation
    %
    % For the low-pressure/vacuum-side calculation:
    % CO2 permeate partial pressure is taken approximately zero.
    %
    % N2 sees the permeate pressure.

    yCO2 = 0;
    yN2  = 1;
    % Partial pressures

    pCO2_h = x * P_feed;
    pN2_h  = xN2 * P_feed;

    pCO2_l = yCO2 * P_perm;
    pN2_l  = yN2 * P_perm;


    % CO2 facilitated permeance

    pCO2_safe = max(pCO2_h,1e-10);

    PCO2 = P0_CO2 .* ...
        ( ...
        1 + eta .* ...
        ( ...
        sqrt(1 + pstar./(pCO2_safe/1000)) - 1 ...
        ) ...
        );


    % CO2 flux


    JCO2 = PCO2 .* ...
        (pCO2_h - pCO2_l);

    JN2 = P_N2 .* ...
        (pN2_h - pN2_l);


    % Total flux

    Jtotal = JCO2 + JN2;


    % Feed material balance
    %
    % dn/dz = -lambda * Jtotal
    dn_dz = -lambda * Jtotal;
    % CO2 mole-fraction balance
    %
    % dx/dz =
    % lambda*(x*Jtotal - JCO2)/n

    dx_dz = lambda * (x*Jtotal - JCO2) / max(n,1e-12);
    dYdz = [
        dn_dz
        dx_dz
    ];

end

% RECOVERY ERROR

function error_value = recovery_error(lambda,target, P_feed,P_perm,P0_CO2,eta,pstar,P_N2,xCO2_in)


    % Integrate

    options = odeset('RelTol',1e-7,...
        'AbsTol',1e-9,...
        'Events',@stopEvent);


    [~,Y] = ode45(@(z,Y) feedODE(z,Y,lambda,P_feed,P_perm,P0_CO2,eta,pstar,P_N2),[0 1],[1 xCO2_in],options);


    % Outlet

    n_out = Y(end,1);
    x_out = Y(end,2);


    % Recovery

    actual_recovery = 1 - (n_out*x_out)/xCO2_in;


    % Error

    error_value = actual_recovery - target;

end

% STOP INTEGRATION IF NUMERICALLY UNPHYSICAL
function [value,isterminal,direction] = stopEvent(~,Y)

    n = Y(1);
    x = Y(2);

    value = min(n-1e-8,x-1e-8);

    isterminal = 1;

    direction = -1;

end
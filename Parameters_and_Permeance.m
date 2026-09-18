%% Partial pressure
pAh = 0.1:0.1:100;       
n = length(pAh);

%% Temperature
Tc = 57:10:77; 
T=Tc + 273.15;
permeance = zeros(length(T), n);

PA0_l_values = zeros(size(T));
eta_values   = zeros(size(T));
pA_star_values = zeros(size(T));

%% Calculate permeance
for i = 1:length(T)

    PA0_l = arrhenius(1.1e9, 38.42, T(i));
    
    eta = arrhenius(1.57e-2, -9.52, T(i));
    
    pA_star = arrhenius( 2.66e5, 29.30, T(i));

    PA0_l_values(i) = PA0_l;
    eta_values(i) = eta;
    pA_star_values(i) = pA_star;

    % Permeance equation
    permeance(i,:) = PA0_l .*(1 + eta .* (sqrt(1 + pA_star ./ pAh) - 1));

end
figure;

semilogx(pAh, permeance, 'LineWidth', 1.5);

xlabel('Partial Pressure, p_A^h');
ylabel('Permeance, P_A/l');
title('Permeance vs Partial Pressure');

legend('T = 57 K', ...
    'T = 67 K', ...
    'T = 77 K', ...
    'Location', 'best');

grid on;
figure;

plot(1000./T, log(PA0_l_values), '-o', 'LineWidth', 1.5);
hold on;

plot(1000./T, log(abs(eta_values)), '-s', 'LineWidth', 1.5);

plot(1000./T, log(pA_star_values), '-^', 'LineWidth', 1.5);

xlabel('1000/T (K^{-1})');
ylabel('ln(Parameter)');

legend('ln(P_{A,0}/l)', 'ln(\eta)', 'ln(p_A^*)', ...
       'Location', 'best');

grid on;
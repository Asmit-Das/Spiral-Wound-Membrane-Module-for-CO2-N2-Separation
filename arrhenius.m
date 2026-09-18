function k = arrhenius(A, Ea, T)
% Arrhenius equation
% A  = Pre-exponential factor
% Ea = Activation energy (kJ/mol)
% R  = Gas constant (J/mol-K)
% T = Temperature in K
R = 8.314;   % J/(mol-K)

%syms T

k = A * exp(-Ea*1000/(R*T));
end
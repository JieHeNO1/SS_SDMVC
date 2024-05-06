function [M,dM] = OSLangevin(B,C)

% This function is a basic magnetization function with external B field
% (H*mu0) in mT and the C(Concentration) in mg/ml as inputs.

Kb = 1.38e-23; %J/k
T = 298.15; % k
Ms = 100e3; % A/m synomag
D = 50e-9;  % 70nm mean core diameter synomag

Vol = 1/6*pi*D^3;%m^3 
m = Vol*Ms;
B = 1e-3*B;  % mT→T

rho_Fe3O4 = 5.18;	% g/cm^3=g/mL
rho_Fe3O4 = rho_Fe3O4*1e9;	% mg/m^3
mass_Fe3O4 = rho_Fe3O4*Vol; % 单个粒子质量 mg
% C = 1;              % 1 mg/mL
C = C./mass_Fe3O4;   % particles/ml
beta = m/(Kb*T);
Xi = beta.*B;

M = C.*m.*(coth(Xi) - 1./(Xi));
M(Xi==0) = 0;
dM = beta.*C.*m.*(1./(Xi.^2)-1./(sinh(Xi).^2));
% dM(Xi==0) = 1/3.*beta.*C.*m;
end
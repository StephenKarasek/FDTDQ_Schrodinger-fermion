function [DATA] = MaxPotential(SPACE, TIME, SIM, PULSE, ITP, LASER, DATA, DEBUG)

%% Potential Energy Map

CONSTANTS;

% V_eff = zeros(SPACE.N, 1);
V_sys = PotentialCoulomb(SIM, SPACE, PULSE);
% V_eff = V_atom;

% Maximum 
V_max = (1+e0*LASER.E0)*max(V_sys);
V_min = (1+e0*LASER.E0)*min(V_sys);

if abs(V_min)>abs(V_max)*1e1; V_min = -abs(mean(V_sys))*5; end
if abs(V_max)>abs(V_min)*1e1; V_max = abs(mean(V_sys))*5; end

if V_max==0; V_max = 0.1*abs(mean(V_min)); end
if V_min==0; V_min = -0.1*abs(mean(V_max)); end

if (V_min==0)&&(V_max==0); V_max = abs(e0); V_min = e0; end

% if abs(V_min)>abs(V_max)*1e1; V_min = -abs(V_max)*10; end
% if abs(V_max)>abs(V_min)*1e1; V_max = abs(V_min)*10; end


DATA.V_.Vsys = V_sys;
DATA.V_.Vmax = V_max;
DATA.V_.Vmin = V_min;





end
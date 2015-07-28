%% CONSTANTS

%% Atomic Units
% Constants
c0 = 2.99792458*10^8;
mu0 = 4*pi*10^-7;
n0 = 1;
eps0 = 1e-9/(36*pi);
eta0 = sqrt(mu0/eps0);
Kb = 1.3806488*10^-23; %J / K


%% Fundamental Atomic Units (Hartree)
hPlanck = 1.0546*10^-34; %J * Sec
hPlanck_full = hPlanck*2*pi;%6.582119*10^-16; %eV * Sec
e0 = -1.60217*10^-19; % C
me = 9.10938291*10^-31; %kg
Ke = 1/(4*pi*eps0);
FineStruct = (e0^2)/(4*pi*eps0*hPlanck*c0); % UNITLESS
ca = 1/FineStruct;
re = Ke*((e0^2)/(me*ca^2));


%% Derived Atomic Units (Hartree)
a0 = (4*pi*eps0*hPlanck^2)/(me*e0^2);
Eh = (me*e0^4)/((4*pi*eps0*hPlanck)^2);
Ta = hPlanck/Eh;
Va = FineStruct*ca;
Fa = Eh/a0;
Tempa = Eh/Kb;
Presa = Eh/(a0^3);
Efa = Eh/(e0*a0);
Edma = e0*a0;

RydbergConst = (me*e0^4)/(8*eps0^2*(hPlanck*2*pi)^3*c0);

%%


% save(
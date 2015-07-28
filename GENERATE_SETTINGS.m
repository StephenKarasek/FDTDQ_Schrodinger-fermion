% function [] = GENERATE_SETTINGS()

clear all;
% clc;
% close all hidden; 



%% SETTINGS
%
% Determines what will be displayed during simulation, what plots will show
% up afterwards, etc...
%
% DEBUG Settings

%**************************************************************************
% GUI
%..........................................................................
% export2wsdlg(checkboxlabels,defaultvariablenames,itemstoexport,title)
%--------------------------------------------------------------------------


%**************************************************************************
% File Name length
%..........................................................................
% >>> If the length of the parameter file name is getting too long, since
% the description is included in the name itself, this shortens it to a
% 4-digit integer number, based on the number of param files within the
% default directory where the params are usually stored.
%--------------------------------------------------------------------------
DEBUG.DetailedFileName = true;
% DEBUG.DetailedFileName = false;
%--------------------------------------------------------------------------
% ...
%--------------------------------------------------------------------------




%% HARDWARE & SOFTWARE

%**************************************************************************
% Parallelization
%--------------------------------------------------------------------------
use_GUI = false;
% DEBUG.use_GUI = 


%**************************************************************************
% GUI
%--------------------------------------------------------------------------
% SETTINGS_gui;




%% VISUALIZATION

%**************************************************************************
% Animate Wavefunction
%--------------------------------------------------------------------------
% 'all', 'wavefunction', 'potential', 'none'

% DEBUG.VIDEO.animate = true;
DEBUG.VIDEO.animate = false;
%..........................................................................
DEBUG.VIDEO.animate_ = 'all';
% DEBUG.VIDEO.animate_ = 'wavefunction';
% DEBUG.VIDEO.animate_ = 'potential';



%**************************************************************************
% Record Movie
%--------------------------------------------------------------------------
% DEBUG.VIDEO.movie_ALL = true;
DEBUG.VIDEO.movie_ALL = false;
%..........................................................................
% DEBUG.VIDEO.movie_WaveFunction = true;
DEBUG.VIDEO.movie_WaveFunction = false;
%..........................................................................
% DEBUG.VIDEO.movie_Potential = true;
DEBUG.VIDEO.movie_Potential = false;



%**************************************************************************
% Data
%--------------------------------------------------------------------------
DEBUG.DATA.store = true;
% DEBUG.DATA.store = false;
%..........................................................................
DEBUG.DATA.save = true;
% DEBUG.DATA.save = false;
%..........................................................................
DEBUG.MEMORY.max = true;
% DEBUG.MEMORY.max = false;



%**************************************************************************
% Graphs
%--------------------------------------------------------------------------
DEBUG.GRAPHS.general = true;
% DEBUG.GRAPHS.general = false;
%..........................................................................
DEBUG.GRAPHS.CentroidPosition = true;
% DEBUG.GRAPHS.CentroidPosition = false;
%..........................................................................
DEBUG.GRAPHS.SpectralAnalysis = true;
% DEBUG.PLOT.SpectralAnalysis = false;
%..........................................................................
DEBUG.GRAPHS.MomentumSpace_VS_Laser = true;
% DEBUG.GRAPHS.MomentumSpace_VS_Laser = false;
%..........................................................................
DEBUG.GRAPHS.InstantaneousIonizationRate = true;
% DEBUG.GRAPHS.InstantaneousIonizationRate = false;
%..........................................................................
DEBUG.GRAPHS.EnergyOutput_VS_Laser = true;
% DEBUG.GRAPHS.EnergyOutput_VS_Laser = false;
%..........................................................................
DEBUG.GRAPHS.PhotonsAbsorbed_VS_Laser = true;
% DEBUG.GRAPHS.PhotonsAbsorbed_VS_Laser = false;
%..........................................................................
DEBUG.GRAPHS.Potential_VS_Laser = true;
% DEBUG.GRAPHS.Potential_VS_Laser = false;
%..........................................................................
DEBUG.GRAPHS.PopulationState = true;
% DEBUG.GRAPHS.PopulationState = false;
%..........................................................................
% DEBUG.GRAPHS.EnergyLevel = true;
DEBUG.GRAPHS.EnergyLevel = false;

%**************************************************************************
% ??? Settings
%
%--------------------------------------------------------------------------
OPTIONS.DisplayResolution_ = [640 480];







%% SIMULATION PARAMETERS

CONSTANTS;


%**************************************************************************
% Dimensionality
%--------------------------------------------------------------------------
% How arrays are calculated; 1D allows for a "cheat" that makes the work go
% a lot faster, but doing it the slow way is actually no slower than the
% higher dimensions. Translation- "slow" method scales fine; "fast"
% method is great but only works on 1D.
SPACE.Dims = 1;
% Dims = 2;
% Dims = 3;



%**************************************************************************
% Spatial Stepping
%--------------------------------------------------------------------------
% SPACE.N = 8192;
% SPACE.N = 4096;
% SPACE.N = 2048;
SPACE.N = 1024;
% SPACE.N = 512;

% SPACE.N = 256;

% SPACE.N = 500;

for n=1:SPACE.Dims; SPACE.Nd(n) = SPACE.N; end

% SPACE.length = a0*1.0e0;
% SPACE.length = a0*2.5e0;
% SPACE.length = a0*5.0e0;
SPACE.length = a0*1.0e1;
% SPACE.length = a0*1.5e1;
% SPACE.length = a0*2.0e1;
% SPACE.length = a0*2.5e1;
% SPACE.length = a0*3.0e1;
% SPACE.length = a0*3.5e1;
% SPACE.length = a0*4.0e1;
% SPACE.length = a0*4.5e1;
% SPACE.length = a0*5.0e1;

% SPACE.R = linspace(-SPACE.length/2,SPACE.length/2, SPACE.N);
% abs(SPACE.R(2)-SPACE.R(1));
SPACE.Axis = linspace(-SPACE.length/2,SPACE.length/2,SPACE.N);
SPACE.delta = abs(SPACE.Axis(2)-SPACE.Axis(1));
SPACE.Origin = 0;



%**************************************************************************
% Temporal Stepping
%--------------------------------------------------------------------------
% Critical Time-Step
TIME.critical = (2*me/hPlanck)*(SPACE.delta)^2;

% Time
% TIME.N = 50e6;
% TIME.N = 25e6; % TIME.save = 2.75e5;
% TIME.N = 22e6; % TIME.save = 2.75e4;
% TIME.N = 18.3e6;
% TIME.N = 13.6e6;

% TIME.N = 11e6;
% TIME.N = 8.16e6;
% TIME.N = 5.5e6;
% TIME.N = 4.4e6;
% TIME.N = 3.90e6;
% TIME.N = 2.0e6;

TIME.N = 1e6;
% TIME.N = 1e5;
% TIME.N = 1e4;

TIME.delta = TIME.critical*5e-2;
TIME.length = TIME.N*TIME.delta;
TIME.save = 5e4;
TIME.saveNum = length(TIME.save:TIME.save:TIME.N);

TIME.FramesPerSec = 10;
% TIME.animate = round(TIME.save*.1*TIME.FramesPerSec);
% TIME.animate = TIME.save;
TIME.animate = 5e3;

% 
% % Time
% t_crit = (2*me/hPlanck)*(delta_x)^2;
% T_len = 1e5;
% t_save = 25;%T_len/T_len;
% t_animate = t_save*100;
% % delta_t = 1*t_crit;
% delta_t = t_crit*1.5e-1;
% Time = linspace(0,T_len,t_save);

% time_dur = length(t_save:t_save:T_len);

% % ACTUAL UNITS
% ACTUAL.dx = delta_x*a0;
% ACTUAL.X_dist = X_dist*a0;
% ACTUAL.space = x*a0;
% %
% ACTUAL.dt = delta_t*Ta;
% ACTUAL.T_len = delta_t*Ta*T_len;
% ACTUAL.time = (1:T_len)*ACTUAL.dt;



%**************************************************************************
% Parameters to pass into Functions
%--------------------------------------------------------------------------
% % Space
% SPACETIME.Dims = Dims;
% SPACETIME.Space = [Xmin Xmax];
% SPACETIME.Nx = N;
% SPACETIME.dx = delta_x;
% 
% % Time
% SPACETIME.Time = [T_len];
% SPACETIME.dt = delta_t;
% 
% %**************************************************************************
% % Zero Padding
% %--------------------------------------------------------------------------
% % SPACE.ZeroPad_Factor = 1;
% SPACE.ZeroPad_Factor = 2;
% % SPACE.ZeroPad_Factor = 3;
% % SPACE.ZeroPad_Factor = 4;
% % SPACE.ZeroPad_Factor = 5;
% 
% % Extent of Added Zeroes
% SPACE.N0 = pow2(ceil(log2(SPACE.N*SPACE.ZeroPad_Factor)));
% % 2^(nextpow2(SPACE.N)+SPACE.ZeroPad_Factor);

% Amount to shave off...
% SPACE.SignalRange = [(-(SPACE.N0-SPACE.N)/2+1) ((SPACE.N0-SPACE.N)/2)]


%% Eigensolver: Imaginary Time Propagation

% Starts out FALSE, toggled if EigenFunction search is set to TRUE
ITP.Function = false; % ITP.Function = true;

% ITP.DEBUGVerbose = true;
ITP.DEBUGVerbose = false;

%**************************************************************************
% Initial Wavefunction
%..........................................................................
% ITP.Type = 'Sine';
% ITP.Type = 'Analytical';
ITP.Type = 'Gaussian';
% ITP.Type = 'Noise';
%--------------------------------------------------------------------------
ITP.InitialPos = 0;
ITP.Width = SPACE.length/10;
ITP.Momentum = 0;
%--------------------------------------------------------------------------
ITP.N_ = 3;
ITP.L_ = 0;
ITP.M_ = 0;


%--------------------------------------------------------------------------

%**************************************************************************
% USED ON 'QUICK' ITP ONLY
%..........................................................................


%__________________________________________________________________________
%:::[InfSq]:::    R=20a0 (dR=2048)
%--------------------------------------------------------------------------
% --> Gaussian
% ITP.ConvergenceRatio = 0;
% ITP.TimeDivisor = 0;

% --> Sine
% ITP.ConvergenceRatio = 0;
% ITP.TimeDivisor = 0;

% --> Noise
% ITP.ConvergenceRatio = 0;
% ITP.TimeDivisor = 0;
%--------------------------------------------------------------------------

%__________________________________________________________________________
%:::[InfSq]:::    R=40a0 (dR=2048)
% Gaussian Ratio.. ~= 
% Sine Ratio...... ~= 
% Noise Ratio..... ~= 
%--------------------------------------------------------------------------
% --> Gaussian
% ITP.ConvergenceRatio = 1e-8;
% ITP.TimeDivisor = 1e0*ITP.ConvergenceRatio;
% ITP.TimeDivisor = 1e-1;

% --> Sine
% ITP.ConvergenceRatio = 1e-7;
% ITP.TimeDivisor = 1e-5;

% --> Noise
% ITP.ConvergenceRatio = 4e-6;
% ITP.TimeDivisor = 2.3e-6;
%--------------------------------------------------------------------------

% ITP.ConvergenceRatio = 5e-5;
% ITP.TimeDivisor = 2e-1*ITP.ConvergenceRatio;
% ITP.TimeDivisor = 2.5e-5;
% ITP.Ratio = ITP.ConvergenceRatio/ITP.TimeDivisor;


%__________________________________________________________________________
%:::[H1]:::    R=40a0 (dR=2048)
%--------------------------------------------------------------------------
% --> Gauss
ITP.ConvergenceRatio = 1e-4;
ITP.TimeDivisor = 5e-3;

% --> Sine
% ITP.ConvergenceRatio = 1e-4;
% ITP.TimeDivisor = 5e-3;

% --> Noise 
% ITP.ConvergenceRatio = 6e-3;
% ITP.TimeDivisor = 2e-3;
%--------------------------------------------------------------------------


%**************************************************************************
% Order
%..........................................................................
ITP.NthOrder = 2;
% ITP.NthOrder = 4;
% ITP.NthOrder = 6;
% ITP.NthOrder = 8;
% ITP.NthOrder = 10;
% ITP.NthOrder = 12;
%--------------------------------------------------------------------------


%**************************************************************************
% States to Find
%..........................................................................
% ITP.EnergyState = [1:6];
ITP.EnergyState = [1:3];
ITP.NumStates = length(ITP.EnergyState);
ITP.FindStates = ITP.NumStates;
ITP.AccurateStates = ITP.FindStates-1;
%--------------------------------------------------------------------------
for n=1:ITP.FindStates
%     ITP.En_(n) = (n*pi)^2/(2*(SPACE.length/a0)^2);
%     (n^2*hPlanck^2*pi)/(2*me*SPACE.length^2)
    ITP.OrbitalNum(n,:) = [n 0 0];
    ITP.OrbitalType(n) = 's';
%     switch ITP.OrbitalNum(n,1)
%         case 
end
%--------------------------------------------------------------------------



%% General Settings
%**************************************************************************
% System Potential
%--------------------------------------------------------------------------
% SIM.PotentialMap = 'H1';
% SIM.PotentialMap = 'H2+';
% SIM.PotentialMap = 'He+';
% SIM.PotentialMap = 'H2O';
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
SIM.PotentialMap = 'Square';
% SIM.PotentialMap = 'QHO';
% SIM.PotentialMap = 'HarmonicOscillator';
% SIM.PotentialMap = 'Vee';
% SIM.PotentialMap = 'Quadratic';
% SIM.PotentialMap = 'Quartic';
%..........................................................................


%**************************************************************************
% (Regularized Coulomb / Soft-Core) Potential
%--------------------------------------------------------------------------
% PULSE.Creg = (1e-3)*a0;
% PULSE.Creg = (5e-3)*a0;
% PULSE.Creg = (1e-2)*a0;
% PULSE.Creg = (5e-2)*a0;
% PULSE.Creg = (1e-1)*a0;
% PULSE.Creg = (5e-1)*a0;
PULSE.Creg = (1e-0)*a0;
% PULSE.Creg = (5e-0)*a0;




%**************************************************************************
% Finite Difference Method
%--------------------------------------------------------------------------
SIM.FDM.Type = 'FDTDQ';
% SIM.FDM.Type = 'GFDTDQ';
% SIM.FDM.Type = 'CrankNicolson';
% SIM.FDM.Type = 'CrankNicolson_NthOrder';
% SIM.FDM.Type = 'RungeKutta';
% SIM.FDM.Type = 'RungeKutta_NthOrder';
% SIM.FDM.Type = 'Numerov';
%..........................................................................


%**************************************************************************
% Boundary Conditions
%--------------------------------------------------------------------------
% SIM.BC.Type = 'Inf';
SIM.BC.Type = 'Dirichlet';
% SIM.BC.Type = 'D_Mask'; %SIM.BC.Type = 'Dirichlet_Mask';
% SIM.BC.Type = 'Neumann';
% SIM.BC.Type = 'Robin';

% SIM.BC.Type = 'Convolution';
% SIM.BC.Type = 'Exact'; % 'ExactStationary', 'ExactTransient';

% SIM.BC.Type = 'ExteriorComplexScaling';
% SIM.BC.Type = 'ComplexAbsorbing';
% SIM.BC.Type = 'RadiationReflection';

% SIM.BC.Type = 'Artificial';
% SIM.BC.Type = 'Transparent';
% SIM.BC.Type = 'DiscreteTransport';
%..........................................................................


%**************************************************************************
% Gauge Choice
%--------------------------------------------------------------------------
SIM.Gauge.Type = 'Coulomb';
% SIM.Gauge.Type = 'Lorenz';

%**************************************************************************
% Gauge Form
%--------------------------------------------------------------------------
SIM.Gauge.Form = 'Length';
% SIM.Gauge.Form = 'Velocity';




%% Finite Difference Scheme
switch SIM.FDM.Type
%**************************************************************************
% FDTQ-Q
    case {'FDTDQ','fdtdq'}
%--------------------------------------------------------------------------
SIM.FDM.Psi_Rorder = 2;
SIM.FDM.Psi_Torder = 1;
SIM.FDM.Veff_Rorder = SIM.FDM.Psi_Rorder;
SIM.FDM.Veff_Torder = SIM.FDM.Psi_Torder;


%**************************************************************************
% Masking Function
    case {'GFDTDQ','gfdtdq'}
%--------------------------------------------------------------------------
SIM.FDM.Psi_Rorder = 2;
SIM.FDM.Psi_Torder = 1;
SIM.FDM.Veff_Rorder = SIM.FDM.Psi_Rorder;
SIM.FDM.Veff_Torder = SIM.FDM.Psi_Torder;


%**************************************************************************
% DEFAULT- Infinite Well
    otherwise
%--------------------------------------------------------------------------
% SIM.FDM.Psi_Rorder = 2;
% SIM.FDM.Psi_Torder = 1;
% SIM.FDM.Veff_Rorder = SIM.FDM.Psi_Rorder;
% SIM.FDM.Veff_Torder = SIM.FDM.Psi_Torder;
error('Finite Difference Scheme undefined!!!')



end



%% Boundary Conditions

switch SIM.BC.Type
%**************************************************************************
% Infinite Boundaries
    case {'Inf'; 'Dirichlet'}
%--------------------------------------------------------------------------
% SIM.BC.Type = 'Inf';
SIM.BC.N = 1;
% SPACE.N_ = SPACE.N-2.*BCs.N;



%**************************************************************************
% Masking Function
    case {'Mask';'Dirichlet_Mask'}
%--------------------------------------------------------------------------
% SIM.BC.Type = 'Mask';
SIM.BC.params_.Mask_param = 1;
% SIM.BC.params_.Mask_param = 0.05;

% Absorbing Region
% SIM.BC.N = round(SPACE.N/128);
% SIM.BC.N = round(SPACE.N/128);

% SIM.BC.N = ceil(2^((log2(SPACE.N)/2)));
SIM.BC.N = ceil(2^((log2(SPACE.N)/3)));


%**************************************************************************
% Exterior Complex Scaling
    case 'ECS'
%--------------------------------------------------------------------------
% SIM.BC.Type = 'ECS';

% ...




%**************************************************************************
% DEFAULT- Infinite Well
    otherwise
%--------------------------------------------------------------------------
% SIM.BC.Type = 'Inf';
SIM.BC.N = 1;
% SPACE.N_ = SPACE.N-2.*BCs.N;




end

SIM.BC.Nratio = (SIM.BC.N/SPACE.N);


%**************************************************************************
% New Spatial Expanse
%--------------------------------------------------------------------------
SPACE.N_ = SPACE.N-2.*SIM.BC.N;



% BCs.length = (SPACE.N/2)-SIM.BC.N;
for n=1:SPACE.Dims
    SIM.BC.R{n,1} = 1:SIM.BC.N;
    SIM.BC.R{n,2} = SPACE.N_-SIM.BC.N+1:SPACE.N_;
end

SPACE.length_ = SPACE.length + 2*SIM.BC.N*SPACE.delta;
SPACE.Axis_ = linspace(-SPACE.length_/2,SPACE.length_/2,SPACE.N);


%**************************************************************************
% Zero Padding
%--------------------------------------------------------------------------
% SPACE.ZeroPad_Factor = 0;
% SPACE.ZeroPad_Factor = 0.5;
% SPACE.ZeroPad_Factor = 1;
% SPACE.ZeroPad_Factor = 1.5;
SPACE.ZeroPad_Factor = 2;
% SPACE.ZeroPad_Factor = 2.5;
% SPACE.ZeroPad_Factor = 3;
% SPACE.ZeroPad_Factor = 3.5;
% SPACE.ZeroPad_Factor = 4;
% SPACE.ZeroPad_Factor = 4.5;
% SPACE.ZeroPad_Factor = 5;
%..........................................................................

% Extent of Added Zeroes
SPACE.N0 = pow2(ceil(log2(SPACE.N)));
% SPACE.N0 = pow2(ceil(log2(SPACE.N_*SPACE.ZeroPad_Factor)));
% 2^(nextpow2(SPACE.N)+SPACE.ZeroPad_Factor);
SPACE.ZeroPad_Ratio = SPACE.N0/SPACE.N;



%% Initial State
%**************************************************************************
% Pulse Parameters
%..........................................................................
% PULSE.Type = 'Gaussian';
% PULSE.Type = 'Fourier-Bessel'; 

PULSE.Type = 'Sine';
% PULSE.Type = 'Dirac';
% PULSE.Type = 'EigenState'; ITP.Function = true;
% PULSE.Type = 'Analytical'; 
%--------------------------------------------------------------------------
PULSE.InitialPos = 0;
% PULSE.Width = a0*(5e01);

% PULSE.Width = a0*(1e-1); PULSE.EnergyState = [1:5];
PULSE.Width = a0*(1e-0); PULSE.EnergyState = [1:5];
% PULSE.Width = a0*(1e1); PULSE.EnergyState = [1:5];

% PULSE.Width = a0*(5e-0);
% PULSE.Width = a0*(1e-0);
% PULSE.Width = a0*(5e-1);
% PULSE.Width = a0*(1e-1);
% PULSE.Width = a0*(5e-2);
% PULSE.Width = a0*(1e-2);

PULSE.Momentum = 0;
% PULSE.Momentum = -1e-3*hPlanck/(me*a0);
%--------------------------------------------------------------------------
PULSE.NumStates = 1; PULSE.EnergyState = 1;
% PULSE.NumStates = 1; PULSE.EnergyState = 2;
% PULSE.NumStates = 1; PULSE.EnergyState = 3;
% PULSE.NumStates = 2; PULSE.EnergyState = [1 2];
% PULSE.NumStates = 2; PULSE.EnergyState = [1 3];
% PULSE.NumStates = 2; PULSE.EnergyState = [2 3];

% PULSE.ViewStates = [1 2];
PULSE.ViewStates = [1 2 3];
% PULSE.ViewStates = [1 2 3 4];
% PULSE.ViewStates = [1 2 3 4 5];
PULSE.NumViewStates = length(PULSE.ViewStates);

% PULSE.CalcStates = [1:(PULSE.ViewStates(end)+0)];
PULSE.CalcStates = [1:(PULSE.ViewStates(end)+1)];
% PULSE.CalcStates = [1:(PULSE.ViewStates(end)+2)];
% PULSE.CalcStates = [1:(PULSE.ViewStates(end)+3)];
% PULSE.CalcStates = [1:6];
% PULSE.CalcStates = [1:5];
% PULSE.CalcStates = [1:4];
% PULSE.CalcStates = [1:3];
PULSE.NumCalcStates = length(PULSE.CalcStates);


PULSE.SpreadStates = (2*(PULSE.NumCalcStates*(SPACE.ZeroPad_Ratio+1))-1);

% PULSE.SpreadStates = (2*SPACE.ZeroPad_Factor+1);
%..........................................................................





%% Laser Source

%**************************************************************************
% Gauge Theory (Laser Coupling)
%--------------------------------------------------------------------------
LASER.GaugeCoupling = 'coulomb';
% LASER.GaugeCoupling = 'lorenz';
%..........................................................................


%**************************************************************************
% Laser Source
%--------------------------------------------------------------------------
% LASER.Type = 'pulse';
LASER.Type = 'CW';
% LASER.Type = 'off';
%..........................................................................


%**************************************************************************
% Carrier-Wave
%--------------------------------------------------------------------------
% LASER.l0 = 800e-9;
% LASER.l0 = 1.9541e-9;
% LASER.l0 = 100e-9;
% LASER.l0 = 121e-9;


%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% 2a0 WELL
%..........................................................................
% 1 Photon
% LASER.l0 = 2.3727e-8;

% 2 Photon
% LASER.l0 = 5.7553e-8;

% 3 Photon
% LASER.l0 = 8.1230e-8;

% 4 Photon
% LASER.l0 = 2.5106e-9;

% 4 Photon
% LASER.l0 = 2.5106e-9;

% 4 Photon
% LASER.l0 = 2.5106e-9;


%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% 5a0 WELL
%..........................................................................
% 1 Photon
% LASER.l0 = 2.3727e-8;

% 2 Photon
% LASER.l0 = 5.7553e-8;

% 3 Photon
% LASER.l0 = 8.1230e-8;

% 4 Photon
% LASER.l0 = 2.5106e-9;

% 4 Photon
% LASER.l0 = 2.5106e-9;

% 4 Photon
% LASER.l0 = 2.5106e-9;



%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% 10a0 WELL
%..........................................................................
% LASER.l0 = 800e-9; LASER.nPhoton = 1; LASER.I0 = 0e0;
% LASER.l0 = 800e-9; LASER.nPhoton = 1; LASER.I0 = 3.5e10;
% LASER.l0 = 800e-9; LASER.nPhoton = 1; LASER.I0 = 3.5e14;
% LASER.l0 = 100e-9; LASER.nPhoton = 1; LASER.I0 = 3.5e12;

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% E1 -> E2
% LASER.l0 = 306.95e-9; LASER.nPhoton = 1; LASER.I0 = 5e10;
% LASER.l0 = 306.95e-9; LASER.nPhoton = 1; LASER.I0 = 2.5e11;
% LASER.l0 = 306.95e-9; LASER.nPhoton = 1; LASER.I0 = 1e12;
LASER.l0 = 306.95e-9; LASER.nPhoton = 1; LASER.I0 = 3.33e12;
% LASER.l0 = 306.95e-9; LASER.nPhoton = 1; LASER.I0 = 1e13;
% LASER.l0 = 306.95e-9; LASER.nPhoton = 1; LASER.I0 = 2.5e13;
% LASER.l0 = 306.95e-9; LASER.nPhoton = 1; LASER.I0 = 5.5e13;
% LASER.l0 = 306.95e-9; LASER.nPhoton = 1; LASER.I0 = 1e14;
% LASER.l0 = 306.95e-9; LASER.nPhoton = 1; LASER.I0 = 3.33e14;


% LASER.l0 = 613.90e-9; LASER.nPhoton = 1; LASER.I0 = 5e10;
% LASER.l0 = 613.90e-9; LASER.nPhoton = 1; LASER.I0 = 2.5e11;
% LASER.l0 = 613.90e-9; LASER.nPhoton = 1; LASER.I0 = 1e12;
% LASER.l0 = 613.90e-9; LASER.nPhoton = 1; LASER.I0 = 3.33e12;
% LASER.l0 = 613.90e-9; LASER.nPhoton = 1; LASER.I0 = 1e13;
% LASER.l0 = 613.90e-9; LASER.nPhoton = 1; LASER.I0 = 2.5e13;
% LASER.l0 = 613.90e-9; LASER.nPhoton = 1; LASER.I0 = 5.5e13;
% LASER.l0 = 613.90e-9; LASER.nPhoton = 1; LASER.I0 = 1e14;
% LASER.l0 = 613.90e-9; LASER.nPhoton = 1; LASER.I0 = 3.33e14;


% LASER.l0 = 460.042e-9; LASER.nPhoton = 1; LASER.I0 = 5e10;
% LASER.l0 = 460.042e-9; LASER.nPhoton = 1; LASER.I0 = 2.5e11;
% LASER.l0 = 460.042e-9; LASER.nPhoton = 1; LASER.I0 = 1e12;
% LASER.l0 = 460.042e-9; LASER.nPhoton = 1; LASER.I0 = 3.33e12;
% LASER.l0 = 460.042e-9; LASER.nPhoton = 1; LASER.I0 = 1e13;
% LASER.l0 = 460.042e-9; LASER.nPhoton = 1; LASER.I0 = 2.5e13;
% LASER.l0 = 460.042e-9; LASER.nPhoton = 1; LASER.I0 = 5.5e13;
% LASER.l0 = 460.042e-9; LASER.nPhoton = 1; LASER.I0 = 1e14;
% LASER.l0 = 460.042e-9; LASER.nPhoton = 1; LASER.I0 = 3.33e14;


% LASER.l0 = 153.4752e-9; LASER.nPhoton = 1; LASER.I0 = 5e10;
% LASER.l0 = 153.4752e-9; LASER.nPhoton = 1; LASER.I0 = 2.5e11;
% LASER.l0 = 153.4752e-9; LASER.nPhoton = 1; LASER.I0 = 1e12;
% LASER.l0 = 153.4752e-9; LASER.nPhoton = 1; LASER.I0 = 3.33e12;
% LASER.l0 = 153.4752e-9; LASER.nPhoton = 1; LASER.I0 = 1e13;
% LASER.l0 = 153.4752e-9; LASER.nPhoton = 1; LASER.I0 = 2.5e13;
% LASER.l0 = 153.4752e-9; LASER.nPhoton = 1; LASER.I0 = 5.5e13;
% LASER.l0 = 153.4752e-9; LASER.nPhoton = 1; LASER.I0 = 1e14;
% LASER.l0 = 153.4752e-9; LASER.nPhoton = 1; LASER.I0 = 3.33e14;





% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% E2 -> E3    *******************
% LASER.l0 = 184.17e-9; LASER.nPhoton = 1; LASER.I0 = (4^0)*0.7175e12;
% LASER.l0 = 184.17e-9; LASER.nPhoton = 1; LASER.I0 = (4^1)*0.7175e12;
% LASER.l0 = 184.17e-9; LASER.nPhoton = 1; LASER.I0 = (4^2)*0.7175e12;
% LASER.l0 = 184.17e-9; LASER.nPhoton = 1; LASER.I0 = (4^3)*0.7175e12;
% LASER.l0 = 184.17e-9; LASER.nPhoton = 1; LASER.I0 = (4^4)*0.7175e12;
% LASER.l0 = 184.17e-9; LASER.nPhoton = 1; LASER.I0 = (4^5)*0.7175e12;


% LASER.l0 = 368.34e-9; LASER.nPhoton = 2; LASER.I0 = (4^0)*2.8698e12;
% LASER.l0 = 368.34e-9; LASER.nPhoton = 2; LASER.I0 = (4^1)*2.8698e12;
% LASER.l0 = 368.34e-9; LASER.nPhoton = 2; LASER.I0 = (4^2)*2.8698e12;
% LASER.l0 = 368.34e-9; LASER.nPhoton = 2; LASER.I0 = (4^3)*2.8698e12;
% LASER.l0 = 368.34e-9; LASER.nPhoton = 2; LASER.I0 = (4^4)*2.8698e12;
% LASER.l0 = 368.34e-9; LASER.nPhoton = 2; LASER.I0 = (4^5)*2.8698e12;


% LASER.l0 = 276.25e-9; LASER.nPhoton = 1; LASER.I0 = (4^0)*0.7175e12;
% LASER.l0 = 276.25e-9; LASER.nPhoton = 1; LASER.I0 = (4^1)*0.7175e12;
% LASER.l0 = 276.25e-9; LASER.nPhoton = 1; LASER.I0 = (4^2)*0.7175e12;
% LASER.l0 = 276.25e-9; LASER.nPhoton = 1; LASER.I0 = (4^3)*0.7175e12;
% LASER.l0 = 276.25e-9; LASER.nPhoton = 1; LASER.I0 = (4^4)*0.7175e12;
% LASER.l0 = 276.25e-9; LASER.nPhoton = 1; LASER.I0 = (4^5)*0.7175e12;


% LASER.l0 = 276.25e-9; LASER.nPhoton = 2; LASER.I0 = (4^0)*2.8698e12;
% LASER.l0 = 276.25e-9; LASER.nPhoton = 2; LASER.I0 = (4^1)*2.8698e12;
% LASER.l0 = 276.25e-9; LASER.nPhoton = 2; LASER.I0 = (4^2)*2.8698e12;
% LASER.l0 = 276.25e-9; LASER.nPhoton = 2; LASER.I0 = (4^3)*2.8698e12;
% LASER.l0 = 276.25e-9; LASER.nPhoton = 2; LASER.I0 = (4^4)*2.8698e12;
% LASER.l0 = 276.25e-9; LASER.nPhoton = 2; LASER.I0 = (4^5)*2.8698e12;

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% 20a0 WELL
..........................................................................
% LASER.l0 = 800e-9; LASER.nPhoton = 1; LASER.I0 = 0e0;

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% E1 -> E2    *******************
%..........................................................................
% LASER.l0 = 1227.8e-9; LASER.nPhoton = 1; LASER.I0 = 1e9;
% LASER.l0 = 1227.8e-9; LASER.nPhoton = 1; LASER.I0 = 5e10;
% LASER.l0 = 1227.8e-9; LASER.nPhoton = 1; LASER.I0 = 2.5e11;
% LASER.l0 = 1227.8e-9; LASER.nPhoton = 1; LASER.I0 = 1e12;
% LASER.l0 = 1227.8e-9; LASER.nPhoton = 1; LASER.I0 = 3.33e12;
% LASER.l0 = 1227.8e-9; LASER.nPhoton = 1; LASER.I0 = 1e13;
% LASER.l0 = 1227.8e-9; LASER.nPhoton = 1; LASER.I0 = 2.5e13;
% LASER.l0 = 1227.8e-9; LASER.nPhoton = 1; LASER.I0 = 5.5e13;
% LASER.l0 = 1227.8e-9; LASER.nPhoton = 1; LASER.I0 = 1e14;
% LASER.l0 = 1227.8e-9; LASER.nPhoton = 1; LASER.I0 = 3.33e14;
%..........................................................................
% LASER.l0 = 1841.7e-9; LASER.nPhoton = 1; LASER.I0 = 5e10;
% LASER.l0 = 1841.7e-9; LASER.nPhoton = 1; LASER.I0 = 2.5e11;
% LASER.l0 = 1841.7e-9; LASER.nPhoton = 1; LASER.I0 = 1e12;
% LASER.l0 = 1841.7e-9; LASER.nPhoton = 1; LASER.I0 = 3.33e12;
% LASER.l0 = 1841.7e-9; LASER.nPhoton = 1; LASER.I0 = 1e13;
% LASER.l0 = 1841.7e-9; LASER.nPhoton = 1; LASER.I0 = 2.5e13;
% LASER.l0 = 1841.7e-9; LASER.nPhoton = 1; LASER.I0 = 5.5e13;
% LASER.l0 = 1841.7e-9; LASER.nPhoton = 1; LASER.I0 = 1e14;
% LASER.l0 = 1841.7e-9; LASER.nPhoton = 1; LASER.I0 = 3.33e14;
%..........................................................................
% LASER.l0 = 2455.6e-9; LASER.nPhoton = 2; LASER.I0 = 5e10;
% LASER.l0 = 2455.6e-9; LASER.nPhoton = 2; LASER.I0 = 2.5e11;
% LASER.l0 = 2455.6e-9; LASER.nPhoton = 2; LASER.I0 = 1e12;
% LASER.l0 = 2455.6e-9; LASER.nPhoton = 2; LASER.I0 = 3.33e12;
% LASER.l0 = 2455.6e-9; LASER.nPhoton = 2; LASER.I0 = 1e13;
% LASER.l0 = 2455.6e-9; LASER.nPhoton = 2; LASER.I0 = 2.5e13;
% LASER.l0 = 2455.6e-9; LASER.nPhoton = 2; LASER.I0 = 5.5e13;
% LASER.l0 = 2455.6e-9; LASER.nPhoton = 2; LASER.I0 = 1e14;
% LASER.l0 = 2455.6e-9; LASER.nPhoton = 2; LASER.I0 = 3.33e14;
%..........................................................................
% LASER.l0 = 613.90e-9; LASER.nPhoton = 1; LASER.I0 = 5e10;
% LASER.l0 = 613.90e-9; LASER.nPhoton = 1; LASER.I0 = 2.5e11;
% LASER.l0 = 613.90e-9; LASER.nPhoton = 1; LASER.I0 = 1e12;
% LASER.l0 = 613.90e-9; LASER.nPhoton = 1; LASER.I0 = 3.33e12;
% LASER.l0 = 613.90e-9; LASER.nPhoton = 1; LASER.I0 = 1e13;
% LASER.l0 = 613.90e-9; LASER.nPhoton = 1; LASER.I0 = 2.5e13;
% LASER.l0 = 613.90e-9; LASER.nPhoton = 1; LASER.I0 = 5.5e13;
% LASER.l0 = 613.90e-9; LASER.nPhoton = 1; LASER.I0 = 1e14;
% LASER.l0 = 613.90e-9; LASER.nPhoton = 1; LASER.I0 = 3.33e14;
%..........................................................................



%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% 30a0 WELL
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% LASER.l0 = 800e-9; LASER.nPhoton = 1; LASER.I0 = 0e0;
% 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% E1 -> E2    *******************
%..........................................................................
% LASER.l0 = 2762.5e-9; LASER.nPhoton = 1; LASER.I0 = 5e10;
% LASER.l0 = 2762.5e-9; LASER.nPhoton = 1; LASER.I0 = 2.5e11;
% LASER.l0 = 2762.5e-9; LASER.nPhoton = 1; LASER.I0 = 1e12;
% LASER.l0 = 2762.5e-9; LASER.nPhoton = 1; LASER.I0 = 3.33e12;
% LASER.l0 = 2762.5e-9; LASER.nPhoton = 1; LASER.I0 = 1e13;
% LASER.l0 = 2762.5e-9; LASER.nPhoton = 1; LASER.I0 = 2.5e13;
% LASER.l0 = 2762.5e-9; LASER.nPhoton = 1; LASER.I0 = 5.5e13;
% LASER.l0 = 2762.5e-9; LASER.nPhoton = 1; LASER.I0 = 1e14;
% LASER.l0 = 2762.5e-9; LASER.nPhoton = 1; LASER.I0 = 3.33e14;
%..........................................................................
% LASER.l0 = 4143.8e-9; LASER.nPhoton = 1; LASER.I0 = 5e10;
% LASER.l0 = 4143.8e-9; LASER.nPhoton = 1; LASER.I0 = 2.5e11;
% LASER.l0 = 4143.8e-9; LASER.nPhoton = 1; LASER.I0 = 1e12;
% LASER.l0 = 4143.8e-9; LASER.nPhoton = 1; LASER.I0 = 3.33e12;
% LASER.l0 = 4143.8e-9; LASER.nPhoton = 1; LASER.I0 = 1e13;
% LASER.l0 = 4143.8e-9; LASER.nPhoton = 1; LASER.I0 = 2.5e13;
% LASER.l0 = 4143.8e-9; LASER.nPhoton = 1; LASER.I0 = 5.5e13;
% LASER.l0 = 4143.8e-9; LASER.nPhoton = 1; LASER.I0 = 1e14;
% LASER.l0 = 4143.8e-9; LASER.nPhoton = 1; LASER.I0 = 3.33e14;
%..........................................................................
% LASER.l0 = 5525.0e-9; LASER.nPhoton = 1; LASER.I0 = 5e10;
% LASER.l0 = 5525.0e-9; LASER.nPhoton = 1; LASER.I0 = 2.5e11;
% LASER.l0 = 5525.0e-9; LASER.nPhoton = 1; LASER.I0 = 1e12;
% LASER.l0 = 5525.0e-9; LASER.nPhoton = 1; LASER.I0 = 3.33e12;
% LASER.l0 = 5525.0e-9; LASER.nPhoton = 1; LASER.I0 = 1e13;
% LASER.l0 = 5525.0e-9; LASER.nPhoton = 1; LASER.I0 = 2.5e13;
% LASER.l0 = 5525.0e-9; LASER.nPhoton = 1; LASER.I0 = 5.5e13;
% LASER.l0 = 5525.0e-9; LASER.nPhoton = 1; LASER.I0 = 1e14;
% LASER.l0 = 5525.0e-9; LASER.nPhoton = 1; LASER.I0 = 3.33e14;
%..........................................................................
% LASER.l0 = 1381.3e-9; LASER.nPhoton = 1; LASER.I0 = 5e10;
% LASER.l0 = 1381.3e-9; LASER.nPhoton = 1; LASER.I0 = 2.5e11;
% LASER.l0 = 1381.3e-9; LASER.nPhoton = 1; LASER.I0 = 1e12;
% LASER.l0 = 1381.3e-9; LASER.nPhoton = 1; LASER.I0 = 3.33e12;
% LASER.l0 = 1381.3e-9; LASER.nPhoton = 1; LASER.I0 = 1e13;
% LASER.l0 = 1381.3e-9; LASER.nPhoton = 1; LASER.I0 = 2.5e13;
% LASER.l0 = 1381.3e-9; LASER.nPhoton = 1; LASER.I0 = 5.5e13;
% LASER.l0 = 1381.3e-9; LASER.nPhoton = 1; LASER.I0 = 1e14;
% LASER.l0 = 1381.3e-9; LASER.nPhoton = 1; LASER.I0 = 3.33e14;
%..........................................................................






%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% 40a0 WELL
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% LASER.l0 = 800e-9; LASER.nPhoton = 1; LASER.I0 = 0e0;
% 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% E1 -> E2    *******************
%..........................................................................
% LASER.l0 = 4911.2e-9; LASER.nPhoton = 1; LASER.I0 = 5e10;
% LASER.l0 = 4911.2e-9; LASER.nPhoton = 1; LASER.I0 = 2.5e11;
% LASER.l0 = 4911.2e-9; LASER.nPhoton = 1; LASER.I0 = 1e12;
% LASER.l0 = 4911.2e-9; LASER.nPhoton = 1; LASER.I0 = 3.33e12;
% LASER.l0 = 4911.2e-9; LASER.nPhoton = 1; LASER.I0 = 1e13;
% LASER.l0 = 4911.2e-9; LASER.nPhoton = 1; LASER.I0 = 2.5e13;
% LASER.l0 = 4911.2e-9; LASER.nPhoton = 1; LASER.I0 = 5.5e13;
% LASER.l0 = 4911.2e-9; LASER.nPhoton = 1; LASER.I0 = 1e14;
% LASER.l0 = 4911.2e-9; LASER.nPhoton = 1; LASER.I0 = 3.33e14;
%..........................................................................
% LASER.l0 = 7366.7e-9; LASER.nPhoton = 1; LASER.I0 = 5e10;
% LASER.l0 = 7366.7e-9; LASER.nPhoton = 1; LASER.I0 = 2.5e11;
% LASER.l0 = 7366.7e-9; LASER.nPhoton = 1; LASER.I0 = 1e12;
% LASER.l0 = 7366.7e-9; LASER.nPhoton = 1; LASER.I0 = 3.33e12;
% LASER.l0 = 7366.7e-9; LASER.nPhoton = 1; LASER.I0 = 1e13;
% LASER.l0 = 7366.7e-9; LASER.nPhoton = 1; LASER.I0 = 2.5e13;
% LASER.l0 = 7366.7e-9; LASER.nPhoton = 1; LASER.I0 = 5.5e13;
% LASER.l0 = 7366.7e-9; LASER.nPhoton = 1; LASER.I0 = 1e14;
% LASER.l0 = 7366.7e-9; LASER.nPhoton = 1; LASER.I0 = 3.33e14;
%..........................................................................
% LASER.l0 = 9822.3e-9; LASER.nPhoton = 1; LASER.I0 = 5e10;
% LASER.l0 = 9822.3e-9; LASER.nPhoton = 1; LASER.I0 = 2.5e11;
% LASER.l0 = 9822.3e-9; LASER.nPhoton = 1; LASER.I0 = 1e12;
% LASER.l0 = 9822.3e-9; LASER.nPhoton = 1; LASER.I0 = 3.33e12;
% LASER.l0 = 9822.3e-9; LASER.nPhoton = 1; LASER.I0 = 1e13;
% LASER.l0 = 9822.3e-9; LASER.nPhoton = 1; LASER.I0 = 2.5e13;
% LASER.l0 = 9822.3e-9; LASER.nPhoton = 1; LASER.I0 = 5.5e13;
% LASER.l0 = 9822.3e-9; LASER.nPhoton = 1; LASER.I0 = 1e14;
% LASER.l0 = 9822.3e-9; LASER.nPhoton = 1; LASER.I0 = 3.33e14;
%..........................................................................
% LASER.l0 = 2455.60e-9; LASER.nPhoton = 1; LASER.I0 = 5e10;
% LASER.l0 = 2455.60e-9; LASER.nPhoton = 1; LASER.I0 = 2.5e11;
% LASER.l0 = 2455.60e-9; LASER.nPhoton = 1; LASER.I0 = 1e12;
% LASER.l0 = 2455.60e-9; LASER.nPhoton = 1; LASER.I0 = 3.33e12;
% LASER.l0 = 2455.60e-9; LASER.nPhoton = 1; LASER.I0 = 1e13;
% LASER.l0 = 2455.60e-9; LASER.nPhoton = 1; LASER.I0 = 2.5e13;
% LASER.l0 = 2455.60e-9; LASER.nPhoton = 1; LASER.I0 = 5.5e13;
% LASER.l0 = 2455.60e-9; LASER.nPhoton = 1; LASER.I0 = 1e14;
% LASER.l0 = 2455.60e-9; LASER.nPhoton = 1; LASER.I0 = 3.33e14;
%..........................................................................






% LASER.l0 = 800e-9; LASER.nPhoton = 1; LASER.I0 = 0e0;
% LASER.l0 = 400e-9; LASER.nPhoton = 2; LASER.I0 = 1e13;
% LASER.l0 = 400e-9; LASER.nPhoton = 2; LASER.I0 = 5e13;



% LASER.l0 = 400e-9; LASER.nPhoton = 2; LASER.I0 = 1e12;
% LASER.l0 = 400e-9; LASER.nPhoton = 2; LASER.I0 = 5e12;
% LASER.l0 = 400e-9; LASER.nPhoton = 2; LASER.I0 = 1e13;
% LASER.l0 = 400e-9; LASER.nPhoton = 2; LASER.I0 = 1.5e13;
% LASER.l0 = 400e-9; LASER.nPhoton = 2; LASER.I0 = 2.5e13;
% LASER.l0 = 400e-9; LASER.nPhoton = 2; LASER.I0 = 5e13;


%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% % Hydrogen-1
% ..........................................................................
% LASER.l0 = 800e-9; LASER.nPhoton = 1; LASER.I0 = 0e0;
% LASER.l0 = 400e-9; LASER.nPhoton = 2; LASER.I0 = 5e12;
% LASER.l0 = 400e-9; LASER.nPhoton = 2; LASER.I0 = 1e13;
% LASER.l0 = 400e-9; LASER.nPhoton = 2; LASER.I0 = 1e12;

% LASER.l0 = 121.5e-9; LASER.nPhoton = 8; LASER.I0 = 1e11;
% LASER.l0 = 121.5e-9; LASER.nPhoton = 8; LASER.I0 = 5e12;

% LASER.l0 = 243.0e-9; LASER.nPhoton = 5; LASER.I0 = 2.5e11;
% LASER.l0 = 243.0e-9; LASER.nPhoton = 5; LASER.I0 = 1.25e13;

% LASER.l0 = 121.5e-9; LASER.nPhoton = 8; LASER.I0 = 7.5e11;
% LASER.l0 = 121.5e-9; LASER.nPhoton = 5; LASER.I0 = 7.5e11;
% LASER.l0 = 182.25e-9; LASER.nPhoton = 8; LASER.I0 = 7.5e11;
% LASER.l0 = 182.25e-9; LASER.nPhoton = 5; LASER.I0 = 7.5e11;
% LASER.l0 = 243.0e-9; LASER.nPhoton = 5; LASER.I0 = 7.5e11;
% LASER.l0 = 303.75e-9; LASER.nPhoton = 5; LASER.I0 = 7.5e11;
% LASER.l0 = 364.5e-9; LASER.nPhoton = 5; LASER.I0 = 7.5e10;

% [H1:{1s}-->{2s}]
% % % % % 
% LASER.l0 = 121.5e-9; LASER.nPhoton = 1; LASER.I0 = 5e10;
% LASER.l0 = 121.5e-9; LASER.nPhoton = 1; LASER.I0 = 2.5e11;
% LASER.l0 = 121.5e-9; LASER.nPhoton = 1; LASER.I0 = 1e12;
% LASER.l0 = 121.5e-9; LASER.nPhoton = 1; LASER.I0 = 3.33e12;
% LASER.l0 = 121.5e-9; LASER.nPhoton = 1; LASER.I0 = 2.5e13;
% LASER.l0 = 121.5e-9; LASER.nPhoton = 1; LASER.I0 = 5.5e13;
% LASER.l0 = 121.5e-9; LASER.nPhoton = 1; LASER.I0 = 1e14;
% LASER.l0 = 121.5e-9; LASER.nPhoton = 1; LASER.I0 = 3.33e14;
% LASER.l0 = 121.5e-9; LASER.nPhoton = 1; LASER.I0 = 2e15;
% LASER.l0 = 121.5e-9; LASER.nPhoton = 1; LASER.I0 = 5.5e15;

% % % % %  
% LASER.l0 = 243.0e-9; LASER.nPhoton = 2; LASER.I0 = 5e10;
% LASER.l0 = 243.0e-9; LASER.nPhoton = 2; LASER.I0 = 2.5e11;
% LASER.l0 = 243.0e-9; LASER.nPhoton = 2; LASER.I0 = 1e12;
% LASER.l0 = 243.0e-9; LASER.nPhoton = 2; LASER.I0 = 3.33e12;
% LASER.l0 = 243.0e-9; LASER.nPhoton = 2; LASER.I0 = 2.5e13;
% LASER.l0 = 243.0e-9; LASER.nPhoton = 2; LASER.I0 = 5.5e13;
% LASER.l0 = 243.0e-9; LASER.nPhoton = 2; LASER.I0 = 1e14;
% LASER.l0 = 243.0e-9; LASER.nPhoton = 2; LASER.I0 = 3.33e14;


% % % % % 
% LASER.l0 = 364.5e-9; LASER.nPhoton = 3; LASER.I0 = 5e10;
% LASER.l0 = 364.5e-9; LASER.nPhoton = 3; LASER.I0 = 2.5e11;
% LASER.l0 = 364.5e-9; LASER.nPhoton = 3; LASER.I0 = 1e12;
% LASER.l0 = 364.5e-9; LASER.nPhoton = 3; LASER.I0 = 3.33e12;
% LASER.l0 = 364.5e-9; LASER.nPhoton = 3; LASER.I0 = 2.5e13;
% LASER.l0 = 364.5e-9; LASER.nPhoton = 3; LASER.I0 = 5.5e13;
% LASER.l0 = 364.5e-9; LASER.nPhoton = 3; LASER.I0 = 1e14;
% LASER.l0 = 364.5e-9; LASER.nPhoton = 3; LASER.I0 = 3.33e14;


% LASER.l0 = 400e-9; LASER.nPhoton = 4; LASER.I0 = 2.5e12;
% % % % % 
% LASER.l0 = 182.25e-9; LASER.nPhoton = 2; LASER.I0 = 5e10;
% LASER.l0 = 182.25e-9; LASER.nPhoton = 2; LASER.I0 = 2.5e11;
% LASER.l0 = 182.25e-9; LASER.nPhoton = 2; LASER.I0 = 1e12;
% LASER.l0 = 182.25e-9; LASER.nPhoton = 1; LASER.I0 = 3.33e12;
% LASER.l0 = 182.25e-9; LASER.nPhoton = 1; LASER.I0 = 2.5e13;
% LASER.l0 = 182.25e-9; LASER.nPhoton = 1; LASER.I0 = 5.5e13;
% LASER.l0 = 182.25e-9; LASER.nPhoton = 1; LASER.I0 = 1e14;
% LASER.l0 = 182.25e-9; LASER.nPhoton = 1; LASER.I0 = 3.33e14;


% % % % % 
% LASER.l0 = 60.75e-9; LASER.nPhoton = 1; LASER.I0 = 5e10;
% LASER.l0 = 60.75e-9; LASER.nPhoton = 1; LASER.I0 = 2.5e11;
% LASER.l0 = 60.75e-9; LASER.nPhoton = 1; LASER.I0 = 1e12;
% LASER.l0 = 60.75e-9; LASER.nPhoton = 1; LASER.I0 = 3.33e12;
% LASER.l0 = 60.75e-9; LASER.nPhoton = 1; LASER.I0 = 2.5e13;
% LASER.l0 = 60.75e-9; LASER.nPhoton = 1; LASER.I0 = 5.5e13;
% LASER.l0 = 60.75e-9; LASER.nPhoton = 1; LASER.I0 = 1e14;
% LASER.l0 = 60.75e-9; LASER.nPhoton = 1; LASER.I0 = 3.33e14;






% [H1:{2s}-->{3s}]
% LASER.l0 = 656.11e-9; LASER.nPhoton = 1; LASER.I0 = 1e12;
% LASER.l0 = 656.11e-9; LASER.nPhoton = 1; LASER.I0 = 5e12;
% LASER.l0 = 656.11e-9; LASER.nPhoton = 1; LASER.I0 = 1e13;
% LASER.l0 = 656.11e-9; LASER.nPhoton = 1; LASER.I0 = 5e13;
% LASER.l0 = 656.11e-9; LASER.nPhoton = 1; LASER.I0 = 1e14;
% LASER.l0 = 656.11e-9; LASER.nPhoton = 1; LASER.I0 = 5e14;

% LASER.l0 = 1312.2e-9; LASER.nPhoton = 2; LASER.I0 = 2.5e11;
% LASER.l0 = 1312.2e-9; LASER.nPhoton = 2; LASER.I0 = 1.25e12;
% LASER.l0 = 1312.2e-9; LASER.nPhoton = 2; LASER.I0 = 6.25e12;
% LASER.l0 = 1312.2e-9; LASER.nPhoton = 2; LASER.I0 = 4.25e13;
% LASER.l0 = 1312.2e-9; LASER.nPhoton = 2; LASER.I0 = 2.25e14;
% LASER.l0 = 1312.2e-9; LASER.nPhoton = 2; LASER.I0 = 1.25e15;

% LASER.l0 = 984.165e-9; LASER.nPhoton = 1; LASER.I0 = 1e12;
% LASER.l0 = 984.165e-9; LASER.nPhoton = 1; LASER.I0 = 5e12;
% LASER.l0 = 984.165e-9; LASER.nPhoton = 1; LASER.I0 = 1e13;
% LASER.l0 = 984.165e-9; LASER.nPhoton = 1; LASER.I0 = 5e13;
% LASER.l0 = 984.165e-9; LASER.nPhoton = 1; LASER.I0 = 1e14;
% LASER.l0 = 984.165e-9; LASER.nPhoton = 1; LASER.I0 = 5e14;

% LASER.l0 = 984.165e-9; LASER.nPhoton = 1; LASER.I0 = 2.5e11;
% LASER.l0 = 984.165e-9; LASER.nPhoton = 1; LASER.I0 = 1.25e12;
% LASER.l0 = 984.165e-9; LASER.nPhoton = 1; LASER.I0 = 6.25e12;
% LASER.l0 = 984.165e-9; LASER.nPhoton = 1; LASER.I0 = 4.25e13;
% LASER.l0 = 984.165e-9; LASER.nPhoton = 1; LASER.I0 = 2.25e14;
% LASER.l0 = 984.165e-9; LASER.nPhoton = 1; LASER.I0 = 1.25e15;


% LASER.l0 = 328.06e-9; LASER.nPhoton = 1; LASER.I0 = 1e12;
% LASER.l0 = 328.06e-9; LASER.nPhoton = 1; LASER.I0 = 5e12;
% LASER.l0 = 328.06e-9; LASER.nPhoton = 1; LASER.I0 = 1e13;
% LASER.l0 = 328.06e-9; LASER.nPhoton = 1; LASER.I0 = 5e13;
% LASER.l0 = 328.06e-9; LASER.nPhoton = 1; LASER.I0 = 1e14;
% LASER.l0 = 328.06e-9; LASER.nPhoton = 1; LASER.I0 = 5e14;






% LASER.l0 = 121.5e-9;    % 1 Photon
% LASER.l0 = 243.0e-9;    % 2 Photon
% LASER.l0 = 486.0e-9;    % 4 Photon
% LASER.l0 = 972.0e-9;    % 8 Photon

% [H1:{2s}-->{3s}]
% LASER.l0 = 656.11e-9;   % 1 Photon
% LASER.l0 = 1312.22e-9;  % 2 Photon
% LASER.l0 = 2624.44e-9;  % 4 Photon
% LASER.l0 = 5248.88e-9;  % 8 Photon

% [H1:{1s}-->{3s}]
% LASER.l0 = 102.52e-9;   % 1 Photon
% LASER.l0 = 205.04e-9;   % 2 Photon
% LASER.l0 = 410.08e-9;   % 4 Photon
% LASER.l0 = 820.16e-9;   % 8 Photon



% LASER.l0 = 12.15e-9;    % 1 Photon

%**************************************************************************
% Chirped-Pulse (CEP)
%--------------------------------------------------------------------------
LASER.Type = 'Continuous-Wave';
% LASER.Type = 'Chirped-Pulse';


LASER.f0 = c0/LASER.l0;
LASER.w0 = 2*pi*LASER.f0;
LASER.Phase = 0;%pi/4;
LASER.T = 1/LASER.f0;


LASER.k0 = 2*pi/LASER.l0;

% LASER.Delay = 10*Ta;%5e-16;%LASER.PulseLength/2;
% LASER.Delay = 10*Ta;
% LASER.Delay = 0.1*TIME.saveNum*(TIME.delta*TIME.N/(TIME.saveNum));

LASER.Delay = 0;

% LASER.PERIODcutoff = 0;
LASER.PERIODcutoff = LASER.T;


% LASER.nPhoton = 1;
LASER.FINALcutoff = LASER.nPhoton*LASER.T;

% LASER.FINALcutoff = 0;
% LASER.FINALcutoff = 1*LASER.T;
% LASER.FINALcutoff = 2*LASER.T;
% LASER.FINALcutoff = 4*LASER.T;
% LASER.FINALcutoff = 8*LASER.T;


% Carrier Wave
LASER.EnvelopeWave = 800e-9;
LASER.EnvelopeFreq = 2*pi*c0/LASER.EnvelopeWave;
LASER.Ncycles = 3;
LASER.PulseLength = LASER.Ncycles*(2*pi/LASER.EnvelopeFreq);

LASER.Offset = 1e-15;%LASER.PulseLength*2;
% LASER.Toffset = round(LASER.offset/TIME.delta);



LASER.Time = LASER.PulseLength+LASER.Offset;

LASER.CEP = pi/8;
LASER.CEOF = 0;


LASER.NumPulses = floor(TIME.length/(LASER.Time)-LASER.Delay/LASER.PulseLength);

LASER.Vgroup = 0;
LASER.Vphase = 0;


%**************************************************************************
% Field Intensity	(W/cm^-2)
%--------------------------------------------------------------------------

% LASER.I0 = 0e00;



%**************************************************************************
% Field Strength (Amplitude & Potential)
%--------------------------------------------------------------------------
LASER.I0 = LASER.I0;
% LASER.I0 = LASER.I0*1e2;

% LASER.E0 = sqrt(2*LASER.I0/(c0*eps0));
% LASER.E0 = sqrt(2*(1e2)*LASER.I0/(c0*eps0));
LASER.E0 = sqrt(2*(1e4)*LASER.I0/(c0*eps0));

% LASER.E0 = sqrt(2*(LASER.I0)/(eps0));
LASER.A0 = LASER.E0/LASER.l0;



% Laser Params
% Laser.lambda0 = 120e-9;
% Laser.lambda0 = 180e-9;
% Laser.lambda0 = 240*10^-9;
% Laser.lambda0 = 300*10^-9;


% Laser.f0 = c0/Laser.lambda0;
% Laser.w0 = 2*pi*Laser.f0;
% Laser.repeat = 80*10^9;
% Laser.T = 1/Laser.f0;


% Amplitude
% if DEBUG.LASER
%     Laser.E0 = 1e7;
%     Laser.E0 = 5e6;
%     Laser.E0 = 1e6;
%     Laser.E0 = 5e5;
% else
%     Laser.E0=0;
% end

% Laser.I0 = 0.5*eps0*c0*Laser.E0^2;

% % Photon Energy
% Laser.photon_E = hPlanck_full*c0/Laser.lambda0;
% % Momentum
% Laser.photon_p = hPlanck_full/Laser.lambda0;


% Strings



%% Observables
% Position-Space

%**************************************************************************
% Potential Energy
%--------------------------------------------------------------------------
DEBUG.Observables_R.Potential_Energy = true;
% DEBUG.Observables_R.Potential_Energy = false;

DEBUG.Observables_R.Kinetic_Energy = true;
% DEBUG.Observables_R.Kinetic_Energy = false;


%**************************************************************************
% Probability Density Function
%--------------------------------------------------------------------------
DEBUG.Observables_R.PDF = true;
% DEBUG.Observables_R.PDF = false;


%**************************************************************************
% Centroid
%--------------------------------------------------------------------------
% DEBUG.Observables_R.PDF = true;
% DEBUG.Observables_R.PDF = false;



% Observables in Momentum & K-Space


%**************************************************************************
% Potential Energy
%--------------------------------------------------------------------------
DEBUG.Observables_F.FreqShift = true;
% DEBUG.Observables_F.FreqShift = false;














%% Stability


% V_eff = zeros(SPACE.N, 1);
V_sys = PotentialCoulomb(SIM, SPACE, PULSE);
% V_eff = V_atom;

% Maximum 
V_max = (1+e0*LASER.E0)*max(V_sys);
V_min = (1+e0*LASER.E0)*min(V_sys);


T_crit = (hPlanck/((hPlanck^2/(me*SPACE.delta^2))+V_max));

if T_crit<=TIME.delta
    Q_str = sprintf('WARNING! Time step is too large for spatial step and maximum potential over grid. Simulation will become unstable.');
    Q_opt1 = sprintf('Halt generation and reduce time step');
    Q_opt2 = sprintf('Continue anyways (probable failure)');
    hCRIT = questdlg(Q_str, 'Simulation Unstable!', Q_opt1, Q_opt2, Q_opt1);
    
    switch hCRIT
        case Q_opt1
        Tcrit_str = sprintf(['\tTime Step >= Critical Time Step...\n\n'...
            '\t-->\t\tGeneration halted. Please reduce time step and try again.\n\n']);
        error(Tcrit_str);
            break;
            
        case Q_opt2
        Tcrit_str = sprintf(['\tTime Step >= Critical Time Step...\n\n'...
            '\t-->\t\tCAUTION! Time step very likely to result in simulation instability.\n\n']);
        warning(Tcrit_str);
            
    end
    
else
    Tcrit_str = sprintf(['\tTime Step >= Critical Time Step...\n\n'...
            '\t-->\t\tPASSED! Simulation will be stable.\n\n']);
        fprintf(Tcrit_str);
end




%% SAVE SETTINGS

% Memory Management
VarMem = TIME.saveNum*SPACE.N^SPACE.Dims*8;
[userview systemview] = memory;
EffMem = systemview.PhysicalMemory.Total;
if (VarMem/EffMem)>.4
            clc;
        txt_warn = sprintf(['\tCANCELING...\n\n'...
            '\t-->\t\tVariables generated will be too large to generate plots.\n\n'...
            '\t-->\t\tReduce the number of saved Space-Time steps.\n\n'...
            '\t-->\t\tTo reattempt operation, please reduce the number of saved Space-Time steps and re-run "GENERATE_SETTINGS()"']);
        error(txt_warn);
%         clearvars;
end
    
% SaveRoot = uigetdir(['../___DATA'],'Please determine root filepath');

DEBUG.prefix = '../___Simulation_Parameter_Files';
if ~exist(DEBUG.prefix, 'dir')
    mkdir(DEBUG.prefix);
end


if DEBUG.DetailedFileName
% Laser Parameters Text
switch LASER.Type
    
    case {'Continuous-Wave', 'CW', 'continuous'}
        
        LaserParameters_txt = ...
    sprintf(['E(Io=' num2str(LASER.I0, '%2.2e') ',wl=' ...
    num2str(LASER.l0, '%2.1e') ',(cyc=' num2str(LASER.nPhoton, '%2.1f') ')']);


    case {'Chirped-Pulse', 'chirped-pulse', 'CEP', 'CPA', 'CP'}
        
        LaserParameters_txt = ...
    sprintf(['(Io=' num2str(LASER.I0, '%2.2e') ',wl=' ...
    num2str(LASER.l0, '%2.1e') ',(cyc=' num2str(LASER.nPhoton, '%2.1f') ')']);
    

    otherwise   % Default = Continuous
        
        LaserParameters_txt = ...
    sprintf(['(Io=' num2str(LASER.I0, '%2.2e') ',wl=' ...
    num2str(LASER.l0, '%2.1e') ',(cyc=' num2str(LASER.nPhoton, '%2.1f') ')']);
        

end
        
        
    
% num2str(LASER.l0, '%2.1e') ',(cyc=' num2str(LASER.nCycle, '%2.1f') ')']);

% LaserParameters_txt = ...
%     sprintf(['Laser[I' num2str(LASER.I0, '%3.2e') '),CEP(' ...
%     num2str(LASER.l0, '%2.1e') '),Ncyc(' num2str(LASER.Ncycles, '%2.0f') '),'...
%     'Op(' num2str(LASER.PERIODcutoff/LASER.T) '),Of(' num2str(LASER.FINALcutoff/LASER.T) ')]']);


switch PULSE.Type
%     case {'Sine'; 'Cosine'}
%         clear en_state; en_state = num2str(PULSE.EnergyState(1));
%         if PULSE.NumStates > 1
%             for n=2:PULSE.NumStates
%                 en_state =...
%                     sprintf([en_state '+' num2str(PULSE.EnergyState(n))]);
%             end
%         end
%         energy_txt = ['Psi0(Sine,En=' en_state ')'];

    case {'Sine'; 'Sin'; 'Cosine'; 'Cos'}
        clear eig_state; eig_state = num2str(PULSE.EnergyState(1));
        if PULSE.NumStates > 1
            for n=2:PULSE.NumStates
                eig_state =...
                    sprintf([eig_state '+' num2str(PULSE.EnergyState(n))]);
            end
        end
        InitialWave_txt = ['Psi0(sin,En=' eig_state ')'];
        
        
    case 'Gaussian'
        InitialWave_txt =...
            ['Psi0(gs,Xw=' num2str(PULSE.Width, '%2.1e') ')'];
%         energy_txt =...
%             ['Eg(Xw=' num2str(PULSE.Width, '%3.2e') ',Ko=' ...
%             num2str(PULSE.Momentum, '%3.2e') ')'];
        
        
        
    case {'EigenState'}
        clear eig_state; eig_state = num2str(PULSE.EnergyState(1));
        if PULSE.NumStates > 1
            for n=2:PULSE.NumStates
                eig_state =...
                    sprintf([eig_state '+' num2str(PULSE.EnergyState(n))]);
            end
        end
        InitialWave_txt = ['Psi0(eigen,(En=' eig_state ')'];
    
    case {'Exact'; 'Analytical'}
        clear eig_state; eig_state = num2str(PULSE.EnergyState(1));
        if PULSE.NumStates > 1
            for n=2:PULSE.NumStates
                eig_state =...
                    sprintf([eig_state '+' num2str(PULSE.EnergyState(n))]);
            end
        end
        InitialWave_txt = ['Psi0(An,En=' eig_state ')'];
    
        
        
    otherwise     
        clc;
        txt_warn = sprintf(['\tCANCELING...\n\n'...
            '\t-->\t\tParameter file for simulation will NOT be created!\n\n'...
            '\t-->\t\tTo reattempt operation, please re-run "GENERATE_SETTINGS()"']);
        warning(txt_warn);
        clearvars;
        
        
end




% Boundary Conditions Text

% if ITP.Function
if length(SIM.BC.Type) >= 5;
    txt_len = 4;
else
    txt_len = length(SIM.BC.Type);
end

BCs_type_txt = ...
    sprintf(['(' SIM.BC.Type(1:txt_len) ')']);%...
txt_region = num2str(SIM.BC.N,'%02.0f');

BCs_txt = ...
    sprintf(['BC[' BCs_type_txt ']']);
% sprintf(['BC[' BCs_type_txt ',Rabs=(' txt_region ')]']);

% end




% Imaginary Time Propagation Text

if ITP.Function
if length(ITP.Type) >= 5;
    txt_len = 5;
else
    txt_len = length(ITP.Type);
end

ITP_txt = ...
    sprintf(['itp(' ITP.Type(1:txt_len) ',(n.l.m='...
    num2str(ITP.N_, '%01.0f') '.'...
    num2str(ITP.L_, '%01.0f') '.'...
    num2str(ITP.M_, '%01.0f') ')']);%...
txt_Nth_Order = num2str(ITP.NthOrder,'%02.0f');
switch txt_Nth_Order(end)
    case '1'
        if ~strcmp(txt_Nth_Order(end-1),'1'); txt_suffix = 'st';
        else txt_suffix = 'th'; end
    case '2'
        if ~strcmp(txt_Nth_Order(end-1),'1'); txt_suffix = 'nd';
        else txt_suffix = 'th'; end        
    case '3'
        if ~strcmp(txt_Nth_Order(end-1),'1'); txt_suffix = 'rd';
        else txt_suffix = 'th'; end        
    otherwise
        txt_suffix = 'th';
end
%     sprintf([ITP_txt ',(N=' num2str(ITP.NthOrder, '%02.0f') '' txt_suffix ')]']);
ITP_txt = ...
sprintf([ITP_txt ',(N=' num2str(ITP.NthOrder, '%02.0f') '' txt_suffix ', ' ...
    num2str(ITP.ConvergenceRatio, '%2.2e') ', '...
    num2str(ITP.TimeDivisor, '%2.2e') ')]']);
else
    ITP_txt = '';
%     ITP_txt = '(itp=Off)';
end


% Fourier Analysis Text
FFT_txt = ...
    sprintf(['(N0=' num2str(SPACE.ZeroPad_Factor, '%01.f') ')']);
% sprintf(['(N0=' num2str(SPACE.ZeroPad_Factor, '%01.1f') 'N)']);

PotentialEnergy_txt = sprintf(['[V(' SIM.PotentialMap ...
        ')]']);
%     '), cr(' num2str(PULSE.Creg/a0, '%1.e') ')]']);

CalcStates_txt = ['Ens(' num2str(PULSE.NumCalcStates,'%1.f') ')'];

Analysis_txt = sprintf(['[' FFT_txt ',' ITP_txt ']']);

% Solution Space Text
SolutionSpace_txt = sprintf(['Sys[X(' num2str(SPACE.Dims, '%01.0f') 'd,'...
    'x=' num2str(SPACE.length/a0, '%2.0f') 'a0,dx=' num2str(SPACE.N) ','...
    'T(t=' num2str(TIME.length, '%3.2e') ',dt=' num2str(TIME.N, '%1.1e') ')]']);

EnergyParams_txt = sprintf(['SRC[' ...
    LaserParameters_txt ',' PotentialEnergy_txt ',' InitialWave_txt ',' CalcStates_txt ',']);


ParamFileName_txt = [SolutionSpace_txt ',' EnergyParams_txt ',' Analysis_txt ',' BCs_txt];

else
    ShortNameNum = num2str(length(dir(fullfile(DEBUG.prefix,'*.mat'))),'%04.0f');
    ParamFileName_txt = ['SIM_parameter_file(' num2str(ShortNameNum) ')'];

end


%% Save

[ChooseParamName, ChooseParamPath] = uiputfile([DEBUG.prefix '/' ParamFileName_txt '.mat']);

% [ChooseParamName, ChooseParamPath] = uiputfile(['../___Simulation_Parameter_Files/*.mat']);
if ~(ChooseParamName==0)
    save([ChooseParamPath '\' ChooseParamName]);
else
    clc;
    txt_warn = sprintf(['\tCANCELING...\n\n'...
        '\t-->\t\tParameter file for simulation will NOT be created!\n\n'...
        '\t-->\t\tTo reattempt operation, please re-run "GENERATE_SETTINGS()"']);
    warning(txt_warn);
    clearvars;
    
end


% SavePath = uisave(['../___Simulation_Parameter_Files'],'Please Choose File Name');














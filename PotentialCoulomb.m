function [V_coulomb] = PotentialCoulomb(SIM, SPACE, PULSE)

% VarName, VarIO, MemIndex)


% [FormatString] = SaveFormat(VarName, VarType, VarSize, FilePath)
% evalstr = SaveFormat('Psi_REAL', 1, [(1e6) (256) (256)], fpath);
%% SaveFormat
% INPUT % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%
% VarName
%	'Psi_REAL'
%
% VarSize
%	[1e4 64 64 64]
%
% VarIO
%	'FDTDQ_file
%
% idx
%	idx.START   = 5e5
%   idx.END     = 1e6
%   idx.length  = 5e5
%
%
% OUTPUT % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%
% FormatString
%	['FDTDQ_file.IO{n}.' FDTDQ_file.vars{n} ' = ' FDTDQ_file.vars{n} ';']
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%


%% Sanity Check

% Input Arguments
narginchk(3,3);

% Output Arguments
nargoutchk(1,1);

%%
% VarName, VarIO, MemIndex


%% Local variable assignment
% SPACE     =       varargin{1};

% if nargin > 1
% WellType     =       varargin{1};
% end


%%
CONSTANTS;

% N = length(SPACE);

% InitializeString = sprintf([VarName ' = zeros(' num2str(VarSize(1))]);
% 
% % N-dimensions
% n_ = length(VarSize);
% for n=2:n_
%     InitializeString = sprintf([InitializeString ', ' num2str(VarSize(n))]);
% end
% 
% InitializeString = sprintf([InitializeString ', ''like'', ' num2str(VarType) ');']);
% 
% eval(InitializeString);




%% Coulomb Potential
% [V_coulomb] = CoulombPotential();

% SPACE_dist = abs(SPACE(end)-SPACE(1));

% Well Width

% well_width = SPACE_dist;
% well_width = 2.5*a0;
% well_width = 10*a0;


% Well Height
% well_height = 10^6;
% well_height = 2.5*a0;
% well_height = 10*a0;


% V_abs = (1/4)*1i*(delta_t/delta_x^2);
% # of Protons
Z = 1;
q = Z*e0^2;
K_coulomb = 1/(4*pi*eps0);
% coulomb_regularized = 3.33e-2*a0;

% R_ = (-SPACE.length/2:SPACE.delta:SPACE.length/2);


%% 

    
switch SIM.PotentialMap

%**************************************************************************
    case {'H_1';'H1'}
%--------------------------------------------------------------------------
% # of Protons
Z = 1;
q = Z*e0^2;

% Regularized Coulomb
%         r_ = sqrt(SPACE.R^2+PULSE.Creg^2);
%         V_coulomb = -q*K_coulomb./r_;

V_coulomb = zeros(1,SPACE.N);
% V_coulomb = zeros(1,SPACE.N_);
for x=2:(SPACE.N-1)
    V_coulomb(x) =...
        -q*K_coulomb/sqrt((SPACE.Axis(x))^2+PULSE.Creg^2);
%     -q*K_coulomb/sqrt((SPACE.Axis_(x))^2+PULSE.Creg^2);
end

% Vcoulomb = zeros(SPACE.N,SPACE.N,SPACE.N);
% for x=1:SPACE.N; for y=1:SPACE.N; for z=1:SPACE.N
% Vcoulomb(x,y,z) = -q*K_coulomb/...
% sqrt((SPACE.R(x))^2+(SPACE.R(y))^2+(SPACE.R(z))^2+PULSE.Creg^2);
% end; end; end

% Vcoulomb = -abs(Vcoulomb./e0);





% REGULARIZED
% if min(min(min(V_t)))<0
%     Vreg = (SPACE.delta/a0);
%     r_ = 1/Vreg;
% 
%     
% for r=1:SPACE.N
%     if r>r_
%         V_t(r) = V_t(r)/abs(e0) + Vreg;
%     else
%         V_t(r) = 0;
%     end
% end   
% else
%     Vreg = 0;
% end








%**************************************************************************
% Square Well
%..........................................................................
    case {'Square'; 'InfSqWell'; 'Inf'}
%--------------------------------------------------------------------------
% Wall.begin = round(SPACE.N/2);
% Wall.R = round(SPACE.N/2);

%         V_coulomb = e0*ones(1,SPACE.N);
V_coulomb = zeros(1,SPACE.N);
% V_coulomb = ones(1,SPACE.N_).*0;
%         Vcoulomb = zeros(SPACE.N,SPACE.N,SPACE.N);
 

%..........................................................................
% Here we have a square well with a floor of (e0), as opposed to zero.
    case 'Square_'
%--------------------------------------------------------------------------
% Wall.begin = round(SPACE.N/2);
% Wall.R = round(SPACE.N/2);

%         V_coulomb = e0*ones(1,SPACE.N);
V_coulomb = ones(1,SPACE.N).*e0;
% V_coulomb = ones(1,SPACE.N_).*e0;
%         Vcoulomb = zeros(SPACE.N,SPACE.N,SPACE.N);


%**************************************************************************
    case {'QHO';'QuantumHarmonicOscillator';'HarmonicOscillator'}
%..........................................................................
% Harmonic Oscillator potentials contain a parameter dependent on the mass
% of the particles interacting within them. So if more particles are
% present within the system, the mass term will 
%--------------------------------------------------------------------------
Nparticles = 1;
SystemMass = Nparticles*me;
V_coulomb = zeros(1,SPACE.N);
w0 = e0/hPlanck;

V_coulomb = 0.5*SystemMass*(w0^2)*(SPACE.Axis.^2);


%**************************************************************************
% None Found...
%--------------------------------------------------------------------------
    otherwise
        warning(sprintf([':\tPOTENTIAL WELL UNDEFINED'...
        '\n\n***\tPlease select potential map and run again.'...
        '\n***\tIf error persists, check & ensure well is formulated correctly.']));
end




%% RETURN
% varname = 'V_coulomb';


% evalin('caller',...
%     ['V_coulomb(1:SPACE.N) = ' num2str(Laser_E) '.*' Laser_coupling ';']);







end
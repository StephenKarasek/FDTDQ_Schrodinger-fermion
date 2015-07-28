function [DATA] = FDM_vars(SPACE, TIME, SIM, PULSE, ITP, LASER, DATA, DEBUG)

%% FDM

% sanity check

% if nargin < 1
%     CMD = '';

%     switch CMD
    
%         case 'define'
            
            
%%
SaveDirectory = DATA.PARAMS_.saveDir;
%**************************************************************************
%   Real Wavefunction
%--------------------------------------------------------------------------
VarName = 'Psi_REAL';
VarSize = [TIME.saveNum SPACE.N];
% VarSize = [TIME.saveNum SPACE.N_];
VarType = 1;
VarLoc = 'DATA.FDTD_';

eval([VarLoc '.' (VarName) ...
    ' = SaveInitialize(SaveDirectory, VarName, VarSize, VarType, VarLoc);']);
fprintf('. ');


%**************************************************************************
%   Imaginary Wavefunction
%--------------------------------------------------------------------------
% Initialize
VarName = 'Psi_IMAG';
VarSize = [TIME.saveNum SPACE.N];
% VarSize = [TIME.saveNum SPACE.N_];
VarType = 1i;
VarLoc = 'DATA.FDTD_';

eval([VarLoc '.' (VarName) ...
    ' = SaveInitialize(SaveDirectory, VarName, VarSize, VarType, VarLoc);']);
fprintf('. ');


%**************************************************************************
%   Probability Density Function
%--------------------------------------------------------------------------
if DEBUG.Observables_R.PDF
% Initialize
VarName = 'PDF_pos';
VarSize = [TIME.saveNum SPACE.N];
% VarSize = [TIME.saveNum SPACE.N_];
VarType = 1;
VarLoc = 'DATA.FDTD_';

eval([VarLoc '.' (VarName) ...
    ' = SaveInitialize(SaveDirectory, VarName, VarSize, VarType, VarLoc);']);
end              
fprintf('. ');


%**************************************************************************
%   Potential Energy
%--------------------------------------------------------------------------
% Initialize
VarName = 'V_pot';
VarSize = [TIME.saveNum SPACE.N];
% VarSize = [TIME.saveNum SPACE.N_];
VarType = 1;
VarLoc = 'DATA.FDTD_';

eval([VarLoc '.' (VarName) ...
    ' = SaveInitialize(SaveDirectory, VarName, VarSize, VarType, VarLoc);']);
fprintf('. ');


%**************************************************************************
%   Atomic Potential
%--------------------------------------------------------------------------
% Initialize
VarName = 'V_atom';
VarSize = [1 SPACE.N];
% VarSize = [1 SPACE.N_];
VarType = 1;
VarLoc = 'DATA.FDTD_';

eval([VarLoc '.' (VarName) ...
    ' = SaveInitialize(SaveDirectory, VarName, VarSize, VarType, VarLoc);']);
fprintf('. ');


%**************************************************************************
%   Interaction Potential
%--------------------------------------------------------------------------
% Initialize
VarName = 'V_laser';
VarSize = [TIME.saveNum SPACE.N];
% VarSize = [TIME.saveNum SPACE.N_];
VarType = 1;
VarLoc = 'DATA.FDTD_';

eval([VarLoc '.' (VarName) ...
    ' = SaveInitialize(SaveDirectory, VarName, VarSize, VarType, VarLoc);']);
fprintf('. ');


%**************************************************************************
%   Interaction Potential
%--------------------------------------------------------------------------
% Initialize
VarName = 'E_laser';
VarSize = [TIME.saveNum 1];
VarType = 1;
VarLoc = 'DATA.FDTD_';

eval([VarLoc '.' (VarName) ...
    ' = SaveInitialize(SaveDirectory, VarName, VarSize, VarType, VarLoc);']);
fprintf('. ');

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 


%**************************************************************************
%   EigenFunction
%--------------------------------------------------------------------------
% Initialize
VarName = 'Psi_eig';
VarSize = [ITP.FindStates SPACE.N];
% VarSize = [ITP.FindStates SPACE.N_];
VarType = 1;
VarLoc = 'DATA.FDTD_';

eval([VarLoc '.' (VarName) ...
    ' = SaveInitialize(SaveDirectory, VarName, VarSize, VarType, VarLoc);']);
fprintf('. ');

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

%**************************************************************************
%   EigenCoefficient
%--------------------------------------------------------------------------
% Initialize
VarName = 'Coeff_eig';
VarSize = [ITP.FindStates 1];
VarType = 1;
VarLoc = 'DATA.FDTD_';

eval([VarLoc '.' (VarName) ...
    ' = SaveInitialize(SaveDirectory, VarName, VarSize, VarType, VarLoc);']);
fprintf('. ');

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

%**************************************************************************
%   EigenEnergy
%--------------------------------------------------------------------------
% Initialize
VarName = 'Energy_eig';
VarSize = [ITP.FindStates 1];
VarType = 1;
VarLoc = 'DATA.FDTD_';

eval([VarLoc '.' (VarName) ...
    ' = SaveInitialize(SaveDirectory, VarName, VarSize, VarType, VarLoc);']);
fprintf('. ');



% 
% %**************************************************************************
% %   Energy State (wave)
% %--------------------------------------------------------------------------
% % Initialize
% VarName = 'V_laser';
% VarSize = [TIME.saveNum SPACE.N];%SPACE.N];
% VarType = 1;
% VarLoc = 'DATA.FDTD_';
% 
% eval([VarLoc '.' (VarName) ...
%     ' = SaveInitialize(SaveDirectory, VarName, VarSize, VarType, VarLoc);']);
% 



end
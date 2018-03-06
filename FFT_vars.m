function [DATA] = FFT_vars(SPACE, TIME, SIM, PULSE, ITP, LASER, DATA, DEBUG)
%%






SaveDirectory = DATA.PARAMS_.saveDir;
%%
%**************************************************************************
%   Real Wavefunction(t,p)
%--------------------------------------------------------------------------
VarName = 'fft_REAL';
% VarSize = [TIME.saveNum SPACE.N];
VarSize = [TIME.saveNum SPACE.N0];
VarType = 1;
VarLoc = 'DATA.FFT_';

eval([VarLoc '.' (VarName) ...
    ' = SaveInitialize(SaveDirectory, VarName, VarSize, VarType, VarLoc);']);
fprintf('. ');


%**************************************************************************
%   Imaginary Wavefunction(t,p)
%--------------------------------------------------------------------------
% Initialize
VarName = 'fft_IMAG';
% VarSize = [TIME.saveNum SPACE.N];
VarSize = [TIME.saveNum SPACE.N0];
VarType = 1i;
VarLoc = 'DATA.FFT_';

eval([VarLoc '.' (VarName) ...
    ' = SaveInitialize(SaveDirectory, VarName, VarSize, VarType, VarLoc);']);
fprintf('. ');


%**************************************************************************
%   Freq Centered Wavefunction(t,k)
if DEBUG.Observables_F.FreqShift
%--------------------------------------------------------------------------
% Initialize
VarName = 'FreqShift';
% VarSize = [TIME.saveNum SPACE.N];
VarSize = [TIME.saveNum SPACE.N0];
VarType = 1;
VarLoc = 'DATA.FFT_';

eval([VarLoc '.' (VarName) ...
    ' = SaveInitialize(SaveDirectory, VarName, VarSize, VarType, VarLoc);']);
fprintf('. ');
end

%**************************************************************************
%   Normalization Constant
%--------------------------------------------------------------------------
% Initialize
VarName = 'NormF_ALL';
VarSize = [TIME.saveNum 1];%SPACE.N];
VarType = 1;
VarLoc = 'DATA.FFT_';

eval([VarLoc '.' (VarName) ...
    ' = SaveInitialize(SaveDirectory, VarName, VarSize, VarType, VarLoc);']);
fprintf('. ');


%**************************************************************************
%   Normalization Constant (first few states)
%--------------------------------------------------------------------------
% Initialize
VarName = 'NormF_CALC';
VarSize = [TIME.saveNum 1];%SPACE.N];
VarType = 1;
VarLoc = 'DATA.FFT_';

eval([VarLoc '.' (VarName) ...
    ' = SaveInitialize(SaveDirectory, VarName, VarSize, VarType, VarLoc);']);
fprintf('. ');

%**************************************************************************
%   Normalization Constant (States in Focus)
%--------------------------------------------------------------------------
% Initialize
VarName = 'NormF_VIEW';
VarSize = [TIME.saveNum 1];%SPACE.N];
VarType = 1;
VarLoc = 'DATA.FFT_';

eval([VarLoc '.' (VarName) ...
    ' = SaveInitialize(SaveDirectory, VarName, VarSize, VarType, VarLoc);']);
fprintf('. ');

%**************************************************************************
%   Kinetic Energy (REAL)
%--------------------------------------------------------------------------
% Initialize
VarName = 'T_kin_REAL';
VarSize = [TIME.saveNum SPACE.N];
% VarSize = [TIME.saveNum SPACE.N_];
VarType = 1;
VarLoc = 'DATA.FFT_';

eval([VarLoc '.' (VarName) ...
    ' = SaveInitialize(SaveDirectory, VarName, VarSize, VarType, VarLoc);']);
fprintf('. ');


%**************************************************************************
%   Kinetic Energy (IMAG)
%--------------------------------------------------------------------------
% Initialize
VarName = 'T_kin_IMAG';
VarSize = [TIME.saveNum SPACE.N];
% VarSize = [TIME.saveNum SPACE.N_];
VarType = 1;
VarLoc = 'DATA.FFT_';

eval([VarLoc '.' (VarName) ...
    ' = SaveInitialize(SaveDirectory, VarName, VarSize, VarType, VarLoc);']);
fprintf('. ');


%**************************************************************************
%   Normalization Constant
%--------------------------------------------------------------------------
% Initialize
VarName = 'NormT';
VarSize = [TIME.saveNum 1];%SPACE.N];
VarType = 1;
VarLoc = 'DATA.FFT_';

eval([VarLoc '.' (VarName) ...
    ' = SaveInitialize(SaveDirectory, VarName, VarSize, VarType, VarLoc);']);
fprintf('. ');


%**************************************************************************
%   Instantaneous Energy (REAL)
%--------------------------------------------------------------------------
% Initialize
VarName = 'En_REAL';
VarSize = [TIME.saveNum 1];%SPACE.N];
VarType = 1;
VarLoc = 'DATA.FFT_';

eval([VarLoc '.' (VarName) ...
    ' = SaveInitialize(SaveDirectory, VarName, VarSize, VarType, VarLoc);']);
fprintf('. ');


%**************************************************************************
%   Instantaneous Energy (REAL)
%--------------------------------------------------------------------------
% Initialize
VarName = 'En_IMAG';
VarSize = [TIME.saveNum 1];%SPACE.N];
VarType = 1;
VarLoc = 'DATA.FFT_';

eval([VarLoc '.' (VarName) ...
    ' = SaveInitialize(SaveDirectory, VarName, VarSize, VarType, VarLoc);']);
fprintf('. ');

%**************************************************************************
%   Instantaneous Energy (REAL)
%--------------------------------------------------------------------------
% Initialize
VarName = 'En_tot';
VarSize = [TIME.saveNum 1];%SPACE.N];
VarType = 1i;
VarLoc = 'DATA.FFT_';

eval([VarLoc '.' (VarName) ...
    ' = SaveInitialize(SaveDirectory, VarName, VarSize, VarType, VarLoc);']);
fprintf('. ');



%**************************************************************************
%   Instantaneous Energy (REAL)
%--------------------------------------------------------------------------
% Initialize
VarName = 'PhotonsAbsorbed';
VarSize = [TIME.saveNum 1];%SPACE.N];
VarType = 1;
VarLoc = 'DATA.FFT_';

eval([VarLoc '.' (VarName) ...
    ' = SaveInitialize(SaveDirectory, VarName, VarSize, VarType, VarLoc);']);
fprintf('. ');


%**************************************************************************
%   Instantaneous Ionization Rate
%--------------------------------------------------------------------------
% Initialize
VarName = 'IIR';
VarSize = [TIME.saveNum 1];%SPACE.N];
VarType = 1;
VarLoc = 'DATA.FFT_';

eval([VarLoc '.' (VarName) ...
    ' = SaveInitialize(SaveDirectory, VarName, VarSize, VarType, VarLoc);']);
fprintf('. ');


%**************************************************************************
%   Population State (normalized to ALL states)
%--------------------------------------------------------------------------
% Initialize
VarName = 'PopulationState';
VarSize = [TIME.saveNum PULSE.NumCalcStates];
% VarSize = [TIME.saveNum PULSE.NumViewStates];
% VarSize = [TIME.saveNum ITP.AccurateStates];%SPACE.N];
VarType = 1;
VarLoc = 'DATA.FFT_';

eval([VarLoc '.' (VarName) ...
    ' = SaveInitialize(SaveDirectory, VarName, VarSize, VarType, VarLoc);']);
fprintf('. ');

%**************************************************************************
%   Population State (normalized to first 5 states)
%--------------------------------------------------------------------------
% Initialize
VarName = 'PopulationStateNorm_CALC';
% VarSize = [TIME.saveNum PULSE.NumCalcStates];
VarSize = [TIME.saveNum PULSE.NumCalcStates];
% VarSize = [TIME.saveNum ITP.AccurateStates];%SPACE.N];
VarType = 1;
VarLoc = 'DATA.FFT_';

eval([VarLoc '.' (VarName) ...
    ' = SaveInitialize(SaveDirectory, VarName, VarSize, VarType, VarLoc);']);
fprintf('. ');

%**************************************************************************
%   Population State (normalized to first 5 states)
%--------------------------------------------------------------------------
% Initialize
VarName = 'PopulationStateNorm_VIEW';
% VarSize = [TIME.saveNum PULSE.NumCalcStates];
VarSize = [TIME.saveNum PULSE.NumViewStates];
% VarSize = [TIME.saveNum ITP.AccurateStates];%SPACE.N];
VarType = 1;
VarLoc = 'DATA.FFT_';

eval([VarLoc '.' (VarName) ...
    ' = SaveInitialize(SaveDirectory, VarName, VarSize, VarType, VarLoc);']);
fprintf('. ');


%**************************************************************************
%   Population State (normalized to first 5 states)
%--------------------------------------------------------------------------
% Initialize
VarName = 'PopulationStateNorm_GIVEN';
% VarSize = [TIME.saveNum PULSE.NumCalcStates];
VarSize = [TIME.saveNum PULSE.SpreadStates];
% VarSize = [TIME.saveNum ITP.AccurateStates];%SPACE.N];
VarType = 1;
VarLoc = 'DATA.FFT_';

eval([VarLoc '.' (VarName) ...
    ' = SaveInitialize(SaveDirectory, VarName, VarSize, VarType, VarLoc);']);
fprintf('. ');


%**************************************************************************
%   Centroid (position)
%--------------------------------------------------------------------------
% Initialize
VarName = 'Centroid_X';
VarSize = [TIME.saveNum 1];%SPACE.N];
VarType = 1;
VarLoc = 'DATA.FFT_';

eval([VarLoc '.' (VarName) ...
    ' = SaveInitialize(SaveDirectory, VarName, VarSize, VarType, VarLoc);']);
fprintf('. ');



end
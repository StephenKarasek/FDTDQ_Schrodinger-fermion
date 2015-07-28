%% MASTER CONTROL PROGRAM
% Atomic Origin of SHG from FDTD and Spectral Analysis of Wavefunction
%__________________________________________________________________________
% AUTHOR:               Stephen Karasek
%__________________________________________________________________________
% AFFILIATION:          Optical Science Laboratory
%                       (Northeastern University)
%__________________________________________________________________________
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%	VERSION                                             DATE
%..........................................................................
%	{01.00}                                             {2013/12/27}
% - Added Leapfrog Method
%..........................................................................
%	{02.00}                                             {2014/01/08}
% - Added Spectral Analysis
%..........................................................................
%	{03.00}                                             {2014/01/30}
% - Added Dirichlet Boundaries
%..........................................................................
%	{04.00}                                             {2014/02/12}
% - Added Coulomb Potential
%..........................................................................
%	{05.00}                                             {2014/02/19}
% - Added Laser Source (CW)
%..........................................................................
%	{06.00}                                             {2014/02/27}
% - Fixed Laser Source (CEP)
%..........................................................................
%	{07.00}                                             {2014/04/03}
% - Saved Graphs
%..........................................................................
%	{08.00}                                             {2014/05/08}
% - Video Recording
%..........................................................................
%	{09.00}                                             {2014/06/04}
% - Vastly expanded Time & Space solution size (save to disk on mem fill)
%..........................................................................
%	{10.00}                                             {2014/06/23}
% - Save Parameters
%..........................................................................
%	{11.00}                                             {2014/07/02}
% - Bug fixes
%..........................................................................
%	{12.00}                                             {2014/07/15}
% - FDTD-Q fixes
%..........................................................................
%	{13.00}                                             {2014/07/29}
% - Bug fixes
%..........................................................................
%	{14.00}                                             {2014/08/23}
% - QFT Exploration
%..........................................................................
%	{15.00}                                             {2014/09/26}
% - Instantaneous Ionization Rate & Calcs
%..........................................................................
%	{16.00}                                             {2014/09/28}
% - Energy Output
%..........................................................................
%	{17.00}                                             {2014/10/04}
% - Bug fixes
%..........................................................................
%	{18.00}                                             {2014/11/17}
% - Population Level (v0.5), ITP 1st draft
%..........................................................................
%	{19.00}                                             {2014/12/04}
% - Results & Analysis generation (OO)
%..........................................................................
%	{20.00}                                             {2014/12/21}
% - Hydrogen Well fixes, Imaginary Time Propagation (v1.0)
%..........................................................................
%	{21.00}                                             {2015/01/07}
% - General bug fixes, Intensity units correction
%..........................................................................
%	{22.00}                                             {2015/01/29}
% - Fourier Analysis fixes, Population Graph accuracy increased
%..........................................................................
%	{23.00}                                             {2015/03/17}
% - Boundary Conditions fixed, Parallelization (pre-built MATLAB functions)
%..........................................................................
%	{24.00}                                             {2015/04/24}
% - Parallelization (MEX-to-C libraries + CUDA Kernels)
%..........................................................................
%	{FINAL}                                             {2015/07/16}
% - Cleanup for general lab use (i.e. non-authors can now understand it)
%..........................................................................
% 
% 
% 
%	Last Updated:::::::::::::::::::::::::::::::::::::	{2015/04/04}
% 
% 
% 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%__________________________________________________________________________
% VERSION:              Rev.023
%__________________________________________________________________________
% (c)   Optical Science Laboratory, Northeastern University 2013
%       LICENSCED AS (CC BY-NC 3.0), non-commercial.
%       Free to use and modify for non-commercial purposes only.
%__________________________________________________________________________



%% Purge Workspace

close all hidden;
clear all;
clc;




%% Settings
%**************************************************************************
% Settings
%--------------------------------------------------------------------------
% GENERATE_SETTINGS

[ParameterFileList,ParameterFilePath] =...
    uigetfile('../___Simulation_Parameter_Files/',...
    'Select Parameter Files to run','MultiSelect','on');
try 
    ParameterFileList = cellstr(ParameterFileList);
catch
    clc;
    txt_warn = sprintf(['\tCANCELING...\n\n'...
        '\t-->\t\tSimulation will NOT be run!\n\n'...
        '\t-->\t\tTo reattempt operation, please re-run:'...
        '\n\n\t\t\t\t"MASTER_CONTROL_PROGRAM()"\n\n']);
    warning(txt_warn);
    clearvars;
    break;
end



for ParamFileNum=1:length(ParameterFileList)
    clearvars -except ParameterFileList ParameterFilePath ParamFileNum
    close all hidden;
% Load Parameters
load([ParameterFilePath ParameterFileList{ParamFileNum}]);

% TIMING
Runtime_CPU.TOTAL = tic;
Runtime_CPU.Setup = tic;
% Create Directory
DATA_STORAGE(SPACE, TIME, SIM);

% TIMING
disp('CPU::: Setup');
toc(Runtime_CPU.Setup);
%% Solution Space Eigensolvers

%**************************************************************************
% Potential Generator
%--------------------------------------------------------------------------

% POTENTIAL_MAP;

%% Finite Difference Methods

%**************************************************************************
% Potential Generator
%--------------------------------------------------------------------------
% INITIAL_STATE_GENERATOR;
% VSYS_run;

%**************************************************************************
% FDTD-Q
%--------------------------------------------------------------------------
Runtime_CPU.FDM = tic;
FDM_run();
% if DEBUG.SYSTEM.use_GPU
%     FDM_parallel;
% else
%     FDM_serial;
% end
disp('CPU::: Finite Difference Methods');
toc(Runtime_CPU.FDM);




%% Spectral Analysis

%**************************************************************************
% FFT
%--------------------------------------------------------------------------
Runtime_CPU.FFT = tic;
FFT_run();
% if DEBUG.SYSTEM.use_GPU
%     FFT_parallel;
% else
%     FFT_serial;
% end
disp('CPU::: Post-Processing');
toc(Runtime_CPU.FFT);



%% Data Processing & Analysis

%**************************************************************************
% Imaginary Time Propagation (Eigensolver)
%..........................................................................
% Run on all saved effective potentials
% ITP_;
%--------------------------------------------------------------------------


% Spectral Analysis
% SPECTRAL_ANALYSIS;


%% Visualization: Plots & Graphs
if DEBUG.GRAPHS.general; GRAPH_; end

disp('CPU::: TOTAL TIME');
toc(Runtime_CPU.TOTAL);

%% Motion Pictures: Recorded Movies
% Record & Save Video
if DEBUG.VIDEO.movie_ALL
    RECORD_MOVIE;
end

if DEBUG.VIDEO.movie_WaveFunction
    RECORD_MOVIE_WaveFunction;
end

if DEBUG.VIDEO.movie_Potential
    RECORD_MOVIE_potentialwell;
end


%% Save Data IO
DATA_FLAG = true;
if DATA_FLAG
    n=1;
    % PARAMS_file.vars{n} = 'time_dur';n=n+1;
    DATA_file.vars{n} = 'FDM';n=n+1;
    DATA_file.vars{n} = 'FFT';n=n+1;
    DATA_file.vars{n} = 'FSF';n=n+1;
    PARAMS_.vars{n} = 1;
    DATA_.fpath = [SaveDirectory '\' 'DATA_.mat'];
    if exist('SaveDirectory')
        DATA_.fpath = [SaveDirectory '\' 'DATA_.mat'];
    else
        DATA_.fpath = [pwd '\' 'DATA_.mat'];
    end
    
    save(DATA_.fpath, '-v7.3');
    DATA_.IO = matfile(DATA_.fpath, 'Writable', true);
end



end


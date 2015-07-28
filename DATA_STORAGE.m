function [DATA] = DATA_STORAGE(SPACE, TIME, SIM, LASER, DEBUG)
%% DATA STORAGE

% Sanity Check
SaveRoot = '../___DATA';
VisualRoot = '../___VISUALIZATION';


if ~exist(SaveRoot, 'dir')
    mkdir(SaveRoot);
end

if ~exist(VisualRoot, 'dir')
    mkdir(VisualRoot);
end


%% Main

%**************************************************************************
% File Text
%--------------------------------------------------------------------------
switch SIM.FDM.Type
    case 'FDTDQ'
        finite_difference.fname = 'FDM(fdtdq)';
        finite_difference.pname = 'FDTD-Q';
    otherwise
        finite_difference.fname = 'FDM';
        finite_difference.pname = 'FDM';
end


switch SIM.BC.Type
    case 'Mask'
        boundary_conditions.fname = 'Mask';
        boundary_conditions.pname = 'Mask';
    case 'Inf'
        boundary_conditions.fname = 'Inf';
        boundary_conditions.pname = 'Inf';
    otherwise
        boundary_conditions.fname = 'BC';
        boundary_conditions.pname = 'BC';
end
        boundary_conditions.param_fname = ['Rabs=' num2str(SIM.BC.N/SPACE.N, '%03.2f')];
        boundary_conditions.param_pname = ['Rabs=' num2str(SIM.BC.N/SPACE.N, '%03.2f')];

switch SIM.PotentialMap
    case 'square'
        well_type.fname = 'Vc(square)';
        well_type.pname = 'Finite Square Well';
    case 'H1'
        well_type.fname = 'Vc(atom)';
        well_type.pname = 'Parabolic Well';
%     case 'infinite'
%         well_type.fname = 'Vc(inf)'; %InfiniteWell
%         well_type.pname = 'Infinite Square Well';
    otherwise
        well_type.fname = 'Vc(0)'; %No Well
        well_type.pname = 'UndefWell';
end


switch LASER.Type
    case 'off'
        input_type.fname = 'VSrc(OFF)';
        input_type.pname = 'V_L_a_s_e_r(0)';
    case 'CW'
        input_type.fname = sprintf('VSrc(wl=%s,I0=%s)',num2str(LASER.l0),num2str(LASER.I0, '%1.2e'));
        input_type.pname = ['V_L_a_s_e_r(\lambda=' num2str(LASER.l0) ', I_0=' num2str(LASER.I0, '%1.2e') ')'];
    case 'pulse'
        input_type.fname = sprintf('VSrc(CE=%s,Ncyc=%s,I0=%s)',num2str(LASER.EnvelopeWave),num2str(LASER.Ncycles),num2str(LASER.I0, '%1.2e'));
        input_type.pname = ['V_L_a_s_e_r(CE_\lambda=' num2str(LASER.EnvelopeWave) ',N_c_y_c=' num2str(LASER.Ncycles) ', I_0=' num2str(LASER.I0, '%1.2e') ')'];
    otherwise
        input_type.fname = 'VSrc(OFF)';
        input_type.pname = 'V_L_a_s_e_r(0)';
end





%**************************************************************************
% Date/Time String
%--------------------------------------------------------------------------
% SaveRoot = uigetdir(['../___DATA'],'Please determine root filepath');


DateTime_str = clock;

if DEBUG.DetailedFileName
NAME_txt = sprintf('(%02.0f-%02.0f-%02.0f---%02.0f-%02.0f-%02.0f)___Data[{%sD}, {t=%s},{r=%s},{%s},{%s},{%s},{%s},(%s)]',...
    DateTime_str(1), DateTime_str(2), DateTime_str(3),...
    DateTime_str(4), DateTime_str(5), round(DateTime_str(6)),...
    num2str(SPACE.Dims, '%1.0f'), num2str(TIME.N, '%3.2e'), num2str(SPACE.length, '%3.2e'),...
    finite_difference.fname, boundary_conditions.fname,...
    boundary_conditions.param_fname, well_type.fname, input_type.fname);
else
    ShortNameNum = num2str(length(dir(fullfile(SaveRoot)))-2,'%04.0f');
    NAME_txt = ['DataStorage_Run(' num2str(ShortNameNum) ')'];
end
SaveDirectory = [SaveRoot '\' NAME_txt];
            
result = mkdir(SaveDirectory);



%%
%**************************************************************************
% NULL Dataset
%--------------------------------------------------------------------------
VarName = 'NULL';
VarSize = [1 1];
VarType = 1;
VarLoc = 'DATA.NULL';

eval([VarLoc '.' (VarName) ...
    ' = SaveInitialize(SaveDirectory, VarName, VarSize, VarType, VarLoc);']);



%% SAVE SETTINGS


% NAME_txt
ShortNameNum = num2str(length(dir(fullfile(DEBUG.prefix,'*.mat'))),'%04.0f');

n=1;
% PARAMS_file.vars{n} = 'time_dur';n=n+1;
DATA.PARAMS_file.vars{n} = 'SPACE';n=n+1;
DATA.PARAMS_file.vars{n} = 'TIME';n=n+1;
DATA.PARAMS_file.vars{n} = 'LASER';n=n+1;
DATA.PARAMS_.vars{n} = 1;
DATA.PARAMS_.IO = cell(1);
if exist('SaveDirectory')
    DATA.PARAMS_.saveDir = SaveDirectory;
    DATA.PARAMS_.saveRoot = SaveRoot;
    DATA.PARAMS_.saveVisual = VisualRoot;
    DATA.PARAMS_.saveName = NAME_txt;
    DATA.PARAMS_.fpath = [SaveDirectory '\' 'PARAMS_.mat'];
else
    DATA.PARAMS_.fpath = [pwd '\' 'PARAMS_.mat'];
end



save(DATA.PARAMS_.fpath, '-v7.3');
DATA.PARAMS_.IO = matfile(DATA.PARAMS_.fpath, 'Writable', true);


end

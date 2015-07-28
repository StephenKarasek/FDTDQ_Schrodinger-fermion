function [DATA] = FDM_run(SPACE, TIME, SIM, PULSE, ITP, LASER, DATA, DEBUG)
% function [varargout] = FDM_run(varargin)
%% FDTD-Q
CONSTANTS;

%% Data Definitions
fprintf('\n\n\t--> BUILDING DATA STRUCTURE... ');

% evalin('caller',FDM_vars(SPACE, TIME, PULSE, ITP, PARAMS_, DEBUG));
[DATA] = FDM_vars(SPACE, TIME, SIM, PULSE, ITP, LASER, DATA, DEBUG);

fprintf('\t\t... DONE!');



%% Concurrent Variables
clear FDTD_Save_Mem FDTD_Load_Mem FDTD_Total_Mem;


%**************************************************************************
% SAVE
% 
% 1)    Calculate Observable
% 
% 2)    Save Observable to File
% 
%--------------------------------------------------------------------------
SaveOrigin = 'DATA.FDTD_';
SaveList = fieldnames(eval(SaveOrigin));

% SAVE Function Handle
grabSaveVar = @(VAR) [SaveOrigin '.' (VAR)];

% List of Data Objects
for n=1:length(SaveList)
    FDTD_Save_Mem{n} = eval(grabSaveVar(SaveList{n}));
end

%
%**************************************************************************
% SAVE
% 
% 1)    Load Observable from File
% 
% 2)    Use Observable to calculate related Observable
% 
%--------------------------------------------------------------------------
% if ITP.Function
%     LoadOrigin = 'DATA.VPOT_';
%     SaveList = fieldnames(eval(LoadOrigin));
%     LoadList = {'EigenFunction'; 'EigenCoefficient'; 'EigenEnergy'; 'EigenSolution'};
%     LoadList = '';
% else
    LoadOrigin = 'DATA.NULL';
    LoadList = '';

% LOAD Function Handle
grabLoadVar = @(VAR) [LoadOrigin '.' (VAR)];

% List of Data Objects
FDTD_Load_Mem = [];




%**************************************************************************
% Total Memory Usage
%--------------------------------------------------------------------------
FDTD_Total_Mem = [FDTD_Save_Mem, FDTD_Load_Mem];



%% MEMORY OPTIMIZATION
% evalstr = ['FDTDQ_file.IO{n}.' FDTDQ_file.vars{n} ' = ' FDTDQ_file.vars{n} ';'];
%     eval(evalstr);
use_GPU = false;

[   DATA.FDTD_MEM.MaxTimeSteps, DATA.FDTD_MEM.PerTimeStep,...
    DATA.FDTD_MEM.SizeToAllocate, DATA.FDTD_MEM.MaxSize_Allocate...
] = ...
    MemAnalysis(use_GPU,FDTD_Total_Mem{:});

%                         VarName,...
%                         VarType,...
%                         VarSize...
% );
% [DATA.FDTD_; DATA.Vpot_]
%                     use_GPU, );
                    
                    
if TIME.saveNum>=DATA.FDTD_MEM.MaxSize_Allocate %.save(1)
    warnstr = sprintf(['\t\tMEMORY OVERFLOW\n'...
        '\n\tVariable sizes are too large for memory allocation .'...
        '\n\tPlease reduce the number of space steps, timesteps, or save'...
        'interval, and try again.']);
    warning(warnstr);
    return;
end


%% MEMORY MANAGEMENT

%**************************************************************************
% Memory Index
%--------------------------------------------------------------------------
%   Calculate index for arrays to be kept in memory
%..........................................................................
fprintf('\n\n\t--> CREATING MEMORY INDEX...\t\t:');
[idx] = createMemIndex(DATA.FDTD_MEM.MaxTimeSteps, TIME.saveNum);
fprintf('\t\t... DONE!\n\n:');


%**************************************************************************
% Initialize Saved Arrays
%--------------------------------------------------------------------------
%   Calculate index for arrays to be kept in memory
%..........................................................................
fprintf('\n\n\t--> PRE-ALLOCATING MEMORY (saves)... ');
for n=1:length(FDTD_Save_Mem)
    [RESULTS] = VarInitialize(idx,FDTD_Save_Mem{n});
end
fprintf('\t\t... DONE!\n\n');


%**************************************************************************
% Load Arrays
%--------------------------------------------------------------------------
%   Calculate index for arrays to be kept in memory
%..........................................................................
fprintf('\n\n\t--> PRE-ALLOCATING MEMORY (loads)... ');
for n=1:length(FDTD_Load_Mem)
    [RESULTS] = LoadFormat(idx,FDTD_Load_Mem{n});
end
fprintf('\t\t... DONE!\n\n');


%% Potential Energy Map

DATA = MaxPotential(SPACE, TIME, SIM, PULSE, ITP, LASER, DATA, DEBUG);

V_sys = DATA.V_.Vsys;
V_max = DATA.V_.Vmax;
V_min = DATA.V_.Vmin;



%% Initial wavefunction

% Varying over Time AND Space
%...

% Varying over Time
%...


if ITP.Function%&&false
    [EigenFunction, EigenCoefficient, EigenEnergy, ImagTimeProp, EigenSolution] =...
    ITP_EigenFinder(SPACE,TIME,SIM,ITP,V_sys,1);

% [EigenFunction, EigenCoefficient, EigenEnergy, ImagTimeProp, EigenSolution] =...
%     ITP_EigenFinder(SPACE,TIME,BCs, SIM,ITP,V_sys,1);


% psi = Psi0_QHO(PULSE,SIM,SPACE,TIME, EigenFunction);    
    psi0 = Psi0_QHO(PULSE,SIM,SPACE,TIME, EigenFunction, EigenCoefficient, EigenEnergy, EigenSolution);    
%     psi = Psi0_QHO(PULSE,SPACE,TIME, Psi0_Norm, Psi0_coeff);    
else
    psi0 = Psi0_QHO(PULSE,SIM,SPACE,TIME);
end


% Local (space-variant & single-point)
R_current = zeros(1,SPACE.N); 
I_current = R_current;
R_next = R_current;
I_next = R_current;

% Propagation Constants
% s1=(hPlanck*TIME.delta)/(me*SPACE.delta^2);
% s2=(TIME.delta)/(hPlanck);
% s=(TIME.delta)/(2*SPACE.delta^2);



%% Setup
% Avg.Max Value of Psi
WaveMax = sqrt(sum(abs(psi0)/2)/SPACE.N);
PdfMax = (WaveMax/2)^2;

% Extract the real and imaginary parts of the wavefunction
% Initialise current real and imaginary parts of psi
R_current=real(psi0);
I_current=imag(psi0);


% BCs
switch SIM.BC.Type
    
    case 'Inf'
R_current(1)=0;%R_next(2);
R_current(SPACE.N)=0;%R_next(N-1);
I_current(1)=0;%R_next(2);
I_current(SPACE.N)=0;%R_next(N-1);

    case 'Mask'
R_current(SIM.BC.R{1})=0;%R_next(2);
R_current(SIM.BC.R{2})=0;%R_next(N-1);
I_current(SIM.BC.R{1})=0;%R_next(2);
I_current(SIM.BC.R{2})=0;%R_next(N-1);

end

% 
% [R_current, I_current] =...
%     FiniteDifferenceMethod(V_eff, SPACE.delta, TIME.delta, DEBUG.FDM, DEBUG.BC,...
%     R_current, I_current, 1);

% t=t+TIME.delta/2;
% I_next(2:SPACE.N-1) = I_current(2:SPACE.N-1) ...
%     + s1*(R_current(3:SPACE.N)-2*R_current(2:SPACE.N-1)+R_current(1:SPACE.N-2))...
%     - s2*V_eff(2:SPACE.N-1).*R_current(2:SPACE.N-1);


% Timeshift
I_next = I_current;
R_next = R_current;



%% Begin FDTD-Q
%##########################################################################
for t = 1:TIME.N;

    
%**************************************************************************
% --> Potential Tilting
% Effects on potential from Laser source
%-------------------------------------------------------------------------- 
% [V_laser] =
% PotentialLaser(LASER, SPACE, TIME, t);
[Laser_E] = PotentialLaser(LASER, SPACE, TIME, t);

% E =  -E0*(sin(pi*t/Laser.T)^2*sin(Laser.w0*t) ...
%     - (pi/(Laser.w0*Laser.T))*sin(cos(Laser.w0*t)*2*pi*t/Laser.T));

switch LASER.GaugeCoupling

    case 'coulomb'
         V_ext = e0*SPACE.Axis_.*Laser_E;
%          V_ext = e0*SPACE.Axis.*Laser_E;
%         Laser_coupling
        
    case 'lorenz'
        V_ext = Laser_E*0;
%     evalin('caller',['V_laser = 0*Laser_E;']);

    otherwise
        V_ext = 0;
end


% E = E0*sin(Laser.w0*t*TIME.delta);
% V_laser = x.*E*e0;
%     V_eff2 = V_ext + V_eff;
    V_eff = V_ext + V_sys;

    


%**************************************************************************
% --> Finite Difference Scheme
%-------------------------------------------------------------------------- 
% Update Equations
[R_next, I_next] = FDM_fdtdq_1D(SPACE, TIME, SIM, PULSE, R_next, I_next, V_eff, R_current, I_current);
% [RESULTS] = FDM_fdtdq_1D_Generate(SPACE, TIME, SIM, PULSE, V_eff, R_next, I_next, R_current, I_current);


%**************************************************************************
% Boundary Conditions
%--------------------------------------------------------------------------
[R_next, I_next] = BC_fdtdq_1D(SPACE, TIME, SIM, PULSE, R_next, I_next, V_eff);
% [RESULTS] = BC_fdtdq_Generate(SPACE, TIME, SIM, PULSE, R_next, I_next, R_current, I_current);


%**************************************************************************
% Update Equations
%--------------------------------------------------------------------------
R_current = R_next;
I_current = I_next;


%//////////////////////////////////////////////////////////////////////////
if rem(t, TIME.save)== 0;
    time_step = t/TIME.save;
    
    
%**************************************************************************
% OBSERVABLES
%--------------------------------------------------------------------------

% Position Wavefunction
Psi_REAL(idx.count,:) = R_next;
Psi_IMAG(idx.count,:) = I_next;

PDF_pos(idx.count,:) = ...
    (conj(Psi_REAL(idx.count,:)+Psi_IMAG(idx.count,:)*1i)...
    .*(Psi_REAL(idx.count,:)+Psi_IMAG(idx.count,:)*1i));
prob_density = PDF_pos(idx.count,:);
% R_next.*R_current+I_next.*I_current;

% Potential
if DEBUG.Observables_R.Potential_Energy
    V_pot(idx.count,:) = V_eff;
%     V_pot2(idx.count,:) = V_eff2;
    V_laser(idx.count,:) = V_ext;
    
%     T_kwave(idx.count,:) = Kwave_;
%     T_kin(idx.count,:) = Kin_;
    E_laser(idx.count) = Laser_E;
end


% Centroid
% centroid(t,:) = sum(Psi_total.*conj(Psi_total).*x);

% KE(t) = 

% PE(t) = 

%**************************************************************************
% ANIMATE
%--------------------------------------------------------------------------
% if false
if (DEBUG.VIDEO.animate)&&(rem(t, TIME.animate)==0)
    
% Psi(R+I), PDF, Vc
    if strcmpi(DEBUG.VIDEO.animate_,'all')
fig0 = animateWavefunction(R_current, I_current, prob_density, V_eff, V_min, V_max, WaveMax,PdfMax,...
    t, TIME.delta, SPACE);
pause(0.001); 
    end
    
% Psi(R+I), PDF
    if strcmpi(DEBUG.VIDEO.animate_,'wavefunction')
fig0 = animateWavefunction_NO_WELL(R_current, I_current, prob_density, V_eff,...
    t, TIME.delta, SPACE);
pause(0.001); 
        
    end
    
% Vc
    if strcmpi(DEBUG.VIDEO.animate_,'potential')
fig0 = animatePotentialWell(V_eff,t, TIME.delta, SPACE);
pause(0.001); 
    end
          

% end
end  
    
    

if idx.count == DATA.FDTD_MEM.MaxTimeSteps
    
    
%**************************************************************************
% SAVE TO DISK
%--------------------------------------------------------------------------

    fprintf('\n\n\t\t\t***COPYING DATA TO DISK***\n\n:'); tic;

for n=1:length(FDTD_Save_Mem)
    
    fprintf('\n%s\t:', getfield(FDTD_Save_Mem{n},'VarName')); tic;
    
    [RESULTS] = SaveFormat(idx, FDTD_Save_Mem{n});

%     eval(SaveString);
    
    fprintf('\t\t%f\tSeconds\n', toc);
end

% fprintf('\n%s\t:', DATA.vars{1}); tic;
% DATA.IO{1}.Psi_REAL(idx.START:idx.END,1:N) = ...
%     Psi_REAL(1:idx.length,1:N);
% fprintf('\t\t%f\tSeconds\n', toc);
%     
% fprintf('\n%s\t:', DATA.vars{2}); tic;
% DATA.IO{2}.Psi_IMAG(idx.START:idx.END,1:N) = ...
%     Psi_IMAG(1:idx.length,1:N);
% fprintf('\t\t%f\tSeconds\n', toc);


% fprintf('\t\t%f\tSeconds\n', toc);

%     for n=1:length(FDTDQ_file.vars)
% fprintf('\n%s\t:', FDTDQ_file.vars{n});
%     switch Dims
%         case 1
%     evalstr = ['FDTDQ_file.IO{n}.' ...
%         FDTDQ_file.vars{n} '(' num2str(idx.START) ':' num2str(idx.END) ','...
%          num2str(1) ':' num2str(FDTDQ_file.info{n}(2)) ') = ' ...
%         FDTDQ_file.vars{n} '(' num2str(1) ':' num2str(idx.length) ','...
%         num2str(1) ':' num2str(FDTDQ_file.info{n}(2)) ');'];
    
%         case 2
            
%         case 3
            
%     end
%     eval(evalstr);
    
% fprintf('\t\t%f\tSeconds\n', toc);
%     end
    
    
%**************************************************************************
% CLEAR MEMORY
%--------------------------------------------------------------------------
fprintf('\n\n\t\t\t***CLEARING MEMORY***\n\n:');% tic;

for n=1:length(FDTD_Total_Mem)

    % Requires variable name search from earlier
    [RESULTS] = VarInitialize(idx,FDTD_Total_Mem{n});
    
    % Run initialization command
%     eval(InitializeString);
end

% Psi_TOTAL(1:idx.length, 1:N) = zeros(1, 'like', [1+1i]);
% Psi_TOTAL(1:idx.length, 1:N) = zeros(1);
% Psi_TOTAL(1:idx.length, 1:N) = zeros(1);
% V_pot(1:idx.length, 1:N) = zeros(1);
    
    
    
    
%     for n=1:length(FDTDQ_file.vars)
        
%     evalstr = [FDTDQ_file.vars{n} '(' num2str(1) ':' num2str(idx.length) ...
%         ',:) = zeros(' num2str(1) ');'];
        
    % Initialize Variable
% switch Dims
    
%     case 1
%     if strcmpi(FDTDQ_file.type{n}, 'real')
%     evalstr = [FDTDQ_file.vars{n} '(' num2str(1) ':' num2str(idx.length) ','...
%         num2str(1) ':' num2str(FDTDQ_file.info{n}(2)) ')'...
%         ' = zeros(' num2str(1) ');'];
%     else
%     evalstr = [FDTDQ_file.vars{n} '(' num2str(1) ':' num2str(idx.length) ','...
%          num2str(1) ':' num2str(FDTDQ_file.info{n}(2)) ')'...
%          ' = zeros(' num2str(1) ', ''like'', [1+1i]);'];
%     end
%     case 2
%     case 3
%   end
%     eval(evalstr);







%**************************************************************************
% RESET INDEX
%--------------------------------------------------------------------------
    idx.START = idx.START+idx.count;
    idx.END = idx.END+idx.count;
    idx.count = 0;
    
    if idx.END>TIME.saveNum; idx.END = TIME.saveNum; end
    
    
    
    end

% Iterate Memory Index
    idx.count = idx.count+1;
    
    
end


show_t = (TIME.N/TIME.save);
if show_t>1e6; show_t = 1e6; end
if show_t<1e4; show_t = 1e4; end
    if ~mod(t, show_t);

%     switch Dims
%         case 1
    fprintf('\tFDM(R-space)\t:::\t\t[%d\tof\t%d]\t=\t<%d*(s)>\n', t, TIME.N, t*TIME.delta);
% fprintf('\tPosition(R-space)\t:::\t\t[%d\tof\t%d])\t=\t<%d*(s)>\n',...
%     time_step*t_save, T_len, time_step*t_save*TIME.delta);
    
    %         
%         case 2
% 	fprintf(['\tFDTD-Q(position)\t:::\t\t[%d\tof\t%d]\t=\t<%d*(s)>\n',...
%         t, T_len, t*TIME.delta]);
%         
%         case 3
%     fprintf(['\tFDTD-Q(position)\t:::\t\t[%d\tof\t%d]\t=\t<%d*(s)>\n',...
%         t, T_len, t*TIME.delta]);
%     
%         otherwise
%     end
% format? '%10.4e\n'
    end


end


if (DEBUG.VIDEO.animate)
    close(fig0);
end

% Come up with a better method to make sure memory index stepping isn't
% stupid.
idx.count = idx.count-1;



%% Save




%**************************************************************************
% SAVE DATA
%--------------------------------------------------------------------------
fprintf('\nSaving FDM (Position Space) Data...\n');
% tic;    


for n=1:length(FDTD_Save_Mem)
    
    fprintf('\n%s\t:', getfield(FDTD_Save_Mem{n},'VarName')); tic;
    
    [RESULTS] = SaveFormat(idx,FDTD_Save_Mem{n});

%     eval(SaveString);
    
    fprintf('\t\t%f\tSeconds\n', toc);
end




%% Cleanup
%**************************************************************************
% CLEAR MEMORY
%--------------------------------------------------------------------------
fprintf('\n\n\t\t\t***CLEARING MEMORY***\n\n:');% tic;

for n=1:length(FDTD_Total_Mem)

    % Requires variable name search from earlier
%     [RESULTS] = VarInitialize(idx,FDTD_Total_Mem{n});
[RESULTS] = CleanUpVars(FDTD_Total_Mem{n});
    
    % Run initialization command
%     eval(InitializeString);
end



fprintf('___________________________________________________________________________\n\n\n');

                        


%% Save PARAMS_FDM_
end


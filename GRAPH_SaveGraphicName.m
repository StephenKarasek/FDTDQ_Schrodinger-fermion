function [TEXT_] = GRAPH_SaveGraphicName(SPACE, TIME, SIM, PULSE, ITP, LASER, DATA, DEBUG)
%% Determine Save Graphic Names
CONSTANTS
% SaveDirectory = DATA.PARAMS_.saveDir;
% load([SaveDirectory '/PARAMS_'], 'PARAMS_');

%% Directory
TEXT_.SaveDirectory = [DATA.PARAMS_.saveDir '\'];
TEXT_.saveVisual = [DATA.PARAMS_.saveVisual];

TEXT_.txt_misc = sprintf([''...
    'FFT0padI(' num2str(SPACE.ZeroPad_Factor) '),'...
    'BC(' SIM.BC.Type ',Nabs=' num2str(SIM.BC.Nratio, '%02.2f') ')']);

TEXT_.txt_pot = ['V(' SIM.PotentialMap '), Creg(' num2str(PULSE.Creg, '%1.e')];
% TEXT_.txt_pot = ['.          V(' SIM.PotentialMap ')'];

%% Pulse Type
switch PULSE.Type
    
%**************************************************************************
% Eigenstate
%--------------------------------------------------------------------------
    case {'Cosine'; 'Cos'; 'Sin'; 'Sine'}
%         TEXT_.txt_pulse = ['_.     \Psi(E:{' num2str(PULSE.EnergyState) '}'];
        TEXT_.txt_pulse = ['\Psi_0(E_n=' num2str(PULSE.EnergyState) ')'];
        
        TEXT_.Estate_ = ['Esin(n=' num2str(PULSE.EnergyState(1)) ];
        if length(PULSE.EnergyState)>1
            for N=2:length(PULSE.EnergyState)
            TEXT_.Estate_ = [TEXT_.Estate_ ',' num2str(PULSE.EnergyState(N)) ];
            end
        end
        TEXT_.Estate_ = [TEXT_.Estate_ ')'];
        
        

%**************************************************************************
% Gaussian
%--------------------------------------------------------------------------
    case 'Gaussian'
%         TEXT_.txt_pulse = ['_.     \Psi(\sigma_' num2str(PULSE.EnergyState) '_-_5)'];
%         TEXT_.txt_pulse = ['\Psi(\sigma_' num2str(PULSE.EnergyState) '_-_5)'];
%         TEXT_.Estate_ = ['Gauss[x_=' num2str(PULSE.Width,'%1.1e') ...
%             ', k_=' num2str(PULSE.Momentum,'%1.1e') ,']'];        

TEXT_.txt_pulse = ['\Psi(E_n=\sigma_1_:_n)'];
TEXT_.Estate_ = ['Eg(Xw=' num2str(PULSE.Width,'%1.1e') ')'];
%..........................................................................


%**************************************************************************
% Eigenstate
%--------------------------------------------------------------------------
    case {'EigenState'; 'ITP'; 'ImagTimeProp'; 'ImaginaryTimePropagation'}
%         TEXT_.txt_pulse = ['_.     \Psi(\sigma_' num2str(PULSE.EnergyState) '_-_5)'];
%         TEXT_.txt_pulse = ['\Psi(\sigma_' num2str(PULSE.EnergyState) '_-_5)'];
        TEXT_.txt_pulse = ['\Psi(E_n_''=' num2str(PULSE.EnergyState) ')'];
        
        TEXT_.Estate_ = ['Eeig(n=' num2str(PULSE.EnergyState(1)) ];
        if length(PULSE.EnergyState)>1
            for N=2:length(PULSE.EnergyState)
            TEXT_.Estate_ = [TEXT_.Estate_ ',' num2str(PULSE.EnergyState(N)) ];
            end
        end
        TEXT_.Estate_ = [TEXT_.Estate_ ')'];
%..........................................................................
        


%**************************************************************************
% Exact
%--------------------------------------------------------------------------
    case {'Analytical'; 'Exact'}
%         TEXT_.txt_pulse = ['_.     \Psi(\sigma_' num2str(PULSE.EnergyState) '_-_5)'];
%         TEXT_.txt_pulse = ['\Psi(\sigma_' num2str(PULSE.EnergyState) '_-_5)'];
        TEXT_.txt_pulse = ['\Psi(E_n_''=' num2str(PULSE.EnergyState) ')'];
        
        TEXT_.Estate_ = ['Eext(n=' num2str(PULSE.EnergyState(1)) ];
        if length(PULSE.EnergyState)>1
            for N=2:length(PULSE.EnergyState)
            TEXT_.Estate_ = [TEXT_.Estate_ ',' num2str(PULSE.EnergyState(N)) ];
            end
        end
        TEXT_.Estate_ = [TEXT_.Estate_ ')'];
%..........................................................................
        

%**************************************************************************
% ...
%--------------------------------------------------------------------------
    otherwise
%         TEXT_.txt_pulse = ['_.          \Psi(\sigma_' num2str(PULSE.EnergyState) '_-_5)'];
%         TEXT_.txt_pulse = ['\Psi(\sigma_' num2str(PULSE.EnergyState) '_-_5)'];
        TEXT_.txt_pulse = ['\Psi(E_n=...)'];
        TEXT_.Estate_ = ['Psi0(...,En=' num2str(PULSE.EnergyState(:)) ];
end

TEXT_.StatePopCalc_ = ['StateCalc(' num2str(PULSE.NumCalcStates, '%01.f') ')'];

% TEXT_.txt_EnPot = {TEXT_.txt_pulse, TEXT_.txt_pot};
TEXT_.txt_EnPot = [TEXT_.txt_pulse ', ' TEXT_.txt_pot];

TEXT_.V0_ = ['V(' SIM.PotentialMap ')'];
TEXT_.EnergyParams = [TEXT_.V0_ ',' TEXT_.Estate_ ',' TEXT_.StatePopCalc_];


%% Space & Time
%//////////////////////////////////////////////////////////////////////////
TEXT_.R0_ = '';
for i=1:length(PULSE.InitialPos)
%     component_ = fprintf('  %s  ', num2str(PULSE.InitialPos(i)))
    TEXT_.R0_ = [TEXT_.R0_ ' ' num2str(PULSE.InitialPos(i))];
end

TEXT_.R0_ = sprintf(['[X(' num2str(SPACE.Dims) 'd,'...
    'x=' num2str(SPACE.length/a0) 'a0,dx=' num2str(SPACE.N) '),'...
    'T(t=' num2str(TIME.length, '%3.2e') 's,dt=' num2str(TIME.N, '%2.1e') ')]']);



%% Laser
%//////////////////////////////////////////////////////////////////////////
switch LASER.Type
    case 'CW'
        LaserInfo_ = ['[' LASER.Type ',' num2str(LASER.l0, '%2.1e') 'm' ...
            'N=' num2str(LASER.nPhoton, '%2.1f'), ',Io=' num2str(LASER.I0, '%2.1e') ']'];
        
    case 'pulse'
        LaserInfo_ = ['[' LASER.Type ',' num2str(LASER.EnvelopeWave, '%2.1e') ','...
            'N=' num2str(LASER.Ncycles, '%1.1g') ',I0=' num2str(LASER.I0, '%2.1e') ']'];
                
    otherwise
        LaserInfo_ = ['NoLaser'];
                
end


I0_ = ['I-' num2str(LASER.I0, '%1.1e')];    

% TEXT_.SaveGraphicName = [LaserInfo_ ',' TEXT_.Estate_ ',' TEXT_.R0_ ' .png'];
% saveas(hfig, TEXT_.SaveGraphicName);
% SaveDirectory



% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Laser Source
% for t=1:TIME.N
% PotentialLaser(LASER.GaugeCoupling, LASER, SPACE, TIME, t, 1);
% end
% if false
% Laser_E = zeros(1,TIME.N/TIME.save);
% Laser_E = zeros(1,TIME.N/TIME.save);



switch LASER.Type
%//////////////////////////////////////////////////////////////////////////
    case {'pulse','CEP', 'chirped-pulse', 'Chirped-Pulse', 'CPA'}
% for t=TIME.save:TIME.save:TIME.N
    
%     Laser_E(t/TIME.save) = PotentialLaser(LASER, SPACE, TIME, t);
%     Laser_E(t/TIME.save) = LASER.E0*sin(LASER.w0*t*TIME.delta/(2*LASER.Ncycles))^2 ...
%     *cos(LASER.w0*t*TIME.delta+LASER.CEP);
    
% end
        TEXT_.txt_Ltype = [LASER.Type '(\lambda_C_E_P=' num2str(LASER.EnvelopeWave, '%2.1e')...
            ', N_c_y_c=' num2str(LASER.Ncycles) ')'];
%                 TEXT_.txt_Ltype = ['_.     ' LASER.Type '(\lambda_C_E_P=' num2str(LASER.EnvelopeWave, '%2.1e')...
%             ', N_c_y_c=' num2str(LASER.Ncycles) ')'];
        TEXT_.txt_Lint = ['Laser:(I_0= ' num2str(LASER.I0, '%2.1e') ')'];
%     ttxt = { TEXT_.txt_Lint, TEXT_.txt_Ltype};
TEXT_.txt_Laser = [TEXT_.txt_Lint ', ' TEXT_.txt_Ltype];
%//////////////////////////////////////////////////////////////////////////
    case {'CW', 'Continuous-Wave', 'continuous', 'continuous-wave'}
%         for t=TIME.save:TIME.save:TIME.N
    
%             Laser_E(t/TIME.save) = PotentialLaser(LASER, SPACE, TIME, t);
            
%             Laser_E(t/TIME.save) = LASER.E0*cos(LASER.w0*t*TIME.delta+LASER.Phase);
%         end
        TEXT_.txt_Ltype = [LASER.Type '(\lambda= ' num2str(LASER.l0, '%2.1e') ')'];
%         TEXT_.txt_Ltype = ['_.     ' LASER.Type '(\lambda= ' num2str(LASER.l0, '%2.1e') ')'];
        TEXT_.txt_Lint = ['Laser:(I_0= ' num2str(LASER.I0, '%2.1e') ')'];
%     ttxt = { TEXT_.txt_Lint, TEXT_.txt_Ltype};
   TEXT_.txt_Laser = [TEXT_.txt_Lint ', ' TEXT_.txt_Ltype];
        
    otherwise
%//////////////////////////////////////////////////////////////////////////
end




% Pulses
if (LASER.PERIODcutoff~=0)
    TEXT_.pulse_gt = ...
    ['\phi_o_f_f=' num2str(2*LASER.PERIODcutoff/LASER.T, '%2.1f') '\pi,inf*pulses'];
%     TEXT_.pulse_txt = [num2str(2*LASER.PERIODcutoff/LASER.T, '%2.1f') '.pi off,inf.pulses'];
TEXT_.pulse_txt = '';
end

if (LASER.FINALcutoff~=0)
    TEXT_.pulse_gt = ['\phi_o_f_f=0\pi,'...
    num2str(LASER.FINALcutoff/LASER.T, '%2.1f') '*pulses'];
%     TEXT_.pulse_txt = ['0.pi off,'num2str(LASER.FINALcutoff/LASER.T, '%2.1f') '.pulses'];
TEXT_.pulse_txt = '';
end

if ((LASER.PERIODcutoff==0)&&(LASER.FINALcutoff==0))
    TEXT_.pulse_gt = ['\phi_o_f_f=0\pi,inf*pulses'];
%     TEXT_.pulse_txt = ['0.pi off,inf.pulses'];
TEXT_.pulse_txt = '';
end

TEXT_.Laser_txt = [LaserInfo_ '.' TEXT_.pulse_txt];

TEXT_.SaveGraphicName = [TEXT_.R0_  ',' TEXT_.EnergyParams ',' TEXT_.Laser_txt ',' TEXT_.txt_misc '.png'];





end
function [DATA] = FFT_run(SPACE, TIME, SIM, PULSE, ITP, LASER, DATA, DEBUG)
%% FOURIER_TRANS

% Must be run AFTER FDTD_Q
CONSTANTS;

%% Save

fprintf('\n\n\t--> BUILDING DATA STRUCTURE... ');

DATA = FFT_vars(SPACE, TIME, SIM, PULSE, ITP, LASER, DATA, DEBUG);

fprintf('\t\t... DONE!');

%% Concurrent Variables
clear FFT_Save_Mem FFT_Load_Mem FFT_Total_Mem;


%**************************************************************************
% SAVE
% 
% 1)    Calculate Observable
% 
% 2)    Save Observable to File
% 
%--------------------------------------------------------------------------
SaveOrigin = 'DATA.FFT_';
SaveList = fieldnames(eval(SaveOrigin));

% SAVE Function Handle
grabSaveVar = @(VAR) [SaveOrigin '.' (VAR)];

% List of Data Objects
for n=1:length(SaveList)
    FFT_Save_Mem{n} = eval(grabSaveVar(SaveList{n}));
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
LoadOrigin = 'DATA.FDTD_';

LoadList = {'Psi_REAL'; 'Psi_IMAG'; 'V_pot'; 'V_laser'; 'V_atom'};

% LOAD Function Handle
grabLoadVar = @(VAR) [LoadOrigin '.' (VAR)];

% List of Data Objects
for n=1:length(LoadList)
    FFT_Load_Mem{n} = eval(grabLoadVar(LoadList{n}));
end


%**************************************************************************
% Total Memory Usage
%--------------------------------------------------------------------------
FFT_Total_Mem = [FFT_Save_Mem, FFT_Load_Mem];




%% MEMORY OPTIMIZATION

use_GPU = false;

[   DATA.FFT_MEM.MaxTimeSteps, DATA.FFT_MEM.PerTimeStep,...
    DATA.FFT_MEM.SizeToAllocate, DATA.FFT_MEM.MaxSize_Allocate...
] = ...
    MemAnalysis(use_GPU,FFT_Total_Mem{:});


%**************************************************************************
% Memory Check
%--------------------------------------------------------------------------
% Is there enough RAM to run [N>=1] spatial step before saving to HDD?
% YES? --> Proceed. NO? --> WARNING(halt)
%
% *** For 1D and reasonably sized 2D spaces this isn't a problem
% *** 3D spaces will provide a greater challenge
%> EXAMPLE... [X,Y,Z]: (2^10) * 8 bits * 2 (real + imag) == 16 GB.
%> Therefore, after EACH spatial step we want to save, all data in memory
% must be saved to disk, cleared, and reloaded to proceed to the next step.
% One simulation of 10^7 time steps would take more than 10 hours...

if TIME.saveNum>=DATA.FFT_MEM.MaxSize_Allocate
    warnstr = sprintf(['\t\tMEMORY OVERFLOW\n'...
        '\n\tVariable sizes are too large for memory allocation.'...
        '\n\tPlease reduce the number of space steps, timesteps, or save'...
        'interval, and try again.']);
    warning(warnstr);
    return;
end



%==========================================================================
%% Memory Index
%..........................................................................

[idx] = createMemIndex(DATA.FFT_MEM.MaxTimeSteps, TIME.saveNum);



%% Concurrent Data Arrays
%**************************************************************************
% Initialize Saved Arrays
%--------------------------------------------------------------------------
fprintf('\n\n\t--> Pre-allocating memory (saves)...\t\t:\n\n');

for n=1:length(FFT_Save_Mem)
    [RESULTS] = VarInitialize(idx,FFT_Save_Mem{n});
end


%**************************************************************************
% Load Vars
%--------------------------------------------------------------------------
fprintf('\n\n\t--> Pre-allocating memory (loads)...\t\t:\n\n');

for n=1:length(FFT_Load_Mem)
    [RESULTS] = LoadFormat(idx,FFT_Load_Mem{n});
end


%==========================================================================
%% Fourier Transform
%..........................................................................


%**************************************************************************
% FFT
%--------------------------------------------------------------------------
% tic
for t=1:TIME.saveNum
    
% FFT --> freq; divide by ...   sqrt(2*pi*hPlanck) = sqrt(hPlanck_full)
% fft_TOTAL = fft((Psi_REAL(idx.count,:)+1i*Psi_IMAG(idx.count,:)), SPACE.N0);
%     *(1/(sqrt(2*pi)*hPlanck));

SignalRange = round(SPACE.N0/2)+[(-(SPACE.N0-SPACE.N)/2+1) ((SPACE.N0-SPACE.N)/2)];

fft_REAL_ = real(fft((Psi_REAL(idx.count,:)+1i*Psi_IMAG(idx.count,:)), SPACE.N0));
fft_REAL(idx.count,:) = real(fft_REAL_);
% fft_REAL(idx.count,:) = real(fft_REAL_(SignalRange(1):SignalRange(2)));

fft_IMAG_ = imag(fft((Psi_REAL(idx.count,:)+1i*Psi_IMAG(idx.count,:)), SPACE.N0));
fft_IMAG(idx.count,:) = imag(fft_IMAG_);
% fft_IMAG(idx.count,:) = imag(fft_IMAG_(SignalRange(1):SignalRange(2)));

% Normalization Factor (all states)
NormF_ALL(idx.count) = (sum(sum(sum(abs(fftshift(...
    fft_REAL(idx.count,:)+1i*fft_IMAG(idx.count,:),2))))));

% Normalization Factor (states of interest)
Ncalc_ = -(PULSE.CalcStates(end)*SPACE.ZeroPad_Factor):(PULSE.CalcStates(end)*SPACE.ZeroPad_Factor);
NormF_CALC(idx.count) = (sum(sum(sum(abs(fftshift(...
    fft_REAL(idx.count,(round(SPACE.N0/2)+1)+Ncalc_)+...
    1i*fft_IMAG(idx.count,(round(SPACE.N0/2)+1)+Ncalc_),2))))));

% Normalization Factor (key states of focus
Nview_ = PULSE.CalcStates(end)*(-SPACE.ZeroPad_Factor:PULSE.CalcStates(end));
NormF_VIEW(idx.count) = (sum(sum(sum(abs(fftshift(...
    fft_REAL(idx.count,(round(SPACE.N0/2)+1)+Nview_)+...
    1i*fft_IMAG(idx.count,(round(SPACE.N0/2)+1)+Nview_),2))))));

% Frequency Shifted FFT
FreqShift_ = (fftshift(fft_REAL(idx.count,:)+1i*fft_IMAG(idx.count,:),2))./abs(NormF_ALL(idx.count));
FreqShift(idx.count,:) = abs(FreqShift_);
% FreqShift(idx.count,:) = abs(sqrt(conj(FreqShift_).*(FreqShift_)));



% FFT --> angular frequency
% fft_TOTAL(1+idx.count,:) = fft((Psi_TOTAL(1+idx.count,:)), SPACE.N0);


% sumNorm = 0;
% for n_ = 0:4
%     sumNorm = sumNorm + (sum(sum(sum(abs(...
%     fftshift(fft_REAL(idx.count,SPACE.N/2+1+n_)...
%     +1i*fft_IMAG(idx.count,SPACE.N/2+1-n_),2))))));
% end
% NormF_mid(idx.count) = sumNorm;


% for n=(1+SPACE.ZeroPad_Factor):(SPACE.N0-SPACE.ZeroPad_Factor)
%    FreqShift_(n) = FreqShift_(n+SPACE.ZeroPad_Factor);
% end
% for n=1:SPACE.ZeroPad_Factor; FreqShift_(n) = 0; end
% for n=(SPACE.N0-SPACE.ZeroPad_Factor):SPACE.N0; FreqShift_(n) = 0; end



% Kinetic Energy
[Kin_, KE_En_grad] = KE_gradient(SPACE, TIME, Psi_REAL(idx.count,:), Psi_IMAG(idx.count,:));
T_kin_REAL(idx.count,:) = real(Kin_);
T_kin_IMAG(idx.count,:) = imag(Kin_);

% T_kin_REAL(idx.count,:) = real(sqrt(hPlanck/(2*pi))*...
%     ((fftshift(fft_REAL(idx.count,:)+1i*fft_IMAG(idx.count,:))).^2)/(2*me));
% T_kin_IMAG(idx.count,:) = imag(sqrt(hPlanck/(2*pi))*...
%     ((fftshift(fft_REAL(idx.count,:)+1i*fft_IMAG(idx.count,:))).^2)/(2*me));


% Instantaneous Energy
En_REAL(idx.count) = real(...
    (SPACE.delta*trapz(conj(Psi_REAL(idx.count,:)+1i*Psi_IMAG(idx.count,:)).*...
    (T_kin_REAL(idx.count,:) + 1i*T_kin_IMAG(idx.count,:)+V_pot(idx.count,:)).*...
    ((Psi_REAL(idx.count,:)+1i*Psi_IMAG(idx.count,:)))))...
    ./(SPACE.delta*trapz(conj((Psi_REAL(idx.count,:)+1i*Psi_IMAG(idx.count,:))).*...
    ((Psi_REAL(idx.count,:)+1i*Psi_IMAG(idx.count,:))))));

En_IMAG(idx.count) = imag(...
    (SPACE.delta*trapz(conj(Psi_REAL(idx.count,:)+1i*Psi_IMAG(idx.count,:)).*...
    (T_kin_REAL(idx.count,:) + 1i*T_kin_IMAG(idx.count,:)+V_pot(idx.count,:)).*...
    ((Psi_REAL(idx.count,:)+1i*Psi_IMAG(idx.count,:)))))...
    ./(SPACE.delta*trapz(conj((Psi_REAL(idx.count,:)+1i*Psi_IMAG(idx.count,:))).*...
    ((Psi_REAL(idx.count,:)+1i*Psi_IMAG(idx.count,:))))));

En_tot(idx.count) = En_REAL(idx.count) + 1i*En_IMAG(idx.count);


% Normalization Factor
NormT(idx.count) = sum((abs(Psi_REAL(idx.count,:)+1i*Psi_IMAG(idx.count,:))).^2);
% NormT(idx.count) = norm(Psi_REAL(idx.count,:)+1i*Psi_IMAG(idx.count,:),2);

% Photons Absorbed
PhotonsAbsorbed(idx.count) = En_REAL(idx.count)/(hPlanck*LASER.w0);

% Instantaneous Ionization Rate
% IIR(idx.count) = -(2*En_IMAG(idx.count)/NormT(idx.count));
IIR(idx.count) = (-(2*En_IMAG(idx.count))/e0)/NormT(idx.count);

% Centroid (position)
% Centroid_X(idx.count) = trapz(conj(Psi_REAL(idx.count,:)+1*Psi_IMAG(idx.count,:)).*...
%         SPACE.Axis.*(Psi_REAL(idx.count,:)+1*Psi_IMAG(idx.count,:)))...
%         /trapz(conj(Psi_REAL(idx.count,:)+1*Psi_IMAG(idx.count,:))...
%         .*(Psi_REAL(idx.count,:)+1*Psi_IMAG(idx.count,:)));
Centroid_X(idx.count) = trapz(conj(Psi_REAL(idx.count,:)+1*Psi_IMAG(idx.count,:)).*...
        SPACE.Axis_.*(Psi_REAL(idx.count,:)+1*Psi_IMAG(idx.count,:)))...
        /trapz(conj(Psi_REAL(idx.count,:)+1*Psi_IMAG(idx.count,:))...
        .*(Psi_REAL(idx.count,:)+1*Psi_IMAG(idx.count,:)));




%**************************************************************************
% Population State
%--------------------------------------------------------------------------
% >>>   Population of states function here works really only for simplistic
%       systems, such as square wells or harmonic oscillators.
% 
% >>>   In next version (once next paper is released), hydrogen population
%       will be included.
% 
%..........................................................................

SpreadStates_radius = floor(PULSE.SpreadStates/2);
Pop_rng = ...
    ((round(SPACE.N/2)+1)-SpreadStates_radius):...
    ((round(SPACE.N/2)+1)+SpreadStates_radius);
Pop_Freq = FreqShift_(Pop_rng);

[PopulationState(idx.count,:),PopulationStateNorm_GIVEN(idx.count,:),...
PopulationStateNorm_CALC(idx.count,:),PopulationStateNorm_VIEW(idx.count,:)] ...
= StatePopulations(SIM, PULSE, SPACE, Pop_Freq);

if idx.count == DATA.FFT_MEM.MaxTimeSteps

%**************************************************************************
% SAVE TO DISK
%--------------------------------------------------------------------------    
    fprintf('\n\n\t\t\t***COPYING DATA TO DISK***\n\n:'); tic;

for n=1:length(FFT_Save_Mem)
    
    fprintf('\n%s\t:', getfield(FFT_Save_Mem{n},'VarName')); tic;
    
    [RESULTS] = SaveFormat(idx, FFT_Save_Mem{n});

%     eval(SaveString);
    
    fprintf('\t\t%f\tSeconds\n', toc);
end


%**************************************************************************
% CLEAR MEMORY
%--------------------------------------------------------------------------
fprintf('\n\n\t\t\t***CLEARING MEMORY***\n\n:');% tic;

for n=1:length(FFT_Total_Mem)
    [RESULTS] = VarInitialize(idx,FFT_Total_Mem{n});
end


%**************************************************************************
% RESET INDEX
%--------------------------------------------------------------------------
    idx.START = idx.START+idx.count;
    idx.END = idx.END+idx.count;
    idx.count = 1;
    
    if idx.END>TIME.saveNum; idx.END = TIME.saveNum; end
    

%**************************************************************************
% LOAD FROM DISK
%--------------------------------------------------------------------------

for n=1:length(FFT_Load_Mem)

    % Requires variable name search from earlier
    [RESULT] = LoadFormat(idx,FFT_Load_Mem{n});
    
end
    
end
 
% Iterate Memory Index
idx.count = idx.count+1;
    
   
    
show_t = (TIME.saveNum/TIME.save);
if (show_t>1e6)
    show_t = 1e6;
end
if (show_t<1e3)
    show_t = 1e4;
end
%     if ~mod(t,show_t)
        if ~mod(t,1e2)

fprintf('\tMomentum(P-space)\t:::\t\t[%d\tof\t%d])\t=\t<%d*(s)>\n',...
    t*TIME.save, TIME.N, t*TIME.save*TIME.delta);
        
    end
    
   

end


idx.count = idx.count-1;



%% Save & Cleanup



%**************************************************************************
% SAVE DATA
%--------------------------------------------------------------------------
fprintf('\nSaving FDTD-Q (Momentum Space) Data...\n');
% tic;    


for n=1:length(FFT_Save_Mem)
    
    fprintf('\n%s\t:', getfield(FFT_Save_Mem{n},'VarName')); tic;
    
    [RESULTS] = SaveFormat(idx,FFT_Save_Mem{n});
    
    fprintf('\t\t%f\tSeconds\n', toc);
end



%**************************************************************************
% CLEAR MEMORY
%--------------------------------------------------------------------------
fprintf('\n\n\t\t\t***CLEARING MEMORY***\n\n:');% tic;

for n=1:length(FFT_Total_Mem)
    [RESULTS] = CleanUpVars(FFT_Total_Mem{n});
end




fprintf('___________________________________________________________________________\n\n\n');


end
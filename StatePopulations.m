function [PopulationState_, PopulationStateNorm_GIVEN,...
    PopulationStateNorm_CALC, PopulationStateNorm_VIEW]...
    =    StatePopulations(SIM, PULSE, SPACE, Pop_Freq)
% 
% 
%==========================================================================
%% DESCROPTION






%% Sanity Check

% Input Arguments
narginchk(1,inf);

% Output Arguments
nargoutchk(0,4);




%% Local Variables
CONSTANTS;

% PopulationState_ = zeros(1,PULSE.NumCalcStates);
PopulationState_CALC = zeros(1,PULSE.NumCalcStates);
PopulationState_VIEW = zeros(1,PULSE.NumViewStates);
PopulationState_GIVEN = zeros(1,length(Pop_Freq));

PopulationStateNorm_CALC = zeros(1,PULSE.NumCalcStates);
PopulationStateNorm_VIEW = zeros(1,PULSE.NumViewStates);
PopulationStateNorm_GIVEN = zeros(1,length(Pop_Freq));


%% Population of States

% Center-point Index (state 1 frequency coefficient)
% CtrPnt = floor(length(Pop_Freq)/2);
CtrPnt = floor(length(Pop_Freq)/2);
Estate = sum(PULSE.EnergyState)/length(PULSE.EnergyState);

ZeroPad_ = (SPACE.N0/SPACE.N);

for n=1:PULSE.NumCalcStates
    N = n-1;
%     C_=Estate*((N^2+Estate)/(N*Estate));
%     C_=Estate*((n+Estate)/(Estate));
    
C_=1;
    
%     N_ = (2*(N-1)):(2*N-1);
    
    
    if n==1
        
        switch SIM.PotentialMap
            case 'Square'
            % State Index (state N frequency coefficient)
%     N_ = (N-1)*[-ZeroPad_ ZeroPad_];
%     N_ = N*[-ZeroPad_ ZeroPad_];
%     N_ = ZeroPad_*[(1-N) (1+N)];

% 	N_ = [  (1-N)*ZeroPad_...
%             ZeroPad_...
%             (N+1)*ZeroPad_];
% N_ = (N-1)*round(ZeroPad_);


%     N_ = ZeroPad_*[(1-N) 1 (1+N)];
    N_ = [(1-N) 1 (1+N)];


            case 'H1'
%     N_ = (ZeroPad_)*(1+[(1-N) (1+N)]);
%     N_ = (ZeroPad_)*[(1-N) (1+N)];
%     N_ = (ZeroPad_)*[(1-N) 1 (1+N)];
    N_ = [(1-N) 1 (1+N)];
                
                
            otherwise
            N_ = [  (-N+1)*ZeroPad_...
            ZeroPad_...
            (N+1)*ZeroPad_];
        end
%         StateCount_ = sum(abs(Pop_Freq(CtrPnt+N_)));
% 	StateCount_ = sum(abs(Pop_Freq(CtrPnt+[0 0])));
%     StateCount_ = 1.5*sum(abs(Pop_Freq(CtrPnt)));
% StateCount_ = C1_*sum(abs(Pop_Freq(CtrPnt + N_)));

%     StateCount_ = 1.5*sum(abs(Pop_Freq(CtrPnt+N_)));
%         StateCount_ = C1_*sum(Pop_Freq(CtrPnt+N_));
%         StateCount_ = C_*sum(Pop_Freq(CtrPnt+N_));
%         StateCount_ = C_*sum((Pop_Freq(CtrPnt+N_)));
        StateCount_ = C_*sum(abs(Pop_Freq(CtrPnt+N_)));

%         
    else
        
        switch SIM.PotentialMap
            
            case 'Square'
%     N_ = [(ZeroPad_-N) (ZeroPad_+N)];
        
%                 N_ = round(ZeroPad_)*[(1-N) (1+N)];
%                 N_ = round(ZeroPad_)*[(1-N) 1 (1+N)];
                N_ = [(1-N) 1 (1+N)];
        
    
            case 'H1'
%                 N_ = round(ZeroPad_)*(1+[(1-N) (N+1)]);
%                 N_ = round(ZeroPad_)*[(-N) (N)];
%                 N_ = round(ZeroPad_)*[(1-N) (1+N)];
%                 N_ = round(ZeroPad_)*[(1-N) 1 (1+N)];
                N_ = [(1-N) 1 (1+N)];
                
            otherwise
%                 N_ = round(ZeroPad_)*[(1-N) 1 (1+N)];
        end
    
%         StateCount_ = C_*sum(Pop_Freq(CtrPnt+N_));
%         StateCount_ = C_*sum(Pop_Freq(CtrPnt+N_));
%         StateCount_ = C_*sum((Pop_Freq(CtrPnt + N_)));
        StateCount_ = C_*sum(abs(Pop_Freq(CtrPnt + N_)));
    end
    
    PopulationState_(n) = StateCount_;
% 	PopulationState_(N) = sum(abs(Pop_Freq(CtrPnt+N_)));
    
end

% Ncalc_ = -PULSE.NumCalcStates:PULSE.NumCalcStates;
% Nview_ = -PULSE.NumViewStates:PULSE.NumViewStates;
% NormFactor_CALC = sum(abs(PopulationState_(CtrPnt+Ncalc_)));
% NormFactor_VIEW = sum(abs(PopulationState_(CtrPnt+Nview_)));

NormFactor_GIVEN = sum(abs(Pop_Freq));
NormFactor_CALC = sum(abs(PopulationState_(1:PULSE.NumCalcStates)));
NormFactor_VIEW = sum(abs(PopulationState_(1:PULSE.NumViewStates)));    

    if (NormFactor_GIVEN<(1e-6)); NormFactor_GIVEN = 1e-6; end
    if (NormFactor_CALC<(1e-6)); NormFactor_CALC = 1e-6; end
    if (NormFactor_VIEW<(1e-6)); NormFactor_VIEW = 1e-6; end
    
for M=1:length(Pop_Freq)
    PopulationStateNorm_GIVEN(M) = abs(Pop_Freq(M))/NormFactor_GIVEN;
end
for M=1:PULSE.NumViewStates
    PopulationStateNorm_VIEW(M) = PopulationState_(M)/NormFactor_VIEW;
end
for M=1:PULSE.NumCalcStates
    PopulationStateNorm_CALC(M) = PopulationState_(M)/NormFactor_CALC;
end

%% OUTPUT

% Popluation of States: Calculated
% (i.e. ALL the ones we used)
% PopulationState_;

% Population of States: Displayed
% (i.e. the ones we care about)
% PopulationStateNorm_;



end
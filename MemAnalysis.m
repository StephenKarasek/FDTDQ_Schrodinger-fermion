function [varargout] = MemAnalysis(useGPU, varargin)
%%
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
% VarName, VarType, VarSize, Space_Steps, Dims, Time_Steps





% MemAnalysis(FDTDQ_file.vars, FDTDQ_file.vartype, FDTDQ_file.VarSize, N, process_dims, T_len, true);

% useGPU = false;
% VarName = {'Psi_REAL', 'Psi_IMAG', 'V_pot'};
% VarType = {1, 1, 1};
% tsize = [1e5];
% rsize = [256 256];
% VarSize = {[1e4 rsize], [1e4 rsize] ,[1e4 rsize]};



%% Sanity Check


% Input Arguments
narginchk(2,inf);

% Output Arguments
nargoutchk(4,4);



%% Var assignment

% useGPU      = varargin{1};
% Parallel Processing


for n=1:length(varargin)

    VarName{n}	= getfield(varargin{n}, 'VarName');
    VarSize{n}	= getfield(varargin{n}, 'VarSize');
    VarType{n}  = getfield(varargin{n}, 'VarType');

end




%% Local Variables

[userview systemview] = memory;
% eff_mem = systemview.PhysicalMemory.Available*.7;
% eff_mem = systemview.PhysicalMemory.Total*.01;

eff_mem = systemview.PhysicalMemory.Available*.5;

REAL_num = 8;
IMAG_num = 16;

% SpaceBlock = SpaceSteps^Dims;

% Quantity of memory required per time step
memWeight           = zeros(length(VarName),1);
spaceStep           = zeros(length(VarName),1);
MaxSize             = zeros(length(VarName),1);
timeStep            = zeros(length(VarName),1);
SizeToAllocate      = zeros(length(VarName),1);



%% Memory consumed by Variables

for n=1:length(VarName)
    
    memWeight(n) = REAL_num;
    
    
% switch isreal(VarType(n))
%     case 1
%         memWeight(n) = REAL_num;
%         
% 	case 0
%         memWeight(n) = IMAG_num;
%         
%     otherwise% assume worst case scenario
%         memWeight(n) = IMAG_num;
        
% end

% Spatial Dimensions of current; 
% centroids are 1-dimensional, stretching only over time
% otherwise, should stretch over each spatial dimension in solution space
Dims = length(VarSize{n});

% Contribution to size over space, per unit time
spaceStep(n) = 1;
for d=0:Dims-2
    spaceStep(n) = spaceStep(n)*VarSize{n}(2+d);
end
% spaceStep(n) = VarSize{n}(2:Dims);


% Time Steps to run
timeStep(n) = VarSize{n}(1);

% Total size of fully allocated variable
SizeToAllocate(n) = timeStep(n)*((spaceStep(n)*memWeight(n)));

% Maximum time steps that may be allocated at once, based on available memory
MaxSize(n) = eff_mem/((spaceStep(n)));
% *memWeight(n)

end


% Memory consumed per Time-Step
PerTimeStep = sum(spaceStep.*memWeight);

% Maximum time steps runnable, growing across all variables
MaxTimeSteps = floor(eff_mem/PerTimeStep);

% 
MaxSize_Allocate = max(MaxSize);



%% Outputs
varargout{1} = MaxTimeSteps;
varargout{2} = PerTimeStep;
varargout{3} = SizeToAllocate;
varargout{4} = MaxSize;




end


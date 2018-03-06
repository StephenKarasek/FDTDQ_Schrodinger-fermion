function [RESULTS] = CleanUpVars(varargin)

% [Success] = VarInitialize(idx, VarName, VarSize, VarType, VarIO);

%% I/O
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%|||||||||||||||||||||||||||||   INPUTS   |||||||||||||||||||||||||||||||||
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%
% VarName:          'Psi_Real
% 
% VarIO
% 
% MemIndex1
% 
% 
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%|||||||||||||||||||||||||||||   OUTPUTS   ||||||||||||||||||||||||||||||||
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%
% VarName
% 
% VarIO
% 
% MemIndex
% 
% 
% 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%% Sanity Check

% Input Arguments
narginchk(1,inf); %inf

% Output Arguments
nargoutchk(0,1);




%% Local variable assignment

switch nargin>1
    
    case 0
        Var{1}  = getfield(varargin{1},'VarName');

    case 1
        
        for n=length(varargin)
           Var{n} = varargin{n};            
        end
end


% VarType{n}  = getfield(varargin{n},'VarType');


%% Initialize Var
try
    for n=1:length(varargin)
        eval(['clear ' Var{n} ';']);
        RESULTS(n) = true;
    end

catch
    RESULTS(n) = false;

end

%% Return
varargout{1} = RESULTS;






end
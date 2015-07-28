function [varargout] = VarInitialize(idx, varargin)
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
narginchk(2,2); %inf

% Output Arguments
nargoutchk(0,1);




%% Local variable assignment

for n=length(varargin)
    VarName  = getfield(varargin{n},'VarName');
    VarSize  = getfield(varargin{n},'VarSize');
    VarType  = getfield(varargin{n},'VarType');
end


% VarType{n}  = getfield(varargin{n},'VarType');


%% Initialize Var

% % Temporal Dimensions
% Initialize_str = sprintf([VarName ' = zeros(' num2str(VarSize(1))]);
% 
% 
% % Spatial Dimensions
% n_ = length(VarSize);
% for n=2:n_
%     Initialize_str = sprintf([Initialize_str ', ' num2str(VarSize(n))]);
% end
% 
% % Check Type
% Initialize_str = sprintf([Initialize_str ', ''like'', ' num2str(VarType) ');']);
% 
% % Run Initialization
% eval(Initialize_str);




% for m=1:length(varargin)
%**************************************************************************
% IO-stream saving TO
%--------------------------------------------------------------------------
% InitializeString = sprintf([VarName{m} ' = zeros(' num2str(idx.length)]);
if VarSize(1)>idx.length
    InitializeString = sprintf([VarName ' = zeros(' num2str(idx.length)]);
else
    InitializeString = sprintf([VarName ' = zeros(' num2str(VarSize(1))]);
end

% N-dimensions
n_ = length(VarSize);
% n_ = length(VarSize{m});

for n=2:n_
    InitializeString = sprintf([InitializeString ', ' num2str(VarSize(n))]);
    
%     InitializeString = sprintf([InitializeString ', ' num2str(VarSize{m}(n))]);
end

InitializeString = sprintf([InitializeString ', ''like'', 1);']);
% InitializeString = sprintf([InitializeString ', ''like'', ' num2str(VarType{m}) ');']);


try
    evalin('caller', InitializeString);
    RESULTS = true;
catch
    RESULTS = false;
end

% end

%% Return
varargout{1} = RESULTS;






end
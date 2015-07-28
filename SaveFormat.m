function [varargout] = SaveFormat(idx, varargin)

% [FormatString] = SaveFormat(VarName, VarType, VarSize, FilePath)
% evalstr = SaveFormat('Psi_REAL', 1, [(1e6) (256) (256)], fpath);
%% SaveFormat
% INPUT % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%
% VarName
%	'Psi_REAL'
%
% VarSize
%	[1e4 64 64 64]
%
% VarIO
%	'FDTDQ_file
%
% idx
%	idx.START   = 5e5
%   idx.END     = 1e6
%   idx.length  = 5e5
%
%
% OUTPUT % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%
% FormatString
%	['FDTDQ_file.IO{n}.' FDTDQ_file.vars{n} ' = ' FDTDQ_file.vars{n} ';']
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%


%% Sanity Check

% Input Arguments
narginchk(2,2); %inf

% Output Arguments
nargoutchk(0,1);


%% Variable Assignment

for n=length(varargin)
    VarName     = getfield(varargin{n},'VarName');
    VarSize     = getfield(varargin{n},'VarSize');
    VarIO       = getfield(varargin{n},'VarIO');
    VarLoc      = getfield(varargin{n},'VarLoc');
end


%%


%     fprintf('\n%s\t:', FDTDQ_file.vars{n});


% for m=1:length(varargin)
%**************************************************************************
% SAVE data to storage source...
%--------------------------------------------------------------------------
FormatString = ['.' VarName '(' num2str(idx.START) ':' num2str(idx.END)];
% FormatString = ['.' VarName '(' num2str(idx.START) ':' num2str(idx.START+idx.count-1)];

% N-dimensions
% n_ = length(VarSize{m});
n_ = length(VarSize);

for n=2:n_
    FormatString = sprintf([FormatString ', ' num2str(1) ':' num2str(VarSize(n))]);
    
%     FormatString = sprintf([FormatString ', ' num2str(1) ':' num2str(VarSize{m}(n))]);
end

FormatString = sprintf([FormatString ') = ']);

%**************************************************************************
% Use this Variable as the source
%--------------------------------------------------------------------------
% FormatString = sprintf([FormatString VarName{m} '(' num2str(1) ':' num2str(idx.count)]);
FormatString = sprintf([FormatString VarName '(' num2str(1) ':' num2str(idx.count)]);

% N-dimensions
% n_ = length(VarSize{m});
n_ = length(VarSize);

for n=2:n_
    FormatString = sprintf([FormatString ', ' num2str(1) ':' num2str(VarSize(n))]);
    
%     FormatString = sprintf([FormatString ', ' num2str(1) ':' num2str(VarSize{m}(n))]);
end
FormatString = sprintf([FormatString ');']);



%**************************************************************************
% Run
%--------------------------------------------------------------------------
SAVE_FormatString = [VarLoc '.' VarName '.VarIO' FormatString];
try
    evalin('caller', SAVE_FormatString);
    RESULTS = true;
catch
    RESULTS = false;
end


% end

%%

%     fprintf('\t\t%f\tSeconds\n', toc);
% eval(['VarIO' FormatString]);





%% RETURN

varargout{1} = RESULTS;


end





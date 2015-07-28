function [varargout] = LoadFormat(idx, varargin)
% VarName, VarIO, MemIndex)


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

%%
% VarName, VarIO, MemIndex


%% Local variable assignment
% idx;%     =       varargin{1};

for n=length(varargin)
    VarName      = getfield(varargin{n},'VarName');
    VarSize      = getfield(varargin{n},'VarSize');
    VarIO        = getfield(varargin{n},'VarIO');
    VarLoc       = getfield(varargin{n},'VarLoc');
end

%%


if (idx.length > VarSize(1))
    idx.length = VarSize(1);
    idx.END = idx.length;
end

% for m=1:length(varargin)
%**************************************************************************
% Overwrite onto Variable...
%--------------------------------------------------------------------------
FormatString = [VarName '(' num2str(1) ':' num2str(idx.length)];
% FormatString = [VarName{m} '(' num2str(1) ':' num2str(idx.length)];

% N-dimensions
n_ = length(VarSize);
% n_ = length(VarSize{m});
for n=2:n_
    FormatString = sprintf([FormatString ', ' num2str(1) ':' num2str(VarSize(n))]);
%     FormatString = sprintf([FormatString ', ' num2str(1) ':' num2str(VarSize{m}(n))]);
end

FormatString_to = sprintf([FormatString ') = ']);


%**************************************************************************
% LOAD from storage source...
%--------------------------------------------------------------------------
FormatString = ['.' VarName '(' num2str(idx.START) ':' num2str(idx.END)];
% FormatString = [FormatString 'VarIO.' VarName '(' num2str(idx.START) ':' num2str(idx.END)];

% N-dimensions
n_ = length(VarSize);
for n=2:n_
    FormatString = sprintf([FormatString ', ' num2str(1) ':' num2str(VarSize(n))]);
end

FormatString_from = sprintf([FormatString ');']);


% DATA.FDTD_.Psi_IMAG.VarIO.Psi_IMAG(1:33455, 1:128);

%**************************************************************************
% Run
%--------------------------------------------------------------------------
LOAD_FormatString = [FormatString_to VarLoc '.' VarName '.VarIO' FormatString_from];
try
    evalin('caller', LOAD_FormatString)
    RESULTS = true;
catch
    RESULTS = false;
end





%%
% varargout = 0;
%     fprintf('\t\t%f\tSeconds\n', toc);
% eval(['VarIO' FormatString]);
% LOAD_FormatString = sprintf(['VarIO' FormatString]);
varargout{1} = RESULTS;





end
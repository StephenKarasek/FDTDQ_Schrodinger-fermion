function [idx] = createMemIndex(varargin)
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
%   idx.loops   = 0;
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
narginchk(2,2);

% Output Arguments
nargoutchk(1,1);


%% Local variable assignment
MaxTimeSteps    =	varargin{1};
Time_N          =	varargin{2};


%%
fprintf('\n\n\t--> CREATING MEMORY INDEX... ');


idx.START = 1;
if MaxTimeSteps < Time_N
    idx.length = idx.START + MaxTimeSteps;
else
    idx.length = Time_N;
end
idx.END = idx.START + idx.length - 1;
idx.count = 1;


fprintf('\t\t... DONE!\n\n:');




%% RETURN

% idx.loops

% varargout{1} = 0;




end
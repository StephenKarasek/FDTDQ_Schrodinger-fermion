function [VarOut] = SaveInitialize(SaveDirectory, varargin)
% function [varargout] = SaveInitialize(SaveDirectory, VarName, VarSize, VarOrigin)
% SaveInitialize(VarName, VarType, VarSize, SaveDirectory)
%% SaveInitialize
% INPUT % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%
% VarName           = 'Psi_REAL'
%
% VarType           = 1
% 
% VarSize           = [1e4 64 64]
% 
% SaveDirectory     = 'D:\...'
% 
%
% OUTPUT % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%
% VarIO             = matlab.io.MatFile
% 
% Filepath          = 'D:\...'
% 
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%


%% Sanity Check

% Input Arguments
narginchk(2,5);

% Output Arguments
nargoutchk(0,1);



%% Variable Assignment

switch nargin
    
    case 2
        VarName     = getfield(varargin,'VarName');
        VarSize     = getfield(varargin,'VarSize');
        VarType     = getfield(varargin,'VarType');
        VarLoc      = getfield(varargin,'VarLoc');
        
    case 5
        VarName     = varargin{1};
        VarSize     = varargin{2};
        VarType     = varargin{3};
        VarLoc      = varargin{4};
        
    otherwise
        error('Incorrect number of arguments... Please check inputs and try again.');

end


%% Initialize Variable
% fprintf('\n\n\t--> BUILDING DATA STRUCTURES...\t\t:');

[CREATEstr] = CreateVar(VarName, VarSize);
eval(CREATEstr);

% InitializeString = sprintf([VarName ' = zeros(' num2str(VarSize(1))]);
% 
% % N-dimensions
% n_ = length(VarSize);
% for n=2:n_
%     InitializeString = sprintf([InitializeString ', ' num2str(VarSize(n))]);
% end
% 
% InitializeString = sprintf([InitializeString ', ''like'', ''1'');']);
% % InitializeString = sprintf([InitializeString ', ''like'', ' num2str(VarType) ');']);
% 
% 
% eval(InitializeString);
    


%% Matfile

[IOstr, SaveMATstr, FilePath_str] = MatFile(SaveDirectory, VarName);
eval(SaveMATstr);
VarIO = eval(IOstr);
% assignin('base', [], IOstr

% % Create MatFile
% FilePath_str = [SaveDirectory '\' VarName '.mat'];
% save(FilePath_str, VarName, '-v7.3');
% 
% % Create IO Interface
% VarIO = eval(['matfile(''' FilePath_str ''')']);
% 
% % Force writable
% VarIO.Properties.Writable = true;


%% Assign Variable w/ IO Interface

% Assign_str = sprintf(['VarIO.' VarName '_(' num2str(1) ':' num2str(VarSize(1))]);
% 
% % N-dimensions
% n_ = length(VarSize);
% for n=2:n_
%     Assign_str = sprintf([Assign_str ', ' num2str(1) ':' num2str(VarSize(n))]);
% end
% 
% Assign_str = sprintf([Assign_str ') = ' VarName '(' num2str(1) ':' num2str(VarSize(1))]);
% 
% % N-dimensions
% n_ = length(VarSize);
% for n=2:n_
%     Assign_str = sprintf([Assign_str ', ' num2str(1) ':' num2str(VarSize(n))]);
% end
% Assign_str = sprintf([Assign_str ');']);
% 
% 
% % Assign Variable
% eval(Assign_str);

%% Clear Memory
eval(['clear ' VarName ';']);



%% Output
% varargout{1} = eval(['VarIO.' VarName '_']);

% switch nargout
    
%     case 1
VarOut.('VarIO')    = VarIO;
VarOut.('FilePath') = FilePath_str;
VarOut.('VarName')  = VarName;
VarOut.('VarSize')  = VarSize;
VarOut.('VarType')  = VarType;
VarOut.('VarLoc')   = VarLoc;
        
% varargout{2} = FilePath_str;

% VarIO = 0;
% FilePath = FilePath_str;



end





%% Initialize Variable
function [CREATEstr] = CreateVar(VarName, VarSize)


CREATEstr = sprintf([VarName ' = zeros(', num2str(VarSize(1))]);

% N-dimensions
n_ = length(VarSize);
for n=2:n_
    CREATEstr = sprintf([CREATEstr ', ' num2str(VarSize(n))]);
end

CREATEstr = sprintf([CREATEstr ', ''like'', 1);']);
% InitializeString = sprintf([InitializeString ', ''like'', ' num2str(VarType) ');']);


% evalin('base', InitializeString);



end






%% Matfile
function [IOstr, SaveMATstr, FilePath_str] = MatFile(SaveDirectory, VarName)
% Create MatFile
FilePath_str = [SaveDirectory '\' VarName '.mat'];
SaveMATstr = ['save(''' FilePath_str ''', ''' VarName ''', ''-v7.3'');'];

% Create IO Interface
IOstr = sprintf(['matfile(''%s'', ''Writable'', %s)'], FilePath_str, 'true');

% Hstr = sprintf(['matfile(''' FilePath_str ''', ''Writable'', %s)'], 'true')

% VarIO = eval(sprintf([VarOrigin '.' VarName '.VarIO = ' Hstr]));

% INstr = sprintf(['matfile(''' FilePath_str ''', ' num2str(1) ')']);
% VarIO = eval('matfile(''' FilePath_str ''')']);

% Force writable
% VarIO.Properties.Writable = true;





end
function [SPACE_] = BoundaryConditions(SPACE, TIME, SIM, PULSE, varargin)

% 


%% Sanity Check

% Input Arguments
narginchk(5,6);

% Output Arguments
nargoutchk(1,1);

switch nargin
    case 5
    PsiR = varargin{1};
    
    case 6
    PsiI = varargin{2};
    
    otherwise
        error('Incorrect number of inputs! Please check function call and try again.');
end



%% Local variable assignment

CONSTANTS;




%%

switch SIM.FDM

    case 'FDTDQ'

        
    case 'GFDTDQ'
%         Space = 0;


    otherwise
%         Space = 0;
end



%% RETURN




end

































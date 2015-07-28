function [R_new, I_new, R_zai, I_zai, R_old, I_old] =...
    FDM_fdtdq_1D(SPACE, TIME, SIM, PULSE, varargin)
% function [varargout] = BC_fdtdq(SPACE, TIME, SIM, PULSE, varargin)
% R_new, I_new, R_zai, I_zai, R_old, I_old)

%% Boundary Conditions for FDTD-Q scheme


%% Sanity Check

% Input Arguments
narginchk(7,9);

% Output Arguments
nargoutchk(2,2);
% nargoutchk(1,1);

switch nargin
    case 7
% Future Steps
    R_new = varargin{1};
    I_new = varargin{2};
    V_new = varargin{3};
    
    case 9
% Future Steps
	R_new = varargin{1};
    I_new = varargin{2};
    V_new = varargin{3};
% Current Steps
    R_zai = varargin{4};
    I_zai = varargin{5};
%     V_zai = varargin{6};
    
    case 10
% Future Steps
	R_new = varargin{1};
    I_new = varargin{2};
    V_new = varargin{3};
% Current Steps
    R_zai = varargin{4};
    I_zai = varargin{5};
    V_zai = varargin{6};

%     % Past Steps
%     R_old = varargin{7};
%     I_old = varargin{8};
%     V_old = varargin{9};
    
    otherwise
        error('Incorrect number of inputs! Please check function call and try again.');
end



%% Local variable assignment

CONSTANTS;


%% Dimensional Operator
switch SPACE.Dims
    case 1
        N_ = 2:SPACE.N-1;
    
    case 2
        N_ = [2:SPACE.N-1; 2:SPACE.N-1];
        
    case 3
        N_ = [2:SPACE.N-1; 2:SPACE.N-1; 2:SPACE.N-1];
        
    otherwise
   error('Dimension parameters not defined properly. Please fix and try again.');
end


%% Main Function

switch SIM.FDM.Type

    
% %%
%**************************************************************************
    case {'FDTDQ','fdtdq'}
%--------------------------------------------------------------------------
% Boundaries
R_new(N_) = R_zai(N_) ...
    - ((hPlanck*TIME.delta)/(me*SPACE.delta^2))...
    *(I_zai(N_+1)-2*I_zai(N_+0)+I_zai(N_-1))...
    + ((TIME.delta)/(hPlanck))*V_new(N_).*I_zai(N_);

% R_new(1) = 0; R_new(SPACE.N) = 0;
% R_zai = R_new;

% R_new(1) = 0; R_new(SPACE.N)=0;
% I_new(1)=0; I_new(SPACE.N)=0;
%..........................................................................
% Update Equations
% R_zai = R_new; I_zai = I_new;
%..........................................................................

I_new(N_) = I_zai(N_) ...
    + ((hPlanck*TIME.delta)/(me*SPACE.delta^2))...
    *(R_new(N_+1)-2*R_new(N_+0)+R_new(N_-1))...
    - ((TIME.delta)/(hPlanck))*V_new(N_).*R_new(N_);

% I_new(N_) = I_zai(N_) ...
%     + ((hPlanck*TIME.delta)/(me*SPACE.delta^2))...
%     *(R_zai(N_+1)-2*R_zai(N_+0)+R_zai(N_-1))...
%     - ((TIME.delta)/(hPlanck))*V_new(N_).*R_zai(N_);

% I_new(1) = 0; I_new(SPACE.N) = 0;
% I_zai = I_new;

%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
% Output
% varargout{1} = R_new;
% varargout{2} = I_new;
%//////////////////////////////////////////////////////////////////////////


    otherwise
%         Space = 0;
end



%% RETURN




end

































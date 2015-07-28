function [R_new, I_new, R_zai, I_zai, R_old, I_old] =...
    BC_fdtdq_1D(SPACE, TIME, SIM, PULSE, varargin)
% function [varargout] = BC_fdtdq(SPACE, TIME, SIM, PULSE, varargin)
% R_new, I_new, R_zai, I_zai, R_old, I_old)

%% Boundary Conditions for FDTD-Q scheme


%% Sanity Check

% Input Arguments
narginchk(6,10);

% Output Arguments
nargoutchk(2,6);
% nargoutchk(1,1);

switch nargin
    case 6
% Future Steps
    R_new = varargin{1};
    I_new = varargin{2};

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
    
    case 10
% Future Steps
	R_new = varargin{1};
    I_new = varargin{2};
    V_new = varargin{3};
% Current Steps
    R_zai = varargin{4};
    I_zai = varargin{5};
    V_zai = varargin{6};
    
% % Past Steps
%     R_old = varargin{5};
%     I_old = varargin{6};
    
    otherwise
        error('Incorrect number of inputs! Please check function call and try again.');
end



%% Local variable assignment

CONSTANTS;




%% Main Function

switch SIM.BC.Type

    
% %%
%**************************************************************************
    case {'Inf';'Dirichlet'}
%--------------------------------------------------------------------------
% Boundaries
R_new(1)=0; R_new(SPACE.N)=0;
I_new(1)=0; I_new(SPACE.N)=0;
%..........................................................................
% Update Equations
% R_zai = R_new; I_zai = I_new;
%..........................................................................


%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
% Output
varargout{1} = R_new;
varargout{2} = I_new;
%//////////////////////////////////////////////////////////////////////////



% %%
%**************************************************************************
    case {'Mask'; 'Dirichlet_Mask'}
%--------------------------------------------------------------------------
% Boundaries
R_new(1)=0; R_new(SPACE.N)=0;
I_new(1)=0; I_new(SPACE.N)=0;
%..........................................................................
% Mask Function
MaskFunc{1} = cos((pi/2)*((SPACE.Axis_(SIM.BC.R{1})/...
	(SPACE.Axis_(SIM.BC.R)-SPACE.Axis_(1)))))...
	.^(SIM.BC.param_.Mask_param);
MaskFunc{2} = cos((pi/2)*((SPACE.Axis_(SIM.BC.R{2})/...
    (SPACE.Axis_(end)-SPACE.Axis_(end-SIM.BC.R)))))...
    .^(SIM.BC.param_.Mask_param);

R_new(SIM.BC.R{1}) = MaskFunc{1}.*R_new(SIM.BC.R{1});
R_new(SIM.BC.R{2}) = MaskFunc{2}.*R_new(SIM.BC.R{2});
I_new(SIM.BC.R{1}) = MaskFunc{1}.*R_new(SIM.BC.R{1});
I_new(SIM.BC.R{2}) = MaskFunc{2}.*R_new(SIM.BC.R{2});
%..........................................................................
% Update Equations
% R_zai = R_new; I_zai = I_new;
%..........................................................................


%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
% Output Assignment
varargout{1} = R_new;
varargout{2} = I_new;
%//////////////////////////////////////////////////////////////////////////

% %%
%**************************************************************************
    case {'Convolution';'Exact';'ExactStationary'}
%--------------------------------------------------------------------------
% Boundaries
R_new(1)=0; R_new(SPACE.N)=0;
I_new(1)=0; I_new(SPACE.N)=0;
%..........................................................................



%..........................................................................
% Update Equations
% R_zai = R_new; I_zai = I_new;
%..........................................................................


%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
% Output Assignment
varargout{1} = R_new; varargout{2} = I_new;
% varargout{3} = R_zai; varargout{4} = I_zai;
% varargout{5} = R_old; varargout{6} = I_old;
%//////////////////////////////////////////////////////////////////////////




    otherwise
%         Space = 0;
end



%% RETURN




end

































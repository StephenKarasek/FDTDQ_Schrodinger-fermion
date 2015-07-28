function [Psi0] = Psi0_QHO(PULSE, SIM, SPACE, TIME, varargin)
%% Initial Wavefunction (@ t = 0), of Infinite Square Well


% % % % % % % % % % % % % For 1D Problems % % % % % % % % % % % % % % % % %
% N = 256;
% Xmin = -a0; Xmax = a0;
% Space = linspace(Xmin,Xmax,N);
% InitialPosition = [(Xmin-Xmax)/2];
% nX = 1;

% Psi0 = Psi0_QHO(Space, InitialPosition, nX);


% % % % % % % % % % % % % For 3D Problems % % % % % % % % % % % % % % % % %
% N = 32;
% Xmin = -a0; Xmax = a0;
% Ymin = -a0; Ymax = a0;
% Zmin = -a0; Zmax = a0;
% Space = [linspace(Xmin,Xmax,N);
%     linspace(Ymin,Ymax,N);...
%     linspace(Zmin,Zmax,N)]
% InitialPosition = [(Xmin-Xmax)/2; (Ymin-Ymax)/2; (Zmin-Zmax)/2];
% nX = 1; nY = 1; nZ = 1;
% 
% Psi0 = Psi0_QHO(Space, InitialPosition, nX, nY, nZ);



%                       ***     1D Image        ***
% close all hidden; plot(Psi0);
%
%                       ***     2D Image        ***
% close all hidden; imagesc(Psi0);
%
%                       ***     3D Image        ***
% close all hidden; hfig = slice(Psi0, floor(N/2), N, floor(N/2));
% set(hfig(:), 'LineStyle', 'none');



%%
% DESCRIPTION % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% 
% Generates initial wavefunction, with different options for simulations of 
% Infinite Square Wells, i.e. "Particle-in-a-Box" problems.
% 
% Example; Instead of excitations using Gaussian Wavepackets, we already
% know individual states before we start. Good for concept demos.
% 
% Perhaps most useful for proofs, as a litmus test, and first step
% 
% 
% INPUT % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%
% delta         = 1e-12;
%
% Spatial step distance; for simplicity spatial derivatives are the same
% for each dimension present in problem.
%
% 
% X             = 100;
%
% Distance in spatial steps, until infinite well boundary begins
%
% 
% nX        = [1]   ... initial wavefunction is wavefunction @ state
%
% Defines enery state(s) to place initial wavefunction in.
% 
%
% OUTPUT % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%
% Psi_0     = 1d, 2d or 3d array in Position space of initial wave function
%
% Since function is run in caller workspace, no return needed.
% 
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %




%% Sanity Check

% Input Arguments
narginchk(4,8);

% assert((~mod(nargin-1,2)),...
%     'InputArgs:LogicError', ...
%     ['ERROR:\tIncorrect number of input arguments.'...
%     '\n\tFor an n-D Particle-in-a-Box, there should be 1+2*n inputs.']);

if nargin>7
    EigenSolution = varargin{4};
end
if nargin>6
    EigenEnergy = varargin{3};
end
if nargin>5
    EigenCoefficient = varargin{2};
end
if nargin>4
    EigenFunction = varargin{1};
end

% Psi0_Norm = ITP_Results.PsiNorm;
% Psi0_coeff = ITP_Results.Coefficients;
% Psi0_Eig = ITP_Results.EigenFunction;

% Output Arguments
nargoutchk(0,1);



%% Variable Check


%**************************************************************************
%   Space
%--------------------------------------------------------------------------
switch SPACE.Dims
    case 1
        Dimension = 1;
        X = SPACE.Axis_;
        x0 =  PULSE.InitialPos(1);% - X(1);% (X(end)-X(1))/2

        
    case 2
        Dimension = 2;
        X = SPACE.Axis_;
        x0 = PULSE.InitialPos(1);
        Y = SPACE.Axis_;
        y0 = PULSE.InitialPos(2);
        
    case 3
        Dimension = 3;
        X = SPACE.Axis_;
        x0 = PULSE.InitialPos(1);
        Y = SPACE.Axis_;
        y0 = PULSE.InitialPos(2);
        Z = SPACE.Axis_;
        z0 = PULSE.InitialPos(3);
           
    otherwise
    error('SolutionSpace:SpaceTime_defs', ...
        ['ERROR:\tSpace/Time size definitions are incorrect.'...
        '\n\n\t\tPlease check definitions and try again.']);
end



% 
% %**************************************************************************
% %   Energy States
% %--------------------------------------------------------------------------
% msgIdent0 = sprintf('EnergyState:Invalid');
% msgString0 = sprintf(['ERROR:\tEnergy state cannot be ''0'' or undefined'...
%     'in an Infinite Square Well, without breaking uncertainty principle.\n\n']);
% 
% switch Dimension
%     case 1
%         try
%         nX = EnergyState(1);
%         catch
%             error(msgIdent0, msgString0);
%         end
%         assert((nX~=0),msgIdent0,msgString0);
% %         
% %     case 2
% %         try
% %         nX = varargin{1};
% %         nY = varargin{2};
% %         catch
% %             error(msgIdent0, msgString0);
% %         end
% %         assert((nX+nY~=0),msgIdent0,msgString0);
% %         
% %     case 3
% %         try
% %         nX = varargin{1};
% %         nY = varargin{2};
% %         nZ = varargin{3};
% %         catch
% %             error(msgIdent0, msgString0);
% %         end
% %         assert((nX+nY+nZ~=0),msgIdent0,msgString0);
%    
%     otherwise
%         assert((nX~=0),msgIdent0,msgString0);
% end
% 




%% Local Variables

CONSTANTS;


% Excited State(s)
switch PULSE.Type
    
    case {'Gaussian'; 'Fourier-Bessel'}
        
        % Pulse Width
        k_width = PULSE.Width;
        
        % Initial Momentum
        k_0 = PULSE.Momentum;
        
        % UNUSED: desired energy state.
        % For use with Imaginary Time Propagation
%         E_desired = varargin{3};

switch Dimension
    case 1
        [psi0] = Gauss_1(X, x0, k_width, k_0);
        
        [Psi0] = NormalizeWave1(X, psi0);
%         Psi0 = psi0;
    case 2
        [psi0] = Gauss_2(X, x0, Y, y0, k_width, k_0);
        
        [Psi0] = NormalizeWave2(X, Y, psi0);
        
    case 3
        [psi0] = Gauss_3(X, x0, Y, y0, Z, z0, k_width, k_0);
        
        [Psi0] = NormalizeWave3(X, Y, Z, psi0);
        
end




        
    case {'Sine';'Cos'}
        
        
switch Dimension
    case 1
%//////////////////////////////////////////////////////////////////////////
%         switch PULSE.NumStates
%             case 1
%         nX = PULSE.EnergyState(1);
        nX = PULSE.EnergyState;
        [Psi0] = ExcitedState_1d(X, x0, nX, TIME.delta);
                
%             case 2
%         nX(1) = PULSE.EnergyState(1);
%         nX(2) = PULSE.EnergyState(2);
%         nX = PULSE.EnergyState;
%         [Psi0] = ExcitedState_1d(X, x0, nX, TIME.delta);
                
%             otherwise
% error('ERROR! No more than two states allowed! (for now...)');
%         end
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
%                 
%     case 2
% %//////////////////////////////////////////////////////////////////////////
%         switch PULSE.NumStates
%             case 1
%         nX = PULSE.EnergyState(1);
%         nY = 1;
%         [Psi0] = ExcitedState_2d(X, x0, nX, Y, y0, nY, TIME.delta);
%                 
%             case 2
%         nX = PULSE.EnergyState(1);
%         nY = PULSE.EnergyState(2);
%         [Psi0] = ExcitedState_2d(X, x0, nX, Y, y0, nY, TIME.delta);
%                 
%             otherwise
% error('ERROR! No more than two states allowed! (for now...)');
%         end
% %\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
%                 
%         
%     case 3
% %//////////////////////////////////////////////////////////////////////////
%         switch PULSE.NumStates
%             case 1
%         nX = PULSE.EnergyState(1);
%         nY = 1;
%         nZ = 1;
%         [Psi0] = ExcitedState_3d(X, x0, nX, Y, y0, nY, Z, z0, nZ, TIME.delta);
%                 
%             case 2
%         nX = PULSE.EnergyState(1);
%         nY = PULSE.EnergyState(2);
%         nZ = 1;
%         [Psi0] = ExcitedState_3d(X, x0, nX, Y, y0, nY, Z, z0, nZ, TIME.delta);
%                 
%             case 3
%         nX = PULSE.EnergyState(1);
%         nY = PULSE.EnergyState(2);
%         nZ = PULSE.EnergyState(3);
%         [Psi0] = ExcitedState_3d(X, x0, nX, Y, y0, nY, Z, z0, nZ, TIME.delta);
%         
%             otherwise
% error('ERROR! No more than 3 states allowed! (for now...)');
%         end
% %\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
%         
end



    case 'Noise'
        [psi0] = rand(size(SPACE.N,1),SPACE.N);
%         [psi0] = rand(size(SPACE.N_,1),SPACE.N_);
        
        [Psi0] = NormalizeWave1(X, psi0);
        
        
        

    case 'EigenState'
%         if nargin<6
            [Psi0] = EigenStateInitialize(PULSE, EigenFunction, SPACE.Dims, X, x0);
%             EigenStateInitialize(PULSE, EigenFunctionSolution, SPACE.Dims, X, x0);
                
%         else
%             [Psi0] = ...
%                 EigenStateInitialize_RAW(PULSE, EigenFunction,EigenCoefficient,SPACE.Dims, X, x0);
%         end
        %         switch nargin
%             case 4
%                 [Psi0] = EigenStateInitialize(X, x0, PULSE, Psi0_Norm);
%             case 5
%                 [Psi0] = EigenStateInitialize(X, x0, PULSE, Psi0_Eig, Psi0_coeff);
%             case 6
%                 [Psi0] = EigenStateInitialize(X, x0, PULSE, Psi0_Norm, Psi0_Eig);
%         end
        

    case 'Analytical'
        [Psi0] = AnalyticalState(PULSE, SIM, SPACE.Dims, TIME.delta, X, x0);
        

    otherwise
            error('SolutionSpace:SpaceTime_defs', ...
        ['ERROR:\tPulse Type Undefined.'...
        '\n\n\t\tPlease check definitions and try again.']);
end
      

%% Outputs
% varargout{1} = Psi0;
% close all hidden;
% plot(Psi0);

end






%% Normalize Wavefunction
function [Psi0_] = NormalizeWave1(X, psi0)
dx = abs(X(1)-X(2));
Lx = abs(X(end)-X(1));
% psi0_all = sum(psi0_,1);

% sum(((conj(psi0).*psi0).^2)*dx);

% normFactor = (sum(conj(psi0).*psi0.*SPACE));
% psi0/(1/abs(1-(trapz(dx*(conj(psi0).*psi0)))));


% (abs(1-(trapz(SPACE.delta*(conj(psi0).*psi0).^2))))

% (abs(1-(trapz(dx*(conj(psi0).*psi0).^2))))

% normFactor = abs(sum(X.*(conj(psi0).*psi0).^2));
% normFactor = sqrt(sum((psi0.^2)));
normFactor = sqrt(sum(abs(psi0.^2)));

% normFactor = (trapz(X,conj(psi0).*psi0));
% normFactor = (trapz(conj(psi0).*psi0));

% normFactor = (trapz(dx*(conj(psi0).*psi0)))/Lx;
% normFactor = (trapz((conj(psi0).*psi0)))/(size(psi0,2));
% normFactor = dx*(sum(conj(psi0).*psi0))/Lx;

% normFactor = (trapz(X,conj(psi0).*psi0)*dx/Lx)*(dx/Lx);
% normFactor = (trapz(X,conj(psi0).*psi0)*dx/Lx);
% normFactor = trapz(psi0);

Psi0_ = psi0./(normFactor);
% Psi0_ = psi0.*(normFactor);

% Psi0_ = 

end



% 2D Normalization
function [Psi0_] = NormalizeWave2(X, Y, psi0)
dx = abs(X(1)-X(2));
dy = abs(Y(1)-Y(2));


trapz(((conj(psi0).*psi0).^2)*dx)

% normFactor = (sum(conj(psi0).*psi0.*SPACE));
% psi0/(1/abs(1-(trapz(dx*(conj(psi0).*psi0)))));
normFactor = (abs(1-(trapz(dx*(conj(psi0).*psi0)))));

Psi0_ = psi0./normFactor;

end





%% Gaussian Wavepacket

% Generate Gaussian
function [psi0_] = Gauss_1(X, x0, k_width, k_0)
CONSTANTS;
% x0 = S_0(1);
dx = abs(X(1)-X(2));
Lx = (X(end)-X(1))/a0;

for x=1:length(X)
psi0_(x) = Lx*exp(-(X(x)-x0).^2/(4*k_width^2)).*exp(1i*k_0*x);

end

end





%% Excited State(s):    1-Dimension
function [Psi0] = ExcitedState_1d(X, x0, Nx, dt)
% if ~isempty(varargin);
%     t_ = varargin{1};
% else
    t_ = 0;
% end

CONSTANTS;
Dim = 1;

% Lengths
dx = abs(X(1)-X(2));
Lx = (X(end)-X(1)); %/a0

% Eigen Energies
for n=1:length(Nx)
    E_(n) = (hPlanck^2*(Nx(n))^2*pi^2)/(2*me*Lx^2);
end

% Wave Vector
Kx = E_*pi/(Lx);

% Wave Function
psi0_ = zeros(length(Nx), length(X));
% psi0_ = linspace(0,0,(length(X)));
% psi0_1x = psi0_;
% psi0_2x = psi0_;

% if length(Nx)>1
%     psi0_2x = linspace(0,0,(length(X)));
% end

% Freq
Wnx = (pi^2*hPlanck*Nx)/(2*Lx^2*me);
% *exp(-1i*Wnx*1)

% Normalization: Simple or Hermitian
A = sqrt((2^Dim)/(Lx)); % 1
% ((pi^-0.25)/(sqrt((2^Enx)*(factorial(Enx)))));
% for n=1:length(Enx)
% HEnx(n) = HermitePolynomial(Enx(n), X./a0);
% end
    
% Operator
% switch length(Nx)
%     case 1
for n=1:length(Nx)

for x=1:length(X)%2:length(X)-1
%     if mod(Nx(n),2)
        psi0_(n,x) = sin(Nx(n)*pi*(X(x)-X(1)+x0)/Lx)*exp(-1i*E_(1)*t_*dt/hPlanck);
%         psi0_1x = A*sin(Enx(1)*pi*(X(x)+x0)/Lx)*exp(-1i*E_(1)*dt/hPlanck);
%     else
%         psi0_(n,x) = cos(Nx(n)*pi*(X(x)-X(1)+x0)/Lx)*exp(-1i*E_(1)*t_*dt/hPlanck);
%     end
    
end

% Boundary Conditions
psi0_(n,1) = 0;
psi0_(n,end) = 0;

end

% Contribution from X component
% if n>1
    psi0 = sqrt(max(Nx))^-1*sum(psi0_,1);
% else
%     psi0 = sum(psi0_,1);
% end

    
    
%     case 2
% for x=1:length(X)%2:length(X)-1
%     if mod(Enx(1),2)
%         psi0_1x(x) = sin(Kx(1)*(X(x)+x0))*exp(-1i*E_(1)*dt/hPlanck);
%         psi0_1x(x) = sin(Nx(1)*pi*x*dx/Lx)*exp(-1i*E_(1)*t_*dt/hPlanck);
%     else
%         psi0_1x(x) = cos(Enx(1)*pi*x*dx/Lx)*exp(-1i*E_(1)*0/hPlanck);
%     end
    
%     if mod(Enx(2),2)
%         psi0_2x(x) = sin(Nx(2)*pi*x*dx/Lx)*exp(-1i*E_(2)*t_*dt/hPlanck);
%     else
%         psi0_2x(x) = cos(Enx(2)*pi*x*dx/Lx)*exp(-1i*E_(2)*0/hPlanck);
%     end
    
    % Contribution from X component
%     psi0_1 = psi0_1x*exp(-1i*E_(1)*dt/hPlanck);
%     psi0_2 = psi0_2x*exp(-1i*E_(2)*dt/hPlanck);
% end
%     psi0_ = (psi0_1x+psi0_2x)./sqrt(2);
% psi0_ = (psi0_1x+psi0_2x);    


% end
        
        
        
% end
% psi0(x) = exp(-(((X(x)-x0)/a0)^2)/2);

% Normalization & Eigenfunctions
% Psi0 = psi0;
Psi0 = NormalizeWave1(X,psi0);



% Psi0 = A*psi0;  
% Psi0 = A*HEnx.*psi0;

% Normalization
%    Psi0 = psi0*(1/abs(1-trapz(dz, dz*trapz(dy,(trapz(dx,dx*(conj(psi0).*psi0),1)),2)),3));
%     
end
    

% 
% %% Excited State(s):    2-Dimensions
% function [Psi0] = ExcitedState_2d(X, x0, Enx, Y, y0, Eny, dt)
% CONSTANTS;
% Dim = 2;
% 
% % Lengths
% dx = abs(X(1)-X(2));
% Lx = (X(end)-X(1)); %/a0
% Ly = (Y(end)-Y(1)); %/a0
% 
% % Wave Vector
% Kx = Enx*pi/(Lx);
% Ky = Eny*pi/(Ly);
% 
% % Wave Function
% psi0 = zeros(length(X), length(Y));
% %meshgrid(1:(length(X)),1:(length(Y)));
% 
% % Freq
% Wnx = (pi^2*hPlanck*Enx)/(2*Lx^2*me);
% Wny = (pi^2*hPlanck*Eny)/(2*Ly^2*me);
% % *exp(-1i*Wnx*1)
% 
% % Normalization: Simple or Hermitian
% A = sqrt((2^Dimension)/(Lx*Ly));
% % ((pi^-0.25)/(sqrt((2^Enx)*(factorial(Enx)))));
% % HEnx = HermitePolynomial(Enx, X./a0);
% % HEny = HermitePolynomial(Eny, Y./a0);
%     
% % Operator
% for x=1:length(X)%2:length(X)-1
% for y=1:length(Y)
%     
%     % Contribution from X component
%     switch mod(Enx,2)
%         case 0
%             psi0_x = sin(Kx*X(x)+x0);
% 
%         otherwise
%             psi0_x = cos(Kx*X(x)+x0);
%     end
%     
%     % Contribution from Y component
%     switch mod(Eny,2)
%         case 0
%             psi0_y = sin(Ky*Y(y)+y0);
% 
%         otherwise
%             psi0_y = cos(Ky*Y(y)+y0);
%     end
%     
%     % Total Wavefunction
%     psi0(x,y) = psi0_x*psi0_y;
%     %     psi0(n) = A*exp(-(((X(n)-x0)/a0)^2))*exp(-1i*Wnx*1); %*(L/a0)    
% end
% end
% % psi0(x) = exp(-(((X(x)-x0)/a0)^2)/2);
% 
% % Normalization & Eigenfunctions
% Psi0 = psi0;  
% % Psi0 = A*psi0;  
% % Psi0 = A*HEnx.*psi0;
% 
% 
% % Normalization
% %    Psi0 = psi0*(1/abs(1-trapz(dz, dz*trapz(dy,(trapz(dx,dx*(conj(psi0).*psi0),1)),2)),3));
% %     
% 
% end
%     
% %% Excited State(s):    2-Dimensions
% function [] = ExcitedState_3d(X, x0, Enx, Y, y0, Eny, Z, z0, Enz)
% CONSTANTS;
% Dim = 2;
% 
% % Lengths
% dx = abs(X(1)-X(2));
% Lx = (X(end)-X(1)); %/a0
% Ly = (Y(end)-Y(1)); %/a0
% Lz = (Z(end)-Z(1)); %/a0
% 
% % Wave Vector
% Kx = Enx*pi/(Lx);
% Ky = Eny*pi/(Ly);
% Kz = Eny*pi/(Lz);
% 
% % Wave Function
% psi0 = zeros(length(X), length(Y), length(Z));
% %meshgrid(1:(length(X)),1:(length(Y)));
% 
% % Freq
% Wnx = (pi^2*hPlanck*Enx)/(2*Lx^2*me);
% Wny = (pi^2*hPlanck*Eny)/(2*Ly^2*me);
% Wnz = (pi^2*hPlanck*Enz)/(2*Lz^2*me);
% % *exp(-1i*Wnx*1)
% 
% % Normalization: Simple or Hermitian
% A = sqrt((2^Dimension)/(Lx*Ly*Lz));
% % ((pi^-0.25)/(sqrt((2^Enx)*(factorial(Enx)))));
% HEnx = HermitePolynomial(Enx, X./a0);
% HEny = HermitePolynomial(Eny, Y./a0);
% HEnz = HermitePolynomial(Enz, Z./a0);
%     
% % Operator
% for x=1:length(X)%2:length(X)-1
% for y=1:length(Y)
% for z=1:length(Z)
%     
%     % Contribution from X component
%     switch mod(Enx,2)
%         case 0
%             psi0_x = sin(Kx*X(x)+x0);
% 
%         case 1
%             psi0_x = cos(Kx*X(x)+x0);
%     end
%     
%     % Contribution from Y component
%     switch mod(Eny,2)
%         case 0
%             psi0_y = sin(Ky*Y(y)+y0);
% 
%         case 1
%             psi0_y = cos(Ky*Y(y)+y0);
%     end
%     
%     % Contribution from Z component
%     switch mod(Enz,2)
%             case 0
%             psi0_z = sin(Kz*Z(z)+z0);
% 
%             case 1
%             psi0_z = cos(Kz*Z(z)+z0);
%     end
%     
%     % Total Wavefunction
%     psi0(x,y,z) = psi0_x*psi0_y*psi0_z;
%     %     psi0(n) = A*exp(-(((X(n)-x0)/a0)^2))*exp(-1i*Wnx*1); %*(L/a0)    
% end
% end
% end
% % psi0(x) = exp(-(((X(x)-x0)/a0)^2)/2);
% 
% % Normalization & Eigenfunctions
% Psi0 = psi0;  
% % Psi0 = A*psi0;  
% % Psi0 = A*HEnx.*psi0;
% 
% 
% % Normalization
% %    Psi0 = psi0*(1/abs(1-trapz(dz, dz*trapz(dy,(trapz(dx,dx*(conj(psi0).*psi0),1)),2)),3));
% %     
% end






%% Hermite Polynomial
function [Hn_] = HermitePolynomial(En, r)
% for now, just store top 10!

Hn_ = 0;

switch En
    case 0
        Hn = 1;
                
    case 1
        Hn = [2 0];
        
    case 2
        Hn = [4 0 2];
        
    case 3
        Hn = [8 0 12 0];
        
    case 4
        Hn = [16 0 -48 0 12];
        
    case 5
        Hn = [32 0 -160 0 120 0];
        
    case 6
        Hn = [64 0 -480 0 720 0 -120];
        
    case 7
        Hn = [128 0 -1344 0 3360 0 -1680 0];
        
    case 8
        Hn = [256 0 -3584 0 13440 0 -13440 0 1680];
        
    case 9
        Hn = [512 0 -9216 0 48384 0 -80640 0 30240 0];
        
    case 10
        Hn = [1024 0 -23040 0 161280 0 -403200 0 302400 0 -30240];
        
    otherwise
        msgID = sprintf('EnergyState:OutOfBounds');
        msgErr = sprintf(['ERROR:\t\ [En > 10]',...
            '\n\t\tEnergy State may not be larger than 10']);
        error(msgID, msgErr);
        
end


for n=1:length(Hn)
    Hn_ = Hn_ + Hn(n).*((r).^(n-1));
end

        
end



%% Eigenstate Setup
function [Psi0] = EigenStateInitialize(PULSE, EigenFunction, Dims, X, x0)
% if nargin>4
%     Psi0_coeff = varargin{1};
% else
%     Psi0_coeff(1:PULSE.NumStates) = 1;
% end
    Psi0_ = zeros(1,length(X));
    CoeffCheck = 0;
    for n=1:PULSE.NumStates
        Psi0_ = Psi0_ + NormalizeWave1(X,EigenFunction(PULSE.EnergyState(n),:));
        
%         CoeffCheck = CoeffCheck + Psi0_coeff(PULSE.EnergyState(n));
    end
    
    % Check for coefficient normalization
%     if abs(sqrt(sum(CoeffCheck.^2))-1)<(.05);
%         Psi0_eig_ = Psi0_;
%     else
% Dims = numel(size(Psi0_eig))-1;
Dims = 1;
%         switch Dims
%     case 1
        [Psi0] = NormalizeWave1(X, Psi0_);
%     case 2
%         [Psi0_eig_] = NormalizeWave2(X, Psi0_);
%     case 3
%         [Psi0_eig_] = NormalizeWave3(X, Psi0_);
%     otherwise
%         end
%     end
end



%% Normalized Wavefunctions with Coefficients Setup
function [Psi0_Norm] = EigenStateInitialize_RAW(PULSE, EigenFunction, EigenCoefficient, Dims, varargin)
switch Dims
    case 1
        X = varargin{1};
        x0 = varargin{2};
        
    case 2
        X = varargin{1};
        x0 = varargin{2};
        Y = varargin{3};
        y0 = varargin{4};
        
    case 3
        X = varargin{1};
        x0 = varargin{2};
        Y = varargin{3};
        y0 = varargin{4};
        Z = varargin{5};
        z0 = varargin{6};
        
    otherwise
        error('SolutionSpace:SpaceDimensionality', ...
        ['ERROR:\tDimensionality of Space not defined.'...
        '\n\n\t\tPlease check definitions and try again.']);
end

    Psi0_ = zeros(1,length(X));
    CoeffCheck = 0;
    for n=1:PULSE.NumStates
        Psi0_ = Psi0_ + NormalizeWave1(X,EigenFunction(n,:)*EigenCoefficient(n));
        
        CoeffCheck = CoeffCheck + EigenCoefficient(PULSE.EnergyState(n));
    end
    
    % Check for coefficient normalization
%     if abs(sqrt(sum(CoeffCheck.^2))-1)<(.05);
%         Psi0_eig_ = Psi0_;
%     else
% Dims = numel(size(Psi0_eig))-1;
Dims = 1;
%         switch Dims
%     case 1
        [Psi0_Norm] = NormalizeWave1(X, Psi0_);
%     case 2
%         [Psi0_eig_] = NormalizeWave2(X, Psi0_);
%     case 3
%         [Psi0_eig_] = NormalizeWave3(X, Psi0_);
%     otherwise
%         end
%     end
end


%% Analytical State
function [Psi0] = AnalyticalState(PULSE, SIM, Dims, dt, varargin)
CONSTANTS;

switch Dims
    case 1
        X = varargin{1};
        x0 = varargin{2};
        
    case 2
        X = varargin{1};
        x0 = varargin{2};
        Y = varargin{3};
        y0 = varargin{4};
        
    case 3
        X = varargin{1};
        x0 = varargin{2};
        Y = varargin{3};
        y0 = varargin{4};
        Z = varargin{5};
        z0 = varargin{6};
        
    otherwise
        error('SolutionSpace:SpaceDimensionality', ...
        ['ERROR:\tDimensionality of Space not defined.'...
        '\n\n\t\tPlease check definitions and try again.']);
end

switch SIM.PotentialMap
    case {'H1'; 'H_1'}
switch Dims
    case 1
        psi0_ = zeros(size(X));
%         psi0 = zeros(size(X));
        psi0 = zeros(PULSE.NumStates, length(X));
        Psi0_ = zeros(size(X));
        R = zeros(size(X));
        %                 F = zeros(size(X));
        %                 P = zeros(size(X));
        l = 0; m = 0;

for N = 1:PULSE.NumStates
    n = PULSE.EnergyState(N);
    switch n
        case 1
        for r=2:length(X)-1
            R(r) = (2/(a0^1.5))...
                *exp(-abs(X(r+x0))/(n*a0));
        end        
        F = 1/sqrt(2*pi);
        P = 1/sqrt(2);
        
        case 2
        for r=2:length(X)-1
            R(r) = (1/(2*sqrt(2)*a0^1.5))...
                *(2-(X(r)/a0))...
                *exp(-abs(X(r+x0))/(n*a0));
        end
        F = 1/sqrt(2*pi);
        P = 1/sqrt(2);

        case 3
        for r=2:length(X)-1
            R(r) = (2/(81*sqrt(3)*a0^1.5))...
                *(27-18*(X(r)/a0)+2*(X(r)/a0)^2)...
                *exp(-abs(X(r+x0))/(n*a0));
        end
        F = 1/sqrt(2*pi);
        P = 1/sqrt(2);
        
        
        otherwise
            error('InitialPulse:EnergyState', ...
            ['ERROR:\tInitial Energy State is invalid.'...
            '\n\n\t\tPlease check definitions and try again.']);
        
    end
        
        for r=2:length(X)-1
%             psi0_(r) = X(r)^2*abs(R(r))^2*F*P;
            psi0_(r) = abs(X(r))*R(r)*F*P;
        end
        
%         psi0(1) = 0; psi0(end) = 0;
        psi0(N,:) = NormalizeWave1(X,psi0_);
        Psi0_ = Psi0_ + psi0(N,:);
end

        
    case 2
        
    
    case 3
    
        
    otherwise
        error('SolutionSpace:SpaceDimensionality', ...
        ['ERROR:\tDimensionality of Space not defined.'...
        '\n\n\t\tPlease check definitions and try again.']);
end
%             X(x)^2;

    case {'Square'}
        Psi0_ = ExcitedState_1d(X, x0, PULSE.EnergyState, dt);
        
        
%**************************************************************************
	case {'QHO'; 'HarmonicOscillator'}
%--------------------------------------------------------------------------
Psi0_ = zeros(size(X));
psi0 = zeros(PULSE.NumStates, length(X));

for n=1:length(PULSE.NumStates)
    Hp = HermitePolynomial(PULSE.EnergyState(n), sqrt(me*e0/hPlanck^2).*X);
    
    psi0(n,:) = sqrt((2^PULSE.EnergyState(n)).*factorial(n)).*...
        ((me*(abs(e0)/hPlanck)/(pi*hPlanck))^(0.25)).*...
        exp(-me*(abs(e0)/hPlanck).*(X.^2)./(2*hPlanck)).*Hp;

    Psi0_ = Psi0_ + psi0;
end    
        
    otherwise
        error('SolutionSpace:PotentialMap_type', ...
        ['ERROR:\tPotential Map Type Undefined.'...
        '\n\n\t\tPlease check definitions and try again.']);
end



% Normalize
switch Dims
    case 1
        [Psi0] = NormalizeWave1(X,Psi0_);
        
    case 2
        
        
    case 3
        
        
    otherwise


end

end
function [varargout] = ITP_EigenFinder(SPACE, TIME, SIM, ITP, Vo, t)
%% Description



%% Sanity Check


%% Local Variables
CONSTANTS;


% ITP_operator = exp(V_t/2).*exp(1/2).*exp(V_t/2);

% a_ = (1-(0.5*dt_i.*V_t/e0))./(1+(0.5*dt_i.*V_t/e0));
% b_ = (1)./(1+(0.5*dt_i.*V_t/e0));


%% Regularize Atomic Potential

% if min(min(min(Vo)))<0
%     Vreg = abs(min(min(min(Vo))))./abs(e0);
% else
%     Vreg = 0;
% end
% V_t = (Vo./abs(e0)+Vreg);
V_t = Vo;
Vreg = 0;

%% Imaginary Time Propagation
%**************************************************************************
% Initial Wave Function
%--------------------------------------------------------------------------
Psi0_ = Psi0_QHO(ITP,SIM, SPACE,TIME);
% Psi0_ = 

% SIM.BC (real)
% Psi0_(1) = 0; Psi0_(end) = 0;


% R_current(1)=0;%R_next(2);
% R_current(SPACE.N_)=0;%R_next(N-1);
% I_current(1)=0;%R_next(2);
% I_current(SPACE.N_)=0;%R_next(N-1);


% SAVE WAVE FUNCTIONS
Psi_ALL = zeros(ITP.NumStates, SPACE.N);
Psi_next = zeros(1, SPACE.N);
% PsiI_ALL = zeros(ITP.NumStates, SPACE.N_);

% DEBUG stuff
% V_t = V_eff(t,:);

%% ITP Setup
%##########################################################################

%**************************************************************************
% Potential --> ITP_EigenFinder
%--------------------------------------------------------------------------
% for t=1:1%length(TIME.Axis)
% V_t = V_eff(t,:);
% for t=TIME.save:TIME.save:TIME.N
%     V_t = V_eff(t,:);
% end

% Propagation Constants
% s1=(TIME.delta)/(SPACE.delta^2);
% s2=(TIME.delta)/(hPlanck);
%0.5*TIME.delta*hPlanck/(me).*V_t)


% k-Space Vector (for Plane Wave Expansion)
% k_ = (1/sqrt(2*pi))*fftaxisshift(fftaxis(SPACE.Axis));

% V_tk = fft(V_t);
% V_ti = ifft(V_tk);

% Imaginary Time Pamateers
dx = (SPACE.length/a0)/SPACE.N;
% dx = (SPACE.length/a0)/SPACE.N_;
dt_i = (1/sqrt(1/dx)^2)*ITP.TimeDivisor;%TIME.delta/Ta;
C_goal = ITP.ConvergenceRatio;

for n=1:ITP.FindStates
fprintf('\t\tITP(i*E-Space)::: [En{%02.0f of %02.0f}]\n\n', n, ITP.FindStates);

    if n>1
%         removeComponent = 0;
        Psi0(n,:) = Psi0(n-1,:) - Coeff_(n-1).*PsiNorm(n-1,:);%.*exp(-EnergyState(n-1));
    else
        Psi0(n,:) = Psi0_;
    end
    
% Boundary COnditions
    Psi_prev = Psi0(n,:);
    Psi_prev(1) = 0; Psi_prev(SPACE.N) = 0;
%     Psi_prev(1) = 0; Psi_prev(SPACE.N_) = 0;
   
% Set ImagTime Step
% dt_i = 1e0*TIME.delta;
% V_ = V_t;
% T_t = abs((1/sqrt(SPACE.length))*ifft((k_.^2/2).*(fft(Psi0(n,:)))));
% T_op = @(Psi) (diff(Psi,2,2).*(hPlanck^2/(2*me)));

% Define ImagTime Evolution Operator
% ITP_operator = ...
%     expm1(-dt_i*V_t).*...
%     expm1(-dt_i*T_t/(1));

%.*...
%     expm1(-dt_i*V_t/(2));
% expm1(dt_i*T_op(Psi_prev)/(1)).*...

%//////////////////////////////////////////////////////////////////////////
% for i_t=1:ITP.TimeSteps(1)

CONVERGENCE = false; En_(1) = 0; ConvStep = 1; flag1 = false; flag2 = false;
while ~CONVERGENCE
    
% Propagation Coeff_icients
% a_ = (1-(0.5*dt_i*V_t))./(1+(0.5*dt_i*V_t));
% b_ = (1./(1+(0.5*dt_i*V_t)));

% a_ = (1-(0.5*dt_i*V_t/(hPlanck^2/me)))./(1+(0.5*dt_i*V_t/(hPlanck^2/me)));
% b_ = (1./(1+(0.5*dt_i*V_t/(hPlanck^2/me))));

a_ = (1-(0.5*dt_i*V_t/e0))./(1+(0.5*dt_i*V_t/e0));
b_ = (1./(1+(0.5*dt_i*V_t/e0)));

%..........................................................................
% REAL Wave Function
% switch mod(ConvStep,2)
%     case 1
%         nN = 1;
%     case 0
%         nN = -1;
% end
Psi_next(2:SPACE.N-1) = a_(2:SPACE.N-1).*Psi_prev(2:SPACE.N-1)...
    +(1)^ConvStep*b_(2:SPACE.N-1).*dt_i/(2*dx^2)...
    .*(Psi_prev(3:SPACE.N)-2*Psi_prev(2:SPACE.N-1)+Psi_prev(1:SPACE.N-2));

switch SIM.BC.Type
    case 'Inf'
Psi_next(1)=0;%R_next(2);
Psi_next(SPACE.N)=0;%R_next(N-1);
    case 'Mask'
%         SIM.BC.params = Mask_param = 1;
        MaskFunc{1} = cos((pi/2)*((SPACE.Axis_(SIM.BC.R{1})/...
            (SPACE.Axis_(SIM.BC.N)-SPACE.Axis_(1)))))...
            .^(SIM.BC.params_.Mask_param);
        MaskFunc{2} = cos((pi/2)*((SPACE.Axis_(SIM.BC.R{2})/...
            (SPACE.Axis_(end)-SPACE.Axis_(end-SIM.BC.N)))))...
            .^(SIM.BC.params_.Mask_param);
        
Psi_next(SIM.BC.R{1}) = MaskFunc{1}.*Psi_next(SIM.BC.R{1});
Psi_next(SIM.BC.R{2}) = MaskFunc{2}.*Psi_next(SIM.BC.R{2});
        
end

% .*dt_i/(2*dx^2)

% (TIME.delta/2*SPACE.delta).*
% ITP_operator.*Psi_prev;%Psi_prev(2:end-1);

% % Boundary conditions
% Psi_next(1)=0;%R_next(2);
% Psi_next(SPACE.N_)=0;%R_next(N-1);


%..........................................................................
% % Display Time
% show_t = (TIME.N/TIME.save);
% if show_t>1e6; show_t = 1e6; end; if show_t<1e4; show_t = 1e4; end
%         if ~mod(i_t, show_t);
% fprintf('\tITP(En{%d} of %d)\t:::\t\t[%d\tof\t%5.0f]\t=\t<%d*i(s)>\n', ...
%     n, ITP.FindStates, i_t, ITP.TimeSteps(1)/n, i_t*TIME.delta);
% % exp(-TIME.delta*V_t/2).*exp(-T_).*exp(-TIME.delta*V_t/2).*PsiR_;
% 
% pdf_ = abs((conj(Psi_next)+Psi_next).*(Psi_next+Psi_next));
% PDF_(n,i_t/show_t,:) = pdf_;
% 
% animateWavefunction(real(Psi_next),imag(Psi_next),pdf_,V_t,...
%     i_t,TIME.delta,SPACE);
% 
%         end
  
% fft(Psi_next);


if n>1
    E0_(n) = E0_(n-1) + dEn(n-1);
%     E0_(n) = dEn(n-1);
%     E0_(n) = EnergyState(n-1);
else
	E0_(n) = 0;
end
% Energy Convergence
En = sum((V_t(2:SPACE.N-1)/e0).*conj(Psi_next(2:SPACE.N-1))...
    .*Psi_next(2:SPACE.N-1)-(Psi_next(2:SPACE.N-1)./((2*dx^2)))...
    .*(Psi_next(3:SPACE.N)-2*Psi_prev(2:SPACE.N-1)+Psi_next(1:SPACE.N-2)))...
    ./(sum((conj(Psi_next(2:SPACE.N-1)).*Psi_next(2:SPACE.N-1))));

% En = sum((V_t(2:SPACE.N_-1)/e0).*conj(Psi_next(2:SPACE.N_-1)).*Psi_next(2:SPACE.N_-1)...
%     - (Psi_next(2:SPACE.N_-1)./((2*dx^2)))...
%     .*(Psi_next(3:SPACE.N_)-2*Psi_prev(2:SPACE.N_-1)+Psi_next(1:SPACE.N_-2)))...
%     ./(sum((conj(Psi_next(2:SPACE.N_-1)).*Psi_next(2:SPACE.N_-1))));


% E0_(n) + 
if isnan(En)
    error('ERROR, Energy must be real number!');
end
% En = dx*trapz((V_t(2:SPACE.N_-1)/e0).*conj(Psi_next(2:SPACE.N_-1)).*Psi_next(2:SPACE.N_-1)...
%     - (Psi_next(2:SPACE.N_-1)./((2*dx^2)))...
%     .*(Psi_next(3:SPACE.N_)-2*Psi_prev(2:SPACE.N_-1)+Psi_next(1:SPACE.N_-2)))...
%     ./(dx*trapz((conj(Psi_next(2:SPACE.N_-1)).*Psi_next(2:SPACE.N_-1))));
      
% ENERGY.En(n) = ...
%     (trapz(V_t.*conj(PsiT).*(PsiT) - (PsiT).*psi_...
%     .*(TIME.delta*hPlanck/(2*me*SPACE.delta^2))))...
%     ./ (trapz(conj(PsiT).*(PsiT)));


% Convergence Check
C_current = abs(1-En/En_);
if C_current<C_goal; CONVERGENCE = true; end

% DISPLACE
% if false
if (C_goal/C_current)>1; C_current = C_goal; end
if ~mod(ConvStep,5e4)||CONVERGENCE;
    fprintf('\tITP(En{%02.0f of %02.0f}--> %2.2d (eV)):::\t[(~%05.2f%%) Converged...]\t@ Time = <%5.4e*(s) + %5.4e*(i*s)>\n', ...
    n, ITP.FindStates, En+E0_(n), 100*C_goal/C_current, t*TIME.delta, ConvStep*dt_i*Ta);

% figure(1); plot(SPACE.Axis, Psi_prev); title(['ITP of \Psi_ ' num2str(n) ' ']);

% if ((C_goal/C_current)>(C_goal/C_past))&&flag1
%     if flag2; CONVERGENCE = true; end
%     flag2 = true;
% end

% if ((C_goal/C_current)<(C_goal/C_past))
%     flag1 = true;
% else
%     flag1 = false;
% end

end
% end


% Timeshift
En_ = En;
Psi_prev=Psi_next;
ConvStep = ConvStep+1;

C_past = C_current;

end
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
% end
PsiT = Psi_next;
Psi_F(n,:) = (PsiT);

%**************************************************************************
% Eigenenergy States
%--------------------------------------------------------------------------
dEn(n) = En;
t_i = ConvStep*dt_i*Ta;

if n>1
% for n_=1:n
%     EnergyState(n) =  -dEn(n-1) + dEn(n);
%     EnergyState(n) =  dEn(n);
    
EnergyState(n) =  sum(dEn(1:n-1)) + dEn(n);
% EnergyState(n-1) + % end
else
%     EnergyState(n) = -dEn(n);
    EnergyState(n) = dEn(n);
    
end
% Estate = -EnergyState - max(EnergyState);
% Estate = -EnergyState;

fprintf('\n\n\n\t***\t[CONVERGENCE]:::(E{%1.0f%s} --> %2.2d eV),  Time = <0*(s) + %5.4e*(i*s)>\n\n\n', ...
    n, ITP.OrbitalType(n), EnergyState(n), ConvStep*dt_i*Ta);


%**************************************************************************
% Normalized Eigenfunction
%--------------------------------------------------------------------------
% NormFactor(n) = dx*sqrt(sum((PsiT).^2))/(SPACE.length/a0);

% NormFactor(n) = dx*trapz(abs(PsiT).^2)/(SPACE.length/a0);

% if n>1
%     norm_ = 0;
%     for n_=1:n-1
%     norm_ = norm_ + (PsiNorm(n_,:).*Psi_F(n,:)).*PsiNorm(n_,:);
%     end
%     NormFactor(n) = sqrt(sum((PsiT - norm_).^2));
%     PsiNorm(n,:) = (PsiT - norm_);
% else
%     NormFactor(n) = sqrt(sum((PsiT-0).^2));
%     PsiNorm(n,:) = (PsiT-0)./NormFactor(n);
% end
 




if n>1
    sumP = 0;
    for m=1:n-1
    sumP = sumP + (PsiNorm(m,:).*(PsiNorm(m,:)).*Psi_F(n,:));
    end
    
    NormFactor(n) = sqrt(sum((PsiT-sumP).^2));
    PsiNorm(n,:) = (PsiT-sumP)./NormFactor(n);
    
else
    NormFactor(n) = sqrt(sum((PsiT-0).^2));
    PsiNorm(n,:) = (PsiT-0)./NormFactor(n);
end

% (dx.*trapz(PsiT.^2)/(SPACE.length/a0))
% trapz(conj(PsiT).*PsiT)*dx;

% PsiNorm(n,:) = (PsiT)./(trapz(conj(PsiT).*(PsiT))/(SPACE.N_));
% PsiNorm(n,:) = (PsiT)./(dx.*trapz(abs(PsiT).^2)/(SPACE.length/a0));
%/SPACE.length);%);

% PsiNorm(n,:) = (PsiT).*(sum(SPACE.delta*(conj(PsiT).*PsiT)/SPACE.length));
% PsiNorm(n,:) = (PsiT)./(cumtrapz((conj(PsiT).*PsiT)));
%/SPACE.length);
% PsiNorm(n,:) = (PsiT)./(trapz(conj(PsiT).*(PsiT)));%/SPACE.length);%);
% PsiNorm(n,:) = (PsiR_next)./abs(trapz(conj(PsiR_next+1i*PsiI_next).*(PsiR_next+1i*PsiI_next)));
% PsiNorm(n,:) = (PsiR_next+1i*PsiI_next)./...
%     trapz(SPACE.delta*conj(PsiR_next+1i*PsiI_next).*(PsiR_next+1i*PsiI_next));

% SPACE.delta*
% Coeff_(n) = sum(conj(PsiT(n,:)).*PsiNorm(n,:));%/SPACE.length;
% Coeff_(n) = 1/trapz(conj(Psi0(n,:)).*PsiNorm(n,:));


%**************************************************************************
% Eigenfunction Coefficient
%--------------------------------------------------------------------------
% Coeff_(n) = dx*sum(conj(Psi0(n,:)).*PsiNorm(n,:))/(SPACE.length/a0);
% Coeff_(n) = sqrt(sum(conj(Psi0(n,:)).*PsiNorm(n,:)));
% Coeff_(n) = dx*trapz(conj(Psi0(n,:)).*PsiNorm(n,:))/(SPACE.length/a0);
% Coeff_(n) = trapz(conj(Psi0(n,:)).*PsiNorm(n,:))/(SPACE.length/a0);
% Coeff_(n) = trapz(conj(Psi0(n,:)).*PsiNorm(n,:));
% Coeff_(n) = trapz(conj(PsiNorm(n,:)).*Psi0(n,:));
% Coeff_(n) = sqrt(sum(PsiNorm(n,:).^2));
Coeff_(n) = (abs(sum(conj(Psi0(n,:)).*PsiNorm(n,:))));


% ENERGY.En(n) = 0;
EigenEnergy(n) = EnergyState(n);
ConvSteps(n) = ConvStep;

% for n=1:ITP.FindStates



% % Display State
% show_t = (TIME.N/TIME.save); if show_t>10000; show_t = 10000; end
%         if ~mod(t*TIME.save, show_t);
% fprintf('\tITP(t --> i*t)\t:::\t\t[%d\tof\t%d]\t=\t<%d*(s)>\n', t_, TIME.N, t_*TIME.delta);
%         end
% end

end


% ITP_Results.dx = dx;
% ITP_Results.Lx = dx*SPACE.N_;

% ITP_Results.dt_i = dt_i;
ImagTimeProp = ConvSteps.*dt_i*Ta;

% ITP_Results.PsiNorm = PsiNorm;
Coefficients = Coeff_/sum(Coeff_);

% for n=1:ITP.FindStates
%     Psi_sln(n,:) = Coefficients(n).*PsiNorm(n,:)...
%         .*exp(-EigenEnergy(n)*ImagTimeProp(n));
% end

Psi_sln = PsiNorm;




%% OUTPUT
% varargout{1} = ITP_Results;
varargout{1} = PsiNorm;
varargout{2} = Coefficients;
varargout{3} = EigenEnergy;
varargout{4} = ImagTimeProp;
varargout{5} = Psi_sln;


% if nargout > 1
%     varargout{2} = Psi_sln;
% end












end



% switch ITP.FDM
%     
%     case 'FDTD-Q'
% %         FDTD_Q();
%         
%     case 'Gram-Schmidt'
% %         GramSchmidt
%         
%     case 'Lowdin'
%         
%         
%     otherwise
%     error('EnergyStateProperties:ImaginaryTimePropagation_FDM', ...
%         ['ERROR:\tImaginary Time Propagation Eigenfinder Method incorrect.'...
%         '\n\n\t\tPlease check definitions and try again.']);
% end
  














%% ttt


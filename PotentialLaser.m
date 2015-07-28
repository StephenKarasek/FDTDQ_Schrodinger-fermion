function [Laser_E] = PotentialLaser(LASER, SPACE, TIME, t, varargin)

% 


%% Sanity Check

% Input Arguments
narginchk(4,5);

% Output Arguments
% nargoutchk(0);

% %%
% if nargin == 5
%     t0 = 0;
% end
% 
% if nargin == 6
%     t0 = varargin{1};
% end
% 
% if nargin == 7
%     t1(1) = varargin{1};
%     t2(2) = varargin{2};
% end

    


%% Local variable assignment

CONSTANTS;


%% 

switch LASER.Type
    case {'pulse'; 'Pulse'; 'CEP'; 'CarrierEnvelopePulse'; 'ChirpedPulse'; 'CPA'}
        t_ = t*TIME.delta;
        nPulse = floor((t_-LASER.Delay)/LASER.Time);
        if nPulse<0; nPulse=0; end
        nt_0 = LASER.Time*nPulse+LASER.Delay;
        nt_1 = LASER.Time*nPulse+LASER.PulseLength+LASER.Delay;
        t0 = t_-nt_0;
        
        if (t_>nt_0)&&(t_<=nt_1)
            E_ = LASER.E0;
            %         A_ = LASER.A0;
        else
            E_ = 0;
            %         A_ = 0;
        end
        
        E_ = sin(LASER.w0*t_/(2*LASER.Ncycles))^2;
        Laser_E = LASER.E0*E_*cos(LASER.w0*t*TIME.delta+LASER.CEP);
        
        
        
        
    case {'CW'; 'continuous'; 'Continuous-Wave'; 'ContinuousWave'}
        t_ = t*TIME.delta;
%         nPulse = floor((t_-LASER.Delay)/LASER.Time)
        
        
        
        CycleSTART = LASER.Delay + floor((t_+LASER.Phase)/LASER.T)*LASER.T;%...
%             + LASER.Delay;
        CycleFINISH = LASER.Delay + floor((t_+LASER.Phase+LASER.PERIODcutoff)/LASER.T);%*LASER.T...
%             + LASER.Delay;
        
%         Laser_E = LASER.E0*cos(LASER.w0*t*TIME.delta+LASER.Phase);
%         Laser_E = LASER.E0*cos(LASER.w0*(t*TIME.delta + LASER.Phase - LASER.Delay));        

Laser_E = LASER.E0*sin(LASER.w0*(t*TIME.delta + LASER.Phase - LASER.Delay));        
% Laser_E = LASER.E0*cos(LASER.w0*(t*TIME.delta + LASER.Phase - LASER.Delay));        


        if (t_<CycleSTART)&&(t_>CycleFINISH)
            Laser_E = 0;
        end
            
        if (LASER.FINALcutoff~=0)&&(t_>LASER.FINALcutoff+LASER.Delay)
            Laser_E = 0;
        end
       
        if (t_<LASER.Delay)
            Laser_E = 0;
        end
        
        
    otherwise
        %                 E_ = 0;
        Laser_E = 0;
end
% evalin('caller', ['Laser_E = ' num2str(Laser_E) ';']);


%%

% switch GAUGE.Type

%     case 'coulomb'
        
%         Laser_coupling = ['e0.*SPACE.Axis.*' num2str(Laser_E)];
%                 'Laser_coupling = e0.*SPACE.Axis.*Laser_E';
%         evalin('caller',['V_laser = e0.*SPACE.Axis.*Laser_E;']);
        
%         
%     case 'lorenz'
%         Laser_E = 0;
%     evalin('caller',['V_laser = 0*Laser_E;']);

%     otherwise
%         Laser_E = 0;
% end



%% RETURN
varname = 'V_laser';

if nargin > 5
    t_delay = varargin{1};
evalin('caller',['E_l(t/' num2str(t_delay) ') = ' num2str(Laser_E) ';']);    
else
% evalin('caller',['V_laser = ' Laser_coupling ';']);
end







end
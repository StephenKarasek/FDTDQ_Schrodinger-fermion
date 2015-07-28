function [varargout] = animateWavefunction(R_current, I_current, prob_density, Veff, Vmin, Vmax,WaveMax,PdfMax, t, dt, SPACE)
% pie = animateWavefunction();

%% I/O
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%|||||||||||||||||||||||||||||   INPUTS   |||||||||||||||||||||||||||||||||
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%
% VarName:          'Psi_Real
% 
% VarIO
% 
% MemIndex
% 
% 
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%|||||||||||||||||||||||||||||   OUTPUTS   ||||||||||||||||||||||||||||||||
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%
% VarName
% 
% VarIO
% 
% MemIndex
% 
% 
% 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%% Sanity Check

% Input Arguments
narginchk(5,15);
% narginchk(10,10);

% Output Arguments
nargoutchk(0,1);




%% Local variable assignment
% n=1;
% R_current       =   varargin{n}; n=n+1;
% I_current       =   varargin{n}; n=n+1;
% prob_density    =	varargin{n}; n=n+1;
% V               =	varargin{n}; n=n+1;
% t               =   varargin{n}; n=n+1;
% dt              =   varargin{n}; n=n+1;
% x               =   varargin{n}; n=n+1;
% 
% switch nargin
%     case 7
%            
%     case 8
%         y       =   varargin{n}; n=n+1;
%     
%     case 9
%         y       =   varargin{n}; n=n+1;
%         z       =   varargin{n}; n=n+1;
%         
%     otherwise
%         error('Please check inputs and try again.');
% end


% Atomic Constants
CONSTANTS;



%% Locals

Scale_Factor = 2.5;

% Vmin = min(min(min(Veff)));
% Vmin = max(max(max(abs(Veff))));
% Wavemin = sum(max(max(max((R_current+1i*I_current).*(R_current-1i*I_current)))));
% Wavemin = Scale_Factor*sqrt(sum(((R_current+1i*I_current).*(R_current-1i*I_current)).^2));
% Wavemin = Scale_Factor*sqrt(sum(abs((R_current+1i*I_current).^2)));
% Wavemin = Scale_Factor*sqrt(abs(sum((R_current+1i*I_current).^2)));

% PdfMax = min(min(min(prob_density)));
% PdfMax = Scale_Factor*(Wavemin/Scale_Factor)^2;

% PdfMax = (WaveMax);


%% Animate

try
    if ~fig0; fig0 = figure(1); end
catch
    fig0 = figure(1);
    set(fig0,'Renderer','zbuffer');
end
set(fig0, 'Position', [400 100 1120 640]);
pause (0.001);
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
% Real & Imaginary Wavefunctions
% if debug.PLOT.Wave_RI
    subplot(3,1,1);
%     subplot(3,3,[1 2 3]);
    haxes = gca;%get(gcf,'CurrentAxes');
    wave_parts = plot(SPACE.Axis_, I_current, SPACE.Axis_, R_current);
    set(wave_parts(1), 'Color', [0 0 1], 'LineWidth', 1);
    set(wave_parts(2), 'Color', [1 0 0], 'LineWidth', 1);
    
    set(haxes,'FontAngle', 'italic');
%     axis([SPACE.Axis_(1) SPACE.Axis_(end) -Wavemin Wavemin]);
    axis([SPACE.Axis_(1) SPACE.Axis_(end) -WaveMax WaveMax]);

% SCALE_factor = 1e3;
%     axis([SPACE.Axis_(1) SPACE.Axis_(end) -SCALE_factor SCALE_factor]);
    
title(['Time = [' num2str(t*dt, '%10.4e\n') ' (sec)]'], ...
    'FontWeight','bold', 'FontSize',16, 'FontAngle', 'normal');

xlabel('','VerticalAlignment','cap',...
    'HorizontalAlignment','center',...
    'FontWeight','normal', 'FontSize',14, 'FontAngle', 'normal');

ylabel('\psi_I_m_a_g, \psi_R_e_a_l','VerticalAlignment','bottom',...
    'HorizontalAlignment','center',...
    'FontWeight','normal', 'FontSize',14, 'FontAngle', 'normal');

%     legend1 = legend(wave_parts,['\psi_I_m_a_g'], ['\Psi_R_e_a_l']);
%     set(legend1,'Orientation','vertical',  'Location', 'NorthEast');
%     pos = get(legend1, 'OuterPosition');
%     pos(3) = 0.05;
%     set(legend1,'OuterPosition',[pos(1) pos(2) 0.05 pos(4)]);
%     set(legend1,'OuterPosition',[pos(1) pos(2) 0.05 pos(4)], 'Location', 'NorthEastOutside');
%     'horizontal');
%     hline1 = 
% end


%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
% Probability Density Function
% if debug.PLOT.PDF
    subplot(3,1,2);
%     subplot(3,3,[4 5 6]);
    haxes = gca;%get(gcf,'CurrentAxes');
    total_wave_plot = area(SPACE.Axis_, prob_density);
    set(total_wave_plot, 'FaceColor', [0.5 0 0.5], 'LineWidth', 1);
    
    set(haxes,'FontAngle', 'italic');
    axis([SPACE.Axis_(1) SPACE.Axis_(end) -0.01*PdfMax PdfMax]);


%     axis([SPACE.Axis_(1) SPACE.Axis_(end) -SCALE_factor^2/100 SCALE_factor^2]);
    
%     title(['PDF @ time = [' num2str(t*dt) ' (sec)]'], ...
%     'FontWeight','bold', 'FontSize',16, 'FontAngle', 'normal');

xlabel('','VerticalAlignment','cap',...
    'HorizontalAlignment','center',...
    'FontWeight','normal', 'FontSize',14, 'FontAngle', 'normal');

ylabel('|\Psi^*\Psi| ^2','VerticalAlignment','bottom',...
    'HorizontalAlignment','center',...
    'FontWeight','normal', 'FontSize',14, 'FontAngle', 'normal');
    
%         legend2 = legend(total_wave,['|\Psi^*\Psi|^2']);
%     set(legend2,'Orientation','vertical', 'Location', 'NorthEast');
%     pos = get(legend2, 'OuterPosition');
%     pos(3) = 0.05;
%     set(legend2,'OuterPosition',[pos(1) pos(2) 0.05 pos(4)], 'Location', 'NorthEastOutside');
%     hline1 = 
% end

%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
% Potential
% if debug.PLOT.Potential
    subplot(3,1,3);
    %     subplot(3,3, [7 8 9]);
    haxes = gca;%get(gcf,'CurrentAxes');
    plot_potential = area(SPACE.Axis_, Veff/abs(e0), 'BaseValue',Vmin/abs(e0));
    set(plot_potential, 'FaceColor', [0 0.7 0], 'LineWidth', 1);
    if Vmin==0; Vmin=e0; end
    set(gca,'FontAngle', 'italic', 'XGrid', 'on', 'YGrid', 'on');
    axis([SPACE.Axis_(1) SPACE.Axis_(end) Vmin/abs(e0) Vmax/abs(e0)]);
    
% title('Effective Potential', ...
%     'FontWeight','bold', 'FontSize',16, 'FontAngle', 'normal');

xlabel('X-Distance (m)','VerticalAlignment','cap',...
    'HorizontalAlignment','center',...
    'FontWeight','normal', 'FontSize',14, 'FontAngle', 'normal');

ylabel('V_e_f_f (eV)','VerticalAlignment','bottom',...
    'HorizontalAlignment','center',...
    'FontWeight','normal', 'FontSize',14, 'FontAngle', 'normal');

%     legend3 = legend(plot_potential, 'V_e_f_f');
%     set(legend3,'Orientation','vertical',  'Location', 'NorthEast');
%     pos = get(legend3, 'OuterPosition');
%     pos(3) = 0.05;
%     set(legend3,'OuterPosition',[pos(1) pos(2) 0.05 pos(4)]);
%     set(legend3,'OuterPosition',[pos(1) pos(2) 0.06 pos(4)], 'Location', 'NorthEastOutside');
%     'horizontal');
%     hline1 = 
% end


%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
drawnow;




%% Return
varargout{1} = fig0;






end
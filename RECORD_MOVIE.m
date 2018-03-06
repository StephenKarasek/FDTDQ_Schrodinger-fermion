function [DATA] = RECORD_MOVIE(SPACE, TIME, SIM, PULSE, ITP, LASER, DATA, DEBUG)
%% RECORD_MOVIE


close all hidden;

CONSTANTS;


%% Load in vars


% Psi_REAL = real(Psi_TOTAL);
% Psi_IMAG = imag(Psi_TOTAL);

% PDF = (conj(Psi_REAL+Psi_IMAG*1i).*(Psi_REAL+Psi_IMAG*1i));
% PDF = (Psi_REAL.^2.*(1i*Psi_IMAG).^2);


    load([DATA.PARAMS_.saveDir '/PARAMS_']);
    load([DATA.PARAMS_.saveDir '/PDF_pos'], 'PDF_pos');
    load([DATA.PARAMS_.saveDir '/Psi_REAL'], 'Psi_REAL');
    load([DATA.PARAMS_.saveDir '/Psi_IMAG'], 'Psi_IMAG');
    load([DATA.PARAMS_.saveDir '/V_pot'], 'V_pot');

% PDF_pos = (conj(Psi_REAL+Psi_IMAG*1i).*(Psi_REAL+Psi_IMAG*1i));


DATA = MaxPotential(SPACE, TIME, SIM, PULSE, ITP, LASER, DATA, DEBUG);

% Vsys = DATA.V_.Vsys;
Vmax = DATA.V_.Vmax;
Vmin = DATA.V_.Vmin;


WaveMax = sqrt(sum(abs(sqrt(PDF_pos(1,:)))/2)/SPACE.N);
PdfMax = (WaveMax/2)^2;


% Vmin = V_min;
% Vmax = V_max;
% Vmin = min(min(min(V_pot)));
% Vmax = abs(Vmin);
%  Psi_REAL = Psi_REAL/5;
%  Psi_IMAG = Psi_IMAG/5;


%% Determine Name
    
txt_pot = ['V(' SIM.PotentialMap ')'];
% txt_pot = ['.          V(' SIM.PotentialMap ')'];

switch PULSE.Type
    case {'Sine'}
%         txt_pulse = ['_.     \Psi(E:{' num2str(PULSE.EnergyState) '}'];
        txt_pulse = ['\Psi(E:' num2str(PULSE.EnergyState) ')'];

    case 'Gaussian'
%         txt_pulse = ['_.     \Psi(\sigma_' num2str(PULSE.EnergyState) '_-_5)'];
        txt_pulse = ['\Psi(\sigma_' num2str(PULSE.EnergyState) '_-_5)'];
        
    otherwise
%         txt_pulse = ['_.          \Psi(\sigma_' num2str(PULSE.EnergyState) '_-_5)'];
        txt_pulse = ['\Psi(\sigma_' num2str(PULSE.EnergyState) '_-_5)'];
end
% txt_EnPot = {txt_pulse, txt_pot};
txt_EnPot = [txt_pulse ', ' txt_pot];


%//////////////////////////////////////////////////////////////////////////
switch PULSE.Type
    case {'Cosine', 'Sine'}
        Estate_ = ['Sin(n=' num2str(PULSE.EnergyState(1))];
        if length(PULSE.EnergyState)>1
            Estate_ = [Estate_ ',' num2str(PULSE.EnergyState(2)) ];
        end
        Estate_ = [Estate_ ')'];
        
    case 'Gaussian'
        Estate_ = ['Gauss[x_=' num2str(PULSE.Width,'%1.1e') ...
            ', k_=' num2str(PULSE.Momentum,'%1.1e') ,']'];
                
    case {'EigenState'; 'Eigenstate'}
        Estate_ = ['Eigen(n=' num2str(PULSE.EnergyState(1))];
        if length(PULSE.EnergyState)>1
            Estate_ = [Estate_ ',' num2str(PULSE.EnergyState(2)) ];
        end
        Estate_ = [Estate_ ')'];
                        
    case {'Exact'; 'Analytical'}
        Estate_ = ['Exact(n=' num2str(PULSE.EnergyState(1))];
        if length(PULSE.EnergyState)>1
            Estate_ = [Estate_ ',' num2str(PULSE.EnergyState(2)) ];
        end
        Estate_ = [Estate_ ')'];
                        
end

%//////////////////////////////////////////////////////////////////////////
r0_ = '';
for i=1:length(PULSE.InitialPos)
%     component_ = fprintf('  %s  ', num2str(PULSE.InitialPos(i)))
    r0_ = [r0_ ' ' num2str(PULSE.InitialPos(i))];
end

R0_ = ['X_=(' num2str(SPACE.N, '%2.1e') 'dx.' num2str(SPACE.length, '%2.1e') ...
    '), T_=(' num2str(TIME.N, '%2.1e') '), V=' SIM.PotentialMap ')'];



%//////////////////////////////////////////////////////////////////////////
switch LASER.Type
    case 'CW'
        LaserInfo_ = ['Laser[' LASER.Type ',' num2str(LASER.l0, '%2.1e') 'nm' ...
            ', I0=' num2str(LASER.I0, '%2.1e') ']'];
        
    case 'pulse'
        LaserInfo_ = ['Laser[' LASER.Type ',' num2str(LASER.EnvelopeWave, '%2.1e') ','...
            'Ncyc=' num2str(LASER.Ncycles, '%1.1g') ', I0=' num2str(LASER.I0*1e4, '%2.1e') ']'];
                
    otherwise
        LaserInfo_ = ['NoLaser'];
                
end

I0_ = ['I-' num2str(LASER.I0*1e4, '%1.1e')]; 

switch LASER.Type
%//////////////////////////////////////////////////////////////////////////
    case 'pulse'
for t=TIME.save:TIME.save:TIME.N
    Laser_E(t/TIME.save) = LASER.E0*sin(LASER.w0*t*TIME.delta/(2*LASER.Ncycles))^2 ...
    *cos(LASER.w0*t*TIME.delta+LASER.CEP);
    
end
        for t=TIME.save:TIME.save:TIME.N
    Laser_E(t/TIME.save) = LASER.E0*cos(LASER.w0*t*TIME.delta+LASER.Phase);
        end
        txt_Ltype = [LASER.Type '(\lambda_C_E_P=' num2str(LASER.EnvelopeWave, '%2.1e')...
            ', N_c_y_c=' num2str(LASER.Ncycles) ')'];
%                 txt_Ltype = ['_.     ' LASER.Type '(\lambda_C_E_P=' num2str(LASER.EnvelopeWave, '%2.1e')...
%             ', N_c_y_c=' num2str(LASER.Ncycles) ')'];
        txt_Lint = ['Laser:(I_0= ' num2str(LASER.I0*1e4, '%2.1e') ')'];
%     ttxt = { txt_Lint, txt_Ltype};
txt_Laser = [txt_Lint ', ' txt_Ltype];
%//////////////////////////////////////////////////////////////////////////
    case 'CW'
        for t=TIME.save:TIME.save:TIME.N
    Laser_E(t/TIME.save) = LASER.E0*cos(LASER.w0*t*TIME.delta+LASER.Phase);
        end
        txt_Ltype = [LASER.Type '(\lambda= ' num2str(LASER.l0, '%2.1e') ')'];
%         txt_Ltype = ['_.     ' LASER.Type '(\lambda= ' num2str(LASER.l0, '%2.1e') ')'];
        txt_Lint = ['Laser:(I_0= ' num2str(LASER.I0*1e4, '%2.1e') ')'];
%     ttxt = { txt_Lint, txt_Ltype};
    txt_Laser = [txt_Lint ', ' txt_Ltype];
        
    otherwise
%//////////////////////////////////////////////////////////////////////////
end

%/////////////////////////////////////////////////////////////////////////
BC_ = sprintf([SIM.BC.Type ',r=(' num2str(SIM.BC.N/SPACE.N,'%02.2f') ')']);


if DEBUG.DetailedFileName
    SaveVideoName = [LaserInfo_ ',' Estate_ ',' R0_ ',' BC_ ' .png'];
else
%     PARAMS_.saveVisual;
% DATA.PARAMS_.saveDir '\' ...
   SaveVideoName = ...
       ['Vid[' SIM.PotentialMap ';' num2str(LASER.I0, '%2.1e') ...
       '(W.cm-2);' num2str(LASER.l0, '%2.1e') '(nm)]'];
%    ['SIM[V(' SIM.PotentialMap '),Io(' num2str(LASER.I0, '%2.1e') ...
%        '),wL(' num2str(LASER.l0, '%2.1e') ')]'];
end





%% Record Movie
% Movie Parameters
mov_length = length((1:TIME.save:TIME.N+1));
mov(mov_length) = struct('cdata', [], 'colormap', []);

% Display Resolution
OPTIONS.DisplayResolution = [720 540];
OPTIONS.DisplayResolution_ = OPTIONS.DisplayResolution*1;
% OPTIONS.DisplayResolution_ = OPTIONS.DisplayResolution*1.5;
% OPTIONS.DisplayResolution_ = OPTIONS.DisplayResolution*2;

% Compression Type
OPTIONS.CompressionType_ = 'Motion JPEG AVI';

% File Name
timestamp = datestr(clock);
timestamp = regexprep(timestamp, '(\s)', '_');
timestamp = regexprep(timestamp, '([:])', '-');
MovieFilePathName = [DATA.PARAMS_.saveDir '\' '(' timestamp ')_' SaveVideoName '.avi'];
% MovieFilePathName = [DATA.PARAMS_.saveDir '\' SaveVideoName '.avi'];




% if DEBUG.VIDEO.save
% I/O
writerObj = VideoWriter(MovieFilePathName, [OPTIONS.CompressionType_]);
%         writerObj = VideoWriter(MovieFileName, 'Uncompressed AVI');
writerObj.FrameRate = 20;
open(writerObj);


fig0 = figure(1); show = 1;
set(fig0, 'Position', [100 100 OPTIONS.DisplayResolution_]);

%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
% Real & Imaginary Wavefunctions
% if DEBUG.PLOT.Wave_RI   
    subplot(3,1,1);
%     subplot(6,1,[1 2]);
    WAVEPLOT.haxes = gca;%get(gcf,'CurrentAxes');
    WAVEPLOT.wave_parts = plot(WAVEPLOT.haxes, SPACE.Axis_, Psi_IMAG(1,:), SPACE.Axis_, Psi_REAL(1,:));
%     WAVEPLOT.wave_parts = plot(WAVEPLOT.haxes, SPACE.Axis, Psi_IMAG(1,:), SPACE.Axis, Psi_REAL(1,:));
    set(WAVEPLOT.wave_parts(1), 'Color', [0 0 1], 'LineWidth', 1);
    set(WAVEPLOT.wave_parts(2), 'Color', [1 0 0], 'LineWidth', 1);
    
%     Rmax = ceil(max(max(max(max(abs(Psi_IMAG)))))*1e4)/1e4;
%     Imax = ceil(max(max(max(max(abs(Psi_IMAG)))))*1e4)/1e4;
%     if Rmax>Imax
%         WaveMax = 2*Rmax;
%     else
%         WaveMax = 2*Imax;
%     end
    
        
    
    set(WAVEPLOT.haxes,'FontAngle', 'italic');
    axis([SPACE.Axis(1) SPACE.Axis(end) -WaveMax WaveMax]);
%     axis([SPACE.Axis(1) SPACE.Axis(end) -1 1]);
    
% title(['Time = [' num2str(1*TIME.delta) ' (sec)]'], ...
%     'FontWeight','bold', 'FontSize',16, 'FontAngle', 'normal');

% xlabel('Distance (m)','VerticalAlignment','cap',...
%     'HorizontalAlignment','center',...
%     'FontWeight','normal', 'FontSize',14, 'FontAngle', 'normal');

% ylabel('\psi_I_m_a_g, \psi_R_e_a_l','VerticalAlignment','bottom',...
%     'HorizontalAlignment','center',...
%     'FontWeight','normal', 'FontSize',14, 'FontAngle', 'normal');


% end

%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
% Probability Density Function
% if DEBUG.PLOT.PDF
    subplot(3,1,2);
%     subplot(6,1,[3 4]);
    PDFPLOT.haxes = gca;%get(gcf,'CurrentAxes');
    PDFPLOT.total_wave = area(PDFPLOT.haxes,SPACE.Axis_, PDF_pos(1,:));
    set(PDFPLOT.total_wave, 'FaceColor', [0.5 0 0.5], 'LineWidth', 1);
    

%     PdfMax_a = max([Rmax Imax])^2;
%     PdfMax_e = ceil(max(max(max(max(abs(PDF_pos))))));
%     if (PdfMax_e)<PdfMax_a
%         PdfMax = 0.5*PdfMax_e;
%     else
%         PdfMax = 0.5*PdfMax_a;
%     end
   
    set(PDFPLOT.haxes,'FontAngle', 'italic');
    axis([SPACE.Axis(1) SPACE.Axis(end) -0.02*PdfMax PdfMax]);
    
%     axis([SPACE.Axis(1) SPACE.Axis(end) 0 1]);
    
%     title(['PDF @ time = [' num2str(time_step*TIME.delta) ' (sec)]'], ...
%     'FontWeight','bold', 'FontSize',16, 'FontAngle', 'normal');

% xlabel('','VerticalAlignment','cap',...
%     'HorizontalAlignment','center',...
%     'FontWeight','normal', 'FontSize',14, 'FontAngle', 'normal');

% ylabel('|\Psi^*\Psi| ^2','VerticalAlignment','bottom',...
%     'HorizontalAlignment','center',...
%     'FontWeight','normal', 'FontSize',14, 'FontAngle', 'normal');
    
%         legend2 = legend(total_wave,['|\Psi^*\Psi|^2']);
%     set(legend2,'Orientation','vertical', 'Location', 'NorthEast');
%     pos = get(legend2, 'OuterPosition');
%     pos(3) = 0.05;
%     set(legend2,'OuterPosition',[pos(1) pos(2) 0.05 pos(4)], 'Location', 'NorthEastOutside');
%     hline1 = 
% end

%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
% Potential
% if DEBUG.PLOT.Potential
    subplot(3,1,3);
    %     subplot(6,1, [5 6]);
    POTENTIALPLOT.haxes = gca;%get(gcf,'CurrentAxes');
    POTENTIALPLOT.potential = area(SPACE.Axis_, V_pot(1,:)/abs(e0), 'BaseValue',Vmin/abs(e0));
    set(POTENTIALPLOT.potential, 'FaceColor', [0 0.7 0], 'LineWidth', 1);
    
     
    XLim_ = [SPACE.Axis(1) SPACE.Axis(end)];
YLim_pot = [-abs(0.5*Vmin/e0) abs(0.25*abs(Vmin)/e0)];
% YLim_pot = [-abs(0.5*Vmin/e0) abs(0.5*abs(Vmin)/e0)];
% YLim_pot = [-abs(0.5*Vmin/e0) abs(0.5*Vmax/e0)];
    
    set(POTENTIALPLOT.haxes,'FontAngle', 'italic', 'XGrid', 'on', 'YGrid', 'on');
    axis([XLim_ YLim_pot]);
    
% title('Effective Potential', ...
%     'FontWeight','bold', 'FontSize',16, 'FontAngle', 'normal');

% xlabel('X-Distance (m)','VerticalAlignment','cap',...
%     'HorizontalAlignment','center',...
%     'FontWeight','normal', 'FontSize',14, 'FontAngle', 'normal');

% ylabel('V_e_f_f (eV)','VerticalAlignment','bottom',...
%     'HorizontalAlignment','center',...
%     'FontWeight','normal', 'FontSize',14, 'FontAngle', 'normal');

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



set(WAVEPLOT.haxes, 'nextplot', 'replacechildren');
set(PDFPLOT.haxes, 'nextplot', 'replacechildren');
set(POTENTIALPLOT.haxes, 'nextplot', 'replacechildren');
set(fig0,'Renderer','zbuffer');
 for time_step = 1:TIME.saveNum;
%--------------------------------------------------------------------------
    cla;
    WAVEPLOT.wave_parts = plot(WAVEPLOT.haxes, SPACE.Axis_, Psi_IMAG(time_step,:), SPACE.Axis_, Psi_REAL(time_step,:));
    set(WAVEPLOT.wave_parts(1), 'Color', [0 0 1], 'LineWidth', 1);
    set(WAVEPLOT.wave_parts(2), 'Color', [1 0 0], 'LineWidth', 1);
    set(WAVEPLOT.haxes,'FontAngle', 'italic', 'XGrid', 'on', 'YGrid', 'on');
    set(WAVEPLOT.haxes,'YLim', [-WaveMax WaveMax],'XLim', [SPACE.Axis_(1) SPACE.Axis_(end)] );
%     set(WAVEPLOT.haxes,'YLim', [-1 1],'XLim', [SPACE.Axis(1) SPACE.Axis(end)] );
    
    set(get(WAVEPLOT.haxes, 'Title'), 'String', ...
        ['Time = [' num2str(time_step*TIME.save*TIME.delta, '%10.4e\n') ' (sec)]'], ...
        'FontWeight','bold', 'FontSize',16, 'FontAngle', 'normal');
    
    set(get(WAVEPLOT.haxes, 'Ylabel'), 'String', '\psi_I_m_a_g, \psi_R_e_a_l',...
        'VerticalAlignment','bottom','HorizontalAlignment','center',...
        'FontWeight','bold', 'FontSize',14, 'FontAngle', 'normal');
    

%--------------------------------------------------------------------------
    cla;
    PDFPLOT.total_wave = area(PDFPLOT.haxes,SPACE.Axis_, PDF_pos(time_step,:));
    set(PDFPLOT.total_wave, 'FaceColor', [0.5 0 0.5], 'LineWidth', 1);
    set(PDFPLOT.haxes,'FontAngle', 'italic', 'XGrid', 'on', 'YGrid', 'on');
    set(PDFPLOT.haxes,'YLim', [0 PdfMax],'XLim', [SPACE.Axis_(1) SPACE.Axis_(end)] );
%     set(PDFPLOT.haxes,'YLim', [0 1],'XLim', [SPACE.Axis(1) SPACE.Axis(end)] );
    
    set(get(PDFPLOT.haxes, 'Ylabel'), 'String', '|\Psi^*\Psi| ^2',...
        'VerticalAlignment','bottom','HorizontalAlignment','center',...
        'FontWeight','bold', 'FontSize',14, 'FontAngle', 'normal');
    
%--------------------------------------------------------------------------
    cla;
    POTENTIALPLOT.potential = ...
    area(POTENTIALPLOT.haxes,SPACE.Axis_, V_pot(time_step,:)/abs(e0),...
    'BaseValue',-abs(Vmin/e0));
% area(POTENTIALPLOT.haxes,SPACE.Axis, V_pot(time_step,:)/abs(e0),...
%     'BaseValue',min(V_pot(time_step,:))/abs(e0));

% XLim_ = [SPACE.Axis(1) SPACE.Axis(end)];
% YLim_pot = [-abs(0.05*Vmin/e0) abs(0.05*Vmax/e0)];

    set(POTENTIALPLOT.potential, 'FaceColor', [0 0.7 0], 'LineWidth', 1);
    set(POTENTIALPLOT.haxes,'FontAngle', 'italic', 'XGrid', 'on', 'YGrid', 'on');
%     set(POTENTIALPLOT.haxes,'XLim', [x(1) x(end)],...
%         'YLim', [min(V_pot(time_step,:))/abs(e0) abs(max(V_pot(time_step,:))/abs(e0))]);
    set(POTENTIALPLOT.haxes,'XLim', XLim_,'YLim',YLim_pot);

    
    set(get(POTENTIALPLOT.haxes, 'Ylabel'), 'String', 'V_e_f_f (eV)',...
        'VerticalAlignment','bottom','HorizontalAlignment','center',...
        'FontWeight','bold', 'FontSize',12, 'FontAngle', 'normal');

    set(get(POTENTIALPLOT.haxes, 'Xlabel'), 'String', 'X-Distance (m)',...
        'VerticalAlignment','top', 'HorizontalAlignment','center',...
        'FontWeight','bold', 'FontSize',14, 'FontAngle', 'normal');
%--------------------------------------------------------------------------

% drawnow;
pause(0.001);

try
    frame = getframe(fig0);
    writeVideo(writerObj,frame);
catch
    warning('MOVIE NOT SAVED!!!');
   close(writerObj);
end

end



close(writerObj)



%%
close(fig0);



end

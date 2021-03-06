function GRAPH_EnergyOutput_VS_Laser(TEXT_, SIM, taxis, xaxis)
%%
close all hidden;
CONSTANTS;
load([TEXT_.SaveDirectory '/En_REAL'], 'En_REAL');
load([TEXT_.SaveDirectory '/E_laser'], 'E_laser');

%%

EoutLaser_fig = figure;
subplot(3,3,4:9);
% sqrt(trapz(TIME.conj(En_REAL+1i*En_IMAG).*(En_REAL+1i*En_IMAG)))
EoutLaser_plt = plot(taxis, En_REAL/(e0), taxis, zeros(size(En_REAL)));


set(EoutLaser_plt(1), 'color', [1 0.5 0])
% set(Eplt(2), 'color', [0 0 1])
set(EoutLaser_plt(1), 'linewidth', 1.5, 'linestyle', ':', 'color', [1 0.5 0]);
set(EoutLaser_plt(2), 'linewidth', 0.5, 'linestyle', '--', 'color', [0 0 0.5]);

% 
% %
% N = size(FreqShift,2)/2;
% M_space = zeros(1,N*2);
% 
% for n=0:(size(FreqShift,2)-1)/2
%     M_space(n+1) = sum(FreqShift(:,N+n) + FreqShift(:,N-n));
% end

% 
% % Energy State
% En_map = M_space./max(M_space);
% Ph_map = M_space;
% 
% 
% % Energy Output
% for n=1:15
%     E_STATES(n) = (n^2*(hPlanck*2*pi)^2)/(8*me*SPACE.length^2);
% end
% n_ = n;
% 
% E_output = En_map(1:n_).*E_STATES/abs(e0);
% 
% E_photons = Ph_map(1:n_);
% % ./E_STATES
% 
% 
% % Plot Energy Output
% En_fig = figure(2);
% subplot(2,1,1);
% Eplt = plot(E_output);
haxes = get(EoutLaser_fig,'CurrentAxes');
set(haxes,'FontSize', 10, 'FontAngle', 'italic', 'YDir', 'normal',...
    'XGrid', 'on', 'YGrid', 'on');
ylabel('\DeltaE_o_u_t (eV), ',...
    'VerticalAlignment','bottom','HorizontalAlignment','center',...
    'FontWeight','normal', 'FontSize',12, 'FontAngle', 'normal');
xlim([taxis(1) taxis(end)]);
% title(['Population of States: [E_1, E_2]'],...
%     'FontWeight','bold', 'FontSize',14, 'FontAngle', 'normal');
xlabel('Time (s)',...
    'VerticalAlignment','cap','HorizontalAlignment','center',...
    'FontWeight','normal', 'FontSize',12, 'FontAngle', 'normal');






% PLOT Laser
% fig2 = figure(2);
subplot(3,3,1:3);

plot_CEP=plot(taxis,E_laser, 'linewidth', 2);
% plot_CEP=area(E_laser,taxis)

set(gca,'FontSize', 10, 'FontAngle', 'italic', 'XGrid', 'on', 'YGrid', 'on');

ylim([-abs(max(E_laser))-1 abs(max(E_laser))+1]);
xlim([taxis(1) taxis(end)]);

% set(plot_CEP, 'FaceColor', [0 0.7 0]);
set(plot_CEP, 'Color', [0 0.7 0]);

title([TEXT_.txt_EnPot ', ' TEXT_.txt_Laser],... % ',' pulse_txt
    'FontWeight','bold', 'FontSize',14, 'FontAngle', 'normal');%,'Color', [0 0 1]);
xlabel('Time (s)',...
    'VerticalAlignment','cap','HorizontalAlignment','center',...
    'FontWeight','normal', 'FontSize',12, 'FontAngle', 'normal');
xlim([taxis(1) taxis(end)]);
ylabel('E_e_x_t (V/cm)',...
    'VerticalAlignment','bottom','HorizontalAlignment','center',...
    'FontWeight','normal', 'FontSize',12, 'FontAngle', 'normal');

%%
set(EoutLaser_fig, 'Position', [100 100 900 700]);
% , [589   217   743   656]);
saveas(EoutLaser_fig, [TEXT_.saveVisual '\{E_Real}.' TEXT_.SaveGraphicName]);


% set(nPfig, 'Position', [589   217   743   656]);
% saveas(nPfig, ['PHOTONS_' SaveGraphicName]);


clear En_real E_laser;



end
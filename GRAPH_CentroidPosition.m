function GRAPH_CentroidPosition(TEXT_, SIM, taxis, xaxis)
%%
close all hidden;
CONSTANTS;
load([TEXT_.SaveDirectory '/Centroid_X'], 'Centroid_X');
load([TEXT_.SaveDirectory '/E_laser'], 'E_laser');

%%

Cenfig = figure;
% imagesc(En_);
% figure;
subplot(3,3,4:9);
Cenplt = plot(taxis, Centroid_X, taxis, zeros(size(Centroid_X)));
% title('E(1)');
% legend('E_1', 'E_2');

% set(Cenplt(1), 'color', [0.7 0 0.7])
% set(Cenplt(2), 'color', [0 0 1])

set(Cenplt(1), 'linewidth', 1.5, 'linestyle', ':', 'color', [0.7 0 0.7]);
set(Cenplt(2), 'linewidth', 0.5, 'linestyle', '--', 'color', [0.5 0 0.5]);

haxes = get(Cenfig,'CurrentAxes');
set(haxes,'FontSize', 10, 'FontAngle', 'italic', 'YDir', 'normal',...
    'XGrid', 'on', 'YGrid', 'on');
% title(['Centroid (position) VS Time'],...
%     'FontWeight','bold', 'FontSize',14, 'FontAngle', 'normal');
ylabel('Centroid, <X> (m_.)',...
    'VerticalAlignment','bottom','HorizontalAlignment','center',...
    'FontWeight','normal', 'FontSize',12, 'FontAngle', 'normal');
xlim([taxis(1) taxis(end)]);
% set(gca,'XTickLabel',{' '});
% xlim([0 100]);
xlabel('Time (s)',...
    'VerticalAlignment','cap','HorizontalAlignment','center',...
    'FontWeight','normal', 'FontSize',12, 'FontAngle', 'normal');
% 

axis([taxis(1) taxis(end) -max(abs(Centroid_X)) max(abs(Centroid_X))]);

% 
% % Plot Photon Output
% % nPfig = figure(2);
% subplot(2,1,2);
% nPplt = plot(E_photons);
% set(nPplt, 'color', [1 0 0]);
% haxes = get(En_fig,'CurrentAxes');
% set(haxes,'FontSize', 10, 'FontAngle', 'italic', 'YDir', 'normal',...
%     'XGrid', 'on', 'YGrid', 'on');
% title([num2str(TEXT_.txt_EnPot, '%2.1e') ', ' TEXT_.txt_Laser ],...
%     'FontWeight','bold', 'FontSize',14, 'FontAngle', 'normal');
% xlabel('Harmonic Order (n)',...
%     'VerticalAlignment','cap','HorizontalAlignment','center',...
%     'FontWeight','normal', 'FontSize',12, 'FontAngle', 'normal');
% xlim([1 n_]);
% ylabel('Total # of Photons',...
%     'VerticalAlignment','bottom','HorizontalAlignment','center',...
%     'FontWeight','normal', 'FontSize',12, 'FontAngle', 'normal');


%%
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

% set(hfig, 'Position', [[993 136 560 656]]);

% saveas(fig2, ['SRC_' SaveGraphicName]);
% end



%%
% set(M_src_fig, 'Position', [993 136 560 656]);
% set(Cenfig, 'Position', [589   217   743   656]);
set(Cenfig, 'Position', [100 100 900 700]);
%, [98 	62  	1089	592]);

saveas(Cenfig, [TEXT_.saveVisual '\{Centroid_Pos}.' TEXT_.SaveGraphicName]);

clear Centroid_X E_laser;
% set(Cenfig, 'Position', [589   217   743   656]);
% saveas(Cenfig, ['ENERGY_' SaveGraphicName]);


% set(nPfig, 'Position', [589   217   743   656]);
% saveas(nPfig, ['PHOTONS_' SaveGraphicName]);





end
function GRAPH_SpectralAnalysis(TEXT_, SIM, PULSE, taxis, xaxis, faxis)
%% Spectral Analysis
close all hidden;
CONSTANTS;
load([TEXT_.SaveDirectory '/FreqShift'], 'FreqShift');
load([TEXT_.SaveDirectory '/E_laser'], 'E_laser');

%%
% taxis = TIME.delta*(1:size(FreqShift,1))*TIME.save;
% xaxis = SPACE.delta*(1:size(FreqShift,2));

dF = abs(faxis(end/2+1)-faxis(end/2));
ndF = 36; % 20


% FreqNorm = FreqShift./max(max(max(max(FreqShift))));
% E_out = FreqShift.^2;%/(2*me);

% FreqNorm = FreqShift(:,1:end).*faxis;
% faxis_label = faxis/faxis(end/2)
% faxis_label = faxis/abs(faxis(end)-faxis(1));
% Es = faxis/faxis(2);

% harmonic_level = -(length(faxis_label)/2-1):length(faxis_label)/2;
% harmonic_level = -(length(faxis)/2-1):length(faxis)/2;

%% Plot
M_src_fig = figure(1);
% PLOT Momentum Space Map
% subplot(1,5,3:5);
subplot(3,3,4:9);
M_plt = imagesc(taxis, faxis, rot90(FreqShift));
% M_src_fig= gcf;
% hfig = imagesc(Es, taxis, FreqShift);
% XTickLabel(num2str(Es));
% haxes = get(gcf,'CurrentAxes');

M_axes = get(M_plt, 'parent');

% set(gca,'FontSize', 10, 'FontAngle', 'italic', 'YGrid', 'on');
set(M_axes,'FontSize', 10, 'FontAngle', 'italic', 'YDir', 'normal',...
    'XGrid', 'on', 'YGrid', 'on');



% title(['|\psi(t,p)|^2 (PDF)'],... 'HorizontalAlignment', 'right',...
%     'VerticalAlignment','bottom','HorizontalAlignment','center',...
%     'FontWeight','normal', 'FontSize',12, 'FontAngle', 'normal');


ylabel(M_axes,['Freq. (k_.)'],...
    'VerticalAlignment','bottom','HorizontalAlignment','center',...
    'FontWeight','normal', 'FontSize',12, 'FontAngle', 'normal');

xlabel(M_axes,'|\psi(t,k)|^2',...
    'VerticalAlignment','middle','HorizontalAlignment','center',...
    'FontWeight','normal', 'FontSize',12, 'FontAngle', 'normal');
set(M_axes,'XTickLabel',{' '});
%\Deltat


% xlim([-1*10^10 1*10^10]);

% ylim(M_axes,(dF)*[-ndF/4 ndF/4]);
% ylim(M_axes,(dF)*[-ndF/3 ndF/3]);
% ylim(M_axes,(dF)*[-ndF/2 ndF/2]);
ylim(M_axes, dF*[-ndF ndF]);

% rng = abs(-ndF:ndF)+1;

rng = -(ndF-mod(ndF,((length(get(M_axes,'XTickLabel'))-1)/2))):...
    floor(ndF/((length(get(M_axes,'XTickLabel'))-1)/2)):...
    (ndF-mod(ndF,((length(get(M_axes,'XTickLabel'))-1)/2)));
harmonic_rng = abs(rng)+1;

% rng = abs(-((length(get(haxes,'XTickLabel'))-1)/2):...
%     ((length(get(haxes,'XTickLabel'))-1)/2))+1;


% set(haxes, 'XTickLabel', harmonic_rng);
% cbar = colorbar;
% caxis(haxes, 'CLim', [0 1]);
caxis(M_axes,[0 0.5]);
cmap = colormap('jet');
% % cmap(1,:) = [1 1 1]; cmap(1,:) = [0.95 0.95 0.95];
colormap(M_axes, cmap);% clabels = {0:0.25:1};
% lcbar = lcolorbar(clabels);
%,  'CLim', [0 1]);

% set(cbar, 'location', 'northoutside');
% set(haxes, 'XTick', [df*(-5:5)]);

% get(haxes,'XTickLabel')
% (length(get(haxes,'XTickLabel'))-1)/2

% set(haxes, 'XMinorGrid', 'on')
% set(haxes, 'XMinorTick', 'on')

% set(haxes, 'XTickLabel', [round(Es(end/2 - 3)):round(Es(end/2 + 3))]);
% set(haxes, 'YTickLabel', [T_len/10:T_len/10:T_len]);


% xrng = linspace(min(xlim),max(xlim),59);
% xrng = min(xlim):abs(faxis(2)-faxis(1)):max(xlim)
% abs(faxis(2)-faxis(1))

% gap = xlim/abs(faxis(4)-faxis(1));
% harmonicLvl.min = min(ceil((gap)));
% harmonicLvl.max = max(floor((gap)));
% harmonicLvl.rng = harmonicLvl.min:2:harmonicLvl.max;


% set(haxes, 'XTick', xrng);
% set(haxes, 'XTickLabel', harmonicLvl.rng);


% 
% % Determine Name
% set(M_src_fig, 'Position', [993 136 560 656]);


%%
% PLOT Laser
% fig2 = figure(2);
subplot(3,3,1:3);
plot_CEP=plot(taxis,(rot90(E_laser)), 'linewidth', 2);
% plot_CEP=area(E_laser,taxis)
CEP_axes = get(plot_CEP, 'parent');
set(CEP_axes,'FontSize', 10, 'FontAngle', 'italic', 'XGrid', 'on', 'YGrid', 'on');

ylim(CEP_axes,[-abs(max(E_laser))-1 abs(max(E_laser))+1]);
xlim(CEP_axes,[taxis(1) taxis(end)]);

% set(plot_CEP, 'FaceColor', [0 0.7 0]);
set(plot_CEP, 'Color', [0 0.7 0]);

title(CEP_axes,[TEXT_.txt_EnPot ', ' TEXT_.txt_Laser],... % ',' pulse_txt
    'FontWeight','bold', 'FontSize',14, 'FontAngle', 'normal');%,'Color', [0 0 1]);
xlabel(CEP_axes,'Time (s)',...
    'VerticalAlignment','cap','HorizontalAlignment','center',...
    'FontWeight','normal', 'FontSize',12, 'FontAngle', 'normal');
ylabel(CEP_axes,'E_e_x_t (V/cm)',...
    'VerticalAlignment','bottom','HorizontalAlignment','center',...
    'FontWeight','normal', 'FontSize',12, 'FontAngle', 'normal');

% set(hfig, 'Position', [[993 136 560 656]]);

% saveas(fig2, ['SRC_' SaveGraphicName]);
% end

%%
% set(M_src_fig, 'Position', [993 136 560 656]);
% set(M_src_fig, 'Position', [603   173   1026   656]);
set(M_src_fig, 'Position', [100 100 900 700]);
% 'Position', [100 100 900 700]);

saveas(M_src_fig, [TEXT_.saveVisual '\{FFT+SRC}.' TEXT_.SaveGraphicName]);

clear FreqShift;

end
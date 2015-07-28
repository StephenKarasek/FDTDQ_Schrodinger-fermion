function GRAPH_SpectralAnalysis_PopulationLevel(TEXT_, SIM, taxis, xaxis, faxis, varargin)
%% Spectral Analysis
close all hidden;
CONSTANTS;
load([TEXT_.SaveDirectory '/E_laser'], 'E_laser');
load([TEXT_.SaveDirectory '/FreqShift'], 'FreqShift');
load([TEXT_.SaveDirectory '/PopulationState'], 'PopulationState');

%% input
if nargin > 5
    
end


%%
% taxis = TIME.delta*(1:size(FreqShift,1))*TIME.save;
% xaxis = SPACE.delta*(1:size(FreqShift,2));

dF = abs(faxis(end/2+1)-faxis(end/2));
ndF = 50;

% FreqNorm = FreqShift./max(max(max(max(FreqShift))));
% E_out = FreqShift.^2;%/(2*me);

% FreqNorm = FreqShift(:,1:end).*faxis;
% faxis_label = faxis/faxis(end/2)
% faxis_label = faxis/abs(faxis(end)-faxis(1));
% Es = faxis/faxis(2);

% harmonic_level = -(length(faxis_label)/2-1):length(faxis_label)/2;
% harmonic_level = -(length(faxis)/2-1):length(faxis)/2;

%% Plot
SpectralPopStatefig = figure(1);
% PLOT Momentum Space Map
% subplot(8,1,3:5);
% subplot(6,1,2:4);
% subplot(5,1,2:3);
subplot(10,1,4:7);


M_plt = imagesc(taxis, faxis, rot90(FreqShift));
% M_src_fig= gcf;
% hfig = imagesc(Es, taxis, FreqShift);
% XTickLabel(num2str(Es));
haxes = get(gcf,'CurrentAxes');
% set(gca,'FontSize', 10, 'FontAngle', 'italic', 'YGrid', 'on');
set(haxes,'FontSize', 10, 'FontAngle', 'italic', 'YDir', 'normal',...
    'XGrid', 'on', 'YGrid', 'on');
% set(haxes, 'XTickLabel


% title(['|\psi(t,p)|^2 (PDF)'],... 'HorizontalAlignment', 'right',...
%     'VerticalAlignment','bottom','HorizontalAlignment','center',...
%     'FontWeight','normal', 'FontSize',12, 'FontAngle', 'normal');


ylabel(['Freq. Response (k_.)'],...
    'VerticalAlignment','bottom','HorizontalAlignment','center',...
    'FontWeight','normal', 'FontSize',12, 'FontAngle', 'normal');

% xlabel('|\psi(t,k)|^2',...
%     'VerticalAlignment','middle','HorizontalAlignment','center',...
%     'FontWeight','normal', 'FontSize',12, 'FontAngle', 'normal');
set(gca,'XTickLabel',{' '});
%\Deltat


% xlim([-1*10^10 1*10^10]);

ylim((dF)*[-ndF/2 ndF/2]);
% rng = abs(-ndF:ndF)+1;

rng = -(ndF-mod(ndF,((length(get(haxes,'XTickLabel'))-1)/2))):...
    floor(ndF/((length(get(haxes,'XTickLabel'))-1)/2)):...
    (ndF-mod(ndF,((length(get(haxes,'XTickLabel'))-1)/2)));
harmonic_rng = abs(rng)+1;

% rng = abs(-((length(get(haxes,'XTickLabel'))-1)/2):...
%     ((length(get(haxes,'XTickLabel'))-1)/2))+1;


% set(haxes, 'XTickLabel', harmonic_rng);
% cbar = colorbar;
% caxis(haxes, 'CLim', [0 1]);
caxis([0 0.5]);
cmap = colormap('jet');
% % cmap(1,:) = [1 1 1];
cmap(1,:) = [.95 .95 .95];
colormap(cmap);

cbar = colorbar('location', 'southoutside', 'FontSize', 10, 'FontAngle', 'italic');

% cbar = colorbar('location', 'north', 'FontSize', 10, 'FontAngle', 'italic');
% clabels = {' '};

% clabels = {0:0.25:1};
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
subplot(10,1,8:10);
% subplot(8,1,6:8);
% subplot(6,1,5:6);
% subplot(5,1,4:5);
% imagesc(En_);
% figure;
% Eplt = plot(taxis, E_well(:,1), taxis, E_well(:,2), taxis, E_psi);
% title('E(1)');

EnSt0 = [0.0   0.7   0.0]; % Green
EnSt1 = [0.7   0.0   0.0]; % Red
EnSt2 = [0.0   0.0   0.7]; % Blue
EnSt3 = [1.0   0.7   0.0]; % Orange
EnSt4 = [0.0   0.0   0.7]; % Yellow

switch size(PopulationState,2)
    case 1
        Eplt = plot(taxis, PopulationState(:,1));
        hleg = legend('E_0');%, 'E_n->inf');
        
%         set(Eplt(1), 'color', [0 0.7 0]);
%         set(Eplt(2), 'color', [0 0 0.7]);
            set(Eplt(1), 'linewidth',1.5, 'linestyle', ':','color', EnSt0);
%                 'marker', 'o', 'markerfacecolor', [0 0.7 0]) ;       
%             set(Eplt(1), 'linewidth',1.5, 'marker', 'o', 'markerfacecolor', [0 0 0.7]) ;       
        
    case 2
        Eplt = plot(taxis, PopulationState(:,1),...
            taxis, PopulationState(:,2));
        hleg = legend('E_0', 'E_1');
        
%         set(Eplt(1), 'color', [0 0.7 0]);
%         set(Eplt(2), 'color', [0.7 0.7 0]);
%         set(Eplt(3), 'color', [0 0 0.7]);
        set(Eplt(1), 'linewidth',1.5, 'linestyle', ':','color', EnSt0);
%             'marker', 'o', 'markerfacecolor', [0 0.7 0]) ;       
        set(Eplt(2), 'linewidth',1.5, 'linestyle', ':','color', EnSt1);
%             'marker', 'o', 'markerfacecolor', [0.7 0.7 0]) ;       
%         set(Eplt(1), 'linewidth',1.5, 'marker', 'o', 'markerfacecolor', [0 0.7 0]) ;       

    case 3
        Eplt = plot(taxis, PopulationState(:,1),...
            taxis, PopulationState(:,2),...
            taxis, PopulationState(:,3));
        hleg = legend('E_0', 'E_1', 'E_2');
        
        set(Eplt(1), 'linewidth',1.5, 'linestyle', ':','color', EnSt0);
% set(Eplt(1), 'linewidth',1.5, 'linestyle', '--','color', [0 0.7 0],'marker', 'o', 'markerfacecolor', [0 0.7 0]);
        set(Eplt(2), 'linewidth',1.5, 'linestyle', ':','color', EnSt1);
%             'marker', 'o', 'markerfacecolor', [0.7 0.7 0]);
        set(Eplt(3), 'linewidth',1.5, 'linestyle', ':','color', EnSt2);
%             'marker', 'o', 'markerfacecolor', [0.7 0 0]);
%         set(Eplt(4), 'color', [0 0 0.7]);
        

    case 4
        Eplt = plot(taxis, PopulationState(:,1),...
            taxis, PopulationState(:,2),...
            taxis, PopulationState(:,3),...
            taxis, PopulationState(:,4));
        hleg = legend('E_0', 'E_1', 'E_2', 'E_3');
        
        set(Eplt(1), 'linewidth',1.5, 'linestyle', ':','color', EnSt1);
%             'marker', 'o', 'markerfacecolor', [0 0.7 0]);
        set(Eplt(2), 'linewidth',1.5, 'linestyle', ':','color', EnSt1);
%             'marker', 'o', 'markerfacecolor', [0.7 0.7 0]);
        set(Eplt(3), 'linewidth',1.5, 'linestyle', ':','color', EnSt2);
%             'marker', 'o', 'markerfacecolor', [0.7 0 0]);
        set(Eplt(4), 'linewidth', 1.5, 'linestyle', ':','color', EnSt3);
%             'marker', 'o', 'markerfacecolor', [1 0.7 0]);
%         set(Eplt(5), 'color', [0 0 0.7]);
        

        

    otherwise
end

haxes = get(SpectralPopStatefig,'CurrentAxes');
set(haxes,'FontSize', 10, 'FontAngle', 'italic', 'YDir', 'normal',...
    'XGrid', 'on', 'YGrid', 'on');
% title(['Population of States: [E_0, E_1, E_2, E_3]'],...
%     'FontWeight','bold', 'FontSize',14, 'FontAngle', 'normal');

% set(hleg,'Orientation', 'horizontal', 'Location', 'NorthEast');
% xlabel('Time (s)',...

% xlabel('Energy State Population (normalized)',...
% xlabel('',...
%     'VerticalAlignment','middle','HorizontalAlignment','center',...
%     'FontWeight','normal', 'FontSize',12, 'FontAngle', 'normal');

% set(hleg,'Orientation', 'horizontal', 'Location', 'SouthOutside', 'FontSize', 10, 'FontAngle', 'italic');

set(hleg,'Orientation', 'horizontal', 'Location', 'NorthEast', 'FontSize', 10, 'FontAngle', 'italic');

xlim([taxis(1) taxis(end)]);
set(gca,'XTickLabel',{' '});

% ylim([0 100]);
ylabel('Population Level',...
    'VerticalAlignment','bottom','HorizontalAlignment','center',...
    'FontWeight','normal', 'FontSize',12, 'FontAngle', 'normal');
ylim([0 1.1]);


%%
% PLOT Laser
% fig2 = figure(2);
subplot(10,1,1:2);
% subplot(6,1,1);
% subplot(8,1,1:2);
plot_CEP=plot(taxis,(rot90(E_laser)), 'linewidth', 2);
% plot_CEP=area(E_laser,taxis)

set(gca,'FontSize', 10, 'FontAngle', 'italic', 'XGrid', 'on', 'YGrid', 'on');

ylim([-abs(max(E_laser))-1 abs(max(E_laser))+1]);
xlim([taxis(1) taxis(end)]);

% set(plot_CEP, 'FaceColor', [0 0.7 0]);
set(plot_CEP, 'Color', [0 0.7 0]);

% title([TEXT_.txt_EnPot ', ' TEXT_.txt_Laser],... % ',' pulse_txt
%     'FontWeight','bold', 'FontSize',14, 'FontAngle', 'normal');%,'Color', [0 0 1]);
xlabel('Time (s)',...
    'VerticalAlignment','cap','HorizontalAlignment','center',...
    'FontWeight','normal', 'FontSize',12, 'FontAngle', 'normal');

% ylabel('I_0 (V/cm^2 )',...
ylabel('E_e_x_t (V/m)',...
    'VerticalAlignment','bottom','HorizontalAlignment','center',...
    'FontWeight','normal', 'FontSize',12, 'FontAngle', 'normal');

% set(hfig, 'Position', [[993 136 560 656]]);

% saveas(fig2, ['SRC_' SaveGraphicName]);
% end

%%
% set(M_src_fig, 'Position', [993 136 560 656]);
% set(M_src_fig, 'Position', [603   173   1026   656]);
set(SpectralPopStatefig, 'Position', [100 100 900 900]);
% 'Position', [100 100 900 700]);
%%
saveas(SpectralPopStatefig, ['FFT+POP+SRC_' TEXT_.SaveGraphicName]);

clear FreqShift;

end
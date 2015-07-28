function GRAPH_PopulationStateNorm_GIVEN(TEXT_, SIM, PULSE, taxis, xaxis)
%% 
close all hidden;
CONSTANTS;
% load([TEXT_.SaveDirectory '/PopulationStateNorm'], 'PopulationStateNorm');
load([TEXT_.SaveDirectory '/PopulationStateNorm_GIVEN']);
load([TEXT_.SaveDirectory '/E_laser'], 'E_laser');

%%
PopStatefig = figure;

subplot(3,3,4:9);
% imagesc(En_);
% figure;
% Eplt = plot(taxis, E_well(:,1), taxis, E_well(:,2), taxis, E_psi);
% title('E(1)');

% EnSt0 = [0.0   0.7   0.0]; % Green
% EnSt1 = [0.0   0.0   0.7]; % Blue
% EnSt2 = [0.7   0.0   0.0]; % Red
% EnSt3 = [0.7   0.7   0.3]; % Yellow


EnSt3 = [0.0   0.7   0.0]; % Green
EnSt1 = [0.7   0.0   0.0]; % Red
EnSt2 = [0.7   0.7   0.3]; % Yellow
EnSt0 = [0.0   0.0   0.7]; % Blue
EnSt_LineWidth = 2;

Estates_ = PULSE.NumViewStates;
% Estates_ = size(PopulationStateNorm,2);
%     Estates_ = 2;
switch Estates_
    case 1
        Eplt = plot(taxis, PopulationStateNorm_GIVEN(:,1));
        hleg = legend('E_0');%, 'E_n->inf');
        
%         set(Eplt(1), 'color', [0 0.7 0]);
%         set(Eplt(2), 'color', [0 0 0.7]);
            set(Eplt(1), 'linewidth',EnSt_LineWidth, 'linestyle', ':','color', EnSt0);
%                 'marker', 'o', 'markerfacecolor', [0 0.7 0]) ;       
%             set(Eplt(1), 'linewidth',EnSt_LineWidth, 'marker', 'o', 'markerfacecolor', [0 0 0.7]) ;       
        
    case 2
        Eplt = plot(taxis, PopulationStateNorm_GIVEN(:,1),...
            taxis, PopulationStateNorm_GIVEN(:,2));
        hleg = legend('E_0', 'E_1');
        
%         set(Eplt(1), 'color', [0 0.7 0]);
%         set(Eplt(2), 'color', [0.7 0.7 0]);
%         set(Eplt(3), 'color', [0 0 0.7]);
        set(Eplt(1), 'linewidth',EnSt_LineWidth, 'linestyle', ':','color', EnSt0);
%             'marker', 'o', 'markerfacecolor', [0 0.7 0]) ;       
        set(Eplt(2), 'linewidth',EnSt_LineWidth, 'linestyle', ':','color', EnSt1);
%             'marker', 'o', 'markerfacecolor', [0.7 0.7 0]) ;       
%         set(Eplt(1), 'linewidth',EnSt_LineWidth, 'marker', 'o', 'markerfacecolor', [0 0.7 0]) ;       

    case 3
        Eplt = plot(taxis, PopulationStateNorm_GIVEN(:,1),...
            taxis, PopulationStateNorm_GIVEN(:,2),...
            taxis, PopulationStateNorm_GIVEN(:,3));
        hleg = legend('E_0', 'E_1', 'E_2');
        
        set(Eplt(1), 'linewidth',EnSt_LineWidth, 'linestyle', ':','color', EnSt0);
% set(Eplt(1), 'linewidth',EnSt_LineWidth, 'linestyle', '--','color', [0 0.7 0],'marker', 'o', 'markerfacecolor', [0 0.7 0]);
        set(Eplt(2), 'linewidth',EnSt_LineWidth, 'linestyle', ':','color', EnSt1);
%             'marker', 'o', 'markerfacecolor', [0.7 0.7 0]);
        set(Eplt(3), 'linewidth',EnSt_LineWidth, 'linestyle', ':','color', EnSt2);
%             'marker', 'o', 'markerfacecolor', [0.7 0 0]);
%         set(Eplt(4), 'color', [0 0 0.7]);
        

    case 4
        Eplt = plot(taxis, PopulationStateNorm_GIVEN(:,1),...
            taxis, PopulationStateNorm_GIVEN(:,2),...
            taxis, PopulationStateNorm_GIVEN(:,3),...
            taxis, PopulationStateNorm_GIVEN(:,4));
        hleg = legend('E_0', 'E_1', 'E_2', 'E_3');
        
        set(Eplt(1), 'linewidth',EnSt_LineWidth, 'linestyle', ':','color', EnSt0);
%             'marker', 'o', 'markerfacecolor', [0 0.7 0]);
        set(Eplt(2), 'linewidth',EnSt_LineWidth, 'linestyle', ':','color', EnSt1);
%             'marker', 'o', 'markerfacecolor', [0.7 0.7 0]);
        set(Eplt(3), 'linewidth',EnSt_LineWidth, 'linestyle', ':','color', EnSt2);
%             'marker', 'o', 'markerfacecolor', [0.7 0 0]);
        set(Eplt(4), 'linewidth', 1.5, 'linestyle', ':','color', EnSt3);
%             'marker', 'o', 'markerfacecolor', [1 0.7 0]);
%         set(Eplt(5), 'color', [0 0 0.7]);
        

    case 5
        Eplt = plot(taxis, PopulationStateNorm_GIVEN(:,1),...
            taxis, PopulationStateNorm_GIVEN(:,2),...
            taxis, PopulationStateNorm_GIVEN(:,3),...
            taxis, PopulationStateNorm_GIVEN(:,4));%,...
%             taxis, PopulationStateNorm_GIVEN(:,5));
        hleg = legend('E_0', 'E_1', 'E_2', 'E_3');%, 'E_4');
        
        set(Eplt(1), 'linewidth',EnSt_LineWidth, 'linestyle', ':','color', EnSt0);
%             'marker', 'o', 'markerfacecolor', [0 0.7 0]);
        set(Eplt(2), 'linewidth',EnSt_LineWidth, 'linestyle', ':','color', EnSt1);
%             'marker', 'o', 'markerfacecolor', [0.7 0.7 0]);
        set(Eplt(3), 'linewidth',EnSt_LineWidth, 'linestyle', ':','color', EnSt2);
%             'marker', 'o', 'markerfacecolor', [0.7 0 0]);
        set(Eplt(4), 'linewidth', 1.5, 'linestyle', ':','color', EnSt3);
%             'marker', 'o', 'markerfacecolor', [1 0.7 0]);
%         set(Eplt(5), 'linewidth', 1.5, 'linestyle', ':','color', EnSt4);
        

    otherwise
end

E_axes_ = get(Eplt,'parent');
if iscell(E_axes_)
    E_axes = E_axes_{1};
else
    E_axes = E_axes_;
end

% haxes = get(PopStatefig,'CurrentAxes');
set(E_axes,'FontSize', 10, 'FontAngle', 'italic', 'YDir', 'normal',...
    'XGrid', 'on', 'YGrid', 'on');
% title(['Population of States: [E_0, E_1, E_2, E_3]'],...
%     'FontWeight','bold', 'FontSize',14, 'FontAngle', 'normal');

% set(hleg,'Orientation', 'horizontal', 'Location', 'NorthEast');
% xlabel('Time (s)',...

% xlabel(E_axes,'Energy State Population',...
%     'VerticalAlignment','middle','HorizontalAlignment','center',...
%     'FontWeight','normal', 'FontSize',12, 'FontAngle', 'normal');
set(hleg,'Orientation', 'horizontal', 'Location', 'SouthOutside', ...
    'FontSize', 10, 'FontAngle', 'italic');

xlim(E_axes,[taxis(1) taxis(end)]);
set(E_axes,'XTickLabel',{' '});

% ylim([0 100]);
ylabel(E_axes,'Population Level',...
    'VerticalAlignment','bottom','HorizontalAlignment','center',...
    'FontWeight','normal', 'FontSize',12, 'FontAngle', 'normal');
ylim(E_axes,[0 1.1]);
% ylim(E_axes,[-1.1 1.1]);

% set(hleg,'Orientation', 'horizontal', 'Location', 'North');
% ,'VerticalAlignment', 'bottom');


%%
% PLOT Laser
% fig2 = figure(2);
subplot(3,3,1:3);

plot_CEP=plot(taxis,E_laser, 'linewidth', 2);
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
xlim(CEP_axes,[taxis(1) taxis(end)]);
ylabel(CEP_axes,'E_e_x_t (V/cm)',...
    'VerticalAlignment','bottom','HorizontalAlignment','center',...
    'FontWeight','normal', 'FontSize',12, 'FontAngle', 'normal');


%%
% set(PopStatefig, 'Position', [589   217   743   656]);
set(PopStatefig, 'Position', [100 100 900 700]);
saveas(PopStatefig, [TEXT_.saveVisual '\{PopLvl_Norm}.' TEXT_.SaveGraphicName]);


% set(nPfig, 'Position', [589   217   743   656]);
% saveas(nPfig, ['PHOTONS_' SaveGraphicName]);

clear PopulationStateNorm E_laser;

end
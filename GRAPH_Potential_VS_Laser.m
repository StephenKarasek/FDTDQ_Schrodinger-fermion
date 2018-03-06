function [] = GRAPH_Potential_VS_Laser(TEXT_, SIM, taxis, xaxis)
%% Potential vs Laser
close all hidden;
CONSTANTS;
load([TEXT_.SaveDirectory '/V_pot'], 'V_pot');
load([TEXT_.SaveDirectory '/E_laser'], 'E_laser');


%% Axes Range
% V_max = max(max(max(max(V_pot/abs(e0)))));
% V_min = min(min(min(min(V_pot/abs(e0)))));
% rng = V_max - V_min;

V_lims = [-100 100];



%% Plot
Potential_tilt_fig = figure;

subplot(3,3,4:9);
P_plt = imagesc(taxis, xaxis, rot90(V_pot/abs(e0)));
P_axes = get(P_plt, 'parent');
set(P_axes,'FontSize', 10, 'FontAngle', 'italic', 'YDir', 'normal',...
    'XGrid', 'on', 'YGrid', 'on');

switch SIM.PotentialMap
    case {'H_1'; 'H1'}
caxis(P_axes, [V_lims(1) V_lims(2)]); cmap = colormap('jet');
% cmap(1,:) = [1 1 1]; cmap(end,:) = [1 1 1];

    case {'H_2'; 'H2';'H_2+'; 'H2+'}
caxis(P_axes, [V_lims(1) V_lims(2)]); cmap = colormap('jet');
% cmap(1,:) = [1 1 1]; cmap(end,:) = [1 1 1];

    case {'InfiniteSquareWell';'SquareWell';'InfSqWell'; 'Square';'Square_'}
caxis(P_axes, [V_lims(1) V_lims(2)]); cmap = colormap('jet');
% cmap(end/2-1,:) = [1 1 1]; cmap(end/2+1,:) = [1 1 1];

    otherwise
end

hcmap = colormap(P_axes, cmap);
cbar = colorbar('location', 'southoutside', 'FontSize', 10, 'FontAngle', 'italic');


% CLIMS(1) = -100;
% CLIMS(2) = 100;

% CLIMS = [-0.2*abs(min(min(min(V_atom/abs(e0)))))...
%   0.1*abs(max(max(max(V_pot/abs(e0)))))];

% if CLIMS(1)<-100; CLIMS(1) = -100; end
% if CLIMS(2)>100; CLIMS(2) = 100; end

% set(gca,  'CLim', CLIMS);

% txt_pot = ['V(' SIM.PotentialMap ')'];
% cbar = colorbar;
% set(cbar, 'location', 'southoutside','FontSize', 10, 'FontAngle', 'italic');
xlim(P_axes, [taxis(1) taxis(end)]);

% title(['Effective Potential (eV)'],...
%     'VerticalAlignment','bottom',...%'HorizontalAlignment', 'right',...
%     'FontWeight','bold', 'FontSize',14, 'FontAngle', 'normal');

xlabel(P_axes, 'V_e_f_f (eV)',...
    'VerticalAlignment','cap','HorizontalAlignment','center',...
    'FontWeight','normal', 'FontSize',12, 'FontAngle', 'normal');
set(P_axes,'XTickLabel',{' '});
ylabel(P_axes, 'Position (m_.)',...
    'VerticalAlignment','bottom','HorizontalAlignment','center',...
    'FontWeight','normal', 'FontSize',12, 'FontAngle', 'normal');

%%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Laser Plot
% PLOT Laser
% fig2 = figure(2);
subplot(3,3,1:3);
plot_CEP=plot(taxis,E_laser, 'linewidth', 2);
CEP_axes = get(plot_CEP, 'parent');

set(CEP_axes,'FontSize', 10, 'FontAngle', 'italic', 'XGrid', 'on', 'YGrid', 'on');
ylim(CEP_axes,[-abs(max(E_laser))-1 abs(max(E_laser))+1]);
xlim(CEP_axes,[taxis(1) taxis(end)]);

set(plot_CEP, 'Color', [0 0.7 0]);

title(CEP_axes,[TEXT_.txt_EnPot ', ' TEXT_.txt_Laser],... % ',' pulse_txt
    'FontWeight','bold', 'FontSize',14, 'FontAngle', 'normal');%,'Color', [0 0 1]);
xlabel(CEP_axes,'Time (s)',...
    'VerticalAlignment','cap','HorizontalAlignment','center',...
    'FontWeight','normal', 'FontSize',12, 'FontAngle', 'normal');
ylabel(CEP_axes,'E_e_x_t (V/cm)',...
    'VerticalAlignment','bottom','HorizontalAlignment','center',...
    'FontWeight','normal', 'FontSize',12, 'FontAngle', 'normal');




%%
% Determine Name
% set(Potential_tilt_fig, 'Position', [603   173   1026   656]);
set(Potential_tilt_fig, 'Position', [100 100 900 700]);
saveas(Potential_tilt_fig, [TEXT_.saveVisual '\{V_tilt}.' TEXT_.SaveGraphicName]);

clear V_pot E_laser;

end
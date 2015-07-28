function [] = GRAPH_(SPACE, TIME, SIM, PULSE, ITP, LASER, DATA, DEBUG)

%% GRAPHS
CONSTANTS;

%% Space-Time Axis & Wave Vector Labels
taxis = (TIME.save:TIME.save:TIME.N).*TIME.delta;
xaxis = SPACE.Axis;

xNlength = SPACE.length*SPACE.ZeroPad_Factor;
xNaxis = linspace(-xNlength/2,xNlength/2,SPACE.N0);

faxis = fftaxisshift(fftaxis(xaxis));
faxis_ = fftaxisshift(fftaxis(xNaxis));
% fs = fftaxis(xaxis);

if (SPACE.ZeroPad_Factor~=0)&&(sum(~isnan(faxis_)))&&(sum(~isinf(faxis_)))
    fNaxis = faxis_;
else
    fNaxis = faxis;
end

%
dF = abs(fNaxis(end/2+1)-fNaxis(end/2));
ndF = 12*(1+(SPACE.N0/SPACE.N));



%% Determine Save Graphic Names
[TEXT_] = GRAPH_SaveGraphicName(SPACE, TIME, SIM, PULSE, ITP, LASER, DATA, DEBUG);
% [TEXT_] = GRAPH_SaveGraphicName(SaveDirectory, SIM, SPACE, TIME, LASER, PULSE, ITP);


%% Potential Well
if DEBUG.GRAPHS.Potential_VS_Laser
    GRAPH_Potential_VS_Laser(TEXT_, SIM, taxis, xaxis);
end


%% FFT Graph
if DEBUG.GRAPHS.SpectralAnalysis
    GRAPH_MomentumResponse_PopulationLevelNorm_CALC(TEXT_, SIM, PULSE, taxis, xaxis, fNaxis, dF, ndF);
    GRAPH_MomentumResponse_PopulationLevelNorm_VIEW(TEXT_, SIM, PULSE, taxis, xaxis, fNaxis, dF, ndF);
    GRAPH_SpectralAnalysis(TEXT_, SIM, PULSE, taxis, xaxis, fNaxis);
end


%% Population States
if DEBUG.GRAPHS.PopulationState
    GRAPH_PopulationState(TEXT_, SIM, PULSE, taxis, xaxis);
    GRAPH_PopulationStateNorm(TEXT_, SIM, PULSE, taxis, xaxis);
end


%% Centroid
if DEBUG.GRAPHS.CentroidPosition
%     GRAPH_CentroidPosition(TEXT_, SIM, taxis, xaxis);
end


%% Energy Output
if DEBUG.GRAPHS.EnergyOutput_VS_Laser
%     GRAPH_EnergyOutput_VS_Laser(TEXT_, SIM, taxis, xaxis);
end


%% Absorbed Photons
if DEBUG.GRAPHS.PhotonsAbsorbed_VS_Laser
%     GRAPH_PhotonsAbsorbed_VS_Laser(TEXT_, SIM, taxis, xaxis);
end


%% Instantaneous Ionization Rate
if DEBUG.GRAPHS.InstantaneousIonizationRate
%     GRAPH_InstantaneousIonizationRate(TEXT_, SIM, taxis, xaxis);
end


%%
close all hidden;

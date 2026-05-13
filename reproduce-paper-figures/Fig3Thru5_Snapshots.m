%%
% Adjust the frequency below to reproduce Figures 3-5
%
% To recreate these datasets from scratch, in the Julia code run: 
%   1) run_simulation_2D.jl 
%   2) data_snapshot.jl 
%
% Need to run this for freq=2, 16, and 64 Hz and n_on_fiber = 1
% These parameters can be adjusted in PackagesAndParameters.jl
%%
clear; close all; clc;

dataDir = './compressed-datasets/snapshots/';

freq = 2; % Figure 3
% freq = 16; % Figure 4
% freq = 64; % Figure 5
vLines = 1;
modelType='nonlinear';
colOpts = {'start','mid','end'};
rowOpts = [1 2];

if freq == 2
    tauVals = [0.6, 0.85, 1.00; 59.1, 59.35, 59.5];
    firingEvents = [0.6, 59.10];
elseif freq == 16
    tauVals = [0.6, 0.631, 0.65; 59.85, 59.881, 59.9];
    firingEvents = [0.6, 59.85];
elseif freq == 64
    tauVals = [0.6, 0.608, 0.613; 59.975, 59.983, 59.988];
    firingEvents = [0.6, 59.975];
end

%%
f = figure(2);
% Use get(gcf,'Position') to figure out what to put here
for ii = 1:length(rowOpts)

    XCenter = load(strcat(dataDir,'varLocX_snapshot.txt'));
    XSpatial = load(strcat(dataDir,'xVec_snapshot.txt'));

    for jj = 1:length(colOpts)

        % Load the data
        midstring = sprintf('%d_%dline_%s',freq,vLines,modelType);
        endString = sprintf('%s%d',colOpts{jj},rowOpts(ii));

        UPeak = load(strcat(dataDir,...
            sprintf('/scatterU%s_%s.txt',midstring,endString)));
        UDataEq = load(strcat(dataDir,...
            sprintf('scatterUEq%s_%s.txt',midstring,endString)));
        USpatial = load(strcat(dataDir,...
            sprintf('Uxy%s_%s.txt',midstring,endString)));
        
        % Plot the data
        subaxis(2,3,(ii-1)*3+jj,'SpacingVert',0.12,'SpacingHoriz',0.06); hold on;
        % subplot(2,3,(ii-1)*3+jj); hold on;
        Title = sprintf('$\\tau$=%.3f s, firing event: %.3f s',...
            tauVals(ii,jj),firingEvents(ii));
        title(Title,'Interpreter','latex')
        plot(XSpatial,USpatial,'linewidth',1.5)
        plot(XCenter,UPeak,'.','markersize',10)
        for kk = 1:5
            plot(XCenter(kk)+[-2, 2],UDataEq(kk)*[1 1],'k-','linewidth',1.5)
        end
        ylim([-0.1 60])
        xlim([-42 42])
        xLabel = sprintf('$X$ ($\\mu$m)');
        xlabel(xLabel,'Interpreter','latex')
        set(gca,'fontsize',12)
        grid on
        yLabel = sprintf('$\\mathcal{U}(\\tau,X,0)$ [mol/$\\mu$m$^2$]$\\cdot 10^{-22}$');
        ylabel(yLabel,'Interpreter','latex')
    end
end

% set(gcf, 'Position', [-49        1432        1166         685]);

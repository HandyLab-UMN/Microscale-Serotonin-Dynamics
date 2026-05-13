%%
% Uses compressed data to recreate Figure 8
%
% To recreate these datasets from scratch, in the Julia code run (in 
% addition to previous figs):  
%   1) run_simuation_2D_lowerbound
% for freq=2 and n_on_fiber = 1
% 
% Adjust the parameters in PackagesAndParameters.jl
%%
clear; close all; clc;

dataDir = 'compressed-datasets/';

freq = 2;
vLines = 1;
colorScheme = [255, 1, 255; 29, 144, 255;205, 93, 93;...
    255, 187, 62; 128, 0, 128]/255;

%% Load the data
varray = load(strcat(dataDir,sprintf('v_array_freq%d_%dline.txt',freq,vLines)));
muave = load(strcat(dataDir,sprintf('mu_ave%d_%dline_nonlinear.txt',freq,vLines)));
tarray = load(strcat(dataDir,'tarray60s.txt'));
kR = 300/4^2; % time conversion
tarray = tarray/kR;

%% Recreate Figure 8
f = figure(1);
for jj =1:2
    subplot(1,2,jj); hold on;
    for ii = 1:5
        h(ii) = plot(tarray, muave(:,ii)- varray(:,ii),'color',...
            colorScheme(ii,:),'linewidth',1.5);
    end

    if jj == 1
        xlim([2 10])
        legend(h,{'$j=1$','$j=2$','$j=3$','$j=4$','$j=5$'},'Interpreter','latex')
    else
        xlim([57.5 59.3])
    end

    ylabel('E$^f_\tau[\mu_j(\tau)]$ - E$^f_\tau[v_j]$','Interpreter','latex')
    xlabel('$\tau [$s]','Interpreter','latex')
    grid on
end
%%
% Uses compressed data to recreate Figure 10
%
% To recreate these datasets from scratch, in the Julia code run (in 
% addition to previous figs):    
%   1) run_simuation_2D_threefreq
%   2) run_simuation_2D_lowerbound_threefreq
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
muarray = load(strcat(dataDir,'mu_array_16Hz_64Hz_2Hz.txt'));
varray = load(strcat(dataDir,'v_array_16Hz_64Hz_2Hz.txt'));
tarray = load(strcat(dataDir,sprintf('tarray6s.txt')));
kR = 300/4^2; % time conversion
tarray = tarray/kR;

lastidx = length(tarray);
%%
f = figure(1); hold on;%clf; hold on;
h =[];
for ii = 1:5
    h(ii) = plot(tarray,muarray(:,ii),'color',colorScheme(ii,:));
    plot(tarray,varray(1:lastidx,ii),'--','color',colorScheme(ii,:),'linewidth',1.5)
end
legend(h,{'$j=1$','$j=2$','$j=3$','$j=4$','$j=5$'},'Interpreter','latex')
xlabel('\tau [s]')
ylabel('\mu_j(\tau) [mol]\cdot 10^{-22}')
set(gca,'fontsize',16)
xlim([0 6])
ylim([-0.1 20])
grid on
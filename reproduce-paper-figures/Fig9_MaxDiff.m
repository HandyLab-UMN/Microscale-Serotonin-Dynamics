%%
% Uses compressed data to recreate Figure 9 and Table 4 
%
% To recreate these datasets from scratch, in the Julia code run (in 
% addition to previous figs):  
%   1) data_spikemax.jl
% for freq=2 and vLines = 1
%
% Requires the use of outside software to insert subplots
%%
clear; close all; clc;

dataDir = 'compressed-datasets/';

freq = [2 16 64];
vLines = 1;
colorScheme = [255, 1, 255; 29, 144, 255;205, 93, 93;...
    255, 187, 62; 128, 0, 128]/255;

tarray = load(strcat(dataDir,sprintf('tarray60s_full.txt')));
kR = 300/4^2; % time conversion
tarray = tarray/kR;

R_0 = 0.3;
h =1;
conversion = 100/(pi*R_0^2*h);
%%
% The subsets allows the system to equlibrate for 
indexStart = find(tarray>1.6022,1);

figure(1); clf; hold on;
for jj = 1:3

    % Load the data
    max_mu = load(strcat(dataDir,sprintf('max_mu%d_%dline.txt',freq(jj),vLines)));
    max_approx = load(strcat(dataDir,sprintf('max_approx%d_%dline.txt',freq(jj),vLines)));

    for ii = 1:5
        subplot(3,2,(jj-1)*2+1); hold on;
        h(ii) = plot(tarray, max_mu(:,ii),'linewidth',1.5,'Color',colorScheme(ii,:));
        plot(tarray, max_approx(ii)*ones(length(tarray),1),'-.','linewidth',1.5,...
            'Color',colorScheme(ii,:));
        ylabel('$\max\{\mu_j\}^\tau_0$','Interpreter','latex')
        xlabel('$\tau$ [s]','Interpreter','latex')


        if jj == 3
            ylim([0 25])
        end

        subplot(3,2,(jj-1)*2+2); hold on;
        plot(tarray(1:end), max_mu(1:end,ii) - max_approx(ii),'linewidth',1.5, ...
            'Color',colorScheme(ii,:));
        ylabel('$\max\{\mu_j\}^\tau_0-v_{max,j}$','Interpreter','latex')


        fprintf('%0.2f, %0.2f, %0.2f\n',max_mu(end,ii)*conversion,max_approx(ii)*conversion,...
            max_mu(end,ii)*conversion-max_approx(ii)*conversion)
    end

    if jj == 1 
        subplot(3,2,1)
        %legend(h,{'$j=1$','$j=2$','$j=3$','$j=4$','$j=5$'},'Interpreter','latex')
        xlabel('$\tau$ [s]','Interpreter','latex')
    end


    fprintf('--------\n')
end

%%


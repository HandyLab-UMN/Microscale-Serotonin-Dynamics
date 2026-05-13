%%
% Uses compressed data to recreate Figures 6 and 7
%
% To recreate these datasets from scratch, in the Julia code run (in 
% addition to previous figs):
%   2) data_LBSS_singlefreq.jl
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
muarray = load(strcat(dataDir,sprintf('mu_array_freq%d_%dline_nonlinear.txt',freq,vLines)));
muave = load(strcat(dataDir,sprintf('mu_ave%d_%dline_nonlinear.txt',freq,vLines)));
mueq = load(strcat(dataDir,sprintf('mu_eq%d_%dline_nonlinear.txt',freq,vLines)));

tarray = load(strcat(dataDir,sprintf('tarray60s.txt')));
kR = 300/4^2; % time conversion
tarray = tarray/kR;

%% Recreate Figure 6
figure(1); clf;
for jj = 1:2
    subplot(1,2,jj); hold on;
    for ii = 1:5
        h(ii) = plot(tarray,muarray(:,ii),'color',colorScheme(ii,:),...
            'linewidth',1.5);
    end

    if jj == 1
        xlim([0 2])
        legend(h,{'$j=1$','$j=2$','$j=3$','$j=4$','$j=5$'},'Interpreter','latex')
    else
        xlim([57.5 59.5])
        legend off
    end
    ylim([-0.1 12])
    grid on
    set(gca,'fontsize',16)
    xlabel('$\tau\; [$s]','Interpreter','latex')
    ylabel('$\mu_j(\tau)\;[$mol$]\cdot 10^{-22}$','Interpreter','latex')
end

%% Recreate Figure 7
f = figure(2); clf; h =[];
for jj = 1:2

    for ii = 1:5

        if ii < 4
            subplot(3,2,2*(ii-1)+jj); hold on;
        elseif ii == 4
            subplot(3,2,2*(2-1)+jj); hold on;
        elseif ii == 5
            subplot(3,2,2*(1-1)+jj); hold on;
        end

        h(1) = plot(tarray,muarray(:,ii),'color',colorScheme(ii,:),...
            'linewidth',1.5);
        h(2) = plot(tarray,mueq(ii)*ones(length(tarray),1),'--','color',colorScheme(ii,:),...
            'linewidth',1.5);
        h(3) = plot(tarray,muave(:,ii),'-.','color',colorScheme(ii,:),...
           'linewidth',1.5);

        if ii == 5
            legend(h,{'$\mu_j$','E$_\infty^{\rm{f}}[\nu_j]$','E$_t^{\rm{f}}[\mu_j]$'},'Interpreter','latex')
            legend AutoUpdate off
        end

        if ii ==1 || ii==5
            title(sprintf('j=1 and j =5'))
            ylim([0 1])
            yticks([0 0.5 1])
        elseif ii == 2 || ii == 4
            title(sprintf('j=2 and j =4'))
        elseif ii == 3
            title(sprintf('j=3'))
        end

        if (ii == 2 || ii == 4) && jj == 2    
            ylim([0.025 0.045])
            yticks([0.03 0.04])
        elseif (ii == 2 || ii == 4) && jj == 1    
            ylim([0 0.04])
            yticks([0 0.02 0.04])
        elseif ii == 3 && jj == 2
            ylim([0.019 0.0212])
            yticks([0.019 0.02 0.021])
       elseif ii == 3 && jj == 1
            ylim([0 2*1e-2])
            yticks([0 0.01 0.02])
        end

        if jj == 1
            xlim([0 2])
        else
            xlim([57.5 59.5])
            legend off
        end
        grid on
        set(gca,'fontsize',16)
        xlabel('$\tau\; [$s]','Interpreter','latex')
        ylabel('$\mu_j(\tau)\;[$mol$]\cdot 10^{-22}$','Interpreter','latex')
    end
end
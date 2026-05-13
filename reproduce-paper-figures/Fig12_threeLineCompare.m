%%
% Uses compressed data to recreate Figure 12
%
% To recreate these datasets from scratch, in the Julia code run (in 
% addition to previous figs):  
%   1) run_simulation_2D.jl 
%   2) data_LBSS_singlefreq.jl
%   3) data_ECS_LBSS_data.jl
% for freq=2 and n_on_fiber = 1 and 3
%
% Adjust the parameters in PackagesAndParameters.jl
%%
clear; close all; clc;

% dataDir = '/Users/ghandy/Library/CloudStorage/Box-Box/Serotonin_Spatial_Concentrations_Proj/datasets/';

dataDir = 'compressed-datasets/';

freq = 2;
vLines = 3;
savefigOpt = 0;

colorScheme = [255, 1, 255; 29, 144, 255;205, 93, 93;...
    255, 187, 62; 128, 0, 128]/255;

%% Load the data
fprintf('Loading the data\n')
endString = sprintf('%d',freq);
muarray = load(strcat(dataDir,sprintf('mu_array_freq%d_%dline_nonlinear.txt',freq,vLines)));
muave = load(strcat(dataDir,sprintf('mu_ave%d_%dline_nonlinear.txt',freq,vLines)));
mueq = load(strcat(dataDir,sprintf('mu_eq%d_%dline_nonlinear.txt',freq,vLines)));

tarray = load(strcat(dataDir,'tarray60s.txt'));
kR = 300/4^2; % time conversion
tarray = tarray/kR;

Xarray = load(strcat(dataDir,'xarray2D.txt'));
Yarray = load(strcat(dataDir,'yarray2D.txt'));
U = load(strcat(dataDir,sprintf('U2D_freq%d_%dlines_nonlinear.txt',freq,vLines)));

% Load the data from the 1-D trial
muarray1 = load(strcat(dataDir,sprintf('mu_array_freq%d_%dline_nonlinear.txt',freq,1)));
U_oneLine = load(strcat(dataDir,sprintf('U2D_freq%d_%dlines_nonlinear_odd.txt',freq,vLines)));

fprintf('Creating the plots\n')
%%
indcsCell{1} = 1:6271;
indcsCell{2} = 179660:179660+6271;

f = figure(1);
for jj = 1:2
    indcs = indcsCell{jj};

    subplot(2,3,jj); hold on;
    for ii = 1:5
        h(ii) = plot(tarray,muarray(:,ii+5),'color',colorScheme(ii,:),...
            'linewidth',1.5);
    end
    if jj == 1
        xlim([0 2])
        legend(h,{'j=1','j=2','j=3','j=4','j=5'})
    else
        xlim([57.5 59.5])
        legend off
    end
    ylim([-0.1 12])
    grid on
    set(gca,'fontsize',12)
    xlabel('$\tau\; [$s]','Interpreter','latex')
    ylabel('$\mu_{j+5}^{(3)}(\tau)\;[$mol$]\cdot 10^{-22}$','Interpreter','latex')
    
end

%%
for jj = 1:2
    indcs = indcsCell{jj};

    subplot(2,3,jj+3); hold on;
    for ii = 1:5
        h(ii) = plot(tarray,muarray(:,ii+5)-muarray1(:,ii),'color',colorScheme(ii,:),...
            'linewidth',1.5);
    end
    if jj == 1
        xlim([0 2])
    else
        xlim([57.5 59.5])
        legend off
    end
    ylim([-0.01 0.95])
    grid on
    set(gca,'fontsize',12)
    xlabel('$\tau\; [$s]','Interpreter','latex')
    ylabel('$\mu_{j+5}^{(3)}(\tau)-\mu_j^{(1)}(\tau)\;[$mol$]\cdot 10^{-22}$','Interpreter','latex')
end

%%
[XX, YY] = meshgrid(Xarray,Yarray);
downsample = 1:10:1501; % downsame the heatmap

f2 = figure(2); 
surf(XX(downsample,downsample),YY(downsample,downsample),U,'EdgeColor','interp')
view(2)
cSave1 = clim;
axis off
exportgraphics(f2,"ECS_SS.png",ContentType="image")
close(f2)

f3 = figure(3); 
surf(XX(downsample,downsample),YY(downsample,downsample),U-U_oneLine,'EdgeColor','interp')
view(2)
cSave2 = clim;
axis off
exportgraphics(f3,"ECS_SS_Diff.png",ContentType="image")
close(f3)
%%

figure(1);
subplot(2,3,3);
image([-60 60], [-60 60],imread('ECS_SS.png'));
ylabel('$Y\; [\mu$m]','Interpreter','latex')
xlabel('$X\; [\mu$m]','Interpreter','latex')
set(gca,'fontsize',12)
title('$\mathcal{U}_e^{(3)}(X,Y) [$mol/$\mu$m$^2]\cdot 10^{-22}$','Interpreter','latex')
yticks([-50 0 50])
clim(cSave1)
colorbar
axis([-60 60 -60 60])
axis square

subplot(2,3,6);
image([-60 60], [-60 60],imread('ECS_SS_Diff.png'));
ylabel('$Y\; [\mu$m]','Interpreter','latex')
xlabel('$X\; [\mu$m]','Interpreter','latex')
set(gca,'fontsize',12)
title('$(\mathcal{U}_e^{(3)}-\mathcal{U}_e^{(1)})(X,Y) [$mol/$\mu$m$^2]\cdot 10^{-22}$','Interpreter','latex')
yticks([-50 0 50])
clim(cSave2)
colorbar
axis([-60 60 -60 60])
axis square

%%
if savefigOpt == 1
    fprintf('Temporarily paused to adjust figure\n')
    pause()
    exportgraphics(f,"test4.pdf",ContentType="vector")
end

%%

f = figure(2); clf; h =[];
for jj = 1:2
    indcs = indcsCell{jj};
    for ii = 1:5
        
        if ii < 4
            subplot(3,2,2*(ii-1)+jj); hold on;
        elseif ii == 4
            subplot(3,2,2*(2-1)+jj); hold on;
        elseif ii == 5
            subplot(3,2,2*(1-1)+jj); hold on;
        end

        h(1) = plot(tarray,muarray(:,ii+5),'color',colorScheme(ii,:),...
            'linewidth',1.5);
        h(2) = plot(tarray,mueq(ii+5)*ones(length(tarray),1),'--','color',colorScheme(ii,:),...
            'linewidth',1.5);
        h(3) = plot(tarray,muave(:,ii+5),'-.','color',colorScheme(ii,:),...
           'linewidth',1.5);
        

        if ii == 5
            legend(h,{'$\mu_j$','E$_\infty^{\rm{f}}[\nu_j]$','E$_t^{\rm{f}}[\mu_j]$'},'Interpreter','latex')
            legend AutoUpdate off
        end
        
        if ii ==1 || ii==5
            title('$j$=6 and $j$ =10','Interpreter','latex')
            ylim([0 1])
            yticks([0 0.5 1])
        elseif ii == 2 || ii == 4
            title('$j$=7 and $j$ =9','Interpreter','latex')
        elseif ii == 3
            title('$j$=8','Interpreter','latex')
        end

        if (ii == 2 || ii == 4) && jj == 2    
            ylim([0.06 0.122])
            yticks([0.08 0.12])
        elseif (ii == 2 || ii == 4) && jj == 1    
            ylim([0 0.1])
            yticks([0 0.04 0.08])
        elseif ii == 3 && jj == 2
            ylim([0.045 0.052])
            yticks([0.045 0.05])
       elseif ii == 3 && jj == 1
            ylim([0 5*1e-2])
            yticks([0 0.02 0.04])
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
        ylabel('$\mu_j^{(3)}(\tau)\;[$mol$]\cdot 10^{-22}$','Interpreter','latex')
    end

  
end

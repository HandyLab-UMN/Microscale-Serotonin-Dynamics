%%
% Uses compressed data to recreate Figure 13
%
% To recreate these datasets from scratch, in the Julia code run (in 
% addition to previous figs):  
%   1) run_simulation_2D.jl 
%   2) data_LBSS_singlefreq.jl
%   3) data_ECS_LBSS_data.jl
% for freq=2, n_on_fiber = 3 and modelType = 'linear' and 'nonlinear'
%
% Adjust the parameters in PackagesAndParameters.jl
%%
clear; close all; clc;

dataDir = 'compressed-datasets/';

freq = 2;
vLines = 3;
savefigOpt = 0;

colorScheme = [255, 1, 255; 29, 144, 255;205, 93, 93;...
    255, 187, 62; 128, 0, 128]/255;

fprintf('Loading the data\n')
endString = sprintf('%d',freq);
muarray = load(strcat(dataDir,sprintf('mu_array_freq%d_%dline_linear.txt',freq,vLines)));

tarray = load(strcat(dataDir,'tarray60s.txt'));
kR = 300/4^2; % time conversion
tarray = tarray/kR;

muarray_nl = load(strcat(dataDir,sprintf('mu_array_freq%d_%dline_nonlinear.txt',freq,vLines)));

U_l = load(strcat(dataDir,sprintf('U2D_freq%d_%dlines_linear.txt',freq,vLines)));
U_nl = load(strcat(dataDir,sprintf('U2D_freq%d_%dlines_nonlinear.txt',freq,vLines)));

Xarray = load(strcat(dataDir,'xarray2D.txt'));
Yarray = load(strcat(dataDir,'yarray2D.txt'));

fprintf('Creating the plots\n')
%%
f=figure(1);
for jj = 1:2

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
    ylabel('$\mu_{L,j+5}^{(3)}(\tau)\;[$mol$]\cdot 10^{-22}$','Interpreter','latex')

end

%%
for jj = 1:2
    subplot(2,3,jj+3); hold on;
    for ii = 1:5
        h(ii) = plot(tarray,muarray(:,ii+5)-muarray_nl(:,ii+5),'color',colorScheme(ii,:),...
            'linewidth',1.5);
    end
    if jj == 1
        xlim([0 2])
    else
        xlim([57.5 59.5])
        legend off
    end
    grid on
    set(gca,'fontsize',12)
    xlabel('$\tau\; [$s]','Interpreter','latex')
    ylabel('$\mu_{L,j+5}^{(3)}(\tau)-\mu_j^{(3)}(\tau)\;[$mol$]\cdot 10^{-22}$','Interpreter','latex')

end

%% Create these subplots first, save as PNG then insert into plot
% Helps with save the whole figure as a PDF
[XX, YY] = meshgrid(Xarray,Yarray);
downsample = 1:10:1501; % downsame the heatmap a bit
f2 = figure(2); 
surf(XX(downsample,downsample),YY(downsample,downsample),U_l,'EdgeColor','interp')
view(2)
axis off
cSave1 = clim;
exportgraphics(f2,"ECS_SS.png",ContentType="image")
close(f2);

f3 = figure(3); 
surf(XX(downsample,downsample),YY(downsample,downsample),U_l-U_nl,'EdgeColor','interp')
view(2)
axis off
cSave2 = clim;
exportgraphics(f3,"ECS_SS_Diff.png",ContentType="image")
close(f3);

%%
figure(1);
subplot(2,3,3);
image([-60 60], [-60 60],imread('ECS_SS.png'));
ylabel('$Y\; [\mu$m]','Interpreter','latex')
xlabel('$X\; [\mu$m]','Interpreter','latex')
set(gca,'fontsize',12)
title('$\mathcal{U}_{L,e}^{(3)}(X,Y) [$mol/$\mu$m$^2]\cdot 10^{-22}$','Interpreter','latex')
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
title('$(\mathcal{U}_{L,e}^{(3)}-\mathcal{U}_e^{(3)})(X,Y) [$mol/$\mu$m$^2]\cdot 10^{-22}$','Interpreter','latex')
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
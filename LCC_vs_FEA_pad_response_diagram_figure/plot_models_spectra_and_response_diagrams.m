close all;
clearvars;
clc;

%% load results

load('DOEs_FEA_ellipsoid_G40P80G40_W10xH20.mat'); % FEA results
d_vec_LUT = d_vec;
dR_vec_LUT = dR_vec;
load('DOEs_LCC_voxel_ellipsoid_G40P80G40_W10xH20.mat'); % LCC voxel results
load('DOEs_LCC_slice_ellipsoid_G40P80G40_W10xH20.mat'); % LCC slice results

addpath('../functions/');

%% compute R + R/C results

pad_to_par_fact = 0.1;
C_eq = c * lp * wc * pad_to_par_fact;

Z_none_RRC = 2 * rho * lg / (wc * hc) + 1./( 1/(rho * lp / (wc*hc)) + 1i*2*pi*f_vec*C_eq );
Z_gaps_RRC = nan([length(d_vec), length(f_vec)]);
Z_pads_RRC = nan([length(d_vec), length(f_vec)]);
for i_d = 1:length(d_vec)
    d_cell = d_vec(i_d);
    dR = max([0, interp1(d_vec_LUT, dR_vec_LUT, d_cell, "linear", "extrap")]); % look up DeltaR in gap from COMSOL dataset
    Z_gaps_RRC(i_d,:) = 2 * rho * lg / (wc * hc) + 1./( 1/(rho * lp / (wc*hc)) + 1i*2*pi*f_vec*C_eq ) + dR;
    Z_pads_RRC(i_d,:) = 2 * rho * lg / (wc * hc) + 1./( 1/(rho * lp / (wc*hc) + dR) + 1i*2*pi*f_vec*C_eq );
end

%% plot R + R/C and LCC / FEA spectra

ind_plot = 12;

% plot R + R/C spectra
fig1 = figure(1); clf; hold on;
plotBodeNyquist(f_vec/1e3, Z_none_RRC/1e6, '-');
plotBodeNyquist(f_vec/1e3, Z_gaps_RRC(ind_plot,:)/1e6, '-');
plotBodeNyquist(f_vec/1e3, Z_pads_RRC(ind_plot,:)/1e6, '-');

subplot(2,2,1);
title('R + R/C model')
axis([0.1,1e3,0.24,0.56]);
grid off;
set(gca, xtick=[0.1,1,10,100,1000], xticklabels=["0.1","1","10","100","1000"]);
set(gca, yscale="linear", ytick=0.25:0.05:0.55, yticklabels=["","0.3","","0.4","","0.5",""]);
xlabel('f [kHz]');
ylabel('|Z| [M\Omega]');
legend(["No cell", "Cell in gap", "Cell over pad"], location="sw");
subplot(2,2,3);
axis([0.1,1e3,-24,0]);
grid off;
set(gca, xtick=[0.1,1,10,100,1000], xticklabels=["0.1","1","10","100","1000"]);
set(gca, ytick=-20:5:0, yticklabels=["-20","","-10","","0"], xaxislocation="origin");
xlabel("");
subplot(1,2,2);
axis([0.24,0.56,0,0.16]);
grid off;
set(gca, xtick=0.25:0.05:0.55, xticklabels=["","0.3","","0.4","","0.5"]);
set(gca, ytick=0:0.05:0.15, yticklabels=["0","","0.1",""]);
xlabel('Re(Z) [M\Omega]');
ylabel('-Im(Z) [M\Omega]');

savename = 'output/RRC_spectra';
saveas(gcf, [savename '.fig']);
saveas(gcf, [savename '.png']);
saveas(gcf, [savename '.pdf']);

% plot LCC and FEA spectra
fig2 = figure(2); clf; hold on;
plotBodeNyquist(f_vec/1e3, does_voxels(ind_plot).Znone/1e6, '-');
plotBodeNyquist(f_vec/1e3, does_voxels(ind_plot).Zgap/1e6, '-');
plotBodeNyquist(f_vec/1e3, does_voxels(ind_plot).Zpad/1e6, '-');
subplot(2,2,1); set(gca, ColorOrderIndex=1); subplot(2,2,3); set(gca, ColorOrderIndex=1); subplot(1,2,2); set(gca, ColorOrderIndex=1); 
plotBodeNyquist(f_vec/1e3, does(ind_plot).Znone/1e6, '--');
plotBodeNyquist(f_vec/1e3, does(ind_plot).Zgap/1e6, '--');
plotBodeNyquist(f_vec/1e3, does(ind_plot).Zpad/1e6, '--');

subplot(2,2,1);
title('LCC model / FEA')
axis([0.1,1e3,0.24,0.56]);
grid off;
set(gca, xtick=[0.1,1,10,100,1000], xticklabels=["0.1","1","10","100","1000"]);
set(gca, yscale="linear", ytick=0.25:0.05:0.55, yticklabels=["","0.3","","0.4","","0.5",""]);
xlabel('f [kHz]');
ylabel('|Z| [M\Omega]');
legend(["No cell", "Cell in gap", "Cell over pad"], location="sw");
subplot(2,2,3);
axis([0.1,1e3,-24,0]);
grid off;
set(gca, xtick=[0.1,1,10,100,1000], xticklabels=["0.1","1","10","100","1000"]);
set(gca, ytick=-20:5:0, yticklabels=["-20","","-10","","0"], xaxislocation="origin");
xlabel("");
subplot(1,2,2);
axis([0.24,0.56,0,0.16]);
grid off;
set(gca, xtick=0.25:0.05:0.55, xticklabels=["","0.3","","0.4","","0.5"]);
set(gca, ytick=0:0.05:0.15, yticklabels=["0","","0.1",""]);
xlabel('Re(Z) [M\Omega]');
ylabel('-Im(Z) [M\Omega]');

savename = 'output/LCC_FEA_spectra';
saveas(gcf, [savename '.fig']);
saveas(gcf, [savename '.png']);
saveas(gcf, [savename '.pdf']);

%% plot model comparison spectra

ind_plot = 12;

% plot Znone spectra
fig3 = figure(3); clf; hold on;
plotBodeNyquist(f_vec/1e3, Z_none_RRC/1e6, '-');
plotBodeNyquist(f_vec/1e3, Z_none_slices/1e6, '--');
plotBodeNyquist(f_vec/1e3, does_voxels(ind_plot).Znone/1e6, '-');
plotBodeNyquist(f_vec/1e3, does(ind_plot).Znone/1e6, '--');

subplot(2,2,1);
title('Z_{none}')
axis([0.1,1e3,0.24,0.56]);
grid off;
set(gca, xtick=[0.1,1,10,100,1000], xticklabels=["0.1","1","10","100","1000"]);
set(gca, yscale="linear", ytick=0.25:0.05:0.55, yticklabels=["","0.3","","0.4","","0.5",""]);
xlabel('f [kHz]');
ylabel('|Z| [M\Omega]');
legend(["R + R/C","LCC slice","LCC voxel","FEA"], location="sw");
subplot(2,2,3);
axis([0.1,1e3,-24,0]);
grid off;
set(gca, xtick=[0.1,1,10,100,1000], xticklabels=["0.1","1","10","100","1000"]);
set(gca, ytick=-20:5:0, yticklabels=["-20","","-10","","0"], xaxislocation="origin");
xlabel("");
subplot(1,2,2);
axis([0.24,0.56,0,0.16]);
grid off;
set(gca, xtick=0.25:0.05:0.55, xticklabels=["","0.3","","0.4","","0.5"]);
set(gca, ytick=0:0.05:0.15, yticklabels=["0","","0.1",""]);
xlabel('Re(Z) [M\Omega]');
ylabel('-Im(Z) [M\Omega]');

savename = 'output/model_comparison_spectra_Znone';
saveas(gcf, [savename '.fig']);
saveas(gcf, [savename '.png']);
saveas(gcf, [savename '.pdf']);

% plot Zgap spectra
fig4 = figure(4); clf; hold on;
plotBodeNyquist(f_vec/1e3, Z_gaps_RRC(ind_plot,:)/1e6, '-');
plotBodeNyquist(f_vec/1e3, Z_gaps_slices(ind_plot,:)/1e6, '--');
plotBodeNyquist(f_vec/1e3, does_voxels(ind_plot).Zgap/1e6, '-');
plotBodeNyquist(f_vec/1e3, does(ind_plot).Zgap/1e6, '--');

subplot(2,2,1);
title('Z_{gap}')
axis([0.1,1e3,0.24,0.56]);
grid off;
set(gca, xtick=[0.1,1,10,100,1000], xticklabels=["0.1","1","10","100","1000"]);
set(gca, yscale="linear", ytick=0.25:0.05:0.55, yticklabels=["","0.3","","0.4","","0.5",""]);
xlabel('f [kHz]');
ylabel('|Z| [M\Omega]');
legend(["R + R/C","LCC slice","LCC voxel","FEA"], location="sw");
subplot(2,2,3);
axis([0.1,1e3,-24,0]);
grid off;
set(gca, xtick=[0.1,1,10,100,1000], xticklabels=["0.1","1","10","100","1000"]);
set(gca, ytick=-20:5:0, yticklabels=["-20","","-10","","0"], xaxislocation="origin");
xlabel("");
subplot(1,2,2);
axis([0.24,0.56,0,0.16]);
grid off;
set(gca, xtick=0.25:0.05:0.55, xticklabels=["","0.3","","0.4","","0.5"]);
set(gca, ytick=0:0.05:0.15, yticklabels=["0","","0.1",""]);
xlabel('Re(Z) [M\Omega]');
ylabel('-Im(Z) [M\Omega]');

savename = 'output/model_comparison_spectra_Zgap';
saveas(gcf, [savename '.fig']);
saveas(gcf, [savename '.png']);
saveas(gcf, [savename '.pdf']);

% plot Zpad spectra
fig5 = figure(5); clf; hold on;
plotBodeNyquist(f_vec/1e3, Z_pads_RRC(ind_plot,:)/1e6, '-');
plotBodeNyquist(f_vec/1e3, Z_pads_slices(ind_plot,:)/1e6, '--');
plotBodeNyquist(f_vec/1e3, does_voxels(ind_plot).Zpad/1e6, '-');
plotBodeNyquist(f_vec/1e3, does(ind_plot).Zpad/1e6, '--');

subplot(2,2,1);
title('Z_{pad}')
axis([0.1,1e3,0.24,0.56]);
grid off;
set(gca, xtick=[0.1,1,10,100,1000], xticklabels=["0.1","1","10","100","1000"]);
set(gca, yscale="linear", ytick=0.25:0.05:0.55, yticklabels=["","0.3","","0.4","","0.5",""]);
xlabel('f [kHz]');
ylabel('|Z| [M\Omega]');
legend(["R + R/C","LCC slice","LCC voxel","FEA"], location="sw");
subplot(2,2,3);
axis([0.1,1e3,-24,0]);
grid off;
set(gca, xtick=[0.1,1,10,100,1000], xticklabels=["0.1","1","10","100","1000"]);
set(gca, ytick=-20:5:0, yticklabels=["-20","","-10","","0"], xaxislocation="origin");
xlabel("");
subplot(1,2,2);
axis([0.24,0.56,0,0.16]);
grid off;
set(gca, xtick=0.25:0.05:0.55, xticklabels=["","0.3","","0.4","","0.5"]);
set(gca, ytick=0:0.05:0.15, yticklabels=["0","","0.1",""]);
xlabel('Re(Z) [M\Omega]');
ylabel('-Im(Z) [M\Omega]');

savename = 'output/model_comparison_spectra_Zpad';
saveas(gcf, [savename '.fig']);
saveas(gcf, [savename '.png']);
saveas(gcf, [savename '.pdf']);

%% plot pad vs. gap response diagrams

Z_none_comsol = [does(1).Znone].';
Z_gaps_comsol = [does(:).Zgap].';
Z_pads_comsol = [does(:).Zpad].';

X_sim_dR = (dR_vec(:)./abs(Z_none_comsol))*100;

Z_none_spice = does_voxels(1).Znone;
Z_gaps_spice = reshape([does_voxels(:).Zgap].', length(f_vec), []).';
Z_pads_spice = reshape([does_voxels(:).Zpad].', length(f_vec), []).';

X_SPICE = (abs(Z_gaps_spice)./abs(Z_none_spice)-1)*100;
Y_SPICE = (abs(Z_pads_spice)./abs(Z_none_spice)-1)*100;
X_SPICE_dR = (dR_vec(:)./abs(Z_none_spice))*100;

fig6 = figure(6); clf; hold on;
scalar = 20;
init_angles_deg = 180/pi*atan2(Y_sim(1,:)*scalar, X_sim(1,:));
init_angles_deg(f_vec > 200e3) = nan;
angle_0 = 90;
[~, ind_freq_1] = min(abs( init_angles_deg - angle_0*0.9 )); % slope 36 deg
[~, ind_freq_2] = min(abs( init_angles_deg - angle_0*0.8 )); % slope 27 deg
[~, ind_freq_3] = min(abs( init_angles_deg - angle_0*0.7 )); % slope 18 deg
[~, ind_freq_4] = min(abs( init_angles_deg - angle_0*0.45 )); % slope 9 deg
[~, ind_freq_5] = min(abs( init_angles_deg - 0 )); % slope 0 deg
[min_slope, ind_freq_6] = min(init_angles_deg); % most negative slope
[~, ind_freq_7] = min(abs( init_angles_deg - min_slope/2)); % half of min slope

f_targets_comsol = f_vec([ind_freq_1, ind_freq_2, ind_freq_3, ind_freq_4, ind_freq_5, ind_freq_6, ind_freq_7]);
f_targets_spice = f_targets_comsol;

% plot SPICE results
subplot(1,2,1); hold on;
plot([0; X_SPICE_dR(:,1)] / scalar, [0; Y_SPICE(:,1)], 'ks-', DisplayName=sprintf('LCC: %0.1fkHz', f_vec(1)/1e3));
subplot(1,2,2); hold on;
plot([0; X_SPICE(:,1)] / scalar, [0; Y_SPICE(:,1)], 'ks-', DisplayName=sprintf('LCC: %0.1fkHz', f_vec(1)/1e3));
subplot(1,2,1); set(gca, ColorOrderIndex=1); subplot(1,2,2); set(gca, ColorOrderIndex=1);
for i = 1:length(f_targets_spice)
    [~, ind] = min(abs(f_vec-f_targets_spice(i)));
    subplot(1,2,1); hold on;
    plot([0; X_SPICE_dR(:,ind)] / scalar, [0; Y_SPICE(:,ind)], 's-', DisplayName=sprintf('LCC: %0.1fkHz', f_vec(ind)/1e3));
    subplot(1,2,2); hold on;
    plot([0; X_SPICE(:,ind)] / scalar, [0; Y_SPICE(:,ind)], 's-', DisplayName=sprintf('LCC: %0.1fkHz', f_vec(ind)/1e3));
end

% plot COMSOL results
subplot(1,2,1); hold on;
plot([0; X_sim_dR(:,1)] / scalar, [0; Y_sim(:,1)], 'kd-', DisplayName=sprintf('FEA: %0.1fkHz', f_vec(1)/1e3));
subplot(1,2,2); hold on;
plot([0; X_sim(:,1)] / scalar, [0; Y_sim(:,1)], 'kd-', DisplayName=sprintf('FEA: %0.1fkHz', f_vec(1)/1e3));
subplot(1,2,1); set(gca, ColorOrderIndex=1); subplot(1,2,2); set(gca, ColorOrderIndex=1);
for i = 1:length(f_targets_comsol)
    [~, ind] = min(abs(f_vec-f_targets_comsol(i)));
    subplot(1,2,1); hold on;
    plot([0; X_sim_dR(:,ind)] / scalar, [0; Y_sim(:,ind)], 'd-', DisplayName=sprintf('FEA: %0.1fkHz', f_vec(ind)/1e3));
    subplot(1,2,2); hold on;
    plot([0; X_sim(:,ind)] / scalar, [0; Y_sim(:,ind)], 'd-', DisplayName=sprintf('FEA: %0.1fkHz', f_vec(ind)/1e3));
end

subplot(1,2,1);
axis equal;
xl = 20;
set(gca, xlim=[0, xl/scalar]);
set(gca, ylim=[-1,1]);
set(gca, XTick=(0:5:xl)./scalar);
set(gca, XTickLabels=get(gca, 'XTick')*scalar);
set(gca, XAxisLocation="origin");
set(gca, YTick=-1:0.25:1);
xlabel("\DeltaR / |Z_{none}| (%)");
ylabel("( |Z_{pad}| - |Z_{none}| ) / |Z_{none}| (%)");
legend(location="eastoutside");
subplot(1,2,2);
axis equal;
set(gca, xlim=[0, xl/scalar]);
set(gca, ylim=[-1,1]);
set(gca, XTick=(0:5:xl)./scalar);
set(gca, XTickLabels=get(gca, 'XTick')*scalar);
set(gca, YTick=-1:0.25:1);
set(gca, XAxisLocation="origin");
xlabel("( |Z_{gap}| - |Z_{none}| ) / |Z_{none}| (%)");
ylabel("( |Z_{pad}| - |Z_{none}| ) / |Z_{none}| (%)");
legend(location="eastoutside");

set(gcf, WindowState="maximized");

savename = 'output/model_comparison_response_diagram';
saveas(gcf, [savename '.fig']);
saveas(gcf, [savename '.png']);
saveas(gcf, [savename '.pdf']);

%% save results

save('plot_models_spectra_and_response_diagrams.mat');

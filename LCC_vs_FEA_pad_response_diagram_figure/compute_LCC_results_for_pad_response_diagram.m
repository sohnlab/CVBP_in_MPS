close all;
clearvars;
clc;

%% load data

load('DOEs_FEA_ellipsoid_G40P80G40_W10xH20.mat');
d_vec_LUT = d_vec;
dR_vec_LUT = dR_vec;

addpath('../functions/');

%% set model parameters

l_vec = [lg, lp, lg];
c = 0.2;
f_vec = [0.001, logspace(1,6,501), 1.885e6];
node_num_vec = [nan, 2, nan];
x_cell_none = -Inf;
x_cell_gap = l_vec(1)/2;
x_cell_pad = l_vec(1) + l_vec(2)/2;

%% compute LCC slice results

n_steps = 250;

Z_none_slices = SPICE_VOL_MPS_RCPE_slices(l_vec, wc, hc, node_num_vec, rho, c, 1, 0, x_cell_none, 0, f_vec, n_steps);
Z_gaps_slices = nan([length(d_vec), length(f_vec)]);
Z_pads_slices = nan([length(d_vec), length(f_vec)]);
tic;
count = 0;
for i_d = 1:length(d_vec)
    d_cell = d_vec(i_d);
    dR = max([0, interp1(d_vec_LUT, dR_vec_LUT, d_cell, "linear", "extrap")]); % look up DeltaR in gap from COMSOL dataset
    Z_gaps_slices(i_d,:) = SPICE_VOL_MPS_RCPE_slices(l_vec, wc, hc, node_num_vec, rho, c, 1, d_cell, x_cell_gap, dR, f_vec, n_steps);
    Z_pads_slices(i_d,:) = SPICE_VOL_MPS_RCPE_slices(l_vec, wc, hc, node_num_vec, rho, c, 1, d_cell, x_cell_pad, dR, f_vec, n_steps);
    
    count = count + 1;
    fprintf('Time elapsed: %0.0f min\nTime remaining: %0.0f min\n', toc/60, toc/60/count*length(d_vec)-toc/60);
end

%% save LCC slice results

save('DOEs_LCC_slice_ellipsoid_G40P80G40_W10xH20.mat', 'Z_none_slices', 'Z_gaps_slices', 'Z_pads_slices', 'd_vec', 'dR_vec', 'lg', 'lp', 'wc', 'hc', 'rho', 'c', 'f_vec');

%% compute LCC voxel results

dxyz = 20/11 .* 1e-6;

[Z_none, result_none] = SPICE_VOL_MPS_RCPE_voxels(l_vec, wc, hc, node_num_vec, rho, c, 1, 0, x_cell_none, 0, f_vec, dxyz, true);
tic;
count = 0;
for i = 1:length(d_vec)
    [Z_gap, result_gap] = SPICE_VOL_MPS_RCPE_voxels(l_vec, wc, hc, node_num_vec, rho, c, 1, d_vec(i), x_cell_gap, dR_vec(i), f_vec, dxyz, true);
    [Z_pad, result_pad] = SPICE_VOL_MPS_RCPE_voxels(l_vec, wc, hc, node_num_vec, rho, c, 1, d_vec(i), x_cell_pad, dR_vec(i), f_vec, dxyz, true);
    doe.name = sprintf('G40-P80-G40_W10xH20_%0.1fum-cell_%0.1fum-voxel', d_vec(i)*1e6, dxyz*1e6);
    doe.diam_um = d_vec(i)*1e6;
    doe.voxel_um = dxyz*1e6;
    doe.rho = rho;
    doe.Cs = c;
    doe.Znone = Z_none;
    doe.x_grid_vec = result_none.x_grid_vec;
    doe.y_grid_vec = result_none.y_grid_vec;
    doe.V_grid_mat_none = result_none.V_grid_mat;
    doe.V_grid_mat_gap = result_gap.V_grid_mat;
    doe.V_grid_mat_pad = result_pad.V_grid_mat;
    doe.Zgap = Z_gap;
    doe.Zpad = Z_pad;
    doe.freq = f_vec;
    doe.ratio_pad_none = abs(Z_pad)./abs(Z_none);
    doe.ratio_gap_none = abs(Z_gap)./abs(Z_none);
    doe.ratio_dip_peak = (1-abs(Z_pad)./abs(Z_none))./(1-abs(Z_gap)./abs(Z_none));
    does_voxels(i) = doe;

    count = count + 1;
    fprintf('Time elapsed: %0.0f min\nTime remaining: %0.0f min\n', toc/60, toc/60/count*length(d_vec)-toc/60);
end

%% save LCC voxel results

save('DOEs_LCC_voxel_ellipsoid_G40P80G40_W10xH20.mat', 'does_voxels', 'd_vec', 'dR_vec', 'lg', 'lp', 'wc', 'hc', 'rho', 'c', 'f_vec');

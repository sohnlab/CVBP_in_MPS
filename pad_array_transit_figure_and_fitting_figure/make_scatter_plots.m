dbstop if error

close all;
clearvars;
clc;

%% run this code block once to collect all detections into one struct

p = 'analyzed_detections/';
dr = dir(p);
fns = {dr.name};
mask = contains(fns, '.mat');
fns = fns(mask);
D = [];
for i = 1:length(fns)
    fn = fns{i};
    load([p fn]);
    [good_detections.fn] = deal(fn);
    D = [D, good_detections];
end
save('pad_array_detections.mat', 'D', '-V7.3');

%% extract data from struct into arrays for ease of plotting

Z_GAPS = []; Z_PADS = []; Z_NONE = []; RMSE = []; RMSE_ = []; MAE = []; MAE_ = []; LINF = []; LINF_ = []; PAD_DURS = []; GAP_DURS = [];
for i = 1:length(D)
    %%% FIX MISTAKES
    D(i).fitting.resid = -D(i).fitting.resid; % FIX MISTAKE
    D(i).fitting.z_mat = D(i).fitting.z_fit + D(i).fitting.resid;
    resid = D(i).fitting.resid;
    z_mat = D(i).fitting.z_mat;
    D(i).fitting.rmse = sqrt(mean(mean(resid.^2)));
    D(i).fitting.mae = mean(mean(abs(resid)));
    D(i).fitting.Linf = max(max(abs(resid)));
    D(i).fitting.rmse_ = D(i).fitting.rmse / sqrt(mean(mean(z_mat.^2)));
    D(i).fitting.mae_ = D(i).fitting.mae / mean(mean(abs(z_mat)));
    D(i).fitting.Linf_ = D(i).fitting.Linf / max(max(abs(z_mat)));
    %%% fix b_vec
    Z_none_left = D(i).b_vec;
    Z_none_right = D(i).b_vec;
    for i_f = 1:size(D(i).Z_dcml, 2)
        zout1 = D(i).Z_dcml(round(length(D(i).Z_dcml)*0.2+1:length(D(i).Z_dcml)*0.4),i_f);
        zout2 = D(i).Z_dcml(round(length(D(i).Z_dcml)*0.6+1:length(D(i).Z_dcml)*0.8),i_f);
        Z_none_left(i_f) = robust_median(zout1);
        Z_none_right(i_f) = robust_median(zout2);
    end
    D(i).Z_none_left = Z_none_left;
    D(i).Z_none_right = Z_none_right;
    %%%
    d = D(i);
    Z_GAPS(:,:,i) = d.Z_gaps;
    Z_PADS(:,:,i) = d.Z_pads;
    Z_NONE_LEFT(1,:,i) = d.Z_none_left;
    Z_NONE_RIGHT(1,:,i) = d.Z_none_right;
    RMSE(i) = d.fitting.rmse;
    RMSE_(i) = d.fitting.rmse;
    MAE(i) = d.fitting.mae;
    MAE_(i) = d.fitting.mae_;
    LINF(i) = d.fitting.Linf;
    LINF_(i) = d.fitting.Linf_;
    PAD_DURS(:,i) = d.edge_times_fit(3:2:end) - d.edge_times_fit(2:2:end-1);
    GAP_DURS(:,i) = d.edge_times_fit(2:2:end) - d.edge_times_fit(1:2:end-1);
end

save('pad_array_detections.mat', 'D', 'Z_GAPS', 'Z_PADS', 'Z_NONE_LEFT', 'Z_NONE_RIGHT', 'RMSE', 'RMSE_', 'MAE', 'MAE_', 'LINF', 'LINF_', 'GAP_DURS', 'PAD_DURS', '-V7.3');

%% load data

load('pad_array_detections.mat');

%% plot scatter plot    

min_gap_dur = 1e-3*2;
min_pad_dur = 0.025e-3*2;
max_rmse_ = prctile(RMSE_, 99.99);
max_mae_ = prctile(MAE_, 99.99);
max_Linf_ = prctile(LINF_, 99.99);
max_none_diff = 0.01;

mask_gap_dur = min(GAP_DURS, [], 1) >= min_gap_dur;
mask_pad_dur = min(PAD_DURS, [], 1) >= min_pad_dur;
mask_rmse_ = RMSE_ <= max_rmse_;
mask_mae_ = MAE_ <= max_mae_;
mask_Linf_ = LINF_ <= max_Linf_;
Z_NONE_AVG = abs(Z_NONE_LEFT)/2+abs(Z_NONE_RIGHT)/2;
mask_none_diff = squeeze(all(abs(abs(Z_NONE_LEFT)-abs(Z_NONE_RIGHT))./Z_NONE_AVG <= max_none_diff, 2)).';

mask = contains({D.fn}, '79padarr');
% mask = contains({D.fn}, '07padarr');
mask = [mask; mask_gap_dur];
mask = [mask; mask_pad_dur];
mask = [mask; mask_rmse_];
mask = [mask; mask_mae_];
mask = [mask; mask_Linf_];
mask = [mask; mask_none_diff];
mask = all(mask, 1);

f_vec = D(1).freq;
f_mask = f_vec < 400e3;

X_data_all = ( (abs(Z_GAPS(1:end-1,f_mask,:))/2+abs(Z_GAPS(2:end,f_mask,:))/2) ./ abs(Z_NONE_AVG(1,f_mask,:)) - 1 ) * 100;
Y_data_all = ( abs(Z_PADS(:,f_mask,:)) ./ abs(Z_NONE_AVG(1,f_mask,:)) - 1 ) * 100;

mask = mask & squeeze(all(all(Y_data_all < X_data_all*1.25, 2), 1)).' & squeeze(all(all(Y_data_all > X_data_all*-0.25, 2), 1)).';
sum(mask)

figure(1); clf; hold on;
plotAll(Z_NONE_AVG(1,f_mask,mask), Z_GAPS(:,f_mask,mask), Z_PADS(:,f_mask,mask), f_vec(f_mask), '.');

subplot(2,4,8);
title('narrower pad array');
shg;

%% try forward modeling

l_chan = [40, 10, 40, 20, 40, 40, 40, 80, 40, 160, 40, 320, 40, 640, 40] * 1e-6;
pad_mult = 1;
pad_add = 0e-6;
l_chan(2:2:end-1) = l_chan(2:2:end-1) * pad_mult + pad_add;
l_extra = sum(l_chan) - 1590e-6;
l_chan(1:2:end) = l_chan(1:2:end) - l_extra / length(l_chan(1:2:end));
if any(l_chan < 0)
    warning('objective_function: l_chan cannot be negative!');
    l_chan(l_chan < 0) = 1e-9;
    l_chan = l_chan ./ sum(l_chan) .* 1590e-6;
    l_chan(l_chan < 0) = 1e-9;
    l_chan = l_chan ./ sum(l_chan) .* 1590e-6;
    l_chan(l_chan < 0) = 1e-9;
    l_chan = l_chan ./ sum(l_chan) .* 1590e-6;
end
w = 20.0e-6;
h = 21.5e-6;
rho = 0.625;
q_vec = [1, nan, nan, nan, nan, 0.1, nan, 0.1, nan, 0.133, nan, 0.0665, nan, 0.0128, nan, 0.0147, nan, 0.0304, nan, nan, nan, nan, 1];
a_vec = [1, nan, nan, nan, nan, 1.0, nan, 1.0, nan, 0.887, nan, 0.900, nan, 0.999, nan, 0.962, nan, 0.896, nan, nan, nan, nan, 1];
load('../LCC_vs_FEA_pad_response_diagram_figure/DOEs_FEA_ellipsoid_G40P80G40_W10xH20.mat', 'd_vec', 'dR_vec');
d_vec_LUT = d_vec;
dR_vec_LUT = dR_vec;

n_diams = 5;
d_vec = linspace(0, (21.5e-6)^3, n_diams).^(1/3);
d_vec(1) = [];
dR_vec = interp1(d_vec_LUT, dR_vec_LUT, d_vec, "linear", "extrap"); % interp from LUT to get dR_vec for the d_vec above

%% compute forward model

tic;
[Z_NONE_sim, Z_GAPS_sim, Z_PADS_sim, X_sim, Y_sim] = pad_array_vol_forward_model(l_chan, w, h, rho, q_vec, a_vec, f_vec(f_mask), d_vec, dR_vec);
toc

%% plot forward model results

figure(2); clf; hold on;
set(gcf, 'WindowState', 'maximized');
plotAll(Z_NONE_AVG(1,f_mask,mask), Z_GAPS(:,f_mask,mask), Z_PADS(:,f_mask,mask), f_vec(f_mask), '.');
for i = 1:8, subplot(2,4,i); set(gca, ColorOrderIndex=1); end
plotAll(Z_NONE_sim(1,f_mask,:), Z_GAPS_sim(:,f_mask,:), Z_PADS_sim(:,f_mask,:), f_vec(f_mask), 'o-');
shg;

%% try fitting

l_chan_init = [40, 10, 40, 20, 40, 40, 40, 80, 40, 160, 40, 320, 40, 640, 40] .* 1e-6;

n_diams = 5;
d_vec = linspace(0, (21.5e-6)^3, n_diams).^(1/3);
d_vec(1) = [];
dR_vec = interp1(d_vec_LUT, dR_vec_LUT, d_vec, "linear", "extrap"); % interp from LUT to get dR_vec for the d_vec above

pad_mult_init = 1;
pad_add_init = 0e-6;
w_init = 20.0e-6;
h_init = 21.5e-6;
rho_init = 0.625;
% q_vec_init = [0.0362    0.9986    0.1971    0.1051    0.0933    0.1969    0.0616    0.0727    0.0362];
% a_vec_init = [0.6808    0.9475    1.0000    0.9575    0.9222    1.0000    0.8873    0.8671    0.6808];
q_vec_init = [1, 0.1, 0.1, 0.133, 0.0665, 0.0128, 0.0147, 0.0304, 1];
a_vec_init = [1, 1.0, 1.0, 0.887, 0.900, 0.999, 0.962, 0.896, 1];
% q_vec_init(:) = 0.2;
% a_vec_init(:) = 1;
dR_mult_init = 1;
params_init = [pad_mult_init, pad_add_init, w_init, h_init, rho_init, q_vec_init, a_vec_init, dR_mult_init];
lb(1) = 1; ub(1) = 1;
lb(2) = 0e-6; ub(2) = 0e-6;
lb(3) = 1*w_init; ub(3) = 1*w_init;
lb(4) = 1*h_init; ub(4) = 1*h_init;
lb(5) = 1*rho_init; ub(5) = 1*rho_init;
lb(6:14) = 0; ub(6:14) = Inf;
lb(15:23) = 0.75; ub(15:23) = 1;
lb(24) = 1*dR_mult_init; ub(24) = 1*dR_mult_init;

% params_fit = params_init;
max_iter = 500;

tic;
[Z_none_hat, Z_gaps_hat, Z_pads_hat, X_hat, Y_hat, ab, params_fit, rmse, mae, Linf, rmse_, mae_, Linf_] = ...
    fit_pads_XY(Z_NONE_AVG(:,f_mask,mask), Z_GAPS(:,f_mask,mask), Z_PADS(:,f_mask,mask), f_vec(f_mask), d_vec, dR_vec, params_fit, lb, ub, max_iter);
toc

%% plot scatter plots with fit on top

f_vec_more = logspace(log10(2.5e1),log10(250e1),25);
[Z_NONE_sim_more, Z_GAPS_sim_more, Z_PADS_sim_more, X_sim_more, Y_sim_more] = pad_array_vol_forward_model([1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 40, 640, 40] .* 1e-6, 20e-6, 20e-6, 0.625, q_vec*0+0.4, a_vec*0+1, f_vec_more, d_vec, dR_vec);

X_data = ( (abs(Z_GAPS(1:end-1,f_mask,mask))/2+abs(Z_GAPS(2:end,f_mask,mask))/2) ./ abs(Z_NONE_AVG(1,f_mask,mask)) - 1 ) * 100;
Y_data = ( abs(Z_PADS(:,f_mask,mask)) ./ abs(Z_NONE_AVG(1,f_mask,mask)) - 1 ) * 100;

figure(5); clf; hold on; set(gcf, 'WindowState', 'Maximized');
plotAll(Z_NONE_AVG(:,f_mask,mask), Z_GAPS(:,f_mask,mask), Z_PADS(:,f_mask,mask), f_vec(f_mask), '.');
for i = 1:8, subplot(2,4,i); set(gca, ColorOrderIndex=1); end
plotAll(Z_none_hat, Z_gaps_hat, Z_pads_hat, f_vec(f_mask), 'o');
for i = 1:8, subplot(2,4,i); set(gca, ColorOrderIndex=1); end
for i_p = 1:size(X_data,1)
    for i_f = 1:size(X_data,2)
        X_poly(i_p, i_f, :) = linspace(0, max([max(X_hat(i_p, i_f, :)), max(X_data(i_p, i_f, :))]), 101);
    end
end
Y_poly = ab(:,:,1) .* X_poly.^2 + ab(:,:,2) .* X_poly;
plotAllXY(X_poly, Y_poly, f_vec(f_mask), '-');
for i = 1:8, subplot(2,4,i); set(gca, ColorOrderIndex=1); end
plotAllXY(X_sim_more, Y_sim_more, f_vec_more, 's--');
str = sprintf('w = %0.2f um, h = %0.2f um, rho = %0.4f Ohm*m, dR_mult = %0.4f\nq_vec = [ %s ] Ss^a\na_vec = [ %s ]\nN = %d, RMSE = %0.4f (%0.2f%%), MAE = %0.4f (%0.2f%%), Linf = %0.4f (%0.2f%%)', ...
               params_fit(3)*1e6, params_fit(4)*1e6, params_fit(5), params_fit(24), ...
               num2str(round(params_fit(7:13),4)), num2str(round(params_fit(16:22),4)), sum(mask), rmse, rmse_*100, mae, mae_*100, Linf, Linf_*100);
subplot(2,4,8); title(str, 'interpreter', 'none');
axis([0,3.5,-1,3.5]);
shg;

for ii = 1:8
    subplot(2,4,ii);
    set(gca, XAxisLocation="origin");
    set(gca, XTick=[0,1,2,3]); set(gca, YTick=[-1,0,1,2,3]);
end

savename = 'output/pads_fit_qaonly_18p79padarr_dev001b_cells_L1.fig';
saveas(gcf, [savename, '.fig']);
saveas(gcf, [savename, '.png']);
saveas(gcf, [savename, '.pdf']);
save([savename, '.mat'], 'Z_none_hat', 'Z_gaps_hat', 'Z_pads_hat', 'X_hat', 'Y_hat', 'ab', 'params_fit', 'rmse', 'mae', 'Linf', 'rmse_', 'mae_', 'Linf_', 'mask', 'f_mask', 'd_vec', 'params_init', 'lb', 'ub');

close all;
clearvars;
clc;

%% fit Cs and dc to match pad-to-gap ratio

model = mphload('../comsol_models/GPG_ellipsoid_no_electrodes.mph');

wc = 10e-6;
hc = 18e-6;
Cs_init = 0.20;
dc_init = 12e-6;
params_init = [Cs_init, dc_init];
target = [5.17, 5.26, 5.19, 3.42, 4.41, -0.71, 1.52, -0.02];
target = target(2:2:end) ./ target(1:2:end-1);
target = [0.052, target]; % add gap response 1 as target

weights = [4,1,1,1,1]; weights = weights ./ sum(weights);
fun = @(x) obj_func(x, model, target, wc, hc) .* weights;
lb = [Cs_init*0.9, dc_init*0.9];
ub = [Cs_init*1.1, dc_init*1.1];
[params_fit, res_norm, res, exit_flag, out, lam, jac] = lsqnonlin(fun, params_init, lb, ub);

[err, rmse_init, mae_init, Linf_init, target_hat] = obj_func(params_fit, model, target, wc, hc) % double check the response values

%% simulate transit

model = mphload('../comsol_models/GPG_ellipsoid_no_electrodes.mph');

x_vec = (-110e-6 : 1e-6 : 0e-6).';
% x_vec = [-110e-6, -70e-6, 0e-6].';
% Cs = params_fit(1);
Cs = 0.18; % [F/m^2]
% dc = params_fit(2);
dc = 11.5e-6; % [um]
dc_ = dc;

model.param.set('Cs', sprintf('%0.9f [F/m^2]', Cs));
model.param.set('dc', sprintf('%0.9f [um]', dc_*1e6));

f_vec = [1e3, 10e3, 50e3, 400e3];
Z_mat = nan(length(x_vec), length(f_vec));

tic;
count = 0;
for i = 1:length(x_vec)
    xc = x_vec(i);
    model.param.set('xc', sprintf('%0.9f [um]', xc*1e6));
    
    succeeded = false;
    while ~succeeded
        try
            model.study('std1').run;
            fprintf('Simulation succeeded!\n')
            succeeded = true;
        catch
            fprintf('Simulation failed: trying again with shrunken cell!\n');
            disp(mphgetexpressions(model.param));
            dc_ = dc_ * 0.99999;
            [ac_, bc_, cc_] = find_abc(100e-6, dc_, 0, 200e-6, wc, hc);
            model.param.set('dc', sprintf('%0.9f [um]', dc_*1e6));
            model.param.set('ac', sprintf('%0.9f [um]', ac_*1e6));
            model.param.set('bc', sprintf('%0.9f [um]', bc_*1e6));
            model.param.set('cc', sprintf('%0.9f [um]', cc_*1e6));
        end
    end
    Z = mphglobal(model, '(ec.V0_1-ec.V0_2)/ec.I0_1');
    Z_mat(i,:) = Z;

    count = count + 1;
    fprintf('Time elapsed: %0.0f min\nTime remaining: %0.0f min\n', toc/60, toc/60/count*length(x_vec)-toc/60);
end

save('simulated_transit_and_params_fit.mat', 'params_fit', 'Cs', 'dc', 'x_vec', 'f_vec', 'Z_mat');

%% plot impedance

x_vec_ = (-250:250).' .* 1e-6;
Z_mat_ = [zeros(140,4) + Z_mat(1,:); Z_mat; flipud(Z_mat(1:end-1,:)); zeros(140,4) + Z_mat(1,:)];
Z_baseline_ = zeros(size(x_vec_)) + Z_mat_(1,:);
rangeZ = range(abs(Z_mat_)/1e3, 1);
xmin = -250;
xmax = 250;
ymin = min(abs(Z_mat_)/1e3, [], 1) - rangeZ*0.1;
ymax = max(abs(Z_mat_)/1e3, [], 1) + rangeZ*0.1;

fig1 = figure; clf;
for i = 1:4
    ax1(i) = subplot(4,1,i); hold on;
    plot(x_vec_*1e6, abs(Z_baseline_(:,i))/1e3, 'k-');
    set(gca, ColorOrderIndex=1);
    plot(x_vec_*1e6, abs(Z_mat_(:,i))/1e3);
    axis([xmin, xmax, ymin(i), ymax(i)]);
    ylabel('|Z| [k\Omega]');
    legend(f_vec(i)/1e3 + "kHz");
end
xlabel('Horizontal cell position [um]');
linkaxes(ax1, 'x');
shg;

savename = 'output/simulated_GPG_transit';
saveas(gcf, [savename '.fig']);
saveas(gcf, [savename '.png']);
saveas(gcf, [savename '.pdf']);

%% plot normalized impedance

Z_norm_ = abs(Z_mat_)./abs(Z_baseline_) - 1;
rangeZ = range(abs(Z_norm_)*100, 1);
xmin = -250;
xmax = 250;
ymin = min(Z_norm_*100, [], 1) - rangeZ*0.1;
ymax = max(Z_norm_*100, [], 1) + rangeZ*0.1;

fig2 = figure; clf;
for i = 1:4
    ax2(i) = subplot(4,1,i); hold on;
    plot(x_vec_*1e6, Z_norm_(:,i)*100);
    axis([xmin, xmax, ymin(i), ymax(i)]);
    ylabel('\Delta|Z|/|Z| (%)');
    legend(f_vec(i)/1e3 + "kHz");
end
xlabel('Horizontal cell position [um]');
linkaxes(ax2, 'x');
shg;

savename = 'output/simulated_GPG_transit_norm';
saveas(gcf, [savename '.fig']);
saveas(gcf, [savename '.png']);
saveas(gcf, [savename '.pdf']);

%% functions

function [err, rmse, mae, Linf, target_hat] = obj_func(params, model, target, wc, hc) % target is [pr1/gr1, pr2/gr2, pr3/gr3, pr4/gr4]

Cs = params(1);
dc = params(2);

[Z_none, Z_gap, Z_pad] = simulate_GPG(model, Cs, dc, wc, hc);
gr1 = (abs(Z_gap(1)) - abs(Z_none(1))) ./ abs(Z_none(1));
pr1 = (abs(Z_pad(1)) - abs(Z_none(1))) ./ abs(Z_none(1));
gr2 = (abs(Z_gap(2)) - abs(Z_none(2))) ./ abs(Z_none(2));
pr2 = (abs(Z_pad(2)) - abs(Z_none(2))) ./ abs(Z_none(2));
gr3 = (abs(Z_gap(3)) - abs(Z_none(3))) ./ abs(Z_none(3));
pr3 = (abs(Z_pad(3)) - abs(Z_none(3))) ./ abs(Z_none(3));
gr4 = (abs(Z_gap(4)) - abs(Z_none(4))) ./ abs(Z_none(4));
pr4 = (abs(Z_pad(4)) - abs(Z_none(4))) ./ abs(Z_none(4));
target_hat = [gr1, pr1/gr1, pr2/gr2, pr3/gr3, pr4/gr4];
err = (target - target_hat) ./ [0.052 / 4, 1, 0.66, 0.16, 0.013];
rmse = sqrt(mean(err.^2));
mae = mean(abs(err));
Linf = max(abs(err));
fprintf('Cs, dc = \n');
disp(params .* [1, 1e6]);
fprintf('target_hat = \n');
disp(target_hat);
fprintf('target = \n');
disp(target);
fprintf('rmse = %0.4f, mae = %0.4f, Linf = %0.4f\n', rmse, mae, Linf);

end

function [Z_none, Z_gap, Z_pad] = simulate_GPG(model, Cs, dc, wc, hc)

dc_ = dc;
[ac_, bc_, cc_] = find_abc(100e-6, dc_, 0, 200e-6, wc, hc);

% set params
model.param.set('Cs', sprintf('%0.9f [F/m^2]', Cs));
model.param.set('dc', sprintf('%0.9f [um]', dc_*1e6));
model.param.set('ac', sprintf('%0.9f [um]', ac_*1e6));
model.param.set('bc', sprintf('%0.9f [um]', bc_*1e6));
model.param.set('cc', sprintf('%0.9f [um]', cc_*1e6));

% cell outside
model.param.set('xc', 'xout');
succeeded = false;
while ~succeeded
    try
        model.study('std1').run;
        fprintf('Simulation succeeded!\n')
        succeeded = true;
    catch
        fprintf('Simulation failed: trying again with shrunken cell!\n');
        disp(mphgetexpressions(model.param));
        dc_ = dc_ * 0.99999;
        [ac_, bc_, cc_] = find_abc(100e-6, dc_, 0, 200e-6, wc, hc);
        model.param.set('dc', sprintf('%0.9f [um]', dc_*1e6));
        model.param.set('ac', sprintf('%0.9f [um]', ac_*1e6));
        model.param.set('bc', sprintf('%0.9f [um]', bc_*1e6));
        model.param.set('cc', sprintf('%0.9f [um]', cc_*1e6));
    end
end
Z_none = mphglobal(model, '(ec.V0_1-ec.V0_2)/ec.I0_1');

% cell in gap
model.param.set('xc', 'xgap');
succeeded = false;
while ~succeeded
    try
        model.study('std1').run;
        fprintf('Simulation succeeded!\n')
        succeeded = true;
    catch
        fprintf('Simulation failed: trying again with shrunken cell!\n');
        disp(mphgetexpressions(model.param));
        dc_ = dc_ * 0.99999;
        [ac_, bc_, cc_] = find_abc(100e-6, dc_, 0, 200e-6, wc, hc);
        model.param.set('dc', sprintf('%0.9f [um]', dc_*1e6));
        model.param.set('ac', sprintf('%0.9f [um]', ac_*1e6));
        model.param.set('bc', sprintf('%0.9f [um]', bc_*1e6));
        model.param.set('cc', sprintf('%0.9f [um]', cc_*1e6));
    end
end
Z_gap = mphglobal(model, '(ec.V0_1-ec.V0_2)/ec.I0_1');

% cell in pad
model.param.set('xc', 'xpad');
succeeded = false;
while ~succeeded
    try
        model.study('std1').run;
        fprintf('Simulation succeeded!\n')
        succeeded = true;
    catch
        fprintf('Simulation failed: trying again with shrunken cell!\n');
        disp(mphgetexpressions(model.param));
        dc_ = dc_ * 0.99999;
        [ac_, bc_, cc_] = find_abc(100e-6, dc_, 0, 200e-6, wc, hc);
        model.param.set('dc', sprintf('%0.9f [um]', dc_*1e6));
        model.param.set('ac', sprintf('%0.9f [um]', ac_*1e6));
        model.param.set('bc', sprintf('%0.9f [um]', bc_*1e6));
        model.param.set('cc', sprintf('%0.9f [um]', cc_*1e6));
    end
end
Z_pad = mphglobal(model, '(ec.V0_1-ec.V0_2)/ec.I0_1');

end

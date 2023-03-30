close all;
clearvars;
clc;

%% load model

model = mphload('../comsol_models/GPG_ellipsoid_no_electrodes.mph');

%% DOEs

wc = 10e-6; % channel width [m]
hc = 20e-6; % channne height [m]
lg = 40e-6; % gap length [m]
lp = 80e-6; % pad length [m]
Cs = 0.2; % specific capacitance [F/m^2]
does = struct([]);
d_vec = linspace(0, lg*wc*hc, 41).^(1/3); d_vec(1) = [];
d_vec = round(d_vec, 10);
for i = 1:length(d_vec)
    does(i).name = sprintf('%0.4f um free diameter ellipsoid in 10 x 20 cross section', d_vec(i)*1e6);
end

%% iterate over DOEs

tic;
count = 0;
for i = 1:length(does)

    %% change params

    model.param.set('lp', sprintf('%0.4f [um]', lp*1e6));
    model.param.set('lg', sprintf('%0.4f [um]', lg*1e6));
    model.param.set('wc', sprintf('%0.4f [um]', wc*1e6));
    model.param.set('hc', sprintf('%0.4f [um]', hc*1e6));
    model.param.set('Cs', sprintf('%0.4f [F/m^2]', Cs));
    
    dc_ = d_vec(i);
    [ac_, bc_, cc_] = find_abc(100e-6, dc_, 0, 200e-6, wc, hc);
    model.param.set('dc', sprintf('%0.9f [um]', dc_*1e6));
    model.param.set('ac', sprintf('%0.9f [um]', ac_*1e6));
    model.param.set('bc', sprintf('%0.9f [um]', bc_*1e6));
    model.param.set('cc', sprintf('%0.9f [um]', cc_*1e6));
    
    params = mphgetexpressions(model.param);
    does(i).params = params;
    does(i).diam_um = d_vec(i)*1e6;
    does(i).shape = 'ellipsoid';
    does(i).vol_pL = 4/3*pi*((does(i).diam_um)*1e-6/2)^3*1e12;
    
    freq = [0.001; logspace(1,6,501).'; 1.885e6]; % make sure .mph file has the same frequency list
    does(i).freq = freq;
    
    %% run study to compute Zpad
    
    dc_ = d_vec(i);
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
    Z = mphglobal(model, '(ec.V0_1-ec.V0_2)/ec.I0_1');
    does(i).Zpad = Z;
    
    %% repeat study for Zgap
    
    dc_ = d_vec(i);
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
    Z = mphglobal(model, '(ec.V0_1-ec.V0_2)/ec.I0_1');
    does(i).Zgap = Z;
    
    count = count + 1;
    fprintf('Time elapsed: %0.0f min\nTime remaining: %0.0f min\n', toc/60, toc/60/count*length(does)-toc/60);
    
end

%% repeat study for Znone

model.param.set('xc', '-lc/2 - 0.02 [um]');
model.param.set('ac', '0.02 [um]'); % make
model.param.set('bc', '0.02 [um]'); % cell
model.param.set('cc', '0.02 [um]'); % small
model.study('std1').run;
Z = mphglobal(model, '(ec.V0_1-ec.V0_2)/ec.I0_1');
[does(:).Znone] = deal(Z);

%% compute ratios

for i = 1:length(does)
    does(i).ratio_pad_none = abs(does(i).Zpad) ./ abs(does(i).Znone);
    does(i).ratio_gap_none = abs(does(i).Zgap) ./ abs(does(i).Znone);
    does(i).ratio_dip_peak = (does(i).ratio_pad_none-1) ./ (does(i).ratio_gap_none-1);
end

%% plot simulated dR/R vs d/D on log-log plot, compare to Deblois-Bean+Smythe and Maxwell

rho = 0.625;
D_eff = sqrt(wc*hc/pi)*2; % guess at effective diameter
L = lg + lp + lg; % length of channel
R_DC = rho * L / (wc * hc); % analytical DC resistance of channel
for i = 1:length(does)
    dR_vec(i) = real(does(i).Zgap(1))-real(does(i).Znone(1));
    d = does(i).diam_um*1e-6;
    if d <= min([wc,hc])
        dRoR_vec_smythe(i) = d^3/(D_eff^2*L)/(1-0.8*(d/D_eff)^3);
    else
        [a_cell, b_cell, c_cell] = find_abc(L/2, d, 0, L, wc, hc);
        d_eff = sqrt(b_cell*c_cell);
        dRoR_vec_smythe(i) = d_eff^3/(D_eff^2*L)/(1-0.8*(d_eff/D_eff)^3) * a_cell/d;
    end
    dRoR_vec_vol(i) = d^3*pi*rho/(6*wc^2*hc^2)/R_DC;
end

dRoR_vec_sim = dR_vec/R_DC;
doD_vec = d_vec/D_eff;

figure; clf; hold on;
plot(doD_vec, dRoR_vec_vol, DisplayName="Maxwell");
plot(doD_vec, dRoR_vec_smythe, DisplayName="Deblois-Bean+Smythe");
plot(doD_vec, dRoR_vec_sim, DisplayName="COMSOL");
set(gca, xscale="log", yscale="log");
xlabel("d/D");
ylabel("\DeltaR/R")
legend(location="northwest");
title('Ellipsoid in long rectangular channel');
grid on; grid minor;
shg;

figure; clf; hold on;
dVoV = 4/3*pi*(d_vec/2).^3./(L*wc*hc);
plot(dVoV, dRoR_vec_vol, DisplayName="Maxwell");
plot(dVoV, dRoR_vec_smythe, DisplayName="Deblois-Bean+Smythe");
plot(dVoV, dRoR_vec_sim, DisplayName="COMSOL");
set(gca, xscale="linear", yscale="linear");
xlabel("\DeltaV/V");
ylabel("\DeltaR/R")
legend(location="northwest");
title('Ellipsoid in long rectangular channel');
grid on; grid minor;
shg;

%% compute gap reponse and pad response values and plot diagram

for i = 1:length(does)
    X_sim_dR(i,:) = ((real(does(i).Zgap(1))-real(does(i).Znone(1)))./abs(does(i).Znone(:)))*100; % dR/|Z|
    X_sim(i,:) = (abs(does(i).Zgap(:))./abs(does(i).Znone(:))-1)*100;
    Y_sim(i,:) = (abs(does(i).Zpad(:))./abs(does(i).Znone(:))-1)*100;
end
figure; clf; hold on;
inds_plot = 1:5:length(does(i).freq);
plot(X_sim(:,inds_plot), Y_sim(:,inds_plot), '.-');
% plot(X_sim_dR(:,inds_plot), Y_sim(:,inds_plot), '.--');
axis equal;

%% save results

save('DOEs_FEA_ellipsoid_G40P80G40_W10xH20.mat', 'does', 'rho', 'wc', 'hc', 'lg', 'lp', 'D_eff', 'L', 'R_DC', 'dRoR_vec_sim', 'dRoR_vec_smythe', 'dRoR_vec_vol', 'dVoV', 'd_vec', 'dR_vec', 'doD_vec', 'X_sim_dR', 'X_sim', 'Y_sim');

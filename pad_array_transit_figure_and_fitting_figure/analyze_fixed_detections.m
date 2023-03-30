% extract the gap and pad response values

close all;
clearvars;
dbstop if error;

%%

p = '../fixed_detections/';
out_p = '../analyzed_detections/';
d = dir(p);
fns = {d.name};
mask = contains(fns, '.mat');
fns = fns(mask);

plotting = [false, false, true];
% plotting = [false, false, false];

%%

for i = length(fns)
    fn = fns{i};
    fprintf('Loading data (%d of %d) from %s\n', i, length(fns), fn);
    load([p fn]);

    %%
    
    for j = 1:length(good_detections)
        
        fprintf('Analyzing detection %d of %d\n', j, length(good_detections));

        gd = good_detections(j);
        Z_dc = gd.Z_dc;
        t_vec = gd.t_vec;
        t_center = gd.edge_center_time;
        t_start = gd.edge_start_time;
        t_stop = gd.edge_stop_time;
        t_transit = gd.edge_transit_time;
        t_edge = 0.25 * t_transit;
        sig_inds = t_vec >= t_start - t_edge & t_vec <= t_stop + t_edge;
        Z_sig = Z_dc(sig_inds,:);
        t_vec = t_vec(sig_inds);
        t_vec_norm = (t_vec - t_center) / t_transit;

        [Z_w, level, wname, npc] = wavelet_pca_denoise(Z_sig);
        
        tt = gd.edge_transit_time;
        len_mean = max([floor(tt * 10 / 1590 * 0.5 * fr), 1]);
        if mod(len_mean, 2) == 0, len_mean = len_mean - 1; end
        len_med = max([floor(tt * 10 / 1590 * 1 * fr), 1]);
        if mod(len_med, 2) == 0, len_med = len_med - 1; end
        Z_wm = movmean(Z_w, len_mean, 1);
        Z_wmm = movmedian(Z_wm, len_med, 1);
        
        f_vec = gd.freq;
        for k = 1:length(f_vec)
            b_vec(k) = robust_median(Z_dc(:,k));
        end
        
        if plotting(1)
            figure(1); clf;
            leg1 = cell(0);
            plotAll(t_vec_norm, abs(Z_sig), f_vec); leg1{end+1} = 'comb filtered';
            plotAll(t_vec_norm, abs(Z_w), f_vec); leg1{end+1} = '+ WMSPCA';
            plotAll(t_vec_norm, abs(Z_wm), f_vec); leg1{end+1} = '+ movmean';
            plotAll(t_vec_norm, abs(Z_wmm), f_vec); leg1{end+1} = '+ movmedian';
            legend(leg1);
            shg;
        end

        if plotting(2)
            figure(2); clf;
            leg2 = cell(0);
            plotAllRes(t_vec_norm, abs(Z_w-Z_sig), f_vec); leg2{end+1} = sprintf('after WMSPCA: J = %0.2f nats', negentropy(abs(Z_w-Z_sig)));
            plotAllRes(t_vec_norm, abs(Z_wm-Z_sig), f_vec); leg2{end+1} = sprintf('after movmean: J = %0.2f nats', negentropy(abs(Z_wm-Z_sig)));
            plotAllRes(t_vec_norm, abs(Z_wmm-Z_sig), f_vec); leg2{end+1} = sprintf('after movmedian: J = %0.2f nats', negentropy(abs(Z_wmm-Z_sig)));
            legend(leg2);
            shg;
        end

        % DO FITTING

%         z_vec = abs(Z_wmm(:,end-1) - RM_vec(end-1)); % use the magnitude of the good freq for fitting
        f_inds = 1:6;
        K = length(f_inds)*2; % number of channels
        z_mat = [real(Z_wmm(:,f_inds)), imag(Z_wmm(:,f_inds))]; % stack all channels of Z
        bl_vec = [real(b_vec(f_inds)), imag(b_vec(f_inds))]; % stack all channels of robust median
        z_mat = (z_mat - bl_vec) ./ repmat(abs(b_vec(f_inds)), [1,2]);
        t_vec = gd.t_vec(sig_inds);
        tc_init = gd.edge_center_time;
        tt_init = gd.edge_transit_time * 0.975;
        fext_init = 0.1;
        fmis_init = 0;
        facc_init = 1;
        amps_init = range(z_mat) .* sign(median(z_mat));
        w_init = 25e-6;
        h_init = 30e-6;
        dR_init = max(real(Z_wmm(t_vec_norm > -0.5 & t_vec_norm < 0.5, 1))) - real(b_vec(1));
        rho_init = 0.625;
        cpad_init = 0.04;
%         cpad_init = 5.53e-4;
        params_init = [tc_init, tt_init, fext_init, fmis_init, facc_init, amps_init, w_init, h_init, dR_init, rho_init, cpad_init];
        lb = nan(size(params_init)); ub = nan(size(params_init));
        lb(1) = tc_init - 0.1 * tt_init; ub(1) = tc_init + 0.1 * tt_init;
        lb(2) = tt_init * 0.9; ub(2) = tt_init * 1.1;
%         lb(3) = 0; ub(3) = 2 * tr_init;
        lb(3) = 0; ub(3) = Inf;
%         lb(4) = -0.1; ub(4) = 0.1;
        lb(4) = 0; ub(4) = 0;
        lb(5) = 0.5; ub(5) = 2;
%         lb(5) = 1; ub(5) = 1;
        lb(6:6+K-1) = -Inf; ub(6:6+K-1) = Inf;
        lb(6+K) = 0.5 * w_init; ub(6+K) = 2 * w_init;
        lb(7+K) = 0.5 * h_init; ub(7+K) = 2 * h_init;
        lb(8+K) = 1 * dR_init; ub(8+K) = 1 * dR_init;
%         lb(9+K) = 0.5 * rho_init; ub(9+K) = 2 * rho_init;
        lb(9+K) = rho_init; ub(9+K) = rho_init;
        lb(10+K) = 0.01 * cpad_init; ub(10+K) = 100 * cpad_init;
%         lb(10+K) = cpad_init; ub(10+K) = cpad_init;

        [z_fit, edge_times_fit, params_fit, rmse, mae, Linf, rmse_, mae_, Linf_] = fit_edge_times(t_vec, z_mat, f_vec(f_inds), params_init, lb, ub);
        Z_fit = z_fit .* repmat(abs(b_vec(f_inds)), [1,2]) + bl_vec;
        Z_fit = Z_fit(:,1:size(Z_fit,2)/2) + 1i .* Z_fit(:,size(Z_fit,2)/2+1:end);
        %%
        if plotting(3)
            figure(3); clf;
            plotAll(t_vec_norm, abs(Z_wmm), f_vec);
            for ii = 1:length(f_inds)
                fi = f_inds(ii);
                figure(3); subplot(2,4,fi); hold on;
                plot(t_vec_norm, abs(Z_fit(:,ii))/1e3, 'k-', DisplayName="fit");
            end
            shg;
        end
        %%

        good_detections(j).Z_sig = Z_sig;
        good_detections(j).Z_w = Z_w;
        good_detections(j).Z_wm = Z_wm;
        good_detections(j).Z_wmm = Z_wmm;
        good_detections(j).b_vec = b_vec;
        good_detections(j).t_sig = t_vec;

        denoise_opt.level = level;
        denoise_opt.wname = wname;
        denoise_opt.npc = npc;
        denoise_opt.len_mean = len_mean;
        denoise_opt.len_med = len_med;
        good_detections(j).denoise_opt = denoise_opt;
        
        fitting.params_init = params_init;
        fitting.lb = lb;
        fitting.ub = ub;
        fitting.z_fit = z_fit;
        fitting.Z_fit = Z_fit;
        fitting.resid = z_fit - z_mat;
        fitting.params_fit = params_fit;
        fitting.rmse = rmse;
        fitting.mae = mae;
        fitting.Linf = Linf;
        fitting.rmse = rmse_;
        fitting.mae = mae_;
        fitting.Linf = Linf_;
        good_detections(j).fitting = fitting;

        pf = params_fit;
        K = length(pf) - 10;
        pfs.tc = pf(1);
        pfs.tt = pf(2);
        pfs.fext = pf(3);
        pfs.fmis = pf(4);
        pfs.facc = pf(5);
        pfs.amps = pf(6:6+K-1);
        pfs.w = pf(6+K);
        pfs.h = pf(7+K);
        pfs.dR = pf(8+K);
        pfs.rho = pf(9+K);
        pfs.cpad = pf(10+K);
        good_detections(j).params_fit_struct = pfs;
        
        good_detections(j).edge_times_fit = edge_times_fit;

    end

end

%% sample the gap and pad impedances

for i = 1:length(fns)
    fn = fns{i};
    load([out_p, fn]);
    
    for j = 1:length(good_detections)
        gd = good_detections(j);
        edge_times = gd.edge_times_fit;
        dur_pads = edge_times(3:2:end)-edge_times(2:2:end-1);
        n_pads = length(dur_pads);
        min_pad_dur = min(dur_pads);
        Z_wmm = gd.Z_wmm;
        b_vec = gd.b_vec;
        Z_sig = gd.Z_sig;
        t_sig = gd.t_sig;
        t_sig_norm = (t_sig - gd.edge_center_time) / gd.edge_transit_time;
        fr = round(1/(t_sig(2) - t_sig(1)));
        t_sig_dif = (t_sig(2:end) + t_sig(1:end-1)) / 2;
        t_sig_dif_norm = (t_sig_dif - gd.edge_center_time) / gd.edge_transit_time;

        % find falling edge and rising edge of each pad
        fine_edge_times = [];
        for i_pad = 1:n_pads
            dur_pad_i = dur_pads(i_pad);    
            len_filt_i = max([1, floor(2 * dur_pad_i * fr)]);
            if mod(len_filt_i, 2) == 0, len_filt_i = len_filt_i + 1; end
            Z_filt_i = movmean(Z_wmm, len_filt_i, 1);
            dif_sig_i = abs(diff(Z_filt_i(:,end-1) - b_vec(end-1))) .* sign(diff(abs(Z_filt_i(:,end-1) - b_vec(end-1))));
            [pn_i, tn_i, wn_i, prn_i] = findpeaks(-dif_sig_i, t_sig_dif);
            [pp_i, tp_i, wp_i, prp_i] = findpeaks(dif_sig_i, t_sig_dif);
            n_peaks = min([15, length(pn_i)]);
            [~, order_n] = sort(prn_i, "descend");
            pn_i = -pn_i(order_n(1:n_peaks));
            tn_i = tn_i(order_n(1:n_peaks));
            n_peaks = min([15, length(pp_i)]);
            [~, order_p] = sort(prp_i, "descend");
            pp_i = pp_i(order_p(1:n_peaks));
            tp_i = tp_i(order_p(1:n_peaks));
            [~, i_tf] = min(abs(tn_i-edge_times(i_pad*2)));
            tf_i = tn_i(i_tf); % time of falling edge for ith pad
            [~, i_tr] = min(abs(tp_i-edge_times(i_pad*2+1)));
            tr_i = tp_i(i_tr); % time of falling edge for ith pad
            fine_edge_times = [fine_edge_times, tf_i, tr_i];
        end
        fine_edge_times_norm = (fine_edge_times - gd.edge_center_time) / gd.edge_transit_time;  

        % sample the gap impedance points
        dur_gaps = edge_times(2:2:end)-edge_times(1:2:end-1);
        t_gaps = (edge_times(2:2:end) + edge_times(1:2:end-1)) / 2; % sample times for gap centers
        t_gaps(1) = edge_times(2) - dur_gaps(2)/2; % account for stalling at the start
        n_gaps = length(dur_gaps);
        Z_gaps = [];
        for i_gap = 1:n_gaps
            dur_gap_i = dur_gaps(i_gap);
            len_filt_i = max([1, floor(dur_gap_i / 3 * fr)]);
            if mod(len_filt_i, 2) == 0, len_filt_i = len_filt_i - 1; end
            Z_filt_i = movmean(Z_sig, len_filt_i, 1);
            [~, i_t_gap] = min(abs(t_sig - t_gaps(i_gap)));
            Z_gap_i = Z_filt_i(i_t_gap,:);
            Z_gaps = [Z_gaps; Z_gap_i];
        end

        % sample the pad impedance points
        t_pads = (edge_times(3:2:end) + edge_times(2:2:end-1)) / 2; % sample times for pad centers
        Z_pads = [];
        for i_pad = 1:n_pads
            dur_gap_i = dur_gaps(i_pad+1); % use gap length to filter pad response
            len_filt_i = max([1, floor(dur_gap_i / 3 * fr)]);
            if mod(len_filt_i, 2) == 0, len_filt_i = len_filt_i - 1; end
            Z_filt_i = movmean(Z_sig, len_filt_i, 1);
            [~, i_t_pad] = min(abs(t_sig - t_pads(i_pad)));
            Z_pad_i = Z_filt_i(i_t_pad,:);
            Z_pads = [Z_pads; Z_pad_i];
        end

        good_detections(j).fine_edge_times = fine_edge_times;
        good_detections(j).Z_gaps = Z_gaps;
        good_detections(j).Z_pads = Z_pads;
        good_detections(j).t_gaps = t_gaps;
        good_detections(j).t_pads = t_pads;

    end

    save([out_p, fn], 'good_detections');

end

%% functions

function plotAll(t_vec, Z_plot, f_vec)
    yr = max(range(Z_plot/1e3))*1.2;
    ax(8) = subplot(2,4,8); hold on;
    Z_mean = mean(Z_plot,2);
    RM_mean = robust_median(Z_mean);
%     plot(t_vec, 100*(Z_mean - RM_mean)/RM_mean, LineWidth=1);
    plot(t_vec, Z_mean/1e3);
    xlabel('Normalized time');
%     ylabel('( |Z| - |Z_0| ) / |Z_0| [%]');
    ylabel('|Z| [k\Omega');
    title("Average of all frequencies");
    set(gca, xlim=[-0.75,0.75]);
    axi = axis;
    ym = (axi(3)+axi(4))/2;
    set(gca, ylim=[ym-yr/2, ym+yr/2]);
    for i = 1:7
        ax(i) = subplot(2,4,i); hold on;
        Z_curr = Z_plot(:,i);
        RM_curr = robust_median(Z_curr);
%         plot(t_vec, 100*(Z_curr - RM_curr)/RM_curr);
        plot(t_vec, Z_curr/1e3);
        xlabel('Normalized time');
%         ylabel('( |Z| - |Z_0| ) / |Z_0| [%]');
        ylabel('|Z| [k\Omega');
        title(f_vec(i)/1e3 + " kHz");
        set(gca, xlim=[-0.75,0.75]);
        axi = axis;
        ym = (axi(3)+axi(4))/2;
        set(gca, ylim=[ym-yr/2, ym+yr/2]);
    end
    linkaxes(ax, 'x');
end

function plotAllRes(t_vec, res, f_vec)
    ax(8) = subplot(2,4,8); hold on;
    res_mean = mean(res,2);
    plot(t_vec, res_mean, LineWidth=1);
    xlabel('Normalized time');
    ylabel('Residual [\Omega]');
    title("Average of all frequencies");
    set(gca, xlim=[-0.75,0.75]);
    for i = 1:7
        ax(i) = subplot(2,4,i); hold on;
        res_curr = res(:,i);
        plot(t_vec, res_curr, LineWidth=1);
        xlabel('Normalized time');
        ylabel('Residual [\Omega]');
        title(f_vec(i)/1e3 + " kHz");
        set(gca, xlim=[-0.75,0.75]);
    end
    linkaxes(ax, 'x');
end

%%

function [rmse, mae, Linf, J] = foms(x, x_orig)
    rmse = RMSE(x, x_orig);
    mae = MAE(x, x_orig);
    Linf = max(abs(x-x_orig), [], "all");
    J = negentropy(x - x_orig);
end

function rmse = RMSE(x, x_orig)
    rmse = sqrt(mean((x - x_orig).^2, "all")) / sqrt(mean((x_orig).^2, "all"));
end

function mae = MAE(x, x_orig)
    mae = mean(abs(x - x_orig), "all") / mean(abs(x_orig), "all");
end

function J = negentropy(x)
    Hx = differential_entropy(x, [], 'KL');
    Sigma = cov(x);
    Hmvn = mvndifent(Sigma);
    J = Hmvn - Hx;
end

function H = mvndifent(Sigma)
    [n, n_] = size(Sigma);
    if n ~= n_, error('Sigma must be a square matrix!'); end
    H = n/2*log(2*pi) + log(det(Sigma))/2 + n/2;
end
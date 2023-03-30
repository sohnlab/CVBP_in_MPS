% new find padarr signals
close all;
clearvars;
clc;

dbstop if error;

%%

%%%
% plotting = true;
plotting = false;
%%%

curr_path = [pwd '/'];
in_path = 'bw1000Hz/';
out_path = 'detections/';
in_dir = dir(in_path);
fns = {in_dir.name};
pattern = {'run'};
mask = contains(fns, pattern);
fns = fns(mask);
fns = fliplr(fns);

bw_lpf = 250; % [Hz]
lambda = 1e10;

for i = 1:length(fns)
    fn = fns{i};
    if exist([out_path fn], 'file')
        fprintf('Skipping %s\n', fn);
        continue;
    else
        fprintf('Loading data from %s\n', fn);
        cd(in_path);
        load(fn);
        cd(curr_path);
    end
    f_inds_merge = freq_vec < 400e3;

    Z = Z_comp;
    fr_orig = fr;

    detections = struct([]);

    t_start = tr(1);
    t_stop = tr(end);
    Z = Z(tr >= t_start & tr <= t_stop, :);

    fr_temp = 10 * bw_lpf;
    m = mod(fr_temp, 500);
    fr = fr_temp - m + 500*(m>0); % new sampling rate
    
    Z_d = unpad_data(resample(pad_data(Z - Z(1,:), fr_orig), fr, fr_orig), fr) + Z(1,:); % resample
    
%     %%% reject huge outliers
%     mm = movmean(abs(mean(Z_d,2)), 1 * fr);
%     mult = 100;
%     max_iter = 5;
%     outlier_mask = false(size(mm));
%     for iter = 1:max_iter
%         med = median(mm(~outlier_mask));
%         outlier_mask = mm > mult * med;
%     end
%     if plotting
%         figure(99); clf; hold on;
%         plot(abs(Z_d));
%         plot(outlier_mask*range(abs(Z_d))+med);
%         set(gca,'yscale','log');
%     end
%     if sum(outlier_mask) > 0
%         fprintf('Removing outliers from data: %0.1f%% outliers\n', sum(outlier_mask)/length(outlier_mask));
%         Z_d = Z_d(~outlier_mask);
%     end
%     %%%
    
    tr = (0:length(Z_d)-1).' / fr; % time vector
    tr_comp = (0:length(Z_comp)-1).' / fr_orig;
    Z_dc = unpad_data(comb_filter_keep_DC(pad_data(abs(Z_d) - abs(Z_d(1,:)), fr), 500, fr), fr) + abs(Z_d(1,:)); % remove 60Hz and harmonics with comb filter

    Fs    = fr;
    N     = 10;       % Order
    Fc    = bw_lpf;      % Cutoff Frequency
    flag  = 'scale';  % Sampling Flag
    Alpha = 2.5;      % Window Parameter
    win = gausswin(N+1, Alpha);
    b  = fir1(N, Fc/(Fs/2), 'low', win, flag).';
    Z_dcl = unpad_data(filtfilt(b, 1, pad_data(Z_dc - Z_dc(1,:), fr)), fr) + Z_dc(1,:);
    Z_merge = abs(mean(Z_dcl(:,f_inds_merge), 2));

    sig_dur_vec = fliplr(logspace(log10(0.01),log10(10),751));
    sig_dur_vec(sig_dur_vec > 5 | sig_dur_vec < 0.05) = [];
    
    baseline_merge = [];
    baseline = nan(size(Z_dcl));

    for sig_dur = sig_dur_vec

        fprintf('Searching for signals with duration %0.1f ms\n', sig_dur*1000);
        gap_dur = sig_dur * 20/1590; % [s]
        pad_dur = sig_dur * 640/1590; % [s]
        len_movmax = round(pad_dur * 1.1 * fr);
        len_movmin = round(sig_dur * 1.1 * fr);
        len_movmed = round(gap_dur * 1.1 * fr);
        
        Z_merge_m = movmedian(Z_merge, len_movmed); % median filter
        ub = movmax(Z_merge_m, len_movmax); % compute upper bound of signal envelope
        lb = movmin(Z_merge_m, len_movmin); % compute lower bound of signal envelope
        env = ub - lb; % signal envelope
        env_thresh = auto_threshold(env, 0, Inf); % compute automatic threshold
        noise_mask = env < env_thresh; % mask of just noise
        Z_noise = Z_merge;
        Z_noise(~noise_mask) = NaN;
        for j = 1:size(Z_dcl,2)
            try
                baseline(:,j) = smooth_masked_fit(abs(Z_dcl(:,j)), noise_mask, lambda);
            catch
                if all(isnan(baseline(:,j)))
                    baseline(:,j) = median(abs(Z_dcl(:,j)));
                    fprintf('WARNING: baseline estimation failed! Using median as baseline\n');
                else
                    fprintf('WARNING: baseline estimation failed! Using last estimated baseline\n');
                end
            end
        end
        try
            baseline_merge = smooth_masked_fit(Z_merge, noise_mask, lambda);
        catch
            if isempty(baseline_merge)
                baseline_merge = zeros(size(Z_merge)) + median(Z_merge);
                fprintf('WARNING: baseline estimation failed! Using median as baseline\n');
            else
                fprintf('WARNING: baseline estimation failed! Using last estimated baseline\n');
            end
        end
        
        %% detection
        
        box_ca = cell(0);
        box_ca{1} = [zeros(1,20-1),0.5, ones(1,1), zeros(1,0)+1, ones(1,1), zeros(1,1)+1, ones(1,1), zeros(1,1)+1, ones(1,1), ...
                    zeros(1,2)+1, ones(1,1), zeros(1,4)+0.9, ones(1,1), 0.8,zeros(1,8-2)+0.6,0.8, ...
                    ones(1,1), 0.7667,0.5333,zeros(1,16-4)+0.3,0.5333,0.7667, ones(1,1), 0.5,zeros(1,20-1)].';
        box_ca{2} = [zeros(1,20-1),0.5, ones(1,1), zeros(1,0)+1, ones(1,1), zeros(1,1)+1, ones(1,1), zeros(1,1)+1, ones(1,1), ...
                    zeros(1,2)+0.9, ones(1,1), zeros(1,4)+0.6, ones(1,1), 0.65,zeros(1,8-2)+0.3,0.65, ...
                    ones(1,1), 0.6667,0.3333,zeros(1,16-4)+0,0.3333,0.6667, ones(1,1), 0.5,zeros(1,20-1)].';
        box_ca{3} = [zeros(1,20-1),0.5, ones(1,1), zeros(1,0)+1, ones(1,1), zeros(1,1)+1, ones(1,1), zeros(1,1)+0.9, ones(1,1), ...
                    zeros(1,2)+0.6, ones(1,1), zeros(1,4)+0.3, ones(1,1), 0.5,zeros(1,8-2)+0,0.5, ...
                    ones(1,1), 0.6667,0.3333,zeros(1,16-4)+0,0.3333,0.6667, ones(1,1), 0.5,zeros(1,20-1)].';
        for ind_ca = 1:length(box_ca)
            box = box_ca{ind_ca};
            box_n = imresize(box, [round(sig_dur * fr * 2), 1], 'nearest');
            box_l = imresize(box, [round(sig_dur * fr * 2), 1], 'bilinear');
            box = box_l; box(box_n==1) = 1;
            box = box - mean(box);
            box = box / norm(box);
            box_ca{ind_ca} = box;
        end
        nxc_mat = nan(length(box)+length(Z_merge)-1, length(box_ca));
        xc_mat = nan(length(box)+length(Z_merge)-1, length(box_ca));
        for ind_ca = 1:length(box_ca)
            [nxc, xc] = normxcorr(box_ca{ind_ca}, Z_merge - baseline_merge);
            nxc_mat(:,ind_ca) = nxc;
            xc_mat(:,ind_ca) = xc;
        end
        [~, ix_vec] = max(nxc_mat, [], 2);
        for ind_nxc = 1:length(nxc_mat)
            nxc(ind_nxc) = nxc_mat(ind_nxc, ix_vec(ind_nxc));
            xc(ind_nxc) = xc_mat(ind_nxc, ix_vec(ind_nxc));
        end
        lag = (0:length(nxc)-1).' / fr - (length(box)-1)/2 / fr;

        xcnxc = nxc .* xc .* sign(xc);
        xcnxc_thresh = auto_threshold(xcnxc, 0, 5) * 0.1;
        peaks_thresh = xcnxc_thresh;
        [peaks, indices] = findpeaks(xcnxc, 'MinPeakHeight', xcnxc_thresh, 'MinPeakDistance', sig_dur * fr);
        [~, order] = sort(peaks, 'descend');
        peaks = peaks(order);
        indices = indices(order);

        if plotting
            xmin = tr(1)+tr(end)/10;
            xmax = tr(end)-tr(end)/10;

            figure(1);
            ax(1) = subplot(4,1,1); cla; hold on;
%             plot(tr, abs(mean(Z_d,2)), 'DisplayName', 'Downsampled');
            plot(tr, abs(mean(Z_dc(:,f_inds_merge),2)), 'DisplayName', '500Hz comb filtered');
            plot(tr, abs(mean(Z_dcl(:,f_inds_merge),2)), 'DisplayName', 'Lowpass filtered');
%             plot(tr, Z_dclm, 'DisplayName', 'Median filtered');
            plot(tr, ub, 'DisplayName', 'Envelope max');
            plot(tr, lb, 'DisplayName', 'Envelope min');
            plot(tr, Z_noise, 'linewidth', 2.5, 'DisplayName', 'Noise only');
            plot(tr, baseline_merge, 'DisplayName', 'Baseline estimate');
            legend();
            ylabel('|Z| [Ohm]');
    
            ax(2) = subplot(4,1,2); cla; hold on;
            plot(tr, env, 'DisplayName', 'Signal envelope');
            plot([tr(1); tr(end)], [env_thresh; env_thresh], 'DisplayName', 'Envelope threshold');
            legend();
            ylabel('Signal envelope [Ohm]');
            xlabel('Time [s]');
    
            ax(3) = subplot(4,1,3); cla; hold on;
            plot(tr, Z_merge - baseline_merge, 'DisplayName', 'Baseline subtracted');
            if isfield(detections, 'time') && ~isempty([detections.time])
                plot([detections.time], zeros(size([detections.time])), 'r.', 'DisplayName', 'Detections');
            end
            ylabel('\Delta|Z| [Ohm]');
            xlabel('Time [s]')
            legend();

            ax(4) = subplot(4,1,4); cla; hold on;
%             % plot(lag, xc1 ./ max(xc1), 'DisplayName', 'Cross-correlation');
%             % plot(lag, nxc1, 'DisplayName', 'Normalized cross-correlation');
%             plot(lag, xcnxc1, 'DisplayName', 'Combined correlation 1');
%             %         plot([lag(1); lag(end)], [xcnxc1_thresh; xcnxc1_thresh], 'DisplayName', 'Correlation threshold');
%             plot(lag(inds1), pks1, 'k.', 'DisplayName', 'Candidate peaks 1');
%             plot(lag(indices1), peaks1, 'ko', 'DisplayName', 'Detected peaks 1');
%             plot([lag(1); lag(end)], [peaks1_thresh; peaks1_thresh], 'DisplayName', 'Min peak height 1');

            set(gca,'ColorOrderIndex',1);
            plot(lag, xcnxc, '-', 'DisplayName', 'Combined correlation');
%             plot(lag(inds1), pks1, 'kx', 'DisplayName', 'Candidate peaks');
            plot(lag(indices), peaks, 'kx', 'DisplayName', 'Detected peaks');
            plot([lag(1); lag(end)], [peaks_thresh; peaks_thresh], '--', 'DisplayName', 'Min peak height');

            xlabel('Lag [s]');
            ylabel('Cross-correlation');
            legend();
            linkaxes(ax, 'x');
            set(gca, 'xlim', [xmin, xmax]);
            shg;
        end
    
        cheap_amp_thresh = 1.0;
        nxc_thresh = 0.75;
        max_detections = 200;

        count = 0;
        count_skip = 0;
        count_rep = 0;
        count_dup = 0;
        count_amp = 0;
        count_nxc = 0;
        for j = 1:length(peaks)
            if count >= max_detections
                break;
            end
            if ~isempty(detections) && any(abs([detections.time] - lag(indices(j))) < 0.25*(sig_dur/2 + [detections.sig_dur]./2)) % too close to existing detection
                [~, best_matching_ind] = min(abs([detections.time] - lag(indices(j))));
                if nxc(indices(j)) > detections(best_matching_ind).nxc && xc(indices(j)) > detections(best_matching_ind).xc
                    % delete existing detection
                    detections(best_matching_ind) = [];
                    count = count - 1;
                    count_rep = count_rep + 1;
                elseif nxc(indices(j)) < detections(best_matching_ind).nxc && xc(indices(j)) < detections(best_matching_ind).xc
                    % skip current detection
                    count_skip = count_skip + 1;
                    continue;
                else
                    % add both and sort it out later...
                    count_dup = count_dup + 1;
                end
            end
            detection.time = lag(indices(j));
            detection.corr = peaks(j);
            detection.xc = xc(indices(j));
            detection.nxc = nxc(indices(j));
            detection.filt_type = 'box';
            detection.freq = freq_vec;
            detection.bw_lpf = bw_lpf;
            detection.sig_dur = sig_dur;
            detection.time_mask = (tr >= lag(indices(j)) - sig_dur) & (tr <= lag(indices(j)) + sig_dur);
            detection.time_vec = tr(detection.time_mask);
            detection.Z_d = Z_d(detection.time_mask,:);
            detection.Z_dc = Z_dc(detection.time_mask,:);
            detection.Z_dcl = Z_dcl(detection.time_mask,:);
            detection.Z_merge = Z_merge(detection.time_mask,:);
            detection.baseline_merge = baseline_merge(detection.time_mask);
            detection.baseline = baseline(detection.time_mask,:);
            detection.cheap_amp = max(detection.Z_merge) ./ median(detection.baseline_merge);
            if detection.cheap_amp < cheap_amp_thresh
                count_amp = count_amp + 1;
                continue;
            elseif detection.nxc < nxc_thresh
                count_nxc = count_nxc + 1;
                continue;
            end
            if isempty(detections)
                detections = detection;
            else
                detections(end+1) = detection;
            end
            count = count + 1;
        end
        fprintf('Detected %d box type events\n', count);
        fprintf('Skipped (due to coincidence) %d box type events\n', count_skip);
        fprintf('Replaced (despite coincidence) %d box type events\n', count_rep);
        fprintf('Saved as duplicate to sort out later %d box type events\n', count_dup);
        fprintf('Cheap amp threshold failed on %d box type events\n', count_amp);
        fprintf('Norm xcorr threshold failed on %d box type events\n', count_nxc);
        fprintf('Total detected events so far: %d\n', length(detections));
    end
    fprintf('Saving detections to %s\n', fn);
    save([out_path fn], 'detections', 'Z_merge', 'Z_d', 'Z_dc', 'Z_dcl', 'bw_lpf', 'fr', 'tr', 'sig_dur_vec', '-V7.3');

end

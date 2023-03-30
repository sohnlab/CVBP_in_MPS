close all;
clearvars;
clc;

%%

p = 'good_detections/';
out_p = 'fixed_detections/';
d = dir(p);
fns = {d.name};
mask = contains(fns, '.mat');
% mask = mask & contains(fns, '');
fns = fns(mask);

for i = 1:length(fns)

    fn = fns{i};
    if exist([out_p fn], 'file')
        fprintf('Skipping %s because it already exists in the fixed_detections folder!\n', fn)
        continue;
    end
    load([p fn]);

    Z_dc = unpad_data(comb_filter_keep_DC(pad_data(Z_d - Z_d(1,:), fr), 500, fr), fr) + Z_d(1,:);

    %%

    for j = 1:length(good_detections)
        
        dj = good_detections(j);
        if isfield(dj, 'edge_center_time') && ~isempty(dj.edge_center_time)
            fprintf('Skipping detection #%d because it already has the new timestamps!\n', j);
            continue;
        end
        x_start = dj.time - dj.sig_dur * 2;
        x_stop  = dj.time + dj.sig_dur * 5;
        x_left  = dj.time - dj.sig_dur/2;
        x_right = dj.time + dj.sig_dur/2;
    
        chunk_mask = tr >= x_start & tr <= x_stop;
        tr_chunk = tr(chunk_mask);
        Z_chunk = abs(Z_dc(chunk_mask,end-1));
        dur_med = dj.sig_dur * 40/1590/2;
        len_med = max([floor(dur_med * fr), 1]);
        Z_chunk = movmedian(Z_chunk, len_med, 1);
        Z_chunk = movmean(Z_chunk, len_med, 1);
        Z_baseline = robust_median(Z_chunk(tr_chunk < x_left | (tr_chunk > x_right & tr_chunk < x_right + (x_left - x_start))));
        Z_chunk = Z_chunk - Z_baseline;
        Z_chunk = Z_chunk / Z_baseline;
        d_Z_chunk = diff(Z_chunk);
        d_Z_chunk = movmedian(d_Z_chunk, 3*len_med, 1);
        d_Z_chunk = movmean(d_Z_chunk, 3*len_med, 1);
        tr_chunk_d = tr_chunk(1:end-1);
        d_mult = 25;
        mask_p = tr_chunk_d >= (dj.time - dj.sig_dur) & tr_chunk_d <= dj.time;
        [pks_p, t_p] = findpeaks(d_Z_chunk(mask_p), tr_chunk_d(mask_p), NPeaks=9, SortStr="descend");
        [t_p, ord_p] = sort(t_p, "ascend"); pks_p = pks_p(ord_p);
        mask_n = tr_chunk_d >= dj.time & tr_chunk_d <= (dj.time + dj.sig_dur*5);
        [pks_n, t_n] = findpeaks(-d_Z_chunk(mask_n), tr_chunk_d(mask_n), NPeaks=49, SortStr="descend");
        [t_n, ord_n] = sort(t_n, "ascend"); pks_n = -pks_n(ord_n);
        [~, ind_start_orig] = min(abs(t_p - (dj.time-dj.sig_dur/2)));
        [~, ind_stop_orig] = min(abs(t_n - (dj.time+dj.sig_dur/2)));
        
        figure(2); clf;
        subplot(2,1,1); hold on;
        plot(tr_chunk, Z_chunk, 'k-', Linewidth=2);
        set(gca, xlim=[x_start, x_stop]);
        axi = axis;
        plot(tr_chunk(1:end-1), d_Z_chunk * d_mult, 'k--');
        title(sprintf('%d of %d', j, length(good_detections)));
        for pn = 1:length(t_p)
            if pn == ind_start_orig
                lw = 2;
            else
                lw = 0.5;
            end
            plot(t_p(pn), pks_p(pn) * d_mult, 'b.');
            text(t_p(pn), axi(4), num2str(pn), Color='b', FontSize=16);
            plot([t_p(pn), t_p(pn)], axi(3:4), 'b-', Linewidth=lw);
        end
        for pn = 1:length(t_n)
            if pn == ind_stop_orig
                lw = 2;
            else
                lw = 0.5;
            end
            plot(t_n(pn), pks_n(pn) * d_mult, 'r.');
            text(t_n(pn), axi(3), num2str(pn), Color='r', FontSize=16);
            plot([t_n(pn), t_n(pn)], axi(3:4), 'r-', Linewidth=lw);
        end
        ylabel('( |Z| - |Z_0| ) /  |Z_0|');
        plot([x_left - dj.sig_dur/5, x_left - dj.sig_dur/5], axi(3:4), 'k-', Linewidth=2);
        plot([x_right + dj.sig_dur/2, x_right + dj.sig_dur/2], axi(3:4), 'k-', Linewidth=2);
        axis(axi);
    
        subplot(2,1,2); hold on;
        plot(tr_chunk, Z_chunk, 'k-', Linewidth=2);
        set(gca, xlim=[x_left - dj.sig_dur/5, x_right + dj.sig_dur/2]);
        axi = axis;
        plot(tr_chunk(1:end-1), d_Z_chunk * d_mult, 'k--');
        for pn = 1:length(t_p)
            if pn == ind_start_orig
                lw = 2;
            else
                lw = 0.5;
            end
            plot(t_p(pn), pks_p(pn) * d_mult, 'b.');
            text(t_p(pn), axi(4), num2str(pn), Color='b', FontSize=16);
            plot([t_p(pn), t_p(pn)], axi(3:4), 'b-', Linewidth=lw);
        end
        for pn = 1:length(t_n)
            if pn == ind_stop_orig
                lw = 2;
            else
                lw = 0.5;
            end
            plot(t_n(pn), pks_n(pn) * d_mult, 'r.');
            text(t_n(pn), axi(3), num2str(pn), Color='r', FontSize=16);
            plot([t_n(pn), t_n(pn)], axi(3:4), 'r-', Linewidth=lw);
        end
        ylabel('( |Z| - |Z_0| ) /  |Z_0|');
        axis(axi);

        str = input('Input 2 or 3 digits to select start time and end time.\nOtherwise default to thick lines.\n', 's');
        if ~isempty(str)
            ind_start = str2double(str(1));
            ind_stop = str2double(str(2:end));
        else
            ind_start = ind_start_orig;
            ind_stop = ind_stop_orig;
        end
        t_start = t_p(ind_start);
        t_stop = t_n(ind_stop);
        t_center = (t_start + t_stop) / 2;
        t_transit = t_stop - t_start;
        t_mask = tr >= t_center - t_transit*2.5 & tr <= t_center + t_transit*2.5;
        good_detections(j).Z_d = Z_d(t_mask,:);
        good_detections(j).Z_dc = Z_dc(t_mask,:);
        dur_filt = t_transit * 40 / 1590 / 2;
        len_filt = max([floor(dur_med * fr), 1]);
        Z_dcm = movmedian(Z_dc(t_mask,:), len_filt, 1);
        Z_dcml = movmean(Z_dcm, len_filt, 1);
        Z_dcl = movmean(Z_dc(t_mask,:), len_filt, 1);
        Z_dclm = movmedian(Z_dcl, len_filt, 1);
        good_detections(j).t_vec = tr(t_mask);
        good_detections(j).Z_dcm = Z_dcm;
        good_detections(j).Z_dcml = Z_dcml;
        good_detections(j).Z_dcl = Z_dcl;
        good_detections(j).Z_dclm = Z_dclm;
        good_detections(j).edge_center_time = t_center;
        good_detections(j).edge_start_time = t_start;
        good_detections(j).edge_stop_time = t_stop;
        good_detections(j).edge_transit_time = t_transit;

    end
    
    %% plot all signals to see if they are aligned well

    figure(99); clf; hold on;
    for ix = 1:length(good_detections)
        gd = good_detections(ix);
        if isfield(gd, 'edge_center_time') && ~isempty(gd.edge_center_time)
            Z_plot = abs(mean(gd.Z_dcml(:,end-1),2));
            plot((gd.t_vec - gd.edge_center_time) / gd.edge_transit_time, Z_plot - robust_median(Z_plot));
        end
    end
    axi = axis;
    plot([-0.5,-0.5], axi(3:4), 'k-');
    plot([0.5,0.5], axi(3:4), 'k-');

    %% save
    
    fprintf('Saving %s\n', fn);
    save([out_p fn], 'good_detections', 'good_mask', 'fr', 'tr', 'Z_merge', 'Z_d', 'Z_dc', 'min_nxc', 'min_cheap_amp', 'min_new_rel_amp', 'min_dur');
end

%% run this code block to remove deprecated fields from good_detections struct

% out_d = dir(out_p);
% out_fns = {out_d.name};
% mask = contains(out_fns, '.mat');
% out_fns = out_fns(mask);
% for i = 1:length(out_fns)
%     out_fn = out_fns{i};
%     load([out_p out_fn]);
%     good_detections = rmfield(good_detections, {'time','filt_type','bw_lpf','sig_dur','time_mask','time_vec','Z_merge','baseline_merge','baseline','cheap_amp'});
%     save([out_p out_fn], 'good_detections', 'good_mask', 'fr', 'tr', 'Z_merge', 'Z_d', 'Z_dc', 'min_nxc', 'min_cheap_amp', 'min_new_rel_amp', 'min_dur');
% end

%%

%% OLD CODE
% %     bw_lpf2 = 250;
%     Z_dc = unpad_data(comb_filter_keep_DC(pad_data(Z_d - Z_d(1,:), fr), 500, fr), fr) + Z_d(1,:);
%     Z_dc2 = unpad_data(comb_filter_keep_DC(pad_data(Z_d - Z_d(1,:), fr), 500, fr), fr) + Z_d(1,:);
% %     Z_dcl2 = unpad_data(lowpass(pad_data(Z_dc2 - Z_dc2(1,:), fr), bw_lpf2, fr), fr) + Z_dc2(1,:);
%     %
%     [~, order] = sort([good_detections.cheap_amp], "descend");
%     d1 = good_detections(order(1));
%     xstart = d1.time - d1.sig_dur*1;
%     xstop = d1.time + d1.sig_dur*1;
%     xleft = d1.time - d1.sig_dur/2;
%     xright = d1.time + d1.sig_dur/2;
% 
%     dur_med = d1.sig_dur * 40/1590/2;
%     len_med = floor(dur_med * fr);
% 
%     Z_dc2 = movmedian(Z_dc2, len_med, 1);
%     Z_dc2 = movmean(Z_dc2, len_med, 1);
% %     Z_dcl2 = movmean(Z_dc2, len_med, 1);
%     Z0 = median(Z_dcl2(tr >= xstart & tr < xleft,:), 1);
% %     Z_chunk = Z_dcl2(tr >= xstart & tr < xleft, :);
% %     [~, big_ind] = max(abs(mean(Z_chunk,2)));
% %     Z1 = Z_chunk(big_ind,:);
% %     DZ = Z1 - Z0;
% %     Z_dcl2 = Z_dcl2 .* conj(DZ./abs(DZ));
% %     Z0 = median(Z_dcl2(tr >= xstart & tr < xleft,:), 1);
% %     Z_dcl2 = Z_dcl2 .* exp(1i*-pi/2);
% %     Z0(:) = 0;
% %     ph0 = angle(Z0);
% %     Z_diff = Z_dcl2 - Z0;
% %     Z_plot = movmedian(Z_diff, len_med, 1);
% 
%     figure(1); clf;
%     ax(1) = subplot(2,2,1); hold on;
%     plot(tr, movmedian(real(Z_dcl2) - real(Z0), len_med, 1)./real(Z0));
%     set(gca, xlim=[xstart, xstop]);
%     axi = axis;
%     plot([xleft, xleft], axi(3:4), 'k');
%     plot([xright, xright], axi(3:4), 'k');
%     xlabel('Time [s]');
%     ylabel('Re(Z) - Re(Z_0) [Re(Z_0)]');
%     legend(d1.freq/1e3 + " kHz");
% 
%     ax(3) = subplot(2,2,3); hold on;
%     plot(tr, movmedian(imag(Z_dcl2) - imag(Z0), len_med, 1)./imag(Z0));
%     set(gca, xlim=[xstart, xstop]);
%     axi = axis;
%     plot([xleft, xleft], axi(3:4), 'k');
%     plot([xright, xright], axi(3:4), 'k');
%     xlabel('Time [s]');
%     ylabel('Im(Z) - Im(Z_0) [Im(Z_0)]');
%     legend(d1.freq/1e3 + " kHz");
% 
% %     ax(2) = subplot(2,3,2); hold on;
% %     plot(tr, real(movmedian(Z_dcl2 - Z0, len_med, 1))/1e3);
% %     set(gca, xlim=[xstart, xstop]);
% %     axi = axis;
% %     plot([xleft, xleft], axi(3:4), 'k');
% %     plot([xright, xright], axi(3:4), 'k');
% %     xlabel('Time [s]');
% %     ylabel('Re(Z - Z_0) [k\Omega]');
% %     legend(d1.freq/1e3 + " kHz");
% % 
% %     ax(5) = subplot(2,3,5); hold on;
% %     plot(tr, imag(movmedian(Z_dcl2 - Z0, len_med, 1))/1e3);
% %     set(gca, xlim=[xstart, xstop]);
% %     axi = axis;
% %     plot([xleft, xleft], axi(3:4), 'k');
% %     plot([xright, xright], axi(3:4), 'k');
% %     xlabel('Time [s]');
% %     ylabel('Im(Z - Z_0) [k\Omega]');
% %     legend(d1.freq/1e3 + " kHz");
% 
%     ax(2) = subplot(2,2,2); hold on;
% %     plot(tr, abs(Z_plot)/1e3);
%     plot(tr, movmedian(abs(Z_dcl2) - abs(Z0), len_med, 1)./abs(Z0));
%     set(gca, xlim=[xstart, xstop]);
%     axi = axis;
%     plot([xleft, xleft], axi(3:4), 'k');
%     plot([xright, xright], axi(3:4), 'k');
%     xlabel('Time [s]');
%     ylabel('|Z| - |Z_0| [|Z_0|]');
%     legend(d1.freq/1e3 + " kHz");
% 
%     ax(4) = subplot(2,2,4); hold on;
%     plot(tr, movmedian(angle(Z_dcl2) - angle(Z0), len_med, 1)./angle(Z0));
%     set(gca, xlim=[xstart, xstop]);
%     axi = axis;
%     plot([xleft, xleft], axi(3:4), 'k');
%     plot([xright, xright], axi(3:4), 'k');
%     xlabel('Time [s]');
%     ylabel('\angleZ - \angleZ_0 [\angleZ_0]');
%     legend(d1.freq/1e3 + " kHz");
%     linkaxes(ax, 'x');
%     pause(0);
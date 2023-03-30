% sort through detections

close all;
clearvars;
clc;

%% select matfile

p = 'detections/';
fn = '';
load([p fn]);

%% visualize detections

counts = zeros(size(sig_dur_vec));

% detections
min_nxc = 0.75;
min_cheap_amp = 1.0;
filt_type = 'box';
R = 10;
C = 10;
figure(1); clf;
for d = detections
    if d.nxc < min_nxc || ~strcmp(d.filt_type, filt_type) || d.cheap_amp < min_cheap_amp
        continue;
    end
    n = max([floor(find(sig_dur_vec == d.sig_dur)/6),1]);
    subplot(R, C, n); hold on;
    plot(d.time_vec - mean(d.time_vec), abs(d.Z_merge) - median(d.baseline_merge)); axis tight;
    counts(n) = counts(n) + 1;
end
for n = 1:length(sig_dur_vec)
    subplot(R, C, max([floor(n/6),1])); hold on;
    title(sprintf('%d at %0.1f ms', counts(n), sig_dur_vec(n)*1e3));
end

%% sort detections by criterion

[~, order] = sort([detections.nxc], 'descend');
detections = detections(order);

%% reset labels

good = nan(size(detections));
good_mask = nan(size(tr));

%% set parameters and start manual pruning

min_cheap_amp = 1.0001;
min_new_rel_amp = 1.0;
min_nxc = 0.75;
min_dur = 0.05;
filt_type = 'box';
for i = 1:length(detections)
    if ~isnan(good(i))
        continue; % already labeled this candidate detection
    end
    d = detections(i);
    if ~strcmp(d.filt_type, filt_type)
        continue; % not currently looking for this signal type
    end
    if d.sig_dur < min_dur
        continue; % candidate detection below sig_dur threshold
    end
    if any(interp1(tr, double(good_mask), (d.time-d.sig_dur/2):1/fr:(d.time+d.sig_dur/2), 'nearest') == 1)
        continue; % already detected nearby signal
    end
    if any(interp1(tr, double(good_mask), (d.time-d.sig_dur/2):1/fr:(d.time+d.sig_dur/2), 'nearest') == -1)
        continue; % already rejected nearby signal
    end
    if d.cheap_amp < min_cheap_amp
        continue; % candidate detection below amplitude threshold
    end
    if d.nxc < min_nxc
        continue; % candidate detection below nxc threshold
    end
    Z_m = movmedian(abs(Z_dc), round(d.sig_dur*40/1590/2*fr));
    wide_time_mask = tr >= d.time - d.sig_dur*3 & tr <= d.time + d.sig_dur*3;
    narrow_time_mask = tr >= d.time - d.sig_dur*0.75 & tr <= d.time + d.sig_dur*0.75;
    mm = max(Z_m(narrow_time_mask));
    mmm = min(Z_m(narrow_time_mask));
    new_rel_amp = mm/mmm;
    if new_rel_amp < min_new_rel_amp
        continue; % candidate detection below new rel amp threshold
    end
    disp(i);
    figure(99); clf; hold on;
    plot(tr(wide_time_mask), abs(Z_dc(wide_time_mask)), 'linewidth', 1.5);
    plot(tr(wide_time_mask), Z_m(wide_time_mask), 'linewidth', 2.5);
%     plot(d.time_vec, d.baseline, 'linewidth', 2.5);
    axis tight; a = axis;
    plot([d.time, d.time], a(3:4), 'k');
    plot([d.time, d.time] - d.sig_dur/2, a(3:4), 'k');
    plot([d.time, d.time] + d.sig_dur/2, a(3:4), 'k');
    title(sprintf('nxc = %0.2f%%, dur = %0.1f ms, new rel amp = %0.1f%%', d.nxc*100, d.sig_dur*1e3, (new_rel_amp-1)*100));
    str = input('press y to mark good or u to undo last\n', 's');
    switch str
        case 'y'
            good(i) = 1;
            good_mask(tr > d.time - d.sig_dur & tr < d.time + d.sig_dur) = 1;
            last_i = i;
        case 'n'
            good(i) = 0;
            good_mask(tr > d.time - d.sig_dur & tr < d.time + d.sig_dur) = -1;
            last_i = i;
        case 'u'
            good(last_i) = 1-good(last_i);
            if good(last_i)
                good_mask(tr > detections(last_i).time - detections(last_i).sig_dur/2 & tr < detections(last_i).time + detections(last_i).sig_dur/2) = 1;
            else
                good_mask(tr > detections(last_i).time - detections(last_i).sig_dur/2 & tr < detections(last_i).time + detections(last_i).sig_dur/2) = 0;
            end
            str2 = input('reversed previous label. press y to mark good\n', 's');
            if strcmp(str2, 'y')
                good(i) = 1;
                good_mask(tr > d.time - d.sig_dur & tr < d.time + d.sig_dur) = 1;
                last_i = i;
            elseif strcmp(str2, 'n')
                good(i) = 0;
                good_mask(tr > d.time - d.sig_dur & tr < d.time + d.sig_dur) = -1;
                last_i = i;
            else
                good(i) = 0;
                good_mask(tr > d.time - d.sig_dur & tr < d.time + d.sig_dur) = 0;
                last_i = i;
            end
        otherwise
            good(i) = 0;
            good_mask(tr > d.time - d.sig_dur/2 & tr < d.time + d.sig_dur/2) = 0;
            last_i = i;
    end
end

%% plot good mask over data

f_inds_merge = detections(1).freq < 400e3;
Z_dc_merge = abs(mean(Z_dc(:,f_inds_merge),2));
Z_dcl_merge = abs(mean(Z_dcl(:,f_inds_merge),2));
Z_dcl_merge_good = Z_dcl_merge;
Z_dcl_merge_good(good_mask~=1) = NaN;
Z_dcl_merge_bad = Z_dcl_merge;
Z_dcl_merge_bad(good_mask~=-1) = NaN;

figure(98); clf; hold on;
% plot(tr, Z_dc_merge, 'b');
plot(tr, Z_dcl_merge);
plot(tr, Z_dcl_merge_bad);
plot(tr, Z_dcl_merge_good);

shg;

%% find nearest detection to time point

t =69.5
[~, torder] = sort([detections.time]);
tind = find([detections(torder).time] > t, 100, 'first');
ind = torder(tind);
i = find(strcmp({detections(ind).filt_type}, 'box'), 1, 'first');
d = detections(ind(i))
figure(98); hold on;
plot(d.time_vec, d.Z_merge, 'linewidth', 2.5);
shg;

%% run this code block if the current detection is good

good(ind(i)) = 1;
good_mask(tr > d.time - d.sig_dur & tr < d.time + d.sig_dur) = 1;

%% save

good_detections = detections(good==1);
out_p = 'good_detections/';
save([out_p fn], 'good_detections', 'good_mask', 'fr', 'tr', 'Z_merge', 'Z_d', 'Z_dc', 'Z_dcl', 'bw_lpf', 'min_nxc', 'min_cheap_amp', 'min_new_rel_amp', 'min_dur');

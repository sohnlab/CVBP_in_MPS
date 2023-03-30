close all;
clearvars;
clc;

%% load data

load('pad_array_detections.mat');

addpath('../functions');

%% plot a smaller signal

fig1 = figure; clf; hold on;
ind = 174;
D(ind)
for i_f = 1:size(D(ind).Z_wmm,2)
    rm_sig = D(ind).Z_dcml(:,i_f);
    rm_sig(round(1:end*1/5)) = nan; rm_sig(round(end*2/5+1:end*3/5)) = nan; rm_sig(round(end*4/5+1:end)) = nan;
    rm_sig(isnan(rm_sig)) = [];
    new_b_vec(i_f) = robust_median(abs(rm_sig));
end
Z_plot = D(ind).Z_sig(:,1:end-1);
for i_f = 1:size(Z_plot,2)
    Z_plot(:, i_f) = wden(abs(Z_plot(:,i_f)), 'modwtsqtwolog', 'h', 'mln', fix(log2(length(Z_plot(:,i_f)))), 'db20');
end
filt_len = floor(min(PAD_DURS(:,ind))*2.5e3*0.25);
Z_plot = movmedian(Z_plot, filt_len, 1);
Z_plot = movmean(Z_plot, filt_len, 1);
for j = 1:6
    subplot(6,2,(j-1)*2+1);
    plot(D(ind).t_sig, (abs(Z_plot(:,j))./abs(new_b_vec(j))-1)*100, '-');
    title(sprintf('%0.1f kHz', D(ind).freq(j)/1e3));
    set(gca, XAxisLocation="origin");
    set(gca, YLim=[-0.3,1.2]);
    set(gca, XLim=D(ind).t_sig([1,end]));
    xlabel('Time [s]');
    ylabel('\Delta|Z|/|Z| (%)');
end
shg;

% plot a bigger signal

ind = 169;
D(ind)
for i_f = 1:size(D(ind).Z_wmm,2)
    rm_sig = D(ind).Z_dcml(:,i_f);
    rm_sig(round(1:end*1/5)) = nan; rm_sig(round(end*2/5+1:end*3/5)) = nan; rm_sig(round(end*4/5+1:end)) = nan;
    rm_sig(isnan(rm_sig)) = [];
    new_b_vec(i_f) = robust_median(abs(rm_sig));
end
Z_plot = D(ind).Z_sig(:,1:end-1);
for i_f = 1:size(Z_plot,2)
    Z_plot(:, i_f) = wden(abs(Z_plot(:,i_f)), 'modwtsqtwolog', 'h', 'mln', fix(log2(length(Z_plot(:,i_f)))), 'db20');
end
filt_len = floor(min(PAD_DURS(:,ind))*2.5e3*0.25);
Z_plot = movmedian(Z_plot, filt_len, 1);
Z_plot = movmean(Z_plot, filt_len, 1);
for j = 1:6
    subplot(6,2,(j-1)*2+2);
    plot(D(ind).t_sig, (abs(Z_plot(:,j))./abs(new_b_vec(j))-1)*100, '-');
    title(sprintf('%0.1f kHz', D(ind).freq(j)/1e3));
    set(gca, XAxisLocation="origin");
    set(gca, YLim=[-0.9,3.6]);
    set(gca, XLim=D(ind).t_sig([1,end]));
    xlabel('Time [s]');
    ylabel('\Delta|Z|/|Z| (%)');
end

shg;

savename = 'output/pad_array_transit_examples_figure';
saveas(gcf, [savename '.fig']);
saveas(gcf, [savename '.png']);
saveas(gcf, [savename '.pdf']);

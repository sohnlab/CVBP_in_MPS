close all;
clearvars;
clc;

%% load data

load('bw1000Hz/G60-P80-G60_dev005a_MCF-7_in_1XPBS+2percFBS_run19_0-10mbar_1Vpp_4sines_sawatleast20good.mat');

addpath('../functions/');

%% filter data and compute baseline

t_mask = tr >= 78.3 & tr < 78.98;

% Z_filt = comb_filter_keep_DC(Z_mat(t_mask,:), 1000, fr, 10);
Z_filt = movmean(Z_filt, 0.001*fr, 1);
Z_baseline = nan(size(Z_filt));
for i = 1:4
    Z_baseline(:,i) = robust_median(Z_mat(t_mask,i));
end

%% plot impedance

fig1 = figure; clf;
for i = 1:4
    ax1(i) = subplot(4,1,i); hold on;
    plot(tr(t_mask), abs(Z_filt(:,i)/1e3));
    plot(tr(t_mask), abs(Z_baseline(:,i)/1e3));
    ylabel('|Z| [k\Omega]');
    legend([freq_vec(i)/1e3 + "kHz", "baseline"]);
    minZ = min(abs(Z_filt(:,i)/1e3));
    maxZ = max(abs(Z_filt(:,i)/1e3));
    rangeZ = range(abs(Z_filt(:,i)/1e3));
    set(gca, YLim=[minZ-rangeZ*0.1, maxZ+rangeZ*0.1]);
    set(gca, XLim=[78.625, 78.653]);
    set(gca, XTick=[78.63, 78.64, 78.65]);
end
xlabel('Time [s]');
linkaxes(ax1, 'x');
shg;

savename = 'output/single_pad_MPS_example';
saveas(gcf, [savename '.fig']);
saveas(gcf, [savename '.png']);
saveas(gcf, [savename '.pdf']);

%% plot normalized impedance

fig2 = figure; clf;
Z_norm = abs(Z_filt)./abs(Z_baseline) - 1;
% Z_norm = real(Z_filt)./real(Z_baseline) - 1;
% Z_norm = imag(Z_filt)./imag(Z_baseline) - 1;
for i = 1:4
    ax2(i) = subplot(4,1,i); hold on;
    plot(tr(t_mask), Z_norm(:,i)*100);
    ylabel('\Delta|Z|/|Z| (%)');
    legend(freq_vec(i)/1e3 + "kHz");
    minZ = min(Z_norm(:,i)*100);
    maxZ = max(Z_norm(:,i)*100);
    rangeZ = range(Z_norm(:,i)*100);
    set(gca, YLim=[minZ-rangeZ*0.1, maxZ+rangeZ*0.1]);
    set(gca, XAxisLocation='origin');
    set(gca, XLim=[78.625, 78.6525]);
    set(gca, XTick=[78.63, 78.64, 78.65]);
end
xlabel('Time [s]');
linkaxes(ax2, 'x');
shg;

savename = 'output/single_pad_MPS_example_normalized';
saveas(gcf, [savename '.fig']);
saveas(gcf, [savename '.png']);
saveas(gcf, [savename '.pdf']);

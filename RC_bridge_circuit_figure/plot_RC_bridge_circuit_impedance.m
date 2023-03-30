close all;
clearvars;
clc;

%% load data

Ra1 = 10.01e3;
Ra2 = 9.999e3;
Cb1 = 9.984e-9;
Cb2 = 9.92e-9;

Ra = mean([Ra1,Ra2]);
Cb = mean([Cb1,Cb2]);

f_Z = readmatrix('RC_bridge_no_e_good_compensation.xlsx', Range="A18:E118");
f = f_Z(:,1);
Z_open = f_Z(:,4) + 1i .* f_Z(:,5);

f_Z = readmatrix('RC_bridge_with_e_good_compensation.xlsx', Range="A18:E118");
% f = freq_Z(:,1);
Z_closed = f_Z(:,4) + 1i .* f_Z(:,5);

x = 2*pi*f * Ra * Cb;
z_open = Z_open / Ra;
z_closed = Z_closed / Ra;

z_open_th = 0.5 * (Ra + 1./(1i*2*pi*f*Cb)) / Ra;
z_closed_th = (1/Ra + (1i*2*pi*f*Cb)).^(-1) * 2 / Ra;

%% plot

fig1 = figure; clf; hold on;
plot(x, abs(z_open), '.', 'markersize', 20);
plot(x, abs(z_closed), '.', 'markersize', 20);
set(gca,'xscale','log','yscale','log');
set(gca, 'ColorOrderIndex', 1);
plot(x, abs(z_open_th), '-', 'linewidth', 2);
plot(x, abs(z_closed_th), '-', 'linewidth', 2);
plot([2-sqrt(3), 2-sqrt(3)], [0.2, 5], 'k-');
plot([2+sqrt(3), 2+sqrt(3)], [0.2, 5], 'k-');

axis([0.1,10,0.2,5]);
xticks([0.1,1,10]);
xticklabels({'0.1','1','10'});
yticks([0.25,0.5,1,2,4]);
xlabel('\omegaRC');
ylabel('|Z|/R');

shg;

savename = 'output/RC_bridge_circuit_impedance';
saveas(gcf, [savename '.fig']);
saveas(gcf, [savename '.png']);
saveas(gcf, [savename '.pdf']);

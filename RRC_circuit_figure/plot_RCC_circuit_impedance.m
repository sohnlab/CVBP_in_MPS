close all;
clearvars;
clc;

%% load data

Ra1 = 10.01e3;
Ra2 = 9.999e3;
Cb1 = 9.984e-9;
% Cb2 = 9.92e-9;

Ra = mean([Ra1,Ra2]);
% Cb = mean([Cb1,Cb2]);
Cb = Cb1;

f_Z = readmatrix('RRC_no_e_good_compensation.xlsx', Range="A18:E118");
f = f_Z(:,1);
Z_open = f_Z(:,4) + 1i .* f_Z(:,5);

f_Z = readmatrix('RRC_with_e_good_compensation.xlsx', Range="A18:E118");
% f = freq_Z(:,1);
Z_closed = f_Z(:,4) + 1i .* f_Z(:,5);

x = 2*pi*f * Ra * Cb;
z_open = Z_open / Ra;
z_closed = Z_closed / Ra;

z_open_th = (Ra + 1./(1i*2*pi*f*Cb)) / Ra;
z_closed_th = ((1/Ra + (1i*2*pi*f*Cb)).^(-1) + Ra) / Ra;

%% plot

fig1 = figure; clf; hold on;
plot(x, abs(z_open), '.', 'markersize', 20);
plot(x, abs(z_closed), '.', 'markersize', 20);
set(gca,'xscale','log','yscale','log');
set(gca, 'ColorOrderIndex', 1);
plot(x, abs(z_open_th), '-', 'linewidth', 2);
plot(x, abs(z_closed_th), '-', 'linewidth', 2);
plot([1/sqrt(2), 1/sqrt(2)], [0.9, 5.7], 'k-');

axis([0.2,10,0.9,5.7]);
xticks([0.1,1,10]);
xticklabels({'1','10'});
yticks([1,2,3,4]);
xlabel('\omegaRC');
ylabel('|Z|/R');

shg;

savename = 'output/RRC_circuit_impedance';
saveas(gcf, [savename '.fig']);
saveas(gcf, [savename '.png']);
saveas(gcf, [savename '.pdf']);

function plotBodeNyquist(f, Z_mat, c, varargin)

f = f(:); % force column
if size(Z_mat,1) == 1 && size(Z_mat,2) > 1 % row vector
    Z_mat = Z_mat.'; % force column
    onedim = true;
elseif size(Z_mat,2) == 1 && size(Z_mat,1) > 1 % column vector
    onedim = true;
else
    onedim = false;
end

if nargin == 2, c = '.-'; end

mag_mean = nanmean(abs(Z_mat), 2);
mag_std = nanstd(abs(Z_mat), [], 2);
ph_mean = nanmean(angle(Z_mat), 2)/pi*180;
ph_std = nanstd(angle(Z_mat), [], 2)/pi*180;

a(1) = subplot(2,2,1); hold on;
if onedim
    plot(f, mag_mean, c);
else
    errorbar(f, mag_mean, mag_std, c);
end
ylabel('|Z| [\Omega]');
set(gca, 'xscale', 'log', 'yscale', 'log');
axis tight; grid on; grid minor;

a(2) = subplot(2,2,3); hold on;
if onedim
    plot(f, ph_mean, c);
else
    errorbar(f, ph_mean, ph_std, c);
end
ylabel('\angleZ [\circ]');
xlabel('f [Hz]');
set(gca, 'xscale', 'log', 'yscale', 'linear');
axis tight; grid on; grid minor;
set(gca, 'ylim', [-90, 0]);

linkaxes(a, 'x');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

real_mean = nanmean(real(Z_mat), 2);
real_std = nanstd(real(Z_mat), [], 2);
imag_mean = nanmean(imag(Z_mat), 2);
imag_std = nanstd(imag(Z_mat), [], 2);

subplot(1,2,2); hold on;

if all(real_std==0)
    plot(real_mean, -imag_mean, c);
else
    errorbarxy(real_mean, -imag_mean, real_std, imag_std, c);
end
ylabel('-Im(Z) [\Omega]');
xlabel('Re(Z) [\Omega]');
grid on; grid minor;
axis equal;
set(gca, 'XAxisLocation', 'origin', 'YAxisLocation', 'origin');
    
end

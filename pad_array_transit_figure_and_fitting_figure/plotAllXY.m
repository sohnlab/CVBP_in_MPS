function plotAllXY(X_PLOT, Y_PLOT, f_vec, plotstr)

if nargin == 3, plotstr = '.'; end

l_vec = [10, 20, 40, 80, 160, 320, 640] * 1e-6;

x_min = min(X_PLOT, [], "all");
x_max = max(X_PLOT, [], "all");
y_min = min(Y_PLOT, [], "all");
y_max = max(Y_PLOT, [], "all");

ax(8) = subplot(2, 4, 8); hold on;
plot(nan(2,length(f_vec)), nan(2,length(f_vec)), '.');
plot(nan(2,2), nan(2,2), 'k-');
legend([f_vec/1e3 + " kHz", "y = 0", "y = x"]);

% get aspect ratio and set x_min, x_max, y_min, y_max
xr_minmax = x_max - x_min;
yr_minmax = y_max - y_min;
axis equal;
axi = axis;
ar_xy = (axi(2) - axi(1)) / (axi(4) - axi(3)); % aspect ratio of subplot (XR/YR)
ar_xy_minmax = xr_minmax / yr_minmax; % current aspect ratio from data
if ar_xy_minmax < ar_xy % current x limits need to stretch
    xm = (x_min + x_max) / 2;
    xr_new = xr_minmax / ar_xy_minmax * ar_xy;
    x_min = xm - xr_new/2;
    x_max = xm + xr_new/2;
elseif ar_xy_minmax > ar_xy % current y limits need to stretch
    ym = (y_min + y_max) / 2;
    yr_new = yr_minmax * ar_xy_minmax / ar_xy;
    y_min = ym - yr_new/2;
    y_max = ym + yr_new/2;
end
axis([x_min, x_max, y_min, y_max]);

for i_pad = 1:length(l_vec)
    ax(i_pad) = subplot(2, 4, i_pad); hold on;
    for i_freq = 1:length(f_vec)
        plot(squeeze(X_PLOT(i_pad, i_freq, :)), squeeze(Y_PLOT(i_pad, i_freq, :)), plotstr);
    end
    plot([0,x_max], [0,0], 'k-');
    plot([0,y_max], [0,y_max], 'k-');
    axis equal;
    axis([x_min, x_max, y_min, y_max]);
    xlabel('Gap Response (%)');
    ylabel('Pad Response (%)');
    title(sprintf('%d um pad', l_vec(i_pad)*1e6));
end
linkaxes(ax, 'xy');

end
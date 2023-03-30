function [a, b, c] = find_abc(x_cell, d_cell, x_l, l_vec, w_vec, h)

plotting = false;
% plotting = true;

% NOTE: b is the minor axis, NOT the minor semi-axis
ellipse_width = @(x, b) real( 2.*sqrt( ( 1 - (x - x_cell).^2 ./ (sqrt(d_cell^3/b)./2).^2 ) .* (b./2).^2 ) );
channel_width = @(x) channel_width_func(x, x_l, l_vec, w_vec);

b = d_cell; % initial value for b is the diameter
a = sqrt(d_cell^3/b);

x_corners = x_l + cumsum(l_vec(1:end-1));

if plotting
    x_left = x_cell - a/2; % left edge of cell
    x_right = x_cell + a/2; % right edge of cell
    x_vec = linspace(x_left, x_right, 25); % vector of x values to check for collision
    x_vec = unique(sort([x_vec, x_corners, x_corners-eps, x_corners+eps]), 'stable');
else
    x_vec = unique(sort([x_cell, x_corners, x_corners-eps, x_corners+eps]), 'stable');
end

channel_width_vec = cell2mat(arrayfun(channel_width, x_vec, 'UniformOutput', false));
ellipse_width_vec = ellipse_width(x_vec, b);

if plotting
    figure; clf; hold on;
    plot(x_vec, channel_width_vec/2, 'k.');
    plot(x_vec, -channel_width_vec/2, 'k.');
    plot(x_vec, ellipse_width_vec/2, 'b.');
    plot(x_vec, -ellipse_width_vec/2, 'b.');
    axis equal;
    shg;
end

count = 0;
while any(ellipse_width_vec > channel_width_vec) && count < 10 % collision detected

    [~, ind_max] = max(ellipse_width_vec - channel_width_vec);
    mult_b = channel_width_vec(ind_max) / ellipse_width_vec(ind_max); % multiply ellipse b by this value to fit

    b = b * mult_b;
    b = floor(b*1e9)/1e9;
    a = sqrt(d_cell^3/b);
    ellipse_width = @(x, b) real( 2.*sqrt( ( 1 - (x - x_cell).^2 ./ (sqrt(d_cell^3/b)./2).^2 ) .* (b./2).^2 ) );

    if plotting
        x_left = x_cell - a/2; % left edge of cell
        x_right = x_cell + a/2; % right edge of cell
        x_vec = linspace(x_left, x_right, 25); % vector of x values to check for collision
        x_vec = unique(sort([x_vec, x_corners, x_corners-eps, x_corners+eps]), 'stable');
    else
        x_vec = unique(sort([x_cell, x_corners, x_corners-eps, x_corners+eps]), 'stable');
    end

    channel_width_vec = cell2mat(arrayfun(channel_width, x_vec, 'UniformOutput', false));
    ellipse_width_vec = ellipse_width(x_vec, b);

    if plotting
        figure; clf; hold on;
        plot(x_vec, channel_width_vec/2, 'k.-');
        plot(x_vec, -channel_width_vec/2, 'k.-');
        plot(x_vec, ellipse_width_vec/2, 'b.-');
        plot(x_vec, -ellipse_width_vec/2, 'b.-');
        axis equal;
        shg;
    end

    count = count + 1;

end

%% now check for collision with ceiling/floor (NOTE: only scalar height is supported)
if a > h % collision, need to update to a, b, c ellipsoid
    c = h; % set ellipsoid z-axis to channel height
    b = sqrt(d_cell^3/c);
    a = b;
    ellipse_width_abc = @(x, b, c) real( 2.*sqrt( ( 1 - (x - x_cell).^2 ./ ((d_cell^3/b/c)./2).^2 ) .* (b./2).^2 ) );

    if plotting
        x_left = x_cell - a/2; % left edge of cell
        x_right = x_cell + a/2; % right edge of cell
        x_vec = linspace(x_left, x_right, 25); % vector of x values to check for collision
        x_vec = unique(sort([x_vec, x_corners, x_corners-eps, x_corners+eps]), 'stable');
    else
        x_vec = unique(sort([x_cell, x_corners, x_corners-eps, x_corners+eps]), 'stable');
    end

    channel_width_vec = cell2mat(arrayfun(channel_width, x_vec, 'UniformOutput', false));
    ellipse_width_vec = ellipse_width_abc(x_vec, b, c);

    count = 0;
    while any(ellipse_width_vec > channel_width_vec) && count < 10 % collision detected

        [~, ind_max] = max(ellipse_width_vec - channel_width_vec);
        mult_b = channel_width_vec(ind_max) / ellipse_width_vec(ind_max); % multiply ellipse b by this value to fit

        b = b * mult_b;
        b = floor(b*1e9)/1e9;
        a = d_cell^3/b/c;

        if plotting
            x_left = x_cell - a/2; % left edge of cell
            x_right = x_cell + a/2; % right edge of cell
            x_vec = linspace(x_left, x_right, 25); % vector of x values to check for collision
            x_vec = unique(sort([x_vec, x_corners, x_corners-eps, x_corners+eps]), 'stable');
        else
            x_vec = unique(sort([x_cell, x_corners, x_corners-eps, x_corners+eps]), 'stable');
        end

        channel_width_vec = cell2mat(arrayfun(channel_width, x_vec, 'UniformOutput', false));
        ellipse_width_vec = ellipse_width_abc(x_vec, b, c);

        if plotting
            figure; clf; hold on;
            plot(x_vec, channel_width_vec/2, 'k.-');
            plot(x_vec, -channel_width_vec/2, 'k.-');
            plot(x_vec, ellipse_width_vec/2, 'b.-');
            plot(x_vec, -ellipse_width_vec/2, 'b.-');
            axis equal;
            shg;
        end

        count = count + 1;

    end
else
    c = a;
end

end
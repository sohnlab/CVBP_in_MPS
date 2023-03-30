% forward model to produce simulated scatter plots from SPICE volume-based model

function [Z_none, Z_gaps, Z_pads, X_data, Y_data] = pad_array_vol_forward_model(l_chan, w, h, rho, q_vec, a_vec, f_vec, d_vec, dR_vec)

n_pads = (length(l_chan)-1)/2;
n_freqs = length(f_vec);
n_diams = length(d_vec);

l_elec = 1000e-6;
w_elec = 1000e-6;
l_marg = 25e-6;
w_marg = 1000e-6;
l_ramp = 150e-6;
w_ramp = 700e-6;
l_node = 100e-6;
w_node = 400e-6;
% l_chan = [40, 10, 40, 20, 40, 40, 40, 80, 40, 160, 40, 320, 40, 640, 40] * 1e-6;
w_chan = w + zeros(size(l_chan));

l_vec = [l_elec, l_marg, l_ramp, l_node, l_chan, l_node, l_ramp, l_marg, l_elec];
w_vec = [w_elec, w_marg, w_ramp, w_node, w_chan, w_node, w_ramp, w_marg, w_elec];

node_num_vec = nan(1, length(l_chan));
node_num_vec(2:2:end-1) = 2:((n_pads)+1);
node_num_vec = [1, nan, nan, nan, node_num_vec, nan, nan, nan, 0];

if isscalar(q_vec), q_vec = repmat(q_vec, size(l_vec)); q_vec(isnan(node_num_vec)) = nan; end
if isscalar(a_vec), a_vec = repmat(a_vec, size(l_vec)); a_vec(isnan(node_num_vec)) = nan; end

% compute Z_none of each segment
n_steps = 250;
n_segs = length(node_num_vec);
Z_segs = nan(n_segs, n_freqs);
parfor i = 1:n_segs
    l_ = l_vec(i);
    w_ = w_vec(i);
    q_ = q_vec(i);
    a_ = a_vec(i);
    node_num_ = node_num_vec(i);
    if i == 1 || i == n_segs % this is one of the electrodes
        l_vec_ = [l_, 1e-9];
        w_vec_ = [w_, 1e0];
        node_num_vec_ = [1, nan];
        q_vec_ = [q_, nan];
        a_vec_ = [a_, nan];
        Z_segs(i,:) = SPICE_VOL_MPS_RCPE_slices(l_vec_, w_vec_, h, node_num_vec_, rho, q_vec_, a_vec_, 0, -Inf, 0, f_vec, n_steps);
    elseif ~isnan(node_num_) % this is a pad segment
        l_vec_ = [1e-9, l_, 1e-9];
        w_vec_ = [1e0, w_, 1e0];
        node_num_vec_ = [nan, 2, nan];
        q_vec_ = [nan, q_, nan];
        a_vec_ = [nan, a_, nan];
        Z_segs(i,:) = SPICE_VOL_MPS_RCPE_slices(l_vec_, w_vec_, h, node_num_vec_, rho, q_vec_, a_vec_, 0, -Inf, 0, f_vec, n_steps);
    else % this is a gap segment or a resistive segment outside of the channel
        Z_segs(i,:) = rho * l_ / (w_ * h); % no need to use SPICE model
    end
end

% compute impedance when there is no particle in the channel
Z_none = sum(Z_segs, 1);

count = 0;
% tic;

% for each diameter, recalculate Z for intersecting segments
x_vec_edges = [0, cumsum(l_vec)]; % segment edge locations
x_vec = l_elec + l_marg + l_ramp + l_node + cumsum(l_chan) - l_chan/2; % particle locations
n_locs = length(x_vec);
N = n_diams * n_locs; % total number of circuit computations to do
Z_mat = nan(N, n_freqs);
fprintf('pad_array_vol_forward_model: computing %d circuits\n', N);
fprintf('%s|\n', char(' ' + zeros(1,N-1)));
ixid_mat = [kron(1:n_locs, ones(1, n_diams)); repmat(1:n_diams, [1, n_locs])];
xd_mat = [kron(x_vec, ones(1, n_diams)); repmat(d_vec, [1, n_locs])];
parfor ixid = 1:N
    d_ = xd_mat(2, ixid);
    x_ = xd_mat(1, ixid);
    [a_cell, b_cell, c_cell] = find_abc(x_, d_, 0, l_vec, w_vec, h);
    w_ = channel_width_func(x_, 0, l_vec, w_vec);
    h_ = h;
    dR_ = interp1(d_vec, dR_vec, d_, "linear", "extrap");
    
    xl_ = x_ - a_cell / 2; % left edge of particle
    xr_ = x_ + a_cell / 2; % right edge of particle
    mask_ = xl_ < x_vec_edges(2:end) & xr_ > x_vec_edges(1:end-1); % mask of intersecting segments
    if any(mask_) % there is at least one intersecting segment
        % compute impedance of intersecting segments with particle
        if mask_(1) && mask_(end) % somehow both voltage electrodes are intersecting
            Z_mat(ixid,:) = SPICE_VOL_MPS_RCPE_slices(l_vec, w_vec, h, node_num_vec, rho, q_vec, a_vec, d_, x_, dR_, f_vec, n_steps);
        elseif mask_(1) % left voltage electrode is intersecting
            l_vec_ = [l_vec(mask_), 1e-9];
            w_vec_ = [w_vec(mask_), 1e0];
            node_num_vec_ = [node_num_vec(mask_), nan];
            q_vec_ = [q_vec(mask_), nan];
            a_vec_ = [a_vec(mask_), nan];
            Z_mask = SPICE_VOL_MPS_RCPE_slices(l_vec_, w_vec_, h, node_num_vec_, rho, q_vec_, a_vec_, d_, x_, dR_, f_vec, n_steps);
            Z_mat(ixid,:) = sum(Z_segs(~mask_, :), 1) + Z_mask;
        elseif mask_(end) % right voltage electrode is intersecting
            l_vec_ = [1e-9, l_vec(mask_)];
            w_vec_ = [1e0, w_vec(mask_)];
            node_num_vec_ = [nan, node_num_vec(mask_)];
            q_vec_ = [nan, q_vec(mask_)];
            a_vec_ = [nan, a_vec(mask_)];
            x_off = x_vec_edges(find(mask_, 1, "first")) - l_vec_(1); % account for x offset
            Z_mask = SPICE_VOL_MPS_RCPE_slices(l_vec_, w_vec_, h, node_num_vec_, rho, q_vec_, a_vec_, d_, x_-x_off, dR_, f_vec, n_steps);
            Z_mat(ixid,:) = sum(Z_segs(~mask_, :), 1) + Z_mask;
        elseif all(isnan(node_num_vec(mask_))) % only intersecting one or more gaps
            Z_mat(ixid,:) = Z_none + dR_; % skip SPICE calculation
        else % neither voltage electrode is intersecting
            l_vec_ = [1e-9, l_vec(mask_), 1e-9];
            w_vec_ = [1e0, w_vec(mask_), 1e0];
            node_num_vec_ = [nan, node_num_vec(mask_), nan];
            q_vec_ = [nan, q_vec(mask_), nan];
            a_vec_ = [nan, a_vec(mask_), nan];
            x_off = x_vec_edges(find(mask_, 1, "first")) - l_vec_(1); % account for x offset
            Z_mask = SPICE_VOL_MPS_RCPE_slices(l_vec_, w_vec_, h, node_num_vec_, rho, q_vec_, a_vec_, d_, x_-x_off, dR_, f_vec, n_steps);
            Z_mat(ixid,:) = sum(Z_segs(~mask_, :), 1) + Z_mask;
        end
    else % somehow there are no intersecting segments
        Z_mat(ixid,:) = Z_none;
    end
    count = count + 1;
%         t = toc;
%         fprintf('Computed Z from SPICE model (%d of %d)\nTime elapsed: %0.0f min\nTime remaining: %0.0f min\n', ... 
%                 count, length(d_vec)*length(x_vec), t/60, t/count*(length(d_vec)*length(x_vec))/60-t/60);
    fprintf('.');
end
fprintf('\n');

% reshape Z_mat
Z_mat = permute(reshape(Z_mat, [n_diams, n_locs, n_freqs]), [2 3 1]);

Z_gaps = Z_mat(1:2:end,:,:);
Z_pads = Z_mat(2:2:end-1,:,:);

X_data = ((0.5*abs(Z_gaps(1:end-1,:,:)) + 0.5*abs(Z_gaps(2:end,:,:))) ./ abs(Z_none) - 1) * 100;
Y_data = (abs(Z_pads) ./ abs(Z_none) - 1) * 100;

end
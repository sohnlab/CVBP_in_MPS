% forward model for template matching

function [template, edge_times] = forward_model(t_vec, tc, tt, fext, fmis, facc, amps, w, h, f_vec, dR, rho, cpad)

% parameters:
%   tc    - center time
%   tt    - transit time (from edge to edge)
%   fext  - extra time fraction for high signal at the start (as fraction of the first gap)
%   fmis  - misalignment fraction (add to first gap and subtract from last gap)
%   facc  - acceleration factor (velocity ratio = v_end / v_start)
%   amps  - amplitude scalar of each channel (size [1,K])
%   w     - channel width [m]
%   h     - channel height [m]
%   f_vec - vector of frequencies [Hz] (size [1,K/2])
%   dR    - deltaR of praticle [Ohm]
%   rho   - resistivity [Ohm*m]
%   cpad  - specific capacitance of pad [F/m^2]

% hard-code rise time because it doesn't need to be a parameter
tr = tt * 10 / 1590 / 3;

% construct nominal edges times [0, 1]
l_vec = [40, 10, 40, 20, 40, 40, 40, 80, 40, 160, 40, 320, 40, 640, 40].' .* 1e-6; % [m]
seg_times = l_vec.'; % from pad array geometry
% gap_inds = 1:2:length(seg_times); % which indices of seg_times correspond to gaps

% use simple pad array model to calculate gap and pad responses
[Z0_vec, Z_mat] = simple_pad_array_model(l_vec, w, h, f_vec, dR, rho, cpad);
resp_mat = Z_mat - Z0_vec;
resp_mat = [real(resp_mat), imag(resp_mat)]; % stack real and imag channels
resp_mat = resp_mat ./ repmat(abs(Z0_vec), [1,2]); % scale to have maximum amplitude of 1

seg_times(1) = seg_times(1) + fmis*seg_times(1); % lengthen first gap
seg_times(end) = seg_times(end) - fmis*seg_times(1); % shorten last gap by same amount
seg_times(1) = seg_times(1) + fext*seg_times(1); % lengthen first gap
seg_times = seg_times / sum(seg_times); % normalized to sum to 1
edge_times = [0, cumsum(seg_times)];
if facc ~=1
    acc_warp = @(t) (sqrt((facc-1)*(1+facc)*t + 1) - 1) / (facc-1);
    edge_times = acc_warp(edge_times);
end
% convert to real edge times
edge_times = edge_times * tt + (tc - tt/2);

% solve for logistic piece parameters
dur_edge = tr * 3; % define duration of edge as 3X rise time
den = 5 * sqrt(7/11) / 4; % modified logistic denominator
off = (44 - 5*sqrt(77)) / 88; % modified logistic vertical offset
a = 2 * log( (8 + sqrt(77)) / 2 ) / tr; % logistic parameter

edge_func = @(left, right, t_mid, t) ((1 ./ (1 + exp(-a .* (t - t_mid))) - off) ./ den) * (right - left) + left; % logistic function that rises from left to right

% construct signal
template = zeros(length(t_vec), length(amps));
prev_amps = zeros(size(amps));
for i = 1:length(edge_times)
    if i == length(edge_times) % this is the last edge
        curr_amps = zeros(size(amps));
    else
        curr_amps = amps .* resp_mat(i,:);
    end
    et = edge_times(i);
    edge_mask = t_vec >= (et - dur_edge/2) & t_vec <= (et + dur_edge/2);
    template(edge_mask,:) = edge_func(prev_amps, curr_amps, et, t_vec(edge_mask));
    rest_mask = t_vec > (et + dur_edge/2);
    template(rest_mask,:) = repmat(curr_amps, [sum(rest_mask), 1]);
    prev_amps = curr_amps;
end

end
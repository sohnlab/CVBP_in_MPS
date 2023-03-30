function [Z0_vec, Z_mat] = simple_pad_array_model(l_vec, w, h, f_vec, dR, rho, cpad)
% inputs:
%   l_vec - size [n_seg, 1] vector of segment lengths [m]
%   w     - channel width [m]
%   h     - channel height [m]
%   f_vec - size [1, n_freq] vector of frequencies [Hz]
%   dR    - deltaR caused by particle [Ohm]
%   rho   - resistivity [Ohm*m]
%   cpad  - specific capacitance of pad [F/m^2]

% outputs:
%   Z0_vec - size [1, n_freq] complex array of baseline impedance
%   Z_mat  - size [n_seg, n_freq] complex array of impedances for when particle is in each segment

n_seg = size(l_vec, 1);
n_freq = size(f_vec, 2);

Z0_vec = sum(rho * l_vec(1:2:end) / (w*h), 1) + sum(( w*h ./ (rho * l_vec(2:2:end)) + 1i * 2*pi*f_vec * cpad * w .* l_vec(2:2:end) ).^-1, 1);

Z_mat = nan(n_seg, n_freq);
for i = 1:n_seg
    if mod(i, 2) == 1 % this is a gap
        Z_mat(i, :) = Z0_vec + dR;
    else % this is a pad
        Z_mat(i, :) = Z0_vec - ( w*h ./ (rho * l_vec(i)) + 1i * 2*pi*f_vec * cpad * w .* l_vec(i) ).^-1 ...
                 + ( 1 ./ (rho * l_vec(i) ./ (w*h) + dR) + 1i * 2*pi*f_vec * cpad * w .* l_vec(i) ).^-1;
    end
end

end
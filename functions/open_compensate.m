function Z_comp_vec = open_compensate(Z_raw_vec, freq_vec, Z_open, freq_open, Z_short, freq_short)
    Z_open_vec = interp_Z(freq_open, Z_open, freq_vec);
    Z_short_vec = interp_Z(freq_short, Z_short, freq_vec);
    Z_comp_vec = (Z_open_vec-Z_short_vec) .* (Z_raw_vec-Z_short_vec) ./ (Z_open_vec-Z_raw_vec);
end

function Z_new = interp_Z(f, Z, f_new)
    Z_new_abs = exp(interp1(log(f), log(abs(Z)), log(f_new), [], 'extrap'));
    Z_new_ang = interp1(log(f), angle(Z), log(f_new), [], 'extrap');
    Z_new = Z_new_abs .* exp(1i .* Z_new_ang);
end

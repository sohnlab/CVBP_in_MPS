function Z_comp_vec = load_compensate(Z_raw_vec, freq_vec, Z_open, freq_open, Z_short, freq_short, Z_load, freq_load, Z_ref)
    Z_open_vec = interp_Z(freq_open, Z_open, freq_vec);
    Z_short_vec = interp_Z(freq_short, Z_short, freq_vec);
    Z_load_vec = interp_Z(freq_load, Z_load, freq_vec);
    Z_comp_vec = Z_ref .* (Z_raw_vec-Z_short_vec) .* (Z_open_vec-Z_load_vec) ...
                     ./ ( (Z_open_vec-Z_raw_vec) .* (Z_load_vec-Z_short_vec) );
end

function Z_new = interp_Z(f, Z, f_new)
    Z_new_abs = exp(interp1(log(f), log(abs(Z)), log(f_new), [], 'extrap'));
    Z_new_ang = interp1(log(f), angle(Z), log(f_new), [], 'extrap');
    Z_new = Z_new_abs .* exp(1i .* Z_new_ang);
end

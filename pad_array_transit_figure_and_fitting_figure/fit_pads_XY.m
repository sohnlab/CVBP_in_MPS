function [Z_none_hat, Z_gaps_hat, Z_pads_hat, X_hat, Y_hat, ab, params_fit, rmse, mae, Linf, rmse_, mae_, Linf_] = fit_pads_XY(Z_none, Z_gaps, Z_pads, f_vec, d_vec, dR_vec, params_init, lb, ub, max_fun_evals)

if nargin == 8, max_fun_evals = 100; end

X_data = ((0.5*abs(Z_gaps(1:end-1,:,:)) + 0.5*abs(Z_gaps(2:end,:,:))) ./ abs(Z_none) - 1) * 100;
Y_data = (abs(Z_pads) ./ abs(Z_none) - 1) * 100;

fun_L2 = @(x) objective_function(x, X_data, Y_data, f_vec, d_vec, dR_vec);
fun_L1 = @(x) sqrt(abs(fun_L2(x)));

[params_fit, resid_norm, fake_resid, exit_flag, out, lam, jac] = lsqnonlin(fun_L2, params_init, lb, ub, optimset('display', 'off', 'MaxFunEvals', max_fun_evals));

[resid, self_resid, Z_none_hat, Z_gaps_hat, Z_pads_hat, X_hat, Y_hat, ab] = objective_function(params_fit, X_data, Y_data, f_vec, d_vec, dR_vec);

rmse = sqrt(mean(resid.^2));
rmse_ = rmse / sqrt(mean(self_resid.^2));
mae = mean(abs(resid));
mae_ = mae / mean(abs(self_resid));
Linf = max(abs(resid));
Linf_ = Linf / max(abs(self_resid));

end

function [resid, self_resid, Z_none_hat, Z_gaps_hat, Z_pads_hat, X_hat, Y_hat, ab] = objective_function(params, X_data, Y_data, f_vec, d_vec, dR_vec)

l_chan = [40, 10, 40, 20, 40, 40, 40, 80, 40, 160, 40, 320, 40, 640, 40] * 1e-6;
pad_mult = params(1);
pad_add = params(2);
l_chan(2:2:end-1) = l_chan(2:2:end-1) * pad_mult + pad_add;
l_extra = sum(l_chan) - 1590e-6;
l_chan(1:2:end) = l_chan(1:2:end) - l_extra / length(l_chan(1:2:end));
if any(l_chan < 0)
    warning('objective_function: l_chan cannot be negative!');
    l_chan(l_chan < 0) = 1e-9;
    l_chan = l_chan ./ sum(l_chan) .* 1590e-6;
    l_chan(l_chan < 0) = 1e-9;
    l_chan = l_chan ./ sum(l_chan) .* 1590e-6;
    l_chan(l_chan < 0) = 1e-9;
    l_chan = l_chan ./ sum(l_chan) .* 1590e-6;
end
w = params(3);
h = params(4);
rho = params(5);
q_vec = nan(1, 23); q_vec([1,6:2:18,23]) = 0;
q_vec(~isnan(q_vec)) = params(6:14);
a_vec = nan(1, 23); a_vec([1,6:2:18,23]) = 0;
a_vec(~isnan(a_vec)) = params(15:23);
dR_mult = params(24);

[Z_none_hat, Z_gaps_hat, Z_pads_hat, X_hat, Y_hat] = pad_array_vol_forward_model(l_chan, w, h, rho, q_vec, a_vec, f_vec, d_vec, dR_vec);

ab = nan(size(X_data,1), size(X_data,2), 2);
fo = fitoptions('Method', 'NonlinearLeastSquares',...
                'Lower', [-Inf, -Inf],...
                'Upper', [Inf, Inf],...
                'StartPoint', [0, 0]);
ft = fittype('a*x^2+b*x', 'options', fo);
for i_p = 1:size(ab,1)
    for i_f = 1:size(ab,2)
        F = fit(squeeze(X_hat(i_p, i_f, :)), squeeze(Y_hat(i_p, i_f, :)), ft);
        ab(i_p, i_f, 1) = F.a; % store fitted params
        ab(i_p, i_f, 2) = F.b; % store fitted params
    end
end

Y_poly = ab(:,:,1) .* X_data.^2 + ab(:,:,2) .* X_data;
Y_poly(X_data < 0) = Y_data(X_data < 0); % to remove influence of negative X values on fit
resid = Y_data(:) - Y_poly(:);
self_resid = Y_data(:);

fprintf('Objective function:\npad_mult = %0.4f, pad_add = %0.0f, l_chan = \n', pad_mult, pad_add*1e6);
disp(round(l_chan*1e6));
fprintf('w = %0.2f, h = %0.2f, rho = %0.4f, dR_mult = %0.4f, q = \n', w*1e6, h*1e6, rho, dR_mult);
disp(params(6:14));
fprintf('a = \n');
disp(params(15:23));
fprintf('X_hat max = %0.3f, Y_hat min = %0.3f\n', max(X_hat, [], "all"), min(Y_hat, [], "all"));
rmse = sqrt(mean(resid.^2));
rmse_ = rmse / sqrt(mean(self_resid.^2));
mae = mean(abs(resid));
mae_ = mae / mean(abs(self_resid));
Linf = max(abs(resid));
Linf_ = Linf / max(abs(self_resid));
fprintf('RMSE = %0.3f (%0.2f%%), MAE = %0.3f (%0.2f%%), Linf = %0.3f (%0.2f%%)\n', rmse, rmse_*100, mae, mae_*100, Linf, Linf_*100);

end
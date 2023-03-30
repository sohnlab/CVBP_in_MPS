function [z_fit, edge_times_fit, params_fit, rmse, mae, Linf, rmse_, mae_, Linf_] = fit_edge_times(t_vec, z_mat, f_vec, params_init, lb, ub)

fun_L2 = @(x) fitting_objective_function(x, t_vec, z_mat, f_vec);
fun_L1 = @(x) sqrt(abs(fun_L2(x)));

[params_fit, resid_norm, fake_resid, exit_flag, out, lam, jac] = lsqnonlin(fun_L1, params_init, lb, ub, optimset('display', 'off'));

P = length(params_init);
K = P - 10;
tc_fit = params_fit(1);
tt_fit = params_fit(2);
fext_fit = params_fit(3);
fmis_fit = params_fit(4);
facc_fit = params_fit(5);
amps_fit = params_fit(6:6+K-1);
w_fit = params_fit(6+K);
h_fit = params_fit(7+K);
dR_fit = params_fit(8+K);
rho_fit = params_fit(9+K);
cpad_fit = params_fit(10+K);

[z_fit, edge_times_fit] = forward_model(t_vec, tc_fit, tt_fit, fext_fit, fmis_fit, facc_fit, amps_fit, w_fit, h_fit, f_vec, dR_fit, rho_fit, cpad_fit);
resid = z_mat - z_fit;
% resid = fake_resid.^2;
rmse = sqrt(mean(mean(resid.^2)));
mae = mean(mean(abs(resid)));
Linf = max(max(abs(resid)));
rmse_ = rmse / sqrt(mean(mean(z_mat.^2)));
mae_ = mae / mean(mean(abs(z_mat)));
Linf_ = Linf / max(max(abs(z_mat)));

end

function res = fitting_objective_function(params, t_vec, z_mat, f_vec)

P = length(params);
K = P - 10; % total number of channels

tc = params(1);
tt = params(2);
fext = params(3);
fmis = params(4);
facc = params(5);
amps = params(6:6+K-1);
w = params(6+K);
h = params(7+K);
dR = params(8+K);
rho = params(9+K);
cpad = params(10+K);

z_hat = forward_model_3(t_vec, tc, tt, fext, fmis, facc, amps, w, h, f_vec, dR, rho, cpad);

res = z_mat - z_hat;

end
close all;
clearvars;
clc;

addpath('../functions');

%%

in_dir = 'bw1000Hz/';
d = dir(in_dir);
fns = {d.name};
mask = contains(fns, '.mat');
fns = fns(mask);

for i = 1:length(fns)
    fn = fns{i};
    comp_dir = '../sweeps/';
    fn_open = '24p07padarr_dev001a_0p2perctween.mat';
    fn_load = '1MOhmdut_par_24p07padarr_dev001a_1MOhmgain_LOAD.mat';
    load([comp_dir fn_open], 'Z_vec', 'freq_vec');
    Z_open = Z_vec;
    f_open = freq_vec;
    load([comp_dir fn_load], 'Z_vec', 'freq_vec');
    Z_load = Z_vec;
    f_load = freq_vec;
    Z_ref = 1e6;
    load([in_dir fn]);
    Z_comp = load_compensate(Z_mat, freq_vec, Z_open, f_open, 0*Z_open, f_open, Z_load, f_load, Z_ref);
    
    save([in_dir fn], 'Z_mat', 'R_ref', 'freq_vec', 'phi_vec', 'fr', 'tr', 'sampleRate', 'Vp_vec', 'bw_nom', 'bw_vec', 'Z_comp', 'Z_open', 'f_open', 'Z_load', 'f_load', '-V7.3');
end

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
    fn_open = 'G60-P80-G60_dev005a_2perctween.mat';
    load([comp_dir fn_open], 'Z_vec', 'freq_vec');
    Z_open = Z_vec;
    f_open = freq_vec;
    load([in_dir fn]);
    Z_comp = open_compensate(Z_mat, freq_vec, Z_open, f_open, 0*Z_open, f_open);
    
   save([out_dir fn], 'Z_mat', 'R_ref', 'freq_vec', 'phi_vec', 'fr', 'tr', 'sampleRate', 'Vp_vec', 'bw_nom', 'bw_vec', 'Z_comp', 'Z_open', 'f_open', '-V7.3');
end

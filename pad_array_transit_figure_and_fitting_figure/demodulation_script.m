close all;
clearvars;
clc;

addpath('../functions');

%% iterate over raw data files

in_dir = 'raw_data/';
d = dir(in_dir);
fns = {d.name};
mask = contains(fns, '.mat');
fns = fns(mask);

bw_nom = 1e3;
out_dir = sprintf('bw%dHz/', bw_nom);

for i_fn = 1:length(fns)
    fn = fns{i_fn};
    if ~contains([in_dir fn], '_run') || exist([out_dir fn '.mat'], 'file')
        continue;
    end
    fprintf('Loading %s\n', fn);
    load([in_dir fn]);
    if ~exist('Vp_vec', 'var')
        Vp_vec = Vp + zeros(size(freq_vec));
    end
    V_TIA = data(:);
    clear data;
    m = mean(V_TIA);
        
    fr = min(bw_nom * 5, sampleRate);
    M = floor(sampleRate / fr); % downsampling factor
    fr = sampleRate / M;
    if fr < 1000
        fr = 1000;
        M = floor(sampleRate / fr);
        fr = sampleRate / M;
    end
    
    V_IQ_vec = nan(1, length(freq_vec));
    V_TIA_IQ_vec = nan(round(length(V_TIA)/sampleRate*fr), length(freq_vec));
    Z_mat = nan(round(length(V_TIA)/sampleRate*fr), length(freq_vec));
    bw_vec = nan(1, length(freq_vec));
    tr = (0:(size(Z_mat, 1)-1)).' ./ fr;
    for i = 1:length(freq_vec)
        bw_vec(i) = min([bw_nom, 0.9*freq_vec(i)]);
        V_IQ_vec(i) = Vp_vec(i) .* exp(1i*phi_vec(i));
        V_TIA_IQ_vec(:,i) = iq_demod(V_TIA, freq_vec(i), sampleRate, fr, bw_vec(i));
        f0 = mod(min([min(freq_vec), sampleRate/2-max(freq_vec)]), 1e3);
        Z_mat(:,i) = -R_ref .* V_IQ_vec(i) ./ comb_filter_keep_DC(V_TIA_IQ_vec(:,i), f0, fr, f0/100); % comb filter
    end

    save([out_dir fn], 'Z_mat', 'R_ref', 'freq_vec', 'phi_vec', 'fr', 'tr', 'sampleRate', 'Vp_vec', 'bw_nom', 'bw_vec', '-V7.3');
end

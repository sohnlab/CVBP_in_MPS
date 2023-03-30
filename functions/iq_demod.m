function [xIQ, tr, fr] = iq_demod(x, f0, fs, fr, bw)

if size(x,1) == 1, x = x.'; end

t = (0:length(x)-1).'/fs; % time vector [s]

mix = 2 * exp(-1i*2*pi*f0*t);
clear t;

xIQ = x .* mix;
clear x;
clear mix;

M = floor(fs / fr); % downsampling factor
fr = fs / M;
if fr < 1000
    fr = 1000;
    M = floor(fs / fr);
    fr = fs / M;
end

m = M;
while mod(m, 10)==0
    xIQ = resample(xIQ - xIQ(1,:), 1, 10) + xIQ(1,:);
    m = m / 10;
end
if m > 1
    xIQ = resample(xIQ - xIQ(1,:), 1, m) + xIQ(1,:);
end

xIQ = lowpass(xIQ - xIQ(1,:), bw, fr, Steepness=0.99) + xIQ(1,:);
tr = (0:length(xIQ)-1).' ./ fr;

end

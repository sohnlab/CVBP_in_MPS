function width = channel_width_func(x, x_l, l_vec, w_vec)
w_vec = [w_vec, Inf];
width = w_vec(find(x > [x_l, x_l + cumsum(l_vec)], 1, 'last'));
if isempty(width), width = Inf; end
end
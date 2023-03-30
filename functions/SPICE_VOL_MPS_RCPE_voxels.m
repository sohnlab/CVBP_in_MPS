function [Z, result, dRoR_o_dVoV_fit] = SPICE_VOL_MPS_RCPE_voxels(l_vec, w_vec, h_vec, node_num_vec, rho, q_vec, a_vec, d_cell_vec, x_cell_vec, dR, f_vec, dxyz, save_all_vars, dRoR_o_dVoV_init)
% NOTE: in this function, the voltage electrodes are considered pads

%% handle inputs

if nargin < 12
    dxyz = 1e-6; % default deltax, deltay, deltaz
end
if isscalar(dxyz)
    dx = dxyz;
    dy = dxyz;
    dz = dxyz;
elseif length(dxyz)==2
    dx = dxyz(1);
    dy = dxyz(2);
    dz = dy;
elseif length(dxyz)==3
    dx = dxyz(1);
    dy = dxyz(2);
    dz = dxyz(3);
else
    error('SPICE_VOL_MPS_RCPE_voxels: invalid input for dxyz');
end
if nargin < 13
    save_all_vars = false;
end
n_segs = length(l_vec); % number of different channel segments
if isempty(node_num_vec)
    node_num_vec = nan(1, n_segs);
    node_num_vec(1) = 1;
    node_num_vec(end) = 0;
    fprintf('SPICE_VOL_MPS_RCPE_voxels: node_num_vec is empty! assuming no metal pads in channel other than electrodes (NPS)!\n');
end
if ~isscalar(w_vec), error('SPICE_VOL_MPS_RCPE_voxels: different widths not currently supported!'); end
if isscalar(w_vec), w_vec = repmat(w_vec, size(l_vec)); end
if ~isscalar(h_vec), error('SPICE_VOL_MPS_RCPE_voxels: different heights not currently supported!'); end
if isscalar(h_vec), h_vec = repmat(h_vec, size(l_vec)); end
if isscalar(q_vec), q_vec = repmat(q_vec, size(l_vec)); end
if isscalar(a_vec), a_vec = repmat(a_vec, size(l_vec)); end
x_vec_seg_edges = [0, cumsum(l_vec)]; % vector of segment edge x locations

% n_blocks_x = max([1, round(sum(l_vec)/dx)]);
n_blocks_y = max([1, round(h_vec(1)/dy)]);
n_blocks_z = max([1, round(w_vec(1)/dz)]);
n_steps_vec = round(l_vec ./ dx); n_steps_vec(n_steps_vec<1) = 1; % n_steps per segment
ind_blocks = []; % which segment index does each block correspond to?
for i = 1:n_segs
    ind_blocks = [ind_blocks, repmat(i, [1, n_steps_vec(i)])];
end

%% set paths and open netlist file

if ispc
    spicePath = 'C:\Program Files\LTC\LTspiceXVII\XVIIx64.exe'; % This is the path to your LT Spice installation
    filePath = 'C:\Users\alan_dong\Desktop\'; % This is the path where you want the netlist, functions and simulation outputs to reside
    netPath = strrep(filePath, '\', '\\'); % This is same as file path but all \ are replaced by \\ because of some Matlab issue in writing Netlist.
elseif isunix
    spicePath = 'XVIIx64.exe'; % This is the path to your LT Spice installation
    filePath = '/home/alan/.wine/drive_c/Program Files/LTC/LTspiceXVII/'; % This is the path where you want the netlist, functions and simulation outputs to reside
    netPath = filePath;
end
fileName = tempname('/');
fileName = fileName(2:end);

netName = sprintf('%s%s.net', netPath, fileName);
if exist(netName, 'file')
    delete(netName);
end
fileID = fopen(netName, 'w');
fprintf(fileID, '* %s.asc\n', fileName); % print path to schematic

%% calculate cell shape (ellipsoid [a, b, c] a.k.a. [x, y, z] diameters)

% first see if there's a cell or not
no_cell = true;
if ~isempty(x_cell_vec) && ~isempty(d_cell_vec)
    % check if cell is outside channel
    for i_c = 1:length(d_cell_vec) % iterate over cells
        d_cell = d_cell_vec(i_c);
        x_cell = x_cell_vec(i_c);
        [a_cell, b_cell, c_cell] = find_abc(x_cell, d_cell, 0, l_vec, w_vec, h_vec(1));
        a_cell_vec(i_c) = a_cell;
        b_cell_vec(i_c) = b_cell;
        c_cell_vec(i_c) = c_cell;
        x_cell_left = x_cell - a_cell/2;
        x_cell_right = x_cell + a_cell/2;
        x_cell_left_vec(i_c) = x_cell_left;
        mask_cell = x_cell_left < x_vec_seg_edges(2:end) & x_cell_right > x_vec_seg_edges(1:end-1); % which segments cell is intersecting
        if any(mask_cell)
            no_cell = false;
        end
    end
end

%% separate each pad segment into n_steps blocks

% compute length of each block, also width and height
l_vec_blocks = l_vec(ind_blocks) ./ n_steps_vec(ind_blocks);
w_vec_blocks = w_vec(ind_blocks);
h_vec_blocks = h_vec(ind_blocks) / n_blocks_y;

node_num_vec_blocks_bottom = node_num_vec(ind_blocks);
q_vec_blocks = q_vec(ind_blocks);
a_vec_blocks = a_vec(ind_blocks);

V_vec_blocks = l_vec_blocks .* w_vec_blocks .* h_vec_blocks; % volume of each block
V_mat_quarterblocks = kron(repmat(V_vec_blocks, [n_blocks_y,1]), ones(2,2)/4);

% compute partial volume of cells in each block and assign dR per block
x_grid_vec = [0, cumsum(l_vec_blocks)];
y_grid_vec = linspace(h_vec(1), 0, n_blocks_y+1).'; % layers
z_grid_vec = reshape(linspace(0, w_vec(1), n_blocks_z+1), 1, 1, []); % depths
% compare min feature size with voxel size, want at least 2 blocks in each dimension of channel/cells
min_dim_x = min([sum(l_vec), a_cell_vec]);
min_dim_y = min([w_vec, c_cell_vec]);
min_dim_z = min([h_vec, b_cell_vec]);
subdiv = ceil(2*max([dx/min_dim_x, dy/min_dim_y, dz/min_dim_z]));
if ~no_cell
    dV_mat_quarterblocks = zeros((numel(y_grid_vec)-1)*2, (numel(x_grid_vec)-1)*2);
    for i_c = 1:length(d_cell_vec) % iterate over cells
        a_cell = a_cell_vec(i_c);
        b_cell = b_cell_vec(i_c);
        c_cell = c_cell_vec(i_c);
        dV_mat_quarterblocks_i = compute_dV_per_quarterblock(x_grid_vec, y_grid_vec, z_grid_vec, x_cell, h_vec(1)/2, w_vec(1)/2, a_cell, b_cell, c_cell, subdiv);
        dV_mat_quarterblocks = dV_mat_quarterblocks + dV_mat_quarterblocks_i;
    end
    dVoV_quarterblocks = dV_mat_quarterblocks ./ V_mat_quarterblocks;
end

%% create grid points and assign node_num

node_num_curr = max(node_num_vec) + 1;
node_num_grid_mat = nan(length(y_grid_vec) + 1, length(x_grid_vec)); % bottom row is for pads' node_num
for r = 1:(length(y_grid_vec)+1)
    for c = 1:length(x_grid_vec)
        if isnan(node_num_vec(1)) && c == 1 % connect leftmost gap to node 1
            node_num_grid_mat(r,c) = 1;
        elseif isnan(node_num_vec(end)) && c == length(x_grid_vec) % connect rightmost gap to node 0
            node_num_grid_mat(r,c) = 0;
        elseif r == (length(y_grid_vec)+1) % bottom metal layer
            if c == length(x_grid_vec) || isnan(node_num_vec_blocks_bottom(c)) % right edge of pad
                node_num_grid_mat(r,c) = node_num_vec_blocks_bottom(c-1);
            else
                node_num_grid_mat(r,c) = node_num_vec_blocks_bottom(c);
            end
        else
            node_num_grid_mat(r,c) = node_num_curr;
            node_num_curr = node_num_curr + 1;
        end
    end
end

%% if there is a cell, solve for the appropriate dRoR_o_dVoV to use to get to the target dR

if ~no_cell
    resistance_func = @(dRoR_o_dVoV) compute_DC_resistance(dRoR_o_dVoV, netPath, spicePath, filePath, x_grid_vec, y_grid_vec, node_num_grid_mat, w_vec_blocks, h_vec_blocks, rho, dV_mat_quarterblocks, V_mat_quarterblocks);
%     R_baseline = resistance_func(0); % baseline resistance
    R_baseline = rho * sum(l_vec ./ w_vec ./ h_vec);
    obj_func = @(dRoR_o_dVoV) (R_baseline + dR) - resistance_func(dRoR_o_dVoV);
    if nargin < 14
        dRoR_o_dVoV_init = dR / R_baseline * sum(V_mat_quarterblocks, "all") / sum(dV_mat_quarterblocks, "all"); % guess a value, this will almost certainly fall short of the target dR
    end
    dRoR_o_dVoV_fit = my_fzero(obj_func, 0, dR, dRoR_o_dVoV_init);
else
    dRoR_o_dVoV_fit = [];
end

%% construct components and print to netlist file

Iin = 1; % Iin [A]
fprintf(fileID, 'Iin %d %d AC %d\n', 0, 1, Iin); % current source always between nodes 0 and 1

for r = 1:length(y_grid_vec)
    for c = 1:length(x_grid_vec)
        if r == 1 || r == length(y_grid_vec)
            mult_horiz = 2;
        else
            mult_horiz = 1;
        end
        if c == 1 || c == length(x_grid_vec)
            mult_vert = 2;
        else
            mult_vert = 1;
        end
        if c < length(x_grid_vec) % connect to right node with resistor
            R = rho * (x_grid_vec(c+1) - x_grid_vec(c)) / (w_vec_blocks(c) * h_vec_blocks(c)) * mult_horiz;
            if ~no_cell
                if r == 1 % top row
                    r_quarterblocks = 2*(r-1) + 1;
                    c_quarterblocks = 2*(c-1) + (1:2);
                    dV_i = sum(dV_mat_quarterblocks(r_quarterblocks, c_quarterblocks), "all");
                    V_i = sum(V_mat_quarterblocks(r_quarterblocks, c_quarterblocks), "all");
                    R = R * (1 + dV_i/V_i * dRoR_o_dVoV_fit);
                elseif r == length(y_grid_vec) % bottom row
                    r_quarterblocks = 2*(r-1) + 0;
                    c_quarterblocks = 2*(c-1) + (1:2);
                    dV_i = sum(dV_mat_quarterblocks(r_quarterblocks, c_quarterblocks), "all");
                    V_i = sum(V_mat_quarterblocks(r_quarterblocks, c_quarterblocks), "all");
                    R = R * (1 + dV_i/V_i * dRoR_o_dVoV_fit);
                else
                    r_quarterblocks = 2*(r-1) + (0:1);
                    c_quarterblocks = 2*(c-1) + (1:2);
                    dV_i = sum(dV_mat_quarterblocks(r_quarterblocks, c_quarterblocks), "all");
                    V_i = sum(V_mat_quarterblocks(r_quarterblocks, c_quarterblocks), "all");
                    R = R * (1 + dV_i/V_i * dRoR_o_dVoV_fit);
                end
            end
            fprintf(fileID, 'R_%d_%d %d %d %0.9g\n', node_num_grid_mat(r,c), node_num_grid_mat(r,c+1), node_num_grid_mat(r,c), node_num_grid_mat(r,c+1), R);
        end
        if r < length(y_grid_vec) && node_num_grid_mat(r,c) ~= node_num_grid_mat(r+1,c) % connect to lower node with resistor
            if c == 1
                w_avg = w_vec_blocks(c);
            elseif c == length(x_grid_vec)
                w_avg = w_vec_blocks(c-1);
            else
                w_avg = (w_vec_blocks(c-1) + w_vec_blocks(c))/2;
            end
            R = rho * (y_grid_vec(r) - y_grid_vec(r+1)) / (w_avg * h_vec_blocks(c)) * mult_vert;
            if ~no_cell
                if c == 1 % leftmost column
                    r_quarterblocks = 2*(r-1) + (1:2);
                    c_quarterblocks = 2*(c-1) + 1;
                    dV_i = sum(dV_mat_quarterblocks(r_quarterblocks, c_quarterblocks), "all");
                    V_i = sum(V_mat_quarterblocks(r_quarterblocks, c_quarterblocks), "all");
                    R = R * (1 + dV_i/V_i * dRoR_o_dVoV);
                elseif c == length(x_grid_vec) % rightmost column
                    r_quarterblocks = 2*(r-1) + (1:2);
                    c_quarterblocks = 2*(c-1) + 0;
                    dV_i = sum(dV_mat_quarterblocks(r_quarterblocks, c_quarterblocks), "all");
                    V_i = sum(V_mat_quarterblocks(r_quarterblocks, c_quarterblocks), "all");
                    R = R * (1 + dV_i/V_i * dRoR_o_dVoV_fit);
                else
                    r_quarterblocks = 2*(r-1) + (1:2);
                    c_quarterblocks = 2*(c-1) + (0:1);
                    dV_i = sum(dV_mat_quarterblocks(r_quarterblocks, c_quarterblocks), "all");
                    V_i = sum(V_mat_quarterblocks(r_quarterblocks, c_quarterblocks), "all");
                    R = R * (1 + dV_i/V_i * dRoR_o_dVoV_fit);
                end
            end
            fprintf(fileID, 'R_%d_%d %d %d %0.9g\n', node_num_grid_mat(r,c), node_num_grid_mat(r+1,c), node_num_grid_mat(r,c), node_num_grid_mat(r+1,c), R);
        elseif ~isnan(node_num_grid_mat(r+1,c)) && node_num_grid_mat(r,c) ~= node_num_grid_mat(r+1,c) % connect to metal layer node with CPE
            if c == length(x_grid_vec) || node_num_grid_mat(r+1,c) ~= node_num_grid_mat(r+1,c+1) % right pad edge
                q = q_vec_blocks(c-1);
                a = a_vec_blocks(c-1);
            else
                q = q_vec_blocks(c);
                a = a_vec_blocks(c);
            end
            if c == 1 || node_num_grid_mat(r+1,c-1) ~= node_num_grid_mat(r+1,c) % left pad edge
                Q = q * (x_grid_vec(c+1) - x_grid_vec(c))/2 * w_vec_blocks(c);
            elseif c == length(x_grid_vec) || node_num_grid_mat(r+1,c) ~= node_num_grid_mat(r+1,c+1) % right pad edge
                Q = q * (x_grid_vec(c) - x_grid_vec(c-1))/2 * w_vec_blocks(c-1);
            else
                Q = q * (x_grid_vec(c) - x_grid_vec(c-1)) * w_vec_blocks(c-1);
            end
            fprintf(fileID, 'R_%d_%d %d %d R=1 Laplace=%0.9g*s^%0.9g\n', node_num_grid_mat(r,c), node_num_grid_mat(r+1,c), node_num_grid_mat(r,c), node_num_grid_mat(r+1,c), Q, a);
        end
    end
end

%% print ac analysis string and end strings and close file

fprintf(fileID, '.ac list ');
fprintf(fileID, '%g ', f_vec);
fprintf(fileID, '\n');
fprintf(fileID, '.backanno\n.end\n');
fclose(fileID);

%% run analysis

if save_all_vars
    result = simulateModel(spicePath, fileName, filePath, 'all');
    Z = result.variable_mat(strcmp(result.variable_name_list, 'V(1)'), :);
    V_grid_mat = zeros(size(node_num_grid_mat,1), size(node_num_grid_mat,2), length(f_vec));
    for r = 1:size(node_num_grid_mat, 1)
        for c = 1:size(node_num_grid_mat, 2)
            node_num = node_num_grid_mat(r,c);
            if node_num > 0
                V_grid_mat(r,c,:) = result.variable_mat(strcmp(result.variable_name_list, ['V(' num2str(node_num) ')']),:);
            end
        end
    end
    result.x_grid_vec = x_grid_vec;
    result.y_grid_vec = y_grid_vec;
    result.node_num_grid_mat = node_num_grid_mat;
    result.V_grid_mat = V_grid_mat;
else
    result = simulateModel(spicePath, fileName, filePath, 1);
    Z = result.variable_mat;
end

end

%% functions

function R = compute_DC_resistance(dRoR_o_dVoV, netPath, spicePath, filePath, x_grid_vec, y_grid_vec, node_num_grid_mat, w_vec_blocks, h_vec_blocks, rho, dV_mat_quarterblocks, V_mat_quarterblocks)

fileName = tempname('/');
fileName = fileName(2:end);

netName = sprintf('%s%s.net', netPath, fileName);
if exist(netName, 'file')
    delete(netName);
end
fileID = fopen(netName, 'w');
fprintf(fileID, '* %s.asc\n', fileName); % print path to schematic

Iin = 1; % Iin [A]
fprintf(fileID, 'Iin %d %d DC %d\n', 0, 1, Iin); % current source always between nodes 0 and 1

for r = 1:length(y_grid_vec)
    for c = 1:length(x_grid_vec)
        if r == 1 || r == length(y_grid_vec)
            mult_horiz = 2;
        else
            mult_horiz = 1;
        end
        if c == 1 || c == length(x_grid_vec)
            mult_vert = 2;
        else
            mult_vert = 1;
        end
        if c < length(x_grid_vec) % connect to right node with resistor
            R = rho * (x_grid_vec(c+1) - x_grid_vec(c)) / (w_vec_blocks(c) * h_vec_blocks(c)) * mult_horiz;
            if r == 1 % top row
                r_quarterblocks = 2*(r-1) + 1;
                c_quarterblocks = 2*(c-1) + (1:2);
                dV_i = sum(dV_mat_quarterblocks(r_quarterblocks, c_quarterblocks), "all");
                V_i = sum(V_mat_quarterblocks(r_quarterblocks, c_quarterblocks), "all");
                R = R * (1 + dV_i/V_i * dRoR_o_dVoV);
            elseif r == length(y_grid_vec) % bottom row
                r_quarterblocks = 2*(r-1) + 0;
                c_quarterblocks = 2*(c-1) + (1:2);
                dV_i = sum(dV_mat_quarterblocks(r_quarterblocks, c_quarterblocks), "all");
                V_i = sum(V_mat_quarterblocks(r_quarterblocks, c_quarterblocks), "all");
                R = R * (1 + dV_i/V_i * dRoR_o_dVoV);
            else
                r_quarterblocks = 2*(r-1) + (0:1);
                c_quarterblocks = 2*(c-1) + (1:2);
                dV_i = sum(dV_mat_quarterblocks(r_quarterblocks, c_quarterblocks), "all");
                V_i = sum(V_mat_quarterblocks(r_quarterblocks, c_quarterblocks), "all");
                R = R * (1 + dV_i/V_i * dRoR_o_dVoV);
            end
            fprintf(fileID, 'R_%d_%d %d %d %0.9g\n', node_num_grid_mat(r,c), node_num_grid_mat(r,c+1), node_num_grid_mat(r,c), node_num_grid_mat(r,c+1), R);
        end
        if r < length(y_grid_vec) && node_num_grid_mat(r,c) ~= node_num_grid_mat(r+1,c) % connect to lower node with resistor
            if c == 1
                w_avg = w_vec_blocks(c);
            elseif c == length(x_grid_vec)
                w_avg = w_vec_blocks(c-1);
            else
                w_avg = (w_vec_blocks(c-1) + w_vec_blocks(c))/2;
            end
            R = rho * (y_grid_vec(r) - y_grid_vec(r+1)) / (w_avg * h_vec_blocks(c)) * mult_vert;
            if c == 1 % leftmost column
                r_quarterblocks = 2*(r-1) + (1:2);
                c_quarterblocks = 2*(c-1) + 1;
                dV_i = sum(dV_mat_quarterblocks(r_quarterblocks, c_quarterblocks), "all");
                V_i = sum(V_mat_quarterblocks(r_quarterblocks, c_quarterblocks), "all");
                R = R * (1 + dV_i/V_i * dRoR_o_dVoV);
            elseif c == length(x_grid_vec) % rightmost column
                r_quarterblocks = 2*(r-1) + (1:2);
                c_quarterblocks = 2*(c-1) + 0;
                dV_i = sum(dV_mat_quarterblocks(r_quarterblocks, c_quarterblocks), "all");
                V_i = sum(V_mat_quarterblocks(r_quarterblocks, c_quarterblocks), "all");
                R = R * (1 + dV_i/V_i * dRoR_o_dVoV);
            else
                r_quarterblocks = 2*(r-1) + (1:2);
                c_quarterblocks = 2*(c-1) + (0:1);
                dV_i = sum(dV_mat_quarterblocks(r_quarterblocks, c_quarterblocks), "all");
                V_i = sum(V_mat_quarterblocks(r_quarterblocks, c_quarterblocks), "all");
                R = R * (1 + dV_i/V_i * dRoR_o_dVoV);
            end
            fprintf(fileID, 'R_%d_%d %d %d %0.9g\n', node_num_grid_mat(r,c), node_num_grid_mat(r+1,c), node_num_grid_mat(r,c), node_num_grid_mat(r+1,c), R);
        end
    end
end

% print dc analysis string and end strings and close file
fprintf(fileID, '.op\n');
fprintf(fileID, '.backanno\n.end\n');
fclose(fileID);

result = simulateModel(spicePath, fileName, filePath, 1);
R = result.variable_mat(1);

fprintf('dRoR_o_dVoV = %0.9g -> R = %0.9g\n', dRoR_o_dVoV, R);

end

%%

function dV = compute_dV_per_quarterblock(x_grid_vec, y_grid_vec, z_grid_vec, x_cell, y_cell, z_cell, a_cell, b_cell, c_cell, subdiv)

if nargin < 10
    subdiv = 512;
end

x_grid_vec_new = [reshape((reshape(x_grid_vec(1:end-1), [], 1) + diff(x_grid_vec(:)).*(0:subdiv*2-1)./(subdiv*2)).', [], 1); x_grid_vec(end)].';
y_grid_vec_new = [reshape((reshape(y_grid_vec(1:end-1), [], 1) + diff(y_grid_vec(:)).*(0:subdiv*2-1)./(subdiv*2)).', [], 1); y_grid_vec(end)];
z_grid_vec_new = reshape([reshape((reshape(z_grid_vec(1:end-1), [], 1) + diff(z_grid_vec(:)).*(0:subdiv*2-1)./(subdiv*2)).', [], 1); z_grid_vec(end)], 1, 1, []);

[part_vol_subdiv, ~] = compute_volume_fraction(x_grid_vec_new, y_grid_vec_new, z_grid_vec_new, x_cell, y_cell, z_cell, a_cell, b_cell, c_cell);

part_vol_flat = sum(part_vol_subdiv, 3);
part_vol = 0;
for i = 0:subdiv-1
    for j = 0:subdiv-1
        part_vol = part_vol + part_vol_flat((1:subdiv:end-subdiv+1) + i, (1:subdiv:end-subdiv+1) + j);
    end
end

dV = sum(part_vol, 3);

end

%%

function x_zero = my_fzero(func, x0, y0, x1, y_tol_mult, max_eval)

if nargin < 6, max_eval = 5; end
if nargin < 5, y_tol_mult = 1e-3; end
y_tol = y_tol_mult .* y0;
y1 = func(x1);
if abs(y1) <= y_tol, x_zero = x1; return; end
x2 = interp1([y0, y1], [x0, x1], 0, "linear", "extrap");
y2 = func(x2);
if abs(y2) <= y_tol, x_zero = x2; return; end
x_vec = [x0, x1, x2];
y_vec = [y0, y1, y2];
while all(abs(y_vec) > y_tol) && length(y_vec) < max_eval
    pi = polyfit(x_vec(end-2:end), y_vec(end-2:end), 2);
    ri = roots(pi);
    [~, ind_closer] = min(abs(ri-x0));
    x_vec = [x_vec, ri(ind_closer)];
    y_vec = [y_vec, func(x_vec(end))];
end
[~, ind_best] = min(abs(y_vec));
x_zero = x_vec(ind_best);

end

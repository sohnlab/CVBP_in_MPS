function [Z, result] = SPICE_VOL_MPS_RCPE_slices(l_vec, w_vec, h_vec, node_num_vec, rho, q_vec, a_vec, d_vec, x_vec, dR, f_vec, n_steps)
% NOTE: in this function, the voltage electrodes are considered pads

%% handle inputs
if nargin == 11
    n_steps = 250; % default n_steps
end
n_segs = length(l_vec); % number of different channel segments
if isempty(node_num_vec)
    node_num_vec = nan(1, n_segs);
    node_num_vec(1) = 1;
    node_num_vec(end) = 0;
    fprintf('generate_MPS_netlist: node_num_vec is empty! assuming no metal pads in channel other than electrodes (NPS)!\n');
end
if isscalar(w_vec), w_vec = repmat(w_vec, size(l_vec)); end
if isscalar(h_vec), h_vec = repmat(h_vec, size(l_vec)); end
if isscalar(q_vec), q_vec = repmat(q_vec, size(l_vec)); end
if isscalar(a_vec), a_vec = repmat(a_vec, size(l_vec)); end
x_vec_seg_edges = [0, cumsum(l_vec)]; % vector of segment edge x locations

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
if ~isempty(x_vec) && ~isempty(d_vec)
    % check if cell is outside channel
    for i_c = 1:length(d_vec) % iterate over cells
        d_cell = d_vec(i_c);
        x_cell = x_vec(i_c);
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
mask_ispad = ~isnan(node_num_vec); % which segments are pads
inds_ispad = find(mask_ispad); % indices of pad segments
ind_blocks = sort([1:length(node_num_vec), repmat(inds_ispad, [1, n_steps-1])], "ascend"); % which segment index does each block correspond to?
n_pads = length(inds_ispad); % number of pads
n_gaps = n_segs - n_pads; % number of gaps
n_blocks = n_gaps + n_pads * n_steps; % total number of blocks

% create vector of center node number and bottom node number for each block
node_num_vec_blocks_bottom = node_num_vec(ind_blocks); % node number of bottom of each block
if isnan(node_num_vec(1)) && ~isnan(node_num_vec(end)) % leftmost segment is gap, first slice should be node 1
    node_num_vec_slices = [1, max([1, max(node_num_vec)]) + (1:n_blocks)]; % node number of each slice
elseif ~isnan(node_num_vec(1)) && isnan(node_num_vec(end)) % rightmost segment is gap, last slice should be node 0
    node_num_vec_slices = [max([1, max(node_num_vec)]) + (1:n_blocks), 0]; % node number of each slice
elseif isnan(node_num_vec(1)) && isnan(node_num_vec(end)) % both leftmost segment and rightmost segment are gap
    node_num_vec_slices = [1, max([1, max(node_num_vec)]) + (1:n_blocks-1), 0]; % node number of each slice
else
    node_num_vec_slices = max([1, max(node_num_vec)]) + (1:n_blocks+1); % node number of each slice
end

q_vec_blocks = q_vec(ind_blocks);
a_vec_blocks = a_vec(ind_blocks);

% compute length of each block, also width and height
l_vec_blocks = l_vec;
l_vec_blocks(mask_ispad) = l_vec_blocks(mask_ispad) / n_steps;
l_vec_blocks = l_vec_blocks(ind_blocks);
w_vec_blocks = w_vec(ind_blocks);
h_vec_blocks = h_vec(ind_blocks);
V_vec_blocks = l_vec_blocks .* w_vec_blocks .* h_vec_blocks; % volume of each block

% compute resistance and CPE parameter of each halfblock and/or block
R_vec_blocks = rho .* l_vec_blocks ./ (w_vec_blocks .* h_vec_blocks);
if ~no_cell
    x_vec_block_edges = [0, cumsum(l_vec_blocks)];
    V_ell_vec_blocks = zeros(size(V_vec_blocks));
    for i_c = 1:length(d_vec) % iterate over cells
        a_cell = a_cell_vec(i_c);
        b_cell = b_cell_vec(i_c);
        c_cell = c_cell_vec(i_c);
        h_ell_vec_blocks_i = x_vec_block_edges(2:end) - x_cell_left_vec(i_c);
        ind_last = find(h_ell_vec_blocks_i > a_cell, 1, "first");
        h_ell_vec_blocks_i(ind_last) = a_cell;
        h_ell_vec_blocks_i(ind_last+1:end) = 0;
        h_ell_vec_blocks_i(h_ell_vec_blocks_i <= 0) = 0; % exclude blocks that don't intersect the ellipsoid
        V_ell_vec_blocks_i = diff([0, pi*b_cell*c_cell/(3*a_cell^2).*h_ell_vec_blocks_i.^2.*(3/2*a_cell-h_ell_vec_blocks_i)]); % volume of left half ellipsoid contained in each halfblock
        V_ell_vec_blocks_i(V_ell_vec_blocks_i < 0) = 0; % volume of cell contained in each halfblock
        V_ell_vec_blocks = V_ell_vec_blocks + V_ell_vec_blocks_i; % add to total volume of all cells contained in each halfblock
    end
    dR_vec_blocks = V_ell_vec_blocks ./ V_vec_blocks;
    dR_vec_blocks = dR_vec_blocks ./ sum(dR_vec_blocks) .* dR;
    R_vec_blocks = R_vec_blocks + dR_vec_blocks;
end
Q_vec_blocks = q_vec_blocks .* l_vec_blocks .* w_vec_blocks; % CPE parameter between each center node and its bottom node

%% construct components and print to netlist file

Iin = 1; % Iin [A]
fprintf(fileID, 'Iin %d %d AC %d\n', 0, 1, Iin); % current source always between nodes 0 and 1

for i = 1:length(node_num_vec_slices) % iterate through the slices
    if (i<length(node_num_vec_slices) && ~isnan(node_num_vec_blocks_bottom(i))) || (i>1 && ~isnan(node_num_vec_blocks_bottom(i-1))) % pad present
        if i == 1 || (isnan(node_num_vec_blocks_bottom(i-1)) && ~isnan(node_num_vec_blocks_bottom(i))) % left edge of pad
            fprintf(fileID, 'R_%d_%d %d %d R=1 Laplace=%0.9g*s^%0.9g\n', node_num_vec_slices(i), node_num_vec_blocks_bottom(i), node_num_vec_slices(i), node_num_vec_blocks_bottom(i), 0.5 * Q_vec_blocks(i), a_vec_blocks(i));
        elseif i == length(node_num_vec_slices) || (~isnan(node_num_vec_blocks_bottom(i-1)) && isnan(node_num_vec_blocks_bottom(i))) % right edge of pad
            fprintf(fileID, 'R_%d_%d %d %d R=1 Laplace=%0.9g*s^%0.9g\n', node_num_vec_slices(i), node_num_vec_blocks_bottom(i-1), node_num_vec_slices(i), node_num_vec_blocks_bottom(i-1), 0.5 * Q_vec_blocks(i-1), a_vec_blocks(i-1));
        elseif ~isnan(node_num_vec_blocks_bottom(i-1)) && ~isnan(node_num_vec_blocks_bottom(i)) % inside of pad
            fprintf(fileID, 'R_%d_%d %d %d R=1 Laplace=%0.9g*s^%0.9g\n', node_num_vec_slices(i), node_num_vec_blocks_bottom(i), node_num_vec_slices(i), node_num_vec_blocks_bottom(i), 0.5 * Q_vec_blocks(i) + 0.5 * Q_vec_blocks(i-1), a_vec_blocks(i));
        end
    end
    if i < length(node_num_vec_slices) % not last slice, connect to node to the right with a resistor
        fprintf(fileID, 'R_%d_%d %d %d %0.9g\n', node_num_vec_slices(i), node_num_vec_slices(i+1), node_num_vec_slices(i), node_num_vec_slices(i+1), R_vec_blocks(i));
    end
end

%% print ac analysis string and end strings and close file
fprintf(fileID, '.ac list ');
fprintf(fileID, '%g ', f_vec);
fprintf(fileID, '\n');
fprintf(fileID, '.backanno\n.end\n');
fclose('all');

%% run analysis
result = simulateModel(spicePath, fileName, filePath, 1);
Z = result.variable_mat;

end

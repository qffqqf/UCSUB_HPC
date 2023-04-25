function [coords_disc, nodes_disc, idx_UCs, nodes_disc_BC, nodes_disc_FR] = get_disconnected_nodes(N_UC_x, N_UC_y, UC_nodes, mesh_data, uc_size, forc_loc, resp_loc, checkplots)

%% Determine real coordinates of the N_UC_x*N_UC_y disconnected UCs
% Coordinates per disconnected UC
n_nodes_UC = size(mesh_data.nd_data,1);
[yy,xx] = meshgrid(1:N_UC_y,1:N_UC_x);
coords_disc.UC = cell(1,N_UC_x*N_UC_y);
L_UC_x = uc_size(1);
L_UC_y = uc_size(2);
for i = 1:N_UC_x*N_UC_y
    local_ids = 1:n_nodes_UC;
    coords_disc.UC{i}(:,1) = local_ids + (i-1)*n_nodes_UC;
    coords_disc.UC{i}(:,2) = mesh_data.nd_data(:,1) + (xx(i)-1)*L_UC_x;
    coords_disc.UC{i}(:,3) = mesh_data.nd_data(:,2) + (yy(i)-1)*L_UC_y;
    coords_disc.UC{i}(:,4) = mesh_data.nd_data(:,3);
end
% Coordinates of all disconnected UCs
coords_disc.all = cat(1, coords_disc.UC{:});

%% Determine all original nodes of the N_UC_x*N_UC_y disconnected UCs
% Disconnected nodes per UC
nodes_disc.UC = cell(1,N_UC_x*N_UC_y);
for i = 1:N_UC_x*N_UC_y
    nodes_disc.UC{i}.I = UC_nodes.I + (i-1)*n_nodes_UC;
    nodes_disc.UC{i}.L = UC_nodes.L + (i-1)*n_nodes_UC;
    nodes_disc.UC{i}.R = UC_nodes.R + (i-1)*n_nodes_UC;
    nodes_disc.UC{i}.B = UC_nodes.B + (i-1)*n_nodes_UC;
    nodes_disc.UC{i}.T = UC_nodes.T + (i-1)*n_nodes_UC;
    nodes_disc.UC{i}.BL = UC_nodes.BL + (i-1)*n_nodes_UC;
    nodes_disc.UC{i}.BR = UC_nodes.BR + (i-1)*n_nodes_UC;
    nodes_disc.UC{i}.TR = UC_nodes.TR + (i-1)*n_nodes_UC;
    nodes_disc.UC{i}.TL = UC_nodes.TL + (i-1)*n_nodes_UC;
    nodes_disc.UC{i}.all = sort([nodes_disc.UC{i}.I, nodes_disc.UC{i}.L, ...
                                 nodes_disc.UC{i}.R, nodes_disc.UC{i}.B, nodes_disc.UC{i}.T, ...
                                 nodes_disc.UC{i}.BL, nodes_disc.UC{i}.BR, nodes_disc.UC{i}.TR, ...
                                 nodes_disc.UC{i}.TL]);
end

% All disconnected nodes
nodes_disc.all = [1:N_UC_x*N_UC_y*n_nodes_UC];
% All disconnected interior nodes (for visualization)
nodes_disc.interior = [];
for i = 1:N_UC_x*N_UC_y
    nodes_disc.interior = [nodes_disc.interior, nodes_disc.UC{i}.I];
end
% All disconnected interface nodes (for visualization)
nodes_disc.interface = setdiff(nodes_disc.all, nodes_disc.interior);

%% Boundary nodes
idx_UCs = reshape([1:N_UC_x*N_UC_y], N_UC_x,N_UC_y);
idx_UCs_B = idx_UCs(:,1);
idx_UCs_T = idx_UCs(:,end);
idx_UCs_L = idx_UCs(1,:);
idx_UCs_R = idx_UCs(end,:);

nodes_disc_BC.B = [];
for i = 1:length(idx_UCs_B(:))
    nodes_disc_BC.B = [nodes_disc_BC.B, nodes_disc.UC{idx_UCs_B(i)}.BL, ...
                       nodes_disc.UC{idx_UCs_B(i)}.B, nodes_disc.UC{idx_UCs_B(i)}.BR];
end
nodes_disc_BC.T = [];
for i = 1:length(idx_UCs_T(:))
    nodes_disc_BC.T = [nodes_disc_BC.T, nodes_disc.UC{idx_UCs_T(i)}.TL, ...
                       nodes_disc.UC{idx_UCs_T(i)}.T, nodes_disc.UC{idx_UCs_T(i)}.TR];
end
nodes_disc_BC.L = [];
for i = 1:length(idx_UCs_L(:))
    nodes_disc_BC.L = [nodes_disc_BC.L, nodes_disc.UC{idx_UCs_L(i)}.BL, ...
                       nodes_disc.UC{idx_UCs_L(i)}.L, nodes_disc.UC{idx_UCs_L(i)}.TL];
end
nodes_disc_BC.R = [];
for i = 1:length(idx_UCs_R(:))
    nodes_disc_BC.R = [nodes_disc_BC.R, nodes_disc.UC{idx_UCs_R(i)}.BR, ...
                       nodes_disc.UC{idx_UCs_R(i)}.R, nodes_disc.UC{idx_UCs_R(i)}.TR];
end
nodes_disc_BC.all = unique([nodes_disc_BC.B, nodes_disc_BC.T, ...
                            nodes_disc_BC.L, nodes_disc_BC.R]);

%% Force response locations
% All disconnected force input location nodes
idx_UCs_F = idx_UCs(forc_loc(1),forc_loc(2));
nodes_disc_FR.F = nodes_disc.UC{idx_UCs_F}.all(forc_loc(3));
% All disconnected response location nodes
idx_UCs_response = idx_UCs(resp_loc(1),resp_loc(2));
nodes_disc_FR.R = nodes_disc.UC{idx_UCs_response}.all(resp_loc(3));


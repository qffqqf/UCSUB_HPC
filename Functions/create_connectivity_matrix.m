function [B_connectivity] = create_connectivity_matrix(dofs_disconnected,alldofs_disconnected,N_UC_x,N_UC_y)

%% Unit cells indices
% Unit cells are ordered first in x-direction, then in y-direction,
% this means that :
% the 1st UC is the utter bottom left UC,
% the N_UC_xth UC is the bottom right UC,
% the N_UC_x*(N_UC_y-1)-1th UC is the top left UC and
% the N_UC_x*N_UC_yth UC is the utter top right UC.
idx_UCs = reshape([1:N_UC_x*N_UC_y],N_UC_x,N_UC_y);

% Indices master UCs for x direction (LR-connectivity): the Left UC
idx_UCs_LR = idx_UCs(1:N_UC_x-1,1:N_UC_y);
% Indices master UCs for y direction (BT-connectivity): the Bottom UC
idx_UCs_BT = idx_UCs(1:N_UC_x,1:N_UC_y-1);


%% Initialise the connected dof vectors
dofs_connected = dofs_disconnected;% Per UC
alldofs_connected = alldofs_disconnected;% For the entire structure


%% Determine connected indices of node and DOF groups
% along x-direction (LR-connectivity)
if ~isempty(idx_UCs_LR)
    for i=1:length(idx_UCs_LR(:))
        idx_UC_master = idx_UCs_LR(i);
        idx_UC_slave =  idx_UCs_LR(i)+1;
        
        % Apply connectivity to dofs
        dofs_connected{idx_UC_slave}.BL = dofs_connected{idx_UC_master}.BR;
        dofs_connected{idx_UC_slave}.L = dofs_connected{idx_UC_master}.R;
        dofs_connected{idx_UC_slave}.TL = dofs_connected{idx_UC_master}.TR;
        % Determine novel, connected dof vector
        alldofs_connected([dofs_disconnected{idx_UC_slave}.BL dofs_disconnected{idx_UC_slave}.L dofs_disconnected{idx_UC_slave}.TL]) = alldofs_connected([dofs_connected{idx_UC_slave}.BL dofs_connected{idx_UC_slave}.L dofs_connected{idx_UC_slave}.TL]);
    end
end

% along y-direction (BT-connectivity)
if ~isempty(idx_UCs_BT)
    for i=1:length(idx_UCs_BT(:))
        idx_UC_master = idx_UCs_BT(i);
        idx_UC_slave =  idx_UCs_BT(i)+N_UC_x;
        
        % Apply connectivity to dofs
        dofs_connected{idx_UC_slave}.BL = dofs_connected{idx_UC_master}.TL;
        dofs_connected{idx_UC_slave}.B = dofs_connected{idx_UC_master}.T;
        dofs_connected{idx_UC_slave}.BR = dofs_connected{idx_UC_master}.TR;
        % Determine novel, connected dof vector
        alldofs_connected([dofs_disconnected{idx_UC_slave}.BL dofs_disconnected{idx_UC_slave}.B dofs_disconnected{idx_UC_slave}.BR]) = alldofs_connected([dofs_connected{idx_UC_slave}.BL dofs_connected{idx_UC_slave}.B dofs_connected{idx_UC_slave}.BR]);
    end
end

% Determine novel dof numbers for the alldofs_connected vector
% (renumbered to correspond to the amount of dofs in the final assembly)
[~,~,alldofs_connected_renumbered] = unique(alldofs_connected);

%% Set up connectivity matrix
B_connectivity = sparse(alldofs_disconnected,alldofs_connected_renumbered,1);
end
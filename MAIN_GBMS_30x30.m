clear all
% close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SETUP  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath([pwd,'/Functions'])
addpath([pwd,'/UCs/bareplates'])
addpath([pwd,'/UCs/vbplates/tet'])
addpath([cd,'/UCs/vbplates/hex'])
% tStart = cputime;
tStart = tic;
%% UC model setup
% Assuming renumbered nodes [1:max_number_nodes], corresponding to
% the matrices!
filename = 'uclvb_spcoarse';
import_option = 'mat'; %'mat' or 'op4' (note that op4 requires binary reader)
n_dofnodes = 3; % Amount of DOFs per node: shells = 6; solids = 3; solids & shells: = 6
z_id = 3;
savename = 'vbplate_30x30_gbms';

%% Finite structure setup
% Composition amount of unit cells
N_UC_x = 30;% Amount of UCs in X-direction
N_UC_y = 30;% Amount of UCs in Y-direction

%% Load single UC model (identical UCs)
disp('Loading UC model ...')
load([filename,'.mat']);
% F_UC_FOM = comsol_matrices.forc;
[K_UC_FOM, M_UC_FOM, C_UC_FOM, UC_nodes, UC_coordinates, L_UC_x, L_UC_y] = import_FE_UC3(filename);

%% UC model order reduction setup
reduction_I_tf = 1; % 0 or 1; Interior modal reduction [BMS]
reduction_A_tf = 1; % 0 or 1; Boundary modal reduction [GBMS] (only 1 if reduction_I_tf = 1)
method_ort = 'SVD'; % 'SVD' or 'QR'; orthogonalization method for boundary modes in the [GBMS] method, SVD default
tol_SVD = 1e-15; % Tolerance of SVD to discard superfluous modes from orthogonalized boundary mode sets (cf https://doi.org/10.1016/j.ymssp.2018.05.031 p584)
n_modes_I = 30; % Number of interior fixed-interface mode shapes [BMS] (choose (e.g. 2x max freq of interest), but always inspect results): off-line pretest available
n_modes_A = 50; % Number of boundary mode shapes/boundary [GBMS] (choose (e.g. 3x interior modes), but always inspect results): off-line pretest available
damping = 1e-2;

%% Final assembly MOR
reduction_T_tf = 1; % Final global modal reduction of the assembled system
n_modes_T = 11000; % Number of global modes
static_enrich_tf = 1; % Use static enrichment for input force yes/no

%% Analysis setup
% Eigenmodes analysis
n_modes = 80; % Number of eigenmodes to be calculated for the full assembeld system
solve_eig_tf = 0; % Boolean to switch eigenmode analysis on/off

% Boundary conditions for the full assembled system
BC = 1; % fixed boundary dofs = 1, free boundary dofs = 0

% Forced response analysis frequencies
freq = 0:2:1000; % [Hz]
omega = 2*pi*freq; % [rad/s]
plot_nfreq = 1000; % Which frequency to plot full field response of?
solve_frf_tf = 1; % Boolean to switch FRF analysis on/off

% Single input point force (uz)
force_UC_x = 12; % Which UC in the x-direction 12
force_UC_y = 12; % Which UC in the y-direction 9
picker_F = @(r) norm(r-[L_UC_x,L_UC_y,0])<1e-10;
force_UC_node = pick_nodes(mesh_data, picker_F);
force_UC_node = force_UC_node(1);
% force_UC_node = 2083; % Which node in the specific UC is excited 2083
% force_UC_node = "all";

% % Single response node (uz)
response_UC_x = 18; % Which UC in the x-direction
response_UC_y = 18; % Which UC in the y-direction
picker_R = @(r) norm(r-[L_UC_x,L_UC_y,0])<1e-10;
response_UC_node = pick_nodes(mesh_data, picker_R);
% response_UC_node = 34; % Which node in the specific UC is the uz response assessed of

% Post-processing
checkplots = 0; % 0 or 1, to verify matrix structure

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SINGLE UC MODEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% !!!!! Hard-coded: assign structural damping !!!!!
K_UC_FOM = K_UC_FOM*(1+damping*1i);

% Inspect original UC matrix structure
if checkplots
    figure
    subplot(121)
    spy(K_UC_FOM)
    title('K UC FOM')
    subplot(122)
    spy(M_UC_FOM)
    title('M UC FOM')
end

% Calculate dofs from nodes
if n_dofnodes == 6
    n_nodes = length(K_UC_FOM)/6;
    delete_rows_UC_FOM = zeros(1,3*n_nodes);% Hard delete
    for i=1:n_nodes
        delete_rows_UC_FOM((i-1)*3+1:i*3)=((i*6-2):1:(i*6));
    end
    K_UC_FOM(delete_rows_UC_FOM,:) = [];
    K_UC_FOM(:,delete_rows_UC_FOM) = [];
    M_UC_FOM(delete_rows_UC_FOM,:) = [];
    M_UC_FOM(:,delete_rows_UC_FOM) = [];
    C_UC_FOM(delete_rows_UC_FOM,:) = [];
    C_UC_FOM(:,delete_rows_UC_FOM) = [];
end

UC_dofs = calculate_dof_indices(UC_nodes,mesh_data);
n_nodes_UC = UC_nodes.nNodes;
n_dofs_UC_FOM = length(K_UC_FOM);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SINGLE UC MOR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Apply BMS
if reduction_I_tf
    disp('Interior modal reduction UC ...')
    [M_UC_redI, K_UC_redI, C_UC_redI, UC_dofs_redI, Phi_I, Psi_IA, V_UC_I] = interior_modal_reduction(M_UC_FOM,K_UC_FOM,C_UC_FOM,UC_dofs,n_modes_I);
end

%% Apply GBMS (only applicable if BMS is applied first!)
if reduction_I_tf&&reduction_A_tf % if BMS has already been performed
    % Note here: additional cleaning of the eigenvectors possible by
    % removing zero rows and re-inserting them
    disp('Boundary modal reduction UC ...')
    [M_UC_red,K_UC_red,C_UC_red,UC_dofs_red,V_UC,L_A] = boundary_modal_reduction(M_UC_FOM,K_UC_FOM,C_UC_FOM,UC_dofs,M_UC_redI,K_UC_redI,UC_dofs_redI,n_modes_A,Phi_I,Psi_IA,method_ort,tol_SVD);
elseif reduction_I_tf&&~reduction_A_tf% If only BMS is to be performed
    M_UC_red = M_UC_redI;
    K_UC_red = K_UC_redI;
    C_UC_red = C_UC_redI;
    UC_dofs_red = UC_dofs_redI;
    V_UC = V_UC_I;
elseif ~reduction_I_tf&&~reduction_A_tf
    V_UC = speye(n_dofs_UC_FOM);
elseif ~reduction_I_tf&&reduction_A_tf
    error('Enable interior modal reduction to apply boundary modal reduction!')
end

%% Verify reduced model eigenvalues
if checkplots && reduction_I_tf
    % Check and compare eigensolutions full and reduced
    [f_n_UC_full, f_n_UC_red] = verify_UC_ROMs(K_UC_FOM,M_UC_FOM,K_UC_red,M_UC_red,20);
    rel_error_f_n = abs((f_n_UC_full-f_n_UC_red)./f_n_UC_full);
    if checkplots
        figure
        spy(V_UC)
        title('V UC ROM')
        figure
        subplot(121)
        spy(K_UC_red)
        title('K UC ROM')
        subplot(122)
        spy(M_UC_red)
        title('M UC ROM')
    end
end

%% Continue with FOM or ROM matrices
if reduction_I_tf||reduction_A_tf % UC ROM
    K_UC = K_UC_red;
    M_UC = M_UC_red;
    C_UC = C_UC_red;
    UC_dofs = UC_dofs_red;
    n_dofs_UC = length(K_UC);
else % UC FOM
    K_UC = K_UC_FOM;
    M_UC = M_UC_FOM;
    C_UC = C_UC_FOM;
    n_dofs_UC = n_dofs_UC_FOM;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DISCONNECTED UCs: PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Determine real coordinates of the N_UC_x*N_UC_y disconnected UCs
% Coordinates per disconnected UC
[yy,xx] = meshgrid(1:N_UC_y,1:N_UC_x);
coords_disconnected = cell(1,N_UC_x*N_UC_y);
for i = 1:N_UC_x*N_UC_y
    coords_disconnected{i}(:,1) = UC_coordinates(:,1) + (i-1)*n_nodes_UC;
    coords_disconnected{i}(:,2) = UC_coordinates(:,2) + (xx(i)-1)*L_UC_x;
    coords_disconnected{i}(:,3) = UC_coordinates(:,3) + (yy(i)-1)*L_UC_y;
    coords_disconnected{i}(:,4) = UC_coordinates(:,4);
end
% Coordinates of all disconnected UCs
allcoords_disconnected = cat(1,coords_disconnected{:});

%% Determine all original nodes of the N_UC_x*N_UC_y disconnected UCs
% Disconnected nodes per UC
nodes_disconnected = cell(1,N_UC_x*N_UC_y);
for i = 1:N_UC_x*N_UC_y
    nodes_disconnected{i}.I = UC_nodes.I + (i-1)*n_nodes_UC;
    nodes_disconnected{i}.L = UC_nodes.L + (i-1)*n_nodes_UC;
    nodes_disconnected{i}.R = UC_nodes.R + (i-1)*n_nodes_UC;
    nodes_disconnected{i}.B = UC_nodes.B + (i-1)*n_nodes_UC;
    nodes_disconnected{i}.T = UC_nodes.T + (i-1)*n_nodes_UC;
    nodes_disconnected{i}.BL = UC_nodes.BL + (i-1)*n_nodes_UC;
    nodes_disconnected{i}.BR = UC_nodes.BR + (i-1)*n_nodes_UC;
    nodes_disconnected{i}.TR = UC_nodes.TR + (i-1)*n_nodes_UC;
    nodes_disconnected{i}.TL = UC_nodes.TL + (i-1)*n_nodes_UC;
    nodes_disconnected{i}.all = sort([nodes_disconnected{i}.I nodes_disconnected{i}.L nodes_disconnected{i}.R nodes_disconnected{i}.B nodes_disconnected{i}.T nodes_disconnected{i}.BL nodes_disconnected{i}.BR nodes_disconnected{i}.TR nodes_disconnected{i}.TL]);
end

% All disconnected nodes
allnodes_disconnected = [1:N_UC_x*N_UC_y*n_nodes_UC];% For the entire structure

% All disconnected interior nodes (for visualization)
allnodes_disconnected_interior = [];
for i = 1:N_UC_x*N_UC_y
    allnodes_disconnected_interior = [allnodes_disconnected_interior nodes_disconnected{i}.I];
end

% All disconnected interface nodes (for visualization)
allnodes_disconnected_interface = setdiff(allnodes_disconnected,allnodes_disconnected_interior);

% All disconnected outer boundary nodes (for boundary condition application)
% Corresponding to connectivity rule:
% Indices of all UCs sorted (from low to high first along x-direction, next
% along y-direction)
idx_UCs = reshape([1:N_UC_x*N_UC_y],N_UC_x,N_UC_y);
% Indices of bottom row of UCs
idx_UCs_B = idx_UCs(:,1);
% Indices of top row of UCs
idx_UCs_T = idx_UCs(:,end);
% Indices of left row of UCs
idx_UCs_L = idx_UCs(1,:);
% Indices of right row of UCs
idx_UCs_R = idx_UCs(end,:);

allnodes_disconnected_BC_B = [];
for i = 1:length(idx_UCs_B(:))
    allnodes_disconnected_BC_B = [allnodes_disconnected_BC_B nodes_disconnected{idx_UCs_B(i)}.BL nodes_disconnected{idx_UCs_B(i)}.B nodes_disconnected{idx_UCs_B(i)}.BR];
end
allnodes_disconnected_BC_T = [];
for i = 1:length(idx_UCs_T(:))
    allnodes_disconnected_BC_T = [allnodes_disconnected_BC_T nodes_disconnected{idx_UCs_T(i)}.TL nodes_disconnected{idx_UCs_T(i)}.T nodes_disconnected{idx_UCs_T(i)}.TR];
end
allnodes_disconnected_BC_L = [];
for i = 1:length(idx_UCs_L(:))
    allnodes_disconnected_BC_L = [allnodes_disconnected_BC_L nodes_disconnected{idx_UCs_L(i)}.BL nodes_disconnected{idx_UCs_L(i)}.L nodes_disconnected{idx_UCs_L(i)}.TL];
end
allnodes_disconnected_BC_R = [];
for i = 1:length(idx_UCs_R(:))
    allnodes_disconnected_BC_R = [allnodes_disconnected_BC_R nodes_disconnected{idx_UCs_R(i)}.BR nodes_disconnected{idx_UCs_R(i)}.R nodes_disconnected{idx_UCs_R(i)}.TR];
end
allnodes_disconnected_BC = unique([allnodes_disconnected_BC_B allnodes_disconnected_BC_T allnodes_disconnected_BC_L allnodes_disconnected_BC_R]);

% All disconnected force input location nodes
idx_UCs_F = idx_UCs(force_UC_x,force_UC_y);
if force_UC_node == 'all'
    allnodes_disconnected_F = nodes_disconnected{idx_UCs_F}.all(:);
else
    allnodes_disconnected_F = nodes_disconnected{idx_UCs_F}.all(force_UC_node);
end

% All disconnected response location nodes
idx_UCs_response = idx_UCs(response_UC_x(:),response_UC_y(:));
for i = 1:length(idx_UCs_response(:))
    allnodes_disconnected_response(i) = nodes_disconnected{idx_UCs_response(i)}.all(response_UC_node);
end

if checkplots
    % Visualize the nodes
    figure
    scatter3(allcoords_disconnected(allnodes_disconnected_interior,2),allcoords_disconnected(allnodes_disconnected_interior,3),allcoords_disconnected(allnodes_disconnected_interior,4),2,'linewidth',0.5)
    hold on
    scatter3(allcoords_disconnected(allnodes_disconnected_interface,2),allcoords_disconnected(allnodes_disconnected_interface,3),allcoords_disconnected(allnodes_disconnected_interface,4),9,'+','linewidth',0.5)
    scatter3(allcoords_disconnected(allnodes_disconnected_BC,2),allcoords_disconnected(allnodes_disconnected_BC,3),allcoords_disconnected(allnodes_disconnected_BC,4),9,'s','linewidth',0.5)
    scatter3(allcoords_disconnected(allnodes_disconnected_F,2),allcoords_disconnected(allnodes_disconnected_F,3),allcoords_disconnected(allnodes_disconnected_F,4),20,'v','linewidth',2)
    scatter3(allcoords_disconnected(allnodes_disconnected_response,2),allcoords_disconnected(allnodes_disconnected_response,3),allcoords_disconnected(allnodes_disconnected_response,4),20,'^','linewidth',2)
    hold off
    legend('Interior nodes','Interface nodes','Outer boundary','Input node','Response node')
    xlabel('X [m]');ylabel('Y [m]');zlabel('Z [m]')
    axis equal
    view([0 0 1])
end

%% Determine all dofs of the disconnected UC models (full or reduced)
% Dofs per UC
dofs_disconnected = cell(1,N_UC_x*N_UC_y);
for i = 1:N_UC_x*N_UC_y
    dofs_disconnected{i}.I = UC_dofs.I + (i-1)*n_dofs_UC;
    dofs_disconnected{i}.L = UC_dofs.L + (i-1)*n_dofs_UC;
    dofs_disconnected{i}.R = UC_dofs.R + (i-1)*n_dofs_UC;
    dofs_disconnected{i}.B = UC_dofs.B + (i-1)*n_dofs_UC;
    dofs_disconnected{i}.T = UC_dofs.T + (i-1)*n_dofs_UC;
    dofs_disconnected{i}.BL = UC_dofs.BL + (i-1)*n_dofs_UC;
    dofs_disconnected{i}.BR = UC_dofs.BR + (i-1)*n_dofs_UC;
    dofs_disconnected{i}.TR = UC_dofs.TR + (i-1)*n_dofs_UC;
    dofs_disconnected{i}.TL = UC_dofs.TL + (i-1)*n_dofs_UC;
end

% All disconnected dofs
alldofs_disconnected = [1:N_UC_x*N_UC_y*n_dofs_UC];% For the entire structure
alldofs_disconnected_FOM = [1:N_UC_x*N_UC_y*n_dofs_UC_FOM];% For the entire structure

% All disconnected outer boundary dofs (for boundary condition application)
alldofs_disconnected_BC_B = [];
for i = 1:length(idx_UCs_B(:))
    alldofs_disconnected_BC_B = [alldofs_disconnected_BC_B dofs_disconnected{idx_UCs_B(i)}.BL dofs_disconnected{idx_UCs_B(i)}.B dofs_disconnected{idx_UCs_B(i)}.BR];
end
alldofs_disconnected_BC_T = [];
for i = 1:length(idx_UCs_T(:))
    alldofs_disconnected_BC_T = [alldofs_disconnected_BC_T dofs_disconnected{idx_UCs_T(i)}.TL dofs_disconnected{idx_UCs_T(i)}.T dofs_disconnected{idx_UCs_T(i)}.TR];
end
alldofs_disconnected_BC_L = [];
for i = 1:length(idx_UCs_L(:))
    alldofs_disconnected_BC_L = [alldofs_disconnected_BC_L dofs_disconnected{idx_UCs_L(i)}.BL dofs_disconnected{idx_UCs_L(i)}.L dofs_disconnected{idx_UCs_L(i)}.TL];
end
alldofs_disconnected_BC_R = [];
for i = 1:length(idx_UCs_R(:))
    alldofs_disconnected_BC_R = [alldofs_disconnected_BC_R dofs_disconnected{idx_UCs_R(i)}.BR dofs_disconnected{idx_UCs_R(i)}.R dofs_disconnected{idx_UCs_R(i)}.TR];
end
alldofs_disconnected_BC = unique([alldofs_disconnected_BC_B alldofs_disconnected_BC_T alldofs_disconnected_BC_L alldofs_disconnected_BC_R]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FINITE SYSTEM ASSEMBLY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Finite structure matrix - disconnected UCs
% (Sparse) block diagonal matrix collecting the subsystem matrices on the
% diagonal
K_disconnected = kron(speye(N_UC_x*N_UC_y),K_UC);
M_disconnected = kron(speye(N_UC_x*N_UC_y),M_UC);
C_disconnected = kron(speye(N_UC_x*N_UC_y),C_UC);

% Column vector collecting the input force contribution: project with
% reduction matrices
F_disconnected_FOM = zeros(length(alldofs_disconnected_FOM),1);
if force_UC_node == 'all'
    alledof_uc = (allnodes_disconnected_F-1)*n_dofnodes;
    alldofs_disconnected_F = sort([alledof_uc+1, alledof_uc+2, alledof_uc+3]);
    F_disconnected_FOM(alldofs_disconnected_F) = F_UC_FOM;
else
    alldofs_disconnected_F = (allnodes_disconnected_F-1)*n_dofnodes+z_id;
    F_disconnected_FOM(alldofs_disconnected_F) = 1;
end
if reduction_I_tf
    disp('Force projection UCs ...')
    % input_selector = speye(N_UC_x*N_UC_y);
    input_UC_selector = spalloc(N_UC_x*N_UC_y,N_UC_x*N_UC_y,length(idx_UCs_F(:)));
    for i=1:length(idx_UCs_F(:))
        input_UC_selector(idx_UCs_F(i),idx_UCs_F(i))=1;
    end
    V = kron(input_UC_selector,V_UC);
    F_disconnected = V.'*F_disconnected_FOM;
else
    F_disconnected = F_disconnected_FOM;
end
if checkplots
    figure
    subplot(121)
    spy(K_disconnected)
    title('K Disconnected system')
    subplot(122)
    spy(M_disconnected)
    title('M Disconnected system')
end

%% Create connectivity matrix
% Non-square boolean matrix which connects the DOFs between unassembled and
% assembled system
disp('Creating connectivity matrix ...')
[B_connectivity] = create_connectivity_matrix(dofs_disconnected,alldofs_disconnected,N_UC_x,N_UC_y);
if checkplots
    figure
    spy(B_connectivity)
end

%% Finite structure matrix - assembly
disp('Assembling system matrices ...')
K_connected = B_connectivity.'*K_disconnected*B_connectivity;
M_connected = B_connectivity.'*M_disconnected*B_connectivity;
C_connected = B_connectivity.'*C_disconnected*B_connectivity;
F_connected = B_connectivity.'*F_disconnected;

if checkplots
    figure
    subplot(121)
    spy(K_connected)
    title('K Assembled system')
    subplot(122)
    spy(M_connected)
    title('M Assembled system')
end

%% Apply boundary conditions
disp('Applying BCs ...')
% Direct eliminiation of boundary DOFs of the assembled system
BC_disconnected = ones(length(alldofs_disconnected),1);
BC_disconnected(alldofs_disconnected_BC) = ~BC;
BC_disconnected = sparse(1:length(alldofs_disconnected),1:length(alldofs_disconnected),BC_disconnected);

BC_connected = B_connectivity.'*BC_disconnected;
dofs_delete_BC = find(~any(BC_connected,2)).';

% Delete zero rows/columns to avoid singularities
rows_delete = [find(~any(K_connected,2)).', dofs_delete_BC];
rows_keep = setdiff(1:length(B_connectivity(1,:)),rows_delete);
K_connected(rows_delete,:) = [];
K_connected(:,rows_delete) = [];
M_connected(rows_delete,:) = [];
M_connected(:,rows_delete) = [];
C_connected(rows_delete,:) = [];
C_connected(:,rows_delete) = [];
F_connected(rows_delete) = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SOLVE FINITE PLATE EIGENMODES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if solve_eig_tf
    Eigenmodes = zeros(length(B_connectivity(1,:)),n_modes);
    disp('Solving eigenmodes ...')
    [Phi,lambda] = eigs(K_connected,M_connected,n_modes,0.0001);%0.0001
    
    Eigenmodes(rows_keep,:) = Phi;
    omega_eig = sqrt(diag(lambda));
    freq_eig = omega_eig/(2*pi);
    
    %
    %     %% Post process eigenmodes
    %     % z-displacement
    %     plot_nmode = 59;
    %     scaling = 1e-4;%1
    %
    %     % Convert modal results to physical results
    %     if reduction_I_tf||reduction_A_tf
    %         Eigenmodes_physical = kron(speye(N_UC_x*N_UC_y),V_UC)*B_connectivity*Eigenmodes;
    %     else
    %         Eigenmodes_physical = B_connectivity*Eigenmodes;
    %     end
    %
    %     % Obtain ux, uy and uz components
    %     ux = Eigenmodes_physical(1:n_dofnodes:end,plot_nmode);
    %     uy = Eigenmodes_physical(2:n_dofnodes:end,plot_nmode);
    %     uz = Eigenmodes_physical(3:n_dofnodes:end,plot_nmode);
    %
    %     figure
    %     scatter3(allcoords_disconnected(:,2)+scaling*ux(:),allcoords_disconnected(:,3)+scaling*uy(:),allcoords_disconnected(:,4)+scaling*uz(:),5,scaling*abs((ux(:).^2+uy(:).^2+uz(:).^2).^0.5))
    %     xlabel('X [m]');ylabel('Y [m]');zlabel('Z [m]')
    %     % axis equal
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SOLVE FINITE PLATE FORCED RESPONSE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if solve_frf_tf
    if reduction_T_tf
        %% USING GLOBAL MODAL REDUCTION
        % Calculate real eigenmodes of assembly
        disp('Calculating assembly eigenmodes ...')
        [Phi_T,lambda_T] = eigs(real(K_connected),real(M_connected),n_modes_T,0.0001);%0.0001
        eigfreqs = (diag(lambda_T).^0.5)/(2*pi);
        f_max = eigfreqs(end);
        
        if static_enrich_tf
            % Calculate static enrichment vector
            disp('Calculating static enrichment vector ...')
            Psi_T = K_connected\F_connected;
            % Mass normalization
            Psi_T = Psi_T/(Psi_T.'*M_connected*Psi_T)^0.5;
            % Compose basis
            Phi_T = [Phi_T,Psi_T];
            % % SVD (Problematic for large basis vectors + does not improve)
            % [Psi_T_U,Phi_T_S,~] = svd(Phi_T);
            % Psi_T = Psi_T_U(:,(diag(Phi_T_S)/Phi_T_S(1,1))>tol_SVD);
        end
        
        % Project assembly system matrices and loads
        disp('Projecting assembly system matrices ...')
        K_connected_red = Phi_T.'*K_connected*Phi_T;
        M_connected_red = Phi_T.'*M_connected*Phi_T;
        C_connected_red = Phi_T.'*C_connected*Phi_T;
        F_connected_red = Phi_T.'*F_connected;
        
        % Direct FRF calculation
        Response_connected_red = zeros(length(K_connected_red(:,1)),length(freq));
        disp('Solving FRF ...')
        for i=1:length(freq)
            i
            % Dynamic stiffness matrix
            Z_connected_red = K_connected_red+1i*omega(i)*C_connected_red-omega(i)^2*M_connected_red;
            % Solve
            Response_connected_red(:,i) =Z_connected_red\F_connected_red;
        end
        % Back-projection
        Response_connected = zeros(length(B_connectivity(1,:)),length(freq));
        disp('Back-projection modal assembly results ...')
        Response_connected(rows_keep,:) = Phi_T*Response_connected_red;
    else
        %% USING ASSEMBLY WITHOUT ADDITIONAL MODAL REDUCTION
        Response_connected = zeros(length(B_connectivity(1,:)),length(freq));
        disp('Solving FRF ...')
        for i=1:length(freq)
	    i
            % Dynamic stiffness matrix
            Z_connected = K_connected+1i*omega(i)*C_connected-omega(i)^2*M_connected;
            % Solve
            Response_connected(rows_keep,i)= Z_connected\F_connected;
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% POST PROCESSING FINITE PLATE FORCED RESPONSE %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Post process forced response field
% % z-displacement
% scaling =1e2;
%
% % Convert modal results to physical results
% if reduction_I_tf||reduction_A_tf
%     Response_physical = kron(speye(N_UC_x*N_UC_y),V_UC)*B_connectivity*Response_connected(:,plot_nfreq);
% else
%     Response_physical = B_connectivity*Response_connected(:,plot_nfreq);
% end
%
% % Obtain ux, uy and uz components
% ux = Response_physical(1:n_dofnodes:end);
% uy = Response_physical(2:n_dofnodes:end);
% uz = Response_physical(3:n_dofnodes:end);
%
% figure
% scatter3(allcoords_disconnected(:,2)+scaling*ux(:),allcoords_disconnected(:,3)+scaling*uy(:),allcoords_disconnected(:,4)+scaling*uz(:),5,scaling*abs((ux(:).^2+uy(:).^2+uz(:).^2).^0.5))
% xlabel('X [mm]');ylabel('Y [mm]');zlabel('Z [mm]')
%
% figure
% scatter3(allcoords_disconnected(:,2)+scaling*ux(:),allcoords_disconnected(:,3)+scaling*uy(:),allcoords_disconnected(:,4)+scaling*uz(:),5,log10(scaling*abs((ux(:).^2+uy(:).^2+uz(:).^2).^0.5)))
% xlabel('X [mm]');ylabel('Y [mm]');zlabel('Z [mm]')
% axis equal

%% Post process FRF plot uz component
alldofs_disconnected_response = (allnodes_disconnected_response-1)*n_dofnodes+z_id; % uz_response
output_selection = spalloc(length(alldofs_disconnected_response),length(alldofs_disconnected_FOM),length(alldofs_disconnected_response));
for i = 1:length(alldofs_disconnected_response)
    output_selection(i,alldofs_disconnected_response(i)) = 1;
end
if reduction_I_tf||reduction_A_tf
    disp('Back-projection UC results ...')
    output_UC_selector = spalloc(N_UC_x*N_UC_y,N_UC_x*N_UC_y,length(idx_UCs_response(:)));
    for i = 1:length(idx_UCs_response(:))
        output_UC_selector(idx_UCs_response(i),idx_UCs_response(i))=1;
    end
    V = kron(output_UC_selector,V_UC);
    uz_gbms = output_selection*V*B_connectivity*Response_connected;
else
    uz_gbms = output_selection*B_connectivity*Response_connected;
end

figure
semilogy(freq,abs(uz_gbms), 'Linewidth', 2)
title('Point to point FRF')
xlabel('Frequency [Hz]');ylabel('|u_z| [m]')
set(gca, 'FontSize', 20)

% timing = cputime - tStart
timing = toc(tStart);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SAVING RESULTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save([savename,'.mat'],'freq','uz_gbms','timing')

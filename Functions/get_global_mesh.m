function mesh_data_global = get_global_mesh(N_UC_x, N_UC_y, mesh_data, coords_disc)

%% nodal coords
mesh_data_global.nd_data = coords_disc.all(:,2:4);
mesh_data_global.element_type = mesh_data.element_type;
nNode = size(mesh_data_global.nd_data,1);
mesh_data_global.nv_data = [];
mesh_data_global.nDOF = 3*nNode;

%% elements
ele_disc = cell(1,N_UC_x*N_UC_y);
n_nd_UC = size(mesh_data.nd_data,1);
n_val_UC = sum(mesh_data.nv_data>0,'all');
existing_id = find(mesh_data.nv_data>0);
mesh_data_global.me_data = containers.Map;
mesh_data_global.me_data("Plate") = [];
mesh_data_global.me_data("Air") = [];

for i = 1:N_UC_x*N_UC_y
    if i<N_UC_x*N_UC_y/2+1
        ele_disc{i} = mesh_data.en_data + (i-1)*n_nd_UC;
        mesh_data_global.me_data("Plate") = [mesh_data_global.me_data("Plate");...
                                         mesh_data.me_data("Plate") + (i-1)*n_nd_UC];
        mesh_data_global.me_data("Air") = [mesh_data_global.me_data("Air");...
                                       mesh_data.me_data("Air") + (i-1)*n_nd_UC];
    end
    nv_data = mesh_data.nv_data;
    nv_data(existing_id) = nv_data(existing_id) + (i-1)*n_val_UC;
    mesh_data_global.nv_data = [mesh_data_global.nv_data; nv_data];
    
end
mesh_data_global.en_data = cat(1, ele_disc{:});

%% surf elements
% ele_disc = cell(1,N_UC_x*N_UC_y);
% for i = 1:N_UC_x*N_UC_y
%     ele_disc{i} = mesh_data.en_data_surf + (i-1)*n_nd_UC;
% end
% mesh_data_global.en_data_surf = cat(1, ele_disc{:});

%% set type
% mesh_data_global.element_type = mesh_data.element_type;
clc
clear all;
addpath("C:\modeling\scom_code\fem\solid_tet");
addpath("C:\modeling\scom_code\fem\solid_tri");
addpath("C:\modeling\scom_code\fem\comsol");
addpath("C:\modeling\experiments\UC_GBMS_2D_Substructuring_FiniteStructureResponse\UCs\vbplates\hex");

%% read plate mesh
element_type = "C3D27";
file_name = "vbplate_coarse_minus";
output = "vbplate_coarse_minus";
mat_type = ["PMMA", "Aluminum"];
surf_set = [];
mesh_data = read_mesh(file_name, element_type, mat_type, surf_set);
alpha = 1;
display_surf_mesh(mesh_data, "all", 0, 10, alpha)
display_mesh(mesh_data, "all", 0, 10, alpha)

%% read matrices
omega = 1e3;
comsol_matrices.stfK = read_matrices(file_name, 0,0,0);
Kc = read_matrices(file_name, omega/2/pi,0,0);
comsol_matrices.mass = (comsol_matrices.stfK - Kc)/omega^2;

%% set element type
nGauss = 3;
interp_data = generate_tet10_interp(nGauss);
interp_data_surf = generate_tri6_interp(nGauss);

%% Set material type
mat_data = containers.Map;
Aluminum_param.density = 2.7e3;
Aluminum_param.youngs_modulus = 7e10;
Aluminum_param.poissons_ratio = 0.3;
mat_data("Aluminum") = Aluminum_param;

PMMA_param.density = 1.188e3;
PMMA_param.youngs_modulus = 4.85e9;
PMMA_param.poissons_ratio = 0.31;
mat_data("PMMA") = PMMA_param;

save(strcat(output, ".mat"), 'mesh_data', 'interp_data', 'interp_data_surf', 'mat_data', 'comsol_matrices');
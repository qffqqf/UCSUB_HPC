function [Dofs] = calculate_dof_indices(Nodes,mesh_data)
% Collect interior dofs in I
Dofs.I=nonzeros(mesh_data.nv_data(Nodes.I,:))';
% Separate boundary dofs
Dofs.L=nonzeros(mesh_data.nv_data(Nodes.L,:))';
Dofs.R=nonzeros(mesh_data.nv_data(Nodes.R,:))';
Dofs.B=nonzeros(mesh_data.nv_data(Nodes.B,:))';
Dofs.T=nonzeros(mesh_data.nv_data(Nodes.T,:))';
Dofs.BL=nonzeros(mesh_data.nv_data(Nodes.BL,:))';
Dofs.BR=nonzeros(mesh_data.nv_data(Nodes.BR,:))';
Dofs.TL=nonzeros(mesh_data.nv_data(Nodes.TL,:))';
Dofs.TR=nonzeros(mesh_data.nv_data(Nodes.TR,:))';
% Collect all boundary dofs in A
Dofs.A = [Dofs.L,Dofs.R,Dofs.B,Dofs.T,Dofs.BL,Dofs.BR,Dofs.TR,Dofs.TL];
Dofs.nDOF = numel([Dofs.I, Dofs.A]);
end
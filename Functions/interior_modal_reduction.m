function [M_redI,K_redI,C_redI,Dofs_redI,Phi_I,Psi_IA,V] = interior_modal_reduction(M,K,C,Dofs,n_modes_I)
%% Verify if n_modes_I is not larger than Dofs.I
if n_modes_I>length(Dofs.I)
    error('n_modes_I is larger than amount of interior dofs!')
end

%% Determine matrix partitions
K_II = real(K(Dofs.I,Dofs.I));
M_II = real(M(Dofs.I,Dofs.I));

% Sorting now [I L R B T BL BR TR TL]
K_IA = real(K(Dofs.I,Dofs.A));
% M_IA = M(Dofs.I,Dofs.A);

%% Calculate interior fixed-interface mode shaps to describe interior DOFs
% Delete potential zero rows/columns (to avoid singularities)
delete_rows_I = ~any(K_II,2);
K_II(delete_rows_I,:) = [];
K_II(:,delete_rows_I) = [];
M_II(delete_rows_I,:) = [];
M_II(:,delete_rows_I) = [];
keep_rows_I = setdiff(1:length(Dofs.I),find(delete_rows_I));

% Calculate interior normal modes
Phi_I = zeros(length(Dofs.I),n_modes_I); %initialize
% Non-symmetric EVP
[Phi_I(keep_rows_I,:),eigval] = eigs(K_II,M_II,n_modes_I,0);
% Mass normalization
% Phi_I(keep_rows_I,:) =  Phi_I(keep_rows_I,:)/(Phi_I(keep_rows_I,:).'*M_II*Phi_I(keep_rows_I,:))^0.5;
% Symmetric EVP
% K_II_tildedot = K_II*M_II^(-0.5)
% K_II_tilde = M_II^(-0.5)*K_II*M_II^(-0.5);
% [Phi_I(keep_rows_I,:),~] = eigs(K_II_tilde,n_modes_I,0);
% Phi(keep_rows_I,:)= M_II^(-0.5)*Phi_I(keep_rows_I,:);%


%% Calculate static constraint modes to describe boundary DOFs
% Delete corresponding zero rows/columns (to avoid singularities)
K_IA(delete_rows_I,:) = [];
% delete_columns_A = ~any(K_IA,1); % does not change/help...
% keep_columns_A = setdiff(1:length(Dofs.A),find(delete_columns_A));
% K_IA(:,delete_columns_A) = [];

% Calculate static constraint modes
Psi_IA = zeros(length(Dofs.I),length(Dofs.A));
Psi_IA(keep_rows_I,:) = -K_II\K_IA;%
% Psi_IA(keep_rows_I,keep_columns_A) = -K_II\K_IA;% does not
% change/help...

%% Reorder Dofs indices according to new Dofs.I (n_modes_I) and Dofs.A
Dofs_redI.I = [1:n_modes_I];
Dofs_redI.L = [n_modes_I+1:n_modes_I+length(Dofs.L)];
Dofs_redI.R = [Dofs_redI.L(end)+1:Dofs_redI.L(end)+length(Dofs.R)];
Dofs_redI.B = [Dofs_redI.R(end)+1:Dofs_redI.R(end)+length(Dofs.B)];
Dofs_redI.T = [Dofs_redI.B(end)+1:Dofs_redI.B(end)+length(Dofs.T)];
Dofs_redI.BL = [Dofs_redI.T(end)+1:Dofs_redI.T(end)+length(Dofs.BL)];
Dofs_redI.BR = [Dofs_redI.BL(end)+1:Dofs_redI.BL(end)+length(Dofs.BR)];
Dofs_redI.TR = [Dofs_redI.BR(end)+1:Dofs_redI.BR(end)+length(Dofs.TR)];
Dofs_redI.TL = [Dofs_redI.TR(end)+1:Dofs_redI.TR(end)+length(Dofs.TL)];
Dofs_redI.A = [Dofs_redI.L,Dofs_redI.R,Dofs_redI.B,Dofs_redI.T,Dofs_redI.BL,Dofs_redI.BR,Dofs_redI.TR,Dofs_redI.TL];

%% Determine modal transformation matrix
% Initialize
V = zeros(length(Dofs.I)+length(Dofs.A),n_modes_I+length(Dofs.A));

% Assign to corresponding DOFs
V(Dofs.I,Dofs_redI.I) = Phi_I;
V(Dofs.I,Dofs_redI.A) = Psi_IA;
V(Dofs.A,Dofs_redI.A) = eye(length(Dofs.A));

% V = full(V);

%% Reduce original UC matrices
M_redI = V.'*M*V;
K_redI = V.'*K*V;
C_redI = V.'*C*V;

%% Symmetrize
M_redI = (M_redI+M_redI.')/2;
K_redI = (K_redI+K_redI.')/2;
C_redI = (C_redI+C_redI.')/2;

end
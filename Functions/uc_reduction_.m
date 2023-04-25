function [M_redI,K_redI, Dofs_redI, V] = uc_reduction_(M,K, Dofs, n_modes, tol_fp, reduction_pgd)

%% Verify
if n_modes > length(Dofs.A)
    error('n_modes_pgd is larger than amount of interface dofs!')
end

%% Determine the input vectors
F = sparse(Dofs.nDOF, numel(Dofs.A));
F(Dofs.A, 1:numel(Dofs.A)) = speye(numel(Dofs.A));

%% Calculate PGD modes
if reduction_pgd
    disp("Solve PGD...")
    figure;
    
    nFP = 40;
    tol_resi = 1e-60;
    freq_pgd = [linspace(50,1000,2000)'];
    
    dv2 = (2*pi*freq_pgd).^2;
    dv0 = ones(size(freq_pgd));
    g0 = ones(size(F,2),1);
    
    freq_func = [];
    freq_dev2 = [];
    spat_func = [];
    mult_func = [];
    
    stat_forc = [];
    iner_forc = [];
    residual = [];
    
    tic
    for iRank = 1:n_modes
        iRank
        flag = 0;
        h_vec = rand(size(freq_pgd));
        h_vec(1) = 100;
        g_vec = rand(size(F,2),1);
        for iter = 1:nFP
            % compute u
            h_vec_old = h_vec;
            H_u = (g_vec'*g_vec)* (h_vec'*h_vec*K - h_vec'*(h_vec.*dv2)*M);
            if isempty(freq_func)
                G_w = h_vec'*dv0* (F*conj(g_vec));
            else
                G_w = h_vec'*dv0* (F*conj(g_vec)) - stat_forc*((h_vec'*freq_func).*(g_vec'*mult_func)).' + iner_forc*((h_vec'*freq_dev2).*(g_vec'*mult_func)).';
            end
            u_vec = H_u\G_w;
            u_vec = u_vec/norm(u_vec+eps);
            % compute g
            g_vec_old = g_vec;
            D_g = h_vec'*h_vec*K - h_vec'*(h_vec.*dv2)*M;
            D_red = u_vec'*D_g*u_vec;
            F_red = u_vec'* (h_vec'*dv0)*F;
            H_g = sparse(1:numel(g_vec), 1:numel(g_vec), D_red*g0);
            if isempty(freq_func)
                G_g = F_red.';
            else
                G_g = F_red.' - mult_func*((h_vec'*freq_func).*(u_vec'*stat_forc)).' + mult_func*((h_vec'*freq_dev2).*(u_vec'*iner_forc)).';
            end
            g_vec = H_g\G_g;
            g_vec = g_vec/norm(g_vec+eps);
            % compute f
            u_vec_old = u_vec;
            K_red = (g_vec'*g_vec)* u_vec'*K*u_vec;
            M_red = (g_vec'*g_vec)* u_vec'*M*u_vec;
            F_red = u_vec'*(F*conj(g_vec));
            dv2 = (2*pi*freq_pgd).^2;
            dv0 = ones(size(freq_pgd));
            H_h = sparse(1:numel(freq_pgd), 1:numel(freq_pgd), dv0*K_red - dv2*M_red);
            if isempty(freq_func)
                G_h = dv0* F_red;
            else
                G_h = dv0* F_red - freq_func*((u_vec'*stat_forc).*(g_vec'*mult_func)).' + freq_dev2*((u_vec'*iner_forc).*(g_vec'*mult_func)).';
            end
            h_vec = H_h\G_h;
            %% check error
            semilogy(abs(freq_pgd), abs(h_vec))
            drawnow;
            error = norm(h_vec_old - h_vec)/norm(h_vec+eps) + norm(u_vec_old - u_vec)/norm(u_vec+eps) + norm(g_vec_old - g_vec)/norm(g_vec+eps)
            if error < tol_fp
                flag = 1;
                break;
            end
        end
        if flag == 0
            disp('warning: fixed point iteration does not converge!')
            freq_func = [freq_func, h_vec];
            freq_dev2 = [freq_dev2, dv2.*h_vec];
            spat_func = [spat_func, u_vec];
            mult_func = [mult_func, g_vec];
            stat_forc = [stat_forc, K*u_vec];
            iner_forc = [iner_forc, M*u_vec];
            continue;
        else
            freq_pgd(1:iRank)
            freq_func = [freq_func, h_vec];
            freq_dev2 = [freq_dev2, dv2.*h_vec];
            spat_func = [spat_func, u_vec];
            mult_func = [mult_func, g_vec];
            stat_forc = [stat_forc, K*u_vec];
            iner_forc = [iner_forc, M*u_vec];
            resi_norm = 1;
            residual = [residual, abs(resi_norm)/(normest(F)*(dv0'*dv0)*(g0'*g0))];
            if residual(end) < tol_resi
                break;
            end
        end
    end
else
    [spat_func,~] = eigs(real(K+K')/2,real(M+M')/2,n_modes,'smallestabs');
end

%% Calculate mode sets I, LR, BT, BL-BR-TR-TL
[Phi_I, nModesI] =  extract_modes(spat_func, Dofs.I);
[Phi_L, nModesL] =  extract_modes(spat_func, Dofs.L);
[Phi_R, nModesR] =  extract_modes(spat_func, Dofs.R);
[Phi_B, nModesB] =  extract_modes(spat_func, Dofs.B);
[Phi_T, nModesT] =  extract_modes(spat_func, Dofs.T);
[Phi_BL, nModesBL] =  extract_modes(spat_func, Dofs.BL);
[Phi_BR, nModesBR] =  extract_modes(spat_func, Dofs.BR);
[Phi_TR, nModesTR] =  extract_modes(spat_func, Dofs.TR);
[Phi_TL, nModesTL] =  extract_modes(spat_func, Dofs.TL);


%% Reorder Dofs indices according to new Dofs.I (n_modes_I) and Dofs.A
Dofs_redI.I = [1:nModesI];
Dofs_redI.L = Dofs_redI.I(end) + [1:nModesL];
Dofs_redI.R = Dofs_redI.L(end) + [1:nModesR];
Dofs_redI.B = Dofs_redI.R(end) + [1:nModesB];
Dofs_redI.T = Dofs_redI.B(end) + [1:nModesT];
Dofs_redI.BL = Dofs_redI.T(end) + [1:nModesBL];
Dofs_redI.BR = Dofs_redI.BL(end) + [1:nModesBR];
Dofs_redI.TR = Dofs_redI.BR(end) + [1:nModesTR];
Dofs_redI.TL = Dofs_redI.TR(end) + [1:nModesTL];
Dofs_redI.A = [Dofs_redI.L, Dofs_redI.R, Dofs_redI.B, Dofs_redI.T, Dofs_redI.BL, Dofs_redI.BR, Dofs_redI.TR, Dofs_redI.TL];
Dofs_redI.nDOF = numel([Dofs_redI.I, Dofs_redI.A]);

%% Create boundary transformation matrix
V = sparse(Dofs.nDOF, Dofs_redI.nDOF);
%I
V(Dofs.I, Dofs_redI.I) = Phi_I;
%LR
V(Dofs.L, Dofs_redI.L) = Phi_L*inv(Phi_L'*Phi_L);
V(Dofs.R, Dofs_redI.R) = Phi_R*inv(Phi_L'*Phi_R);
%BT
V(Dofs.B, Dofs_redI.B) = Phi_B*inv(Phi_B'*Phi_B);
V(Dofs.T, Dofs_redI.T) = Phi_T*inv(Phi_B'*Phi_T);
%BL-BR-TR-TL
V(Dofs.BL, Dofs_redI.BL) = Phi_BL;
V(Dofs.BR, Dofs_redI.BR) = Phi_BR;
V(Dofs.TR, Dofs_redI.TR) = Phi_TR;
V(Dofs.TL, Dofs_redI.TL) = Phi_TL;
V = clean_zeros(V, 1e15);

%% Reduce original UC matrices
M_redI = V'*M*V;
K_redI = V'*K*V;

%% Symmetrize
M_redI = (M_redI+M_redI')/2;
K_redI = (K_redI+K_redI')/2;
M_redI = clean_zeros(M_redI, 1e15);
K_redI = clean_zeros(K_redI, 1e15);
end

function M = clean_zeros(M, thres)
    scalar = max(max(M));
    M(abs(M) < scalar/thres) = 0;
    M = sparse(M);
end

function Phi_orth = compress_modes(Phi)
    if length(Phi(1,:))>=length(Phi(:,1))
        Phi_orth = eye(length(Phi(:,1)));
    else % orthogonalize modeset
        [U,S,~] = svd(Phi);
        Phi_orth = U(:,(diag(S)/S(1,1))>1e-15);
    end
end

function [Phi, nModes] =  extract_modes(spat_func, dofs)
    Phi = spat_func(dofs,:);
    Phi = compress_modes(Phi);
    nModes = size(Phi,2);
end
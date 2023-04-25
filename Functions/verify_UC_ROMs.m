function [f_n_full,f_n_red] = verify_UC_ROMs(K_UC,M_UC,K_UC_red,M_UC_red,n_modes)

% Full UC model
cols_delete_full = ~any(K_UC,1);
K_UC(cols_delete_full,:) = [];
K_UC(:,cols_delete_full) = [];
M_UC(cols_delete_full,:) = [];
M_UC(:,cols_delete_full) = [];
[~,eigs_full] = eigs(K_UC,M_UC,n_modes,0.001);
f_n_full = (diag(eigs_full).^0.5)/(2*pi);

% Reduced UC model
cols_delete_red = ~any(K_UC_red,1);
K_UC_red(cols_delete_red,:) = [];
K_UC_red(:,cols_delete_red) = [];
M_UC_red(cols_delete_red,:) = [];
M_UC_red(:,cols_delete_red) = [];
[~,eigs_red] = eigs(K_UC_red,M_UC_red,n_modes,0.001);
f_n_red = (diag(eigs_red).^0.5)/(2*pi);

%% Plot eigensolutions
figure
subplot(121)
scatter(1:length(eigs_full),real(f_n_full), 300*ones(size(f_n_full)), 'Linewidth', 5)
hold on
scatter(1:length(eigs_full),real(f_n_red), 300*ones(size(f_n_full)), '*', 'Linewidth', 5)
hold off
xlabel('Number');ylabel('Re(f_{eig})')
legend('FOM UC','ROM UC','Location','NorthWest')
ylim([0,3e4])
set(gca, 'FontSize', 40)

subplot(122)
scatter(1:length(eigs_full),imag(f_n_full), 300*ones(size(f_n_full)), 'Linewidth', 5)
hold on
scatter(1:length(eigs_full),imag(f_n_red), 300*ones(size(f_n_full)), '*', 'Linewidth', 5)
hold off
xlabel('Number');ylabel('Im(f_{eig})')
legend('FOM UC','ROM UC','Location','NorthWest')
set(gca, 'FontSize', 40)
ylim([0,3e4])

figure
% scatter(1:length(eigs_full), abs(f_n_full-f_n_red)./abs(f_n_full+1), 300*ones(size(f_n_full)), 'Linewidth', 5)
semilogy(1:length(eigs_full), abs(f_n_full-f_n_red)./abs(f_n_full+1), '*', 'Linewidth', 10)
xlabel('Number');ylabel('Relative error')
title('Relative error of eigen frequencies')
set(gca, 'FontSize', 20)
xlim([7,length(eigs_full)])
ylim([1e-6,1e-2])
end
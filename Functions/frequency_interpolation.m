function uz_range = frequency_interpolation(omega_LR, omega_range, uz_sample)

%% split frequency
[omega_LR, index_sort] = sort(omega_LR);
uz_sample = uz_sample(index_sort);
omega_L = omega_LR(1:2:end);
omega_R = omega_LR(2:2:end);
uz_L = uz_sample(1:2:end);
uz_R = uz_sample(2:2:end);

%% Construct matrices
N = numel(omega_L);
omega_i = omega_L.' * ones(1,N);
omega_j = ones(N,1) * omega_R;
uz_i = uz_L.' * ones(1,N);
uz_j = ones(N,1) * uz_R;
K = (omega_i.^2.*uz_i - omega_j.^2.*uz_j)./ (omega_i.^2 - omega_j.^2);
M = (uz_i - uz_j)./ (omega_i.^2 - omega_j.^2);
cond(K)
cond(M)
uz_range = [];

%% frequency loop
for iFreq = 1:numel(omega_range)
    omega = omega_range(iFreq);
    resp = uz_R* ((K - omega^2*M)\uz_L.');
    uz_range = [uz_range, resp];
end

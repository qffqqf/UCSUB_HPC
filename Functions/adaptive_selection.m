function [new_index, error_max] = adaptive_selection(omega_LR_index, uz1, uz0, nSamp, nStep)

%% Get errors
error = abs(uz1-uz0)./abs(uz1);
error(omega_LR_index) = 0;
[error_n, new_index_] = maxk(error,nSamp*3);
error_max = error_n(1);

new_index_ = [new_index_, randsample(numel(error),50)'];
new_index = zeros(1,nSamp);
jSamp = 0;
for iSamp = 1:numel(new_index_)
    if  sum(abs([omega_LR_index, new_index(1:jSamp)] - new_index_(iSamp)) < nStep) == 0
        jSamp = jSamp + 1;
        new_index(jSamp) = new_index_(iSamp);
    end
    if jSamp == nSamp
        break
    end
end


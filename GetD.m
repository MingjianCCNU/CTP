function [D] = GetD(Cov,V_adj,N_mode)
%% variance of the nullifiers
% var(p_i)
D = diag(Cov((1:2:2*N_mode) + 1, (1:2:2*N_mode) + 1));

for n_d = 1:N_mode
    %     var(x_i) and cov(p_i,x_j)
    for n_i = 1:N_mode
        D(n_d) = D(n_d) + V_adj(n_d, n_i) * Cov(2 * n_i - 1, 2 * n_i - 1)...
            - 2 * V_adj(n_d, n_i) * Cov(2 * n_d, 2 * n_i - 1);
    end
    %     cov(x_i,x_j)
    for n_i = 1:N_mode
        for n_j = n_i+1:N_mode
            D(n_d) = D(n_d) + 2 * V_adj(n_d, n_i) * V_adj(n_d, n_j) *...
                Cov(2 * n_i - 1, 2 * n_j - 1);
        end
    end
end


% D = D./[3; 3; 3; 3];
% for i = 1:2*N_mode
%    if Cov_ini(i, i) < 0
%       D = 2 * ones(1,N_mode);
%       break
%    end
% end
end
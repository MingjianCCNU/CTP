function log_neg = gaussian_entanglement(covariance_matrix)
    % 计算高斯态的log-negativity
    %
    % 参数:
    % covariance_matrix: 4x4协方差矩阵，顺序为[x1, p1, x2, p2]
    %
    % 返回:
    % log_neg: log-negativity值
    
    % 提取子矩阵
    A = covariance_matrix(1:2, 1:2);
    B = covariance_matrix(3:4, 3:4);
    C = covariance_matrix(1:2, 3:4);
    
    % 计算Δ和行列式
    delta = det(A) + det(B) - 2 * det(C);
    det_sigma = det(covariance_matrix);
    
    % 计算部分转置后的最小辛特征值
    nu_tilde_minus = sqrt((delta - sqrt(delta^2 - 4 * det_sigma)) / 2);
    
    % 计算log-negativity
    log_neg = max(-log2(nu_tilde_minus));
end
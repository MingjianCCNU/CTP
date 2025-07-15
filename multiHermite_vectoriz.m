function result = multiHermite_vectoriz(n, x, M)
% MULTIHERMITE_VECTORIZED 向量化计算多维厄米特多项式
%   n: 多重指标向量 [n1, n2, ..., nd]
%   x: 变量张量，维度为 [d, N, N]
%   M: 协方差张量，维度为 [d, d, N, N]
%   返回值: N×N 矩阵，每个元素对应一组 (x(:,:,i,j), M(:,:,i,j)) 的计算结果

    % 检查输入维度
    d = length(n);
    [d_x, N1, N2] = size(x);
    [d_M1, d_M2, N_M1, N_M2] = size(M);
    
    assert(d_x == d, 'x的第一维必须与n的长度一致');
    assert(d_M1 == d && d_M2 == d, 'M的前两维必须为d×d');
    assert(N1 == N2 && N1 == N_M1 && N1 == N_M2, 'x和M的后两维必须为N×N');
    
    % 初始化结果矩阵
    result = zeros(N1, N2);
    
    % 递归计算每个位置的多项式值
    result(:) = compute_multiHermite(n, x, M);
end

function H = compute_multiHermite(n, x, M)
% 递归计算多维厄米特多项式（处理张量输入）

    % 递归终止条件: 全零指标
    if all(n == 0)
        [~, N1, N2] = size(x);
        H = ones(N1, N2);
        return;
    end
    
    % 选择第一个非零的k
    k = find(n > 0, 1, 'first');
    
    % 构造 n - e_k
    n_minus_ek = n;
    n_minus_ek(k) = n(k) - 1;
    n = n_minus_ek;
    
    % 递归计算 H_{n - e_k}(x)
    H_prev = compute_multiHermite(n_minus_ek, x, M);
    
    % 计算主导项: x_k * H_{n - e_k}(x)
    leading_term = reshape(x(k,:,:), size(H_prev)) .* H_prev;
    
    % 计算求和项: 2 * sum_j M_{kj} * n_j * H_{n - e_j}(x)
    sum_term = zeros(size(H_prev));
    
    for j = 1:length(n)
        if n(j) > 0
            % 构造 n - e_j
            n_minus_ej = n;
            n_minus_ej(j) = n(j) - 1;
            
            % 递归计算 H_{n - e_j}(x)
            H_prev_ej = compute_multiHermite(n_minus_ej, x, M);
            
            % 获取 M_{kj} 张量
            M_kj = reshape(M(k,j,:,:), size(H_prev_ej));
            
            sum_term = sum_term + 2 * M_kj .* n(j) * H_prev_ej;
        end
    end
    
    % 最终结果
    H = leading_term + sum_term;
end
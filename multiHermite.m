function result = multiHermite(n, x, M)
% MULTIHERMITE 计算多维厄米特多项式 H_n(x)
%   n: 多重指标向量 [n1, n2, ..., nd]
%   x: 变量向量 [x1, x2, ..., xd]
%   M: 协方差矩阵 (d×d)
%   返回值: H_n(x) 的值

    % 检查输入维度
    d = length(n);
    assert(length(x) == d, 'x的维度必须与n一致');
    assert(size(M,1) == d && size(M,2) == d, 'M必须是d×d的矩阵');
    
    % 初始化记忆化缓存
    persistent cache;
    if isempty(cache)
        cache = containers.Map();
    end
    
    % 生成当前指标的唯一键
    key = sprintf('%d_', n);
    
    % 检查缓存
    if isKey(cache, key)
        result = cache(key);
        return;
    end
    
    % 递归终止条件: 全零指标
    if all(n == 0)
        result = 1;
        cache(key) = result;
        return;
    end
    
    % 选择第一个非零的k
    k = find(n > 0, 1, 'first');
    
    % 构造 n - e_k
    n_minus_ek = n;
    % disp(n_minus_ek)
    n_minus_ek(k) = n(k) - 1;
    n = n_minus_ek;
    
    % 递归计算 H_{n - e_k}(x)
    H_prev = multiHermite(n_minus_ek, x, M);
    
    % 计算主导项: x_k * H_{n - e_k}(x)
    leading_term = x(k) * H_prev;
    
    % 计算求和项: 2 * sum_j M_{kj} * n_j * H_{n - e_j}(x)
    sum_term = 0;
    for j = 1:d
        if n(j) > 0
            % 构造 n - e_j
            n_minus_ej = n;
            n_minus_ej(j) = n(j) - 1;
            
            % 递归计算 H_{n - e_j}(x)
            H_prev_ej = multiHermite(n_minus_ej, x, M);
            
            sum_term = sum_term + 2 * M(k,j) * n(j) * H_prev_ej;
        end
    end
    
    % 最终结果
    result = leading_term + sum_term;
    
    % 存入缓存
    cache(key) = result;
end